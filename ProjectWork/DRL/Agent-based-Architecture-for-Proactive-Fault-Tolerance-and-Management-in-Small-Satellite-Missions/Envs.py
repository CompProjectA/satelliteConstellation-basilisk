from Basilisk.architecture import messaging
import numpy as np
from Basilisk import __path__
from Basilisk.simulation import reactionWheelStateEffector, simpleBattery, spacecraft, extForceTorque, ReactionWheelPower
from Basilisk.utilities import RigidBodyKinematics, unitTestSupport, macros, simIncludeRW, SimulationBaseClass, vizSupport
import gymnasium as gym
from gymnasium import spaces
import  os
import pandas as pd  # Import pandas to handle the data logging

bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])


class BasiliskModel:
    def __init__(self, I: list, ref_MRP: np.ndarray, torque_mode: str = "wheel", timestep: int = macros.sec2nano(0.01)) -> None:
        self.ref_MRP = ref_MRP
        self.torque_mode = torque_mode
        simulationTime = macros.sec2nano(8000)
        self.timeData=None
        self.step_count = 0
        self.fault_time = 2000#np.random.randint(low=2000, high=7000)
        self.wheel_num = 0#np.random.randint(low=0, high=3)
        self.faulty=True
        # Initialize the simulation environment
        self.scSim = SimulationBaseClass.SimBaseClass()
        task_name = "sim_task"
        process_name = "sim_process"
        dyn_process = self.scSim.CreateNewProcess(process_name)
        dyn_process.addTask(self.scSim.CreateNewTask(task_name, timestep))
        
        # Initialize the spacecraft object
        self.scObject = spacecraft.Spacecraft()
        self.scObject.hub.r_BcB_B = [[0.0], [0.0], [0.0]]
        self.scObject.hub.IHubPntBc_B = unitTestSupport.np2EigenMatrix3d(I)
        self.scObject.hub.sigma_BNInit = np.zeros(3, dtype=np.float32)
        self.scObject.hub.omega_BN_BInit = np.zeros(3, dtype=np.float32)
        self.scSim.AddModelToTask(task_name, self.scObject, None)

        # Initialize the reaction wheels if using "wheel" torque mode
        if torque_mode == "wheel":
            self.map_matrix = np.array(
                [
                    [1, 0, 0, np.sqrt(3) / 3],
                    [0, 1, 0, np.sqrt(3) / 3],
                    [0, 0, 1, np.sqrt(3) / 3],
                ]
            )
            self.rwFactory = simIncludeRW.rwFactory()
            varRWModel = messaging.BalancedWheels
            RW1 = self.rwFactory.create("Honeywell_HR16", [1, 0, 0], maxMomentum=50.0, Omega=0.0, RWModel=varRWModel)
            RW2 = self.rwFactory.create("Honeywell_HR16", [0, 1, 0], maxMomentum=50.0, Omega=0.0, RWModel=varRWModel)
            RW3 = self.rwFactory.create("Honeywell_HR16", [0, 0, 1], maxMomentum=50.0, Omega=0.0, RWModel=varRWModel)
            RW4 = self.rwFactory.create(
                "Honeywell_HR16", [np.sqrt(3) / 3, np.sqrt(3) / 3, np.sqrt(3) / 3], maxMomentum=50.0, Omega=0.0, RWModel=varRWModel
            )
            self.rwStateEffector = reactionWheelStateEffector.ReactionWheelStateEffector()
            self.rwStateEffector.ModelTag = "RW_cluster"
            self.rwFactory.addToSpacecraft(self.scObject.ModelTag, self.rwStateEffector, self.scObject)
            self.scSim.AddModelToTask(task_name, self.rwStateEffector, None)
            self.wheelTorqueMsg = messaging.ArrayMotorTorqueMsg()
            self.rwStateEffector.rwMotorCmdInMsg.subscribeTo(self.wheelTorqueMsg)
            
            # 1- Initialize reaction wheel power modules
            self.rwPowerList = []
            for c in range(4):
                powerRW = ReactionWheelPower.ReactionWheelPower()
                powerRW.ModelTag = self.scObject.ModelTag + "RWPower" + str(c)
                powerRW.basePowerNeed = 5.0  # baseline power draw, Watts
                powerRW.rwStateInMsg.subscribeTo(self.rwStateEffector.rwOutMsgs[c])
                self.scSim.AddModelToTask(task_name, powerRW)
                self.rwPowerList.append(powerRW)
            # 2-  create battery module
            battery = simpleBattery.SimpleBattery()
            battery.ModelTag =  self.scObject.ModelTag
            battery.storageCapacity = 300000  # W-s
            battery.storedCharge_Init = battery.storageCapacity * 0.8  # 20% depletion
            self.scSim.AddModelToTask(task_name, battery)
            # 3- connect RW power to the battery module
            for c in range(4):
                battery.addPowerNodeToModel(self.rwPowerList[c].nodePowerOutMsg)
           # 4- To log the RW information, the following code is used:
            numDataPoints = 100000
            samplingTime = unitTestSupport.samplingTime(simulationTime, timestep, numDataPoints)
            print(samplingTime)
            self.rwOutLog = []
            self.rwPowLog = []
            for c in range(4):
                self.rwOutLog.append(self.rwStateEffector.rwOutMsgs[c].recorder(samplingTime))
                self.rwPowLog.append(self.rwPowerList[c].nodePowerOutMsg.recorder(samplingTime))
                self.scSim.AddModelToTask(task_name, self.rwOutLog[-1])
                self.scSim.AddModelToTask(task_name, self.rwPowLog[-1])
            self.batPowLog = battery.batPowerOutMsg.recorder(samplingTime)
            self.scSim.AddModelToTask(task_name, self.batPowLog)
            self.rwSpeedLog = self.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
            self.scSim.AddModelToTask(task_name, self.rwSpeedLog)
        # Initialize external force/torque object for direct torque application
        self.extFTObject = extForceTorque.ExtForceTorque()
        self.extFTObject.extTorquePntB_B = [[0.0], [0.0], [0.0]]
        self.scObject.addDynamicEffector(self.extFTObject)

        # Initialize the spacecraft state
        self.scObject.SelfInit()
        self.scObject.initializeDynamics()
        vizSupport.enableUnityVisualization(self.scSim, task_name, self.scObject,  saveFile=fileName, rwEffectorList=self.rwStateEffector)
        self.scSim.InitializeSimulation()

        # Initialize state variables
        self.cur_MRP = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        self.cur_error_MRP = RigidBodyKinematics.subMRP(self.ref_MRP, self.cur_MRP)  # Q = subMRP(Q1,Q2) from Q2 to Q1.
        self.cur_error_angle = 4 * np.arctan(np.linalg.norm(self.cur_error_MRP))
        self.cur_omega = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        self.cur_omega_dot = np.array([0.0, 0.0, 0.0], dtype=np.float32)

        # History tracking
        self.MRP_history = [self.cur_MRP]
        self.error_MRP_history = [self.cur_error_MRP]
        self.error_angle_history = [self.cur_error_angle]
        self.omega_history = [self.cur_omega]
        self.omega_dot_history = [self.cur_omega_dot]
        self.rw_power_history = [[] for _ in range(4)]  # Initialize power history for each RW

        self.dataRW = [[] for _ in range(4)]
        self.dataRwPower = [[] for _ in range(4)]
       
        self.x=0


       

    def step(self, cur_nano_second: int, torque: np.ndarray):
        self.step_count+=1

        
        if self.torque_mode == "wheel":
           
            wheelTorqueBuffer = messaging.ArrayMotorTorqueMsgPayload()
            
            wheelTorqueBuffer.motorTorque = torque
            
            self.wheelTorqueMsg.write(wheelTorqueBuffer)
            
        elif self.torque_mode == "axis":
            self.extFTObject.extTorquePntB_B = torque
        else:
            raise ValueError("torque_mode must be 'wheel' or 'axis'")
        
        #if cur_nano_second>=3:
         #   torque[0]=0
    
        

        self.scSim.ConfigureStopTime(cur_nano_second)
        self.scSim.ExecuteSimulation()

        state = self.scObject.scStateOutMsg.read()
        self.cur_MRP = np.array(state.sigma_BN, dtype=np.float32)
        self.cur_error_MRP = RigidBodyKinematics.subMRP(self.ref_MRP, self.cur_MRP)
        self.cur_error_angle = 4 * np.arctan(np.linalg.norm(self.cur_error_MRP))
        self.cur_omega = np.array(state.omega_BN_B, dtype=np.float32)
        self.cur_omega_dot = np.array(state.omegaDot_BN_B, dtype=np.float32)
       
        self.MRP_history.append(self.cur_MRP)
        self.error_MRP_history.append(self.cur_error_MRP)
        self.error_angle_history.append(self.cur_error_angle)
        self.omega_history.append(self.cur_omega)
        self.omega_dot_history.append(self.cur_omega_dot)
        
        self.dataOmegaRW = self.rwSpeedLog.wheelSpeeds
        #print('length self.dataOmegaRW = ', len(self.dataOmegaRW))
        #for c in range(0, 4):
        #    self.dataRW[c].append(self.rwOutLog[c].u_current)
        #    self.dataRwPower[c].append(self.rwPowLog[c].netPower)
            #print(self.rwPowLog[c].netPower)
        batteryStorageLog = self.batPowLog.storageLevel
        currentNetPower = self.batPowLog.currentNetPower
        self.timeData = self.rwOutLog[0].times() 
        #print(len(self.timeData))
        #print('batteryStorageLog ' , batteryStorageLog)
        #print('currentNetPower ' , currentNetPower)
        #self.x+=1
       # while self.x>3:
        #    x=0
        # Log RW power consumption
        #for c in range(0,4):
            #self.rw_power_history[c].append(self.rwPowLog[c].netPower[0])
            
            #print(self.rwPowLog[c].netPower[0])

        return self.cur_MRP, self.cur_error_MRP, self.cur_error_angle, self.cur_omega, self.cur_omega_dot

    
import Tools


class BasiliskEnv(gym.Env):
    metadata = {"render_modes": ["console"]}

    def __init__(self, render_mode="console", faulty: bool = False, torque_mode: str = "wheel"):
        super(BasiliskEnv, self).__init__()
        self.faulty = faulty
        self.torque_mode = torque_mode
        self.render_mode = render_mode
        self.model: BasiliskModel = None
        self.ref_MRP = np.array([0.0, 0.0, 0.0])
        self.reward = 0
        self.action = np.zeros(4, dtype=np.float32)
        self.observation = np.zeros(6, dtype=np.float32)
        self.step_count = 0
        self.fault_time = 0
        self.wheel_num = 0
        self.torque_history = [np.array([0, 0, 0,0])]

        # Define observation space as a Dict space
        self.observation_space = spaces.Dict({
            "observation": spaces.Box(low=-np.inf, high=np.inf, shape=(6,), dtype=np.float32),
            "desired_goal": spaces.Box(low=-np.inf, high=np.inf, shape=(3,), dtype=np.float32),
            "achieved_goal": spaces.Box(low=-np.inf, high=np.inf, shape=(3,), dtype=np.float32)
        })

        if self.torque_mode == "wheel":
            self.action_space = spaces.Box(low=-1, high=1, shape=(4,), dtype=np.float32)
        elif self.torque_mode == "axis":
            self.action_space = spaces.Box(low=-1, high=1, shape=(3,), dtype=np.float32)
        else:
            raise ValueError("torque_mode must be 'wheel' or 'axis'")

        self.reset()

    def reset(self, seed=None, options=None):
        super().reset(seed=seed, options=options)
        self.ref_MRP = RigidBodyKinematics.euler1232MRP(Tools.random_euler())
        self.model = BasiliskModel(I=[0.025, 0, 0, 0, 0.05, 0, 0, 0, 0.065], ref_MRP=self.ref_MRP, torque_mode=self.torque_mode)
        self.reward = 0
        self.action = np.zeros(4, dtype=np.float32)  # Four-wheel torque
        self.observation = np.zeros(6, dtype=np.float32)  # Error MRP, angular velocity
        self.step_count = 0
        self.fault_time =  3000#np.random.randint(low=0, high=3000)
        self.wheel_num =  0#np.random.randint(low=-1, high=4)

        self.observation = np.hstack([self.model.cur_error_MRP, self.model.cur_omega], dtype=np.float32)

        # Returning a dict for HER
        #return {
          #  "observation": self.observation,
         #   "desired_goal": self.ref_MRP,
        #    "achieved_goal": self.model.cur_MRP
       # }, {}
        return self.observation, {}

    def step(self, action):
        terminated = False
        truncated = False
        self.step_count += 1
        self.action = 0.01 * action  

        if self.faulty:
            if self.step_count > self.fault_time and self.wheel_num != -1:
                self.action[self.wheel_num] = 0
                # wheelSpeeds Message
                rwSpeedMessage = messaging.RWSpeedMsgPayload()
                Omega = [10.0, 25.0, 50.0, 100.0]  # rad/sec
                rwSpeedMessage.wheelSpeeds = Omega
                rwSpeedInMsg = messaging.RWSpeedMsg().write(rwSpeedMessage)
                #self.model.rwStateEffector.rwSpeedOutMsg.read().wheelSpeeds[0] = 0.0
                #print('length self.model.rwStateEffector.rwSpeedOutMsg.read().wheelSpeeds',len(self.model.rwStateEffector.rwSpeedOutMsg.read().wheelSpeeds))

                #print('self.model.rwStateEffector.rwSpeedOutMsg.read().wheelSpeeds',self.model.rwStateEffector.rwSpeedOutMsg.read().wheelSpeeds)
            else:
                self.action[3]=0.0
                #print('action =', self.action)

        self.model.step(macros.sec2nano((self.step_count) * 0.01), self.action)
        
        self.torque_history.append(self.action)
        self.observation = np.hstack([self.model.cur_error_MRP, self.model.cur_omega], dtype=np.float32)

        self.reward = self.calculate_reward()

        if self.step_count == 6000:
            truncated = True
        if np.abs(self.model.cur_omega[0]) > 1 or np.abs(self.model.cur_omega[1]) > 1 or np.abs(self.model.cur_omega[2]) > 1:
            terminated = True

        info = {}
       # return {
           # "observation": self.observation,
          #  "desired_goal": self.ref_MRP,
         #   "achieved_goal": self.model.cur_MRP
        #}, self.reward, terminated, truncated, info
        return self.observation, self.reward, terminated, truncated, info

    def calculate_reward(self):
        pre_error = self.model.error_angle_history[-2]
        cur_error = self.model.error_angle_history[-1]

        r1 = (10 - cur_error) * (pre_error - cur_error) / np.pi

        r2 = 0
        if np.abs(self.model.cur_omega[0]) > 1 or np.abs(self.model.cur_omega[1]) > 1 or np.abs(self.model.cur_omega[2]) > 1:
            r2 = -1

        r3 = 0
        if cur_error < 0.0043633 and np.linalg.norm(self.model.cur_omega) < 0.001:
            r3 = 1

        reward = r1 + r2 + r3
        # print(f"r1: {r1}, r2: {r2}")
        return reward

    def render(self):
        pass

    def close(self):
        pass





