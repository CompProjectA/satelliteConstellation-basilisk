import os
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.simulation import spacecraft, spacecraftLocation, simSynch,extForceTorque,simpleNav
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion, simIncludeGravBody, unitTestSupport, vizSupport)
from Basilisk import __path__
from Basilisk.architecture import messaging
from Basilisk.utilities import SimulationBaseClass

bskPath = __path__[0]

# Define a simple message payload for strings
from Basilisk.architecture import messaging
from Basilisk.architecture import sysModel

class MessageReceiver(sysModel.SysModel):
    def __init__(self):
        super(MessageReceiver, self).__init__()
        self.inputMsg = messaging.DoubleMsg()
        self.inputMsgName = "satMessage"
    
    def Reset(self, currentTime):
        print("Receiver reset at time:", currentTime)

    def UpdateState(self, currentTime):
        if self.inputMsg.isLinked():
            msg = self.inputMsg()
            print(f"[t={currentTime*1e-9:.2f}s] Debris received message: {msg.data}")

# class SatMessage(messaging.MessagePayload):
#     def __init__(self):
#         super(SatMessage, self).__init__()
#         self.value = 0.0  # Example field

def run(show_plots=True):
    print("Basilisk Path:", bskPath)

    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()

    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

    
    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    earth.isCentralBody = True
    mu = earth.mu
    earth_radius = earth.radEquator
    
    num_clusters = 2
    sats_per_cluster = 4
    total_sats = num_clusters * sats_per_cluster
    satellites = []
    relative_anomalies = [0, 5, -5, 10]  # Leader + close followers
    cluster_RAAC = [0, 0, 90, 270]       # RAAN for each cluster
    cluster_anomaly_offset = [0, 180, 0, 0]  # Offset for anomaly (cluster 2 is opposite side)
    a_base = 7000e3
    mu = 398600.4418 * 1e9
    for i in range(num_clusters):
        for j in range(sats_per_cluster):
            sat = spacecraft.Spacecraft()
            sat.ModelTag = f"Satellite{i+1}"
            satellites.append(sat)
            scSim.AddModelToTask(simTaskName, sat)

    # Add the sun
    # sun = gravFactory.createSun()
    # sun.isCentralBody = False  # Ensure the sun is not the central body for the orbits of the satellites

    for sat in satellites:
        gravFactory.addBodiesTo(sat)
    # extFTObject = extForceTorque.ExtForceTorque()
    # extFTObject.ModelTag = "extTorque"
    # satellites[0].addDynamicEffector(extFTObject)
    # scSim.AddModelToTask(simTaskName, extFTObject)
    # sNavObject = simpleNav.SimpleNav()
    # sNavObject.ModelTag = "SimpleNavigation"
    # sNavObject.scStateInMsg.subscribeTo(satellites[0].scStateOutMsg)
    # scSim.AddModelToTask(simTaskName, sNavObject)
    # Define a base orbit
    rLEO_base = 7000. * 1000  # 7000 km
   
    # Define the angular speed based on the base velocity and position
    # angular_speed=vN_base/ np.linalg.norm(rN_base)
    # print(rN_base)
    # Formation offsets relative to the leader - TIGHT FORMATION
    relative_offsets = [
        np.array([0.0, 0.0, 0.0]),       # Leader (Sat1)
        np.array([200e3, 0.0, 0.0]),     # Child (Sat2) - Very close
        np.array([400e3,0.0 , 0.0]),     # Child (Sat3) - Very close
        np.array([600e3, 0.0, 0.0])    # Child (Sat4) - Very close
    ]

    for i, sat in enumerate(satellites):
        # offset = relative_offsets[i]
        # r_init_N = rN_base + offset
        # sat.hub.r_CN_NInit = r_init_N

        # # Approximate circular orbital speed at this radius
        # if(i==0):
        #     sat.hub.v_CN_NInit = vN_base
        # else:
        #     r_mag = np.linalg.norm(r_init_N)
        #     v_mag = np.sqrt(mu / r_mag)
        #     orbit_normal_N = np.cross(rN_base, vN_base)
        #     orbit_normal_N = orbit_normal_N / np.linalg.norm(orbit_normal_N)
        # # Calculate a velocity vector perpendicular to the position vector
        #     velocity_direction = np.cross(orbit_normal_N, r_init_N)
        #     velocity_direction = velocity_direction / np.linalg.norm(velocity_direction)
        #     v_init_N = velocity_direction * v_mag
        #     sat.hub.v_CN_NInit = v_init_N
        if(i >=0 and i < 4):
            # Set the leader's orbit
            oe = orbitalMotion.ClassicElements()
            oe.a = rLEO_base + i * 200e3 # Different altitudes
            oe.e = 0.0  # Circular
            oe.i = 55 * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = 45 * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            sat.hub.r_CN_NInit = rN
            sat.hub.v_CN_NInit = vN
            sat.hub.sigma_BNInit = [0.0, 0.0, 0.0]
            sat.hub.omega_BN_BInit = [0.0, 0.0, 0.0]
        else:
            # Opposite side 
            oe = orbitalMotion.ClassicElements()
            oe.a = rLEO_base + (i-4) * 200e3
            oe.e = 0.0  # Circular
            oe.i = 55 * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = (45+180) * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            sat.hub.r_CN_NInit = rN
            sat.hub.v_CN_NInit = vN
            sat.hub.sigma_BNInit = [0.0, 0.0, 0.0]
            sat.hub.omega_BN_BInit = [0.0, 0.0, 0.0]
    # Set up communication modules
    print(satellites)
    scLocationModules = []
    accessRecorders = []
    maximum_range = 2000e3  # Increased maximum range to 5000 km
    print(f"Maximum Range: {maximum_range/1000} km")
   
    # leader_pos = satellites[0].scStateOutMsg.read().r_BN_N
    # follower_pos = [sat.scStateOutMsg.read().r_BN_N for sat in satellites[1:]]
    # avg_dir = np.mean([np.array(pos) - np.array(leader_pos) for pos in follower_pos], axis=0)
    # aHat_B = avg_dir / np.linalg.norm(avg_dir)

    scLocation = spacecraftLocation.SpacecraftLocation()
    scLocation.ModelTag = f"CommCheck{i+1}"
    scLocation.rEquator = earth_radius
    scLocation.rPolar = earth_radius * 0.98
    scLocation.aHat_B = [0, 1.0, 0]  # Pointing in the +Y direction
    scLocation.theta = np.radians(30.0)  # Adjust as needed to be narrower
    scLocation.maximumRange = maximum_range
    for i, sat in enumerate(satellites):
        if i == 0:
            scLocation.primaryScStateInMsg.subscribeTo(sat.scStateOutMsg) # Subscribe all scLocations to their respective satellite's state
        else:
            scLocation.addSpacecraftToModel(sat.scStateOutMsg)
        # Create recorders for all relevant pairs
        

    scLocationModules.append(scLocation)
    scSim.AddModelToTask(simTaskName, scLocation)
    # Force initial update of spacecraft location modules
    # for scLocation in scLocationModules:
    #     scLocation.UpdateState(macros.sec2nano(0.0))
    for j in range(total_sats-1):
            recorder_name = f"AccessRecorder_{i+1}_{j+1}"
            recorder = scLocation.accessOutMsgs[j].recorder()
            # sending message to the debris
            # servicerMessage = messaging.DoubleMsg()
            # servicerMessagePayload = messaging.DoubleMsgPayload()
            # servicerMessagePayload.data = 42.0  # Example data
            # servicerMsg = servicerMessage.write(servicerMessagePayload)
            # receiver = MessageReceiver()
            # receiver.inputMsg.subscribeTo(servicerMsg)
            # scSim.AddModelToTask(simTaskName, receiver)
            recorder.ModelTag = recorder_name
            accessRecorders.append(recorder)
            scSim.AddModelToTask(simTaskName, recorder)
    # Visualization
    if vizSupport.vizFound:
        print("Vizard Visualization Found. Enabling live streaming...")
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)
        spriteList = ["satellite"] * 8
        viz = vizSupport.enableUnityVisualization(scSim,
                                                 simTaskName,
                                                 satellites,
                                                 liveStream=True,
                                                 )
        # Add a rotating frame for visualization
        # vizSupport.setCameraTarget(viz, [0, 0, 0])  # Set the camera target
        # vizSupport.setCameraMode(viz, "velocityRelative")  # Use velocity-relative mode
        # vizSupport.addCelestialBody(viz, sunName='Sun',
        #                                r_SB_N=[0, 0, 0],
        #                                )

        for i, scLocation in enumerate(scLocationModules):
            vizSupport.addLocation(
                viz,
                stationName=f"Satellite{i + 1}",
                parentBodyName=f"Satellite{i + 1}",
                r_GP_P=[0, 0, 0],
                gHat_P=[0, 1, 0],
                fieldOfView=np.pi / 4,
                range=scLocation.maximumRange
            )

    else:
        print("Vizard Visualization Module Not Found. Check Basilisk installation.")

    simulationTime = macros.sec2nano(2400.0)  # Simulate for 120 seconds
    scSim.InitializeSimulation()
    print(scLocation.accessOutMsgs[0].recorder().hasAccess)
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()

    # Print debugging information: distances
    print("\nSatellite Distances:")
    for i, scLocation in enumerate(scLocationModules):
        if i == 0:  # Leader
            for j in range(1, num_satellites):
                dist = np.linalg.norm(np.array(satellites[i].hub.r_CN_NInit) - np.array(satellites[j].hub.r_CN_NInit))
                print(f"Sat{i+1} to Sat{j+1}: {dist/1000} km")
        elif i > 0:
             dist = np.linalg.norm(np.array(satellites[i].hub.r_CN_NInit) - np.array(satellites[0].hub.r_CN_NInit))
             print(f"Sat{i+1} to Sat{0+1}: {dist/1000} km")
    print("-" * 30)

    # Plotting access
    if show_plots:
        time_vec = accessRecorders[0].times() * macros.NANO2MIN
        plt.figure(figsize=(10, 6))

        
        for i in range(1,num_satellites):
                    # Plot access from leader to others and from others to leader
                    access_data = accessRecorders[i-1].hasAccess
                    label = f'Sat1 to Sat{i+1}'
                    plt.plot(time_vec, access_data, label=label)

        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader (Sat1) and Children (Sat2-4)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    run(show_plots=True)
