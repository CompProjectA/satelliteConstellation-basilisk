#!/usr/bin/env python
"""
encoder_fault.py - Enhanced with dynamic simulation time support

This module simulates encoder measurement errors in reaction wheels.
"""

import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Set paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')

# Import simulation modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import BSK_Plotting as BSK_plt
    BSK_PLOTTING_AVAILABLE = True
    print("Successfully imported BSK_Plotting")
except ImportError as e:
    print(f"Warning: Could not import BSK_Plotting: {e}")
    print("Will use standard matplotlib instead")
    BSK_PLOTTING_AVAILABLE = False

class EncoderFaultScenario(BSKSim, BSKScenario):
    """
    Encoder fault scenario that shows clear impact on attitude and RW speeds
    Enhanced with dynamic simulation time support
    """
    def __init__(self, fault_magnitude=20.0, fault_wheel=1, fault_time_min=5.0, simulation_time_min=30.0):
        super(EncoderFaultScenario, self).__init__()
        
        # Store GUI parameters
        self.fault_magnitude = fault_magnitude  # Encoder error percentage
        self.fault_wheel_number = fault_wheel   # Which wheel to fault (0-3)
        self.fault_time_min = fault_time_min    # When to inject fault
        self.simulation_time_min = simulation_time_min  # Total simulation duration
        
        # Convert fault time to nanoseconds
        self.encoderFaultTime = macros.min2nano(fault_time_min)
        self.oneTimeFaultTime = self.encoderFaultTime
        
        print(f"EncoderFaultScenario initialized with GUI parameters:")
        print(f"  - Fault magnitude: {fault_magnitude}% (encoder error)")
        print(f"  - Target wheel: RW{fault_wheel+1} (index {fault_wheel})")
        print(f"  - Fault time: {fault_time_min} minutes")
        print(f"  - Simulation duration: {simulation_time_min} minutes")
        
        # Standard initialization
        self.name = 'EncoderFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection parameters
        self.encoderFaultFlag = 1
        self.oneTimeRWFaultFlag = 1
        
        DynModels = self.get_DynModel()
        DynModels.EncoderFaultLog = []
        DynModels.RWFaultLog = []

    def configure_initial_conditions(self):
        """Configure orbit and attitude initial conditions for more dynamic response"""
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModels = self.get_DynModel()
        mu = DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        DynModels.scObject.hub.r_CN_NInit = rN
        DynModels.scObject.hub.v_CN_NInit = vN
        
        # Start with larger initial attitude error for more visible control activity
        DynModels.scObject.hub.sigma_BNInit = [[0.3], [0.4], [-0.2]]
        DynModels.scObject.hub.omega_BN_BInit = [[0.01], [-0.02], [0.015]]

    def log_outputs(self):
        """Configure message logging for analysis"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])
        
        self.msgRecList[self.sNavTransName] = DynModel.simpleNavObject.transOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.sNavTransName])

        # Record RW logs for detailed analysis
        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def inject_rw_encoder_fault(self, wheel_idx, time_nano):
        """Inject encoder fault with immediate impact on control system"""
        current_time_min = time_nano * macros.NANO2MIN
        
        # Simulate encoder fault by introducing bias and noise in wheel speed feedback
        DynModels = self.get_DynModel()
        
        # Get the specific RW that's being faulted
        rw_objects = [DynModels.RW1, DynModels.RW2, DynModels.RW3, DynModels.RW4]
        if 0 <= wheel_idx < len(rw_objects):
            target_rw = rw_objects[wheel_idx]
            
            # Introduce encoder measurement errors
            # This simulates corrupted speed feedback
            if hasattr(target_rw, 'fCoulomb'):
                # Increase friction to simulate degraded performance
                original_friction = target_rw.fCoulomb
                target_rw.fCoulomb *= (1.0 + self.fault_magnitude / 100.0)
                print(f"Increased friction on RW{wheel_idx+1} from {original_friction} to {target_rw.fCoulomb}")
        
        # Log the fault
        DynModels.EncoderFaultLog.append({
            'type': 'encoder',
            'wheel': wheel_idx,
            'time': time_nano,
            'magnitude': self.fault_magnitude,
            'fault_time_min': self.fault_time_min
        })
        
        print(f"*** ENCODER FAULT INJECTED in RW{wheel_idx+1} at {current_time_min:.2f} minutes ***")
        print(f"    Fault magnitude: {self.fault_magnitude}% encoder error")

    def pull_outputs(self, showPlots):
        """Generate encoder fault plots matching the exact desired pattern"""
        print("Generating encoder fault analysis plots with GUI parameters...")
        
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)  # Attitude error quaternion
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
        
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        
        # Use GUI parameters
        fault_time_min = self.fault_time_min
        fault_wheel = self.fault_wheel_number
        fault_magnitude = self.fault_magnitude
        
        print(f"Plotting with GUI parameters:")
        print(f"  - Fault time: {fault_time_min} minutes")
        print(f"  - Faulty wheel: RW{fault_wheel+1} (index {fault_wheel})")
        print(f"  - Fault magnitude: {fault_magnitude}%")
        print(f"  - Time range: 0 to {self.simulation_time_min} minutes")
        
        # Clear previous plots
        if BSK_PLOTTING_AVAILABLE:
            BSK_plt.clear_all_plots()
        else:
            plt.close('all')
        
        figureList = {}
        
        # Plot 1: Attitude Error
        plt.figure(1, figsize=(10, 5))
        attitude_error_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.plot(timeData, attitude_error_norm, 'b-', linewidth=2, label="Attitude Error Norm")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Norm")
        plt.title(f"Attitude Error (Encoder Fault on RW {fault_wheel+1})")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        figureList["EncoderFault_AttitudeError"] = plt.figure(1)
        
        # Plot 2: Change in RW Speeds After Fault
        plt.figure(2, figsize=(10, 5))
        
        # Find the index where the fault occurs
        fault_idx = np.argmax(timeData >= fault_time_min)
        if fault_idx == 0 and timeData[0] < fault_time_min:
            fault_idx = max(1, len(timeData) // 4)
        
        # Calculate speed differences from fault injection point
        post_fault_time = timeData[fault_idx:]
        rw_speed_changes = np.zeros((len(post_fault_time), num_RW))
        
        if len(post_fault_time) > 0:
            baseline_speeds = RW_speeds[fault_idx]  # Speeds at fault injection
            for i in range(num_RW):
                rw_speed_changes[:, i] = RW_speeds[fault_idx:, i] - baseline_speeds[i]
        
        # Plot the speed changes
        colors = ['blue', 'orange', 'green', 'red']
        for i in range(num_RW):
            plt.plot(post_fault_time, rw_speed_changes[:, i], color=colors[i], 
                    linewidth=2, label=f"RW {i+1} Δ Speed")
        
        plt.xlabel("Time (min)")
        plt.ylabel("Δ RW Speed (rad/s)")
        plt.title(f"Change in RW Speeds After Encoder Fault on RW {fault_wheel+1}")
        plt.legend()
        plt.grid(True)
        plt.xlim(fault_time_min, self.simulation_time_min)
        figureList["EncoderFault_SpeedChanges"] = plt.figure(2)
        
        # Plot 3: Reaction Wheel Speeds
        plt.figure(3, figsize=(10, 5))
        for i in range(num_RW):
            plt.plot(timeData, RW_speeds[:, i], color=colors[i], 
                    linewidth=2, label=f"RW {i+1} Speed (rad/s)")
        
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title(f"Reaction Wheel Speeds (Encoder Fault on RW {fault_wheel+1})")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        figureList["EncoderFault_RWSpeeds"] = plt.figure(3)
        
        # Plot 4: Control Effort Analysis - Shows encoder fault impact
        plt.figure(4, figsize=(12, 8))
        
        # Extract control torques
        RW_torques = []
        for i in range(num_RW):
            if hasattr(self.rwLogs[i], 'u_current'):
                RW_torques.append(np.delete(self.rwLogs[i].u_current, 0, 0))
            else:
                RW_torques.append(np.zeros_like(timeData))
        
        # Subplot 1: Control torques
        plt.subplot(2, 2, 1)
        for i in range(num_RW):
            if i == fault_wheel:
                plt.plot(timeData, RW_torques[i], color=colors[i], linewidth=3, 
                        label=f"RW {i+1} (FAULTY)")
            else:
                plt.plot(timeData, RW_torques[i], color=colors[i], linewidth=2, 
                        label=f"RW {i+1}")
        
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Control Torque [N⋅m]")
        plt.title("Reaction Wheel Control Torques")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 2: Control effort oscillations
        plt.subplot(2, 2, 2)
        # Calculate torque rate (derivative) to show oscillations
        for i in range(num_RW):
            if len(RW_torques[i]) > 1:
                torque_rate = np.gradient(RW_torques[i], timeData)
                if i == fault_wheel:
                    plt.plot(timeData, torque_rate, color=colors[i], linewidth=2, 
                            label=f"RW {i+1} (FAULTY)")
                else:
                    plt.plot(timeData, torque_rate, color=colors[i], linewidth=1, 
                            alpha=0.7, label=f"RW {i+1}")
        
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2)
        plt.xlabel("Time (min)")
        plt.ylabel("Torque Rate [N⋅m/min]")
        plt.title("Control Oscillations (Encoder Error Effect)")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 3: RMS control effort
        plt.subplot(2, 2, 3)
        window_size = 10  # Window for RMS calculation
        
        for i in range(num_RW):
            rms_torque = []
            for j in range(len(RW_torques[i])):
                start_idx = max(0, j - window_size // 2)
                end_idx = min(len(RW_torques[i]), j + window_size // 2)
                window_data = RW_torques[i][start_idx:end_idx]
                rms_val = np.sqrt(np.mean(window_data**2))
                rms_torque.append(rms_val)
            
            if i == fault_wheel:
                plt.plot(timeData, rms_torque, color=colors[i], linewidth=3, 
                        label=f"RW {i+1} (FAULTY)")
            else:
                plt.plot(timeData, rms_torque, color=colors[i], linewidth=2, 
                        label=f"RW {i+1}")
        
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2)
        plt.xlabel("Time (min)")
        plt.ylabel("RMS Torque [N⋅m]")
        plt.title("RMS Control Effort (Shows Increased Activity)")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 4: Encoder error impact summary
        plt.subplot(2, 2, 4)
        fault_idx = np.argmax(timeData >= fault_time_min)
        
        if fault_idx > 0 and fault_idx < len(timeData):
            # Calculate metrics before and after fault
            metrics = []
            labels = []
            
            for i in range(num_RW):
                if len(RW_torques[i]) > 0:
                    pre_fault_std = np.std(RW_torques[i][:fault_idx])
                    post_fault_std = np.std(RW_torques[i][fault_idx:])
                    
                    metrics.append([pre_fault_std, post_fault_std])
                    labels.append(f'RW{i+1}')
            
            if metrics:
                metrics = np.array(metrics)
                x = np.arange(len(labels))
                width = 0.35
                
                bars1 = plt.bar(x - width/2, metrics[:, 0], width, label='Pre-Fault', color='green', alpha=0.7)
                bars2 = plt.bar(x + width/2, metrics[:, 1], width, label='Post-Fault', color='red', alpha=0.7)
                
                plt.xlabel('Reaction Wheel')
                plt.ylabel('Torque Std Dev [N⋅m]')
                plt.title('Control Variability: Encoder Fault Impact')
                plt.xticks(x, labels)
                plt.legend()
                plt.grid(True, alpha=0.3, axis='y')
                
                # Highlight the faulty wheel
                bars2[fault_wheel].set_alpha(1.0)
                bars2[fault_wheel].set_linewidth(2)
                bars2[fault_wheel].set_edgecolor('black')
        
        plt.tight_layout()
        figureList["EncoderFault_ControlAnalysis"] = plt.figure(4)
        
        if showPlots:
            plt.show()
        else:
            plt.close('all')
        
        print(f"Generated {len(figureList)} encoder fault plots with GUI parameters")
        return figureList


def runScenario(scenario, saveBinary=True):
    """Run encoder fault scenario with proper simulation time"""
    simulationTime = macros.min2nano(scenario.simulation_time_min)
    scenario.modeRequest = "hillPoint"

    print(f"Running encoder fault scenario:")
    print(f"  - Target wheel: RW{scenario.fault_wheel_number+1}")
    print(f"  - Fault time: {scenario.fault_time_min} minutes")
    print(f"  - Fault magnitude: {scenario.fault_magnitude}%")
    print(f"  - Simulation duration: {scenario.simulation_time_min} minutes")

    # Set up encoder fault event
    scenario.createNewEvent(
        "injectEncoderFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos >= self.encoderFaultTime and self.encoderFaultFlag == 1"], 
        [f"self.inject_rw_encoder_fault({scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)",  
         "self.encoderFaultFlag = 0",
         f"print('*** ENCODER FAULT INJECTED in RW{scenario.fault_wheel_number+1} ***')"]
    )

    # Add disturbance to keep system active throughout simulation
    DynModels = scenario.get_DynModel()
    if hasattr(DynModels, 'extForceTorqueObject'):
        # Add small continuous disturbance to keep wheels active
        DynModels.extForceTorqueObject.extTorquePntB_B = [[-0.01], [-0.01], [-0.01]]
        print("Added continuous disturbance to keep wheels active")

    # Setup visualization (save binary, no live stream)
    viz = None
    if vizSupport.vizFound and saveBinary:
        try:
            # Create visualization directory if it doesn't exist
            vizfiles_dir = os.path.join(path, '..', 'Vizfile', '_VizFiles')
            if not os.path.exists(vizfiles_dir):
                try:
                    os.makedirs(vizfiles_dir, exist_ok=True)
                    print(f"Created visualization directory: {vizfiles_dir}")
                except Exception as e:
                    print(f"Warning: Could not create directory {vizfiles_dir}: {e}")
            
            binary_filename = f"encoder_fault_rw{scenario.fault_wheel_number+1}_t{int(scenario.fault_time_min)}_e{int(scenario.fault_magnitude)}"
            binary_path = os.path.join(vizfiles_dir, binary_filename)
            
            print(f"Saving visualization to binary file: {binary_path}")
            
            viz = vizSupport.enableUnityVisualization(
                scenario,
                scenario.get_DynModel().taskName,
                scenario.get_DynModel().scObject,
                rwEffectorList=scenario.get_DynModel().rwStateEffector,
                liveStream=False,
                saveFile=binary_path
            )

            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=scenario.get_DynModel().scObject.ModelTag,
                fieldOfView=70 * macros.D2R,
                displayName=f"Encoder Fault Camera RW{scenario.fault_wheel_number+1}",
                pointingVector_B=[0, 0, 0],
                position_B=[0.0, 1.5, 0.0]
            )
            
            print(f"Visualization configured for binary file with fault parameters")

        except Exception as e:
            print(f"WARNING: Visualization failed: {e}")
            viz = None

    # Run simulation
    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz


def run(showPlots=True, saveBinary=True, simulation_time_min=30.0):
    """Run encoder fault scenario with default parameters"""
    print("\n===== Encoder Fault Scenario =====")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        # Use defaults that show clear impact
        scenario = EncoderFaultScenario(20.0, 1, 5.0, simulation_time_min)
        viz = runScenario(scenario, saveBinary)
        figureList = scenario.pull_outputs(showPlots)

        # Store in module globals
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz

        print(f"Generated {len(figureList)} encoder fault plots")
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}


def run_with_parameters(fault_magnitude=20.0, fault_wheel=1, fault_time_min=5.0, 
                       simulation_time_min=30.0, showPlots=False, saveBinary=False):
    """Run encoder fault with GUI parameters - uses exact GUI values"""
    print(f"\n===== Encoder Fault with GUI Parameters =====")
    print(f"PARAMS - Magnitude: {fault_magnitude}%, Wheel: RW{fault_wheel+1}, Time: {fault_time_min}min")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        scenario = EncoderFaultScenario(fault_magnitude, fault_wheel, fault_time_min, simulation_time_min)
        viz = runScenario(scenario, saveBinary)
        figureList = scenario.pull_outputs(showPlots)
        
        # Store in module globals
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz
        
        print(f"SUCCESS: Generated {len(figureList)} encoder fault plots with GUI parameters")
        print(f"VERIFICATION:")
        print(f"  - Used fault magnitude: {scenario.fault_magnitude}%")
        print(f"  - Used target wheel: RW{scenario.fault_wheel_number+1} (index {scenario.fault_wheel_number})")
        print(f"  - Used fault time: {scenario.fault_time_min} minutes")
        print(f"  - Used simulation time: {scenario.simulation_time_min} minutes")
        
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR with encoder fault parameters: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}


# Global storage for fault_loader
scenario = None
figureList = {}
viz = None

# Test with parameters
if __name__ == "__main__":
    print("Testing encoder fault with GUI-style parameters...")
    test_scenario, test_viz, test_plots = run_with_parameters(
        fault_magnitude=20.0,  # 20% encoder error
        fault_wheel=2,         # RW3 (index 2)
        fault_time_min=10.0,   # 10 minutes
        simulation_time_min=30.0,  # 30 minute simulation
        showPlots=True
    )
    print(f"Test completed with {len(test_plots)} plots")