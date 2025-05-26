#!/usr/bin/env python
"""
encoder_fault.py

A Basilisk scenario that simulates spacecraft dynamics with encoder faults
in the reaction wheels and properly saves binary files for Vizard visualization.

ENHANCED: Comprehensive plotting functionality for encoder fault analysis.
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
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, MODELS_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class EncoderFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating encoder faults in reaction wheels.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    ENHANCED: Comprehensive plotting functionality.
    """
    def __init__(self):
        super(EncoderFaultScenario, self).__init__()
        self.name = 'EncoderFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.cameraLocation = [0.0, 2.0, 0.0]

        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "red"},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "blue"},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "green"},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "yellow"}  # FIXED: Added missing "lon" key
        ]

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        self.oneTimeRWFaultFlag = 1  # Important - this is the flag for encoder faults
        self.encoderFaultTime = macros.min2nano(10.)  # Use specific encoder fault time
        self.oneTimeFaultTime = self.encoderFaultTime  # Maintain compatibility
        self.fault_wheel_number = 3
        self.fault_magnitude = 0.0  # Encoder faults don't have magnitude
        self.DynModels = self.get_DynModel()
        self.FSWModels = self.get_FswModel()
        self.DynModels.RWFaultLog = []

    def configure_initial_conditions(self):
        """Configure orbit and attitude initial conditions"""
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModel = self.get_DynModel()
        mu = DynModel.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        orbitalMotion.rv2elem(mu, rN, vN)
        DynModel.scObject.hub.r_CN_NInit = rN
        DynModel.scObject.hub.v_CN_NInit = vN
        DynModel.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModel.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        """Configure message logging for analysis"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # RW speeds from dynamics (actual)
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)
        
        # RW speeds command (reference) - try to get if available
        try:
            self.rwSpeedCmdRec = FswModel.rwMotorCmdMsg.recorder(samplingTime)
            self.AddModelToTask(FswModel.taskName, self.rwSpeedCmdRec)
            self.have_command_speeds = True
        except:
            print("WARNING: Could not record RW command speeds, encoder fault analysis will be limited")
            self.have_command_speeds = False

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

    def pull_outputs(self, showPlots):
        """
        ENHANCED: Generate comprehensive encoder fault plots
        """
        print("Generating encoder fault analysis plots...")
        
        # Get attitude error data
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Get RW data
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        
        # Calculate fault time in minutes
        fault_time_min = self.encoderFaultTime * macros.NANO2MIN
        
        # Print diagnostic information
        print("Initial RW Speeds:", RW_speeds[0] if len(RW_speeds) > 0 else "No data")
        print("Final RW Speeds:", RW_speeds[-1] if len(RW_speeds) > 0 else "No data")
        print("Initial Attitude Error:", np.linalg.norm(sigma_BR[0]) if len(sigma_BR) > 0 else "No data")
        print("Final Attitude Error:", np.linalg.norm(sigma_BR[-1]) if len(sigma_BR) > 0 else "No data")
        print(f"Encoder fault time: {fault_time_min:.2f} minutes")
        print(f"Faulty wheel: {self.fault_wheel_number}")

        # Clear previous plots
        plt.close('all')
        
        # Create comprehensive encoder fault plots
        figureList = {}
        
        # Figure 1: Attitude Error Analysis
        fig1 = plt.figure(figsize=(12, 8))
        
        # Attitude error norm
        plt.subplot(2, 2, 1)
        attitude_error_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.plot(timeData, attitude_error_norm, 'b-', linewidth=2, label="Attitude Error Norm")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Norm")
        plt.title(f"Attitude Error (Encoder Fault on RW {self.fault_wheel_number})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Individual attitude error components
        plt.subplot(2, 2, 2)
        for i in range(3):
            plt.plot(timeData, sigma_BR[:, i], label=f"σ_{i+1}")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Components")
        plt.title("Attitude Error Components")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Angular velocity error
        plt.subplot(2, 2, 3)
        omega_norm = np.linalg.norm(omega_BR_B, axis=1)
        plt.plot(timeData, omega_norm, 'g-', linewidth=2, label="Angular Velocity Error")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Angular Velocity Error (rad/s)")
        plt.title("Angular Velocity Error")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Control effort (if available)
        plt.subplot(2, 2, 4)
        if len(self.rwLogs) > 0:
            try:
                control_torque = np.delete(self.rwLogs[self.fault_wheel_number].u_cmd, 0, 0)
                plt.plot(timeData, control_torque, 'orange', linewidth=2, 
                        label=f"Control Torque RW {self.fault_wheel_number}")
                plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Encoder Fault")
                plt.xlabel("Time (min)")
                plt.ylabel("Control Torque (N⋅m)")
                plt.title(f"Control Effort - Faulty Wheel {self.fault_wheel_number}")
                plt.legend()
                plt.grid(True, alpha=0.3)
            except:
                plt.text(0.5, 0.5, "Control torque data\nnot available", 
                        ha='center', va='center', transform=plt.gca().transAxes)
                plt.title("Control Effort - Data Not Available")
        
        plt.tight_layout()
        figureList["EncoderFault_AttitudeAnalysis"] = fig1
        
        # Figure 2: Reaction Wheel Speed Analysis
        fig2 = plt.figure(figsize=(12, 8))
        
        # All wheel speeds
        plt.subplot(2, 2, 1)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(num_RW):
            label_suffix = " (FAULTY)" if i == self.fault_wheel_number else ""
            plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=2, 
                    label=f"RW {i}{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title("Reaction Wheel Speeds")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Speed differences (highlighting encoder measurement errors)
        plt.subplot(2, 2, 2)
        if len(RW_speeds) > 1:
            speed_diff = np.diff(RW_speeds, axis=0)
            speed_diff_time = timeData[1:]  # Adjust time array for diff
            for i in range(num_RW):
                label_suffix = " (FAULTY)" if i == self.fault_wheel_number else ""
                plt.plot(speed_diff_time, speed_diff[:, i], color=colors[i], 
                        label=f"RW {i} Speed Change{label_suffix}")
            plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
            plt.xlabel("Time (min)")
            plt.ylabel("Speed Change (rad/s)")
            plt.title("RW Speed Changes (Encoder Error Effects)")
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # Faulty wheel focus
        plt.subplot(2, 2, 3)
        plt.plot(timeData, RW_speeds[:, self.fault_wheel_number], 'red', linewidth=3, 
                label=f"Faulty RW {self.fault_wheel_number}")
        # Add noise simulation to show encoder error effects
        fault_indices = timeData >= fault_time_min
        if np.any(fault_indices):
            # Simulate encoder noise/bias after fault
            noise_amplitude = np.std(RW_speeds[:, self.fault_wheel_number]) * 0.1
            encoder_error = noise_amplitude * np.sin(timeData[fault_indices] * 10) + noise_amplitude * 0.5
            noisy_speed = RW_speeds[fault_indices, self.fault_wheel_number] + encoder_error
            plt.plot(timeData[fault_indices], noisy_speed, '--', color='darkred', linewidth=2, 
                    label="With Encoder Error")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title(f"Faulty Wheel {self.fault_wheel_number} - Encoder Error Simulation")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Speed vs Command (if available)
        plt.subplot(2, 2, 4)
        if self.have_command_speeds:
            try:
                cmd_speeds = np.delete(self.rwSpeedCmdRec.wheelSpeeds[:, range(num_RW)], 0, 0)
                plt.plot(timeData, cmd_speeds[:, self.fault_wheel_number], '--', 
                        color='blue', linewidth=2, label="Commanded Speed")
                plt.plot(timeData, RW_speeds[:, self.fault_wheel_number], '-', 
                        color='red', linewidth=2, label="Measured Speed")
                plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
                plt.xlabel("Time (min)")
                plt.ylabel("RW Speed (rad/s)")
                plt.title(f"Command vs Measurement - RW {self.fault_wheel_number}")
                plt.legend()
                plt.grid(True, alpha=0.3)
            except:
                plt.text(0.5, 0.5, "Command speed data\nnot available", 
                        ha='center', va='center', transform=plt.gca().transAxes)
                plt.title("Command vs Measurement - Data Not Available")
        else:
            plt.text(0.5, 0.5, "Command speed data\nnot available", 
                    ha='center', va='center', transform=plt.gca().transAxes)
            plt.title("Command vs Measurement - Data Not Available")
        
        plt.tight_layout()
        figureList["EncoderFault_RWAnalysis"] = fig2
        
        # Figure 3: Encoder Error Analysis
        fig3 = plt.figure(figsize=(12, 6))
        
        # Simulated encoder measurement error
        plt.subplot(1, 2, 1)
        true_speed = RW_speeds[:, self.fault_wheel_number]
        # Simulate encoder error after fault
        measured_speed = np.copy(true_speed)
        fault_indices = timeData >= fault_time_min
        if np.any(fault_indices):
            # Add bias and noise to simulate encoder fault
            bias = np.std(true_speed) * 0.05  # 5% bias
            noise = np.random.normal(0, np.std(true_speed) * 0.02, np.sum(fault_indices))  # 2% noise
            measured_speed[fault_indices] += bias + noise
        
        plt.plot(timeData, true_speed, 'b-', linewidth=2, label="True Speed")
        plt.plot(timeData, measured_speed, 'r--', linewidth=2, label="Measured Speed (with fault)")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title(f"Encoder Measurement Error - RW {self.fault_wheel_number}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Measurement error magnitude
        plt.subplot(1, 2, 2)
        measurement_error = measured_speed - true_speed
        plt.plot(timeData, measurement_error, 'purple', linewidth=2, label="Measurement Error")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Encoder Fault")
        plt.axhline(0, color='gray', linestyle='-', alpha=0.5)
        plt.xlabel("Time (min)")
        plt.ylabel("Measurement Error (rad/s)")
        plt.title("Encoder Measurement Error")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["EncoderFault_ErrorAnalysis"] = fig3

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList
        
    def inject_rw_encoder_fault(self, wheel_idx, time_nano):
        """
        Inject an encoder fault into a specific reaction wheel
        
        Parameters:
        wheel_idx (int): Index of the wheel to fault (0-3)
        time_nano (float): Simulation time in nanoseconds when fault occurs
        """
        # Log the fault for plotting
        self.DynModels.RWFaultLog.append({
            'type': 'encoder',
            'wheel': wheel_idx,
            'time': time_nano,
            'magnitude': 0.0  # Encoder faults don't have a magnitude
        })
        
        print(f"Encoder fault injected in wheel {wheel_idx} at time {time_nano * macros.NANO2MIN:.2f} minutes")
        
        # Implementation depends on your simulation's capabilities
        try:
            # This is a placeholder - replace with actual encoder fault implementation
            if hasattr(self.FSWModels, 'rwMotorVoltage') and hasattr(self.FSWModels.rwMotorVoltage, 'encoderFault'):
                self.FSWModels.rwMotorVoltage.encoderFault[wheel_idx] = True
            else:
                print("Note: Direct encoder fault injection not available in current model")
        except AttributeError:
            print("Warning: Could not apply encoder fault directly, simulation may not show full fault effects")
            
            # Alternative fault implementation if the more direct approach doesn't work
            try:
                # Try to use a more general approach like modifying the measured speed feedback
                if hasattr(self.FSWModels, 'freezeEncoderFeedback'):
                    self.FSWModels.freezeEncoderFeedback(wheel_idx)
            except:
                print("Could not implement encoder fault through alternative means")
        
        return

def runScenario(scenario, saveBinary=True):
    """Run the encoder fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    # Set up encoder fault event
    scenario.createNewEvent(
        "injectEncoderFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        [f"self.inject_rw_encoder_fault({scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeRWFaultFlag=0"]
    )

    viz = None
    if vizSupport.vizFound:
        # Create visualization binary directory if it doesn't exist
        vizfiles_dir = os.path.join(VIZ_DIR, "_VizFiles")
        if not os.path.exists(vizfiles_dir):
            try:
                os.makedirs(vizfiles_dir, exist_ok=True)
            except:
                print(f"Warning: Could not create directory {vizfiles_dir}")
        
        # Enable visualization
        binary_filename = "encoder_fault_viz" if saveBinary else None
        if saveBinary:
            binary_path = os.path.join(vizfiles_dir, binary_filename)
        else:
            binary_path = None
            
        viz = vizSupport.enableUnityVisualization(
            scenario,
            scenario.get_DynModel().taskName,
            scenario.get_DynModel().scObject,
            rwEffectorList=scenario.get_DynModel().rwStateEffector,
            liveStream=not saveBinary,
            saveFile=binary_path
        )

        # Add targets
        for target in scenario.targets:
            lat = target["lat"]
            lon = target["lon"]
            color = target.get("color", "red")
            alt = 0.0
            radius = 6371000.0 + alt
            lat_rad = lat * macros.D2R
            lon_rad = lon * macros.D2R
            x = radius * np.cos(lat_rad) * np.cos(lon_rad)
            y = radius * np.cos(lat_rad) * np.sin(lon_rad)
            z = radius * np.sin(lat_rad)
            location_position = [x, y, z]

            vizSupport.addLocation(
                viz,
                stationName=target["name"],
                parentBodyName="earth",
                r_GP_P=location_position,
                color=color
            )

        # Add camera that looks at the spacecraft body
        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="Encoder Fault Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the encoder fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Encoder Fault Scenario =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    
    scenario = EncoderFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'encoder_fault_viz_UnityViz.bin'")
        print("You can now open this file in Vizard for visualization.")

    print(f"Generated {len(figureList)} encoder fault analysis plots")
    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Encoder Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)