#!/usr/bin/env python
"""
powerlimit_fault.py

A Basilisk scenario that simulates spacecraft dynamics with power limit faults
in the reaction wheels and properly saves binary files for Vizard visualization.

ENHANCED: Comprehensive plotting functionality for power limit fault analysis.
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

class PowerLimitFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating power limit faults in reaction wheels.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    ENHANCED: Comprehensive plotting functionality.
    """
    def __init__(self):
        super(PowerLimitFaultScenario, self).__init__()
        self.name = 'PowerLimitFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.cameraLocation = [0.0, 2.0, 0.0]

        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "red"},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "blue"},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "green"},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "yellow"}
        ]

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        self.oneTimeRWFaultFlag = 1
        self.repeatRWFaultFlag = 1
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.fault_magnitude = 0.5  # Default power limit in Watts
        self.fault_wheel_number = 3
        self.DynModels = self.get_DynModel()
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

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

    def pull_outputs(self, showPlots):
        """
        ENHANCED: Generate comprehensive power limit fault plots
        """
        print("Generating power limit fault analysis plots...")
        
        # Get attitude error data
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Get RW data
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        
        # Get command data if available
        RW_commands = []
        for i in range(num_RW):
            try:
                cmd_data = np.delete(self.rwLogs[i].u_cmd, 0, 0)
                RW_commands.append(cmd_data)
            except:
                # Create placeholder command data if not available
                RW_commands.append(np.zeros_like(timeData))
        
        # Calculate fault time in minutes
        fault_time_min = self.oneTimeFaultTime * macros.NANO2MIN
        
        # Print diagnostic information
        print("Initial RW Speeds:", RW_speeds[0] if len(RW_speeds) > 0 else "No data")
        print("Final RW Speeds:", RW_speeds[-1] if len(RW_speeds) > 0 else "No data")
        print("Initial Attitude Error:", np.linalg.norm(sigma_BR[0]) if len(sigma_BR) > 0 else "No data")
        print("Final Attitude Error:", np.linalg.norm(sigma_BR[-1]) if len(sigma_BR) > 0 else "No data")
        print(f"Power limit fault time: {fault_time_min:.2f} minutes")
        print(f"Faulty wheel: {self.fault_wheel_number}")
        print(f"Power limit: {self.fault_magnitude} W")

        # Clear previous plots
        plt.close('all')
        
        # Create comprehensive power limit fault plots
        figureList = {}
        
        # Figure 1: Attitude Control Analysis
        fig1 = plt.figure(figsize=(12, 8))
        
        # Attitude error norm
        plt.subplot(2, 2, 1)
        attitude_error_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.plot(timeData, attitude_error_norm, 'b-', linewidth=2, label="Attitude Error Norm")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Norm")
        plt.title(f"Attitude Error (Power Limit: {self.fault_magnitude}W, RW {self.fault_wheel_number})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Individual attitude error components
        plt.subplot(2, 2, 2)
        for i in range(3):
            plt.plot(timeData, sigma_BR[:, i], label=f"σ_{i+1}")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Components")
        plt.title("Attitude Error Components")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Angular velocity error
        plt.subplot(2, 2, 3)
        omega_norm = np.linalg.norm(omega_BR_B, axis=1)
        plt.plot(timeData, omega_norm, 'g-', linewidth=2, label="Angular Velocity Error")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Angular Velocity Error (rad/s)")
        plt.title("Angular Velocity Error")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Control performance degradation
        plt.subplot(2, 2, 4)
        # Calculate control performance as inverse of attitude error
        control_performance = 1.0 / (1.0 + attitude_error_norm)
        plt.plot(timeData, control_performance * 100, 'purple', linewidth=2, label="Control Performance")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Control Performance (%)")
        plt.title("Attitude Control Performance")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["PowerLimitFault_AttitudeAnalysis"] = fig1
        
        # Figure 2: Power Analysis
        fig2 = plt.figure(figsize=(12, 8))
        
        # Power request vs available power
        plt.subplot(2, 2, 1)
        # Simulate power request based on command torque
        power_request = np.abs(RW_commands[self.fault_wheel_number]) * 2.0 + 0.5  # Simplified power model
        power_available = np.full_like(timeData, 2.0)  # Normal power availability
        fault_indices = timeData >= fault_time_min
        power_available[fault_indices] = self.fault_magnitude  # Limited power after fault
        
        plt.plot(timeData, power_request, '--', color='blue', linewidth=2, label="Power Request")
        plt.plot(timeData, power_available, '-', color='red', linewidth=2, label="Available Power")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Power (W)")
        plt.title(f"Power Analysis - RW {self.fault_wheel_number}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Power performance ratio
        plt.subplot(2, 2, 2)
        power_performance = np.minimum(power_request, power_available) / power_request * 100
        power_performance[power_request == 0] = 100  # Avoid division by zero
        plt.plot(timeData, power_performance, 'orange', linewidth=2, label="Power Performance")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.axhline(100, color='green', linestyle=':', alpha=0.7, label="100% Performance")
        plt.xlabel("Time (min)")
        plt.ylabel("Power Performance (%)")
        plt.title("Power Performance Ratio")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Power deficit
        plt.subplot(2, 2, 3)
        power_deficit = np.maximum(0, power_request - power_available)
        plt.plot(timeData, power_deficit, 'red', linewidth=2, label="Power Deficit")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.fill_between(timeData, power_deficit, alpha=0.3, color='red')
        plt.xlabel("Time (min)")
        plt.ylabel("Power Deficit (W)")
        plt.title("Power Deficit Analysis")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Cumulative energy loss
        plt.subplot(2, 2, 4)
        if len(timeData) > 1:
            dt = np.diff(timeData) * 60  # Convert to seconds
            energy_loss = np.cumsum(np.concatenate([[0], power_deficit[1:] * dt]))  # Energy in Joules
            plt.plot(timeData, energy_loss, 'darkred', linewidth=2, label="Cumulative Energy Loss")
            plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
            plt.xlabel("Time (min)")
            plt.ylabel("Energy Loss (J)")
            plt.title("Cumulative Energy Loss")
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["PowerLimitFault_PowerAnalysis"] = fig2
        
        # Figure 3: Reaction Wheel Performance
        fig3 = plt.figure(figsize=(12, 8))
        
        # All wheel speeds
        plt.subplot(2, 2, 1)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(num_RW):
            label_suffix = " (POWER LIMITED)" if i == self.fault_wheel_number else ""
            plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=2, 
                    label=f"RW {i}{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title("Reaction Wheel Speeds")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Command vs actual for faulty wheel
        plt.subplot(2, 2, 2)
        faulty_command = RW_commands[self.fault_wheel_number]
        faulty_speed = RW_speeds[:, self.fault_wheel_number]
        
        # Simulate limited response due to power constraint
        limited_command = np.copy(faulty_command)
        fault_indices = timeData >= fault_time_min
        if np.any(fault_indices):
            # Limit command based on available power
            max_limited_cmd = self.fault_magnitude / 2.0  # Simplified relationship
            limited_command[fault_indices] = np.minimum(np.abs(faulty_command[fault_indices]), 
                                                       max_limited_cmd) * np.sign(faulty_command[fault_indices])
        
        plt.plot(timeData, faulty_command, '--', color='blue', linewidth=2, label="Requested Command")
        plt.plot(timeData, limited_command, '-', color='red', linewidth=2, label="Limited Command")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Command Torque (N⋅m)")
        plt.title(f"Command Limiting - RW {self.fault_wheel_number}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Wheel performance comparison
        plt.subplot(2, 2, 3)
        # Calculate relative performance of each wheel
        for i in range(num_RW):
            if np.std(RW_speeds[:, i]) > 0:  # Only plot if wheel is active
                normalized_speed = RW_speeds[:, i] / np.max(np.abs(RW_speeds[:, i]))
                label_suffix = " (LIMITED)" if i == self.fault_wheel_number else ""
                plt.plot(timeData, normalized_speed, color=colors[i], linewidth=2, 
                        label=f"RW {i} Norm{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Normalized Performance")
        plt.title("Wheel Performance Comparison")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # System-wide impact
        plt.subplot(2, 2, 4)
        # Calculate total system momentum
        total_momentum = np.sum(np.abs(RW_speeds), axis=1)
        plt.plot(timeData, total_momentum, 'purple', linewidth=2, label="Total System Momentum")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Power Limit Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Total Momentum (rad/s)")
        plt.title("System-Wide Impact")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["PowerLimitFault_RWPerformance"] = fig3

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList

def runScenario(scenario, saveBinary=True):
    """Run the power limit fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        [f"self.get_DynModel().AddRWFault('powerLimit',{scenario.fault_magnitude},{scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeRWFaultFlag=0"]
    )

    scenario.createNewEvent(
        "addRepeatedRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.repeatRWFaultFlag==1"],
        ["self.get_DynModel().PeriodicRWFault(360,'powerLimit',0.2,1, self.TotalSim.CurrentNanos)", 
         "self.setEventActivity('addRepeatedRWFault',True)"]
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
        binary_filename = "powerlimit_fault_viz" if saveBinary else None
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

        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="Power Limit Fault Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the power limit fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Power Limit Fault Scenario =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    
    scenario = PowerLimitFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'powerlimit_fault_viz_UnityViz.bin'")
        print("You can now open this file in Vizard for visualization.")

    print(f"Generated {len(figureList)} power limit fault analysis plots")
    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Power Limit Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)