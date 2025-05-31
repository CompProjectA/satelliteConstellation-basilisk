
#!/usr/bin/env python
"""
friction_fault.py

A Basilisk scenario that simulates spacecraft dynamics with friction faults
in the reaction wheels and properly saves binary files for Vizard visualization.

ENHANCED: Comprehensive plotting functionality for friction fault analysis.
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

class FrictionFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating friction faults in reaction wheels.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    ENHANCED: Comprehensive plotting functionality.
    """
    def __init__(self):
        super(FrictionFaultScenario, self).__init__()
        self.name = 'FrictionFaultScenario'
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
        self.fault_magnitude = 0.0005
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
        ENHANCED: Generate comprehensive friction fault plots
        """
        print("Generating friction fault analysis plots...")
        
        # Get attitude error data
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Get RW data
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        
        # Get friction data
        RW_friction = []
        for i in range(num_RW):
            try:
                friction_data = np.delete(self.rwLogs[i].u_f, 0, 0)
                RW_friction.append(friction_data)
            except:
                # Create placeholder friction data if not available
                RW_friction.append(np.zeros_like(timeData))
        
        # Calculate fault time in minutes
        fault_time_min = self.oneTimeFaultTime * macros.NANO2MIN
        
        # Print diagnostic information
        print("Initial RW Speeds:", RW_speeds[0] if len(RW_speeds) > 0 else "No data")
        print("Final RW Speeds:", RW_speeds[-1] if len(RW_speeds) > 0 else "No data")
        print("Initial Attitude Error:", np.linalg.norm(sigma_BR[0]) if len(sigma_BR) > 0 else "No data")
        print("Final Attitude Error:", np.linalg.norm(sigma_BR[-1]) if len(sigma_BR) > 0 else "No data")
        print(f"Friction fault time: {fault_time_min:.2f} minutes")
        print(f"Faulty wheel: {self.fault_wheel_number}")
        print(f"Fault magnitude: {self.fault_magnitude} N⋅m")

        # Clear previous plots
        plt.close('all')
        
        # Create comprehensive friction fault plots
        figureList = {}
        
        # Figure 1: Attitude Control Analysis
        fig1 = plt.figure(figsize=(12, 8))
        
        # Attitude error norm
        plt.subplot(2, 2, 1)
        attitude_error_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.plot(timeData, attitude_error_norm, 'b-', linewidth=2, label="Attitude Error Norm")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Norm")
        plt.title(f"Attitude Error (Friction Fault, RW {self.fault_wheel_number})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Individual attitude error components
        plt.subplot(2, 2, 2)
        for i in range(3):
            plt.plot(timeData, sigma_BR[:, i], label=f"σ_{i+1}")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Components")
        plt.title("Attitude Error Components")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Angular velocity error
        plt.subplot(2, 2, 3)
        omega_norm = np.linalg.norm(omega_BR_B, axis=1)
        plt.plot(timeData, omega_norm, 'g-', linewidth=2, label="Angular Velocity Error")
        plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Angular Velocity Error (rad/s)")
        plt.title("Angular Velocity Error")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Rate of attitude error change
        plt.subplot(2, 2, 4)
        if len(attitude_error_norm) > 1:
            error_rate = np.diff(attitude_error_norm) / np.diff(timeData)
            error_rate_time = timeData[1:]
            plt.plot(error_rate_time, error_rate, 'purple', linewidth=2, label="Error Rate")
            plt.axvline(fault_time_min, linestyle="--", color="red", linewidth=2, label="Friction Fault")
            plt.xlabel("Time (min)")
            plt.ylabel("Attitude Error Rate (1/min)")
            plt.title("Rate of Attitude Error Change")
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["FrictionFault_AttitudeAnalysis"] = fig1
        
        # Figure 2: Reaction Wheel Analysis
        fig2 = plt.figure(figsize=(12, 8))
        
        # All wheel speeds
        plt.subplot(2, 2, 1)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(num_RW):
            label_suffix = " (FAULTY)" if i == self.fault_wheel_number else ""
            plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=2, 
                    label=f"RW {i}{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title("Reaction Wheel Speeds")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Normalized wheel speeds for comparison
        plt.subplot(2, 2, 2)
        # Normalize by initial values to show relative changes
        initial_speeds = RW_speeds[0, :] if len(RW_speeds) > 0 else np.ones(num_RW)
        initial_speeds[initial_speeds == 0] = 1  # Avoid division by zero
        RW_speeds_norm = RW_speeds / initial_speeds[np.newaxis, :]
        
        for i in range(num_RW):
            label_suffix = " (FAULTY)" if i == self.fault_wheel_number else ""
            plt.plot(timeData, RW_speeds_norm[:, i], color=colors[i], linewidth=2, 
                    label=f"RW {i} Norm{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Normalized RW Speed")
        plt.title("Normalized Reaction Wheel Speeds")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Focus on faulty wheel
        plt.subplot(2, 2, 3)
        plt.plot(timeData, RW_speeds[:, self.fault_wheel_number], 'red', linewidth=3, 
                label=f"Faulty RW {self.fault_wheel_number}")
        # Show the effect of increased friction
        fault_indices = timeData >= fault_time_min
        if np.any(fault_indices):
            # Simulate reduced speed due to friction
            pre_fault_speed = RW_speeds[timeData < fault_time_min, self.fault_wheel_number]
            if len(pre_fault_speed) > 0:
                avg_pre_fault = np.mean(pre_fault_speed[-10:]) if len(pre_fault_speed) >= 10 else np.mean(pre_fault_speed)
                friction_effect = avg_pre_fault * (1 - self.fault_magnitude * 100)  # Simplified friction effect
                plt.axhline(friction_effect, linestyle=':', color='darkred', linewidth=2, 
                           label="Expected with Friction")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("RW Speed (rad/s)")
        plt.title(f"Faulty Wheel {self.fault_wheel_number} Analysis")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Power consumption estimate
        plt.subplot(2, 2, 4)
        # Estimate power based on speed and friction
        base_power = 0.5  # Base power consumption in Watts
        power_consumption = np.full_like(timeData, base_power)
        
        for i in range(num_RW):
            # Add power based on wheel speed
            speed_power = np.abs(RW_speeds[:, i]) * 0.001  # Simplified power model
            power_consumption += speed_power
            
            # Add extra power for faulty wheel due to friction
            if i == self.fault_wheel_number:
                fault_indices = timeData >= fault_time_min
                extra_friction_power = self.fault_magnitude * np.abs(RW_speeds[fault_indices, i]) * 200
                power_consumption[fault_indices] += extra_friction_power
        
        plt.plot(timeData, power_consumption, 'orange', linewidth=2, label="Total Power")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Power Consumption (W)")
        plt.title("Estimated Power Consumption")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["FrictionFault_RWAnalysis"] = fig2
        
        # Figure 3: Friction Analysis
        fig3 = plt.figure(figsize=(12, 8))
        
        # Friction torques for all wheels
        plt.subplot(2, 2, 1)
        for i in range(num_RW):
            label_suffix = " (FAULTY)" if i == self.fault_wheel_number else ""
            plt.plot(timeData, RW_friction[i], color=colors[i], linewidth=2, 
                    label=f"RW {i} Friction{label_suffix}")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Friction Torque (N⋅m)")
        plt.title("Reaction Wheel Friction Torques")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Faulty wheel friction focus
        plt.subplot(2, 2, 2)
        faulty_friction = RW_friction[self.fault_wheel_number]
        plt.plot(timeData, faulty_friction, 'red', linewidth=2, 
                label=f"RW {self.fault_wheel_number} Friction")
        # Show expected friction increase
        baseline_friction = np.mean(faulty_friction[timeData < fault_time_min]) if np.any(timeData < fault_time_min) else 0
        expected_friction = baseline_friction + self.fault_magnitude
        fault_indices = timeData >= fault_time_min
        if np.any(fault_indices):
            plt.axhline(expected_friction, linestyle=':', color='darkred', linewidth=2, 
                       label=f"Expected (+{self.fault_magnitude} N⋅m)")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Friction Torque (N⋅m)")
        plt.title(f"Faulty Wheel {self.fault_wheel_number} Friction")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Friction vs Speed relationship
        plt.subplot(2, 2, 3)
        faulty_speed = RW_speeds[:, self.fault_wheel_number]
        plt.scatter(faulty_speed[timeData < fault_time_min], faulty_friction[timeData < fault_time_min], 
                   c='blue', alpha=0.6, label="Pre-fault", s=20)
        if np.any(fault_indices):
            plt.scatter(faulty_speed[fault_indices], faulty_friction[fault_indices], 
                       c='red', alpha=0.6, label="Post-fault", s=20)
        plt.xlabel("RW Speed (rad/s)")
        plt.ylabel("Friction Torque (N⋅m)")
        plt.title(f"Friction vs Speed - RW {self.fault_wheel_number}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Total system friction
        plt.subplot(2, 2, 4)
        total_friction = np.sum(np.abs(RW_friction), axis=0)
        plt.plot(timeData, total_friction, 'purple', linewidth=2, label="Total System Friction")
        plt.axvline(fault_time_min, linestyle="--", color="black", linewidth=2, label="Friction Fault")
        plt.xlabel("Time (min)")
        plt.ylabel("Total Friction Torque (N⋅m)")
        plt.title("Total System Friction")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        figureList["FrictionFault_FrictionAnalysis"] = fig3

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList

def runScenario(scenario, saveBinary=True):
    """Run the friction fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        [f"self.get_DynModel().AddRWFault('friction',{scenario.fault_magnitude},{scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeRWFaultFlag=0"]
    )

    scenario.createNewEvent(
        "addRepeatedRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.repeatRWFaultFlag==1"],
        ["self.get_DynModel().PeriodicRWFault(360,'friction',0.1,1, self.TotalSim.CurrentNanos)", 
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
        binary_filename = "friction_fault_viz" if saveBinary else None
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
            displayName="Friction Fault Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the friction fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Friction Fault Scenario =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    
    scenario = FrictionFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'friction_fault_viz_UnityViz.bin'")
        print("You can now open this file in Vizard for visualization.")

    print(f"Generated {len(figureList)} friction fault analysis plots")
    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Friction Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)

# Import necessary standard libraries
import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import utilities from Basilisk
from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Determine file path and append relevant directories to system path
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')

# Import simulation modules
from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

# Define the scenario class for RW fault analysis
class scenario_AddRWFault(BSKSim, BSKScenario):
    def __init__(self):
        super(scenario_AddRWFault, self).__init__()
        self.name = 'scenario_AddRWFault'  # Scenario name for identification
        self.msgRecList = {}  # Placeholder for any future message subscriptions
        self.sNavTransName = "sNavTransMsg"  # Navigation translation message name
        self.attGuidName = "attGuidMsg"  # Attitude guidance message name
        self.cameraLocation = [0.0, 3.0, 0.0]  # Camera position for visualization

        # Set the dynamics and flight software models
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        # Set up initial conditions and log desired outputs
        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection configuration
        self.oneTimeRWFaultFlag = 1  # Enable one-time fault
        self.repeatRWFaultFlag = 1  # Enable repeated fault
        self.oneTimeFaultTime = macros.min2nano(10.)  # Inject one-time fault at 10 minutes
        self.get_DynModel().RWFaultLog = []  # Initialize fault log

    def configure_initial_conditions(self):
        # Set classical orbital elements
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0  # Semi-major axis in meters
        oe.e = 0.01  # Eccentricity
        oe.i = 33.3 * macros.D2R  # Inclination (converted to radians)
        oe.Omega = 48.2 * macros.D2R  # RAAN
        oe.omega = 347.8 * macros.D2R  # Argument of perigee
        oe.f = 85.3 * macros.D2R  # True anomaly

        # Get dynamics model and compute position/velocity from orbital elements
        DynModel = self.get_DynModel()
        mu = DynModel.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)

        # Initialize spacecraft position and velocity
        DynModel.scObject.hub.r_CN_NInit = rN
        DynModel.scObject.hub.v_CN_NInit = vN

        # Initialize spacecraft attitude and angular velocity
        DynModel.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModel.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        # Set up message recording for RW speeds and friction torques
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # Record RW speed message
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

        # Record RW output torque messages
        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def pull_outputs(self, showPlots):

        # FSW process outputs, remove first data point as it is before FSW is called
        attErrRec = self.msgRecList[self.attGuidName]

        # B refers to the body frame (attached to the spacecraft).
        # R refers to the reference frame (desired orientation).
        # So, sigma_BR gives the orientation difference between the body and reference frames.

        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)

        # Extract recorded RW data (speed and friction)
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = []
        for i in range(numRW):
            RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))

        # Estimate RW temperatures based on speed and friction
        self.no_cooling, self.with_cooling = self.calculate_temperatures(RW_speeds, RW_friction)

        # Plotting section
        BSK_plt.clear_all_plots()
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Plot RW speeds and friction
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.get_DynModel().RWFaultLog)
        BSK_plt.plot_attitude_error(timeData, sigma_BR)

        # Plot temperatures
        self.plot_rw_temperature(timeData, self.no_cooling, numRW)
        self.plot_rw_C_temperature(timeData, self.with_cooling, numRW)

        # Return or show/save figures
        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["RWSpeeds", "RWFriction","RWTemperatures(c)","RWTemperatures","attitudeErrorNorm" ]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        return figureList

    def calculate_temperatures(self, rw_speeds, rw_friction):
        """Estimate temperatures with and without cooling based on RW friction."""
        numRW = len(rw_friction)
        num_samples = len(rw_speeds)
        no_cooling = []
        with_cooling = []

        T_ambient = 10.0  # Ambient temperature in Celsius

        for rw in range(numRW):
            temp_nc = np.zeros(num_samples)
            temp_c = np.zeros(num_samples)
            temp_nc[0] = T_ambient
            temp_c[0] = T_ambient

            for i in range(1, num_samples):
                omega = rw_speeds[i, rw] * 2 * np.pi / 60
                P_friction = abs(rw_friction[rw][i] * omega) if i < len(rw_friction[rw]) else 0
                temp_rise = P_friction * 0.2
                cooling = 0.005 * (temp_c[i-1] - T_ambient)  # Arbitrary cooling term

                temp_nc[i] = temp_nc[i-1] + temp_rise
                temp_c[i] = temp_c[i-1] + (temp_rise) - cooling

            no_cooling.append(temp_nc)
            with_cooling.append(temp_c)

        return no_cooling, with_cooling




    def plot_rw_temperature(self, timeData, RW_temperatures, numRW):
        """Generate plot of RW temperatures over time."""
        plt.figure()
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('Reaction Wheel Temperatures(Without Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning')
        plt.tight_layout()

    def plot_rw_C_temperature(self, timeData, RW_temperatures, numRW):
        """Generate plot of RW temperatures over time."""
        plt.figure()
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('Reaction Wheel Temperatures(With Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning')
        plt.tight_layout()

    def runScenario(self):
        # Define simulation duration
        simulationTime = macros.min2nano(30.)
        self.modeRequest = "hillPoint"  # Request a specific control mode

        # If Unity viz is available, set up visualization
        if vizSupport.vizFound:
            viz = vizSupport.enableUnityVisualization(
                self,
                self.get_DynModel().taskName,
                self.get_DynModel().scObject,
                rwEffectorList=self.get_DynModel().rwStateEffector,
                liveStream=True,
                saveFile="friction_fault"
            )

            # Set up default camera
            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=self.get_DynModel().scObject.ModelTag,
                fieldOfView=30 * macros.D2R,
                displayName="RW Camera",
                pointingVector_B=[0, 0, 0],
                position_B=self.cameraLocation
            )

        # Run the simulation
        self.InitializeSimulation()
        self.ConfigureStopTime(simulationTime)
        self.ExecuteSimulation()

# Entry point for running the scenario as a script
def run(showPlots=True):
    scenario = scenario_AddRWFault()
    scenario.runScenario()  # Run the simulation
    scenario.pull_outputs(showPlots)  # Plot or save outputs

# Run the script if executed directly
if __name__ == "__main__":
    run(True)
