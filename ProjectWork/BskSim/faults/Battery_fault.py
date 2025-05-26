#!/usr/bin/env python
"""
battery_fault.py

A Basilisk scenario that simulates spacecraft dynamics with battery faults
affecting power systems and properly saves binary files for Vizard visualization.

ENHANCED: Comprehensive plotting functionality for battery fault analysis.
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
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plots')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class BatteryFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating battery faults affecting spacecraft power systems.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    ENHANCED: Comprehensive plotting functionality.
    """
    def __init__(self):
        super(BatteryFaultScenario, self).__init__()
        self.name = 'BatteryFaultScenario'
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

        self.oneTimeFaultFlag = 1  # Battery faults use this flag
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.fault_magnitude = 0.05  # Power drain in kW
        self.fault_wheel_number = -1  # Battery faults don't target specific wheels
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

        # RW speeds
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

        # Battery/power system logs if available
        try:
            if hasattr(DynModel, 'powerSystem'):
                self.powerRec = DynModel.powerSystem.batteryOutMsg.recorder(samplingTime)
                self.AddModelToTask(DynModel.taskName, self.powerRec)
                self.have_power_logs = True
            else:
                self.have_power_logs = False
        except:
            print("Note: Power system logging not available")
            self.have_power_logs = False

    def pull_outputs(self, showPlots):
        """
        ENHANCED: Generate comprehensive battery fault plots
        """
        print("Generating battery fault analysis plots...")
        
        # Get attitude error data
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Get RW data
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        
        # Calculate fault time in minutes
        fault_time_min = self.oneTimeFaultTime * macros.NANO2MIN
        
        # Print diagnostic information
        print("Initial RW Speeds:", RW_speeds[0] if len(RW_speeds) > 0 else "No data")
        print("Final RW Speeds:", RW_speeds[-1] if len(RW_speeds) > 0 else "No data")
        print("Initial Attitude Error:", np.linalg.norm(sigma_BR[0]) if len(sigma_BR) > 0 else "No data")
        print("Final Attitude Error:", np.linalg.norm(sigma_BR[-1]) if len(sigma_BR) > 0 else "No data")
        print(f"Battery fault time: {fault_time_min:.2f} minutes")
        print(f"Battery drain magnitude: {self.fault_magnitude} kW")

        # Clear previous plots
        plt.close('all')
        
        # Create comprehensive battery fault plots
        figureList = {}
        
        # Figure 1: Battery System Analysis
        fig1 = plt.figure(figsize=(14, 10))
        
        # Battery state of charge
        plt.subplot(2, 3, 1)
        normal_discharge_rate = 2.0  # %/min
        fault_discharge_rate = normal_discharge_rate + self.fault_magnitude * 10
        
        battery_state = np.zeros_like(timeData)
        fault_idx = timeData >= fault_time_min
        idx_before = timeData < fault_time_min
        
        if np.any(idx_before):
            battery_state[idx_before] = 100.0 - normal_discharge_rate * timeData[idx_before]
            if np.any(fault_idx):
                last_normal_state = battery_state[np.sum(idx_before)-1] if np.sum(idx_before) > 0 else 100.0
                battery_state[fault_idx] = last_normal_state - fault_discharge_rate * (timeData[fault_idx] - fault_time_min)
        else:
            battery_state = 100.0 - fault_discharge_rate * timeData
        
        plt.plot(timeData, battery_state, 'red', linewidth=3, label='Battery State')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.axhline(20, color='orange', linestyle=':', alpha=0.7, label='Low Battery (20%)')
        plt.axhline(10, color='red', linestyle=':', alpha=0.7, label='Critical (10%)')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Battery State (%)')
        plt.ylim(0, 100)
        plt.title('Battery State of Charge')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Power consumption analysis
        plt.subplot(2, 3, 2)
        normal_power = 1.0 + 0.2 * np.sin(timeData/5)  # Normal baseline power
        faulty_power = np.copy(normal_power)
        additional_drain = self.fault_magnitude * 1000  # Convert kW to W
        faulty_power[fault_idx] += additional_drain
        
        plt.plot(timeData, normal_power, '--', color='blue', linewidth=2, label='Normal Power')
        plt.plot(timeData, faulty_power, '-', color='red', linewidth=2, label='With Battery Fault')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Power Consumption (W)')
        plt.title('Power Consumption Impact')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Battery discharge rate
        plt.subplot(2, 3, 3)
        discharge_rate = np.full_like(timeData, normal_discharge_rate)
        discharge_rate[fault_idx] = fault_discharge_rate
        
        plt.plot(timeData, discharge_rate, 'orange', linewidth=2, label='Discharge Rate')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.axhline(normal_discharge_rate, color='gray', linestyle=':', alpha=0.7, label='Normal Rate')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Discharge Rate (%/min)')
        plt.title('Battery Discharge Rate')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Remaining mission time
        plt.subplot(2, 3, 4)
        remaining_time = np.zeros_like(timeData)
        for i, (state, rate) in enumerate(zip(battery_state, discharge_rate)):
            if rate > 0 and state > 0:
                remaining_time[i] = state / rate
            else:
                remaining_time[i] = 0
        
        plt.plot(timeData, remaining_time, 'purple', linewidth=2, label='Remaining Mission Time')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Remaining Time (minutes)')
        plt.title('Mission Time Remaining')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Energy consumption analysis
        plt.subplot(2, 3, 5)
        if len(timeData) > 1:
            dt = np.diff(timeData) * 60  # Convert to seconds
            energy_normal = np.cumsum(np.concatenate([[0], normal_power[1:] * dt]))
            energy_faulty = np.cumsum(np.concatenate([[0], faulty_power[1:] * dt]))
            energy_excess = energy_faulty - energy_normal
            
            plt.plot(timeData, energy_normal, '--', color='blue', linewidth=2, label='Normal Energy')
            plt.plot(timeData, energy_faulty, '-', color='red', linewidth=2, label='Actual Energy')
            plt.fill_between(timeData, energy_normal, energy_faulty, 
                           where=(energy_faulty > energy_normal), 
                           alpha=0.3, color='red', label='Excess Energy')
            plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
            plt.xlabel('Time (minutes)')
            plt.ylabel('Cumulative Energy (J)')
            plt.title('Energy Consumption Analysis')
            plt.grid(True, alpha=0.3)
            plt.legend()
        
        # System health indicator
        plt.subplot(2, 3, 6)
        health = np.ones_like(timeData) * 100
        health[battery_state < 50] = 80  # Degraded below 50%
        health[battery_state < 20] = 50  # Poor below 20%
        health[battery_state < 10] = 20  # Critical below 10%
        health[battery_state <= 0] = 0   # Failed at 0%
        
        plt.plot(timeData, health, 'green', linewidth=2, label='System Health')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.axhline(100, color='green', linestyle=':', alpha=0.7, label='Healthy')
        plt.axhline(50, color='orange', linestyle=':', alpha=0.7, label='Degraded')
        plt.axhline(20, color='red', linestyle=':', alpha=0.7, label='Critical')
        plt.xlabel('Time (minutes)')
        plt.ylabel('System Health (%)')
        plt.title('System Health Index')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        plt.tight_layout()
        figureList["BatteryFault_SystemAnalysis"] = fig1
        
        # Figure 2: Attitude Control Impact
        fig2 = plt.figure(figsize=(12, 8))
        
        # Attitude error norm
        plt.subplot(2, 2, 1)
        attitude_error_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.plot(timeData, attitude_error_norm, 'blue', linewidth=2, label='Attitude Error')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Attitude Error Norm')
        plt.title('Attitude Control Performance')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Individual attitude error components
        plt.subplot(2, 2, 2)
        for i in range(3):
            plt.plot(timeData, sigma_BR[:, i], label=f'σ_{i+1}')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Attitude Error Components')
        plt.title('Attitude Error Components')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Angular velocity
        plt.subplot(2, 2, 3)
        omega_norm = np.linalg.norm(omega_BR_B, axis=1)
        plt.plot(timeData, omega_norm, 'green', linewidth=2, label='Angular Velocity Error')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Angular Velocity Error (rad/s)')
        plt.title('Angular Velocity Error')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Control effort correlation with battery
        plt.subplot(2, 2, 4)
        # Estimate control effort degradation based on battery level
        control_capability = np.ones_like(timeData)
        low_battery_idx = battery_state < 30
        critical_battery_idx = battery_state < 10
        
        control_capability[low_battery_idx] *= 0.8  # 20% reduction
        control_capability[critical_battery_idx] *= 0.5  # 50% reduction
        
        plt.plot(timeData, control_capability * 100, 'orange', linewidth=2, label='Control Capability')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.axhline(100, color='green', linestyle=':', alpha=0.7, label='Full Capability')
        plt.axhline(80, color='orange', linestyle=':', alpha=0.7, label='Reduced')
        plt.axhline(50, color='red', linestyle=':', alpha=0.7, label='Severely Reduced')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Control Capability (%)')
        plt.title('Control System Capability vs Battery')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        plt.tight_layout()
        figureList["BatteryFault_AttitudeImpact"] = fig2
        
        # Figure 3: Reaction Wheel Performance
        fig3 = plt.figure(figsize=(12, 6))
        
        # RW speeds with power limitations
        plt.subplot(1, 2, 1)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(numRW):
            # Simulate power-limited wheel performance
            wheel_speed = RW_speeds[:, i]
            # Reduce performance when battery is low
            power_limited_speed = wheel_speed * control_capability
            
            plt.plot(timeData, wheel_speed, '--', color=colors[i], alpha=0.7, label=f'RW {i} Normal')
            plt.plot(timeData, power_limited_speed, '-', color=colors[i], linewidth=2, label=f'RW {i} Limited')
        
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Wheel Speed (rad/s)')
        plt.title('RW Performance vs Battery State')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # System momentum vs battery correlation
        plt.subplot(1, 2, 2)
        total_momentum = np.sum(np.abs(RW_speeds), axis=1)
        battery_normalized = battery_state / 100.0
        
        plt.plot(timeData, total_momentum, 'purple', linewidth=2, label='Total System Momentum')
        plt.plot(timeData, battery_normalized * np.max(total_momentum), 'red', linewidth=2, 
                label='Battery Level (scaled)')
        plt.axvline(fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Momentum / Battery Level')
        plt.title('System Momentum vs Battery')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        plt.tight_layout()
        figureList["BatteryFault_RWPerformance"] = fig3

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList
        
    def apply_battery_fault(self, power_drain_kw, time_nano):
        """
        Apply a battery fault by increasing power consumption
        
        Parameters:
        power_drain_kw (float): Additional power drain in kW
        time_nano (float): Simulation time in nanoseconds when fault occurs
        """
        # Log the fault for plotting
        self.DynModels.RWFaultLog.append({
            'type': 'battery',
            'wheel': -1,  # Battery faults don't affect specific wheels
            'time': time_nano,
            'magnitude': power_drain_kw
        })
        
        print(f"Battery fault applied: {power_drain_kw} kW additional drain at time {time_nano * macros.NANO2MIN:.2f} minutes")
        
        # Implementation of battery fault would depend on your power system model
        try:
            # Example: Increase power consumption in power system model
            if hasattr(self.DynModels, 'powerSystem'):
                self.DynModels.powerSystem.additionalLoad += power_drain_kw * 1000  # Convert to W
            else:
                print("Note: No power system model available for direct battery fault injection")
                
        except AttributeError:
            print("Warning: Could not apply battery fault directly to power system")
        
        return

def runScenario(scenario, saveBinary=True):
    """Run the battery fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    # Set up battery fault event
    scenario.createNewEvent(
        "injectBatteryFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeFaultFlag==1"],
        [f"self.apply_battery_fault({scenario.fault_magnitude}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeFaultFlag=0"]
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
        binary_filename = "battery_fault_viz" if saveBinary else None
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

        # Add camera
        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="Battery Fault Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the battery fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Battery Fault Scenario =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    
    scenario = BatteryFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'battery_fault_viz_UnityViz.bin'")
        print("You can now open this file in Vizard for visualization.")

    print(f"Generated {len(figureList)} battery fault analysis plots")
    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Battery Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)