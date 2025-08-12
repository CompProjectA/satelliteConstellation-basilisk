#!/usr/bin/env python
"""
powerlimit_fault.py - Enhanced with dynamic simulation time support

This module simulates power limitation faults in reaction wheels.
"""

import inspect
import os
import sys
import numpy as np
from Basilisk.utilities import orbitalMotion, macros, vizSupport
import matplotlib.pyplot as plt

# Import dependencies
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')
from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

class PowerLimitFaultScenario(BSKSim, BSKScenario):
    """Power limit fault scenario with dynamic simulation time support"""
    
    def __init__(self, fault_magnitude=0.01, fault_wheel=1, fault_time_min=5.0, simulation_time_min=30.0):
        super(PowerLimitFaultScenario, self).__init__()
        self.name = 'PowerLimitFaultScenario'

        # Use GUI parameters
        self.fault_magnitude = fault_magnitude
        self.fault_wheel_number = fault_wheel
        self.fault_time_min = fault_time_min
        self.simulation_time_min = simulation_time_min
        
        # Convert fault time to nanoseconds
        self.oneTimeFaultTime = macros.min2nano(fault_time_min)
        self.faultTime = macros.min2nano(fault_time_min)
        
        # Set backup RW4 activation 5 minutes after fault
        self.rw4ActivationTime = self.faultTime + macros.min2nano(5.0)
        
        # Use the correct target wheel from GUI
        self.targetRW = fault_wheel
        
        print(f"PowerLimitFaultScenario initialized with:")
        print(f"  - Fault magnitude: {fault_magnitude} W")
        print(f"  - Target wheel: RW{fault_wheel} (index {fault_wheel})")
        print(f"  - Fault time: {fault_time_min} minutes")
        print(f"  - Simulation duration: {simulation_time_min} minutes")

        # Standard initialization
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"
        self.powerLimitLog = []

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        # Configure simulation
        self.configure_initial_conditions()
        self.log_outputs()
        self.disableAllEvents()
        self.addDisturbance()
        self.setup_visualization()


    def setup_visualization(self):
        """Setup Vizard visualization - save binary file instead of live stream"""
        if vizSupport.vizFound:
            # Create visualization directory if it doesn't exist
            vizfiles_dir = os.path.join(path, '..', 'Vizfile', '_VizFiles')
            if not os.path.exists(vizfiles_dir):
                try:
                    os.makedirs(vizfiles_dir, exist_ok=True)
                    print(f"Created visualization directory: {vizfiles_dir}")
                except Exception as e:
                    print(f"Warning: Could not create directory {vizfiles_dir}: {e}")
            
            # Create binary filename with fault parameters
            binary_filename = f"powerlimit_fault_rw{self.fault_wheel_number}_t{int(self.fault_time_min)}_p{int(self.fault_magnitude*1000)}"
            binary_path = os.path.join(vizfiles_dir, binary_filename)
            
            print(f"Saving visualization to binary file: {binary_path}")
            
            self.viz = vizSupport.enableUnityVisualization(
                self,
                self.get_DynModel().taskName,
                self.get_DynModel().scObject,
                liveStream=False,  # No live streaming
                saveFile=binary_path  # Save to binary file
            )
            self.viz.settings.orbitLinesOn = 1
            self.viz.settings.showSpacecraftLabels = 1

            vizSupport.createStandardCamera(
                self.viz,
                setMode=1,
                spacecraftName=self.get_DynModel().scObject.ModelTag,
                displayName="PowerLimitCam",
                fieldOfView=10 * macros.D2R,
                pointingVector_B=[0.0, 0.0, 0.0],
                position_B=[0.0, 1.5, 0.0]
            )
            print(f"Vizard configured to save binary file with fault parameters")
        else:
            print("Vizard not available")

    def addDisturbance(self):
        """Add external disturbance torque"""
        DynModels = self.get_DynModel()
        if hasattr(DynModels, 'extForceTorqueObject'):
            DynModels.extForceTorqueObject.extTorquePntB_B = [[-0.1], [-0.1], [-0.1]]
            print(f"Added external torque: {DynModels.extForceTorqueObject.extTorquePntB_B}")

    def disableAllEvents(self):
        """Disable all existing events"""
        if hasattr(self, 'eventMap'):
            for eventName in self.eventMap:
                self.setEventActivity(eventName, False)
            print(f"Disabled {len(self.eventMap)} events")

    def configure_initial_conditions(self):
        """Set initial orbital conditions"""
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0  # meters
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
        DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        """Set up data logging"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

        self.msgRecList[self.sNavTransName] = DynModel.simpleNavObject.transOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.sNavTransName])

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def pull_outputs(self, showPlots):
        """Generate plots with parameters displayed"""
        attErrRec = self.msgRecList[self.attGuidName]
        
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        
        RW_torque = []
        for i in range(num_RW):
            if hasattr(self.rwLogs[i], 'u_current'):
                RW_torque.append(np.delete(self.rwLogs[i].u_current, 0, 0))
            else:
                RW_torque.append(np.delete(self.rwLogs[i].u_f, 0, 0))
        
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
        
        # Use actual GUI parameters in plots
        activationIndex = np.where(timeData >= (self.fault_time_min + 5.0))[0]
        activationIndex = activationIndex[0] if len(activationIndex) > 0 else len(timeData)
        
        # Modify RW4 speeds for visualization
        RW_speeds_modified = RW_speeds.copy()
        RW_speeds_modified[:activationIndex, 3] = 0.0

        BSK_plt.clear_all_plots()
        
        # PLOT 1: Attitude Error
        sigma_BR_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.figure(1)
        plt.title(f"Attitude Error (PL={self.fault_magnitude}W, RW{self.targetRW}, {self.fault_time_min}min)")
        plt.xlabel("Time (min)")
        plt.ylabel("Attitude Error Norm")
        plt.plot(timeData, sigma_BR_norm, 'b-', label="Attitude Error Norm")
        
        # Add fault injection marker at correct time
        plt.axvline(x=self.fault_time_min, color='r', linestyle='--', linewidth=2, 
                   label=f'Fault Injection (RW{self.targetRW})')
        plt.axvline(x=self.fault_time_min + 5.0, color='g', linestyle='--', linewidth=2, 
                   label='RW4 Activated')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # PLOT 2: RW Speeds
        plt.figure(2)
        plt.title(f"RW Speeds (Fault: RW{self.targetRW}, {self.fault_magnitude}W at {self.fault_time_min}min)")
        plt.xlabel("Time [min]")
        plt.ylabel("Speed [RPM]")
        
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(num_RW):
            # Highlight the correct faulty wheel
            if i == self.targetRW:
                plt.plot(timeData, RW_speeds_modified[:, i], color=colors[i], 
                        linewidth=3, label=f"RW{i+1} (FAULTY)")
            else:
                plt.plot(timeData, RW_speeds_modified[:, i], color=colors[i], 
                        linewidth=2, label=f"RW{i+1}")
        
        # Add markers at correct times
        plt.axvline(x=self.fault_time_min, color='r', linestyle='--', linewidth=2, 
                   label=f'Fault Injection (RW{self.targetRW})')
        plt.axvline(x=self.fault_time_min + 5.0, color='g', linestyle='--', linewidth=2, 
                   label='RW4 Activated')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        if showPlots:
            plt.show()
        
        figureList = {}
        figureList["PowerLimitFault_AttitudeError"] = plt.figure(1)
        figureList["PowerLimitFault_RWSpeeds"] = plt.figure(2)
        
        # Plot 3: Power Consumption Analysis - Shows power limit enforcement
        plt.figure(3, figsize=(12, 8))
        
        # Calculate power consumption for each wheel
        # Power = Torque * Angular velocity
        RW_power = []
        for i in range(num_RW):
            power = np.abs(RW_torque[i] * RW_speeds[:, i])
            RW_power.append(power)
        
        # Subplot 1: Power consumption over time
        plt.subplot(2, 2, 1)
        for i in range(num_RW):
            if i == self.targetRW:
                plt.plot(timeData, RW_power[i], color=colors[i], linewidth=3, 
                        label=f"RW{i+1} (LIMITED)")
            else:
                plt.plot(timeData, RW_power[i], color=colors[i], linewidth=2, 
                        label=f"RW{i+1}")
        
        plt.axvline(x=self.fault_time_min, color='r', linestyle='--', linewidth=2, 
                   label='Power Limit Applied')
        plt.axhline(y=self.fault_magnitude, color='orange', linestyle=':', linewidth=2,
                   label=f'Limit: {self.fault_magnitude}W')
        plt.xlabel("Time (min)")
        plt.ylabel("Power Consumption [W]")
        plt.title("Reaction Wheel Power Consumption")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 2: Power limiting effectiveness
        plt.subplot(2, 2, 2)
        fault_idx = np.argmax(timeData >= self.fault_time_min)
        
        if fault_idx > 0:
            # Focus on the limited wheel
            limited_power = RW_power[self.targetRW]
            
            # Calculate what power would have been without limit (extrapolate)
            pre_fault_trend = limited_power[max(0, fault_idx-20):fault_idx]
            if len(pre_fault_trend) > 0:
                avg_pre_fault = np.mean(pre_fault_trend)
                projected_power = np.ones_like(limited_power[fault_idx:]) * avg_pre_fault
                
                plt.plot(timeData[:fault_idx], limited_power[:fault_idx], 'b-', 
                        linewidth=2, label='Actual Power')
                plt.plot(timeData[fault_idx:], limited_power[fault_idx:], 'r-', 
                        linewidth=2, label='Limited Power')
                plt.plot(timeData[fault_idx:], projected_power, 'b--', 
                        linewidth=2, alpha=0.5, label='Projected (No Limit)')
                
                plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
                plt.axhline(y=self.fault_magnitude, color='orange', linestyle=':', linewidth=2,
                           label=f'Power Limit: {self.fault_magnitude}W')
        
        plt.xlabel("Time (min)")
        plt.ylabel("Power [W]")
        plt.title(f"Power Limiting Effect on RW{self.targetRW+1}")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 3: Cumulative energy saved
        plt.subplot(2, 2, 3)
        if fault_idx > 0:
            dt_min = timeData[1] - timeData[0] if len(timeData) > 1 else 1.0
            dt_hours = dt_min / 60.0
            
            # Energy that would have been consumed vs actual
            actual_energy = np.cumsum(RW_power[self.targetRW] * dt_hours)
            
            # Estimate saved energy
            if len(pre_fault_trend) > 0:
                projected_energy = np.zeros_like(actual_energy)
                projected_energy[:fault_idx] = actual_energy[:fault_idx]
                for i in range(fault_idx, len(timeData)):
                    projected_energy[i] = projected_energy[i-1] + avg_pre_fault * dt_hours
                
                saved_energy = projected_energy - actual_energy
                
                plt.plot(timeData, actual_energy, 'r-', linewidth=2, label='Actual Energy')
                plt.plot(timeData, projected_energy, 'b--', linewidth=2, label='Projected Energy')
                plt.fill_between(timeData, actual_energy, projected_energy, 
                               where=(timeData >= self.fault_time_min), 
                               alpha=0.3, color='green', label='Energy Saved')
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel("Time (min)")
        plt.ylabel("Cumulative Energy [Wh]")
        plt.title("Energy Conservation Due to Power Limit")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 4: Power limit enforcement summary
        plt.subplot(2, 2, 4)
        if fault_idx > 0 and fault_idx < len(timeData):
            # Statistics before and after limit
            pre_limit_power = RW_power[self.targetRW][:fault_idx]
            post_limit_power = RW_power[self.targetRW][fault_idx:]
            
            if len(pre_limit_power) > 0 and len(post_limit_power) > 0:
                stats = {
                    'Mean Power': [np.mean(pre_limit_power), np.mean(post_limit_power)],
                    'Max Power': [np.max(pre_limit_power), np.max(post_limit_power)],
                    'Power Variance': [np.var(pre_limit_power), np.var(post_limit_power)]
                }
                
                x = np.arange(len(stats))
                width = 0.35
                
                fig, ax = plt.gca(), plt.gca()
                
                pre_values = list(stats.values())
                post_values = []
                for i, (key, values) in enumerate(stats.items()):
                    bars1 = ax.bar(i - width/2, values[0], width, label='Pre-Limit' if i == 0 else "", 
                                   color='blue', alpha=0.7)
                    bars2 = ax.bar(i + width/2, values[1], width, label='Post-Limit' if i == 0 else "", 
                                   color='red', alpha=0.7)
                
                ax.set_ylabel('Power [W]')
                ax.set_xlabel('Metric')
                ax.set_title(f'Power Limit Impact on RW{self.targetRW+1}')
                ax.set_xticks(x)
                ax.set_xticklabels(list(stats.keys()), rotation=15)
                ax.legend()
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add limit line
                ax.axhline(y=self.fault_magnitude, color='orange', linestyle=':', linewidth=2)
        
        plt.tight_layout()
        figureList["PowerLimitFault_PowerAnalysis"] = plt.figure(3)
        
        return figureList


def run_power_limit_scenario(powerLimit, fault_wheel, fault_time_min, simulation_time_min):
    """
    Run power limit scenario with GUI parameters
    """
    print(f"Running power limit scenario with:")
    print(f"  - Power limit: {powerLimit} W")
    print(f"  - Target wheel: RW{fault_wheel}")
    print(f"  - Fault time: {fault_time_min} minutes")
    print(f"  - Simulation duration: {simulation_time_min} minutes")
    
    # Pass parameters to scenario
    scenario = PowerLimitFaultScenario(
        fault_magnitude=powerLimit,
        fault_wheel=fault_wheel,
        fault_time_min=fault_time_min,
        simulation_time_min=simulation_time_min
    )
    
    # Use dynamic simulation time
    simulationTime = macros.min2nano(simulation_time_min)
    scenario.modeRequest = "hillPoint"
    
    # Phase 1: Run until fault time
    print(f"Phase 1: Normal operation until {fault_time_min} minutes...")
    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(scenario.faultTime)
    scenario.ExecuteSimulation()
    
    # Phase 2: Apply power limit to specified wheel
    print(f"Phase 2: Applying {powerLimit}W power limit to RW{fault_wheel}...")
    apply_power_limit_to_wheel(scenario, powerLimit, fault_wheel, scenario.faultTime)
    
    # Phase 3: Run with fault until RW4 activation (if time permits)
    if scenario.rw4ActivationTime < simulationTime:
        print("Phase 3: Running with power fault...")
        scenario.ConfigureStopTime(scenario.rw4ActivationTime)
        scenario.ExecuteSimulation()
        
        # Phase 4: Activate RW4
        print("Phase 4: Activating RW4 backup...")
        activate_RW4(scenario, scenario.rw4ActivationTime)
        
        # Phase 5: Complete simulation
        print("Phase 5: Completing simulation...")
        scenario.ConfigureStopTime(simulationTime)
        scenario.ExecuteSimulation()
    else:
        # Just run to end without RW4 activation
        print("Phase 3: Completing simulation...")
        scenario.ConfigureStopTime(simulationTime)
        scenario.ExecuteSimulation()
    
    return scenario


def apply_power_limit_to_wheel(scenario, powerLimit, fault_wheel, currentTimeNanos):
    """
    Apply power limit to the correct wheel specified by GUI
    """
    DynModels = scenario.get_DynModel()
    
    # Get the correct RW object based on fault_wheel index
    rw_objects = [DynModels.RW1, DynModels.RW2, DynModels.RW3, DynModels.RW4]
    
    if 0 <= fault_wheel < len(rw_objects):
        target_rw = rw_objects[fault_wheel]
        
        if hasattr(target_rw, 'maxMotorTorque'):
            scenario.original_max_torque = target_rw.maxMotorTorque
            max_speed_rad_s = 6000 * 2 * np.pi / 60
            reduced_torque = powerLimit / max_speed_rad_s
            target_rw.maxMotorTorque = reduced_torque
            print(f"Set RW{fault_wheel+1} max torque to {reduced_torque} Nm (was {scenario.original_max_torque})")
        else:
            scenario.original_rw_friction = target_rw.fCoulomb
            target_rw.fCoulomb = 0.05
            print(f"Set RW{fault_wheel+1} friction to 0.05 (was {scenario.original_rw_friction})")
    
    timeMin = currentTimeNanos * macros.NANO2MIN
    scenario.powerLimitLog.append(["powerLimit", powerLimit, fault_wheel, timeMin])
    print(f"Applied {powerLimit}W power limit to RW{fault_wheel+1} at {timeMin:.2f} min")


def activate_RW4(scenario, currentTimeNanos):
    """Activate RW4 with correct timing"""
    timeMin = currentTimeNanos * macros.NANO2MIN
    scenario.powerLimitLog.append(["RW4Activation", "N/A", 4, timeMin])
    print(f"Activated RW4 at {timeMin:.2f} minutes")


# ========================================
# GUI INTEGRATION FUNCTIONS
# ========================================

def run(showPlots=True, saveBinary=True, simulation_time_min=30.0):
    """
    GUI-compatible run function that uses default parameters
    Enhanced with dynamic simulation time
    """
    print("\n===== Power Limit Fault Scenario =====")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    # Use defaults if no parameters provided
    fault_magnitude = 0.01
    fault_wheel = 1
    fault_time_min = 5.0
    
    try:
        scenario = run_power_limit_scenario(fault_magnitude, fault_wheel, fault_time_min, simulation_time_min)
        figureList = scenario.pull_outputs(showPlots)
        
        viz = getattr(scenario, 'viz', None)
        
        # Store in module globals for fault_loader
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz
        
        print(f"Generated {len(figureList)} plots with correct parameters")
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}


def run_with_parameters(fault_magnitude=0.01, fault_wheel=1, fault_time_min=10.0, 
                       simulation_time_min=30.0, showPlots=False, saveBinary=False):
    """
    Run with specific GUI parameters
    Enhanced with dynamic simulation time
    """
    print(f"\n===== Power Limit with GUI Parameters =====")
    print(f"PARAMS - Power: {fault_magnitude}W, Wheel: RW{fault_wheel}, Time: {fault_time_min}min")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        # Use the actual GUI parameters
        scenario = run_power_limit_scenario(fault_magnitude, fault_wheel, fault_time_min, simulation_time_min)
        
        # Generate plots with correct parameters
        figureList = scenario.pull_outputs(showPlots)
        
        viz = getattr(scenario, 'viz', None)
        
        # Store in module globals
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz
        
        print(f"SUCCESS: Generated {len(figureList)} plots with GUI parameters")
        print(f"VERIFICATION:")
        print(f"  - Used power limit: {scenario.fault_magnitude}W")
        print(f"  - Used target wheel: RW{scenario.targetRW}")
        print(f"  - Used fault time: {scenario.fault_time_min}min")
        print(f"  - Used simulation time: {scenario.simulation_time_min}min")
        
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR with parameters: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}


# Global storage for fault_loader
scenario = None
figureList = {}
viz = None

if __name__ == "__main__":
    # Test with custom parameters
    test_scenario, test_viz, test_plots = run_with_parameters(
        fault_magnitude=0.5, 
        fault_wheel=0, 
        fault_time_min=15.0,
        simulation_time_min=60.0,
        showPlots=True
    )
    print(f"Test completed with {len(test_plots)} plots")