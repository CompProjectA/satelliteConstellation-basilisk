#!/usr/bin/env python
"""
friction_fault.py - Enhanced with dynamic simulation time support

This module simulates friction faults in reaction wheels with BSK_Plotting integration.
"""

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
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import BSK_Plotting as BSK_plt
    BSK_PLOTTING_AVAILABLE = True
    print("Successfully imported BSK_Plotting for friction fault")
except ImportError as e:
    print(f"Warning: Could not import BSK_Plotting: {e}")
    print("Will use standard matplotlib instead")
    BSK_PLOTTING_AVAILABLE = False

class FrictionFaultScenario(BSKSim, BSKScenario):
    """
    Friction fault scenario with GUI parameters and BSK_Plotting support
    Enhanced with dynamic simulation time support
    """
    def __init__(self, fault_magnitude=0.0005, fault_wheel=3, fault_time_min=10.0, simulation_time_min=30.0):
        super(FrictionFaultScenario, self).__init__()
        
        # Store GUI parameters
        self.fault_magnitude = fault_magnitude
        self.fault_wheel_number = fault_wheel
        self.fault_time_min = fault_time_min
        self.simulation_time_min = simulation_time_min
        
        # Convert fault time to nanoseconds
        self.oneTimeFaultTime = macros.min2nano(fault_time_min)
        
        print(f"FrictionFaultScenario initialized with GUI parameters:")
        print(f"  - Fault magnitude: {fault_magnitude} N⋅m")
        print(f"  - Target wheel: RW{fault_wheel}")
        print(f"  - Fault time: {fault_time_min} minutes")
        print(f"  - Simulation duration: {simulation_time_min} minutes")
        
        # Original initialization
        self.name = 'FrictionFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"
        self.cameraLocation = [0.0, 3.0, 0.0]

        # Set the dynamics and flight software models
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        # Set up initial conditions and log desired outputs
        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection configuration
        self.oneTimeRWFaultFlag = 1
        self.repeatRWFaultFlag = 1
        self.get_DynModel().RWFaultLog = []

    def configure_initial_conditions(self):
        """Set initial orbital conditions"""
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
        """Set up message recording for RW speeds and friction torques"""
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

    def calculate_temperatures(self, rw_speeds, rw_friction):
        """Calculate temperature evolution based on friction"""
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
        """Enhanced temperature plot with GUI parameters"""
        plt.figure(figsize=(12, 6))
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            # Highlight faulty wheel
            if idx == self.fault_wheel_number:
                plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                        label=f'RW {idx} (FAULTY)', linewidth=4)
            else:
                plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                        label=f'RW {idx}', linewidth=2)

        # Add fault time marker
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('Temperature [°C]', fontsize=12)
        plt.title(f'RW Temperatures (No Cooling)\nFriction Fault: RW{self.fault_wheel_number}, +{self.fault_magnitude} N⋅m at {self.fault_time_min}min', 
                 fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning (30°C)')
        plt.axhline(y=50, color='red', linestyle='--', alpha=0.7, label='Critical (50°C)')
        plt.xlim(0, self.simulation_time_min)
        plt.tight_layout()

    def plot_rw_C_temperature(self, timeData, RW_temperatures, numRW):
        """Enhanced cooling temperature plot with GUI parameters"""
        plt.figure(figsize=(12, 6))
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            # Highlight faulty wheel
            if idx == self.fault_wheel_number:
                plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                        label=f'RW {idx} (FAULTY)', linewidth=4)
            else:
                plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                        label=f'RW {idx}', linewidth=2)

        # Add fault time marker
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('Temperature [°C]', fontsize=12)
        plt.title(f'RW Temperatures (With Cooling)\nFriction Fault: RW{self.fault_wheel_number}, +{self.fault_magnitude} N⋅m at {self.fault_time_min}min', 
                 fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning (30°C)')
        plt.xlim(0, self.simulation_time_min)
        plt.tight_layout()

    def pull_outputs(self, showPlots):
        """
        Use BSK_Plotting when available, with GUI parameters integrated
        """
        print("Generating friction fault analysis plots with GUI parameters...")

        # FSW process outputs, remove first data point as it is before FSW is called
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)

        # Extract recorded RW data (speed and friction)
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = []
        for i in range(numRW):
            try:
                RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))
            except:
                # Create placeholder if friction data not available
                timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
                RW_friction.append(np.zeros_like(timeData))

        # Estimate RW temperatures based on speed and friction
        self.no_cooling, self.with_cooling = self.calculate_temperatures(RW_speeds, RW_friction)

        # Use BSK_Plotting if available, otherwise matplotlib
        if BSK_PLOTTING_AVAILABLE:
            return self._generate_plots_with_bsk(RW_speeds, RW_friction, sigma_BR, showPlots)
        else:
            return self._generate_plots_with_matplotlib(RW_speeds, RW_friction, sigma_BR, showPlots)

    def _generate_plots_with_bsk(self, RW_speeds, RW_friction, sigma_BR, showPlots):
        """Generate plots using BSK_Plotting"""
        print("Using BSK_Plotting for friction fault plot generation...")
        
        # Use BSK_Plotting
        BSK_plt.clear_all_plots()
        timeData = np.delete(self.msgRecList[self.attGuidName].times(), 0, 0) * macros.NANO2MIN
        numRW = 4

        # Temporarily store fault info for BSK_plt functions to use
        if hasattr(BSK_plt, 'fault_info'):
            BSK_plt.fault_info = {
                'fault_time': self.fault_time_min,
                'fault_wheel': self.fault_wheel_number,
                'fault_magnitude': self.fault_magnitude
            }

        # Use BSK_Plotting functions
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.get_DynModel().RWFaultLog)
        BSK_plt.plot_attitude_error(timeData, sigma_BR)

        # Add custom temperature plots with GUI parameters
        self.plot_rw_temperature(timeData, self.no_cooling, numRW)
        self.plot_rw_C_temperature(timeData, self.with_cooling, numRW)

        # Add custom friction analysis plot with GUI parameters
        self._add_friction_analysis_plot(timeData, RW_speeds, RW_friction, numRW)

        # Return or show/save figures
        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            # Updated figure names to include GUI parameters
            figureNames = [
                f"RWSpeeds_RW{self.fault_wheel_number}Fault", 
                f"RWFriction_Mag{self.fault_magnitude}", 
                f"RWTemperatures_NoCooling_RW{self.fault_wheel_number}", 
                f"RWTemperatures_WithCooling_RW{self.fault_wheel_number}",
                "attitudeErrorNorm",
                f"FrictionAnalysis_RW{self.fault_wheel_number}"
            ]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        print(f"Generated {len(figureList)} friction fault plots using BSK_Plotting")
        return figureList

    def _generate_plots_with_matplotlib(self, RW_speeds, RW_friction, sigma_BR, showPlots):
        """Fallback plotting using standard matplotlib"""
        print("Using matplotlib for friction fault plot generation...")
        
        plt.close('all')
        timeData = np.delete(self.msgRecList[self.attGuidName].times(), 0, 0) * macros.NANO2MIN
        figureList = {}

        # Basic RW speeds plot
        fig1 = plt.figure(figsize=(12, 8))
        colors = ['blue', 'green', 'red', 'orange']
        
        for i in range(4):
            if i == self.fault_wheel_number:
                plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=4, 
                        label=f"RW{i} (FAULTY)")
            else:
                plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=2, 
                        label=f"RW{i}")
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('RW Speed [rad/s]', fontsize=12)
        plt.title(f'Friction Fault RW Speeds\n(RW{self.fault_wheel_number}, +{self.fault_magnitude} N⋅m at {self.fault_time_min}min)', 
                 fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        plt.tight_layout()
        figureList["FrictionFault_RWSpeeds"] = fig1

        # Generate temperature plots
        self.plot_rw_temperature(timeData, self.no_cooling, 4)
        figureList["FrictionFault_RWTemperatures"] = plt.gcf()

        self.plot_rw_C_temperature(timeData, self.with_cooling, 4)
        figureList["FrictionFault_RWTemperatures_Cooling"] = plt.gcf()

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        print(f"Generated {len(figureList)} friction fault plots using matplotlib")
        return figureList

    def _add_friction_analysis_plot(self, timeData, RW_speeds, RW_friction, numRW):
        """Add detailed friction analysis plot with GUI parameters"""
        plt.figure(figsize=(14, 8))
        
        # Friction comparison plot
        plt.subplot(2, 2, 1)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(numRW):
            if i == self.fault_wheel_number:
                plt.plot(timeData, RW_friction[i], color=colors[i], linewidth=4, 
                        label=f"RW{i+1} Friction (FAULTY)")
            else:
                plt.plot(timeData, RW_friction[i], color=colors[i], linewidth=2, 
                        label=f"RW{i+1} Friction")
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.xlabel('Time [min]')
        plt.ylabel('Friction Torque [N⋅m]')
        plt.title(f'Friction Analysis: +{self.fault_magnitude} N⋅m to RW{self.fault_wheel_number+1}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Power consumption due to friction
        plt.subplot(2, 2, 2)
        for i in range(numRW):
            # Power = Torque * Angular velocity
            power_friction = np.abs(RW_friction[i] * RW_speeds[:, i])
            if i == self.fault_wheel_number:
                plt.plot(timeData, power_friction, color=colors[i], linewidth=4, 
                        label=f"RW{i+1} Power (FAULTY)")
            else:
                plt.plot(timeData, power_friction, color=colors[i], linewidth=2, 
                        label=f"RW{i+1} Power")
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('Power Loss [W]')
        plt.title('Power Consumption Due to Friction')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Cumulative energy loss
        plt.subplot(2, 2, 3)
        dt_min = timeData[1] - timeData[0] if len(timeData) > 1 else 1.0
        dt_hours = dt_min / 60.0
        
        for i in range(numRW):
            power_friction = np.abs(RW_friction[i] * RW_speeds[:, i])
            energy_loss = np.cumsum(power_friction * dt_hours)  # Wh
            
            if i == self.fault_wheel_number:
                plt.plot(timeData, energy_loss, color=colors[i], linewidth=4, 
                        label=f"RW{i+1} Energy (FAULTY)")
            else:
                plt.plot(timeData, energy_loss, color=colors[i], linewidth=2, 
                        label=f"RW{i+1} Energy")
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('Cumulative Energy Loss [Wh]')
        plt.title('Total Energy Lost to Friction')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Fault impact analysis
        plt.subplot(2, 2, 4)
        # Calculate the difference in friction after fault
        fault_idx = np.argmax(timeData >= self.fault_time_min)
        
        if fault_idx > 0 and fault_idx < len(timeData):
            pre_fault_friction = np.mean(RW_friction[self.fault_wheel_number][:fault_idx])
            post_fault_friction = np.mean(RW_friction[self.fault_wheel_number][fault_idx:])
            friction_increase = post_fault_friction - pre_fault_friction
            
            # Bar chart showing friction increase
            categories = ['Pre-Fault', 'Post-Fault', 'Increase']
            values = [pre_fault_friction, post_fault_friction, friction_increase]
            bars = plt.bar(categories, values, color=['green', 'red', 'orange'])
            
            # Add value labels with scientific notation for small values
            for bar, value in zip(bars, values):
                height = bar.get_height()
                if abs(value) < 0.01:
                    label_text = f'{value:.2e}'  # Scientific notation for small values
                else:
                    label_text = f'{value:.5f}'
                plt.text(bar.get_x() + bar.get_width()/2., height + abs(height) * 0.01,
                        label_text, ha='center', va='bottom')
            
            plt.ylabel('Friction Torque [N⋅m]')
            plt.title(f'Friction Fault Impact on RW{self.fault_wheel_number+1}')
            plt.grid(True, alpha=0.3, axis='y')
            
            # Set y-axis to show small values properly
            plt.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        
        plt.tight_layout()

    def runScenario(self, saveBinary=True):
        """
        Run scenario with binary file saving (no live viz)
        Using dynamic simulation time
        """
        # Define simulation duration based on parameter
        simulationTime = macros.min2nano(self.simulation_time_min)
        self.modeRequest = "hillPoint"

        print(f"Running friction fault scenario:")
        print(f"  - Fault magnitude: {self.fault_magnitude} N⋅m")
        print(f"  - Target wheel: RW{self.fault_wheel_number}")
        print(f"  - Fault time: {self.fault_time_min} minutes")
        print(f"  - Simulation duration: {self.simulation_time_min} minutes")

        # Add friction fault event with GUI parameters
        self.createNewEvent(
            "addOneTimeRWFault",
            self.get_FswModel().processTasksTimeStep,
            True,
            ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
            [f"self.get_DynModel().AddRWFault('friction',{self.fault_magnitude},{self.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
             "self.oneTimeRWFaultFlag=0"]
        )

        # Setup visualization to save binary file (no live streaming)
        viz = None
        if vizSupport.vizFound and saveBinary:
            try:
                # Create binary filename with fault parameters
                binary_filename = f"friction_fault_rw{self.fault_wheel_number}_t{int(self.fault_time_min)}_f{int(self.fault_magnitude*10000)}"
                
                viz = vizSupport.enableUnityVisualization(
                    self,
                    self.get_DynModel().taskName,
                    self.get_DynModel().scObject,
                    rwEffectorList=self.get_DynModel().rwStateEffector,
                    liveStream=False,      # No live streaming
                    saveFile=binary_filename  # Save to binary file
                )

                # Set up camera
                vizSupport.createStandardCamera(
                    viz,
                    setMode=1,
                    spacecraftName=self.get_DynModel().scObject.ModelTag,
                    fieldOfView=30 * macros.D2R,
                    displayName=f"Friction Fault Camera RW{self.fault_wheel_number}",
                    pointingVector_B=[0, 0, 0],
                    position_B=self.cameraLocation
                )
                
                print(f"Visualization configured for binary file saving")

            except Exception as e:
                print(f"WARNING: Visualization setup failed: {e}")
                viz = None

        # Run the simulation
        self.InitializeSimulation()
        self.ConfigureStopTime(simulationTime)
        self.ExecuteSimulation()
        
        return viz

def run(showPlots=True, saveBinary=True, simulation_time_min=30.0):
    """
    GUI-compatible run function using BSK_Plotting
    Enhanced with dynamic simulation time
    """
    print("\n===== Friction Fault Scenario with BSK_Plotting =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        scenario = FrictionFaultScenario(0.0005, 3, 10.0, simulation_time_min)
        viz = scenario.runScenario(saveBinary)
        figureList = scenario.pull_outputs(showPlots)

        # Store in module globals for fault_loader
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz

        print(f"Generated {len(figureList)} friction fault analysis plots")
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR: Friction fault failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}

def run_with_parameters(fault_magnitude=0.0005, fault_wheel=3, fault_time_min=10.0, 
                       simulation_time_min=30.0, showPlots=False, saveBinary=False):
    """
    Run friction fault with GUI parameters using BSK_Plotting
    Enhanced with dynamic simulation time
    """
    print(f"\n===== Friction Fault with GUI Parameters (BSK_Plotting) =====")
    print(f"PARAMS - Magnitude: {fault_magnitude} N⋅m, Wheel: RW{fault_wheel}, Time: {fault_time_min}min")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        # Use actual GUI parameters including simulation time
        scenario = FrictionFaultScenario(fault_magnitude, fault_wheel, fault_time_min, simulation_time_min)
        
        # Run simulation with parameters
        viz = scenario.runScenario(saveBinary)
        
        # Generate plots with correct parameters
        figureList = scenario.pull_outputs(showPlots)
        
        # Store in module globals for fault_loader
        import sys
        current_module = sys.modules[__name__]
        current_module.scenario = scenario
        current_module.figureList = figureList
        current_module.viz = viz
        
        print(f"SUCCESS: Generated {len(figureList)} friction fault plots with GUI parameters")
        print(f"VERIFICATION:")
        print(f"  - Used fault magnitude: {scenario.fault_magnitude} N⋅m")
        print(f"  - Used target wheel: RW{scenario.fault_wheel_number}")
        print(f"  - Used fault time: {scenario.fault_time_min} minutes")
        print(f"  - Used simulation time: {scenario.simulation_time_min} minutes")
        print(f"  - Used BSK_Plotting: {BSK_PLOTTING_AVAILABLE}")
        
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR with friction fault parameters: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}

# Global storage for fault_loader
scenario = None
figureList = {}
viz = None

if __name__ == "__main__":
    run(True, True, 30.0)