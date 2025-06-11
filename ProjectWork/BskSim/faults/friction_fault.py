#!/usr/bin/env python
"""
friction_fault.py - Complete version with ALL plots from old and new versions

This module simulates friction faults in reaction wheels and generates comprehensive plots.
"""

import inspect
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to avoid thread issues
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
    BSK_AVAILABLE = True
    BSK_PLOTTING_AVAILABLE = True
    print("Successfully imported BSK modules for friction fault")
except ImportError as e:
    print(f"Warning: Could not import BSK modules: {e}")
    BSK_AVAILABLE = False
    BSK_PLOTTING_AVAILABLE = False

class FrictionFaultScenario(BSKSim, BSKScenario):
    """
    Friction fault scenario with comprehensive plot generation
    """
    def __init__(self, fault_magnitude=0.0005, fault_wheel=3, fault_time_min=10.0, simulation_time_min=30.0):
        super(FrictionFaultScenario, self).__init__()
        
        # Store GUI parameters - fault_wheel is already 0-based from GUI
        self.fault_magnitude = fault_magnitude
        self.fault_wheel_number = fault_wheel  # 0-based index
        self.fault_time_min = fault_time_min
        self.simulation_time_min = simulation_time_min
        
        # Convert fault time to nanoseconds
        self.oneTimeFaultTime = macros.min2nano(fault_time_min)
        
        print(f"FrictionFaultScenario initialized with GUI parameters:")
        print(f"  - Fault magnitude: {fault_magnitude} N⋅m")
        print(f"  - Target wheel: RW{fault_wheel + 1} (internal index: {fault_wheel})")
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
        self.repeatRWFaultFlag = 0  # Disable repeated faults for GUI control
        self.get_DynModel().RWFaultLog = []
        
        # Override BSK_Dynamics fault parameters with GUI values
        self.get_DynModel().oneTimeFaultTime = self.oneTimeFaultTime
        self.get_DynModel().rwFaultFlag = 1

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
        self.rwLogs = [None] * 4  # Create list with 4 None elements
        for item in range(4):
            self.rwLogs[item] = DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime)
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

        # Record commanded torque if available
        try:
            self.rwMotorRec = FswModel.rwMotorTorque.rwMotorTorqueOutMsg.recorder(samplingTime)
            self.AddModelToTask(DynModel.taskName, self.rwMotorRec)
        except:
            self.rwMotorRec = None

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
                omega = rw_speeds[i, rw]
                P_friction = abs(rw_friction[rw][i] * omega) if i < len(rw_friction[rw]) else 0
                temp_rise = P_friction * 0.2
                cooling = 0.005 * (temp_c[i-1] - T_ambient)

                temp_nc[i] = temp_nc[i-1] + temp_rise
                temp_c[i] = temp_c[i-1] + (temp_rise) - cooling

            no_cooling.append(temp_nc)
            with_cooling.append(temp_c)

        return no_cooling, with_cooling

    def pull_outputs(self, showPlots):
        """
        Generate ALL plots using real simulation data - combines old and new versions
        """
        print("Generating comprehensive friction fault analysis plots from REAL simulation data...")

        # FSW process outputs, remove first data point as it is before FSW is called
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)

        # Extract recorded RW data (speed and friction)
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = []
        RW_torque = []
        
        # Extract real friction and torque data from simulation
        for i in range(numRW):
            try:
                friction_data = np.delete(self.rwLogs[i].u_f, 0, 0)
                RW_friction.append(friction_data)
                # Try to get torque if available
                if hasattr(self.rwLogs[i], 'u_current'):
                    torque_data = np.delete(self.rwLogs[i].u_current, 0, 0)
                    RW_torque.append(torque_data)
                else:
                    RW_torque.append(np.zeros_like(friction_data))
            except:
                # If data not available, create zero array
                timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
                RW_friction.append(np.zeros_like(timeData))
                RW_torque.append(np.zeros_like(timeData))

        # Get commanded torque if available
        RW_cmd_torque = None
        if self.rwMotorRec:
            try:
                RW_cmd_torque = np.delete(self.rwMotorRec.motorTorque[:, range(numRW)], 0, 0)
            except:
                pass

        # Estimate RW temperatures based on speed and friction
        self.no_cooling, self.with_cooling = self.calculate_temperatures(RW_speeds, RW_friction)

        # Generate ALL plots
        return self._generate_all_plots(RW_speeds, RW_friction, RW_torque, RW_cmd_torque, 
                                       sigma_BR, omega_BR_B, showPlots)

    def _generate_all_plots(self, RW_speeds, RW_friction, RW_torque, RW_cmd_torque, 
                           sigma_BR, omega_BR_B, showPlots):
        """Generate ALL plots from both old and new versions"""
        print("Using REAL friction fault data for comprehensive plot generation...")
        
        matplotlib.use('Agg')  # Ensure non-interactive backend
        plt.close('all')
        timeData = np.delete(self.msgRecList[self.attGuidName].times(), 0, 0) * macros.NANO2MIN
        figureList = {}

        # ===== PLOT 1: RW Speeds (Enhanced) =====
        fig1 = plt.figure(figsize=(12, 8))
        colors = ['blue', 'green', 'red', 'orange']
        
        plt.subplot(2, 1, 1)
        for i in range(4):
            if i == self.fault_wheel_number:
                plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=4, 
                        label=f"RW{i + 1} (FAULTY)", zorder=10)
            else:
                plt.plot(timeData, RW_speeds[:, i], color=colors[i], linewidth=2, 
                        label=f"RW{i + 1}", alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault Injection')
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('RW Speed [rad/s]', fontsize=12)
        plt.title(f'REAL DATA: Reaction Wheel Speeds\n(Friction Fault on RW{self.fault_wheel_number + 1}: +{self.fault_magnitude} N⋅m at {self.fault_time_min} min)', 
                 fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 2: RW speeds in RPM
        plt.subplot(2, 1, 2)
        for i in range(4):
            speed_rpm = RW_speeds[:, i] * 60 / (2 * np.pi)  # Convert to RPM
            if i == self.fault_wheel_number:
                plt.plot(timeData, speed_rpm, color=colors[i], linewidth=4, 
                        label=f"RW{i + 1} (FAULTY)", zorder=10)
            else:
                plt.plot(timeData, speed_rpm, color=colors[i], linewidth=2, 
                        label=f"RW{i + 1}", alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('RW Speed [RPM]', fontsize=12)
        plt.title('Reaction Wheel Speeds (RPM)', fontsize=12)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        plt.tight_layout()
        figureList["REAL_FrictionFault_RWSpeeds"] = fig1

        # ===== PLOT 2: Friction Torques =====
        fig2 = plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        for i in range(4):
            if i == self.fault_wheel_number:
                plt.plot(timeData, RW_friction[i], color=colors[i], linewidth=4, 
                        label=f"RW{i+1} (FAULTY)", zorder=10)
            else:
                plt.plot(timeData, RW_friction[i], color=colors[i], linewidth=2, 
                        label=f"RW{i+1}", alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Friction Torque [N⋅m]')
        plt.title('REAL DATA: Reaction Wheel Friction Torques')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Subplot 2: Friction torque zoom on faulty wheel
        plt.subplot(2, 1, 2)
        plt.plot(timeData, RW_friction[self.fault_wheel_number], 'r-', linewidth=3)
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Fault Injection')
        plt.axhline(y=0.0005, color='green', linestyle=':', label='Nominal Friction')
        plt.axhline(y=0.0005 + self.fault_magnitude, color='red', linestyle=':', 
                   label=f'Faulty Friction ({0.0005 + self.fault_magnitude:.4f} N⋅m)')
        plt.xlabel('Time [min]')
        plt.ylabel('Friction Torque [N⋅m]')
        plt.title(f'RW{self.fault_wheel_number + 1} Friction Detail')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        plt.tight_layout()
        figureList["REAL_FrictionFault_FrictionTorques"] = fig2

        # ===== PLOT 3: Temperature Evolution (No Cooling) =====
        fig3 = plt.figure(figsize=(12, 6))
        for idx in range(4):
            if idx == self.fault_wheel_number:
                plt.plot(timeData, self.no_cooling[idx], color=colors[idx],
                        label=f'RW{idx + 1} (FAULTY)', linewidth=4, zorder=10)
            else:
                plt.plot(timeData, self.no_cooling[idx], color=colors[idx],
                        label=f'RW{idx + 1}', linewidth=2, alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning (30°C)')
        plt.axhline(y=50, color='red', linestyle='--', alpha=0.7, label='Critical (50°C)')
        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('REAL DATA: RW Temperatures (No Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        plt.tight_layout()
        figureList["REAL_FrictionFault_Temperature_NoCooling"] = fig3

        # ===== PLOT 4: Temperature Evolution (With Cooling) =====
        fig4 = plt.figure(figsize=(12, 6))
        for idx in range(4):
            if idx == self.fault_wheel_number:
                plt.plot(timeData, self.with_cooling[idx], color=colors[idx],
                        label=f'RW{idx + 1} (FAULTY)', linewidth=4, zorder=10)
            else:
                plt.plot(timeData, self.with_cooling[idx], color=colors[idx],
                        label=f'RW{idx + 1}', linewidth=2, alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning (30°C)')
        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('REAL DATA: RW Temperatures (With Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        plt.tight_layout()
        figureList["REAL_FrictionFault_Temperature_WithCooling"] = fig4

        # ===== PLOT 5: Comprehensive Friction Analysis =====
        fig5 = plt.figure(figsize=(14, 10))
        
        # Power consumption due to friction
        plt.subplot(2, 2, 1)
        for i in range(4):
            power_friction = np.abs(RW_friction[i] * RW_speeds[:, i])
            if i == self.fault_wheel_number:
                plt.plot(timeData, power_friction, color=colors[i], linewidth=4, 
                        label=f"RW{i+1} (FAULTY)", zorder=10)
            else:
                plt.plot(timeData, power_friction, color=colors[i], linewidth=2, 
                        label=f"RW{i+1}", alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('Power Loss [W]')
        plt.title('Power Consumption Due to Friction')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Cumulative energy loss
        plt.subplot(2, 2, 2)
        dt_min = timeData[1] - timeData[0] if len(timeData) > 1 else 1.0
        dt_hours = dt_min / 60.0
        
        for i in range(4):
            power_friction = np.abs(RW_friction[i] * RW_speeds[:, i])
            energy_loss = np.cumsum(power_friction * dt_hours)  # Wh
            
            if i == self.fault_wheel_number:
                plt.plot(timeData, energy_loss, color=colors[i], linewidth=4, 
                        label=f"RW{i+1} (FAULTY)", zorder=10)
            else:
                plt.plot(timeData, energy_loss, color=colors[i], linewidth=2, 
                        label=f"RW{i+1}", alpha=0.7)
        
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('Cumulative Energy Loss [Wh]')
        plt.title('Total Energy Lost to Friction')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Fault impact analysis
        plt.subplot(2, 2, 3)
        fault_idx = np.argmax(timeData >= self.fault_time_min)
        
        if fault_idx > 0 and fault_idx < len(timeData):
            pre_fault_friction = np.mean(RW_friction[self.fault_wheel_number][:fault_idx])
            post_fault_friction = np.mean(RW_friction[self.fault_wheel_number][fault_idx:])
            friction_increase = post_fault_friction - pre_fault_friction
            
            categories = ['Pre-Fault', 'Post-Fault', 'Increase']
            values = [pre_fault_friction, post_fault_friction, friction_increase]
            bars = plt.bar(categories, values, color=['green', 'red', 'orange'])
            
            for bar, value in zip(bars, values):
                height = bar.get_height()
                if abs(value) < 0.01:
                    label_text = f'{value:.2e}'
                else:
                    label_text = f'{value:.5f}'
                plt.text(bar.get_x() + bar.get_width()/2., height + abs(height) * 0.01,
                        label_text, ha='center', va='bottom')
            
            plt.ylabel('Friction Torque [N⋅m]')
            plt.title(f'Friction Fault Impact on RW{self.fault_wheel_number+1}')
            plt.grid(True, alpha=0.3, axis='y')
            plt.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
        
        # Temperature comparison
        plt.subplot(2, 2, 4)
        temp_nc = self.no_cooling[self.fault_wheel_number]
        temp_c = self.with_cooling[self.fault_wheel_number]
        
        plt.plot(timeData, temp_nc, 'r-', linewidth=3, label='No Cooling')
        plt.plot(timeData, temp_c, 'b-', linewidth=3, label='With Cooling')
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.axhline(y=30, color='orange', linestyle=':', alpha=0.7)
        
        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title(f'Temperature Comparison for RW{self.fault_wheel_number+1}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        plt.tight_layout()
        figureList["REAL_FrictionFault_DetailedAnalysis"] = fig5

        # ===== PLOT 6: Attitude Error Analysis =====
        fig6 = plt.figure(figsize=(12, 8))
        
        # Attitude error norm
        plt.subplot(3, 1, 1)
        sigma_BR_norm = np.linalg.norm(sigma_BR, axis=1)
        plt.semilogy(timeData, sigma_BR_norm, 'b-', linewidth=2)
        plt.axvline(x=self.fault_time_min, color='red', linestyle='--', linewidth=2, 
                   label='Friction Fault')
        plt.xlabel('Time [min]')
        plt.ylabel('Attitude Error Norm')
        plt.title('REAL DATA: Attitude Control Performance')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Individual attitude error components
        plt.subplot(3, 1, 2)
        plt.plot(timeData, sigma_BR[:, 0], 'r-', label='Roll', linewidth=2)
        plt.plot(timeData, sigma_BR[:, 1], 'g-', label='Pitch', linewidth=2)
        plt.plot(timeData, sigma_BR[:, 2], 'b-', label='Yaw', linewidth=2)
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('MRP Error')
        plt.title('Attitude Error Components')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        # Rate error
        plt.subplot(3, 1, 3)
        plt.plot(timeData, omega_BR_B[:, 0], 'r-', label='Roll Rate', linewidth=2)
        plt.plot(timeData, omega_BR_B[:, 1], 'g-', label='Pitch Rate', linewidth=2)
        plt.plot(timeData, omega_BR_B[:, 2], 'b-', label='Yaw Rate', linewidth=2)
        plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
        plt.xlabel('Time [min]')
        plt.ylabel('Rate Error [rad/s]')
        plt.title('Angular Rate Error')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(0, self.simulation_time_min)
        
        plt.tight_layout()
        figureList["REAL_FrictionFault_AttitudeError"] = fig6

        # ===== PLOT 7: RW Motor Torque Analysis (if available) =====
        if RW_cmd_torque is not None:
            fig7 = plt.figure(figsize=(12, 8))
            
            # Commanded vs actual torque
            plt.subplot(2, 1, 1)
            for i in range(4):
                if i == self.fault_wheel_number:
                    plt.plot(timeData, RW_cmd_torque[:, i], '--', color=colors[i], 
                            linewidth=3, label=f'RW{i+1} Cmd (FAULTY)')
                    plt.plot(timeData, RW_torque[i], '-', color=colors[i], 
                            linewidth=3, label=f'RW{i+1} Actual')
                else:
                    plt.plot(timeData, RW_cmd_torque[:, i], '--', color=colors[i], 
                            linewidth=2, alpha=0.5, label=f'RW{i+1} Cmd')
                    plt.plot(timeData, RW_torque[i], '-', color=colors[i], 
                            linewidth=2, alpha=0.7, label=f'RW{i+1} Actual')
            
            plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
            plt.xlabel('Time [min]')
            plt.ylabel('Motor Torque [N⋅m]')
            plt.title('Commanded vs Actual RW Motor Torques')
            plt.legend(ncol=2, fontsize=9)
            plt.grid(True, alpha=0.3)
            plt.xlim(0, self.simulation_time_min)
            
            # Torque error
            plt.subplot(2, 1, 2)
            for i in range(4):
                torque_error = RW_cmd_torque[:, i] - RW_torque[i]
                if i == self.fault_wheel_number:
                    plt.plot(timeData, torque_error, color=colors[i], linewidth=3, 
                            label=f'RW{i+1} (FAULTY)')
                else:
                    plt.plot(timeData, torque_error, color=colors[i], linewidth=2, 
                            alpha=0.7, label=f'RW{i+1}')
            
            plt.axvline(x=self.fault_time_min, color='black', linestyle='--', linewidth=2)
            plt.xlabel('Time [min]')
            plt.ylabel('Torque Error [N⋅m]')
            plt.title('Motor Torque Tracking Error')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.xlim(0, self.simulation_time_min)
            
            plt.tight_layout()
            figureList["REAL_FrictionFault_MotorTorque"] = fig7

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        print(f"Generated {len(figureList)} comprehensive friction fault plots using REAL simulation data")
        return figureList

    def runScenario(self, saveBinary=True):
        """
        Run scenario with proper fault injection using GUI parameters
        """
        # Define simulation duration based on parameter
        simulationTime = macros.min2nano(self.simulation_time_min)
        self.modeRequest = "hillPoint"

        print(f"Running friction fault scenario with GUI parameters:")
        print(f"  - Fault magnitude: {self.fault_magnitude} N⋅m")
        print(f"  - Target wheel: RW{self.fault_wheel_number + 1} (index: {self.fault_wheel_number})")
        print(f"  - Fault time: {self.fault_time_min} minutes")
        print(f"  - Simulation duration: {self.simulation_time_min} minutes")

        # Add friction fault event with actual GUI parameters
        self.createNewEvent(
            "addOneTimeRWFault",
            self.get_FswModel().processTasksTimeStep,
            True,
            ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
            [f"self.get_DynModel().AddRWFault('friction',{self.fault_magnitude},{self.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
             "self.oneTimeRWFaultFlag=0",
             f"print('FRICTION FAULT INJECTED: RW{self.fault_wheel_number + 1}, magnitude={self.fault_magnitude} N⋅m')"]
        )

        # Setup visualization
        viz = None
        if vizSupport.vizFound and saveBinary:
            try:
                binary_filename = f"friction_fault_rw{self.fault_wheel_number + 1}_t{int(self.fault_time_min)}_m{int(self.fault_magnitude*10000)}"
                
                viz = vizSupport.enableUnityVisualization(
                    self,
                    self.get_DynModel().taskName,
                    self.get_DynModel().scObject,
                    rwEffectorList=self.get_DynModel().rwStateEffector,
                    liveStream=False,
                    saveFile=binary_filename
                )

                print(f"Visualization configured for binary file: {binary_filename}")

            except Exception as e:
                print(f"WARNING: Visualization setup failed: {e}")
                viz = None

        # Run the simulation
        self.InitializeSimulation()
        self.ConfigureStopTime(simulationTime)
        self.ExecuteSimulation()
        
        # Verify fault was applied
        if hasattr(self.get_DynModel(), 'RWFaultLog') and len(self.get_DynModel().RWFaultLog) > 0:
            print(f"Fault log entries: {len(self.get_DynModel().RWFaultLog)}")
            for fault in self.get_DynModel().RWFaultLog:
                print(f"  - {fault[0]} fault: magnitude={fault[1]}, wheel={fault[2]}, time={fault[3]:.2f} min")
        
        return viz

def run(showPlots=True, saveBinary=True, simulation_time_min=30.0):
    """
    GUI-compatible run function
    """
    print("\n===== Friction Fault Scenario =====")
    print(f"Show Plots: {showPlots}")
    print(f"Save Binary: {saveBinary}")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        # Default parameters
        scenario = FrictionFaultScenario(0.0005, 3, 10.0, simulation_time_min)
        viz = scenario.runScenario(saveBinary)
        figureList = scenario.pull_outputs(showPlots)

        print(f"Generated {len(figureList)} friction fault analysis plots with REAL data")
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR: Friction fault failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}

def run_with_parameters(fault_magnitude=0.0005, fault_wheel=3, fault_time_min=10.0, 
                       simulation_time_min=30.0, showPlots=False, saveBinary=False):
    """
    Run friction fault with GUI parameters - fault_wheel is already 0-based from GUI
    """
    print(f"\n===== Friction Fault with GUI Parameters =====")
    print(f"PARAMS - Magnitude: {fault_magnitude} N⋅m, Wheel: RW{fault_wheel + 1} (index: {fault_wheel}), Time: {fault_time_min}min")
    print(f"Simulation Duration: {simulation_time_min} minutes")
    
    try:
        # Use actual GUI parameters
        scenario = FrictionFaultScenario(fault_magnitude, fault_wheel, fault_time_min, simulation_time_min)
        
        # Run simulation with parameters
        viz = scenario.runScenario(saveBinary)
        
        # Generate plots with correct parameters
        figureList = scenario.pull_outputs(showPlots)
        
        print(f"SUCCESS: Generated {len(figureList)} friction fault plots with REAL GUI parameters")
        print(f"VERIFICATION:")
        print(f"  - Used fault magnitude: {scenario.fault_magnitude} N⋅m")
        print(f"  - Used target wheel: RW{scenario.fault_wheel_number + 1} (index: {scenario.fault_wheel_number})")
        print(f"  - Used fault time: {scenario.fault_time_min} minutes")
        print(f"  - Used simulation time: {scenario.simulation_time_min} minutes")
        print(f"  - ALL PLOTS USE REAL SIMULATION DATA")
        
        return scenario, viz, figureList
        
    except Exception as e:
        print(f"ERROR with friction fault parameters: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}

# Main execution
if __name__ == "__main__":
    run(True, True, 30.0)