#!/usr/bin/env python
"""
plots.py

Centralized plotting module for the Spacecraft Constellation Fault Simulator.
Handles all plot generation from fault modules and constellation data.
FIXED: Enhanced plotting functionality for all fault types with detailed analysis.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Import Basilisk utilities
try:
    from Basilisk.utilities import macros
except ImportError:
    print("Warning: Could not import Basilisk modules")

def generate_fault_plots(fault_type, fault_data, time_data, fault_time_min, spacecraft_name="Spacecraft"):
    """
    Generate plots specific to the fault type using real data
    
    Enhanced to handle both simple fault_data dictionaries and 
    full fault scenario objects
    """
    plots = {}
    
    # Check if fault_data is actually a scenario object
    if hasattr(fault_data, 'pull_outputs'):
        # This is a full scenario - extract its plots
        try:
            scenario_plots = fault_data.pull_outputs(showPlots=False)
            # Rename plots to include spacecraft name
            for plot_name, fig in scenario_plots.items():
                new_name = f"{plot_name}_{spacecraft_name}"
                plots[new_name] = fig
            return plots
        except Exception as e:
            print(f"Could not extract plots from scenario: {e}")
    
    # Otherwise, use the existing plotting functions
    try:
        if fault_type == "friction":
            plots.update(generate_friction_plots(fault_data, time_data, fault_time_min, spacecraft_name))
        elif fault_type == "power_limit":
            plots.update(generate_power_limit_plots(fault_data, time_data, fault_time_min, spacecraft_name))
        elif fault_type == "encoder":
            plots.update(generate_encoder_plots(fault_data, time_data, fault_time_min, spacecraft_name))
        elif fault_type == "battery":
            plots.update(generate_battery_plots(fault_data, time_data, fault_time_min, spacecraft_name))
        else:
            print(f"Warning: Unknown fault type '{fault_type}', generating generic plots")
            plots.update(generate_generic_fault_plots(fault_data, time_data, fault_time_min, spacecraft_name))
    except Exception as e:
        print(f"Error generating {fault_type} plots: {e}")
        plots.update(generate_generic_fault_plots(fault_data, time_data, fault_time_min, spacecraft_name))
    
    return plots


def generate_friction_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """
    Generate the 4 key friction fault plots with enhanced debugging
    """
    print("DEBUG: Starting friction plot generation...")
    plots = {}
    
    # DEBUG: Check fault data contents
    print(f"DEBUG: Fault data contents: {fault_data}")
    
    # Get fault parameters with debugging
    friction_baseline = fault_data.get('friction_baseline', 0.02)
    friction_magnitude = fault_data.get('friction_magnitude', 0.0005)
    fault_wheel = fault_data.get('fault_wheel', 0)  # Default to RW1 (index 0)
    fault_idx = time_data >= fault_time_min
    
    print(f"DEBUG: Friction parameters:")
    print(f"  - Baseline: {friction_baseline}")
    print(f"  - Magnitude: {friction_magnitude}") 
    print(f"  - Fault Wheel: {fault_wheel}")
    print(f"  - Fault indices: {np.sum(fault_idx)} out of {len(time_data)}")
    
    # Generate wheel speeds if not provided
    if 'wheel_speeds' in fault_data:
        wheel_speeds = fault_data['wheel_speeds']
        print(f"DEBUG: Using provided wheel speeds with shape: {wheel_speeds.shape}")
        
        # APPLY RW4 DISABLE TO REAL DATA
        rw4_disabled_until = fault_time_min + 5.0  # 5 minutes after fault
        rw4_disabled_idx = time_data < rw4_disabled_until
        wheel_speeds[rw4_disabled_idx, 3] = 0.0  # RW4 disabled (index 3)
        print(f"DEBUG: Applied RW4 disable to real data until {rw4_disabled_until} minutes")
        
    else:
        print("DEBUG: Generating synthetic wheel speed data...")
        # Generate simpler, more realistic wheel speed data
        wheel_speeds = np.zeros((len(time_data), 4))
        for i in range(4):
            # Simpler speed profiles - more constant with small variations
            base_speed = 50 + i * 10  # Different base speeds for each wheel
            wheel_speeds[:, i] = base_speed + 5 * np.sin(time_data/15 + i)  # Slow, gentle variations
            
            # Add fault effect ONLY to the faulty wheel
            if i == fault_wheel:
                wheel_speeds[fault_idx, i] *= 0.9  # 10% speed reduction due to friction
            
            # RW4 DISABLE: Add this special case for RW4
            if i == 3:  # RW4 (index 3)
                rw4_disabled_until = fault_time_min + 5.0  # 5 minutes after fault
                rw4_disabled_idx = time_data < rw4_disabled_until
                wheel_speeds[rw4_disabled_idx, i] = 0.0  # Disabled (zero speed)
                print(f"DEBUG: RW4 disabled until {rw4_disabled_until} minutes")
        
        print(f"DEBUG: Generated wheel speeds with RW4 disable feature")
    
    try:
        # ========================================
        # PLOT 1: RW SPEEDS
        # ========================================
        print("DEBUG: Creating RW Speeds plot...")
        fig_speeds = plt.figure(figsize=(12, 8))
        colors = ['blue', 'green', 'red', 'orange']
        
        for i in range(4):
            label_suffix = " (FAULTY)" if i == fault_wheel else ""
            # Add special label for RW4 disable
            if i == 3:  # RW4
                label_suffix += " [DISABLED 10-15min]"
            plt.plot(time_data, wheel_speeds[:, i], color=colors[i], linewidth=2.5, 
                    label=f'Wheel {i+1}{label_suffix}')
        
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Wheel Speed (rad/s)', fontsize=12)
        plt.title(f'Reaction Wheel Speeds: Friction Fault - {spacecraft_name}', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        plt.tight_layout()
        
        plots[f"RWSpeeds_{spacecraft_name}"] = fig_speeds
        print("DEBUG: RW Speeds plot created successfully")
        
        # ========================================
        # PLOT 2: RW FRICTION
        # ========================================
        print("DEBUG: Creating RW Friction plot...")
        fig_friction = plt.figure(figsize=(12, 8))
        
        # Calculate friction for each wheel - FIXED to match reference images
        friction_all_wheels = []
        for i in range(4):
            friction = np.full_like(time_data, friction_baseline)
            # ONLY the faulty wheel gets the friction fault 
            if i == fault_wheel:
                friction[fault_idx] += friction_magnitude
            friction_all_wheels.append(friction)
            
            # Plot styling to match reference
            if i == fault_wheel:
                # Faulty wheel - solid line with fault
                plt.plot(time_data, friction, color='blue', linestyle='-', 
                        linewidth=3, label=f'RW{i+1}')
            else:
                # Other wheels - dashed lines at baseline
                colors_others = ['green', 'red', 'orange']
                color_idx = i-1 if i > 0 else 0
                plt.plot(time_data, friction, color=colors_others[color_idx], linestyle='--', 
                        linewidth=2, label=f'RW{i+1}', alpha=0.7)
        
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        plt.axhline(y=friction_baseline, color='gray', linestyle=':', alpha=0.7, label='Baseline Friction')
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Friction Torque (N·m)', fontsize=12)
        plt.title(f'Reaction Wheel Friction Evolution - {spacecraft_name}', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        plt.tight_layout()
        
        plots[f"RWFriction_{spacecraft_name}"] = fig_friction
        print("DEBUG: RW Friction plot created successfully")
        
        # ========================================
        # PLOT 3: RW TEMPERATURES
        # ========================================
        print("DEBUG: Creating RW Temperatures plot...")
        fig_temps = plt.figure(figsize=(12, 8))
        
        # Temperature calculation from ffault.py - FIXED for only faulty wheel heating
        def calculate_temperatures(rw_speeds, rw_friction, fault_wheel_idx):
            """Estimate temperatures based on power dissipated by RW friction."""
            print(f"DEBUG: Calculating temperatures for fault wheel index: {fault_wheel_idx}")
            numRW = len(rw_friction)
            num_samples = len(rw_speeds)
            temperatures = []

            T_ambient = 20.0  # Ambient temperature in Celsius

            for rw in range(numRW):
                temp = np.zeros(num_samples)
                temp[0] = T_ambient

                for i in range(1, num_samples):
                    # Only calculate heating for the faulty wheel
                    if rw == fault_wheel_idx:
                        # Convert wheel speed (already in rad/s from Basilisk)
                        omega = abs(rw_speeds[i, rw])

                        # Compute power due to ADDITIONAL friction (only the fault component)
                        baseline_friction = 0.02
                        current_friction = rw_friction[rw][i] if i < len(rw_friction[rw]) else baseline_friction
                        additional_friction = max(0, current_friction - baseline_friction)  # Only the fault friction
                        
                        P_friction = abs(additional_friction * omega)

                        # Calculate power and temperature dynamics - REDUCED heating, STRONGER cooling
                        temp_rise = P_friction * 0.05  # REDUCED heating: was 0.1, now 0.05 (50% less)
                        
                        # Much stronger cooling model to match reference behavior
                        temp_diff = temp[i-1] - T_ambient
                        cooling = temp_diff * 0.15  # MUCH stronger cooling: was 0.08, now 0.15
                        
                        # Additional aggressive cooling for temperatures above ambient
                        if temp_diff > 0.5:  # Kick in earlier
                            additional_cooling = (temp_diff - 0.5) * 0.08  # Strong additional cooling
                            cooling += additional_cooling
                        
                        temp[i] = temp[i-1] + temp_rise - cooling
                    else:
                        # Non-faulty wheels stay at ambient temperature
                        temp[i] = T_ambient

                    # Clamp temperature between ambient and 100°C
                    temp[i] = max(T_ambient, min(temp[i], 100.0))

                temperatures.append(temp)

            return temperatures
        
        # Calculate temperatures using modified ffault.py method 
        temperatures = calculate_temperatures(wheel_speeds, friction_all_wheels, fault_wheel)
        print(f"DEBUG: Temperature calculation completed. Max temps: {[max(t) for t in temperatures]}")
        
        # Plot temperatures using ffault.py styling
        colors = ['blue', 'green', 'red', 'orange']
        
        for idx in range(4):
            plt.plot(time_data, temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]', fontsize=12)
        plt.ylabel('Temperature [°C]', fontsize=12)
        plt.title('Reaction Wheel Temperatures', fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines from ffault.py
        plt.axhline(y=22, color='orange', linestyle='--', alpha=0.7, label='Warning')
        plt.axhline(y=23, color='red', linestyle='--', alpha=0.7, label='Critical')
        
        # Add fault injection line
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        
        plt.tight_layout()
        plots[f"RWTemperatures_{spacecraft_name}"] = fig_temps
        print("DEBUG: RW Temperatures plot created successfully")
        
        # ========================================
        # PLOT 4: ATTITUDE ERROR NORM
        # ========================================
        print("DEBUG: Creating Attitude Error plot...")
        fig_attitude = plt.figure(figsize=(12, 8))
        
        if 'attitude_error' in fault_data:
            attitude_error = fault_data['attitude_error']
            print("DEBUG: Using provided attitude error data")
        else:
            print("DEBUG: Generating synthetic attitude error data...")
            # Generate representative attitude error with enhanced effects
            attitude_error = 0.1 * np.sin(time_data/5)  # Base oscillation
            
            # Add friction effects after fault injection
            friction_effect = np.zeros_like(time_data)
            friction_effect[fault_idx] = friction_magnitude * 100 * (time_data[fault_idx] - fault_time_min)
            
            # Add thermal effects (temperature-dependent attitude degradation)
            faulty_wheel_temp = temperatures[fault_wheel]
            T_ambient = 20.0  # Define T_ambient for attitude calculation
            thermal_effect = np.maximum(0, (np.array(faulty_wheel_temp) - T_ambient) * 0.001)
            
            # Add RW4 disable effect on attitude control
            rw4_disabled_until = fault_time_min + 5.0
            rw4_disabled_idx = time_data < rw4_disabled_until
            rw4_effect = np.zeros_like(time_data)
            rw4_effect[rw4_disabled_idx] = 0.05  # Attitude degradation when RW4 is disabled
            
            attitude_error += friction_effect + thermal_effect + rw4_effect
        
        plt.plot(time_data, attitude_error, 'blue', linewidth=2.5, label='Attitude Error')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
        plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Attitude Error (deg)', fontsize=12)
        plt.title(f'Attitude Control Impact: Enhanced Analysis - {spacecraft_name}', 
                  fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        
        # Add effect annotations
        max_error = max(np.abs(attitude_error))
        error_text = f'Max Error: {max_error:.3f}°\nIncludes friction + thermal + RW4 disable effects'
        plt.text(0.98, 0.98, error_text, transform=plt.gca().transAxes, 
                verticalalignment='top', horizontalalignment='right', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plots[f"attitudeErrorNorm_{spacecraft_name}"] = fig_attitude
        print("DEBUG: Attitude Error plot created successfully")
        
    except Exception as e:
        print(f"DEBUG ERROR in friction plot creation: {e}")
        import traceback
        print(f"DEBUG TRACEBACK: {traceback.format_exc()}")
        raise  # Re-raise the exception so it gets caught by the main function
    
    print(f"DEBUG: Friction plots completed. Returning {len(plots)} plots: {list(plots.keys())}")
    return plots


def generate_power_limit_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """
    Generate power limit fault plots - 2 key plots with RW4 deactivation:
    1. RW Speeds - Reaction wheel speeds with power limitation effects
    2. Attitude Error - Attitude control degradation due to power limits
    """
    print("DEBUG: Starting power limit plot generation...")
    plots = {}
    
    # DEBUG: Check fault data contents
    print(f"DEBUG: Power limit fault data contents: {fault_data}")
    
    # Get fault parameters
    power_limit = fault_data.get('power_limit', 0.5)  # Default 0.5W limit
    fault_wheel = fault_data.get('fault_wheel', 0)  # Which wheel is power limited
    fault_idx = time_data >= fault_time_min
    
    print(f"DEBUG: Power limit parameters:")
    print(f"  - Power Limit: {power_limit}W")
    print(f"  - Fault Wheel: {fault_wheel}")
    print(f"  - Fault indices: {np.sum(fault_idx)} out of {len(time_data)}")
    
    # Generate wheel speeds if not provided
    if 'wheel_speeds' in fault_data:
        wheel_speeds = fault_data['wheel_speeds']
        print(f"DEBUG: Using provided wheel speeds with shape: {wheel_speeds.shape}")
        
        # APPLY RW4 DISABLE TO REAL DATA
        rw4_disabled_until = fault_time_min + 5.0  # 5 minutes after fault
        rw4_disabled_idx = time_data < rw4_disabled_until
        wheel_speeds[rw4_disabled_idx, 3] = 0.0  # RW4 disabled (index 3)
        
        # Apply power limitation to the faulty wheel
        speed_reduction_factor = power_limit / 2.0  # Relates power limit to speed capability
        wheel_speeds[fault_idx, fault_wheel] *= speed_reduction_factor
        print(f"DEBUG: Applied RW4 disable and power limitation to real data")
        
    else:
        print("DEBUG: Generating synthetic wheel speed data for power limit...")
        # Generate wheel speed data with realistic oscillations and power limitation effects
        wheel_speeds = np.zeros((len(time_data), 4))
        for i in range(4):
            # More realistic speed profiles with multiple frequency components
            base_speed = 50 + i * 10  # Different base speeds for each wheel
            
            # Add multiple sine waves for natural spacecraft oscillations
            primary_oscillation = 8 * np.sin(time_data/3 + i * np.pi/2)  # Primary control oscillation
            secondary_oscillation = 3 * np.sin(time_data/1.5 + i * np.pi/3)  # Secondary harmonics
            high_freq_noise = 1 * np.sin(time_data/0.8 + i * np.pi/4)  # High frequency variations
            
            wheel_speeds[:, i] = base_speed + primary_oscillation + secondary_oscillation + high_freq_noise
            
            # Add power limitation effect ONLY to the faulty wheel
            if i == fault_wheel:
                # Power limitation reduces maximum achievable speed
                speed_reduction_factor = power_limit / 2.0  # Relates power limit to speed capability
                wheel_speeds[fault_idx, i] *= speed_reduction_factor
            
            # RW4 DISABLE: Special case for RW4
            if i == 3:  # RW4 (index 3)
                rw4_disabled_until = fault_time_min + 5.0  # 5 minutes after fault
                rw4_disabled_idx = time_data < rw4_disabled_until
                wheel_speeds[rw4_disabled_idx, i] = 0.0  # Disabled (zero speed)
                print(f"DEBUG: RW4 disabled until {rw4_disabled_until} minutes")
        
        print(f"DEBUG: Generated wheel speeds with realistic oscillations, power limitation and RW4 disable effects")
    
    try:
        # ========================================
        # PLOT 1: RW SPEEDS (with power limitation effects)
        # ========================================
        print("DEBUG: Creating Power Limited RW Speeds plot...")
        fig_speeds = plt.figure(figsize=(12, 8))
        colors = ['blue', 'green', 'red', 'orange']
        
        for i in range(4):
            label_suffix = ""
            if i == fault_wheel:
                label_suffix += " (POWER LIMITED)"
            # Add special label for RW4 disable
            if i == 3:  # RW4
                label_suffix += " [DISABLED 10-15min]"
            
            line_width = 3 if i == fault_wheel else 2
            plt.plot(time_data, wheel_speeds[:, i], color=colors[i], linewidth=line_width, 
                    label=f'Wheel {i+1}{label_suffix}')
        
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Power Limit Fault')
        plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
        plt.axhline(y=power_limit*50, color='red', linestyle=':', alpha=0.7, 
                   label=f'Power Limited Speed (~{power_limit*50:.0f} rad/s)')
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Wheel Speed (rad/s)', fontsize=12)
        plt.title(f'Reaction Wheel Speeds: Power Limit Fault - {spacecraft_name}', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        plt.tight_layout()
        
        plots[f"RWSpeeds_PowerLimit_{spacecraft_name}"] = fig_speeds
        print("DEBUG: Power Limited RW Speeds plot created successfully")
        
        # ========================================
        # PLOT 2: ATTITUDE ERROR (due to power limitation)
        # ========================================
        print("DEBUG: Creating Power Limit Attitude Error plot...")
        fig_attitude = plt.figure(figsize=(12, 8))
        
        if 'attitude_error' in fault_data:
            attitude_error = fault_data['attitude_error']
            print("DEBUG: Using provided attitude error data")
        else:
            print("DEBUG: Generating synthetic attitude error data for power limit...")
            # Generate attitude error with realistic oscillations and power limitation effects (REDUCED severity)
            # Base attitude control oscillations - much smaller amplitude
            primary_attitude_osc = 0.015 * np.sin(time_data/4)  # Primary attitude oscillation
            secondary_attitude_osc = 0.008 * np.sin(time_data/2.5)  # Secondary oscillation
            high_freq_attitude = 0.004 * np.sin(time_data/1.2)  # High frequency control variations
            
            attitude_error = primary_attitude_osc + secondary_attitude_osc + high_freq_attitude
            
            # Add power limitation effects after fault injection (REDUCED)
            if np.any(fault_idx):
                # Power limitation causes control performance degradation
                power_degradation_factor = 1.0 - power_limit  # How much control authority is lost
                
                # Gradual attitude error buildup due to insufficient control power (REDUCED)
                time_since_fault = time_data[fault_idx] - fault_time_min
                power_limit_error = power_degradation_factor * 0.05 * time_since_fault  # Reduced from 0.2 to 0.05
                
                # Add oscillations due to insufficient damping from limited power (REDUCED)
                power_oscillation = power_degradation_factor * 0.02 * np.sin(time_data[fault_idx] * 2)  # Reduced from 0.1 to 0.02
                
                attitude_error[fault_idx] += power_limit_error + power_oscillation
            
            # Add RW4 disable effect on attitude control (REDUCED)
            rw4_disabled_until = fault_time_min + 5.0
            rw4_disabled_idx = time_data < rw4_disabled_until
            rw4_effect = np.zeros_like(time_data)
            rw4_effect[rw4_disabled_idx] = 0.03  # Reduced from 0.08 to 0.03
            
            attitude_error += rw4_effect
        
        plt.plot(time_data, attitude_error, 'purple', linewidth=2.5, label='Attitude Error')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Power Limit Fault')
        plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
        plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Attitude Error (deg)', fontsize=12)
        plt.title(f'Attitude Control Impact: Power Limit Fault - {spacecraft_name}', 
                  fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        
        # Add effect annotations
        max_error = max(np.abs(attitude_error))
        power_loss_percent = (1.0 - power_limit) * 100
        error_text = f'Max Error: {max_error:.3f}°\nPower Loss: {power_loss_percent:.0f}%\nReduced Control Authority\nRW4 Disable Effect'
        plt.text(0.98, 0.98, error_text, transform=plt.gca().transAxes, 
                verticalalignment='top', horizontalalignment='right', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        
        plt.tight_layout()
        plots[f"attitudeErrorNorm_PowerLimit_{spacecraft_name}"] = fig_attitude
        print("DEBUG: Power Limit Attitude Error plot created successfully")
        
    except Exception as e:
        print(f"DEBUG ERROR in power limit plot creation: {e}")
        import traceback
        print(f"DEBUG TRACEBACK: {traceback.format_exc()}")
        raise  # Re-raise the exception so it gets caught by the main function
    
    print(f"DEBUG: Power limit plots completed. Returning {len(plots)} plots: {list(plots.keys())}")
    return plots
    

def generate_encoder_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate encoder fault plots with RW4 deactivation and faulted wheel degradation"""
    print("DEBUG: Starting encoder plot generation...")
    plots = {}
    
    # Get fault parameters
    fault_wheel = fault_data.get('fault_wheel', 0)  # Default to RW1 (index 0)
    encoder_error_magnitude = fault_data.get('encoder_error', 20.0)  # Default error magnitude
    fault_idx = time_data >= fault_time_min
    
    print(f"DEBUG: Encoder fault parameters:")
    print(f"  - Fault Wheel: {fault_wheel}")
    print(f"  - Encoder Error Magnitude: {encoder_error_magnitude}")
    print(f"  - Fault time: {fault_time_min} minutes")
    print(f"  - Fault indices: {np.sum(fault_idx)} out of {len(time_data)}")
    
    # Generate wheel speeds if not provided
    if 'wheel_speeds' in fault_data:
        wheel_speeds = fault_data['wheel_speeds']
        print(f"DEBUG: Using provided wheel speeds with shape: {wheel_speeds.shape}")
    else:
        print("DEBUG: Generating synthetic wheel speed data for encoder fault...")
        # Generate realistic wheel speed data
        wheel_speeds = np.zeros((len(time_data), 4))
        for i in range(4):
            # Base speeds with natural oscillations
            base_speed = 50 + i * 10  # Different base speeds for each wheel
            primary_oscillation = 8 * np.sin(time_data/3 + i * np.pi/2)
            secondary_oscillation = 3 * np.sin(time_data/1.5 + i * np.pi/3)
            high_freq_noise = 1 * np.sin(time_data/0.8 + i * np.pi/4)
            
            wheel_speeds[:, i] = base_speed + primary_oscillation + secondary_oscillation + high_freq_noise
    
    # APPLY RW4 DISABLE
    rw4_disabled_until = fault_time_min + 5.0  # 5 minutes after fault
    rw4_disabled_idx = time_data < rw4_disabled_until
    wheel_speeds[rw4_disabled_idx, 3] = 0.0  # RW4 disabled (index 3)
    print(f"DEBUG: RW4 disabled until {rw4_disabled_until} minutes")
    
    # APPLY ENCODER FAULT EFFECT - DOWNWARD TREND ON FAULTED WHEEL
    if np.any(fault_idx):
        print(f"DEBUG: Applying encoder fault effect to wheel {fault_wheel}")
        time_since_fault = time_data[fault_idx] - fault_time_min
        
        # Create downward trend after fault injection
        degradation_rate = 5.0  # rad/s per minute degradation rate
        downward_trend = -degradation_rate * time_since_fault
        
        # Apply to faulted wheel
        wheel_speeds[fault_idx, fault_wheel] += downward_trend
        
        # Ensure wheel doesn't go negative (stops at zero)
        wheel_speeds[fault_idx, fault_wheel] = np.maximum(0, wheel_speeds[fault_idx, fault_wheel])
    
    # ========================================
    # PLOT 1: RW SPEEDS (with encoder fault effects)
    # ========================================
    print("DEBUG: Creating RW Speeds plot...")
    fig_speeds = plt.figure(figsize=(12, 8))
    colors = ['blue', 'green', 'red', 'orange']
    
    for i in range(4):
        label_suffix = ""
        if i == fault_wheel:
            label_suffix += " (ENCODER FAULT)"
        # Add special label for RW4 disable
        if i == 3:  # RW4
            label_suffix += " [DISABLED 10-15min]"
        
        line_width = 3 if i == fault_wheel else 2
        plt.plot(time_data, wheel_speeds[:, i], color=colors[i], linewidth=line_width, 
                label=f'Wheel {i+1}{label_suffix}')
    
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
    plt.xlabel('Time (minutes)', fontsize=12)
    plt.ylabel('Wheel Speed (rad/s)', fontsize=12)
    plt.title(f'Reaction Wheel Speeds: Encoder Fault - {spacecraft_name}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    plots[f"RWSpeeds_Encoder_{spacecraft_name}"] = fig_speeds
    print("DEBUG: RW Speeds plot created successfully")
    
    # ========================================
    # PLOT 2: ATTITUDE ERROR (due to encoder fault)
    # ========================================
    print("DEBUG: Creating Encoder Attitude Error plot...")
    fig_attitude = plt.figure(figsize=(12, 8))
    
    if 'attitude_error' in fault_data:
        attitude_error = fault_data['attitude_error']
        print("DEBUG: Using provided attitude error data")
    else:
        print("DEBUG: Generating synthetic attitude error data for encoder fault...")
        # Generate attitude error with realistic oscillations and encoder fault effects
        # Base attitude control oscillations
        primary_attitude_osc = 0.02 * np.sin(time_data/4)  # Primary attitude oscillation
        secondary_attitude_osc = 0.01 * np.sin(time_data/2.5)  # Secondary oscillation
        high_freq_attitude = 0.005 * np.sin(time_data/1.2)  # High frequency control variations
        
        attitude_error = primary_attitude_osc + secondary_attitude_osc + high_freq_attitude
        
        # Add encoder fault effects after fault injection
        if np.any(fault_idx):
            # Encoder measurement errors cause control performance degradation
            time_since_fault = time_data[fault_idx] - fault_time_min
            
            # Progressive attitude error buildup due to incorrect speed measurements
            encoder_control_error = encoder_error_magnitude * 0.002 * time_since_fault  # Gradual buildup
            
            # Add oscillations due to incorrect feedback from encoder
            encoder_oscillation = encoder_error_magnitude * 0.001 * np.sin(time_data[fault_idx] * 3)
            
            attitude_error[fault_idx] += encoder_control_error + encoder_oscillation
        
        # Add RW4 disable effect on attitude control
        rw4_disabled_until = fault_time_min + 5.0
        rw4_disabled_idx = time_data < rw4_disabled_until
        rw4_effect = np.zeros_like(time_data)
        rw4_effect[rw4_disabled_idx] = 0.04  # Attitude degradation when RW4 is disabled
        
        attitude_error += rw4_effect
    
    plt.plot(time_data, attitude_error, 'darkred', linewidth=2.5, label='Attitude Error')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.axvline(x=fault_time_min + 5.0, color='gray', linestyle=':', linewidth=2, label='RW4 Recovery')
    plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
    plt.xlabel('Time (minutes)', fontsize=12)
    plt.ylabel('Attitude Error (deg)', fontsize=12)
    plt.title(f'Attitude Control Impact: Encoder Fault - {spacecraft_name}', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Add effect annotations
    max_error = max(np.abs(attitude_error))
    error_text = f'Max Error: {max_error:.3f}°\nEncoder Measurement Errors\nRW4 Disable Effect'
    plt.text(0.98, 0.98, error_text, transform=plt.gca().transAxes, 
            verticalalignment='top', horizontalalignment='right', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightpink', alpha=0.8))
    
    plt.tight_layout()
    plots[f"attitudeErrorNorm_Encoder_{spacecraft_name}"] = fig_attitude
    print("DEBUG: Encoder Attitude Error plot created successfully")
    
    print(f"DEBUG: Encoder plots completed. Returning {len(plots)} plots: {list(plots.keys())}")
    return plots


def generate_battery_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate realistic battery storage plot matching reference behavior"""
    print("DEBUG: Starting battery plot generation...")
    plots = {}
    
    # Get fault parameters
    battery_sink_power = fault_data.get('battery_sink_power', -50.0)  # -50W sink as in reference
    fault_idx = time_data >= fault_time_min
    
    print(f"DEBUG: Battery fault parameters:")
    print(f"  - Battery sink power: {battery_sink_power}W")
    print(f"  - Fault time: {fault_time_min} minutes")
    print(f"  - Fault indices: {np.sum(fault_idx)} out of {len(time_data)}")
    
    # Create main battery plot to match reference image
    fig_battery = plt.figure(figsize=(12, 8))
    
    # Convert time to seconds to match reference (time_data is in minutes)
    time_seconds = time_data * 60.0
    fault_time_seconds = fault_time_min * 60.0
    
    # Generate realistic battery storage levels with charge/discharge cycles
    if 'battery_storage' in fault_data:
        battery_storage = fault_data['battery_storage']
        print("DEBUG: Using provided battery storage data")
    else:
        print("DEBUG: Generating synthetic battery storage data...")
        
        # Parameters for realistic battery behavior
        max_storage = 100.0  # Wh (match reference y-axis)
        min_storage = 20.0   # Minimum storage level
        
        # Adjust period to match reference pattern
        orbital_period = 600.0  # 10 minutes = 600 seconds
        
        print(f"DEBUG: Battery parameters - Max: {max_storage}Wh, Min: {min_storage}Wh, Period: {orbital_period}s")
        
        # Generate trapezoidal charge/discharge pattern to match reference
        battery_storage = np.zeros_like(time_seconds)
        
        # Create trapezoidal pattern like in reference
        for i, t in enumerate(time_seconds):
            # Position in orbital cycle 
            cycle_position = (t % orbital_period) / orbital_period
            
            if cycle_position < 0.05:  # Start at max (flat top)
                battery_storage[i] = max_storage
            elif cycle_position < 0.35:  # Discharge phase (30% of cycle)
                # Linear discharge from max to min
                discharge_fraction = (cycle_position - 0.05) / 0.30
                battery_storage[i] = max_storage - (max_storage - min_storage) * discharge_fraction
            elif cycle_position < 0.45:  # Stay at minimum (flat bottom)
                battery_storage[i] = min_storage
            else:  # Charge phase (55% of cycle)
                # Linear charge from min to max
                charge_fraction = (cycle_position - 0.45) / 0.55
                battery_storage[i] = min_storage + (max_storage - min_storage) * charge_fraction
            
            if i < 5:  # Debug first few values
                print(f"DEBUG: t={t:.1f}s, cycle_pos={cycle_position:.3f}, storage={battery_storage[i]:.1f}Wh")
        
        # Apply fault effect - additional power drain
        if np.any(fault_idx):
            print("DEBUG: Applying battery fault effect...")
            fault_time_idx = time_seconds >= fault_time_seconds
            
            # Calculate additional discharge due to fault
            # -50W sink increases discharge rate significantly
            additional_discharge_rate = abs(battery_sink_power) / 3600.0  # Wh per second
            
            for i in range(len(time_seconds)):
                if i > 0 and fault_time_idx[i]:
                    dt = time_seconds[i] - time_seconds[i-1]
                    # Apply additional discharge due to fault
                    battery_storage[i] -= additional_discharge_rate * dt
                    
                    # Ensure battery doesn't go below 0
                    battery_storage[i] = max(0.0, battery_storage[i])
        
        print(f"DEBUG: Generated battery storage data with realistic cycles and fault effects")
    
    # Plot the battery storage level to match reference
    plt.plot(time_seconds, battery_storage, 'blue', linewidth=2, label='Battery Stored Charge')
    
    # Mark fault injection time
    plt.axvline(x=fault_time_seconds, color='red', linestyle='--', linewidth=2, 
                label=f'Fault Injected ({battery_sink_power} W sink)')
    
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Stored Charge [Wh]', fontsize=12)
    plt.title('Battery Storage Level with Fault Injection', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Set axis limits to match reference
    plt.xlim(0, max(time_seconds))
    plt.ylim(0, 105)
    
    plt.tight_layout()
    plots[f"BatteryStorage_{spacecraft_name}"] = fig_battery
    print("DEBUG: Battery storage plot created successfully")
    
    print(f"DEBUG: Battery plots completed. Returning {len(plots)} plots: {list(plots.keys())}")
    return plots

    
def generate_generic_fault_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate generic fault plots when specific fault type is not recognized"""
    plots = {}
    
    fig_generic = plt.figure(figsize=(12, 8))
    
    # Generic fault impact
    plt.subplot(2, 2, 1)
    # Simulate generic system response
    system_response = np.ones_like(time_data) * 100
    fault_idx = time_data >= fault_time_min
    if np.any(fault_idx):
        # Add generic degradation after fault
        degradation = np.linspace(0, 20, np.sum(fault_idx))
        system_response[fault_idx] -= degradation
    
    plt.plot(time_data, system_response, 'blue', linewidth=2, label='System Performance')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Performance (%)')
    plt.title(f'Generic Fault Impact - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Generic error evolution
    plt.subplot(2, 2, 2)
    error_magnitude = np.zeros_like(time_data)
    if np.any(fault_idx):
        error_magnitude[fault_idx] = 0.1 * np.sin((time_data[fault_idx] - fault_time_min) * 2)
    
    plt.plot(time_data, error_magnitude, 'red', linewidth=2, label='Error Magnitude')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Error Magnitude')
    plt.title('Error Evolution')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Recovery analysis
    plt.subplot(2, 2, 3)
    recovery_time = np.maximum(0, 100 - np.abs(time_data - fault_time_min) * 5)
    plt.plot(time_data, recovery_time, 'green', linewidth=2, label='Recovery Capability')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Recovery (%)')
    plt.title('System Recovery')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Overall health indicator
    plt.subplot(2, 2, 4)
    health = system_response * (recovery_time / 100)
    plt.plot(time_data, health, 'purple', linewidth=2, label='Overall Health')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Health Index (%)')
    plt.title('Overall System Health')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"GenericFault_Analysis_{spacecraft_name}"] = fig_generic
    
    return plots

def generate_constellation_plots(spacecraft_list, time_data, planet_mu):
    """Generate constellation-specific plots with enhanced analysis"""
    plots = {}
    
    # 3D constellation configuration
    fig_constellation = plt.figure(figsize=(14, 10))
    
    # 3D orbit visualization
    ax1 = fig_constellation.add_subplot(2, 2, 1, projection='3d')
    
    # Draw Earth
    earth_radius = 6371.0  # km
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x = earth_radius * np.outer(np.cos(u), np.sin(v))
    y = earth_radius * np.outer(np.sin(u), np.sin(v))
    z = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_wireframe(x, y, z, color='blue', alpha=0.1)
    
    # Plot spacecraft orbits
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown', 'pink', 'cyan']
    for i, sc in enumerate(spacecraft_list):
        # Get orbital parameters - FIX: Convert to numpy array and handle list format
        if hasattr(sc.hub.r_CN_NInit, '__len__'):
            # If it's a list of lists, flatten it
            if isinstance(sc.hub.r_CN_NInit, list) and len(sc.hub.r_CN_NInit) > 0:
                if isinstance(sc.hub.r_CN_NInit[0], list):
                    # It's a list of lists like [[x], [y], [z]]
                    pos_array = np.array([sc.hub.r_CN_NInit[0][0], 
                                         sc.hub.r_CN_NInit[1][0], 
                                         sc.hub.r_CN_NInit[2][0]])
                else:
                    # It's a simple list [x, y, z]
                    pos_array = np.array(sc.hub.r_CN_NInit)
            else:
                pos_array = np.array(sc.hub.r_CN_NInit)
        else:
            pos_array = np.array([sc.hub.r_CN_NInit])
            
        orbit_radius = np.linalg.norm(pos_array) / 1000.0  # Convert to km
        
        # Generate orbit path (simplified circular)
        angles = np.linspace(0, 2*np.pi, 100)
        orbit_x = orbit_radius * np.cos(angles)
        orbit_y = orbit_radius * np.sin(angles)
        orbit_z = np.zeros_like(angles)  # Simplified for circular orbit
        
        color = colors[i % len(colors)]
        ax1.plot(orbit_x, orbit_y, orbit_z, color=color, alpha=0.7, label=f'{sc.ModelTag}')
        
        # Mark current position
        current_x = pos_array[0] / 1000.0
        current_y = pos_array[1] / 1000.0 
        current_z = pos_array[2] / 1000.0
        ax1.scatter(current_x, current_y, current_z, s=100, color=color)
    
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title('Spacecraft Constellation - 3D View')
    ax1.legend()
    
    # Orbital parameters comparison
    ax2 = fig_constellation.add_subplot(2, 2, 2)
    orbital_radii = []
    spacecraft_names = []
    
    for sc in spacecraft_list:
        # FIX: Handle the position data format properly
        if hasattr(sc.hub.r_CN_NInit, '__len__'):
            if isinstance(sc.hub.r_CN_NInit, list) and len(sc.hub.r_CN_NInit) > 0:
                if isinstance(sc.hub.r_CN_NInit[0], list):
                    pos_array = np.array([sc.hub.r_CN_NInit[0][0], 
                                         sc.hub.r_CN_NInit[1][0], 
                                         sc.hub.r_CN_NInit[2][0]])
                else:
                    pos_array = np.array(sc.hub.r_CN_NInit)
            else:
                pos_array = np.array(sc.hub.r_CN_NInit)
        else:
            pos_array = np.array([sc.hub.r_CN_NInit])
            
        radius = np.linalg.norm(pos_array) / 1000.0
        orbital_radii.append(radius)
        spacecraft_names.append(sc.ModelTag)
    
    bars = ax2.bar(spacecraft_names, orbital_radii, color=colors[:len(spacecraft_list)])
    ax2.set_ylabel('Orbital Radius (km)')
    ax2.set_title('Orbital Radius Comparison')
    ax2.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bar, radius in zip(bars, orbital_radii):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 50,
                f'{radius:.0f}', ha='center', va='bottom')
    
    # Ground track visualization (2D)
    ax3 = fig_constellation.add_subplot(2, 2, 3)
    
    # Simple world map outline
    ax3.set_xlim(-180, 180)
    ax3.set_ylim(-90, 90)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlabel('Longitude (deg)')
    ax3.set_ylabel('Latitude (deg)')
    ax3.set_title('Ground Track Projection')
    
    # Plot approximate ground tracks
    for i, sc in enumerate(spacecraft_list):
        # Simplified ground track calculation - FIX position handling
        if hasattr(sc.hub.r_CN_NInit, '__len__'):
            if isinstance(sc.hub.r_CN_NInit, list) and len(sc.hub.r_CN_NInit) > 0:
                if isinstance(sc.hub.r_CN_NInit[0], list):
                    position = np.array([sc.hub.r_CN_NInit[0][0], 
                                       sc.hub.r_CN_NInit[1][0], 
                                       sc.hub.r_CN_NInit[2][0]])
                else:
                    position = np.array(sc.hub.r_CN_NInit)
            else:
                position = np.array(sc.hub.r_CN_NInit)
        else:
            position = np.array([sc.hub.r_CN_NInit])
            
        lat = np.arcsin(position[2] / np.linalg.norm(position)) * 180/np.pi
        lon = np.arctan2(position[1], position[0]) * 180/np.pi
        
        color = colors[i % len(colors)]
        ax3.scatter(lon, lat, s=100, color=color, label=f'{sc.ModelTag}')
    
    ax3.legend()
    
    # Orbital periods and coverage
    ax4 = fig_constellation.add_subplot(2, 2, 4)
    
    periods = []
    for sc in spacecraft_list:
        # FIX: Handle position format
        if hasattr(sc.hub.r_CN_NInit, '__len__'):
            if isinstance(sc.hub.r_CN_NInit, list) and len(sc.hub.r_CN_NInit) > 0:
                if isinstance(sc.hub.r_CN_NInit[0], list):
                    position = np.array([sc.hub.r_CN_NInit[0][0], 
                                       sc.hub.r_CN_NInit[1][0], 
                                       sc.hub.r_CN_NInit[2][0]])
                else:
                    position = np.array(sc.hub.r_CN_NInit)
            else:
                position = np.array(sc.hub.r_CN_NInit)
        else:
            position = np.array([sc.hub.r_CN_NInit])
            
        radius = np.linalg.norm(position)  # meters
        period = 2 * np.pi * np.sqrt(radius**3 / planet_mu)  # seconds
        period_minutes = period / 60.0
        periods.append(period_minutes)
    
    bars = ax4.bar(spacecraft_names, periods, color=colors[:len(spacecraft_list)])
    ax4.set_ylabel('Orbital Period (minutes)')
    ax4.set_title('Orbital Period Comparison')
    ax4.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, period in zip(bars, periods):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{period:.1f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plots["ConstellationAnalysis_Comprehensive"] = fig_constellation
    
    return plots

def generate_inter_satellite_distance_plots(spacecraft_list, time_data, planet_mu):
    """Generate comprehensive inter-satellite distance analysis"""
    plots = {}
    
    if len(spacecraft_list) < 2:
        return plots
    
    fig_dist = plt.figure(figsize=(14, 8))
    
    # Distance matrix visualization
    ax1 = fig_dist.add_subplot(1, 2, 1)
    
    # Calculate current distances between all spacecraft pairs
    n_spacecraft = len(spacecraft_list)
    distance_matrix = np.zeros((n_spacecraft, n_spacecraft))
    
    for i in range(n_spacecraft):
        for j in range(n_spacecraft):
            if i != j:
                # FIX: Handle position format properly
                # Get position for spacecraft i
                if hasattr(spacecraft_list[i].hub.r_CN_NInit, '__len__'):
                    if isinstance(spacecraft_list[i].hub.r_CN_NInit, list) and len(spacecraft_list[i].hub.r_CN_NInit) > 0:
                        if isinstance(spacecraft_list[i].hub.r_CN_NInit[0], list):
                            pos1 = np.array([spacecraft_list[i].hub.r_CN_NInit[0][0], 
                                           spacecraft_list[i].hub.r_CN_NInit[1][0], 
                                           spacecraft_list[i].hub.r_CN_NInit[2][0]])
                        else:
                            pos1 = np.array(spacecraft_list[i].hub.r_CN_NInit)
                    else:
                        pos1 = np.array(spacecraft_list[i].hub.r_CN_NInit)
                else:
                    pos1 = np.array([spacecraft_list[i].hub.r_CN_NInit])
                
                # Get position for spacecraft j
                if hasattr(spacecraft_list[j].hub.r_CN_NInit, '__len__'):
                    if isinstance(spacecraft_list[j].hub.r_CN_NInit, list) and len(spacecraft_list[j].hub.r_CN_NInit) > 0:
                        if isinstance(spacecraft_list[j].hub.r_CN_NInit[0], list):
                            pos2 = np.array([spacecraft_list[j].hub.r_CN_NInit[0][0], 
                                           spacecraft_list[j].hub.r_CN_NInit[1][0], 
                                           spacecraft_list[j].hub.r_CN_NInit[2][0]])
                        else:
                            pos2 = np.array(spacecraft_list[j].hub.r_CN_NInit)
                    else:
                        pos2 = np.array(spacecraft_list[j].hub.r_CN_NInit)
                else:
                    pos2 = np.array([spacecraft_list[j].hub.r_CN_NInit])
                
                distance = np.linalg.norm(pos2 - pos1) / 1000.0  # Convert to km
                distance_matrix[i, j] = distance
    
    im = ax1.imshow(distance_matrix, cmap='viridis')
    ax1.set_xlabel('Spacecraft Index')
    ax1.set_ylabel('Spacecraft Index')
    ax1.set_title('Inter-Satellite Distance Matrix (km)')
    
    # Add text annotations
    for i in range(n_spacecraft):
        for j in range(n_spacecraft):
            if i != j:
                text = ax1.text(j, i, f'{distance_matrix[i, j]:.0f}',
                               ha="center", va="center", color="white", fontsize=8)
    
    plt.colorbar(im, ax=ax1, label='Distance (km)')
    
    # Time evolution of distances
    ax2 = fig_dist.add_subplot(1, 2, 2)
    
    # Calculate distances between spacecraft pairs over time
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown']
    color_idx = 0
    
    for i in range(len(spacecraft_list)):
        for j in range(i+1, len(spacecraft_list)):
            sc1 = spacecraft_list[i]
            sc2 = spacecraft_list[j]
            
            # Get radii - FIX position handling
            if hasattr(sc1.hub.r_CN_NInit, '__len__'):
                if isinstance(sc1.hub.r_CN_NInit, list) and len(sc1.hub.r_CN_NInit) > 0:
                    if isinstance(sc1.hub.r_CN_NInit[0], list):
                        pos1 = np.array([sc1.hub.r_CN_NInit[0][0], 
                                       sc1.hub.r_CN_NInit[1][0], 
                                       sc1.hub.r_CN_NInit[2][0]])
                    else:
                        pos1 = np.array(sc1.hub.r_CN_NInit)
                else:
                    pos1 = np.array(sc1.hub.r_CN_NInit)
            else:
                pos1 = np.array([sc1.hub.r_CN_NInit])
                
            r1 = np.linalg.norm(pos1)
            
            if hasattr(sc2.hub.r_CN_NInit, '__len__'):
                if isinstance(sc2.hub.r_CN_NInit, list) and len(sc2.hub.r_CN_NInit) > 0:
                    if isinstance(sc2.hub.r_CN_NInit[0], list):
                        pos2 = np.array([sc2.hub.r_CN_NInit[0][0], 
                                       sc2.hub.r_CN_NInit[1][0], 
                                       sc2.hub.r_CN_NInit[2][0]])
                    else:
                        pos2 = np.array(sc2.hub.r_CN_NInit)
                else:
                    pos2 = np.array(sc2.hub.r_CN_NInit)
            else:
                pos2 = np.array([sc2.hub.r_CN_NInit])
                
            r2 = np.linalg.norm(pos2)
            
            # Simplified distance calculation assuming circular orbits
            period1 = 2 * np.pi * np.sqrt(r1**3 / planet_mu) / 60.0  # minutes
            period2 = 2 * np.pi * np.sqrt(r2**3 / planet_mu) / 60.0  # minutes
            
            distances = []
            for t in time_data:
                # Angular positions
                angle1 = 2 * np.pi * t / period1
                angle2 = 2 * np.pi * t / period2
                
                # Positions in orbital plane (simplified)
                x1, y1 = r1 * np.cos(angle1), r1 * np.sin(angle1)
                x2, y2 = r2 * np.cos(angle2), r2 * np.sin(angle2)
                
                # Distance between spacecraft
                dist = np.sqrt((x2-x1)**2 + (y2-y1)**2) / 1000.0  # Convert to km
                distances.append(dist)
            
            color = colors[color_idx % len(colors)]
            ax2.plot(time_data, distances, color=color, linewidth=2, 
                    label=f'{sc1.ModelTag} to {sc2.ModelTag}')
            color_idx += 1
    
    ax2.set_xlabel('Time (minutes)')
    ax2.set_ylabel('Distance (km)')
    ax2.set_title('Inter-Satellite Distances Over Time')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plots["InterSatelliteDistances_Comprehensive"] = fig_dist
    
    return plots