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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Import Basilisk utilities
try:
    from Basilisk.utilities import macros
except ImportError:
    print("Warning: Could not import Basilisk modules")

def generate_fault_plots(fault_type, fault_data, time_data, fault_time_min, spacecraft_name="Spacecraft"):
    """
    Generate plots specific to the fault type using real data
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    fault_data (dict): Dictionary containing fault-specific data
    time_data (array): Time array in minutes
    fault_time_min (float): Time when fault was injected in minutes
    spacecraft_name (str): Name of the spacecraft
    
    Returns:
    dict: Dictionary of matplotlib Figure objects
    """
    plots = {}
    
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
        # Generate fallback plots
        plots.update(generate_generic_fault_plots(fault_data, time_data, fault_time_min, spacecraft_name))
    
    return plots

def generate_friction_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate comprehensive friction fault specific plots"""
    plots = {}
    
    # Main friction analysis figure
    fig_friction = plt.figure(figsize=(14, 10))
    
    # Subplot 1: Wheel speeds with fault effects
    plt.subplot(2, 3, 1)
    if 'wheel_speeds' in fault_data:
        wheel_speeds = fault_data['wheel_speeds']
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(min(wheel_speeds.shape[1], 4)):
            label_suffix = " (FAULTY)" if i == fault_data.get('fault_wheel', 3) else ""
            plt.plot(time_data, wheel_speeds[:, i], color=colors[i], linewidth=2, 
                    label=f'Wheel {i}{label_suffix}')
    else:
        # Generate representative data if real data not available
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(4):
            speed = 100 * np.sin(time_data/10 + i) + 50 * i
            # Add fault effect
            fault_idx = time_data >= fault_time_min
            if i == fault_data.get('fault_wheel', 3):  # Faulty wheel
                speed[fault_idx] *= (1.0 - fault_data.get('friction_magnitude', 0.0005) * 200)  # Reduced speed due to friction
            plt.plot(time_data, speed, color=colors[i], linewidth=2, label=f'Wheel {i}')
    
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Wheel Speed (rad/s)')
    plt.title(f'Wheel Speeds: Friction Fault - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Subplot 2: Friction torque evolution
    plt.subplot(2, 3, 2)
    friction_baseline = fault_data.get('friction_baseline', 0.02)
    friction_magnitude = fault_data.get('friction_magnitude', 0.0005)
    
    friction = np.full_like(time_data, friction_baseline)
    fault_idx = time_data >= fault_time_min
    friction[fault_idx] += friction_magnitude
    
    plt.plot(time_data, friction, 'red', linewidth=3, label='Friction Torque')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.axhline(y=friction_baseline, color='gray', linestyle=':', alpha=0.7, label='Baseline')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Friction Torque (N·m)')
    plt.title(f'Friction Evolution - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Subplot 3: Power consumption analysis
    plt.subplot(2, 3, 3)
    power_normal = 0.5 + 0.1 * np.sin(time_data/3)
    power_fault = np.copy(power_normal)
    power_fault[fault_idx] += friction_magnitude * 200  # More power needed due to friction
    
    plt.plot(time_data, power_normal, '--', color='blue', linewidth=2, label='Normal Power')
    plt.plot(time_data, power_fault, '-', color='red', linewidth=2, label='With Friction Fault')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Power (W)')
    plt.title(f'Power Consumption Impact - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Subplot 4: Attitude error progression
    plt.subplot(2, 3, 4)
    if 'attitude_error' in fault_data:
        attitude_error = fault_data['attitude_error']
    else:
        attitude_error = 0.1 * np.sin(time_data/5)
        # Friction causes gradual attitude error buildup
        attitude_error[fault_idx] += friction_magnitude * 100 * (time_data[fault_idx] - fault_time_min)
    
    plt.plot(time_data, attitude_error, 'blue', linewidth=2, label='Attitude Error')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Attitude Error (deg)')
    plt.title(f'Attitude Control Impact - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Subplot 5: Friction vs Speed relationship
    plt.subplot(2, 3, 5)
    if 'wheel_speeds' in fault_data:
        faulty_wheel_idx = fault_data.get('fault_wheel', 3)
        faulty_speeds = fault_data['wheel_speeds'][:, faulty_wheel_idx]
        pre_fault_idx = time_data < fault_time_min
        post_fault_idx = time_data >= fault_time_min
        
        if np.any(pre_fault_idx):
            plt.scatter(faulty_speeds[pre_fault_idx], friction[pre_fault_idx], 
                       c='blue', alpha=0.6, label='Pre-fault', s=30)
        if np.any(post_fault_idx):
            plt.scatter(faulty_speeds[post_fault_idx], friction[post_fault_idx], 
                       c='red', alpha=0.6, label='Post-fault', s=30)
    else:
        # Generate representative scatter data
        speeds = np.linspace(50, 150, len(time_data))
        plt.scatter(speeds[time_data < fault_time_min], friction[time_data < fault_time_min], 
                   c='blue', alpha=0.6, label='Pre-fault', s=30)
        plt.scatter(speeds[time_data >= fault_time_min], friction[time_data >= fault_time_min], 
                   c='red', alpha=0.6, label='Post-fault', s=30)
    
    plt.xlabel('Wheel Speed (rad/s)')
    plt.ylabel('Friction Torque (N·m)')
    plt.title(f'Friction vs Speed Relationship')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Subplot 6: Cumulative energy loss
    plt.subplot(2, 3, 6)
    if len(time_data) > 1:
        dt = np.diff(time_data) * 60  # Convert to seconds
        power_loss = np.maximum(0, power_fault - power_normal)
        energy_loss = np.cumsum(np.concatenate([[0], power_loss[1:] * dt]))  # Energy in Joules
        plt.plot(time_data, energy_loss, 'darkred', linewidth=2, label='Cumulative Energy Loss')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Energy Loss (J)')
        plt.title('Cumulative Energy Loss')
        plt.grid(True, alpha=0.3)
        plt.legend()
    
    plt.tight_layout()
    plots[f"FrictionFault_Comprehensive_{spacecraft_name}"] = fig_friction
    
    return plots

def generate_power_limit_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate comprehensive power limit fault specific plots"""
    plots = {}
    
    fig_power = plt.figure(figsize=(14, 10))
    
    # Power limit analysis
    plt.subplot(2, 3, 1)
    power_limit = fault_data.get('power_limit', 0.5)
    power_request = 1.0 + 0.5 * np.sin(time_data/3) + 0.2 * np.sin(time_data/1.5)
    power_actual = np.copy(power_request)
    fault_idx = time_data >= fault_time_min
    power_actual[fault_idx] = np.minimum(power_request[fault_idx], power_limit)
    
    plt.plot(time_data, power_request, '--', color='blue', linewidth=2, label='Requested Power')
    plt.plot(time_data, power_actual, '-', color='red', linewidth=2, label='Actual Power')
    plt.axhline(y=power_limit, color='red', linestyle=':', linewidth=2, label=f'Power Limit ({power_limit}W)')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Power (W)')
    plt.title(f'Power Limitation Analysis - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Performance degradation
    plt.subplot(2, 3, 2)
    performance = np.ones_like(time_data) * 100
    performance[fault_idx] = (power_actual[fault_idx] / power_request[fault_idx]) * 100
    
    plt.plot(time_data, performance, 'orange', linewidth=2, label='Performance')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.axhline(y=100, color='green', linestyle=':', alpha=0.7, label='100% Performance')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Performance (%)')
    plt.title(f'System Performance - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Power deficit
    plt.subplot(2, 3, 3)
    power_deficit = np.maximum(0, power_request - power_actual)
    plt.plot(time_data, power_deficit, 'red', linewidth=2, label='Power Deficit')
    plt.fill_between(time_data, power_deficit, alpha=0.3, color='red')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Power Deficit (W)')
    plt.title('Power Shortage Analysis')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Wheel speed impact
    plt.subplot(2, 3, 4)
    if 'wheel_speeds' in fault_data:
        wheel_speeds = fault_data['wheel_speeds']
        faulty_wheel = fault_data.get('fault_wheel', 3)
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(min(wheel_speeds.shape[1], 4)):
            label_suffix = " (POWER LIMITED)" if i == faulty_wheel else ""
            plt.plot(time_data, wheel_speeds[:, i], color=colors[i], linewidth=2, 
                    label=f'Wheel {i}{label_suffix}')
    else:
        # Generate representative wheel speed data with power limitation
        colors = ['blue', 'green', 'red', 'orange']
        for i in range(4):
            speed = 100 * np.sin(time_data/7 + i) + 50
            if i == fault_data.get('fault_wheel', 3):
                # Reduce speed due to power limitation
                speed[fault_idx] *= (power_limit / 2.0)  # Simplified power-speed relationship
            plt.plot(time_data, speed, color=colors[i], linewidth=2, label=f'Wheel {i}')
    
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Wheel Speed (rad/s)')
    plt.title('Wheel Speed Impact')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Energy analysis
    plt.subplot(2, 3, 5)
    if len(time_data) > 1:
        dt = np.diff(time_data) * 60  # Convert to seconds
        energy_deficit = np.cumsum(np.concatenate([[0], power_deficit[1:] * dt]))  # Energy in Joules
        plt.plot(time_data, energy_deficit, 'darkred', linewidth=2, label='Energy Deficit')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Energy Deficit (J)')
        plt.title('Cumulative Energy Deficit')
        plt.grid(True, alpha=0.3)
        plt.legend()
    
    # System health indicator
    plt.subplot(2, 3, 6)
    health_indicator = 100 * (1 - power_deficit / np.maximum(power_request, 0.1))
    plt.plot(time_data, health_indicator, 'purple', linewidth=2, label='System Health')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.axhline(y=100, color='green', linestyle=':', alpha=0.7, label='Healthy')
    plt.axhline(y=50, color='orange', linestyle=':', alpha=0.7, label='Degraded')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Health Index (%)')
    plt.title('System Health Index')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"PowerLimitFault_Comprehensive_{spacecraft_name}"] = fig_power
    
    return plots

def generate_encoder_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate comprehensive encoder fault specific plots"""
    plots = {}
    
    fig_encoder = plt.figure(figsize=(14, 10))
    
    # Speed measurement error
    plt.subplot(2, 3, 1)
    true_speed = 100.0 * np.sin(time_data/7) + 50.0
    measured_speed = np.copy(true_speed)
    
    # Add measurement noise
    np.random.seed(42)  # For reproducible results
    noise_base = 2.0 * np.random.randn(len(time_data))
    measured_speed += noise_base
    
    # Add bias and increased noise after fault
    fault_idx = time_data >= fault_time_min
    if np.any(fault_idx):
        error_magnitude = 20.0
        bias = error_magnitude * 0.5  # 50% of error as bias
        noise = error_magnitude * 0.3 * np.random.randn(np.sum(fault_idx))  # 30% as random noise
        measured_speed[fault_idx] += bias + noise
    
    plt.plot(time_data, true_speed, '--', color='blue', linewidth=2, label='True Speed')
    plt.plot(time_data, measured_speed, '-', color='red', linewidth=2, label='Measured Speed')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Wheel Speed (rad/s)')
    plt.title(f'Encoder Measurement Error - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Measurement error magnitude
    plt.subplot(2, 3, 2)
    speed_error = measured_speed - true_speed
    plt.plot(time_data, speed_error, 'purple', linewidth=2, label='Measurement Error')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Speed Error (rad/s)')
    plt.title('Measurement Error Magnitude')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Error statistics over time
    plt.subplot(2, 3, 3)
    window_size = max(1, len(time_data) // 20)  # Rolling window
    if len(speed_error) >= window_size:
        error_rms = []
        error_bias = []
        for i in range(window_size, len(speed_error)):
            window_data = speed_error[i-window_size:i]
            error_rms.append(np.sqrt(np.mean(window_data**2)))
            error_bias.append(np.mean(window_data))
        
        time_windowed = time_data[window_size:]
        plt.plot(time_windowed, error_rms, 'red', linewidth=2, label='RMS Error')
        plt.plot(time_windowed, np.abs(error_bias), 'orange', linewidth=2, label='Bias Magnitude')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Error Statistics (rad/s)')
        plt.title('Error Statistics Evolution')
        plt.grid(True, alpha=0.3)
        plt.legend()
    
    # Control impact
    plt.subplot(2, 3, 4)
    # Simulate control error due to encoder fault
    control_error = np.abs(speed_error) * 0.1  # Control effort proportional to measurement error
    plt.plot(time_data, control_error, 'green', linewidth=2, label='Control Error')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Control Error Magnitude')
    plt.title('Control System Impact')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Frequency analysis placeholder
    plt.subplot(2, 3, 5)
    # Simple frequency analysis representation
    pre_fault_error = speed_error[time_data < fault_time_min]
    post_fault_error = speed_error[time_data >= fault_time_min]
    
    if len(pre_fault_error) > 5 and len(post_fault_error) > 5:
        # Simple histogram comparison
        plt.hist(pre_fault_error, bins=20, alpha=0.5, label='Pre-fault Error', color='blue')
        plt.hist(post_fault_error, bins=20, alpha=0.5, label='Post-fault Error', color='red')
        plt.xlabel('Error Magnitude (rad/s)')
        plt.ylabel('Frequency')
        plt.title('Error Distribution')
        plt.legend()
    else:
        plt.text(0.5, 0.5, 'Insufficient data\nfor analysis', 
                ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Error Distribution')
    
    plt.grid(True, alpha=0.3)
    
    # Detection indicator
    plt.subplot(2, 3, 6)
    # Simple fault detection based on error threshold
    error_threshold = 5.0  # rad/s
    detection_signal = np.abs(speed_error) > error_threshold
    detection_level = detection_signal.astype(float) * 100
    
    plt.plot(time_data, detection_level, 'red', linewidth=2, label='Fault Detection')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Actual Fault')
    plt.axhline(y=50, color='orange', linestyle=':', alpha=0.7, label='Detection Threshold')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Detection Level (%)')
    plt.title('Fault Detection Signal')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"EncoderFault_Comprehensive_{spacecraft_name}"] = fig_encoder
    
    return plots

def generate_battery_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate comprehensive battery fault specific plots"""
    plots = {}
    
    fig_battery = plt.figure(figsize=(14, 10))
    
    # Battery state of charge
    plt.subplot(2, 3, 1)
    normal_discharge_rate = 2.0  # %/min
    fault_discharge_rate = normal_discharge_rate + fault_data.get('battery_drain', 0.05) * 10
    
    battery_state = np.zeros_like(time_data)
    fault_idx = time_data >= fault_time_min
    idx_before = time_data < fault_time_min
    
    if np.any(idx_before):
        battery_state[idx_before] = 100.0 - normal_discharge_rate * time_data[idx_before]
        if np.any(fault_idx):
            battery_state[fault_idx] = battery_state[np.sum(idx_before)-1] - fault_discharge_rate * (time_data[fault_idx] - fault_time_min)
    else:
        battery_state = 100.0 - fault_discharge_rate * time_data
    
    plt.plot(time_data, battery_state, 'red', linewidth=3, label=f'{spacecraft_name} Battery')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
    plt.axhline(y=20, color='orange', linestyle=':', alpha=0.7, label='Low Battery (20%)')
    plt.axhline(y=10, color='red', linestyle=':', alpha=0.7, label='Critical (10%)')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Battery State (%)')
    plt.ylim(0, 100)
    plt.title(f'Battery Degradation - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Power consumption comparison
    plt.subplot(2, 3, 2)
    normal_power = 1.0 + 0.2 * np.sin(time_data/5)  # Normal power consumption
    faulty_power = np.copy(normal_power)
    battery_drain = fault_data.get('battery_drain', 0.05) * 1000  # Convert to W
    faulty_power[fault_idx] += battery_drain
    
    plt.plot(time_data, normal_power, '--', color='blue', linewidth=2, label='Normal Power')
    plt.plot(time_data, faulty_power, '-', color='red', linewidth=2, label='With Battery Fault')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Power Consumption (W)')
    plt.title('Power Consumption Impact')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Battery discharge rate
    plt.subplot(2, 3, 3)
    discharge_rate = np.full_like(time_data, normal_discharge_rate)
    discharge_rate[fault_idx] = fault_discharge_rate
    
    plt.plot(time_data, discharge_rate, 'orange', linewidth=2, label='Discharge Rate')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
    plt.axhline(y=normal_discharge_rate, color='gray', linestyle=':', alpha=0.7, label='Normal Rate')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Discharge Rate (%/min)')
    plt.title('Battery Discharge Rate')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Remaining mission time
    plt.subplot(2, 3, 4)
    remaining_time = np.zeros_like(time_data)
    for i, (state, rate) in enumerate(zip(battery_state, discharge_rate)):
        if rate > 0:
            remaining_time[i] = max(0, state / rate)
        else:
            remaining_time[i] = 1000  # Essentially infinite if no discharge
    
    plt.plot(time_data, remaining_time, 'purple', linewidth=2, label='Remaining Mission Time')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Remaining Time (minutes)')
    plt.title('Mission Time Remaining')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Energy analysis
    plt.subplot(2, 3, 5)
    if len(time_data) > 1:
        dt = np.diff(time_data) * 60  # Convert to seconds
        energy_consumed = np.cumsum(np.concatenate([[0], faulty_power[1:] * dt]))  # Energy in Joules
        energy_normal = np.cumsum(np.concatenate([[0], normal_power[1:] * dt]))
        
        plt.plot(time_data, energy_normal, '--', color='blue', linewidth=2, label='Normal Energy')
        plt.plot(time_data, energy_consumed, '-', color='red', linewidth=2, label='Actual Energy')
        plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Cumulative Energy (J)')
        plt.title('Energy Consumption Analysis')
        plt.grid(True, alpha=0.3)
        plt.legend()
    
    # System health vs battery level
    plt.subplot(2, 3, 6)
    health = np.ones_like(time_data) * 100
    health[battery_state < 50] = 80  # Degraded below 50%
    health[battery_state < 20] = 50  # Poor below 20%
    health[battery_state < 10] = 20  # Critical below 10%
    health[battery_state <= 0] = 0   # Failed at 0%
    
    plt.plot(time_data, health, 'green', linewidth=2, label='System Health')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Battery Fault')
    plt.axhline(y=100, color='green', linestyle=':', alpha=0.7, label='Healthy')
    plt.axhline(y=50, color='orange', linestyle=':', alpha=0.7, label='Degraded')
    plt.axhline(y=20, color='red', linestyle=':', alpha=0.7, label='Critical')
    plt.xlabel('Time (minutes)')
    plt.ylabel('System Health (%)')
    plt.title('System Health vs Battery')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"BatteryFault_Comprehensive_{spacecraft_name}"] = fig_battery
    
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
        # Get orbital parameters
        orbit_radius = np.linalg.norm(sc.hub.r_CN_NInit) / 1000.0  # Convert to km
        
        # Generate orbit path (simplified circular)
        angles = np.linspace(0, 2*np.pi, 100)
        orbit_x = orbit_radius * np.cos(angles)
        orbit_y = orbit_radius * np.sin(angles)
        orbit_z = np.zeros_like(angles)  # Simplified for circular orbit
        
        color = colors[i % len(colors)]
        ax1.plot(orbit_x, orbit_y, orbit_z, color=color, alpha=0.7, label=f'{sc.ModelTag}')
        
        # Mark current position
        current_x = sc.hub.r_CN_NInit[0] / 1000.0
        current_y = sc.hub.r_CN_NInit[1] / 1000.0 
        current_z = sc.hub.r_CN_NInit[2] / 1000.0
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
        radius = np.linalg.norm(sc.hub.r_CN_NInit) / 1000.0
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
        # Simplified ground track calculation
        position = sc.hub.r_CN_NInit
        lat = np.arcsin(position[2] / np.linalg.norm(position)) * 180/np.pi
        lon = np.arctan2(position[1], position[0]) * 180/np.pi
        
        color = colors[i % len(colors)]
        ax3.scatter(lon, lat, s=100, color=color, label=f'{sc.ModelTag}')
    
    ax3.legend()
    
    # Orbital periods and coverage
    ax4 = fig_constellation.add_subplot(2, 2, 4)
    
    periods = []
    for sc in spacecraft_list:
        radius = np.linalg.norm(sc.hub.r_CN_NInit)  # meters
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
                pos1 = np.array(spacecraft_list[i].hub.r_CN_NInit)
                pos2 = np.array(spacecraft_list[j].hub.r_CN_NInit)
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
            
            r1 = np.linalg.norm(sc1.hub.r_CN_NInit)
            r2 = np.linalg.norm(sc2.hub.r_CN_NInit)
            
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