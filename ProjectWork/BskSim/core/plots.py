#!/usr/bin/env python
"""
plots.py

Centralized plotting module for the Spacecraft Constellation Fault Simulator.
Handles all plot generation from fault modules and constellation data.
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

# Import fault loader for real simulation integration
try:
    from fault_loader import (
        run_scenario,
        extract_fault_data_from_scenario,
        get_available_fault_types,
        is_fault_type_available
    )
    FAULT_LOADER_AVAILABLE = True
    print("Successfully imported fault_loader")
except ImportError as e:
    print(f"Warning: Could not import fault_loader: {e}")
    FAULT_LOADER_AVAILABLE = False

def generate_fault_plots(fault_type, fault_data, time_data, fault_time_min, spacecraft_name="Spacecraft"):
    """
    Generate plots specific to the fault type using real data from fault_loader
    
    Enhanced to use fault_loader.py for real simulation data when available
    """
    plots = {}
    
    # Calculate simulation time from time_data
    simulation_time_min = time_data[-1] if len(time_data) > 0 else 30.0
    
    # If fault_data contains a 'use_real_simulation' flag, run the actual simulation
    if isinstance(fault_data, dict) and fault_data.get('use_real_simulation', False) and FAULT_LOADER_AVAILABLE:
        print(f"Running real {fault_type} simulation for {simulation_time_min} minutes...")
        
        try:
            # Extract simulation parameters
            simulation_params = fault_data.get('simulation_params', {})
            # Add simulation time to parameters
            simulation_params['simulation_time_min'] = simulation_time_min
            
            # Run the actual fault simulation using fault_loader
            print(f"Calling run_scenario for {fault_type} with {simulation_time_min} minute duration...")
            result = run_scenario(
                fault_type,
                showPlots=False,  # Don't show plots, we'll handle them
                saveBinary=False,  # Don't save binary files
                **simulation_params
            )
            
            # Handle return values
            if result is None:
                print(f"run_scenario returned None for {fault_type}")
                scenario, viz, figure_list = None, None, {}
            elif isinstance(result, tuple):
                if len(result) == 3:
                    scenario, viz, figure_list = result
                elif len(result) == 2:
                    scenario, viz = result
                    figure_list = {}
                else:
                    print(f"Unexpected tuple length from run_scenario: {len(result)}")
                    scenario, viz, figure_list = result[0] if len(result) > 0 else None, None, {}
            else:
                # If result is not a tuple, treat it as a scenario object
                print(f"run_scenario returned non-tuple: {type(result)}")
                scenario = result
                viz = None
                figure_list = {}
            
            print(f"Processed run_scenario result - scenario: {scenario is not None}, "
                  f"viz: {viz is not None}, figure_list: {len(figure_list) if figure_list else 0}")
            
            if scenario is not None:
                print(f"Successfully ran {fault_type} simulation for {simulation_time_min} minutes")
                
                # Extract data from the completed scenario
                real_fault_data = extract_fault_data_from_scenario(scenario, fault_type)
                
                # If the scenario has its own plots in figure_list, use those
                if figure_list and len(figure_list) > 0:
                    for plot_name, fig in figure_list.items():
                        new_name = f"REAL_{plot_name}_{spacecraft_name}"  # PREFIX with REAL_
                        plots[new_name] = fig
                    return plots
                
                # Check if scenario has a pull_outputs method
                elif hasattr(scenario, 'pull_outputs'):
                    try:
                        print(f"Attempting to extract plots using pull_outputs...")
                        scenario_plots = scenario.pull_outputs(showPlots=False)
                        if scenario_plots and len(scenario_plots) > 0:
                            print(f"Extracted {len(scenario_plots)} plots using pull_outputs")
                            # Rename plots to include spacecraft name
                            for plot_name, fig in scenario_plots.items():
                                new_name = f"{plot_name}_{spacecraft_name}"
                                plots[new_name] = fig
                            return plots
                        else:
                            print(f"pull_outputs returned empty or None")
                    except Exception as e:
                        print(f"pull_outputs failed: {e}")
                
                # If no plots available from scenario, generate from extracted data
                print(f"Generating plots from extracted real data")
                real_plots = generate_plots_from_real_data(
                    fault_type, real_fault_data, scenario, fault_time_min, spacecraft_name
                )
                if real_plots and len(real_plots) > 0:
                    print(f"Generated {len(real_plots)} plots from real data")
                    return real_plots
                else:
                    print(f"No plots generated from real data, falling back to synthetic")
            else:
                print(f"Failed to run {fault_type} simulation (scenario is None), falling back to synthetic")
                
        except Exception as e:
            print(f"ERROR: Failed to run real simulation for {fault_type}: {e}")
            import traceback
            traceback.print_exc()
            print(f"Falling back to synthetic plots...")
    
    # Check if fault_data is actually a scenario object with pull_outputs method
    if hasattr(fault_data, 'pull_outputs'):
        try:
            print(f"Extracting plots from existing {fault_type} scenario...")
            scenario_plots = fault_data.pull_outputs(showPlots=False)
            # Rename plots to include spacecraft name
            for plot_name, fig in scenario_plots.items():
                new_name = f"{plot_name}_{spacecraft_name}"
                plots[new_name] = fig
            return plots
        except Exception as e:
            print(f"Could not extract plots from scenario: {e}")
    
    # Fall back to synthetic plotting functions
    print(f"Using synthetic plotting for {fault_type}")
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


def generate_plots_from_real_data(fault_type, real_fault_data, scenario, fault_time_min, spacecraft_name):
    """
    Generate plots using real data extracted from a completed simulation scenario
    """
    plots = {}
    
    print(f"Generating plots from real {fault_type} fault data")
    print(f"Available data keys: {list(real_fault_data.keys()) if real_fault_data else 'None'}")
    
    try:
        if fault_type == "friction":
            plots.update(generate_real_friction_plots(real_fault_data, scenario, fault_time_min, spacecraft_name))
        else:
            print(f"Real plotting not implemented for {fault_type}, using synthetic fallback")
            # Create default time data and fall back to synthetic
            time_data = np.linspace(0, 30, 100)
            plots.update(generate_generic_fault_plots(real_fault_data, time_data, fault_time_min, spacecraft_name))
            
    except Exception as e:
        print(f"Error generating real {fault_type} plots: {e}")
        import traceback
        traceback.print_exc()
        # Fall back to generic plots
        time_data = np.linspace(0, 30, 100)
        plots.update(generate_generic_fault_plots(real_fault_data, time_data, fault_time_min, spacecraft_name))
    
    return plots


def generate_real_friction_plots(real_fault_data, scenario, fault_time_min, spacecraft_name):
    """Generate plots from real friction fault simulation data"""
    plots = {}
    
    # Placeholder for real friction plotting implementation
    # This would use the actual data from the simulation
    fig = plt.figure(figsize=(12, 8))
    plt.text(0.5, 0.5, f"Real Friction Fault Analysis\n{spacecraft_name}", 
             ha='center', va='center', transform=plt.gca().transAxes)
    plt.title("Friction Fault Analysis (Real Data)")
    
    plots[f"REAL_FrictionAnalysis_{spacecraft_name}"] = fig
    return plots


def create_fault_config_for_real_simulation(fault_type, fault_params):
    """
    Create a fault configuration dictionary that triggers real simulation
    
    Args:
        fault_type (str): Type of fault to simulate
        fault_params (dict): Parameters for the fault
        
    Returns:
        dict: Configuration that will trigger real simulation in generate_fault_plots
    """
    return {
        'use_real_simulation': True,
        'simulation_params': fault_params,
        'fault_type': fault_type
    }


def generate_friction_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate friction fault plots - synthetic"""
    plots = {}
    
    fig = plt.figure(figsize=(12, 8))
    
    # Friction torque over time
    plt.subplot(2, 2, 1)
    baseline_friction = 0.02
    friction_magnitude = fault_data.get('friction_magnitude', 0.0005)
    
    friction = np.ones_like(time_data) * baseline_friction
    fault_idx = time_data >= fault_time_min
    friction[fault_idx] += friction_magnitude
    
    plt.plot(time_data, friction, 'b-', linewidth=2)
    plt.axvline(x=fault_time_min, color='red', linestyle='--', linewidth=2, label='Fault Injection')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Friction Torque (N·m)')
    plt.title(f'Friction Fault - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"FrictionFault_Analysis_{spacecraft_name}"] = fig
    
    return plots


def generate_power_limit_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate power limit fault plots - synthetic"""
    plots = {}
    
    fig = plt.figure(figsize=(12, 8))
    
    # Power consumption over time
    plt.subplot(2, 2, 1)
    normal_power = 10.0  # Watts
    power_limit = fault_data.get('power_limit', 0.5)
    
    power = np.ones_like(time_data) * normal_power
    fault_idx = time_data >= fault_time_min
    power[fault_idx] = np.minimum(power[fault_idx], power_limit)
    
    plt.plot(time_data, power, 'g-', linewidth=2)
    plt.axvline(x=fault_time_min, color='red', linestyle='--', linewidth=2, label='Power Limit Applied')
    plt.axhline(y=power_limit, color='orange', linestyle=':', label=f'Limit: {power_limit}W')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Power (W)')
    plt.title(f'Power Limit Fault - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"PowerLimitFault_Analysis_{spacecraft_name}"] = fig
    
    return plots


def generate_encoder_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate encoder fault plots - synthetic"""
    plots = {}
    
    fig = plt.figure(figsize=(12, 8))
    
    # Encoder error over time
    plt.subplot(2, 2, 1)
    encoder_error = fault_data.get('encoder_error', 20.0)  # percentage
    
    error = np.zeros_like(time_data)
    fault_idx = time_data >= fault_time_min
    error[fault_idx] = encoder_error * (1 + 0.1 * np.random.randn(np.sum(fault_idx)))
    
    plt.plot(time_data, error, 'r-', linewidth=2)
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Encoder Fault')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Encoder Error (%)')
    plt.title(f'Encoder Fault - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"EncoderFault_Analysis_{spacecraft_name}"] = fig
    
    return plots


def generate_battery_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """Generate battery fault plots - synthetic"""
    plots = {}
    
    fig = plt.figure(figsize=(12, 8))
    
    # Battery state of charge
    plt.subplot(2, 2, 1)
    initial_charge = 50.0  # Wh
    battery_drain = fault_data.get('battery_drain', 50.0)  # W
    
    # Simple discharge model
    charge = np.ones_like(time_data) * initial_charge
    for i in range(1, len(time_data)):
        dt = (time_data[i] - time_data[i-1]) * 60  # Convert to seconds
        if time_data[i] >= fault_time_min:
            charge[i] = max(0, charge[i-1] - battery_drain * dt / 3600)  # Wh
        else:
            charge[i] = charge[i-1] - 10 * dt / 3600  # Normal drain
    
    plt.plot(time_data, charge, 'b-', linewidth=2)
    plt.axvline(x=fault_time_min, color='red', linestyle='--', linewidth=2, label='Battery Fault')
    plt.axhline(y=10, color='orange', linestyle=':', label='Low Battery Warning')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Battery Charge (Wh)')
    plt.title(f'Battery Fault - {spacecraft_name}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plots[f"BatteryFault_Analysis_{spacecraft_name}"] = fig
    
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
        # Get orbital parameters - handle list format
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
        # Handle the position data format properly
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
        # Simplified ground track calculation - handle position
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
        # Handle position format
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
                # Handle position format properly
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
            
            # Get radii - handle position
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