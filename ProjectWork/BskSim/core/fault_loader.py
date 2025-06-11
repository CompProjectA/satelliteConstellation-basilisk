#!/usr/bin/env python
"""
fault_loader.py - Enhanced to pass simulation time to fault modules

This module dynamically loads and runs fault scenarios with proper parameter passing.
"""

import os
import sys
import importlib
import numpy as np
from typing import Optional, Dict, Any, Tuple
from Basilisk.utilities import macros

def safe_import_fault_module(module_name, class_name, fallback_class=None):
    """Safely import a fault module with fallback options"""
    try:
        module = importlib.import_module(f"faults.{module_name}")
        fault_class = getattr(module, class_name)
        run_function = getattr(module, "run", None)
        run_with_parameters = getattr(module, "run_with_parameters", None)
        return fault_class, run_function, run_with_parameters
    except ImportError as e:
        print(f"Warning: Could not import {module_name}: {e}")
        if fallback_class:
            print(f"Using fallback class for {class_name}")
            return fallback_class, None, None
        return None, None, None
    except AttributeError as e:
        print(f"Warning: Could not find {class_name} in {module_name}: {e}")
        return None, None, None

# Import fault modules with enhanced parameter support
fault_modules = {}
run_functions = {}
run_with_parameters_functions = {}
fault_implementations = {}

# Try to import friction fault
friction_class, friction_run, friction_run_with_params = safe_import_fault_module("friction_fault", "FrictionFaultScenario")
if friction_class:
    fault_modules["friction"] = friction_class
    run_functions["friction"] = friction_run
    run_with_parameters_functions["friction"] = friction_run_with_params

# Try to import power limit fault
powerlimit_class, powerlimit_run, powerlimit_run_with_params = safe_import_fault_module("powerlimit_fault", "PowerLimitFaultScenario")
if powerlimit_class:
    fault_modules["power_limit"] = powerlimit_class
    run_functions["power_limit"] = powerlimit_run
    run_with_parameters_functions["power_limit"] = powerlimit_run_with_params

# Try to import encoder fault
encoder_class, encoder_run, encoder_run_with_params = safe_import_fault_module("encoder_fault", "EncoderFaultScenario")
if encoder_class:
    fault_modules["encoder"] = encoder_class
    run_functions["encoder"] = encoder_run
    run_with_parameters_functions["encoder"] = encoder_run_with_params

# Try to import battery fault
battery_class, battery_run, battery_run_with_params = safe_import_fault_module("battery_fault", "BatteryFaultScenario")
if battery_class:
    fault_modules["battery"] = battery_class
    run_functions["battery"] = battery_run
    run_with_parameters_functions["battery"] = battery_run_with_params
    print("Successfully imported battery fault module")
else:
    print("Failed to import battery fault module")

# Fault implementation functions
def apply_friction_fault(dynModel, fault_magnitude, fault_wheel, current_time):
    """Apply friction fault to a reaction wheel"""
    if hasattr(dynModel, 'rwFactory'):
        rw_list = [dynModel.RW1, dynModel.RW2, dynModel.RW3, dynModel.RW4]
        if 0 <= fault_wheel < len(rw_list) and rw_list[fault_wheel]:
            rw_list[fault_wheel].fCoulomb += fault_magnitude
            print(f"Applied friction fault: {fault_magnitude} N⋅m to RW {fault_wheel}")
            return True
    return False

def apply_power_limit_fault(dynModel, power_limit, fault_wheel, current_time):
    """Apply power limit fault to a reaction wheel"""
    if not hasattr(dynModel, 'power_limits'):
        dynModel.power_limits = {}
    dynModel.power_limits[fault_wheel] = power_limit
    print(f"Applied power limit: {power_limit} W to RW {fault_wheel}")
    return True

def apply_encoder_fault(dynModel, fault_magnitude, fault_wheel, current_time):
    """Apply encoder fault to a reaction wheel"""
    if not hasattr(dynModel, 'encoder_faults'):
        dynModel.encoder_faults = {}
    dynModel.encoder_faults[fault_wheel] = True
    print(f"Applied encoder fault to RW {fault_wheel}")
    return True

def apply_battery_fault(dynModel, drain_rate, fault_wheel, current_time):
    """Apply battery drain fault"""
    if not hasattr(dynModel, 'battery_drain'):
        dynModel.battery_drain = 0.0
    dynModel.battery_drain += drain_rate
    print(f"Applied battery drain: {drain_rate} kW")
    return True

# Map fault types to implementation functions
fault_implementations = {
    "friction": apply_friction_fault,
    "power_limit": apply_power_limit_fault,
    "encoder": apply_encoder_fault,
    "battery": apply_battery_fault
}

def apply_fault_to_spacecraft(dynModel, fault_type, fault_magnitude, fault_wheel, current_time):
    """Apply a fault to a spacecraft's reaction wheel system"""
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type in fault_implementations:
        impl_func = fault_implementations[fault_type]
        return impl_func(dynModel, fault_magnitude, fault_wheel, current_time)
    else:
        print(f"Warning: Unknown fault type '{fault_type}'")
        return False

def extract_fault_data_from_scenario(scenario, fault_type):
    """Extract data from a fault scenario for plotting"""
    fault_data = {}
    
    try:
        # Extract basic fault parameters
        if hasattr(scenario, 'fault_wheel_number'):
            fault_data['fault_wheel'] = scenario.fault_wheel_number
        elif hasattr(scenario, 'fault_wheel'):
            fault_data['fault_wheel'] = scenario.fault_wheel
        else:
            fault_data['fault_wheel'] = 3  # Default
            
        if hasattr(scenario, 'oneTimeFaultTime'):
            fault_data['fault_time'] = scenario.oneTimeFaultTime
        else:
            fault_data['fault_time'] = 10.0  # Default
        
        # Extract fault-specific parameters
        if fault_type == "friction":
            if hasattr(scenario, 'fault_magnitude'):
                fault_data['friction_magnitude'] = scenario.fault_magnitude
            else:
                fault_data['friction_magnitude'] = 0.0005  # Default
            fault_data['friction_baseline'] = 0.02
            
        elif fault_type == "power_limit":
            if hasattr(scenario, 'fault_magnitude'):
                fault_data['power_limit'] = scenario.fault_magnitude
            else:
                fault_data['power_limit'] = 0.5  # Default
                
        elif fault_type == "encoder":
            if hasattr(scenario, 'fault_magnitude'):
                fault_data['encoder_error'] = scenario.fault_magnitude
            else:
                fault_data['encoder_error'] = 20.0  # Default
                
        elif fault_type == "battery":
            if hasattr(scenario, 'fault_magnitude'):
                fault_data['battery_drain'] = scenario.fault_magnitude
            else:
                fault_data['battery_drain'] = 0.05  # Default
        
        # Extract wheel speeds if available
        if hasattr(scenario, 'rwSpeedRec') and scenario.rwSpeedRec:
            try:
                if hasattr(scenario.rwSpeedRec, 'wheelSpeeds'):
                    import numpy as np
                    wheel_speeds_raw = scenario.rwSpeedRec.wheelSpeeds
                    if len(wheel_speeds_raw) > 1:
                        fault_data['wheel_speeds'] = np.delete(wheel_speeds_raw, 0, 0)
                        print(f"Extracted wheel speeds with shape: {fault_data['wheel_speeds'].shape}")
            except Exception as e:
                print(f"Could not extract wheel speeds: {e}")
        
        # Extract attitude error if available
        if hasattr(scenario, 'msgRecList'):
            try:
                if hasattr(scenario, 'attGuidName') and scenario.attGuidName in scenario.msgRecList:
                    attErrRec = scenario.msgRecList[scenario.attGuidName]
                    if hasattr(attErrRec, 'sigma_BR'):
                        import numpy as np
                        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
                        fault_data['attitude_error'] = np.linalg.norm(sigma_BR, axis=1)
                        print(f"Extracted attitude error data with {len(fault_data['attitude_error'])} points")
            except Exception as e:
                print(f"Could not extract attitude error: {e}")
        
        print(f"Extracted fault data keys: {list(fault_data.keys())}")
        
    except Exception as e:
        print(f"ERROR: Failed to extract fault data: {e}")
        import traceback
        traceback.print_exc()
    
    return fault_data

# Keep existing dictionaries for backward compatibility
FAULT_SCENARIOS = fault_modules
RUN_FUNCTIONS = run_functions

def get_available_fault_types():
    """Get list of available fault types"""
    return list(FAULT_SCENARIOS.keys())

def is_fault_type_available(fault_type):
    """Check if a specific fault type is available"""
    return fault_type.lower().replace(" ", "_") in FAULT_SCENARIOS

def get_fault_scenario_class(fault_type: str):
    """Get the appropriate scenario class for the given fault type"""
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type not in FAULT_SCENARIOS:
        available_types = list(FAULT_SCENARIOS.keys())
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {available_types}")
    
    return FAULT_SCENARIOS[fault_type]

def create_scenario(fault_type: str, **kwargs):
    """Create a scenario instance for the given fault type"""
    try:
        scenario_class = get_fault_scenario_class(fault_type)
        return scenario_class(**kwargs)
    except Exception as e:
        print(f"Error creating scenario for {fault_type}: {e}")
        return None

def run_scenario_enhanced(fault_type: str, **kwargs):
    """
    Enhanced run_scenario that properly handles battery capacity reduction
    """
    fault_type = fault_type.lower().replace(" ", "_")
    
    # Extract standard parameters
    show_plots = kwargs.get('showPlots', False)
    save_binary = kwargs.get('saveBinary', False)
    
    # Extract simulation parameters
    sim_params = kwargs.get('simulation_params', {})
    fault_magnitude = sim_params.get('fault_magnitude', kwargs.get('fault_magnitude', 0.0005))
    fault_wheel = sim_params.get('fault_wheel', kwargs.get('fault_wheel', 3))
    fault_time_min = sim_params.get('fault_time_min', kwargs.get('fault_time_min', 10.0))
    simulation_time_min = sim_params.get('simulation_time_min', kwargs.get('simulation_time_min', 30.0))
    
    # Special handling for battery faults
    if fault_type == "battery":
        # Check if this is a capacity reduction (magnitude between 10-90)
        if 10 <= fault_magnitude <= 90:
            battery_capacity_percentage = fault_magnitude
            print(f"Detected battery capacity reduction: {battery_capacity_percentage}%")
            
            # Add capacity percentage to kwargs
            kwargs['battery_capacity_percentage'] = battery_capacity_percentage
            sim_params['battery_capacity_percentage'] = battery_capacity_percentage
        else:
            # Regular power drain fault
            print(f"Battery power drain fault: {fault_magnitude}W")
    
    print(f"Running enhanced {fault_type} scenario with parameters:")
    print(f"  - showPlots: {show_plots}")
    print(f"  - saveBinary: {save_binary}")
    print(f"  - fault_magnitude: {fault_magnitude}")
    print(f"  - fault_wheel: {fault_wheel}")
    print(f"  - fault_time_min: {fault_time_min}")
    print(f"  - simulation_time_min: {simulation_time_min}")
    
    # Call the appropriate run function
    if fault_type in run_with_parameters_functions and run_with_parameters_functions[fault_type]:
        print(f"Using run_with_parameters() for {fault_type}")
        
        try:
            run_with_params_func = run_with_parameters_functions[fault_type]
            
            # For battery fault, include capacity percentage if available
            if fault_type == "battery" and 'battery_capacity_percentage' in kwargs:
                result = run_with_params_func(
                    fault_magnitude=fault_magnitude,
                    fault_wheel=fault_wheel,
                    fault_time_min=fault_time_min,
                    simulation_time_min=simulation_time_min,
                    showPlots=show_plots,
                    saveBinary=save_binary,
                    battery_capacity_percentage=kwargs['battery_capacity_percentage']
                )
            else:
                result = run_with_params_func(
                    fault_magnitude=fault_magnitude,
                    fault_wheel=fault_wheel,
                    fault_time_min=fault_time_min,
                    simulation_time_min=simulation_time_min,
                    showPlots=show_plots,
                    saveBinary=save_binary
                )
            
            print(f"run_with_parameters() returned: {type(result)}")
            
            # Handle the result
            if result is None:
                print(f"run_with_parameters returned None")
                return None, None, {}
            elif isinstance(result, tuple):
                if len(result) >= 2:
                    return result[0], result[1], result[2] if len(result) > 2 else {}
                else:
                    return result[0] if len(result) > 0 else None, None, {}
            else:
                return result, None, {}
                
        except Exception as e:
            print(f"ERROR: run_with_parameters failed for {fault_type}: {e}")
            import traceback
            traceback.print_exc()
    
    # Fall back to regular run if enhanced version not available
    return run_scenario(fault_type, **kwargs)

def run_scenario(fault_type: str, **kwargs):
    """
    Run a simulation using run_with_parameters() when available
    
    Enhanced to accept simulation_time_min parameter
    """
    fault_type = fault_type.lower().replace(" ", "_")
    
    # Extract standard parameters
    show_plots = kwargs.get('showPlots', False)
    save_binary = kwargs.get('saveBinary', False)
    
    # Extract simulation parameters for scenario configuration
    sim_params = kwargs.get('simulation_params', {})
    fault_magnitude = sim_params.get('fault_magnitude', kwargs.get('fault_magnitude', 0.0005))
    fault_wheel = sim_params.get('fault_wheel', kwargs.get('fault_wheel', 3))
    fault_time_min = sim_params.get('fault_time_min', kwargs.get('fault_time_min', 10.0))
    
    # Extract simulation duration
    simulation_time_min = sim_params.get('simulation_time_min', kwargs.get('simulation_time_min', 30.0))
    
    print(f"Running {fault_type} with parameters:")
    print(f"  - showPlots: {show_plots}")
    print(f"  - saveBinary: {save_binary}")
    print(f"  - fault_magnitude: {fault_magnitude}")
    print(f"  - fault_wheel: {fault_wheel}")
    print(f"  - fault_time_min: {fault_time_min}")
    print(f"  - simulation_time_min: {simulation_time_min}")
    
    # Check if the enhanced run_with_parameters function is available
    if fault_type in run_with_parameters_functions and run_with_parameters_functions[fault_type]:
        print(f"Using run_with_parameters() for {fault_type}")
        
        try:
            run_with_params_func = run_with_parameters_functions[fault_type]
            
            # Call the enhanced function with all parameters including simulation time
            result = run_with_params_func(
                fault_magnitude=fault_magnitude,
                fault_wheel=fault_wheel,
                fault_time_min=fault_time_min,
                simulation_time_min=simulation_time_min,
                showPlots=show_plots,
                saveBinary=save_binary
            )
            
            print(f"run_with_parameters() returned: {type(result)}")
            
            # Handle the result
            if result is None:
                print(f"run_with_parameters returned None")
                return None, None, {}
            elif isinstance(result, tuple):
                if len(result) >= 2:
                    return result[0], result[1], result[2] if len(result) > 2 else {}
                else:
                    return result[0] if len(result) > 0 else None, None, {}
            else:
                return result, None, {}
                
        except Exception as e:
            print(f"ERROR: run_with_parameters failed for {fault_type}: {e}")
            import traceback
            traceback.print_exc()
            # Fall back to regular run function
    
    # Fall back to regular run function if run_with_parameters is not available
    if fault_type in RUN_FUNCTIONS and RUN_FUNCTIONS[fault_type]:
        print(f"Falling back to regular run() for {fault_type}")
        
        run_function = RUN_FUNCTIONS[fault_type]
        
        try:
            import inspect
            run_sig = inspect.signature(run_function)
            run_params = {}
            
            if 'showPlots' in run_sig.parameters:
                run_params['showPlots'] = show_plots
            if 'saveBinary' in run_sig.parameters:
                run_params['saveBinary'] = save_binary
            if 'simulation_time_min' in run_sig.parameters:
                run_params['simulation_time_min'] = simulation_time_min
                
            result = run_function(**run_params)
            
            # Handle different return formats
            if result is None:
                print(f"Run function returned None for {fault_type}")
                return None, None, {}
            elif isinstance(result, tuple) and len(result) >= 2:
                if len(result) == 2:
                    return result[0], result[1], {}
                else:
                    return result[0], result[1], result[2] if len(result) > 2 else {}
            else:
                return result, None, {}
                
        except Exception as e:
            print(f"ERROR: Failed to run {fault_type} scenario: {e}")
            import traceback
            traceback.print_exc()
            return None, None, {}
    else:
        print(f"No run function available for {fault_type}")
        return None, None, {}


# Print diagnostic information
print("fault_loader.py loaded successfully")
print(f"Available fault types: {list(fault_modules.keys())}")
print(f"Enhanced run_with_parameters available for: {[k for k, v in run_with_parameters_functions.items() if v]}")
print(f"Regular run functions available for: {[k for k, v in run_functions.items() if v]}")