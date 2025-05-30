#!/usr/bin/env python
"""
fault_loader.py - ENHANCED VERSION

Provides a unified interface for loading and accessing different fault types
with proper integration for constellation simulations.
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
        return fault_class, run_function
    except ImportError as e:
        print(f"Warning: Could not import {module_name}: {e}")
        if fallback_class:
            print(f"Using fallback class for {class_name}")
            return fallback_class, None
        return None, None
    except AttributeError as e:
        print(f"Warning: Could not find {class_name} in {module_name}: {e}")
        return None, None

# Import specific fault modules with error handling
fault_modules = {}
run_functions = {}
fault_implementations = {}

# Try to import friction fault
friction_class, friction_run = safe_import_fault_module("friction_fault", "FrictionFaultScenario")
if friction_class:
    fault_modules["friction"] = friction_class
    run_functions["friction"] = friction_run

# Try to import power limit fault
powerlimit_class, powerlimit_run = safe_import_fault_module("powerlimit_fault", "PowerLimitFaultScenario")
if powerlimit_class:
    fault_modules["power_limit"] = powerlimit_class
    run_functions["power_limit"] = powerlimit_run

# Try to import encoder fault
encoder_class, encoder_run = safe_import_fault_module("encoder_fault", "EncoderFaultScenario")
if encoder_class:
    fault_modules["encoder"] = encoder_class
    run_functions["encoder"] = encoder_run

# Try to import battery fault
battery_class, battery_run = safe_import_fault_module("battery_fault", "BatteryFaultScenario")
if battery_class:
    fault_modules["battery"] = battery_class
    run_functions["battery"] = battery_run

# Create fault implementation functions for constellation simulations
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
    # This would require modifying the RW model to support power limits
    # For now, we'll store the limit for plotting purposes
    if not hasattr(dynModel, 'power_limits'):
        dynModel.power_limits = {}
    dynModel.power_limits[fault_wheel] = power_limit
    print(f"Applied power limit: {power_limit} W to RW {fault_wheel}")
    return True

def apply_encoder_fault(dynModel, fault_magnitude, fault_wheel, current_time):
    """Apply encoder fault to a reaction wheel"""
    # This would require modifying the RW model to support encoder errors
    # For now, we'll flag it for plotting purposes
    if not hasattr(dynModel, 'encoder_faults'):
        dynModel.encoder_faults = {}
    dynModel.encoder_faults[fault_wheel] = True
    print(f"Applied encoder fault to RW {fault_wheel}")
    return True

def apply_battery_fault(dynModel, drain_rate, fault_wheel, current_time):
    """Apply battery drain fault"""
    # This would require adding battery models to the spacecraft
    # For now, we'll store the drain rate for plotting
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
    """
    Apply a fault to a spacecraft's reaction wheel system
    
    Parameters:
    dynModel: The dynamics model containing the RW system
    fault_type (str): Type of fault to apply
    fault_magnitude (float): Magnitude/severity of the fault
    fault_wheel (int): Which wheel to fault (0-3)
    current_time (float): Current simulation time in nanoseconds
    
    Returns:
    bool: True if fault was successfully applied
    """
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type in fault_implementations:
        impl_func = fault_implementations[fault_type]
        return impl_func(dynModel, fault_magnitude, fault_wheel, current_time)
    else:
        print(f"Warning: Unknown fault type '{fault_type}'")
        return False

def extract_fault_data_from_scenario(scenario, fault_type):
    """
    Extract data from a fault scenario for plotting
    
    Parameters:
    scenario: The fault scenario object
    fault_type (str): Type of fault
    
    Returns:
    dict: Fault data for plotting
    """
    fault_data = {
        'fault_wheel': getattr(scenario, 'fault_wheel_number', 3),
        'fault_time': getattr(scenario, 'oneTimeFaultTime', macros.min2nano(10.))
    }
    
    if fault_type == "friction":
        fault_data['friction_magnitude'] = getattr(scenario, 'fault_magnitude', 0.0005)
        fault_data['friction_baseline'] = 0.02
    elif fault_type == "power_limit":
        fault_data['power_limit'] = getattr(scenario, 'fault_magnitude', 0.5)
    elif fault_type == "battery":
        fault_data['battery_drain'] = getattr(scenario, 'fault_magnitude', 0.05)
    
    # Extract wheel speeds if available
    if hasattr(scenario, 'rwSpeedRec') and scenario.rwSpeedRec:
        fault_data['wheel_speeds'] = scenario.rwSpeedRec.wheelSpeeds
    
    # Extract attitude error if available
    if hasattr(scenario, 'msgRecList') and 'attGuidName' in scenario.msgRecList:
        attErrRec = scenario.msgRecList[scenario.attGuidName]
        if hasattr(attErrRec, 'sigma_BR'):
            sigma_BR = attErrRec.sigma_BR
            fault_data['attitude_error'] = np.linalg.norm(sigma_BR, axis=1)
    
    return fault_data

# Dictionary mapping fault types to their corresponding scenario classes
FAULT_SCENARIOS = fault_modules

# Dictionary mapping fault types to their run functions  
RUN_FUNCTIONS = run_functions

def get_available_fault_types():
    """Get list of available fault types"""
    return list(FAULT_SCENARIOS.keys())

def is_fault_type_available(fault_type):
    """Check if a specific fault type is available"""
    return fault_type.lower().replace(" ", "_") in FAULT_SCENARIOS

def get_fault_scenario_class(fault_type: str):
    """
    Get the appropriate scenario class for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    
    Returns:
    class: The scenario class for the specified fault type
    """
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type not in FAULT_SCENARIOS:
        available_types = list(FAULT_SCENARIOS.keys())
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {available_types}")
    
    return FAULT_SCENARIOS[fault_type]

def create_scenario(fault_type: str, **kwargs):
    """
    Create a scenario instance for the given fault type
    
    Parameters:
    fault_type (str): Type of fault
    **kwargs: Additional parameters to pass to the scenario constructor
    
    Returns:
    object: An instance of the appropriate scenario class
    """
    try:
        scenario_class = get_fault_scenario_class(fault_type)
        return scenario_class(**kwargs)
    except Exception as e:
        print(f"Error creating scenario for {fault_type}: {e}")
        return None

def run_scenario(fault_type: str, **kwargs):
    """
    Run a simulation for the given fault type
    
    Parameters:
    fault_type (str): Type of fault
    **kwargs: Additional parameters to pass to the run function
    
    Returns:
    tuple: Results from the run function (scenario, viz, figureList)
    """
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type in RUN_FUNCTIONS and RUN_FUNCTIONS[fault_type]:
        run_function = RUN_FUNCTIONS[fault_type]
        return run_function(**kwargs)
    else:
        print(f"No run function available for {fault_type}")
        return None, None, {}

# For backward compatibility
try:
    RWFaultScenario = FAULT_SCENARIOS.get("friction")
    run = RUN_FUNCTIONS.get("friction", lambda: (None, None, {}))
except:
    class RWFaultScenario:
        def __init__(self):
            self.name = "RWFaultScenario_Compatibility"
    
    def run():
        return None, None, {}