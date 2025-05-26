#!/usr/bin/env python
"""
fault_loader.py - ENHANCED VERSION

Provides a unified interface for loading and accessing different fault types
and their associated scenario classes with robust error handling.
"""

import os
import sys
import importlib
from typing import Optional, Dict, Any, Tuple

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

# Try to import battery fault (with fallback to friction if needed)
battery_class, battery_run = safe_import_fault_module("battery_fault", "BatteryFaultScenario", 
                                                     fallback_class=friction_class)
if battery_class:
    fault_modules["battery"] = battery_class
    run_functions["battery"] = battery_run or friction_run

# Create fallback classes if none were imported
if not fault_modules:
    print("Warning: No fault modules could be imported. Creating minimal fallback.")
    
    class FallbackFaultScenario:
        """Minimal fallback fault scenario"""
        def __init__(self):
            self.name = 'FallbackFaultScenario'
            print("Using fallback fault scenario - limited functionality")
        
        def run(self):
            print("Fallback scenario run method called")
            return None, None, {}
    
    # Use fallback for all fault types
    for fault_type in ["friction", "power_limit", "encoder", "battery"]:
        fault_modules[fault_type] = FallbackFaultScenario
        run_functions[fault_type] = lambda: (None, None, {})

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
    # Convert from GUI name format to module format if needed
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type not in FAULT_SCENARIOS:
        available_types = list(FAULT_SCENARIOS.keys())
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {available_types}")
    
    return FAULT_SCENARIOS[fault_type]

def get_fault_run_function(fault_type: str):
    """
    Get the run function for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    
    Returns:
    function: The run function for the specified fault type
    """
    # Convert from GUI name format to module format if needed
    fault_type = fault_type.lower().replace(" ", "_")
    
    if fault_type not in RUN_FUNCTIONS:
        available_types = list(RUN_FUNCTIONS.keys())
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {available_types}")
    
    run_func = RUN_FUNCTIONS[fault_type]
    if run_func is None:
        print(f"Warning: No run function available for {fault_type}")
        return lambda: (None, None, {})  # Return empty function
    
    return run_func

def create_scenario(fault_type: str, **kwargs):
    """
    Create a scenario instance for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    **kwargs: Additional parameters to pass to the scenario constructor
    
    Returns:
    object: An instance of the appropriate scenario class
    """
    try:
        scenario_class = get_fault_scenario_class(fault_type)
        return scenario_class(**kwargs)
    except Exception as e:
        print(f"Error creating scenario for {fault_type}: {e}")
        # Return fallback scenario
        return FAULT_SCENARIOS[list(FAULT_SCENARIOS.keys())[0]]()

def run_scenario(fault_type: str, **kwargs):
    """
    Run a simulation for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    **kwargs: Additional parameters to pass to the run function
    
    Returns:
    tuple: Results from the run function (usually scenario, viz, figureList)
    """
    try:
        run_function = get_fault_run_function(fault_type)
        return run_function(**kwargs)
    except Exception as e:
        print(f"Error running scenario for {fault_type}: {e}")
        return None, None, {}

def test_fault_modules():
    """Test all available fault modules"""
    print("Testing fault modules...")
    print(f"Available fault types: {get_available_fault_types()}")
    
    for fault_type in get_available_fault_types():
        try:
            print(f"Testing {fault_type}...")
            scenario_class = get_fault_scenario_class(fault_type)
            run_function = get_fault_run_function(fault_type)
            
            print(f"  ✓ Scenario class: {scenario_class.__name__}")
            print(f"  ✓ Run function: {'Available' if run_function else 'Not available'}")
            
        except Exception as e:
            print(f"  ✗ Error testing {fault_type}: {e}")

# For backward compatibility
try:
    RWFaultScenario = FAULT_SCENARIOS.get("friction")
    run = RUN_FUNCTIONS.get("friction", lambda: (None, None, {}))
except:
    # Create minimal compatibility objects
    class RWFaultScenario:
        def __init__(self):
            self.name = "RWFaultScenario_Compatibility"
    
    def run():
        print("Compatibility run function called")
        return None, None, {}

if __name__ == "__main__":
    # Test the fault loader when run directly
    test_fault_modules()