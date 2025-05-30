#!/usr/bin/env python
"""
fault_loader.py - FIXED VERSION

Provides a unified interface for loading and accessing different fault types
from the faults folder and their associated scenario classes.
"""

import os
import sys
import importlib
from typing import Optional, Dict, Any, Tuple

# Get the faults directory path
import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
FAULTS_DIR = os.path.join(ROOT_DIR, 'faults')

# Add faults directory to path
sys.path.insert(0, FAULTS_DIR)

def safe_import_fault_module(module_name, class_name, fallback_class=None):
    """Safely import a fault module from the faults folder with fallback options"""
    try:
        # Try importing from faults folder directly
        module = importlib.import_module(module_name)
        fault_class = getattr(module, class_name, None)
        
        if fault_class is None:
            print(f"Warning: Could not find {class_name} in {module_name}")
            return fallback_class, None
            
        # Look for run function
        run_function = getattr(module, "run", None)
        print(f"Successfully imported {class_name} from {module_name}")
        return fault_class, run_function
        
    except ImportError as e:
        print(f"Warning: Could not import {module_name}: {e}")
        if fallback_class:
            print(f"Using fallback class for {class_name}")
            return fallback_class, None
        return None, None
    except Exception as e:
        print(f"Error importing {module_name}: {e}")
        return None, None

# Import specific fault modules with error handling
fault_modules = {}
run_functions = {}

# Try to import friction fault
print("Loading fault modules from:", FAULTS_DIR)
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

# Create fallback classes if none were imported
if not fault_modules:
    print("Warning: No fault modules could be imported. Creating minimal fallback.")
    
    class FallbackFaultScenario:
        """Minimal fallback fault scenario"""
        def __init__(self):
            self.name = 'FallbackFaultScenario'
            self.fault_type = 'unknown'
            self.fault_magnitude = 0.0
            self.fault_wheel_number = 0
            self.fault_time = 10.0
            print("Using fallback fault scenario - limited functionality")
        
        def run(self):
            print("Fallback scenario run method called")
            return None, None, {}
        
        def pull_outputs(self, showPlots):
            print("Fallback pull_outputs called")
            return {}
    
    # Use fallback for all fault types
    for fault_type in ["friction", "power_limit", "encoder", "battery"]:
        fault_modules[fault_type] = FallbackFaultScenario
        run_functions[fault_type] = lambda: (None, None, {})
else:
    print(f"Successfully loaded {len(fault_modules)} fault modules")

# Dictionary mapping fault types to their corresponding scenario classes
FAULT_SCENARIOS = fault_modules

# Dictionary mapping fault types to their run functions  
RUN_FUNCTIONS = run_functions

# For backward compatibility
RWFaultScenario = fault_modules.get("friction", None)

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
        print(f"Unknown fault type: {fault_type}. Available types: {available_types}")
        # Return friction fault as default
        return FAULT_SCENARIOS.get("friction", FallbackFaultScenario)
    
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
        print(f"Unknown fault type: {fault_type}. Available types: {available_types}")
        return lambda: (None, None, {})  # Return empty function
    
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
        scenario = scenario_class()
        
        # Set common parameters if they exist
        if hasattr(scenario, 'fault_magnitude') and 'magnitude' in kwargs:
            scenario.fault_magnitude = kwargs['magnitude']
        if hasattr(scenario, 'fault_wheel_number') and 'wheel' in kwargs:
            scenario.fault_wheel_number = kwargs['wheel']
        if hasattr(scenario, 'fault_time') and 'time' in kwargs:
            scenario.fault_time = kwargs['time']
        if hasattr(scenario, 'oneTimeFaultTime') and 'time' in kwargs:
            scenario.oneTimeFaultTime = kwargs['time'] * 60e9  # Convert minutes to nanoseconds
            
        return scenario
    except Exception as e:
        print(f"Error creating scenario for {fault_type}: {e}")
        # Return fallback scenario
        return FallbackFaultScenario()

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
        if run_function:
            return run_function(**kwargs)
        else:
            print(f"No run function for {fault_type}")
            return None, None, {}
    except Exception as e:
        print(f"Error running scenario for {fault_type}: {e}")
        import traceback
        traceback.print_exc()
        return None, None, {}

def test_fault_modules():
    """Test all available fault modules"""
    print("\nTesting fault modules...")
    print("="*50)
    print(f"Faults directory: {FAULTS_DIR}")
    print(f"Available fault types: {get_available_fault_types()}")
    
    for fault_type in get_available_fault_types():
        try:
            print(f"\nTesting {fault_type}...")
            scenario_class = get_fault_scenario_class(fault_type)
            print(f"  ✓ Scenario class: {scenario_class.__name__}")
            
            # Try to create instance
            scenario = create_scenario(fault_type, magnitude=0.001, wheel=2, time=10.0)
            print(f"  ✓ Instance created: {scenario.name if hasattr(scenario, 'name') else 'OK'}")
            
            # Check run function
            run_function = get_fault_run_function(fault_type)
            print(f"  ✓ Run function: {'Available' if run_function else 'Not available'}")
            
        except Exception as e:
            print(f"  ✗ Error testing {fault_type}: {e}")

def run(showPlots=True, saveBinary=True, fault_type=None):
    """
    Run a fault scenario with the specified type
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    fault_type (str, optional): Type of fault to run
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    # Default to friction if not specified
    if fault_type is None:
        fault_type = "friction"
    
    # Get the run function for this fault type
    run_func = get_fault_run_function(fault_type)
    
    if run_func:
        return run_func(showPlots=showPlots, saveBinary=saveBinary)
    else:
        print(f"No run function available for fault type: {fault_type}")
        return None, None, {}

# Export key classes and functions
__all__ = [
    'FAULT_SCENARIOS',
    'RUN_FUNCTIONS',
    'get_available_fault_types',
    'is_fault_type_available',
    'get_fault_scenario_class',
    'get_fault_run_function',
    'create_scenario',
    'run_scenario',
    'test_fault_modules',
    'run'
]

if __name__ == "__main__":
    # Test the fault loader when run directly
    test_fault_modules()