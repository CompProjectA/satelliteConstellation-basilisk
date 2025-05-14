#!/usr/bin/env python
"""
fault_loader.py - FIXED VERSION

Provides a unified interface for loading and accessing different fault types
and their associated scenario classes.
"""

import os
import sys
import importlib
from typing import Optional, Dict, Any, Tuple

# Import specific fault modules
try:
    # Try to import fault scenario classes from different modules
    from faults.friction_fault import FrictionFaultScenario
    from faults.powerlimit_fault import PowerLimitFaultScenario
    from faults.encoder_fault import EncoderFaultScenario
    
    # Try to import BatteryFaultScenario, but provide a fallback if it fails
    try:
        from faults.battery_fault import BatteryFaultScenario
    except ImportError:
        # Define a fallback BatteryFaultScenario that inherits from FrictionFaultScenario
        print("Warning: battery_fault.py module not found, using FrictionFaultScenario as fallback")
        class BatteryFaultScenario(FrictionFaultScenario):
            """Fallback BatteryFaultScenario that inherits from FrictionFaultScenario"""
            def __init__(self):
                super(BatteryFaultScenario, self).__init__()
                self.name = 'BatteryFaultScenario (Fallback)'
    
    # Fallback import for backward compatibility
    from faults.rw_fault import run as run_rw_fault
except ImportError as e:
    print(f"WARNING: Some fault modules could not be imported: {e}")
    # Define fallback classes if needed

# Dictionary mapping fault types to their corresponding scenario classes
FAULT_SCENARIOS = {
    "friction": FrictionFaultScenario,
    "power_limit": PowerLimitFaultScenario,
    "encoder": EncoderFaultScenario,
    "battery": BatteryFaultScenario
}

# Dictionary mapping fault types to their run functions
try:
    from faults.friction_fault import run as run_friction
    from faults.powerlimit_fault import run as run_powerlimit
    from faults.encoder_fault import run as run_encoder
    
    # Try to import battery run function, but provide a fallback if it fails
    try:
        from faults.battery_fault import run as run_battery
    except ImportError:
        # Use friction run as fallback for battery
        run_battery = run_friction
    
    # Default run functions, if the scenario classes don't provide their own
    RUN_FUNCTIONS = {
        "friction": run_friction,
        "power_limit": run_powerlimit,
        "encoder": run_encoder,
        "battery": run_battery
    }
except ImportError:
    print("WARNING: Some fault run functions could not be imported")
    RUN_FUNCTIONS = {
        "friction": getattr(FrictionFaultScenario, "run", run_rw_fault),
        "power_limit": getattr(PowerLimitFaultScenario, "run", run_rw_fault),
        "encoder": getattr(EncoderFaultScenario, "run", run_rw_fault),
        "battery": getattr(BatteryFaultScenario, "run", run_rw_fault)
    }

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
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {list(FAULT_SCENARIOS.keys())}")
    
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
        raise ValueError(f"Unknown fault type: {fault_type}. Available types: {list(RUN_FUNCTIONS.keys())}")
    
    return RUN_FUNCTIONS[fault_type]

def create_scenario(fault_type: str, **kwargs):
    """
    Create a scenario instance for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    **kwargs: Additional parameters to pass to the scenario constructor
    
    Returns:
    object: An instance of the appropriate scenario class
    """
    scenario_class = get_fault_scenario_class(fault_type)
    return scenario_class(**kwargs)

def run_scenario(fault_type: str, **kwargs):
    """
    Run a simulation for the given fault type
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder', 'battery')
    **kwargs: Additional parameters to pass to the run function
    
    Returns:
    tuple: Results from the run function (usually scenario, viz, figureList)
    """
    run_function = get_fault_run_function(fault_type)
    return run_function(**kwargs)

# For backward compatibility
RWFaultScenario = FrictionFaultScenario
run = run_rw_fault