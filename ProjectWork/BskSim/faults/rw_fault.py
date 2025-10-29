#!/usr/bin/env python
"""
rw_fault.py

A backward-compatible module that provides access to the various reaction wheel fault scenarios.
This serves as a compatibility layer for existing code that relies on the rw_fault.py interface.
"""
import inspect
import os
import sys
import numpy as np

# Fix path resolution to work with new project structure
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import the specific fault modules - this ensures they are available for the fault loader
from .friction_fault import FrictionFaultScenario, run as run_friction
from .powerlimit_fault import PowerLimitFaultScenario, run as run_powerlimit
from .encoder_fault import EncoderFaultScenario, run as run_encoder


# Create a global fault type that can be changed by users or client code
_current_fault_type = "friction"

def set_fault_type(fault_type):
    """
    Set the global fault type for simulations run through this module
    
    Parameters:
    fault_type (str): Type of fault ('friction', 'power_limit', 'encoder')
    """
    global _current_fault_type
    valid_types = {"friction", "power_limit", "encoder"}
    
    if fault_type not in valid_types:
        raise ValueError(f"Invalid fault type '{fault_type}'. Valid options are: {valid_types}")
    
    _current_fault_type = fault_type
    print(f"Fault type set to: {_current_fault_type}")

def get_fault_type():
    """
    Get the current global fault type
    
    Returns:
    str: Current fault type
    """
    return _current_fault_type

# For backward compatibility, define RWFaultScenario as the FrictionFaultScenario
RWFaultScenario = FrictionFaultScenario

def run(showPlots=True, saveBinary=True, fault_type=None):
    """
    Run a fault scenario with the specified type or the current global fault type
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    fault_type (str, optional): Override the global fault type for this run
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    # Use the specified fault type or fall back to the global setting
    active_fault_type = fault_type or _current_fault_type
    
    # Dispatch to the appropriate run function
    if active_fault_type == "friction":
        return run_friction(showPlots, saveBinary)
    elif active_fault_type == "power_limit":
        return run_powerlimit(showPlots, saveBinary)
    elif active_fault_type == "encoder":
        return run_encoder(showPlots, saveBinary)
    else:
        # This shouldn't happen due to the set_fault_type validation,
        # but we add it as a safety check
        print(f"WARNING: Unknown fault type '{active_fault_type}', using friction fault")
        return run_friction(showPlots, saveBinary)

# Allow direct running of this module with a command line fault type selection
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run a Reaction Wheel Fault Scenario")
    parser.add_argument("--fault-type", choices=["friction", "power_limit", "encoder"], 
                        default="friction", help="Type of fault to simulate")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    
    args = parser.parse_args()
    
    # Set the global fault type so it's available to other modules
    set_fault_type(args.fault_type)
    
    # Run the simulation with the specified fault type
    run(not args.no_plots, not args.no_binary, args.fault_type)