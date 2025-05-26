#!/usr/bin/env python
"""
basilisk_checker.py

Utility to check available modules in Basilisk installation and find alternatives
to simIncludeGravBody.
"""
import os
import sys
import importlib
import pkgutil

def explore_package(package_name):
    """Print all available modules in a package and its subpackages"""
    package = importlib.import_module(package_name)
    print(f"Available modules in {package_name}:")
    
    for _, module_name, is_pkg in pkgutil.iter_modules(package.__path__):
        full_module_name = f"{package_name}.{module_name}"
        print(f"  - {module_name}" + (" (package)" if is_pkg else ""))
        
        # Recurse into subpackages
        if is_pkg:
            try:
                explore_package(full_module_name)
            except ImportError as e:
                print(f"    Error exploring {full_module_name}: {e}")

def find_gravity_modules(package_name="Basilisk"):
    """Search for gravity-related modules in Basilisk"""
    package = importlib.import_module(package_name)
    gravity_modules = []
    
    print(f"Searching for gravity-related modules in {package_name}...")
    
    for finder, module_name, is_pkg in pkgutil.iter_modules(package.__path__):
        full_module_name = f"{package_name}.{module_name}"
        
        if any(keyword in module_name.lower() for keyword in ["grav", "earth", "planet", "body"]):
            gravity_modules.append((full_module_name, is_pkg))
            print(f"  Found: {full_module_name}" + (" (package)" if is_pkg else ""))
        
        # Recurse into subpackages
        if is_pkg:
            try:
                sub_modules = find_gravity_modules(full_module_name)
                gravity_modules.extend(sub_modules)
            except ImportError:
                pass
    
    return gravity_modules

def check_import(module_path, class_name=None):
    """Try to import a module or class and print whether it succeeds"""
    try:
        if class_name:
            module = importlib.import_module(module_path)
            getattr(module, class_name)
            print(f"✓ Successfully imported {class_name} from {module_path}")
        else:
            importlib.import_module(module_path)
            print(f"✓ Successfully imported {module_path}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import {module_path}: {e}")
        return False
    except AttributeError as e:
        print(f"✗ {module_path} exists but doesn't have {class_name}: {e}")
        return False

if __name__ == "__main__":
    print("Basilisk Module Checker\n")
    
    # Check basic Basilisk import
    if not check_import("Basilisk"):
        print("Basilisk package not found! Make sure it's installed and in your PYTHONPATH.")
        sys.exit(1)
    
    # Print Basilisk version if available
    try:
        import Basilisk
        print(f"Basilisk version: {getattr(Basilisk, '__version__', 'Unknown')}")
    except:
        print("Could not determine Basilisk version")
    
    print("\n--- Checking for specific gravity modules ---")
    # Check common alternatives for simIncludeGravBody
    modules_to_check = [
        "Basilisk.simulation.simIncludeGravBody",
        "Basilisk.simulation.gravityEffector",
        "Basilisk.utilities.orbitalMotion",
        "Basilisk.architecture.gravityEffector",
        "Basilisk.dynamics.gravBody"
    ]
    
    for module in modules_to_check:
        check_import(module)
    
    print("\n--- Available modules in Basilisk.simulation ---")
    try:
        explore_package("Basilisk.simulation")
    except ImportError as e:
        print(f"Error exploring Basilisk.simulation: {e}")
    
    print("\n--- Searching for gravity-related modules ---")
    gravity_modules = find_gravity_modules()
    
    print("\n--- Conclusion ---")
    print("Based on the available modules, you likely need to update your imports in spacecraft_simulation.py")
    print("Look for modules related to gravity bodies in the lists above.")