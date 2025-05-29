#!/usr/bin/env python
"""
run_fault_simulator.py

A simple runner script to launch the Spacecraft Fault Simulator GUI
with all modules properly initialized and robust error handling.
"""
import os
import sys

# Check required modules
required_packages = [
    "numpy", 
    "matplotlib",
    "Basilisk"
]

def check_requirements():
    """Check if required packages are installed"""
    missing_packages = []
    for package in required_packages:
        try:
            if package == "Basilisk":
                # Special case for Basilisk
                from Basilisk import __path__
                print(f"Found Basilisk at: {__path__[0]}")
            else:
                __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("Warning: The following required packages are missing:")
        for package in missing_packages:
            print(f"  - {package}")
        return False
    return True

def setup_environment():
    """Setup environment variables and paths"""
    # Get current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Check if we're in the core directory and adjust path accordingly
    current_dir_name = os.path.basename(current_dir)
    project_root = current_dir
    if current_dir_name == 'core':
        # If running from core directory, use parent as project root
        project_root = os.path.dirname(current_dir)
    
    # Add main directory and subdirectories to Python path
    sys.path.insert(0, project_root)
    sys.path.insert(0, os.path.join(project_root, 'core'))
    sys.path.insert(0, os.path.join(project_root, 'gui_tab'))
    sys.path.insert(0, os.path.join(project_root, 'faults'))
    
    # Create required directories if they don't exist
    directories_to_create = ['logs', 'plots', 'Vizfile', 'faults']  # Changed from 'plotting' to 'plots'
    
    for dir_name in directories_to_create:
        dir_path = os.path.join(project_root, dir_name)
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path, exist_ok=True)
                print(f"Created directory: {dir_path}")
            except Exception as e:
                print(f"Warning: Could not create directory {dir_path}: {e}")
    
    # Also create _VizFiles subdirectory in Vizfile directory
    viz_files_dir = os.path.join(project_root, 'Vizfile', '_VizFiles')
    if not os.path.exists(viz_files_dir):
        try:
            os.makedirs(viz_files_dir, exist_ok=True)
            print(f"Created directory: {viz_files_dir}")
        except Exception as e:
            print(f"Warning: Could not create directory {viz_files_dir}: {e}")
    
    return project_root

def check_file_structure():
    """Check if essential files exist and report any issues"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = current_dir if os.path.basename(current_dir) != 'core' else os.path.dirname(current_dir)
    
    essential_files = [
        'core/spacecraft_simulator_gui.py',
        'core/spacecraft_simulation.py',
        'core/fault_loader.py',
        'core/plots.py',
        'faults/friction_fault.py',
        'faults/powerlimit_fault.py',
        'faults/encoder_fault.py',
        'faults/battery_fault.py'
    ]
    
    missing_files = []
    for file_path in essential_files:
        full_path = os.path.join(project_root, file_path)
        if not os.path.exists(full_path):
            missing_files.append(file_path)
    
    if missing_files:
        print("Warning: Missing essential files:")
        for file_path in missing_files:
            print(f"  - {file_path}")
        return False
    
    return True


    

def main():
    """Main function to run the simulator"""
    print("Starting Spacecraft Fault Simulator...")
    print("=" * 50)
    
    # Setup environment first
    project_root = setup_environment()
    print(f"Project root: {project_root}")
    
    # Check file structure
    if not check_file_structure():
        print("\nSome essential files are missing.")
        proceed = input("Do you want to proceed anyway? (y/n): ")
        if proceed.lower() != 'y':
            print("Exiting. Please ensure al required files are present.")
            return
    
    
    
    # Check requirements
    print("\nChecking dependencies...")
    check_result = check_requirements()
    if not check_result:
        print("\nWarning: Some dependencies are missing.")
        print("The simulator may not work correctly without all required packages.")
        proceed = input("Do you want to proceed anyway? (y/n): ")
        if proceed.lower() != 'y':
            print("Exiting. Please install the required packages and try again.")
            return
    
    print("\nStarting GUI...")
    
    # Import after path setup to ensure all modules are found
    try:
        from spacecraft_simulator_gui import SatelliteSimulatorApp
        print("Successfully imported GUI module")
    except ImportError as e:
        print(f"Error importing GUI module: {e}")
        print("\nTrying to identify the specific issue...")
        
        # Try importing dependencies step by step
        try:
            import tkinter as tk
            print("✓ tkinter imported successfully")
        except ImportError:
            print("✗ tkinter import failed")
            
        try:
            import numpy as np
            print("✓ numpy imported successfully")
        except ImportError:
            print("✗ numpy import failed")
            
        try:
            import matplotlib.pyplot as plt
            print("✓ matplotlib imported successfully")
        except ImportError:
            print("✗ matplotlib import failed")
            
        try:
            from Basilisk.utilities import macros
            print("✓ Basilisk imported successfully")
        except ImportError:
            print("✗ Basilisk import failed")
        
        # Try importing the main modules individually
        try:
            from spacecraft_simulation import SimulationConfig
            print("✓ spacecraft_simulation imported successfully")
        except ImportError as e:
            print(f"✗ spacecraft_simulation import failed: {e}")
            
        try:
            from fault_loader import get_fault_scenario_class
            print("✓ fault_loader imported successfully")
        except ImportError as e:
            print(f"✗ fault_loader import failed: {e}")
            
        try:
            from plots import generate_fault_plots
            print("✓ plots module imported successfully")
        except ImportError as e:
            print(f"✗ plots module import failed: {e}")
        
        print(f"\nDetailed error: {e}")
        print("Please check the error messages above and fix any missing dependencies.")
        input("Press Enter to exit...")
        return
    
    except Exception as e:
        print(f"Unexpected error during import: {e}")
        print("Please check the console output for more details.")
        input("Press Enter to exit...")
        return
    
    # Import tkinter and create GUI
    try:
        import tkinter as tk
        print("Creating GUI window...")
        root = tk.Tk()
        
        # Set up window properties
        root.title("Spacecraft Constellation Fault Simulator")
        root.geometry("1200x900")
        
        print("Initializing application...")
        app = SatelliteSimulatorApp(root)
        
        print("Starting GUI main loop...")
        print("=" * 50)
        print("GUI is now running. Close the window to exit.")
        
        root.mainloop()
        
    except ImportError as e:
        print(f"Error importing tkinter: {e}")
        print("tkinter is required but not available. Please check your Python installation.")
        input("Press Enter to exit...")
        
    except Exception as e:
        print(f"Error running GUI: {e}")
        import traceback
        print("Full traceback:")
        traceback.print_exc()
        
        # Show error dialog if possible
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()  # Hide the root window
            messagebox.showerror("Startup Error", 
                               f"Failed to start the simulator:\n\n{str(e)}\n\nCheck the console for more details.")
        except:
            pass
        
        input("Press Enter to exit...")

if __name__ == "__main__":
    main()