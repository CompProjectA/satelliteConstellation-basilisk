#!/usr/bin/env python
"""
run_fault_simulator.py

A simple runner script to launch the Spacecraft Fault Simulator GUI
with all modules properly initialized and robust error handling.
"""
import os
import sys

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
    directories_to_create = ['logs', 'plots', 'Vizfile', 'faults']
    
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

def check_basilisk():
    """Check if Basilisk is available"""
    try:
        from Basilisk import __path__
        print(f"Found Basilisk at: {__path__[0]}")
        return True
    except ImportError:
        print("ERROR: Basilisk not found. Please ensure it's installed and in PYTHONPATH")
        return False

def main():
    """Main function to run the simulator"""
    print("Starting Spacecraft Fault Simulator...")
    print("=" * 50)
    
    # Setup environment first
    project_root = setup_environment()
    print(f"Project root: {project_root}")
    
    # Check Basilisk
    print("\nChecking dependencies...")
    if not check_basilisk():
        print("Cannot continue without Basilisk.")
        input("Press Enter to exit...")
        return
    
    print("\nStarting GUI...")
    
    # Import after path setup to ensure all modules are found
    try:
        # This will trigger the fault module checking in spacecraft_simulation.py
        from spacecraft_simulator_gui import SatelliteSimulatorApp
        print("Successfully imported GUI module")
    except ImportError as e:
        print(f"Error importing GUI module: {e}")
        print("\nPlease check that all required files are present.")
        input("Press Enter to exit...")
        return
    except Exception as e:
        print(f"Unexpected error during import: {e}")
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