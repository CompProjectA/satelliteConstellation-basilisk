#!/usr/bin/env python
"""
spacecraft_fault_cli.py

A simplified command-line interface for running Basilisk spacecraft simulations
with reaction wheel fault injection capabilities.
"""
import os
import sys
import inspect
import argparse
import numpy as np
import textwrap
from datetime import datetime

# Fix path resolution to work with new project structure
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
FAULTS_DIR = os.path.join(ROOT_DIR, 'faults')
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')
LOGS_DIR = os.path.join(ROOT_DIR, 'logs')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, FAULTS_DIR, MODELS_DIR, PLOTTING_DIR])

# Import modules
try:
    from Basilisk.utilities import macros, vizSupport
    from faults.rw_fault import RWFaultScenario, run as run_scenario
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class TargetDefinition:
    """Class to hold target location definitions"""
    def __init__(self, name, latitude, longitude, color="red"):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.color = color
        
    @staticmethod
    def from_string(target_str):
        """Parse a target string in format 'name:lat:lon[:color]'"""
        parts = target_str.split(':')
        if len(parts) < 3:
            raise ValueError("Target format should be 'name:latitude:longitude[:color]'")
        
        name = parts[0]
        try:
            lat = float(parts[1])
            lon = float(parts[2])
        except ValueError:
            raise ValueError("Latitude and longitude must be numbers")
            
        color = parts[3] if len(parts) > 3 else "red"
        return TargetDefinition(name, lat, lon, color)
    
    def to_dict(self):
        """Convert to dictionary format for the scenario"""
        return {
            "name": self.name,
            "lat": self.latitude,
            "lon": self.longitude,
            "color": self.color
        }

class SimulationConfig:
    """Configuration class for simulation parameters"""
    def __init__(self):
        self.simulation_time = 30.0  # minutes
        self.fault_magnitude = 0.0005
        self.fault_wheel_number = 3  # 0-indexed wheel number
        self.fault_time = 10.0  # minutes
        self.enable_periodic_fault = False
        self.periodic_fault_interval = 360  # seconds
        self.periodic_fault_magnitude = 0.1
        self.periodic_fault_wheel = 1  # 0-indexed wheel number
        self.binary_filename = "rw_fault_viz"
        self.targets = [
            TargetDefinition("Melbourne", -37.8136, 144.9631, "red"),
            TargetDefinition("New York", 40.71, -74.00, "blue"),
            TargetDefinition("Tokyo", 35.68, 139.77, "green"),
            TargetDefinition("London", 51.51, -0.13, "yellow")
        ]
        self.camera_position = [0.0, 2.0, 0.0]
        self.show_plots = True
        self.save_binary = True
        self.interactive = False
        
    def validate(self):
        """Validate configuration parameters"""
        if self.fault_wheel_number not in range(4):
            raise ValueError("Fault wheel number must be between 0 and 3")
        if self.periodic_fault_wheel not in range(4):
            raise ValueError("Periodic fault wheel number must be between 0 and 3")
        if self.fault_time >= self.simulation_time:
            raise ValueError("Fault time must be less than simulation time")
        if self.fault_magnitude <= 0:
            raise ValueError("Fault magnitude must be positive")
        if self.enable_periodic_fault and self.periodic_fault_magnitude <= 0:
            raise ValueError("Periodic fault magnitude must be positive")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Spacecraft Reaction Wheel Fault Injection CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples:
          # Run with default parameters
          python spacecraft_fault_cli.py
          
          # Run with custom fault parameters
          python spacecraft_fault_cli.py --fault-magnitude 0.001 --fault-wheel 2 --fault-time 15.0
          
          # Enable periodic fault injection
          python spacecraft_fault_cli.py --enable-periodic-fault --periodic-fault-interval 180
          
          # Add custom target locations
          python spacecraft_fault_cli.py --add-target "Sydney:-33.87:151.21:purple" --add-target "Cairo:30.04:31.24:orange"
          
          # Run in interactive mode
          python spacecraft_fault_cli.py -i
        """)
    )
    
    # Simulation parameters
    parser.add_argument("--sim-time", type=float, help="Simulation time in minutes (default: 30.0)")
    parser.add_argument("--no-plots", action="store_true", help="Disable plotting of results")
    parser.add_argument("--no-binary", action="store_true", help="Disable binary file generation")
    parser.add_argument("-o", "--output", help="Binary output filename (default: rw_fault_viz)")
    
    # Fault parameters
    parser.add_argument("--fault-magnitude", type=float, help="Magnitude of the friction fault (default: 0.0005)")
    parser.add_argument("--fault-wheel", type=int, choices=range(4), 
                        help="Wheel number to inject fault into (0-3, default: 3)")
    parser.add_argument("--fault-time", type=float, help="Time to inject fault in minutes (default: 10.0)")
    
    # Periodic fault parameters
    parser.add_argument("--enable-periodic-fault", action="store_true", help="Enable periodic fault injection")
    parser.add_argument("--periodic-fault-interval", type=float, help="Interval between periodic faults in seconds (default: 360)")
    parser.add_argument("--periodic-fault-magnitude", type=float, help="Magnitude of periodic faults (default: 0.1)")
    parser.add_argument("--periodic-fault-wheel", type=int, choices=range(4), 
                        help="Wheel number for periodic faults (0-3, default: 1)")
    
    # Target parameters
    parser.add_argument("--add-target", action="append", help="Add target location in format 'name:lat:lon[:color]'")
    parser.add_argument("--clear-targets", action="store_true", help="Clear default target locations")
    
    # Camera parameters
    parser.add_argument("--camera-position", type=str, help="Camera position in body frame as 'x,y,z' (default: 0,2,0)")
    
    # Mode
    parser.add_argument("-i", "--interactive", action="store_true", help="Run in interactive mode")
    
    return parser.parse_args()

def apply_args_to_config(args, config):
    """Apply parsed arguments to configuration"""
    # Simulation parameters
    if args.sim_time is not None:
        config.simulation_time = args.sim_time
    if args.no_plots:
        config.show_plots = False
    if args.no_binary:
        config.save_binary = False
    if args.output:
        config.binary_filename = args.output
    
    # Fault parameters
    if args.fault_magnitude is not None:
        config.fault_magnitude = args.fault_magnitude
    if args.fault_wheel is not None:
        config.fault_wheel_number = args.fault_wheel
    if args.fault_time is not None:
        config.fault_time = args.fault_time
    
    # Periodic fault parameters
    if args.enable_periodic_fault:
        config.enable_periodic_fault = True
    if args.periodic_fault_interval is not None:
        config.periodic_fault_interval = args.periodic_fault_interval
    if args.periodic_fault_magnitude is not None:
        config.periodic_fault_magnitude = args.periodic_fault_magnitude
    if args.periodic_fault_wheel is not None:
        config.periodic_fault_wheel = args.periodic_fault_wheel
    
    # Target parameters
    if args.clear_targets:
        config.targets = []
    if args.add_target:
        for target_str in args.add_target:
            config.targets.append(TargetDefinition.from_string(target_str))
    
    # Camera parameters
    if args.camera_position:
        try:
            x, y, z = map(float, args.camera_position.split(','))
            config.camera_position = [x, y, z]
        except (ValueError, TypeError):
            print("ERROR: Camera position must be in format 'x,y,z' with numeric values")
            sys.exit(1)
    
    # Mode
    config.interactive = args.interactive

def interactive_mode(config):
    """Run interactive mode to configure simulation parameters"""
    print("\n===== Spacecraft RW Fault Simulation Interactive Mode =====")
    
    # Simulation time
    try:
        sim_time = float(input(f"Simulation time in minutes [{config.simulation_time}]: ") or config.simulation_time)
        config.simulation_time = sim_time
    except ValueError:
        print("Invalid input, using default value")
    
    # Fault parameters
    try:
        config.fault_magnitude = float(input(f"Fault magnitude [{config.fault_magnitude}]: ") or config.fault_magnitude)
    except ValueError:
        print("Invalid input, using default value")
    
    try:
        config.fault_wheel_number = int(input(f"Fault wheel number (0-3) [{config.fault_wheel_number}]: ") or config.fault_wheel_number)
        if config.fault_wheel_number not in range(4):
            raise ValueError("Wheel number must be between 0 and 3")
    except ValueError as e:
        print(f"Invalid input: {e}, using default value")
    
    try:
        config.fault_time = float(input(f"Fault time in minutes [{config.fault_time}]: ") or config.fault_time)
    except ValueError:
        print("Invalid input, using default value")
    
    # Periodic fault
    enable_periodic = input(f"Enable periodic fault (y/n) [{'y' if config.enable_periodic_fault else 'n'}]: ") or ('y' if config.enable_periodic_fault else 'n')
    config.enable_periodic_fault = enable_periodic.lower() == 'y'
    
    if config.enable_periodic_fault:
        try:
            config.periodic_fault_interval = float(input(f"Periodic fault interval in seconds [{config.periodic_fault_interval}]: ") or config.periodic_fault_interval)
        except ValueError:
            print("Invalid input, using default value")
        
        try:
            config.periodic_fault_magnitude = float(input(f"Periodic fault magnitude [{config.periodic_fault_magnitude}]: ") or config.periodic_fault_magnitude)
        except ValueError:
            print("Invalid input, using default value")
        
        try:
            config.periodic_fault_wheel = int(input(f"Periodic fault wheel (0-3) [{config.periodic_fault_wheel}]: ") or config.periodic_fault_wheel)
            if config.periodic_fault_wheel not in range(4):
                raise ValueError("Wheel number must be between 0 and 3")
        except ValueError as e:
            print(f"Invalid input: {e}, using default value")
    
    # Target management
    print("\nCurrent targets:")
    for i, target in enumerate(config.targets):
        print(f"  {i}. {target.name} (lat: {target.latitude}, lon: {target.longitude}, color: {target.color})")
    
    modify_targets = input("Modify targets? (y/n) [n]: ") or 'n'
    if modify_targets.lower() == 'y':
        clear_targets = input("Clear existing targets? (y/n) [n]: ") or 'n'
        if clear_targets.lower() == 'y':
            config.targets = []
        
        add_targets = input("Add targets? (y/n) [y]: ") or 'y'
        if add_targets.lower() == 'y':
            while True:
                target_str = input("Enter target as 'name:lat:lon[:color]' (empty to finish): ")
                if not target_str:
                    break
                try:
                    config.targets.append(TargetDefinition.from_string(target_str))
                except ValueError as e:
                    print(f"Invalid target format: {e}")
    
    # Camera position
    camera_pos_str = input(f"Camera position as 'x,y,z' [{','.join(map(str, config.camera_position))}]: ") or ','.join(map(str, config.camera_position))
    try:
        x, y, z = map(float, camera_pos_str.split(','))
        config.camera_position = [x, y, z]
    except (ValueError, TypeError):
        print("Invalid input, using default value")
    
    # Output options
    config.binary_filename = input(f"Binary output filename [{config.binary_filename}]: ") or config.binary_filename
    
    show_plots = input(f"Show plots (y/n) [{'y' if config.show_plots else 'n'}]: ") or ('y' if config.show_plots else 'n')
    config.show_plots = show_plots.lower() == 'y'
    
    save_binary = input(f"Save binary file (y/n) [{'y' if config.save_binary else 'n'}]: ") or ('y' if config.save_binary else 'n')
    config.save_binary = save_binary.lower() == 'y'

def create_custom_scenario(config):
    """Create and configure a custom RW fault scenario based on the configuration"""
    # Create the base scenario
    scenario = RWFaultScenario()
    
    # Set simulation time
    simulationTime = macros.min2nano(config.simulation_time)
    
    # Set camera position
    scenario.cameraLocation = config.camera_position
    
    # Set targets
    scenario.targets = [target.to_dict() for target in config.targets]
    
    # Set fault parameters
    scenario.oneTimeFaultTime = macros.min2nano(config.fault_time)
    
    # Override the fault event creation
    def customize_scenario(scenario):
        """Customize the scenario with our specific fault configuration"""
        # Set mode request
        scenario.modeRequest = "hillPoint"
        
        # Configure one-time fault event
        scenario.createNewEvent(
            "addOneTimeRWFault",
            scenario.get_FswModel().processTasksTimeStep,
            True,
            ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
            [f"self.get_DynModel().AddRWFault('friction',{config.fault_magnitude},{config.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
             "self.oneTimeRWFaultFlag=0"]
        )
        
        # Configure periodic fault event if enabled
        if config.enable_periodic_fault:
            scenario.createNewEvent(
                "addRepeatedRWFault",
                scenario.get_FswModel().processTasksTimeStep,
                True,
                ["self.repeatRWFaultFlag==1"],
                [f"self.get_DynModel().PeriodicRWFault({config.periodic_fault_interval},'friction',{config.periodic_fault_magnitude},{config.periodic_fault_wheel}, self.TotalSim.CurrentNanos)", 
                 "self.setEventActivity('addRepeatedRWFault',True)"]
            )
        
        # Set up visualization if vizSupport is available
        viz = None
        if vizSupport.vizFound and config.save_binary:
            # Create visualization with RW effector list
            binary_path = os.path.join(VIZ_DIR, config.binary_filename)
            viz = vizSupport.enableUnityVisualization(
                scenario,
                scenario.get_DynModel().taskName,
                scenario.get_DynModel().scObject,
                rwEffectorList=scenario.get_DynModel().rwStateEffector,
                liveStream=False,
                saveFile=binary_path  # <- now uses absolute path in Vizfile
            )

            
            # Add target locations
            for target in scenario.targets:
                lat = target["lat"]
                lon = target["lon"]
                color = target.get("color", "red")
                
                # Convert lat/lon to ECEF coordinates
                alt = 0.0
                radius = 6371000.0 + alt  # Earth radius + altitude
                lat_rad = lat * macros.D2R
                lon_rad = lon * macros.D2R
                x = radius * np.cos(lat_rad) * np.cos(lon_rad)
                y = radius * np.cos(lat_rad) * np.sin(lon_rad)
                z = radius * np.sin(lat_rad)
                location_position = [x, y, z]
                
                # Add location to visualization
                vizSupport.addLocation(
                    viz,
                    stationName=target["name"],
                    parentBodyName="earth",
                    r_GP_P=location_position,
                    color=color
                )
            
            # Add camera that looks at the spacecraft body
            vizSupport.createStandardCamera(
                viz,
                setMode=1,  # Standard camera mode
                spacecraftName=scenario.get_DynModel().scObject.ModelTag,
                fieldOfView=70 * macros.D2R,
                displayName="RW Camera",
                pointingVector_B=[0, 0, 0],  # Look at spacecraft center
                position_B=scenario.cameraLocation  # Camera position in body frame
            )
        
        # Initialize and run the simulation
        scenario.InitializeSimulation()
        scenario.ConfigureStopTime(simulationTime)
        scenario.ExecuteSimulation()
        
        return viz
    
    return scenario, customize_scenario

def run_custom_simulation(config):
    """Run a customized simulation based on the configuration"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(LOGS_DIR, f"sim_results_{timestamp}")
    
    print("\n===== Running Custom RW Fault Scenario =====")
    print(f"Configuration:")
    print(f"  - Simulation time: {config.simulation_time} minutes")
    print(f"  - Fault magnitude: {config.fault_magnitude}")
    print(f"  - Fault wheel: {config.fault_wheel_number}")
    print(f"  - Fault time: {config.fault_time} minutes")
    
    if config.enable_periodic_fault:
        print(f"  - Periodic fault enabled:")
        print(f"    - Interval: {config.periodic_fault_interval} seconds")
        print(f"    - Magnitude: {config.periodic_fault_magnitude}")
        print(f"    - Wheel: {config.periodic_fault_wheel}")
    
    print(f"  - Save binary: {config.save_binary}")
    if config.save_binary:
        print(f"  - Binary filename: {config.binary_filename}")
    
    print(f"  - Show plots: {config.show_plots}")
    print(f"  - Targets: {len(config.targets)} defined")
    
    # Create and run the custom scenario
    scenario, customize_scenario = create_custom_scenario(config)
    viz = customize_scenario(scenario)
    
    # Generate and display/save plots
    figureList = scenario.pull_outputs(config.show_plots)
    
    # Save configuration summary
    try:
        os.makedirs(output_dir, exist_ok=True)
        with open(f"{output_dir}/simulation_summary.txt", "w") as f:
            f.write(f"===== RW Fault Simulation Summary =====\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Simulation time: {config.simulation_time} minutes\n")
            f.write(f"Fault magnitude: {config.fault_magnitude}\n")
            f.write(f"Fault wheel: {config.fault_wheel_number}\n")
            f.write(f"Fault time: {config.fault_time} minutes\n\n")
            
            if config.enable_periodic_fault:
                f.write(f"Periodic fault enabled:\n")
                f.write(f"  Interval: {config.periodic_fault_interval} seconds\n")
                f.write(f"  Magnitude: {config.periodic_fault_magnitude}\n")
                f.write(f"  Wheel: {config.periodic_fault_wheel}\n\n")
            
            f.write(f"Targets:\n")
            for i, target in enumerate(config.targets):
                f.write(f"  {i+1}. {target.name} (lat: {target.latitude}, lon: {target.longitude}, color: {target.color})\n")
            
            f.write(f"\nResults saved to: {output_dir}\n")
            
            if config.save_binary and viz:
                    vizard_path = os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin")
                    if os.path.exists(vizard_path):
                        file_size = os.path.getsize(vizard_path) / (1024*1024)  # Size in MB
                        print(f"\nBinary file created: {vizard_path} ({file_size:.2f} MB)")
                        print("\nTo view this file in Vizard:")
                        print("1. Open Vizard application")
                        print(f"2. Load the binary file: {vizard_path}")
                    else:
                        print(f"\nWarning: Expected binary file {vizard_path} was not created")
        
        print(f"\nSimulation summary saved to {output_dir}/simulation_summary.txt")
    except Exception as e:
        print(f"Warning: Could not save simulation summary: {e}")
    
    if config.save_binary and viz:
        vizard_path = f"./_VizFiles/{config.binary_filename}_UnityViz.bin"
        if os.path.exists(vizard_path):
            file_size = os.path.getsize(vizard_path) / (1024*1024)  # Size in MB
            print(f"\nBinary file created: {vizard_path} ({file_size:.2f} MB)")
            print("\nTo view this file in Vizard:")
            print("1. Open Vizard application")
            print(f"2. Load the binary file: {vizard_path}")
        else:
            print(f"\nWarning: Expected binary file {vizard_path} was not created")
    
    return scenario, viz, figureList, output_dir

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_args()
    
    # Create default configuration
    config = SimulationConfig()
    
    # Apply command line arguments to config
    apply_args_to_config(args, config)
    
    # Run interactive mode if requested
    if config.interactive:
        interactive_mode(config)
    
    try:
        # Validate configuration
        config.validate()
        
        # Run the simulation
        scenario, viz, figureList, output_dir = run_custom_simulation(config)
        
        print("\n===== Simulation Complete =====")
        return 0
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())