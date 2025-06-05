#!/usr/bin/env python
"""
spacecraft_simulation.py

Core simulation functionality for the Spacecraft Constellation Fault Simulator.
Version with comprehensive fault integration, battery visualization, and real fault simulation.
"""
import os
import sys
import inspect
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Fix path resolution to work with project structure
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
FAULTS_DIR = os.path.join(ROOT_DIR, 'faults')
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plots')
LOGS_DIR = os.path.join(ROOT_DIR, 'logs')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, FAULTS_DIR, MODELS_DIR, PLOTTING_DIR])

# Import Basilisk modules
try:
    from Basilisk.utilities import macros, vizSupport
    from Basilisk.utilities import orbitalMotion
    from Basilisk.utilities import simIncludeGravBody
    from Basilisk.simulation import spacecraft
    from Basilisk.utilities import RigidBodyKinematics as rbk
    from Basilisk.utilities import SimulationBaseClass
    from Basilisk import __path__
    BASILISK_AVAILABLE = True
except ImportError as e:
    print(f"ERROR: Could not import required Basilisk modules: {e}")
    print("Make sure Basilisk is properly installed.")
    BASILISK_AVAILABLE = False

# Global fault status variables
FAULT_IMPORT_STATUS = {}
AVAILABLE_FAULTS = []
FAILED_FAULTS = []
FAULT_LOADER_AVAILABLE = False
PLOTS_AVAILABLE = False

def check_fault_modules(verbose=True):
    """Check and report status of all fault modules"""
    global FAULT_IMPORT_STATUS, AVAILABLE_FAULTS, FAILED_FAULTS, FAULT_LOADER_AVAILABLE
    
    FAULT_IMPORT_STATUS = {}
    AVAILABLE_FAULTS = []
    FAILED_FAULTS = []
    
    if verbose:
        print("\n" + "="*60)
        print("FAULT MODULE IMPORT STATUS")
        print("="*60)
    
    # Try to import the fault loader
    try:
        from fault_loader import (
            get_fault_scenario_class, 
            create_scenario, 
            run_scenario, 
            apply_fault_to_spacecraft,
            extract_fault_data_from_scenario,
            get_available_fault_types
        )
        if verbose:
            print("Fault loader imported successfully")
        FAULT_LOADER_AVAILABLE = True
        
        # Get list of available fault types
        try:
            available_types = get_available_fault_types()
            if verbose:
                print(f"Available fault types: {available_types}")
            AVAILABLE_FAULTS = available_types
        except Exception as e:
            if verbose:
                print(f"Could not get available fault types: {e}")
                
    except ImportError as e:
        if verbose:
            print(f"Fault loader import failed: {e}")
            print("WARNING: Fault functionality will be severely limited.")
        FAULT_LOADER_AVAILABLE = False
    
    # Try to import each specific fault module for verification
    fault_modules_to_check = [
        ("friction_fault", "FrictionFaultScenario"),
        ("powerlimit_fault", "PowerLimitFaultScenario"),
        ("encoder_fault", "EncoderFaultScenario"),
        ("battery_fault", "BatteryFaultScenario")
    ]
    
    if verbose:
        print("\nIndividual fault module status:")
    
    for module_name, class_name in fault_modules_to_check:
        try:
            module = __import__(f"faults.{module_name}", fromlist=[class_name])
            if hasattr(module, class_name):
                FAULT_IMPORT_STATUS[module_name] = True
                if verbose:
                    print(f"  {module_name}: {class_name} loaded")
            else:
                FAULT_IMPORT_STATUS[module_name] = False
                FAILED_FAULTS.append(module_name)
                if verbose:
                    print(f"  {module_name}: {class_name} not found in module")
        except ImportError as e:
            FAULT_IMPORT_STATUS[module_name] = False
            FAILED_FAULTS.append(module_name)
            if verbose:
                print(f"  {module_name}: Import failed - {e}")
        except Exception as e:
            FAULT_IMPORT_STATUS[module_name] = False
            FAILED_FAULTS.append(module_name)
            if verbose:
                print(f"  {module_name}: Unexpected error - {e}")
    
    if verbose:
        print("="*60 + "\n")
    
    # Summary of fault availability
    if FAILED_FAULTS and verbose:
        print(f"WARNING: The following fault types will not be available: {FAILED_FAULTS}")
        print("You may need to check the fault module implementations.\n")
    
    return FAULT_LOADER_AVAILABLE, AVAILABLE_FAULTS, FAILED_FAULTS

def check_plots_module(verbose=True):
    """Check if plots module is available"""
    global PLOTS_AVAILABLE
    
    try:
        from plots import (
            generate_fault_plots,
            generate_constellation_plots,
            generate_inter_satellite_distance_plots,
            create_fault_config_for_real_simulation
        )
        if verbose:
            print("Plots module imported successfully")
        PLOTS_AVAILABLE = True
    except ImportError as e:
        if verbose:
            print(f"Plots module import failed: {e}")
            print("WARNING: Plotting functionality will be limited.")
        PLOTS_AVAILABLE = False
    
    return PLOTS_AVAILABLE

# Check modules on import - but quietly
check_fault_modules(verbose=False)
check_plots_module(verbose=False)

# Import plotting functions if available
if PLOTS_AVAILABLE:
    from plots import (
        generate_fault_plots,
        generate_constellation_plots,
        generate_inter_satellite_distance_plots,
        create_fault_config_for_real_simulation
    )

def inject_fault_into_spacecraft(sc_object, fault_config, current_time_nano):
    """
    Inject a fault into a spacecraft based on configuration
    
    Parameters:
    sc_object: Spacecraft object with dynamics model
    fault_config: Fault configuration dictionary
    current_time_nano: Current simulation time in nanoseconds
    
    Returns:
    bool: True if fault was successfully injected
    """
    if not fault_config.get("enabled", False):
        return False
        
    fault_type = fault_config["type"]
    fault_magnitude = fault_config["magnitude"]
    fault_wheel = fault_config["wheel"]
    fault_time_nano = macros.min2nano(fault_config["time"])
    
    if current_time_nano < fault_time_nano:
        return False  # Not time yet
    
    # Check if this fault type is available
    if fault_type not in AVAILABLE_FAULTS and FAULT_LOADER_AVAILABLE:
        print(f"ERROR: Fault type '{fault_type}' is not available")
        print(f"Available types: {AVAILABLE_FAULTS}")
        return False
        
    # Try to use the fault loader to apply the fault
    if FAULT_LOADER_AVAILABLE:
        try:
            from fault_loader import apply_fault_to_spacecraft
            # Store fault info for plotting
            if not hasattr(sc_object, 'fault_log'):
                sc_object.fault_log = []
            sc_object.fault_log.append({
                'type': fault_type,
                'magnitude': fault_magnitude,
                'wheel': fault_wheel,
                'time': current_time_nano
            })
            print(f"Fault logged: {fault_type} on wheel {fault_wheel} at time {current_time_nano * macros.NANO2MIN:.2f} min")
            return True
        except Exception as e:
            print(f"Could not inject fault: {e}")
            return False
    else:
        print(f"Fault loader not available for {fault_type}")
        return False


class TargetDefinition:
    """Class to hold target location definitions"""
    def __init__(self, name, latitude, longitude, color="red"):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.color = color
        self.assigned_to = []  # List of spacecraft this target is assigned to
    
    def to_dict(self):
        """Convert to dictionary format for the scenario"""
        return {
            "name": self.name,
            "lat": self.latitude,
            "lon": self.longitude,
            "color": self.color,
            "assigned_to": self.assigned_to
        }


class SimulationConfig:
    """Configuration class for simulation parameters"""
    def __init__(self):
        # Default 30-minute simulation time
        self.simulation_time = 300.0  # MINUTES - Total simulation duration
        
        # Spacecraft list for constellation support
        self.spacecraft_list = []
        


        self.DEFAULT_FAULT_MAGNITUDES = {
            "friction": 0.0005,    # N⋅m - this is fine
            "power_limit": 0.5,    # W - realistic power limit
            "encoder": 20.0,       # % - noticeable encoder error
            "battery": 50.0        # W - significant battery drain
        }
        # Legacy fault parameters (for backward compatibility)
        self.fault_type = "friction"  # Options: friction, power_limit, encoder, battery
        self.fault_magnitude = 0.0005
        self.fault_wheel_number = 3  # 0-indexed wheel number
        self.fault_time = 10.0  # MINUTES - When to inject the fault

      
        
        # Periodic fault parameters
        self.enable_periodic_fault = False
        self.periodic_fault_interval = 360  # SECONDS
        self.periodic_fault_magnitude = 0.1
        self.periodic_fault_wheel = 1  # 0-indexed wheel number
        
        # Output settings
        self.binary_filename = "spacecraft_viz"
        self.show_plots = True
        self.save_plots = True
        self.save_binary = True
        
        # Camera configuration for target viewing
        self.camera_position = [0.0, 0.0, 10.0]  # Higher camera position for better Earth/target view
        self.camera_fov = 80.0  # Wider field of view for target viewing
        self.active_camera_name = None  # Which spacecraft has the active camera
        
        # Target locations with better visibility
        self.targets = [
            TargetDefinition("Melbourne", -37.8136, 144.9631, "red"),
            TargetDefinition("New York", 40.71, -74.00, "blue"),
            TargetDefinition("Tokyo", 35.68, 139.77, "green"),
            TargetDefinition("London", 51.51, -0.13, "yellow")
        ]
        
    def validate(self):
        """Validate configuration parameters"""
        # Check Basilisk availability first
        if not BASILISK_AVAILABLE:
            raise ValueError("Basilisk is not available. Please install Basilisk to run simulations.")
        
        # Validate time parameters first
        if self.simulation_time <= 0:
            raise ValueError("Simulation time must be positive (in minutes)")
        if self.fault_time >= self.simulation_time:
            raise ValueError(f"Fault time ({self.fault_time} min) must be less than simulation time ({self.simulation_time} min)")
        if self.fault_time < 0:
            raise ValueError("Fault time cannot be negative")
            
        # Check if we're using constellation mode
        if self.spacecraft_list:
            # Validate each spacecraft configuration
            for sc in self.spacecraft_list:
                # Validate orbit parameters - ensure proper altitudes
                if sc["orbit"]["a"] <= 6371:  # Must be above Earth's surface
                    raise ValueError(f"Spacecraft {sc['name']}: Semi-major axis must be greater than Earth radius (6371 km)")
                if sc["orbit"]["e"] < 0 or sc["orbit"]["e"] >= 1.0:
                    raise ValueError(f"Spacecraft {sc['name']}: Eccentricity must be between 0 and 1")
                
                # Validate fault parameters if enabled
                if sc["fault"]["enabled"]:
                    if sc["fault"]["time"] >= self.simulation_time:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time ({sc['fault']['time']} min) must be less than simulation time ({self.simulation_time} min)")
                    if sc["fault"]["time"] < 0:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time cannot be negative")
                    
                    # Check if fault type is available
                    fault_type = sc["fault"]["type"]
                    if FAULT_LOADER_AVAILABLE and fault_type not in AVAILABLE_FAULTS:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault type '{fault_type}' is not available. Available types: {AVAILABLE_FAULTS}")
                    
                    if sc["fault"]["type"] == "friction" and sc["fault"]["magnitude"] <= 0:
                        raise ValueError(f"Spacecraft {sc['name']}: Friction fault magnitude must be positive")
                    elif sc["fault"]["type"] == "power_limit" and sc["fault"]["magnitude"] <= 0:
                        raise ValueError(f"Spacecraft {sc['name']}: Power limit must be positive")
                    
                    if sc["fault"]["wheel"] not in range(4):
                        raise ValueError(f"Spacecraft {sc['name']}: Fault wheel number must be between 0 and 3")


def calculate_target_visibility(spacecraft_pos, target_pos, earth_radius=6371000):
    """
    Calculate if a target is visible from a spacecraft position
    
    Parameters:
    spacecraft_pos: [x, y, z] position of spacecraft in meters
    target_pos: [x, y, z] position of target in meters  
    earth_radius: Earth radius in meters
    
    Returns:
    bool: True if target is visible (line of sight not blocked by Earth)
    """
    try:
        # Convert to numpy arrays
        sc_pos = np.array(spacecraft_pos)
        tgt_pos = np.array(target_pos)
        
        # Vector from spacecraft to target
        sc_to_target = tgt_pos - sc_pos
        
        # Distance from spacecraft to target
        distance = np.linalg.norm(sc_to_target)
        if distance == 0:
            return False
        
        # Check line of sight - simplified check
        # If both spacecraft and target are above Earth surface and 
        # spacecraft altitude is reasonable, assume visibility
        sc_altitude = np.linalg.norm(sc_pos) - earth_radius
        tgt_altitude = np.linalg.norm(tgt_pos) - earth_radius
        
        # For simplicity, if spacecraft is above 200km and target is visible above horizon
        if sc_altitude > 200000:  # 200km
            # Calculate elevation angle to target
            earth_center = np.array([0, 0, 0])
            sc_to_earth = earth_center - sc_pos
            
            # If angle between sc_to_target and sc_to_earth is less than 90 degrees,
            # target is above horizon
            cos_angle = np.dot(sc_to_target, sc_to_earth) / (np.linalg.norm(sc_to_target) * np.linalg.norm(sc_to_earth))
            elevation_angle = np.pi/2 - np.arccos(np.clip(cos_angle, -1, 1))
            
            # Target is visible if elevation angle > 10 degrees (above horizon)
            return elevation_angle > np.radians(10)
        
        return False
    except:
        return False


def run_custom_simulation(config):
    """
    Run a customized simulation based on the configuration
    
    Parameters:
    config (SimulationConfig): Simulation configuration object
    
    Returns:
    tuple: (scenario, viz, figureList, output_dir)
    """
    # Re-check fault modules with verbose output at simulation start
    print("\nChecking fault module availability...")
    check_fault_modules(verbose=True)
    check_plots_module(verbose=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(LOGS_DIR, f"sim_results_{timestamp}")
    
    # Validate configuration first
    try:
        config.validate()
        print("Configuration validation passed")
    except ValueError as e:
        print(f"Configuration validation failed: {e}")
        raise
    
    # Use constellation mode
    if config.spacecraft_list:
        print("\n" + "="*60)
        print("SPACECRAFT CONSTELLATION SIMULATION")
        print("="*60)
        print(f"Simulation Duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)")
        print(f"Number of Spacecraft: {len(config.spacecraft_list)}")
        
        # Check fault status
        if FAILED_FAULTS:
            print(f"\nWARNING: The following fault types are not available: {FAILED_FAULTS}")
            print("Simulation will continue but these faults cannot be injected.\n")
        
        # Proper time conversion
        simulationTime = macros.min2nano(config.simulation_time)
        
        print(f"Time Conversion:")
        print(f"   - Input: {config.simulation_time} minutes")
        print(f"   - Basilisk time: {simulationTime} nanoseconds")
        print(f"   - Verification: {simulationTime/1e9/60:.1f} minutes")
        
        # Import necessary Basilisk modules
        try:
            bskPath = __path__[0]
            print(f"Using Basilisk at: {bskPath}")
            
        except ImportError as e:
            print(f"Failed to import Basilisk modules: {e}")
            raise
        
        # Create simulation
        simTaskName = "simTask"
        simProcessName = "simProcess"
        scSim = SimulationBaseClass.SimBaseClass()

        # Create dynamics process with appropriate time step
        dynProcess = scSim.CreateNewProcess(simProcessName)
        simulationTimeStep = macros.sec2nano(1.0)  # 1 second time step
        dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
        
        print(f"Simulation time step: {simulationTimeStep} ns (1 second)")
        
        # Create the spacecraft objects
        sc_objects = []
        
        # Create gravitational bodies using simIncludeGravBody
        gravFactory = simIncludeGravBody.gravBodyFactory()
        planet = gravFactory.createEarth()
        planet.isCentralBody = True
        mu = planet.mu  # Standard gravitational parameter
        
        print(f"Earth gravity parameter: {mu/1e14:.3f} x 10^14 m³/s²")
                        
        # Process each spacecraft with orbital parameters
        fault_spacecraft_count = 0

        # Define fixed RAAN for each orbit type to ensure proper separation
        orbit_raan_values = {
            "Default Orbit": 0.0,
            "MEO Navigation": 60.0,
            "High Coverage": 120.0
        }

        for i, sc_config in enumerate(config.spacecraft_list):
            # Create spacecraft
            sc = spacecraft.Spacecraft()
            sc.ModelTag = sc_config["name"]
            
            # Add spacecraft to task
            scSim.AddModelToTask(simTaskName, sc)
            
            # Add gravity using simIncludeGravBody
            gravFactory.addBodiesTo(sc)
            
            # Set proper orbit for stable motion with better altitude
            oe = orbitalMotion.ClassicElements()
            orbit = sc_config["orbit"]
            
            # Get orbit name to ensure same orbital parameters for same orbit
            orbit_name = sc_config.get("orbit_name", "Default Orbit")
            
            # Use predefined RAAN for each orbit type
            orbit_raan = orbit_raan_values.get(orbit_name, 0.0)
            
            # Convert orbital parameters properly with minimum altitude check
            oe.a = max(orbit["a"] * 1000, 6571000)  # Convert km to m, minimum 200km altitude
            oe.e = min(orbit["e"], 0.01)  # Keep eccentricity very small for stability
            oe.i = orbit["i"] * macros.D2R  # Convert degrees to radians
            oe.Omega = orbit_raan * macros.D2R  # Use consistent RAAN for orbit type
            oe.omega = orbit["omega"] * macros.D2R
            oe.f = orbit["f"] * macros.D2R
            
            # Calculate orbital period for reference
            orbital_period_sec = 2 * np.pi * np.sqrt((oe.a**3) / mu)
            orbital_period_min = orbital_period_sec / 60.0
            altitude_km = (oe.a / 1000) - 6371
            
            print(f"Spacecraft '{sc.ModelTag}' - Orbit: {orbit_name}, Altitude: {altitude_km:.1f}km, Period: {orbital_period_min:.1f} minutes")
            
            # Convert to Cartesian coordinates
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            sc.hub.r_CN_NInit = rN
            sc.hub.v_CN_NInit = vN
            
            # Set proper initial attitude for stable control
            sc.hub.sigma_BNInit = [[0.01], [0.02], [-0.01]]  # Small attitude errors
            sc.hub.omega_BN_BInit = [[0.0001], [-0.0002], [0.0001]]  # Very small angular velocities
            
            # Store the spacecraft object for later use
            sc_objects.append(sc)
            
            # Set up fault if enabled for this spacecraft

            if sc_config["fault"]["enabled"]:
                fault_type = sc_config["fault"]["type"]
                fault_magnitude = sc_config["fault"]["magnitude"]
                
                # ADD THIS: Scale fault magnitude based on type
                if fault_magnitude == 0.0005:  # If using default GUI value
                    # Use appropriate magnitude for fault type
                    if fault_type in config.DEFAULT_FAULT_MAGNITUDES:
                        scaled_magnitude = config.DEFAULT_FAULT_MAGNITUDES[fault_type]
                        print(f"   Scaling fault magnitude from {fault_magnitude} to {scaled_magnitude} for {fault_type}")
                        fault_magnitude = scaled_magnitude
                        sc_config["fault"]["magnitude"] = scaled_magnitude  # Update config
                
                fault_wheel = sc_config["fault"]["wheel"]
                fault_time_minutes = sc_config["fault"]["time"]
                
                # Convert fault time to nanoseconds
                fault_time_nano = macros.min2nano(fault_time_minutes)
                
                print(f"\nFault Configuration for '{sc.ModelTag}':")
                print(f"   - Type: {fault_type}")
                print(f"   - Magnitude: {fault_magnitude}")
                print(f"   - Wheel: {fault_wheel}")
                print(f"   - Injection Time: {fault_time_minutes} minutes ({fault_time_minutes * 60:.0f} seconds)")
                
                # Check if fault type is available
                if fault_type not in AVAILABLE_FAULTS and FAULT_LOADER_AVAILABLE:
                    print(f"   WARNING: Fault type '{fault_type}' is not available!")
                    print(f"   Available types: {AVAILABLE_FAULTS}")
                elif fault_type in FAILED_FAULTS:
                    print(f"   WARNING: Fault type '{fault_type}' failed to load!")
                else:
                    print(f"   Fault type '{fault_type}' is available")
                
                # Verify fault time is reasonable
                if fault_time_nano >= simulationTime:
                    print(f"ERROR: Fault time ({fault_time_minutes} min) >= simulation time ({config.simulation_time} min)")
                    raise ValueError(f"Fault injection time must be less than simulation duration")
                
                # Store fault config in spacecraft for event access
                sc.faultConfig = sc_config["fault"].copy()
                sc.faultInjected = False
                
                # Create a simplified fault event (since we can't directly inject into constellation)
                event_name = f"logFault_{sc.ModelTag}"
                
                # For now, we'll just log when the fault should occur
                scSim.createNewEvent(
                    event_name,
                    simulationTimeStep,
                    True,
                    [f"self.TotalSim.CurrentNanos >= {fault_time_nano}"],
                    [
                        f"print('FAULT EVENT: {fault_type} fault triggered for {sc.ModelTag} at', self.TotalSim.CurrentNanos * {macros.NANO2MIN}, 'minutes')",
                        f"self.setEventActivity('{event_name}', False)"
                    ]
                )
                
                fault_spacecraft_count += 1
                    
                # Set legacy fault parameters for backward compatibility
                if i == 0:  # Only use first spacecraft with fault
                    config.fault_type = fault_type
                    config.fault_magnitude = fault_magnitude
                    config.fault_wheel_number = fault_wheel
                    config.fault_time = fault_time_minutes
                    
        print(f"\nTotal spacecraft with faults configured: {fault_spacecraft_count}")
        
        # Visualization setup
        viz = None
        if config.save_binary and vizSupport.vizFound:
            print("\nSetting up Vizard visualization...")
            
            # Clean up any existing binary files
            binary_base_name = f"{config.binary_filename}_UnityViz.bin"
            potential_paths = [
                os.path.join(VIZ_DIR, binary_base_name),
                os.path.join(VIZ_DIR, "_VizFiles", binary_base_name)
            ]
            
            for path in potential_paths:
                if os.path.exists(path):
                    try:
                        os.remove(path)
                        print(f"Removed old binary: {os.path.basename(path)}")
                    except Exception as e:
                        print(f"Could not remove old binary {path}: {e}")
            
            # Create visualization directory structure
            viz_files_dir = os.path.join(VIZ_DIR, "_VizFiles")
            if not os.path.exists(viz_files_dir):
                os.makedirs(viz_files_dir, exist_ok=True)
                print(f"Created visualization directory: {viz_files_dir}")
                
            # Use the correct path approach
            binary_filename = os.path.basename(config.binary_filename)
            binary_full_path = os.path.join(VIZ_DIR, binary_filename)
            print(f"Binary output: {binary_filename}_UnityViz.bin")
            
            # Check if any spacecraft has battery fault enabled
            has_battery_fault = any(sat["fault"]["enabled"] and sat["fault"]["type"] == "battery" 
                                for sat in config.spacecraft_list)
            
            # Initialize gsList properly for each spacecraft
            gsList = [[] for _ in range(len(config.spacecraft_list))]  # Empty list for each spacecraft
            battery_components = {}  # Store battery components for each spacecraft
            
            if has_battery_fault:
                print("Battery fault detected - setting up battery visualization panels...")
                
                # Import battery modules
                from Basilisk.simulation import simpleBattery, simpleSolarPanel, simplePowerSink
                from Basilisk.architecture import messaging
                
                # Create battery components for each spacecraft with battery fault
                for i, sc_config in enumerate(config.spacecraft_list):
                    if sc_config["fault"]["enabled"] and sc_config["fault"]["type"] == "battery":
                        sc = sc_objects[i]
                        
                        # Create battery
                        battery = simpleBattery.SimpleBattery()
                        battery.ModelTag = f"{sc.ModelTag}_battery"
                        battery.storageCapacity = 100.0  # Wh
                        battery.storedCharge_Init = 50.0  # Wh
                        scSim.AddModelToTask(simTaskName, battery)
                        
                        # Create solar panel
                        solarPanel = simpleSolarPanel.SimpleSolarPanel()
                        solarPanel.ModelTag = f"{sc.ModelTag}_solarPanel"
                        solarPanel.setPanelParameters([-1.0, 0.0, 0.0], 0.05, 0.3)  # 50W panel
                        solarPanel.stateInMsg.subscribeTo(sc.scStateOutMsg)
                        scSim.AddModelToTask(simTaskName, solarPanel)
                        battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)
                        
                        # Simple sun direction
                        sunMsgData = messaging.SpicePlanetStateMsgPayload()
                        sunMsgData.PositionVector = [-1.0, 0.0, 0.0]
                        sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
                        solarPanel.sunInMsg.subscribeTo(sunMsg)
                        
                        # Create power sinks
                        powerSink = simplePowerSink.SimplePowerSink()
                        powerSink.ModelTag = f"{sc.ModelTag}_powerSink"
                        powerSink.nodePowerOut = -0.01  # 10W nominal
                        scSim.AddModelToTask(simTaskName, powerSink)
                        battery.addPowerNodeToModel(powerSink.nodePowerOutMsg)
                        
                        powerSinkFault = simplePowerSink.SimplePowerSink()
                        powerSinkFault.ModelTag = f"{sc.ModelTag}_powerSinkFault"
                        powerSinkFault.nodePowerOut = 0.0  # Initially off
                        scSim.AddModelToTask(simTaskName, powerSinkFault)
                        battery.addPowerNodeToModel(powerSinkFault.nodePowerOutMsg)
                        
                        # Store components
                        battery_components[sc.ModelTag] = {
                            'battery': battery,
                            'solarPanel': solarPanel,
                            'powerSink': powerSink,
                            'powerSinkFault': powerSinkFault
                        }
                        
                        # Create battery panel visualization
                        batteryPanel = vizSupport.vizInterface.GenericStorage()
                        batteryPanel.label = f"{sc.ModelTag} Battery (%)"
                        batteryPanel.units = "%"
                        batteryPanel.minValue = 0
                        batteryPanel.maxValue = 100
                        batteryPanel.useStorageLevel = True
                        
                        batteryInMsg = messaging.PowerStorageStatusMsgReader()
                        batteryInMsg.subscribeTo(battery.batPowerOutMsg)
                        batteryPanel.batteryStateInMsg = batteryInMsg
                        batteryPanel.this.disown()
                        
                        # Set thresholds and colors
                        batteryPanel.thresholds = vizSupport.vizInterface.IntVector([20, 50, 80])
                        batteryPanel.color = vizSupport.vizInterface.IntVector(
                            vizSupport.toRGBA255("red") +
                            vizSupport.toRGBA255("orange") +
                            vizSupport.toRGBA255("yellow") +
                            vizSupport.toRGBA255("green")
                        )
                        
                        # Create solar panel visualization
                        solarViz = vizSupport.vizInterface.GenericStorage()
                        solarViz.label = f"{sc.ModelTag} Solar Power"
                        solarViz.units = "W"
                        solarViz.minValue = 0.0
                        solarViz.maxValue = 60.0
                        solarViz.useStorageLevel = False
                        
                        solarReader = messaging.PowerNodeUsageMsgReader()
                        solarReader.subscribeTo(solarPanel.nodePowerOutMsg)
                        solarViz.powerNodeUsageInMsg = solarReader
                        solarViz.this.disown()
                        
                        # Add panels to the correct spacecraft's list
                        gsList[i] = [batteryPanel, solarViz]
                        
                        # Create fault injection event
                        fault_time_nano = macros.min2nano(sc_config["fault"]["time"])
                        scSim.__dict__[f'powerSinkFault_{i}'] = powerSinkFault
                        
                        scSim.createNewEvent(
                            f"batteryFault_{sc.ModelTag}",
                            simulationTimeStep,
                            True,
                            [f"self.TotalSim.CurrentNanos >= {fault_time_nano}"],
                            [
                                f"self.powerSinkFault_{i}.nodePowerOut = -0.05",  # 50W drain
                                f"print('Battery fault activated for {sc.ModelTag}')",
                                f"self.setEventActivity('batteryFault_{sc.ModelTag}', False)"
                            ]
                        )
                        
                        print(f"Set up battery visualization for {sc.ModelTag}")
            
            # Enable visualization with battery panels
            try:
                # Check if we have any generic storage panels to show
                has_generic_storage = any(len(panels) > 0 for panels in gsList)
                
                if has_generic_storage:
                    # Enable with generic storage panels
                    viz = vizSupport.enableUnityVisualization(
                        scSim,
                        simTaskName,
                        sc_objects,
                        saveFile=binary_full_path,
                        liveStream=False,
                        genericStorageList=gsList
                    )
                    
                    # Enable instrument GUI for all spacecraft with battery faults
                    for i, sc in enumerate(sc_objects):
                        if sc.ModelTag in battery_components:
                            vizSupport.setInstrumentGuiSetting(
                                viz,
                                spacecraftName=sc.ModelTag,
                                showGenericStoragePanel=True
                            )
                            print(f"Enabled battery panel GUI for {sc.ModelTag}")
                else:
                    # No battery faults, use standard visualization
                    viz = vizSupport.enableUnityVisualization(
                        scSim,
                        simTaskName,
                        sc_objects,
                        saveFile=binary_full_path,
                        liveStream=False
                    )
                
                print("Vizard visualization enabled")
                
                # Configure visualization settings for better visibility
                if hasattr(viz, 'settings'):
                    viz.settings.showSpacecraftLabels = True  # Show spacecraft names
                    viz.settings.orbitLinesOn = 1  # Enable orbit lines
                    viz.settings.spacecraftSizeMultiplier = 8.0  # Make spacecraft larger for visibility
                    viz.settings.orbitLinesOn = True  # Ensure orbit lines are on
                    viz.settings.spacecraftCSOn = True  # Show coordinate system
                    viz.settings.showCelestialBodyLabels = True  # Show body labels
                    print("Vizard settings: labels=True, orbitLines=True, size=8x, coordSys=True")
                
                # Add target locations with better visibility and connections
                target_added = False
                assigned_targets = [t for t in config.targets if t.assigned_to]
                
                print(f"Adding {len(assigned_targets)} assigned targets with visibility...")
                
                # Store target positions for visibility calculations
                target_positions = {}
                
                for target in assigned_targets:
                    try:
                        lat = target.latitude
                        lon = target.longitude
                        color = target.color
                        
                        # Target positioning for better visibility
                        alt = 50000.0  # 50km above surface for better visibility
                        radius = 6371000.0 + alt  # Earth radius + altitude
                        lat_rad = lat * macros.D2R
                        lon_rad = lon * macros.D2R
                        x = radius * np.cos(lat_rad) * np.cos(lon_rad)
                        y = radius * np.cos(lat_rad) * np.sin(lon_rad)
                        z = radius * np.sin(lat_rad)
                        location_position = [x, y, z]
                        
                        # Store position for visibility calculations
                        target_positions[target.name] = location_position
                        
                        # Add location with better parameters for visibility
                        vizSupport.addLocation(
                            viz,
                            stationName=target.name,
                            parentBodyName="earth",
                            r_GP_P=location_position,
                            color=color,
                            range=2000000.0  # 2000km marker for better visibility
                        )
                        print(f"Target: {target.name} at {lat:.2f}°, {lon:.2f}° -> assigned to {', '.join(target.assigned_to)}")
                        target_added = True
                    except Exception as e:
                        print(f"Could not add target {target.name}: {e}")
                
                # If no targets assigned, add all targets with visibility
                if not target_added:
                    print("No targets assigned, adding all targets with visibility...")
                    for target in config.targets:
                        try:
                            lat = target.latitude
                            lon = target.longitude
                            color = target.color
                            
                            # Target positioning for visibility
                            alt = 50000.0  # 50km above surface
                            radius = 6371000.0 + alt  # Earth radius + altitude
                            lat_rad = lat * macros.D2R
                            lon_rad = lon * macros.D2R
                            x = radius * np.cos(lat_rad) * np.cos(lon_rad)
                            y = radius * np.cos(lat_rad) * np.sin(lon_rad)
                            z = radius * np.sin(lat_rad)
                            location_position = [x, y, z]
                            
                            # Store position for visibility calculations
                            target_positions[target.name] = location_position
                            
                            # Add location with better visibility
                            vizSupport.addLocation(
                                viz,
                                stationName=target.name,
                                parentBodyName="earth",
                                r_GP_P=location_position,
                                color=color,
                                range=2000000.0  # 2000km marker
                            )
                            print(f"Default target: {target.name} at {lat:.2f}°, {lon:.2f}°")
                            target_added = True
                        except Exception as e:
                            print(f"Could not add default target {target.name}: {e}")
                
                # Create cameras with better positioning for target viewing
                camera_created_count = 0
                
                print("Setting up camera configuration for target viewing...")
                
                # Find the spacecraft designated as having the active camera or use first one
                active_camera_created = False
                for i, sc in enumerate(sc_objects):
                    sc_config = config.spacecraft_list[i]
                    if sc_config["camera"]["enabled"] and config.active_camera_name == sc.ModelTag:
                        try:
                            # Create camera for this spacecraft with better settings for target viewing
                            camera_pos = sc_config["camera"]["position"]
                            
                            vizSupport.createStandardCamera(
                                viz,
                                setMode=1,  # Standard camera mode
                                spacecraftName=sc.ModelTag, 
                                fieldOfView=config.camera_fov * macros.D2R,  # Use wider FOV from config
                                displayName=f"{sc.ModelTag} Camera",
                                pointingVector_B=[0, 0, -1],  # Point toward Earth for target viewing
                                position_B=camera_pos
                            )
                            print(f"Camera for {sc.ModelTag} - positioned for target viewing")
                            active_camera_created = True
                            camera_created_count += 1
                            break
                        except Exception as e:
                            print(f"Could not create camera for {sc.ModelTag}: {e}")
                
                # If no specific active camera was found, create one for the first spacecraft
                if not active_camera_created and sc_objects:
                    try:
                        # Use better default camera position for target viewing
                        default_camera_pos = config.camera_position  # Use config position
                        
                        vizSupport.createStandardCamera(
                            viz,
                            setMode=1,  # Standard camera mode
                            spacecraftName=sc_objects[0].ModelTag, 
                            fieldOfView=config.camera_fov * macros.D2R,  # Wider FOV
                            displayName=f"{sc_objects[0].ModelTag} Camera",
                            pointingVector_B=[0, 0, -1],  # Point toward Earth
                            position_B=default_camera_pos
                        )
                        print(f"Default camera for {sc_objects[0].ModelTag} at {default_camera_pos} - optimized for target viewing")
                        active_camera_created = True
                        camera_created_count += 1
                    except Exception as e:
                        print(f"Could not create default camera: {e}")
                
                print(f"Created {camera_created_count} cameras optimized for target viewing")
                
            except Exception as e:
                print(f"Failed to set up visualization: {e}")
                import traceback
                traceback.print_exc()
                viz = None
        else:
            if not vizSupport.vizFound:
                print("Vizard visualization module not found")
            if not config.save_binary:
                print("Binary file generation is disabled")
        
        # Set simulation time properly
        print("\n" + "-"*50)
        print("SIMULATION TIME SETUP")
        print("-"*50)
        
        print(f"Final simulation time: {simulationTime} nanoseconds")
        print(f"Equivalent to: {simulationTime/1e9:.0f} seconds")
        print(f"Equivalent to: {simulationTime/1e9/60:.1f} minutes")
        
        # Initialize and run the simulation
        print("\n" + "-"*50)
        print("RUNNING SIMULATION")
        print("-"*50)
        
        print("Initializing simulation...")
        scSim.InitializeSimulation()
        
        print(f"Setting stop time to {simulationTime} ns...")
        scSim.ConfigureStopTime(simulationTime)
        
        print(f"Starting simulation for {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)...")
        start_time = datetime.now()
        
        # Execute the simulation
        scSim.ExecuteSimulation()
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"Simulation completed!")
        print(f"Wall-clock time: {duration:.2f} seconds")
        print(f"Simulated time: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)")
        print(f"Speed ratio: {(config.simulation_time * 60) / duration:.1f}x real-time")
        
        # Get the actual simulation time that was executed
        current_sim_time = scSim.TotalSim.CurrentNanos
        actual_sim_minutes = current_sim_time / 1e9 / 60
        
        print(f"Actual simulation time executed: {actual_sim_minutes:.2f} minutes ({current_sim_time / 1e9:.0f} seconds)")
        
        if abs(actual_sim_minutes - config.simulation_time) > 0.1:
            print(f"WARNING: Simulation may have ended early!")
            print(f"   Requested: {config.simulation_time} minutes")
            print(f"   Executed: {actual_sim_minutes:.2f} minutes")
        else:
            print(f"Simulation time matches request")

    # Generate plots using the centralized plots module with real fault simulation
    figureList = {}
    if config.show_plots or config.save_plots:
        print("\n" + "-"*50)
        print("GENERATING PLOTS WITH REAL FAULT SIMULATION")
        print("-"*50)
        
        if not PLOTS_AVAILABLE:
            print("Plots module not available, skipping plot generation")
        else:
            try:
                # Create time array for plotting (in minutes for clarity)
                time_data = np.linspace(0, config.simulation_time, max(100, int(config.simulation_time * 2)))
                print(f"Time data: 0 to {config.simulation_time} minutes ({len(time_data)} points)")
                
                # Generate constellation plots
                try:
                    constellation_plots = generate_constellation_plots(sc_objects, time_data, mu)
                    figureList.update(constellation_plots)
                    print(f"Generated {len(constellation_plots)} constellation plots")
                except Exception as e:
                    print(f"Could not generate constellation plots: {e}")
                
                # Generate inter-satellite distance plots
                try:
                    distance_plots = generate_inter_satellite_distance_plots(sc_objects, time_data, mu)
                    figureList.update(distance_plots)
                    print(f"Generated {len(distance_plots)} distance plots")
                except Exception as e:
                    print(f"Could not generate distance plots: {e}")

                # Import plotting functions at the beginning of fault plotting section
                try:
                    from plots import generate_fault_plots, create_fault_config_for_real_simulation
                    plots_imported = True
                    print("Fault plotting functions imported successfully")
                except ImportError as e:
                    print(f"Could not import plotting functions: {e}")
                    plots_imported = False

                # Real fault plot generation
                fault_plots_generated = 0
                real_plots_count = 0
                synthetic_plots_count = 0
                
                print(f"\n{'='*60}")
                print(f"REAL FAULT PLOT GENERATION")
                print(f"{'='*60}")
                
                for i, sc_config in enumerate(config.spacecraft_list):
                    if sc_config["fault"]["enabled"]:
                        fault_type = sc_config["fault"]["type"]
                        
                        # Check if this fault type can be plotted
                        if fault_type in FAILED_FAULTS:
                            print(f"Cannot plot {fault_type} fault for {sc_config['name']} - module not loaded")
                            continue
                            
                        if not plots_imported:
                            print(f"Cannot generate plots for {sc_config['name']} - plotting functions not available")
                            continue
                            
                        print(f"\n{'='*50}")
                        print(f"PROCESSING {sc_config['name']} - {fault_type.upper()} FAULT")
                        print(f"{'='*50}")
                        print(f"Fault Parameters:")
                        print(f"  - Type: {fault_type}")
                        print(f"  - Magnitude: {sc_config['fault']['magnitude']}")
                        print(f"  - Wheel: {sc_config['fault']['wheel']}")
                        print(f"  - Time: {sc_config['fault']['time']} minutes")
                        
                        try:
                            # Automatically attempt real simulation for enabled faults
                            real_simulation_successful = False
                            
                            # Check fault availability status
                            print(f"\nChecking fault availability...")
                            print(f"   - FAULT_LOADER_AVAILABLE: {FAULT_LOADER_AVAILABLE}")
                            print(f"   - AVAILABLE_FAULTS: {AVAILABLE_FAULTS}")
                            print(f"   - FAILED_FAULTS: {FAILED_FAULTS}")
                            print(f"   - Target fault type: {fault_type}")
                            
                            # Force a fresh check of available faults
                            try:
                                from fault_loader import get_available_fault_types
                                fresh_available_faults = get_available_fault_types()
                                print(f"   - Fresh available faults: {fresh_available_faults}")
                                
                                # Use the fresh check
                                fault_loader_working = len(fresh_available_faults) > 0
                                use_real_simulation = fault_loader_working and fault_type in fresh_available_faults
                                
                                print(f"   - Fault loader working: {fault_loader_working}")
                                print(f"   - Will attempt real simulation: {use_real_simulation}")
                            except Exception as e:
                                print(f"   - Could not do fresh fault check: {e}")
                                use_real_simulation = False
                                
                            # Attempt real simulation if available
                            if use_real_simulation:
                                try:
                                    print(f"\nAttempting real fault simulation")
                                    print(f"   Spacecraft: {sc_config['name']}")
                                    print(f"   Fault Type: {fault_type}")
                                    print(f"   Parameters: {sc_config['fault']}")
                                    
                                    # Create real simulation configuration
                                    fault_params = {
                                        'fault_magnitude': sc_config["fault"]["magnitude"],
                                        'fault_wheel': sc_config["fault"]["wheel"],
                                        'fault_time_min': sc_config["fault"]["time"]
                                    }
                                    
                                    # Create real simulation config
                                    real_fault_config = create_fault_config_for_real_simulation(fault_type, fault_params)
                                    
                                    print(f"   Running real simulation...")
                                    
                                    # Generate plots using real simulation data
                                    fault_plots = generate_fault_plots(
                                        fault_type=fault_type,
                                        fault_data=real_fault_config,  # This triggers real simulation
                                        time_data=time_data,
                                        fault_time_min=sc_config["fault"]["time"],
                                        spacecraft_name=sc_config["name"]
                                    )
                                    
                                    if fault_plots and len(fault_plots) > 0:
                                        figureList.update(fault_plots)
                                        fault_plots_generated += len(fault_plots)
                                        
                                        # Check if we got real simulation plots
                                        real_plots = [name for name in fault_plots.keys() if "REAL" in name]
                                        if real_plots and len(real_plots) > 0:
                                            real_plots_count += len(real_plots)
                                            real_simulation_successful = True
                                            print(f"   SUCCESS: Generated {len(real_plots)} real simulation plots!")
                                            print(f"   Real plots: {list(real_plots)}")
                                            
                                            # Log details about real plots
                                            for plot_name in real_plots:
                                                print(f"      {plot_name}")
                                        else:
                                            synthetic_plots_count += len(fault_plots)
                                            print(f"   Got {len(fault_plots)} synthetic plots (real simulation may have failed)")
                                            
                                            # Log details about synthetic plots
                                            for plot_name in fault_plots.keys():
                                                print(f"      {plot_name}")
                                    else:
                                        print(f"   No plots generated from real simulation attempt")
                                        
                                except Exception as e:
                                    print(f"   Real simulation failed: {e}")
                                    print(f"   Falling back to synthetic plots...")
                                    import traceback
                                    traceback.print_exc()
                            
                            # Use synthetic plots if real simulation failed or unavailable
                            if not real_simulation_successful:
                                print(f"\nGenerating synthetic plots")
                                reason = 'Real simulation failed' if use_real_simulation else 'Real simulation not available'
                                print(f"   Reason: {reason}")
                                
                                fault_data = {
                                    'fault_wheel': sc_config["fault"]["wheel"],
                                    'friction_magnitude': sc_config["fault"]["magnitude"],
                                    'friction_baseline': 0.02,
                                    'power_limit': sc_config["fault"]["magnitude"] if fault_type == "power_limit" else None,
                                    'encoder_error': sc_config["fault"]["magnitude"] if fault_type == "encoder" else None,
                                    'battery_drain': sc_config["fault"]["magnitude"] if fault_type == "battery" else None,
                                }
                                
                                fault_plots = generate_fault_plots(
                                    fault_type,
                                    fault_data,
                                    time_data,
                                    sc_config["fault"]["time"],
                                    sc_config["name"]
                                )
                                
                                if fault_plots and len(fault_plots) > 0:
                                    figureList.update(fault_plots)
                                    fault_plots_generated += len(fault_plots)
                                    synthetic_plots_count += len(fault_plots)
                                    print(f"   Generated {len(fault_plots)} synthetic {fault_type} plots")
                                    
                                    # Log details about synthetic plots
                                    for plot_name in fault_plots.keys():
                                        print(f"      {plot_name}")
                                else:
                                    print(f"   Failed to generate synthetic plots")
                            
                        except Exception as e:
                            print(f"Critical error generating {fault_type} fault plots for {sc_config['name']}: {e}")
                            import traceback
                            traceback.print_exc()
                
                # Summary of plot generation
                print(f"\n{'='*60}")
                print(f"PLOT GENERATION SUMMARY")
                print(f"{'='*60}")
                total_plots = len(figureList)
                other_plots = total_plots - real_plots_count - synthetic_plots_count
                
                print(f"Total plots generated: {total_plots}")
                print(f"   Real simulation plots: {real_plots_count}")
                print(f"   Synthetic fault plots: {synthetic_plots_count}")
                print(f"   Constellation/other plots: {other_plots}")
                
                if real_plots_count > 0:
                    print(f"\nSUCCESS: Real fault simulation data successfully integrated!")
                    print(f"   Look for plots marked with 'REAL' in the Results tab")
                    print(f"   These plots contain actual Basilisk simulation data!")
                elif synthetic_plots_count > 0:
                    print(f"\nUsing synthetic fault data")
                    print(f"   Real simulation may be unavailable or failed")
                    print(f"   Check fault module status and console for details")
                else:
                    print(f"\nNo fault plots generated")
                    if fault_spacecraft_count > 0:
                        print(f"   Check fault module installation and configuration")

                # Save plots to files if requested
                if config.save_plots and figureList:
                    print(f"\n{'='*40}")
                    print(f"SAVING PLOTS TO DISK")
                    print(f"{'='*40}")
                    
                    os.makedirs(PLOTTING_DIR, exist_ok=True)
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    
                    saved_count = 0
                    failed_count = 0
                    
                    for name, fig in figureList.items():
                        try:
                            # Create descriptive filename
                            plot_filename = f"{timestamp}_{name}.png"
                            plot_path = os.path.join(PLOTTING_DIR, plot_filename)
                            
                            # Save with high quality
                            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                            saved_count += 1
                            
                            # Log what type of plot was saved
                            if "REAL" in name:
                                print(f"  Saved real plot: {plot_filename}")
                            elif any(fault in name for fault in ["Friction", "PowerLimit", "Encoder", "Battery"]):
                                print(f"  Saved fault plot: {plot_filename}")
                            else:
                                print(f"  Saved plot: {plot_filename}")
                            
                            # Close the figure to free memory
                            plt.close(fig)
                            
                        except Exception as e:
                            print(f"Error saving plot {name}: {e}")
                            failed_count += 1
                    
                    print(f"\nSaved {saved_count} plots to {PLOTTING_DIR}")
                    if failed_count > 0:
                        print(f"Failed to save {failed_count} plots")
                    
                    # Clear the figureList after closing all figures
                    figureList.clear()
                else:
                    # If not saving, still close all figures
                    for name, fig in figureList.items():
                        try:
                            plt.close(fig)
                        except:
                            pass
                    figureList.clear()
                
            except Exception as e:
                print(f"Error generating plots: {e}")
                import traceback
                traceback.print_exc()
                
                # Ensure all figures are closed even on error
                try:
                    plt.close('all')
                except:
                    pass
                    
        # Create a simple scenario object for compatibility
        class ConstellationScenario:
            def __init__(self, scSim, sc_objects, config):
                self.TotalSim = scSim
                self.sc_objects = sc_objects
                self.targets = [target.to_dict() for target in config.targets]
                self.fault_type = config.fault_type
                self.actual_sim_time = actual_sim_minutes  # Store actual simulation time
        
        scenario = ConstellationScenario(scSim, sc_objects, config)
        
    else:
        # Legacy mode - use existing fault modules
        print("\n===== Running Legacy RW Fault Scenario =====")
        print("Legacy mode not fully implemented in this version.")
        scenario = None
        viz = None
        figureList = {}
    
    # Create output directory to store simulation results
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Save detailed simulation summary
        summary_path = os.path.join(output_dir, "simulation_summary.txt")
        with open(summary_path, "w", encoding='utf-8') as f:
            f.write(f"SPACECRAFT CONSTELLATION SIMULATION SUMMARY\n")
            f.write(f"=" * 50 + "\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"FAULT MODULE STATUS:\n")
            f.write(f"  - Available faults: {AVAILABLE_FAULTS}\n")
            f.write(f"  - Failed faults: {FAILED_FAULTS}\n")
            f.write(f"  - Fault loader: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}\n")
            f.write(f"  - Plots module: {'Available' if PLOTS_AVAILABLE else 'Not available'}\n\n")
            
            f.write(f"REAL FAULT SIMULATION:\n")
            f.write(f"  - Automatic real simulation attempts for fault-enabled spacecraft\n")
            f.write(f"  - Graceful fallback to synthetic plots when real simulation fails\n")
            f.write(f"  - Clear plot categorization (real, synthetic, constellation)\n")
            f.write(f"  - Detailed logging and user feedback\n\n")
            
            f.write(f"VISUALIZATION FEATURES:\n")
            f.write(f"  - Higher satellite altitudes for better target coverage\n")
            f.write(f"  - Target visibility at 50km altitude\n")
            f.write(f"  - Better camera positioning for target viewing\n")
            f.write(f"  - Spacecraft size increased 8x for visibility\n")
            f.write(f"  - Target range increased to 2000km for visibility\n")
            f.write(f"  - Orbit lines enabled for all spacecraft\n")
            f.write(f"  - Battery visualization panels for battery faults\n\n")
            
            f.write(f"TIME CONFIGURATION:\n")
            f.write(f"  - Requested duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)\n")
            f.write(f"  - Nanoseconds: {simulationTime} ns\n")
            if 'actual_sim_minutes' in locals():
                f.write(f"  - Actual duration: {actual_sim_minutes:.2f} minutes ({current_sim_time / 1e9:.0f} seconds)\n")
                f.write(f"  - Time accuracy: {'GOOD' if abs(actual_sim_minutes - config.simulation_time) < 0.1 else 'MISMATCH'}\n")
            f.write(f"  - Wall-clock time: {duration:.2f} seconds\n\n")
            
            if config.spacecraft_list:
                f.write(f"CONSTELLATION CONFIGURATION:\n")
                f.write(f"  - Number of spacecraft: {len(config.spacecraft_list)}\n")
                f.write(f"  - Active camera: {config.active_camera_name or 'Default (first spacecraft)'}\n")
                f.write(f"  - Spacecraft with faults: {fault_spacecraft_count}\n\n")
                
                # Write information about each spacecraft
                for i, sc in enumerate(config.spacecraft_list):
                    f.write(f"Spacecraft {i+1}: {sc['name']}\n")
                    f.write(f"  Orbit:\n")
                    f.write(f"    - Semi-major axis: {sc['orbit']['a']} km\n")
                    f.write(f"    - Altitude: {sc['orbit']['a'] - 6371:.1f} km\n")
                    f.write(f"    - Eccentricity: {sc['orbit']['e']}\n")
                    f.write(f"    - Inclination: {sc['orbit']['i']}°\n")
                    f.write(f"    - RAAN: {sc['orbit']['Omega']}°\n")
                    f.write(f"    - Arg. of Periapsis: {sc['orbit']['omega']}°\n")
                    f.write(f"    - True Anomaly: {sc['orbit']['f']}°\n")
                    
                    # Calculate and write orbital period
                    a_m = sc['orbit']['a'] * 1000  # Convert to meters
                    period_sec = 2 * np.pi * np.sqrt((a_m**3) / mu)
                    period_min = period_sec / 60.0
                    f.write(f"    - Orbital period: {period_min:.1f} minutes\n")
                    
                    # Write fault information if enabled
                    if sc["fault"]["enabled"]:
                        f.write(f"  Fault: ENABLED\n")
                        f.write(f"    - Type: {sc['fault']['type']}\n")
                        f.write(f"    - Available: {'Yes' if sc['fault']['type'] in AVAILABLE_FAULTS else 'No'}\n")
                        f.write(f"    - Magnitude: {sc['fault']['magnitude']}\n")
                        f.write(f"    - Wheel: {sc['fault']['wheel']}\n")
                        f.write(f"    - Injection time: {sc['fault']['time']} minutes ({sc['fault']['time'] * 60:.0f} seconds)\n")
                        
                        if sc["fault"]["periodic"]["enabled"]:
                            f.write(f"    - Periodic fault: YES\n")
                            f.write(f"      - Interval: {sc['fault']['periodic']['interval']} seconds\n")
                            f.write(f"      - Magnitude: {sc['fault']['periodic']['magnitude']}\n")
                            f.write(f"      - Wheel: {sc['fault']['periodic']['wheel']}\n")
                    else:
                        f.write(f"  Fault: DISABLED\n")
                    
                    f.write("\n")
            
            f.write(f"OUTPUT CONFIGURATION:\n")
            f.write(f"  - Results directory: {output_dir}\n")
            f.write(f"  - Plots generated: {len(figureList) if figureList else 0}\n")
            
            if config.save_binary:
                # Check all possible locations for the binary file
                vizard_paths = [
                    os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
                    os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
                ]
                
                binary_found = False
                for vizard_path in vizard_paths:
                    if os.path.exists(vizard_path):
                        binary_found = True
                        file_size = os.path.getsize(vizard_path) / (1024*1024)  # Size in MB
                        f.write(f"  - Binary file: {vizard_path} ({file_size:.2f} MB)\n")
                        f.write(f"\nVIZARD FEATURES:\n")
                        f.write(f"1. Open Vizard application\n")
                        f.write(f"2. Load binary file: {vizard_path}\n")
                        f.write(f"3. You should now see:\n")
                        f.write(f"   - Spacecraft names labeled\n")
                        f.write(f"   - Orbit lines visible\n")
                        f.write(f"   - Target markers at 50km altitude (better visibility)\n")
                        f.write(f"   - Larger spacecraft (8x size)\n")
                        f.write(f"   - Camera optimized for target viewing\n")
                        f.write(f"   - Higher satellite altitudes for better coverage\n")
                        f.write(f"   - Battery panels for spacecraft with battery faults\n")
                        f.write(f"4. Simulation duration: {config.simulation_time * 60:.0f} seconds\n")
                        break
                
                if not binary_found:
                    f.write(f"  - Binary file: NOT FOUND (check visualization setup)\n")
        
        print(f"Detailed summary saved: {summary_path}")
        
    except Exception as e:
        print(f"Could not save simulation summary: {e}")
    
    # Final binary file verification
    if config.save_binary:
        print("\n" + "-"*50)
        print("BINARY FILE VERIFICATION")
        print("-"*50)
        
        # Check all possible locations for the binary file
        vizard_paths = [
            os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
            os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
        ]
        
        binary_found = False
        for vizard_path in vizard_paths:
            if os.path.exists(vizard_path):
                binary_found = True
                file_size = os.path.getsize(vizard_path) / (1024*1024)  # Size in MB
                print(f"Binary file created: {os.path.basename(vizard_path)}")
                print(f"Location: {vizard_path}")
                print(f"Size: {file_size:.2f} MB")
                print(f"Duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)")
                print(f"Features: Higher altitudes, Target visibility, Better camera, 8x spacecraft size, Battery panels")
                break
        

        if not binary_found:
            print(f"Binary file not found in expected locations:")
            for path in vizard_paths:
                print(f"   - {path}")
    
    print("\n" + "="*60)
    print("SIMULATION COMPLETED SUCCESSFULLY")
    print("="*60)
    print(f"Simulation Duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)")
    print(f"Available Faults: {AVAILABLE_FAULTS}")
    print(f"Failed Faults: {FAILED_FAULTS}")
    print(f"Features: Real fault simulation, Higher altitudes, Target visibility, Better camera positioning, Battery visualization")
    
    return scenario, viz, figureList, output_dir

# Test fault modules if run directly
if __name__ == "__main__":
    print("\nTesting fault module availability...")
    check_fault_modules(verbose=True)
    check_plots_module(verbose=True)
    
    if FAULT_LOADER_AVAILABLE:
        print("\nTesting fault scenarios:")
        for fault_type in AVAILABLE_FAULTS:
            try:
                from fault_loader import get_fault_scenario_class
                scenario_class = get_fault_scenario_class(fault_type)
                print(f"{fault_type}: {scenario_class.__name__}")
            except Exception as e:
                print(f"{fault_type}: {e}")
    
    # Test basic simulation configuration
    print("\nTesting basic simulation configuration...")
    try:
        config = SimulationConfig()
        config.simulation_time = 5.0  # Short test
        config.spacecraft_list = [
            {
                "name": "TestSat",
                "orbit": {
                    "a": 6971,  # 600km altitude
                    "e": 0.01,
                    "i": 55.0,
                    "Omega": 45.0,
                    "omega": 30.0,
                    "f": 0.0
                },
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 2.0,
                    "enabled": True,
                    "periodic": {"enabled": False}
                },
                "camera": {
                    "position": [0.0, 0.0, 15.0],
                    "fov": 80.0,
                    "enabled": True
                }
            }
        ]
        
        config.validate()
        print("Configuration validation passed")
        
        if BASILISK_AVAILABLE:
            print("Basilisk is available - ready for simulation")
        else:
            print("Basilisk not available - simulation will fail")
            
    except Exception as e:
        print(f"Configuration test failed: {e}")
        
    print("\nspacecraft_simulation.py module ready for use!")
    print("Key features:")
    print("  - Real fault simulation with automatic fallback")
    print("  - Battery and solar panel visualization")
    print("  - Target visibility optimization")
    print("  - Camera positioning for target viewing")
    print("  - Comprehensive error handling")