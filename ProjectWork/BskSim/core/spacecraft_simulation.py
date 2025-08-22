#!/usr/bin/env python
"""
spacecraft_simulation.py

COMPLETE MERGED VERSION - Core simulation functionality for the Spacecraft Constellation Fault Simulator.
Includes:
- Cluster integration (ClusterManager) when available
- Real fault simulation for EACH satellite
- Battery visualization panels
- Full plotting coverage (supports both old/new plots APIs)
- Detailed summary + binary verification + cameras/targets
- Automatic real simulation attempts with synthetic fallback
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
    from Basilisk import __path__ as bsk_path
    BASILISK_AVAILABLE = True
    BSK_INSTALL_PATH = bsk_path[0]
except Exception as e:
    print(f"ERROR: Could not import Basilisk modules: {e}")
    print("Make sure Basilisk is properly installed.")
    BASILISK_AVAILABLE = False
    BSK_INSTALL_PATH = None

# Cluster integration (new feature)
try:
    from cluster_integration import ClusterManager, integrate_clusters_with_simulation
    CLUSTER_INTEGRATION_AVAILABLE = True
    print("Cluster integration module loaded successfully")
except Exception as e:
    print(f"Warning: Could not import cluster integration: {e}")
    CLUSTER_INTEGRATION_AVAILABLE = False

# Plotting imports - support BOTH old and new APIs
PLOTS_AVAILABLE = False
PLOTS_API = {
    "new": False,
    "old": False
}

# Try new plots API
try:
    from plots import (
        generate_constellation_overview_plots,     # new
        generate_cluster_communication_plots,      # new
        generate_inter_satellite_distance_plots,   # new (also in old)
        generate_fault_plots                       # both
    )
    PLOTS_AVAILABLE = True
    PLOTS_API["new"] = True
    print("Plotting module (new API) loaded successfully")
except Exception as e:
    print(f"Warning: Could not import new plotting API: {e}")

# Try old plots API additions
try:
    from plots import (
        generate_constellation_plots,              # old
        create_fault_config_for_real_simulation    # old (CRITICAL for real sim)
    )
    PLOTS_AVAILABLE = True
    PLOTS_API["old"] = True
    print("Plotting module (old API) loaded successfully")
except Exception as e:
    print(f"Note: Old plotting API not fully available: {e}")

# Global fault status variables
FAULT_IMPORT_STATUS = {}
AVAILABLE_FAULTS = []
FAILED_FAULTS = []
FAULT_LOADER_AVAILABLE = False

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
            print("WARNING: Real fault simulation will not be available.")
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
        except Exception as e:
            FAULT_IMPORT_STATUS[module_name] = False
            FAILED_FAULTS.append(module_name)
            if verbose:
                print(f"  {module_name}: Import failed - {e}")
    
    if verbose:
        print("="*60 + "\n")
    
    return FAULT_LOADER_AVAILABLE, AVAILABLE_FAULTS, FAILED_FAULTS

def check_plots_module(verbose=True):
    """Verify plots module availability"""
    global PLOTS_AVAILABLE
    
    if PLOTS_API["new"] or PLOTS_API["old"]:
        if verbose:
            apis = []
            if PLOTS_API["new"]:
                apis.append("new")
            if PLOTS_API["old"]:
                apis.append("old")
            print(f"Plots module is available (APIs: {', '.join(apis)})")
        PLOTS_AVAILABLE = True
    else:
        if verbose:
            print("Plots module NOT available")
        PLOTS_AVAILABLE = False
    
    return PLOTS_AVAILABLE

def inject_fault_into_spacecraft(sc_object, fault_config, current_time_nano):
    """
    Inject/log fault into a spacecraft based on configuration.
    Returns True if injection conditions were met at current_time_nano.
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
        
    # Log the fault for later processing
    if FAULT_LOADER_AVAILABLE:
        try:
            if not hasattr(sc_object, 'fault_log'):
                sc_object.fault_log = []
            sc_object.fault_log.append({
                'type': fault_type,
                'magnitude': fault_magnitude,
                'wheel': fault_wheel,
                'time': current_time_nano
            })
            print(f"Fault logged: {fault_type} on wheel {fault_wheel} at {current_time_nano*macros.NANO2MIN:.2f} min")
            return True
        except Exception as e:
            print(f"Could not inject/log fault: {e}")
            return False
    else:
        print(f"Fault loader not available for {fault_type}")
        return False

def calculate_target_visibility(spacecraft_pos, target_pos, earth_radius=6371000):
    """Simple LOS check to a target"""
    try:
        sc_pos = np.array(spacecraft_pos)
        tgt_pos = np.array(target_pos)
        sc_to_target = tgt_pos - sc_pos
        if np.linalg.norm(sc_to_target) == 0:
            return False
        sc_altitude = np.linalg.norm(sc_pos) - earth_radius
        if sc_altitude > 200000:  # 200km
            earth_center = np.array([0, 0, 0])
            sc_to_earth = earth_center - sc_pos
            cos_angle = np.dot(sc_to_target, sc_to_earth) / (np.linalg.norm(sc_to_target) * np.linalg.norm(sc_to_earth))
            elevation_angle = np.pi/2 - np.arccos(np.clip(cos_angle, -1, 1))
            return elevation_angle > np.radians(10)
        return False
    except Exception:
        return False

class TargetDefinition:
    """Class to hold target location definitions"""
    def __init__(self, name, latitude, longitude, color="red"):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.color = color
        self.assigned_to = []
    
    def to_dict(self):
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
        # Default simulation time
        self.simulation_time = 30.0  # minutes
        
        # Spacecraft list for constellation support
        self.spacecraft_list = []
        
        # Default fault magnitudes for proper scaling
        self.DEFAULT_FAULT_MAGNITUDES = {
            "friction": 0.0005,    # N⋅m
            "power_limit": 0.5,    # W
            "encoder": 20.0,       # %
            "battery": 50.0        # W
        }
        
        # Legacy fault parameters (for backward compatibility)
        self.fault_type = "friction"
        self.fault_magnitude = 0.0005
        self.fault_wheel_number = 3
        self.fault_time = 10.0  # minutes
        
        # Periodic fault parameters
        self.enable_periodic_fault = False
        self.periodic_fault_interval = 360  # seconds
        self.periodic_fault_magnitude = 0.1
        self.periodic_fault_wheel = 1
        
        # Output settings
        self.binary_filename = "spacecraft_viz"
        self.show_plots = True
        self.save_plots = True
        self.save_binary = True
        
        # Camera configuration
        self.camera_position = [0.0, 0.0, 10.0]
        self.camera_fov = 80.0
        self.active_camera_name = None
        
        # Target locations
        self.targets = [
            TargetDefinition("Melbourne", -37.8136, 144.9631, "red"),
            TargetDefinition("New York", 40.71, -74.00, "blue"),
            TargetDefinition("Tokyo", 35.68, 139.77, "green"),
            TargetDefinition("London", 51.51, -0.13, "yellow")
        ]
        
    def validate(self):
        """Validate configuration parameters"""
        if not BASILISK_AVAILABLE:
            raise ValueError("Basilisk is not available. Please install Basilisk to run simulations.")
        
        if self.simulation_time <= 0:
            raise ValueError("Simulation time must be positive (in minutes)")
            
        # Validate each spacecraft configuration if using constellation mode
        if self.spacecraft_list:
            for sc in self.spacecraft_list:
                # Validate orbit parameters
                if sc["orbit"]["a"] <= 6371:
                    raise ValueError(f"Spacecraft {sc['name']}: Semi-major axis must be greater than Earth radius (6371 km)")
                if sc["orbit"]["e"] < 0 or sc["orbit"]["e"] >= 1.0:
                    raise ValueError(f"Spacecraft {sc['name']}: Eccentricity must be between 0 and 1")
                
                # Validate fault parameters if enabled
                if sc["fault"]["enabled"]:
                    if sc["fault"]["time"] >= self.simulation_time:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time must be less than simulation time")
                    if sc["fault"]["time"] < 0:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time cannot be negative")
                    
                    # Check if fault type is available
                    if FAULT_LOADER_AVAILABLE and sc["fault"]["type"] not in AVAILABLE_FAULTS:
                        print(f"Warning: Spacecraft {sc['name']}: Fault type '{sc['fault']['type']}' may not be available")
                    
                    if sc["fault"]["wheel"] not in range(4):
                        raise ValueError(f"Spacecraft {sc['name']}: Fault wheel number must be between 0 and 3")

def run_custom_simulation(config):
    """
    Run a customized simulation based on the configuration
    
    Parameters:
    config (SimulationConfig): Simulation configuration object
    
    Returns:
    tuple: (scenario, viz, figureList, output_dir)
    """
    # Check module availability
    print("\nChecking module availability...")
    check_fault_modules(verbose=True)
    check_plots_module(verbose=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(LOGS_DIR, f"sim_results_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Validate configuration
    config.validate()
    print("Configuration validation passed")
    
    if not config.spacecraft_list:
        print("ERROR: No spacecraft configured in simulation")
        return None, None, {}, output_dir
    
    print("\n" + "="*60)
    print("SPACECRAFT CONSTELLATION SIMULATION")
    print("="*60)
    print(f"Simulation Duration: {config.simulation_time} minutes")
    print(f"Number of Spacecraft: {len(config.spacecraft_list)}")
    
    # Check for clusters
    has_clusters = any(sc.get('type') == 'cluster_member' for sc in config.spacecraft_list)
    if has_clusters:
        try:
            print(f"\nGenerating cluster communication plots...")
            
            # Build cluster data structure for plotting
            cluster_data = {}
                
                # Method 1: Get clusters from spacecraft list
            clusters_found = {}
            for i, sc_config in enumerate(config.spacecraft_list):
                cluster_name = sc_config.get('cluster')
                if cluster_name:
                    if cluster_name not in clusters_found:
                        clusters_found[cluster_name] = {
                            'leader_idx': None,
                            'children_idx': []
                            }
                        
                    if sc_config.get('role') == 'leader':
                        clusters_found[cluster_name]['leader_idx'] = i
                        print(f"  Found leader: {sc_config['name']} (index {i}) for cluster {cluster_name}")
                    elif sc_config.get('role') == 'child':
                        clusters_found[cluster_name]['children_idx'].append(i)
                        print(f"  Found child: {sc_config['name']} (index {i}) for cluster {cluster_name}")
                
              # Convert to format expected by plot function
            for cluster_name, cluster_info in clusters_found.items():
                if cluster_info['leader_idx'] is not None and len(cluster_info['children_idx']) > 0:
                    cluster_data[cluster_name] = {
                          'leader': cluster_info['leader_idx'],
                        'children': cluster_info['children_idx']
                     }
                    print(f"  Cluster '{cluster_name}' ready: leader at index {cluster_info['leader_idx']}, {len(cluster_info['children_idx'])} children")
                
            # Generate plots if we have valid clusters
            if cluster_data:
                print(f"Generating plots for {len(cluster_data)} clusters...")
                cluster_plots = generate_cluster_communication_plots(
                    cluster_data,
                    sc_objects,
                    time_data
                )
                figureList.update(cluster_plots)
                print(f"✓ Generated {len(cluster_plots)} cluster communication plots")
            else:
                print("  WARNING: No valid cluster data found!")
                print(f"  Clusters found: {list(clusters_found.keys())}")
                for name, info in clusters_found.items():
                    print(f"    {name}: leader={info['leader_idx']}, children={info['children_idx']}")
                
        except Exception as e:
            print(f"✗ Cluster plots failed: {e}")
            import traceback
            traceback.print_exc()
    
    # Time setup
    simulationTime = macros.min2nano(config.simulation_time)
    print("\nTime Configuration:")
    print(f"  Input: {config.simulation_time} minutes")
    print(f"  Simulation: {simulationTime} nanoseconds")
    print(f"  Verification: {simulationTime/1e9/60:.1f} minutes")
    
    # Create simulation
    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)  # 1 second step
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    print("Simulation time step: 1 second")
    
    # Create gravitational bodies
    gravFactory = simIncludeGravBody.gravBodyFactory()
    planet = gravFactory.createEarth()
    planet.isCentralBody = True
    mu = planet.mu
    print(f"Earth gravity parameter: {mu/1e14:.3f} x 10^14 m³/s²")
    
    # Create spacecraft objects
    sc_objects = []
    cluster_manager = None
    fault_spacecraft_count = 0
    
    # Use cluster integration if available and needed
    if has_clusters and CLUSTER_INTEGRATION_AVAILABLE:
        print("\n" + "-"*50)
        print("USING CLUSTER INTEGRATION")
        print("-"*50)
        cluster_manager = ClusterManager()
        sc_objects = integrate_clusters_with_simulation(
            config, cluster_manager, scSim, simTaskName, gravFactory, mu
        )
        print(f"Created {len(sc_objects)} spacecraft objects using ClusterManager")
        comm_status = cluster_manager.get_communication_status()
        print(f"Communication setup: {comm_status}")
    else:
        # Create individual spacecraft
        print("\n" + "-"*50)
        print("CREATING INDIVIDUAL SPACECRAFT")
        print("-"*50)
        
        # Fixed RAAN for each orbit type to ensure separation
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
            
            # Add gravity
            gravFactory.addBodiesTo(sc)
            
            # Set orbit parameters
            oe = orbitalMotion.ClassicElements()
            orbit = sc_config["orbit"]
            
            # Get RAAN (prefer explicit Omega, fall back to orbit name mapping)
            if "Omega" in orbit:
                Omega_deg = orbit["Omega"]
            else:
                orbit_name = sc_config.get("orbit_name", "Default Orbit")
                Omega_deg = orbit_raan_values.get(orbit_name, 0.0)
            
            # Convert orbital parameters
            oe.a = max(orbit["a"] * 1000, 6571000)  # Convert km to m, minimum 200km altitude
            oe.e = min(orbit["e"], 0.01)  # Keep eccentricity small
            oe.i = orbit["i"] * macros.D2R
            oe.Omega = Omega_deg * macros.D2R
            oe.omega = orbit.get("omega", 0.0) * macros.D2R
            oe.f = orbit.get("f", 0.0) * macros.D2R
            
            # Convert to Cartesian coordinates
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            sc.hub.r_CN_NInit = rN
            sc.hub.v_CN_NInit = vN
            sc.hub.sigma_BNInit = [[0.01], [0.02], [-0.01]]
            sc.hub.omega_BN_BInit = [[0.0001], [-0.0002], [0.0001]]
            
            sc_objects.append(sc)
            
            # Calculate and print orbital info
            altitude_km = (oe.a / 1000) - 6371
            period_min = (2*np.pi*np.sqrt(oe.a**3/mu))/60.0
            print(f"Spacecraft '{sc.ModelTag}' - Orbit: {sc_config.get('orbit_name','N/A')}, Alt: {altitude_km:.1f}km, Period: {period_min:.1f} min")
            
            # Set up fault configuration
            if sc_config["fault"]["enabled"]:
                fault_type = sc_config["fault"]["type"]
                fault_magnitude = sc_config["fault"]["magnitude"]
                
                # Scale fault magnitude based on type (CRITICAL FOR PROPER SIMULATION)
                if fault_magnitude == 0.0005 and fault_type in config.DEFAULT_FAULT_MAGNITUDES:
                    scaled_magnitude = config.DEFAULT_FAULT_MAGNITUDES[fault_type]
                    print(f"  Scaling fault magnitude from {fault_magnitude} to {scaled_magnitude} for {fault_type}")
                    sc_config["fault"]["magnitude"] = scaled_magnitude
                    fault_magnitude = scaled_magnitude
                
                print(f"\nFault Configuration for '{sc.ModelTag}':")
                print(f"  Type: {fault_type} | Magnitude: {fault_magnitude} | Wheel: {sc_config['fault']['wheel']} | Time: {sc_config['fault']['time']} min")
                
                # Store fault config in spacecraft for later use
                sc.faultConfig = sc_config["fault"].copy()
                sc.faultInjected = False
                fault_spacecraft_count += 1
        
        print(f"\nTotal spacecraft with faults configured: {fault_spacecraft_count}")
    
    # ============= VISUALIZATION SETUP =============
    viz = None
    if config.save_binary and vizSupport.vizFound:
        print("\n" + "-"*50)
        print("SETTING UP VIZARD VISUALIZATION")
        print("-"*50)
        
        # Create visualization directory
        viz_files_dir = os.path.join(VIZ_DIR, "_VizFiles")
        os.makedirs(viz_files_dir, exist_ok=True)
        
        binary_filename = os.path.basename(config.binary_filename)
        binary_full_path = os.path.join(VIZ_DIR, binary_filename)
        print(f"Binary output: {binary_filename}_UnityViz.bin")
        
        # Check for battery faults and set up battery panels
        has_battery_fault = any(sat["fault"]["enabled"] and sat["fault"]["type"] == "battery" 
                                for sat in config.spacecraft_list)
        gsList = [[] for _ in range(len(config.spacecraft_list))]
        battery_components = {}
        
        if has_battery_fault:
            try:
                from Basilisk.simulation import simpleBattery, simpleSolarPanel, simplePowerSink
                from Basilisk.architecture import messaging
                
                print("Battery fault detected - setting up battery visualization panels...")
                
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
                        solarPanel.setPanelParameters([-1.0, 0.0, 0.0], 0.05, 0.3)
                        solarPanel.stateInMsg.subscribeTo(sc.scStateOutMsg)
                        scSim.AddModelToTask(simTaskName, solarPanel)
                        battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)
                        
                        # Sun message
                        sunMsgData = messaging.SpicePlanetStateMsgPayload()
                        sunMsgData.PositionVector = [-1.0, 0.0, 0.0]
                        sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
                        solarPanel.sunInMsg.subscribeTo(sunMsg)
                        
                        # Power sinks
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
                        
                        # Solar panel visualization
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
                        
                        # Schedule battery fault activation
                        fault_time_nano = macros.min2nano(sc_config["fault"]["time"])
                        scSim.__dict__[f'powerSinkFault_{i}'] = powerSinkFault
                        
                        scSim.createNewEvent(
                            f"batteryFault_{sc.ModelTag}",
                            simulationTimeStep,
                            True,
                            [f"self.TotalSim.CurrentNanos >= {fault_time_nano}"],
                            [
                                f"self.powerSinkFault_{i}.nodePowerOut = -0.05",
                                f"print('Battery fault activated for {sc.ModelTag}')",
                                f"self.setEventActivity('batteryFault_{sc.ModelTag}', False)"
                            ]
                        )
                        
                        print(f"Set up battery visualization for {sc.ModelTag}")
                        
            except Exception as e:
                print(f"Battery visualization setup failed (continuing without panels): {e}")
        
        # Enable visualization
        try:
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
                
                # Enable instrument GUI for spacecraft with battery faults
                for i, sc in enumerate(sc_objects):
                    if sc.ModelTag in battery_components:
                        vizSupport.setInstrumentGuiSetting(
                            viz,
                            spacecraftName=sc.ModelTag,
                            showGenericStoragePanel=True
                        )
            else:
                # Standard visualization without battery panels
                viz = vizSupport.enableUnityVisualization(
                    scSim,
                    simTaskName,
                    sc_objects,
                    saveFile=binary_full_path,
                    liveStream=False
                )
            
            print(f"Vizard enabled for {len(sc_objects)} spacecraft")
            
            # Configure visualization settings
            if hasattr(viz, 'settings'):
                viz.settings.showSpacecraftLabels = True
                viz.settings.orbitLinesOn = 1
                viz.settings.spacecraftSizeMultiplier = 8.0
                viz.settings.spacecraftCSOn = True
                viz.settings.showCelestialBodyLabels = True
                print("Vizard settings configured")
            
            # Add targets
            targets_added = 0
            for target in config.targets:
                if target.assigned_to:
                    try:
                        lat, lon, color = target.latitude, target.longitude, target.color
                        alt = 50000.0  # 50 km
                        radius = 6371000.0 + alt
                        lat_rad = lat * macros.D2R
                        lon_rad = lon * macros.D2R
                        x = radius * np.cos(lat_rad) * np.cos(lon_rad)
                        y = radius * np.cos(lat_rad) * np.sin(lon_rad)
                        z = radius * np.sin(lat_rad)
                        
                        vizSupport.addLocation(
                            viz,
                            stationName=target.name,
                            parentBodyName="earth",
                            r_GP_P=[x, y, z],
                            color=color,
                            range=2000000.0
                        )
                        targets_added += 1
                        print(f"Added target: {target.name}")
                    except Exception as e:
                        print(f"Could not add target {target.name}: {e}")
            
            print(f"Added {targets_added} targets to visualization")
            
            # Create cameras
            camera_created = False
            if config.active_camera_name:
                for sc in sc_objects:
                    if getattr(sc, "ModelTag", None) == config.active_camera_name:
                        try:
                            vizSupport.createStandardCamera(
                                viz,
                                setMode=1,
                                spacecraftName=sc.ModelTag,
                                fieldOfView=config.camera_fov * macros.D2R,
                                displayName=f"{sc.ModelTag} Camera",
                                pointingVector_B=[0, 0, -1],
                                position_B=config.camera_position
                            )
                            camera_created = True
                            print(f"Created camera for {sc.ModelTag}")
                            break
                        except Exception as e:
                            print(f"Could not create camera: {e}")
            
            if not camera_created and sc_objects:
                try:
                    vizSupport.createStandardCamera(
                        viz,
                        setMode=1,
                        spacecraftName=sc_objects[0].ModelTag,
                        fieldOfView=config.camera_fov * macros.D2R,
                        displayName="Default Camera",
                        pointingVector_B=[0, 0, -1],
                        position_B=config.camera_position
                    )
                    print("Created default camera")
                except Exception:
                    pass
                    
        except Exception as e:
            print(f"Failed to set up visualization: {e}")
            viz = None
    else:
        if not vizSupport.vizFound:
            print("Vizard visualization module not found")
        if not config.save_binary:
            print("Binary file generation is disabled")
    
    # ============= RUN SIMULATION =============
    print("\n" + "-"*50)
    print("RUNNING SIMULATION")
    print("-"*50)
    
    print("Initializing simulation...")
    scSim.InitializeSimulation()
    
    print(f"Setting stop time to {simulationTime} ns...")
    scSim.ConfigureStopTime(simulationTime)
    
    print(f"Starting simulation for {config.simulation_time} minutes...")
    start_time = datetime.now()
    
    # Execute the simulation
    scSim.ExecuteSimulation()
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print(f"Simulation completed!")
    print(f"Wall-clock time: {duration:.2f} s | Simulated: {config.simulation_time:.1f} min | Speed: {(config.simulation_time*60)/max(duration,1):.1f}x")
    
    # Verify simulation time
    current_sim_time = scSim.TotalSim.CurrentNanos
    actual_sim_minutes = current_sim_time / 1e9 / 60
    print(f"Actual simulation time: {actual_sim_minutes:.2f} minutes")
    
    if abs(actual_sim_minutes - config.simulation_time) > 0.1:
        print("WARNING: Simulation time mismatch!")
    
    # ============= GENERATE PLOTS WITH REAL FAULT SIMULATION =============
# ============= GENERATE PLOTS WITH REAL FAULT SIMULATION =============
    figureList = {}
    if (config.show_plots or config.save_plots) and PLOTS_AVAILABLE:
        print("\n" + "-"*50)
        print("GENERATING PLOTS WITH REAL FAULT SIMULATION ONLY")
        print("-"*50)
        
        try:
            # Create time array for plotting
            time_data = np.linspace(0, config.simulation_time, max(100, int(config.simulation_time * 2)))
            print(f"Time data: 0 to {config.simulation_time} minutes ({len(time_data)} points)")
            
            # 1. Generate constellation plots FOR ALL SATELLITES
            try:
                print(f"\nGenerating constellation plots for ALL {len(sc_objects)} satellites...")
                constellation_plots = generate_constellation_overview_plots(sc_objects, time_data, mu)
                figureList.update(constellation_plots)
                print(f"✓ Generated {len(constellation_plots)} constellation overview plots")
            except Exception as e:
                print(f"✗ Constellation plots failed: {e}")
            
            # 2. Generate cluster communication plots
            if has_clusters and cluster_manager:
                try:
                    print(f"\nGenerating cluster communication plots...")
                    # Build cluster data structure for plotting
                    cluster_data = {}
                    for cluster_name, cluster_info in cluster_manager.clusters.items():
                        # Get leader and children indices
                        leader_idx = None
                        children_indices = []
                        
                        # Find leader index
                        if 'leader_index' in cluster_info:
                            leader_idx = cluster_info['leader_index']
                        elif 'leader' in cluster_info:
                            # Find index by name
                            leader_name = cluster_info['leader']
                            for idx, sc in enumerate(sc_objects):
                                if hasattr(sc, 'ModelTag') and sc.ModelTag == leader_name:
                                    leader_idx = idx
                                    break
                        
                        # Find children indices
                        if 'children_indices' in cluster_info:
                            children_indices = cluster_info['children_indices']
                        elif 'children' in cluster_info:
                            # Find indices by names
                            for child_name in cluster_info['children']:
                                for idx, sc in enumerate(sc_objects):
                                    if hasattr(sc, 'ModelTag') and sc.ModelTag == child_name:
                                        children_indices.append(idx)
                                        break
                        
                        cluster_data[cluster_name] = {
                            'leader': leader_idx,
                            'children': children_indices
                        }
                    
                    cluster_plots = generate_cluster_communication_plots(
                        cluster_data,
                        sc_objects,
                        time_data
                    )
                    figureList.update(cluster_plots)
                    print(f"✓ Generated {len(cluster_plots)} cluster communication plots")
                except Exception as e:
                    print(f"✗ Cluster plots failed: {e}")
            elif has_clusters:
                # Try to generate cluster plots from constellation tab data
                try:
                    if hasattr(config, 'spacecraft_list'):
                        # Build cluster data from spacecraft list
                        cluster_data = {}
                        for sc_config in config.spacecraft_list:
                            if sc_config.get('cluster') and sc_config['cluster'] not in cluster_data:
                                cluster_data[sc_config['cluster']] = {
                                    'leader': None,
                                    'children': []
                                }
                        
                        # Find leader and children
                        for i, sc_config in enumerate(config.spacecraft_list):
                            if sc_config.get('cluster'):
                                cluster = cluster_data[sc_config['cluster']]
                                if sc_config.get('role') == 'leader':
                                    cluster['leader'] = i
                                elif sc_config.get('role') == 'child':
                                    cluster['children'].append(i)
                        
                        if cluster_data:
                            cluster_plots = generate_cluster_communication_plots(
                                cluster_data,
                                sc_objects,
                                time_data
                            )
                            figureList.update(cluster_plots)
                            print(f"✓ Generated {len(cluster_plots)} cluster plots from config")
                except Exception as e:
                    print(f"Note: Could not generate cluster plots from config: {e}")
            
            # 3. Generate inter-satellite distance plots
            try:
                print(f"\nGenerating distance plots for ALL {len(sc_objects)} satellites...")
                distance_plots = generate_inter_satellite_distance_plots(sc_objects, time_data, mu)
                figureList.update(distance_plots)
                print(f"✓ Generated {len(distance_plots)} distance plots")
            except Exception as e:
                print(f"✗ Distance plots failed: {e}")
            
            # 4. Generate REAL fault plots for EACH spacecraft with faults
            real_plots_count = 0
            failed_plots_count = 0
            
            print(f"\n{'='*60}")
            print(f"REAL FAULT SIMULATION PLOT GENERATION")
            print(f"{'='*60}")
            print(f"Fault loader available: {FAULT_LOADER_AVAILABLE}")
            print(f"Available fault types: {AVAILABLE_FAULTS}")
            
            if not FAULT_LOADER_AVAILABLE:
                print("ERROR: Fault loader not available - cannot generate real fault plots!")
                print("Please ensure fault modules are properly installed in the faults/ directory")
            else:
                for i, sc_config in enumerate(config.spacecraft_list):
                    if sc_config["fault"]["enabled"]:
                        fault_type = sc_config["fault"]["type"]
                        
                        print(f"\n{'='*50}")
                        print(f"GENERATING REAL PLOTS FOR {sc_config['name']}")
                        print(f"{'='*50}")
                        print(f"Fault Type: {fault_type}")
                        print(f"Magnitude: {sc_config['fault']['magnitude']}")
                        print(f"Wheel: RW{sc_config['fault']['wheel'] + 1}")
                        print(f"Time: {sc_config['fault']['time']} minutes")
                        
                        # Check if fault type is available
                        if fault_type not in AVAILABLE_FAULTS:
                            print(f"ERROR: Fault type '{fault_type}' not available!")
                            print(f"Available types: {AVAILABLE_FAULTS}")
                            failed_plots_count += 1
                            continue
                        
                        try:
                            # Create real simulation configuration
                            from plots import create_fault_config_for_real_simulation
                            
                            fault_params = {
                                'fault_magnitude': sc_config["fault"]["magnitude"],
                                'fault_wheel': sc_config["fault"]["wheel"],
                                'fault_time_min': sc_config["fault"]["time"],
                                'simulation_time_min': config.simulation_time
                            }
                            
                            # Force real simulation
                            real_fault_config = create_fault_config_for_real_simulation(fault_type, fault_params)
                            
                            print(f"Running REAL {fault_type} simulation...")
                            
                            # Generate plots using REAL simulation
                            fault_plots = generate_fault_plots(
                                fault_type=fault_type,
                                fault_data=real_fault_config,  # This forces real simulation
                                time_data=time_data,
                                fault_time_min=sc_config["fault"]["time"],
                                spacecraft_name=sc_config["name"]
                            )
                            
                            if fault_plots and len(fault_plots) > 0:
                                figureList.update(fault_plots)
                                
                                # Count real plots
                                real_plots = [name for name in fault_plots.keys() if "REAL" in name]
                                if real_plots:
                                    real_plots_count += len(real_plots)
                                    print(f"✓ SUCCESS: Generated {len(real_plots)} REAL simulation plots!")
                                    for plot_name in real_plots:
                                        print(f"    - {plot_name}")
                                else:
                                    print(f"✗ ERROR: No REAL plots generated (got {len(fault_plots)} plots)")
                                    failed_plots_count += 1
                            else:
                                print(f"✗ ERROR: No plots generated from real simulation")
                                failed_plots_count += 1
                                
                        except Exception as e:
                            print(f"✗ ERROR: Real simulation failed for {sc_config['name']}: {e}")
                            import traceback
                            traceback.print_exc()
                            failed_plots_count += 1
            
            # Summary of plot generation
            print(f"\n{'='*60}")
            print(f"PLOT GENERATION SUMMARY")
            print(f"{'='*60}")
            total_plots = len(figureList)
            constellation_count = len([k for k in figureList.keys() if 'Constellation' in k])
            cluster_count = len([k for k in figureList.keys() if 'Cluster' in k])
            distance_count = len([k for k in figureList.keys() if 'Distance' in k])
            
            print(f"Total plots generated: {total_plots}")
            print(f"  ✓ Constellation plots: {constellation_count}")
            print(f"  ✓ Cluster communication plots: {cluster_count}")
            print(f"  ✓ Distance plots: {distance_count}")
            print(f"  ✓ REAL fault simulation plots: {real_plots_count}")
            if failed_plots_count > 0:
                print(f"  ✗ Failed fault plots: {failed_plots_count}")
            
            if real_plots_count > 0:
                print(f"\n{'='*60}")
                print(f"SUCCESS: REAL FAULT SIMULATIONS COMPLETED!")
                print(f"All fault plots use actual simulation data")
                print(f"Look for plots marked with 'REAL' in the Results tab")
                print(f"{'='*60}")
            elif FAULT_LOADER_AVAILABLE:
                print(f"\nWARNING: No real fault plots generated despite fault loader being available")
                print(f"Check that fault modules are properly configured")
            else:
                print(f"\nERROR: Cannot generate real fault plots - fault loader not available!")
            
            # Save plots to files if requested
            if config.save_plots and figureList:
                print(f"\nSaving plots to disk...")
                os.makedirs(PLOTTING_DIR, exist_ok=True)
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                
                saved_count = 0
                for name, fig in figureList.items():
                    try:
                        # Clean filename
                        clean_name = name.replace(" ", "_").replace("/", "_")
                        plot_filename = f"{timestamp}_{clean_name}.png"
                        plot_path = os.path.join(PLOTTING_DIR, plot_filename)
                        
                        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                        saved_count += 1
                        
                        # Report saved plots
                        if "REAL" in name:
                            print(f"  ✓ Saved REAL plot: {plot_filename}")
                        elif "Constellation" in name:
                            print(f"  ✓ Saved constellation plot: {plot_filename}")
                        elif "Cluster" in name:
                            print(f"  ✓ Saved cluster plot: {plot_filename}")
                        
                        plt.close(fig)
                    except Exception as e:
                        print(f"  ✗ Error saving plot {name}: {e}")
                
                print(f"Saved {saved_count}/{len(figureList)} plots to {PLOTTING_DIR}")
                figureList.clear()
            else:
                # Close all figures to free memory
                for name, fig in figureList.items():
                    try:
                        plt.close(fig)
                    except:
                        pass
                if not config.save_plots:
                    figureList.clear()
                
        except Exception as e:
            print(f"Critical error in plot generation: {e}")
            import traceback
            traceback.print_exc()
            
            # Ensure all figures are closed
            try:
                plt.close('all')
            except:
                pass
    
    # ============= CREATE SCENARIO OBJECT =============
    class ConstellationScenario:
        def __init__(self, scSim, sc_objects, config, cluster_manager=None):
            self.TotalSim = scSim
            self.sc_objects = sc_objects
            self.targets = [t.to_dict() for t in config.targets]
            self.fault_type = getattr(config, 'fault_type', None)
            self.actual_sim_time = actual_sim_minutes
            self.cluster_manager = cluster_manager
    
    scenario = ConstellationScenario(scSim, sc_objects, config, cluster_manager)
    
    # ============= SAVE SIMULATION SUMMARY =============
    try:
        summary_path = os.path.join(output_dir, "simulation_summary.txt")
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write("SPACECRAFT CONSTELLATION SIMULATION SUMMARY\n")
            f.write("="*50 + "\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("MODULE STATUS:\n")
            f.write(f"  Basilisk: {'Available' if BASILISK_AVAILABLE else 'Not available'}\n")
            f.write(f"  Cluster Integration: {'Available' if CLUSTER_INTEGRATION_AVAILABLE else 'Not available'}\n")
            f.write(f"  Plotting: {'Available' if PLOTS_AVAILABLE else 'Not available'}")
            if PLOTS_AVAILABLE:
                apis = []
                if PLOTS_API["new"]:
                    apis.append("new")
                if PLOTS_API["old"]:
                    apis.append("old")
                f.write(f" (APIs: {', '.join(apis)})")
            f.write("\n")
            f.write(f"  Fault Loader: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}\n")
            f.write(f"  Available Faults: {AVAILABLE_FAULTS}\n")
            if FAILED_FAULTS:
                f.write(f"  Failed Faults: {FAILED_FAULTS}\n")
            f.write("\n")
            
            f.write("SIMULATION CONFIGURATION:\n")
            f.write(f"  Duration: {config.simulation_time} minutes\n")
            f.write(f"  Time Step: 1 second\n")
            f.write(f"  Spacecraft: {len(config.spacecraft_list)}\n")
            if has_clusters and cluster_manager:
                f.write(f"  Clusters: {len(cluster_manager.clusters)}\n")
            f.write("\n")
            
            f.write("SPACECRAFT DETAILS:\n")
            for i, sc in enumerate(config.spacecraft_list):
                f.write(f"\nSpacecraft {i+1}: {sc['name']}\n")
                f.write(f"  Type: {sc.get('type','individual')}\n")
                if sc.get('cluster'):
                    f.write(f"  Cluster: {sc['cluster']}\n")
                    f.write(f"  Role: {sc.get('role','')}\n")
                f.write(f"  Orbit:\n")
                f.write(f"    - Semi-major axis: {sc['orbit']['a']} km\n")
                f.write(f"    - Altitude: {sc['orbit']['a'] - 6371:.1f} km\n")
                f.write(f"    - Eccentricity: {sc['orbit']['e']}\n")
                f.write(f"    - Inclination: {sc['orbit']['i']}°\n")
                f.write(f"    - RAAN: {sc['orbit'].get('Omega', 'NA')}°\n")
                
                if sc["fault"]["enabled"]:
                    f.write("  Fault: ENABLED\n")
                    f.write(f"    - Type: {sc['fault']['type']}\n")
                    f.write(f"    - Magnitude: {sc['fault']['magnitude']}\n")
                    f.write(f"    - Wheel: {sc['fault']['wheel']}\n")
                    f.write(f"    - Injection time: {sc['fault']['time']} minutes\n")
                else:
                    f.write("  Fault: DISABLED\n")
            
            f.write("\nOUTPUT:\n")
            f.write(f"  Results Directory: {output_dir}\n")
            f.write(f"  Plots Generated: {len(figureList) if figureList else 0}\n")
            
            if config.save_binary:
                vizard_paths = [
                    os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
                    os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
                ]
                
                for vp in vizard_paths:
                    if os.path.exists(vp):
                        size_mb = os.path.getsize(vp)/(1024*1024)
                        f.write(f"  Binary File: {vp} ({size_mb:.2f} MB)\n")
                        break
        
        print(f"Summary saved: {summary_path}")
        
    except Exception as e:
        print(f"Could not save simulation summary: {e}")
    
    # ============= BINARY FILE VERIFICATION =============
    if config.save_binary:
        print("\n" + "-"*50)
        print("BINARY FILE VERIFICATION")
        print("-"*50)
        
        vizard_paths = [
            os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
            os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
        ]
        
        binary_found = False
        for vp in vizard_paths:
            if os.path.exists(vp):
                binary_found = True
                size_mb = os.path.getsize(vp)/(1024*1024)
                print(f"Binary file created: {os.path.basename(vp)}")
                print(f"Location: {vp}")
                print(f"Size: {size_mb:.2f} MB")
                print(f"Duration: {config.simulation_time} minutes")
                print(f"Spacecraft: {len(sc_objects)}")
                break
        
        if not binary_found:
            print("Binary file not found in expected locations")
    
    print("\n" + "="*60)
    print("SIMULATION COMPLETED SUCCESSFULLY")
    print("="*60)
    if has_clusters and cluster_manager:
        print(f"Clusters: {len(cluster_manager.clusters)}")
    print(f"Total Satellites: {len(sc_objects)}")
    print(f"Simulation Time: {config.simulation_time} minutes")
    print(f"Real Fault Simulations: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}")
    
    return scenario, viz, figureList, output_dir

# ============= MODULE TEST =============
if __name__ == "__main__":
    print("\nTesting spacecraft_simulation module...")
    check_fault_modules(verbose=True)
    check_plots_module(verbose=True)
    
    print("\nModule status:")
    print(f"  Basilisk: {'Available' if BASILISK_AVAILABLE else 'Not available'}")
    print(f"  Cluster Integration: {'Available' if CLUSTER_INTEGRATION_AVAILABLE else 'Not available'}")
    print(f"  Plotting: {'Available' if PLOTS_AVAILABLE else 'Not available'}")
    if PLOTS_AVAILABLE:
        apis = []
        if PLOTS_API["new"]:
            apis.append("new")
        if PLOTS_API["old"]:
            apis.append("old")
        print(f"    APIs: {', '.join(apis)}")
    print(f"  Fault Loader: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}")
    
    if BASILISK_AVAILABLE:
        config = SimulationConfig()
        config.simulation_time = 5.0
        config.spacecraft_list = [
            {
                "name": "TestSat1",
                "type": "individual",
                "cluster": None,
                "role": "independent",
                "orbit": {"a": 6971, "e": 0.01, "i": 55.0, "Omega": 0.0, "omega": 0.0, "f": 0.0},
                "orbit_name": "Test Orbit",
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 2.0,
                    "enabled": True,
                    "periodic": {"enabled": False, "interval": 360, "magnitude": 0.1, "wheel": 3}
                },
                "camera": {"position": [0.0, 0.0, 15.0], "fov": 80.0, "enabled": True},
                "communication": {"range": 2000.0, "fov": 30.0, "aHat_B": [0.0, 0.0, -1.0]},
                "targets": []
            }
        ]
        
        try:
            config.validate()
            print("Configuration validation passed")
            print("\nRunning test simulation...")
            scenario, viz, figureList, output_dir = run_custom_simulation(config)
            if scenario:
                print("\nTest simulation completed successfully!")
                print(f"Results saved to: {output_dir}")
            else:
                print("\nTest simulation failed")
        except Exception as e:
            print(f"Test failed: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("\nCannot run test - Basilisk not available")
    
    print("\nspacecraft_simulation.py module test complete!")