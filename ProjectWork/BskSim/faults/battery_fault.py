#!/usr/bin/env python
"""
Battery Fault Simulation Module

This module simulates battery degradation faults in spacecraft systems.
It includes power consumption modeling, fault injection, and safe-mode transitions.
"""

import os
import sys
import inspect
import matplotlib.pyplot as plt
import numpy as np
from Basilisk import __path__

# Path setup
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

sys.path.append(path + '/../')
sys.path.append(path + '/../modelsMultiSat')
sys.path.append(path + '/../plottingMultiSat')

bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])

# Import simulation modules
from Basilisk.simulation import spacecraft
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion,
                                simIncludeGravBody, unitTestSupport, vizSupport)
from Basilisk.simulation import simSynch
from Basilisk.simulation import simpleBattery 
from Basilisk.architecture import messaging
from Basilisk.simulation import simpleSolarPanel
from Basilisk.simulation import eclipse
from Basilisk.simulation import simplePowerSink

# Set additional paths
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, MODELS_DIR])

try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)


class BatteryFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating battery degradation faults.
    """
    def __init__(self):
        super(BatteryFaultScenario, self).__init__()
        self.name = 'BatteryFaultScenario'

    def run(self, **kwargs):
        """Run the battery fault scenario with given parameters."""
        return run(**kwargs)


def run(show_plots=False, liveStream=False, timeStep=1.0, orbitCase='LEO', 
        useSphericalHarmonics=False, planetCase='Earth', fault_magnitude=50.0, 
        fault_time_min=15.0, simulation_time_min=30.0, showPlots=None, saveBinary=False, **kwargs):
    """
    Main battery fault simulation function.
    
    Args:
        show_plots: Whether to display plots (legacy parameter)
        showPlots: Whether to display plots (GUI parameter) 
        liveStream: Whether to enable live streaming to Vizard
        timeStep: Simulation time step in seconds
        orbitCase: Orbit type ('LEO', 'GEO', 'GTO')
        useSphericalHarmonics: Whether to use spherical harmonics gravity
        planetCase: Planet ('Earth' or 'Mars')
        fault_magnitude: Fault power consumption in Watts
        fault_time_min: Time to inject fault in minutes
        simulation_time_min: Total simulation duration in minutes
        saveBinary: Whether to save binary files
    
    Returns:
        tuple: (scenario, viz, figureList) - format expected by fault_loader
    """
    
    # Handle parameter name variants and extract from kwargs
    if showPlots is not None:
        show_plots = showPlots
    
    # Extract parameters from kwargs if they're passed that way
    fault_magnitude = kwargs.get('fault_magnitude', fault_magnitude)
    fault_time_min = kwargs.get('fault_time_min', fault_time_min)
    simulation_time_min = kwargs.get('simulation_time_min', simulation_time_min)
    
    print(f"Running battery fault simulation with:")
    print(f"  - Fault magnitude: {fault_magnitude}W")
    print(f"  - Fault time: {fault_time_min} minutes")
    print(f"  - Simulation duration: {simulation_time_min} minutes")
    print(f"  - Orbit: {orbitCase}")
    print(f"  - Planet: {planetCase}")
    
    # Create simulation variable names
    simTaskName = "simTask"
    simProcessName = "simProcess"

    # Create a sim module as an empty container
    scSim = SimulationBaseClass.SimBaseClass()

    # Create the simulation process
    dynProcess = scSim.CreateNewProcess(simProcessName)

    # Create the dynamics task and specify the integration update time
    simulationTimeStep = macros.sec2nano(timeStep)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

    # Initialize spacecraft object and set properties
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = "bskSat"

    # Camera and target setup
    cameraLocation = [0.0, 2.0, 0.0]
    targets = [
        {"name": "Tokyo",          "lat":  35.68,   "lon": 139.77,   "color": "green"},
        {"name": "Sri Lanka",      "lat":   6.9271, "lon":  79.8612, "color": "orange"},
        {"name": "Central Africa", "lat":   1.5333, "lon":  17.6667, "color": "purple"}
    ]

    # Create and configure the battery
    battery = simpleBattery.SimpleBattery()
    battery.ModelTag = "satBattery"
    battery.storageCapacity = 100.0 
    battery.storedCharge_Init = 50.0
    scSim.AddModelToTask(simTaskName, battery)

    # Primary power sink (10 W baseline consumption)
    powerSink = simplePowerSink.SimplePowerSink()
    powerSink.ModelTag = "powerSink"
    powerSink.nodePowerOut = -0.01  # sink 10 W
    scSim.AddModelToTask(simTaskName, powerSink)
    battery.addPowerNodeToModel(powerSink.nodePowerOutMsg)

    # Secondary (fault-injection) load
    powerSinkFault = simplePowerSink.SimplePowerSink()
    powerSinkFault.ModelTag = "powerSinkFault"
    powerSinkFault.nodePowerOut = 0.0  # start off
    scSim.AddModelToTask(simTaskName, powerSinkFault)
    battery.addPowerNodeToModel(powerSinkFault.nodePowerOutMsg)

    # Expose components for event access
    scSim.powerSink = powerSink
    scSim.powerSinkFault = powerSinkFault

    # Record the battery status
    battRec = battery.batPowerOutMsg.recorder(simulationTimeStep)
    scSim.battRec = battRec
    scSim.AddModelToTask(simTaskName, battRec)

    # Calculate fault time in nanoseconds
    faultTime = macros.min2nano(fault_time_min)
    print(f"Fault will be injected at {fault_time_min} minutes ({fault_time_min * 60} seconds)")
    
    # Convert fault magnitude from W to kW (Basilisk units)
    fault_power_kw = fault_magnitude / 1000.0

    # Create fault injection event
    scSim.createNewEvent(
        "powerSinkFaultEvent",        
        simulationTimeStep,
        True,                   
        [f"self.TotalSim.CurrentNanos >= {faultTime}"],
        [
            f"print('*** BATTERY FAULT INJECTED: {fault_magnitude}W additional load ***')",
            f"self.powerSinkFault.nodePowerOut = -{fault_power_kw}",
            "self.setEventActivity('powerSinkFaultEvent', False)"
        ]
    )

    # Safe-mode event (20% battery threshold)
    safeThresh = battery.storageCapacity * 0.2
    cond = (
        f"(len(self.battRec.storageLevel) > 0) and "
        f"(self.battRec.storageLevel[-1] <= {safeThresh})"
    )
    scSim.createNewEvent(
        "batterySafeMode",
        simulationTimeStep,
        True,
        [cond],
        [
            "print('*** ENTERING SAFE MODE: Battery SOC ≤ 20% ***')",
            "self.powerSink.nodePowerOut = 0.0",
            "self.powerSinkFault.nodePowerOut = 0.0",
            "self.setEventActivity('batterySafeMode', False)"
        ]
    )

    # Solar Panel setup
    solarPanel = simpleSolarPanel.SimpleSolarPanel()
    solarPanel.ModelTag = "solarPanel"
    solarPanel.setPanelParameters([-1.0, -10.0, -1.0], 0.00001, 0.0000001)
    solarPanel.stateInMsg.subscribeTo(scObject.scStateOutMsg)
    scSim.AddModelToTask(simTaskName, solarPanel)
    battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)

    # Sun direction setup
    rawSun = np.array([-1.0, -10.0, -1.0])
    sunDir = (rawSun / np.linalg.norm(rawSun)).tolist()
    sunMsgData = messaging.SpicePlanetStateMsgPayload()
    sunMsgData.PositionVector = sunDir
    sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
    solarPanel.sunInMsg.subscribeTo(sunMsg)

    # Visualization Setup
    gsList = []
    
    if liveStream:
        # Battery visualization panel
        batteryPanel = vizSupport.vizInterface.GenericStorage()
        batteryPanel.label = "Battery (%)"
        batteryPanel.units = "%"
        batteryPanel.minValue = 0
        batteryPanel.maxValue = 100
        batteryPanel.useStorageLevel = True
        
        batteryInMsg = messaging.PowerStorageStatusMsgReader()
        batteryInMsg.subscribeTo(battery.batPowerOutMsg)
        batteryPanel.batteryStateInMsg = batteryInMsg
        batteryPanel.this.disown()

        # Color thresholds
        batteryPanel.thresholds = vizSupport.vizInterface.IntVector([20, 50, 80])
        batteryPanel.color = vizSupport.vizInterface.IntVector(
            vizSupport.toRGBA255("red") +
            vizSupport.toRGBA255("orange") +
            vizSupport.toRGBA255("yellow") +
            vizSupport.toRGBA255("green")
        )

        gsList.append([batteryPanel])

    # Add spacecraft object to the simulation process
    scSim.AddModelToTask(simTaskName, scObject)

    # Setup Gravity Body
    gravFactory = simIncludeGravBody.gravBodyFactory()
    if planetCase == 'Mars':
        planet = gravFactory.createMarsBarycenter()
        planet.isCentralBody = True
        if useSphericalHarmonics:
            planet.useSphericalHarmonicsGravityModel(bskPath + '/supportData/LocalGravData/GGM2BData.txt', 100)
    else:  # Earth
        planet = gravFactory.createEarth()
        planet.isCentralBody = True
        if useSphericalHarmonics:
            planet.useSphericalHarmonicsGravityModel(bskPath + '/supportData/LocalGravData/GGM03S-J2-only.txt', 2)
    
    mu = planet.mu
    gravFactory.addBodiesTo(scObject)

    # Setup orbit using classical orbit elements
    oe = orbitalMotion.ClassicElements()
    rLEO = 7000. * 1000      # meters
    rGEO = 42000. * 1000     # meters
    
    if orbitCase == 'GEO':
        oe.a = rGEO
        oe.e = 0.00001
        oe.i = 0.0 * macros.D2R
    elif orbitCase == 'GTO':
        oe.a = (rLEO + rGEO) / 2.0
        oe.e = 1.0 - rLEO / oe.a
        oe.i = 0.0 * macros.D2R
    else:  # LEO case, default
        oe.a = rLEO
        oe.e = 0.0001
        oe.i = 33.3 * macros.D2R
    
    oe.Omega = 48.2 * macros.D2R
    oe.omega = 347.8 * macros.D2R
    oe.f = 85.3 * macros.D2R
    rN, vN = orbitalMotion.elem2rv(mu, oe)
    oe = orbitalMotion.rv2elem(mu, rN, vN)

    # Initialize Spacecraft States
    scObject.hub.r_CN_NInit = rN  # m   - r_BN_N
    scObject.hub.v_CN_NInit = vN  # m/s - v_BN_N

    # Set the simulation time based on the parameter
    simulationTime = macros.sec2nano(simulation_time_min * 60.0)
    print(f"Simulation time set to: {simulation_time_min} minutes ({simulation_time_min * 60} seconds)")

    # Setup data logging
    numDataPoints = min(100, int(simulation_time_min * 2))  # Scale data points with sim time
    samplingTime = unitTestSupport.samplingTime(simulationTime, simulationTimeStep, numDataPoints)
    dataLog = scObject.scStateOutMsg.recorder(samplingTime)
    scSim.AddModelToTask(simTaskName, dataLog)

    # Record the battery data every time step
    batteryLog = battery.batPowerOutMsg.recorder(simulationTimeStep)
    scSim.AddModelToTask(simTaskName, batteryLog)

    # Setup live streaming if enabled
    if liveStream:
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)

        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, scObject,
                                                liveStream=True,
                                                genericStorageList=gsList)

        vizSupport.setInstrumentGuiSetting(viz, 
                                          spacecraftName=scObject.ModelTag,
                                          showGenericStoragePanel=True)

        # Add ground targets
        bodyName = planetCase.lower()
        earthRadius = 6371000.0

        for tgt in targets:
            lat_r = tgt["lat"] * macros.D2R
            lon_r = tgt["lon"] * macros.D2R
            x = earthRadius * np.cos(lat_r) * np.cos(lon_r)
            y = earthRadius * np.cos(lat_r) * np.sin(lon_r)
            z = earthRadius * np.sin(lat_r)
            location_position = [x, y, z]

            vizSupport.addLocation(
                viz,
                stationName=tgt["name"],
                parentBodyName=bodyName,
                r_GP_P=location_position,
                color=tgt["color"],
                range=2000e3 
            )

        # Add camera
        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scObject.ModelTag,
            displayName="BatteryCam",
            fieldOfView=30 * macros.D2R,
            pointingVector_B=[0.0, 0.0, 0.0],
            position_B=cameraLocation
        )

    # Initialize and run simulation
    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()

    # Retrieve logged data
    posData = dataLog.r_BN_N
    velData = dataLog.v_BN_N

    # Create plots
    plt.close("all")
    figureList = {}

    # Plot 1: Inertial position
    plt.figure(1)
    fig = plt.gcf()
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False, style='plain')
    time_minutes = dataLog.times() * macros.NANO2MIN
    for idx in range(3):
        plt.plot(time_minutes, posData[:, idx] / 1000.,
                color=unitTestSupport.getLineColor(idx, 3),
                label='$r_{BN,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [minutes]')
    plt.ylabel('Inertial Position [km]')
    plt.xlim(0, simulation_time_min)
    pltName = fileName + "1" + orbitCase + str(int(useSphericalHarmonics)) + planetCase
    figureList[pltName] = plt.figure(1)

    # Plot 2: Simple orbit plot
    plt.figure(2)
    plt.plot(posData[:, 0] / 1000, posData[:, 1] / 1000, 'b-', linewidth=2)
    plt.xlabel('X Position [km]')
    plt.ylabel('Y Position [km]')
    plt.title('Spacecraft Orbit')
    plt.grid(True)
    plt.axis('equal')
    pltName = fileName + "2" + orbitCase + str(int(useSphericalHarmonics)) + planetCase
    figureList[pltName] = plt.figure(2)

    # Plot 3: Battery storage level with fault injection
    plt.figure(3)
    timeData = batteryLog.times() * macros.NANO2MIN  # Convert to minutes
    storageData = batteryLog.storageLevel                  
    
    plt.plot(timeData, storageData, 'b-', linewidth=2, label='Battery Stored Charge')
    
    # Mark the fault injection moment
    plt.axvline(x=fault_time_min, color='r', linestyle='--', linewidth=2, 
               label=f'Fault Injected ({fault_magnitude}W additional load)')
    
    # Mark safe mode threshold
    safeThresh = battery.storageCapacity * 0.2
    plt.axhline(y=safeThresh, color='orange', linestyle=':', linewidth=2,
               label='Safe Mode Threshold (20%)')
    
    plt.xlabel('Time [minutes]')
    plt.ylabel('Stored Charge [Wh]')
    plt.title(f'Battery Storage Level with {fault_magnitude}W Fault Injection')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, battery.storageCapacity * 1.1)
    plt.xlim(0, simulation_time_min)
    
    # Print timing information
    print(f"Plot time range: {timeData[0]:.1f} to {timeData[-1]:.1f} minutes")
    print(f"Fault injection time: {fault_time_min:.1f} minutes")
    print(f"Simulation duration: {timeData[-1]:.1f} minutes")
    
    pltName = fileName + "3BatteryFault" + orbitCase + str(int(useSphericalHarmonics)) + planetCase
    figureList[pltName] = plt.figure(3)
    
    # Plot 4: Power Balance Analysis - Shows fault functionality
    plt.figure(4)
    
    # Calculate power consumption rate from battery data
    power_consumption = np.zeros_like(timeData)
    for i in range(1, len(timeData)):
        dt_hours = (timeData[i] - timeData[i-1]) / 60.0  # Convert minutes to hours
        if dt_hours > 0:
            dE = storageData[i-1] - storageData[i]  # Energy decrease (Wh)
            power_consumption[i] = dE / dt_hours  # Power (W)
    
    # Simple moving average for smoothing
    window_size = 5
    power_consumption_smooth = np.convolve(power_consumption, np.ones(window_size)/window_size, mode='same')
    
    # Create subplot layout
    plt.subplot(2, 1, 1)
    plt.plot(timeData, power_consumption_smooth, 'r-', linewidth=2, label='Total Power Consumption')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    plt.axhline(y=10, color='green', linestyle=':', alpha=0.7, label='Baseline (10W)')
    plt.axhline(y=10 + fault_magnitude, color='red', linestyle=':', alpha=0.7, 
               label=f'With Fault ({10 + fault_magnitude}W)')
    plt.xlabel('Time [minutes]')
    plt.ylabel('Power Consumption [W]')
    plt.title('Power Consumption Analysis - Battery Fault Verification')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, simulation_time_min)
    plt.ylim(0, max(70, 10 + fault_magnitude + 10))
    
    # Subplot 2: Discharge rate
    plt.subplot(2, 1, 2)
    discharge_rate = np.gradient(storageData, timeData)  # Wh/min
    discharge_rate_w = discharge_rate * 60  # Convert to W
    
    # Simple moving average
    window_size = 5
    discharge_rate_smooth = np.convolve(discharge_rate_w, np.ones(window_size)/window_size, mode='same')
    
    plt.plot(timeData, -discharge_rate_smooth, 'b-', linewidth=2, label='Discharge Rate')
    plt.axvline(x=fault_time_min, color='black', linestyle='--', linewidth=2, label='Fault Injection')
    
    # Find safe mode activation time
    safe_mode_idx = np.where(storageData <= safeThresh)[0]
    if len(safe_mode_idx) > 0:
        safe_mode_time = timeData[safe_mode_idx[0]]
        plt.axvline(x=safe_mode_time, color='orange', linestyle='--', linewidth=2, label='Safe Mode')
    
    plt.xlabel('Time [minutes]')
    plt.ylabel('Discharge Rate [W]')
    plt.title('Battery Discharge Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, simulation_time_min)
    plt.ylim(-10, max(70, 10 + fault_magnitude + 10))
    
    plt.tight_layout()
    pltName = fileName + "4BatteryPowerAnalysis" + orbitCase + str(int(useSphericalHarmonics)) + planetCase
    figureList[pltName] = plt.figure(4)

    if show_plots:
        plt.show()

    # Close plots to avoid memory issues
    plt.close("all")

    print(f"Battery fault returning {len(figureList)} figures: {list(figureList.keys())}")
    
    # Return the format expected by fault_loader: (scenario, viz, figure_list)
    # Create a simple scenario object with the required data
    class BatteryScenarioResult:
        def __init__(self, figure_list):
            self.figure_list = figure_list
            
        def pull_outputs(self, showPlots=False):
            return self.figure_list
    
    scenario = BatteryScenarioResult(figureList)
    return scenario, None, figureList



# --- Add to faults/battery_fault.py ---

def run_with_parameters(fault_magnitude=50.0,
                        fault_wheel=3,               # kept for API parity; not used by battery
                        fault_time_min=15.0,
                        simulation_time_min=30.0,
                        showPlots=False,
                        saveBinary=False,
                        **kwargs):
    """
    GUI-compatible wrapper that tolerates extra kwargs (e.g., battery_capacity_percentage)
    and calls run() with our standard parameters. Returns (scenario, viz, figureList).
    """
    try:
        scenario, viz, figureList = run(
            show_plots=showPlots,        # accept both spellings
            showPlots=showPlots,
            saveBinary=saveBinary,
            fault_magnitude=fault_magnitude,
            fault_time_min=fault_time_min,
            simulation_time_min=simulation_time_min
        )
        print(f"SUCCESS: Generated {len(figureList)} battery plots with GUI parameters")
        return scenario, viz, figureList
    except Exception as e:
        print(f"ERROR with battery parameters: {e}")
        import traceback; traceback.print_exc()
        return None, None, {}


def generate_battery_plots(fault_data, time_data, fault_time_min, spacecraft_name):
    """
    Generate battery fault plots using the standard plotting interface.
    
    This function is called by plots.py with standard parameters.
    """
    import matplotlib.pyplot as plt
    plots = {}
    
    print(f"Generating battery plots for {spacecraft_name}")
    print(f"Fault time: {fault_time_min} minutes")
    print(f"Time data range: {time_data[0]:.1f} to {time_data[-1]:.1f} minutes")
    
    # Extract fault parameters
    if isinstance(fault_data, dict):
        fault_magnitude = fault_data.get('battery_drain', fault_data.get('fault_magnitude', 50.0))
    else:
        fault_magnitude = 50.0  # Default
    
    # Create battery discharge simulation
    fig = plt.figure(figsize=(12, 6))
    
    # Simulate battery behavior
    initial_charge = 50.0  # Starting at 50% as in your simulation
    max_capacity = 100.0
    
    # Normal power consumption (10W baseline)
    baseline_consumption = 10.0  # Watts
    
    # Create time series in seconds
    time_seconds = time_data * 60.0  # Convert minutes to seconds
    battery_level = np.zeros_like(time_seconds)
    battery_level[0] = initial_charge
    
    # Simulate battery discharge
    for i in range(1, len(time_seconds)):
        dt = time_seconds[i] - time_seconds[i-1]  # Time step in seconds
        
        # Power consumption changes at fault time
        if time_data[i] >= fault_time_min:
            total_consumption = baseline_consumption + fault_magnitude  # Additional load from fault
        else:
            total_consumption = baseline_consumption
        
        # Simple discharge model (ignoring solar charging for simplicity)
        energy_consumed = total_consumption * dt / 3600.0  # Convert to Wh (W * hours)
        battery_level[i] = max(0, battery_level[i-1] - energy_consumed)
        
        # Stop at safe mode threshold (20%)
        if battery_level[i] <= max_capacity * 0.2:
            battery_level[i:] = battery_level[i]  # Maintain level (safe mode)
            break
    
    # Plot battery level
    plt.plot(time_data, battery_level, 'b-', linewidth=2, label='Battery Stored Charge')
    
    # Mark fault injection
    plt.axvline(x=fault_time_min, color='r', linestyle='--', linewidth=2, 
               label=f'Fault Injected ({fault_magnitude}W additional load)')
    
    # Mark safe mode threshold
    safe_threshold = max_capacity * 0.2
    plt.axhline(y=safe_threshold, color='orange', linestyle=':', linewidth=2,
               label='Safe Mode Threshold (20%)')
    
    plt.xlabel('Time [minutes]')
    plt.ylabel('Stored Charge [Wh]')
    plt.title(f'Battery Storage Level with {fault_magnitude}W Fault - {spacecraft_name}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, max_capacity * 1.1)
    
    plt.tight_layout()
    plots[f"BatteryFault_StorageLevel_{spacecraft_name}"] = fig
    
    return plots


# Main execution
if __name__ == "__main__":
    run(
        show_plots=True,
        liveStream=True,
        timeStep=1.0,
        orbitCase='LEO',
        useSphericalHarmonics=False,
        planetCase='Earth',
        fault_magnitude=50.0,
        fault_time_min=15.0,
        simulation_time_min=30.0
    )