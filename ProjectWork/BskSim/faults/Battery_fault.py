#!/usr/bin/env python
"""
battery_fault.py

A Basilisk scenario that simulates spacecraft dynamics with battery faults
and properly saves binary files for Vizard visualization.

For use with the spacecraft fault simulator GUI.
"""
import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Set paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import plotting.BSK_Plotting as BSK_plt
    from Basilisk.simulation import simpleBattery, simplePowerSink, simpleSolarPanel
    from Basilisk.architecture import messaging
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class BatteryFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating battery faults in spacecraft.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    """
    def __init__(self):
        super(BatteryFaultScenario, self).__init__()
        self.name = 'BatteryFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.cameraLocation = [0.0, 2.0, 0.0]

        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "red"},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "blue"},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "green"},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "yellow"}
        ]

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        self.oneTimeFaultFlag = 1  # Important - this is the flag for battery faults
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.fault_magnitude = 0.05  # Default power drain in kW (50W)
        self.fault_wheel_number = 0  # Not used for battery faults but needed for GUI compatibility
        self.DynModels = self.get_DynModel()
        self.DynModels.BatteryFaultLog = []
        
        # Battery and power components
        self.battery = None
        self.powerSink = None
        self.solarPanel = None
        self.batteryReader = None
        self.solarReader = None

    def configure_initial_conditions(self):
        """Configure orbit and attitude initial conditions"""
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModel = self.get_DynModel()
        mu = DynModel.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        orbitalMotion.rv2elem(mu, rN, vN)
        DynModel.scObject.hub.r_CN_NInit = rN
        DynModel.scObject.hub.v_CN_NInit = vN
        DynModel.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModel.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def setup_power_system(self, simTaskName):
        """Set up the battery and power components"""
        # Create and configure the battery
        self.battery = simpleBattery.SimpleBattery()
        self.battery.ModelTag = "satBattery"
        self.battery.storageCapacity = 100.0  # Amp-hours
        self.battery.storedCharge_Init = 100.0  # Start with full charge
        self.AddModelToTask(simTaskName, self.battery)
        
        # Create battery status reader for visualization
        self.batteryReader = messaging.PowerStorageStatusMsgReader()
        self.batteryReader.subscribeTo(self.battery.batPowerOutMsg)
        
        # Create main power sink (baseline power consumption)
        self.powerSink = simplePowerSink.SimplePowerSink()
        self.powerSink.ModelTag = "powerSink"
        self.powerSink.nodePowerOut = -0.01  # 10W baseline consumption
        self.AddModelToTask(simTaskName, self.powerSink)
        
        # Connect power sink to battery
        self.battery.addPowerNodeToModel(self.powerSink.nodePowerOutMsg)
        
        # Create solar panel
        self.solarPanel = simpleSolarPanel.SimpleSolarPanel()
        self.solarPanel.ModelTag = "solarPanel"
        self.solarPanel.setPanelParameters([-1.0, -10.0, -1.0], 0.0001, 0.000001)
        self.solarPanel.stateInMsg.subscribeTo(self.DynModels.scObject.scStateOutMsg)
        self.AddModelToTask(simTaskName, self.solarPanel)
        self.battery.addPowerNodeToModel(self.solarPanel.nodePowerOutMsg)
        
        # Configure sun direction for solar panel
        sunDir = np.array([-1.0, -10.0, -1.0])
        sunDir = (sunDir / np.linalg.norm(sunDir)).tolist()
        sunMsgData = messaging.SpicePlanetStateMsgPayload()
        sunMsgData.PositionVector = sunDir
        sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
        self.solarPanel.sunInMsg.subscribeTo(sunMsg)
        
        # Create solar panel power output reader
        self.solarReader = messaging.PowerNodeUsageMsgReader()
        self.solarReader.subscribeTo(self.solarPanel.nodePowerOutMsg)

    def log_outputs(self):
        """Configure message logging for analysis"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # RW speeds from dynamics (actual) - needed for plotting
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

    def log_battery_data(self, samplingTime):
        """Configure battery data logging"""
        if self.battery:
            self.batteryLog = self.battery.batPowerOutMsg.recorder(samplingTime)
            self.AddModelToTask(self.get_DynModel().taskName, self.batteryLog)
        
        if self.solarPanel:
            self.solarLog = self.solarPanel.nodePowerOutMsg.recorder(samplingTime)
            self.AddModelToTask(self.get_DynModel().taskName, self.solarLog)
        
        if self.powerSink:
            self.powerSinkLog = self.powerSink.nodePowerOutMsg.recorder(samplingTime)
            self.AddModelToTask(self.get_DynModel().taskName, self.powerSinkLog)

    def apply_battery_fault(self, fault_magnitude, time_nano):
        """
        Apply a battery fault by increasing power drain
        
        Parameters:
        fault_magnitude (float): Power drain increase in kW
        time_nano (float): Time to apply fault in nanoseconds
        """
        if self.powerSink:
            # Convert kW to normalized power units and make negative (power draw)
            power_draw = -fault_magnitude
            
            # Log the fault event
            self.DynModels.BatteryFaultLog.append({
                'type': 'battery',
                'time': time_nano,
                'magnitude': fault_magnitude
            })
            
            # Set the new power draw
            self.powerSink.nodePowerOut = power_draw
            
            print(f"Battery fault injected: {fault_magnitude} kW power drain at time {time_nano * macros.NANO2MIN:.2f} minutes")
            return True
        else:
            print("ERROR: Power sink not initialized, cannot apply battery fault")
            return False

    def pull_outputs(self, showPlots):
        """Process and plot simulation outputs"""
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)

        # Get attitude logs
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)

        # Clear all existing plots
        plt.close('all')
        
        # Get time data for plots
        timeData = np.delete(self.rwSpeedRec.times(), 0, 0) * macros.NANO2MIN
        fault_time_min = self.oneTimeFaultTime * macros.NANO2MIN
        
        # Generate plots
        # Attitude error plot (figure 1)
        plt.figure(1)
        plt.plot(timeData, np.linalg.norm(sigma_BR, axis=1))
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Attitude Error Norm')
        plt.title('Attitude Error')
        plt.grid(True)
        plt.legend()
        
        # Rate error plot (figure 2)
        plt.figure(2)
        for i in range(3):
            plt.plot(timeData, omega_BR_B[:, i], label=f'Axis {i+1}')
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Angular Rate Error [rad/s]')
        plt.title('Rate Error')
        plt.grid(True)
        plt.legend()
        
        # RW speeds plot (figure 3)
        plt.figure(3)
        for i in range(numRW):
            plt.plot(timeData, RW_speeds[:, i], label=f'RW{i+1}')
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Wheel Speed [rad/s]')
        plt.title('Reaction Wheel Speeds')
        plt.grid(True)
        plt.legend()
        
        # Battery specific plot (figure 4)
        if hasattr(self, 'batteryLog'):
            batteryData = np.delete(self.batteryLog.storageLevel, 0, 0)
            powerData = np.delete(self.batteryLog.powerOutFlow, 0, 0)
            
            plt.figure(4)
            plt.subplot(211)
            plt.plot(timeData, batteryData, 'b-', label='Battery Charge Level')
            plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
            plt.xlabel('Time [min]')
            plt.ylabel('Battery Level [%]')
            plt.title('Battery Charge Status')
            plt.grid(True)
            plt.legend()
            
            plt.subplot(212)
            plt.plot(timeData, powerData, 'g-', label='Power Flow')
            plt.axvline(x=fault_time_min, color='r', linestyle='--')
            plt.axhline(y=-self.fault_magnitude, color='m', linestyle='-.', label='Fault Power Level')
            plt.xlabel('Time [min]')
            plt.ylabel('Power Flow [kW]')
            plt.title('Battery Power Flow')
            plt.grid(True)
            plt.legend()
            
            plt.tight_layout()

        # Create figure list to return
        figureList = {
            "attitudeErrorNorm": plt.figure(1),
            "rateError": plt.figure(2),
            "RWSpeeds": plt.figure(3)
        }
        
        if hasattr(self, 'batteryLog'):
            figureList["BatteryStatus"] = plt.figure(4)

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList

def runScenario(scenario, saveBinary=True):
    """Run the battery fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"
    
    # Set up power system
    scenario.setup_power_system(scenario.get_DynModel().taskName)
    
    # Log battery data
    scenario.log_battery_data(scenario.get_FswModel().processTasksTimeStep)

    # Set up battery fault event
    scenario.createNewEvent(
        "injectBatteryFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeFaultFlag==1"],
        [f"self.apply_battery_fault({scenario.fault_magnitude}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeFaultFlag=0"]
    )

    viz = None
    if vizSupport.vizFound:
        # Create generic storage for battery visualization
        batteryPanel = vizSupport.vizInterface.GenericStorage()
        batteryPanel.label = "Battery (%)"
        batteryPanel.units = "%"
        batteryPanel.minValue = 0
        batteryPanel.maxValue = 100
        batteryPanel.useStorageLevel = True
        batteryPanel.batteryStateInMsg = scenario.batteryReader
        batteryPanel.this.disown()
        batteryPanel.thresholds = vizSupport.vizInterface.IntVector([20, 50, 80])
        batteryPanel.color = vizSupport.vizInterface.IntVector(
            vizSupport.toRGBA255("red") +
            vizSupport.toRGBA255("orange") +
            vizSupport.toRGBA255("yellow") +
            vizSupport.toRGBA255("green")
        )
        
        # Create solar panel visualization
        solarViz = vizSupport.vizInterface.GenericStorage()
        solarViz.label = "Solar Power"
        solarViz.units = "W"
        solarViz.minValue = 0.0
        solarViz.maxValue = 20.0
        solarViz.useStorageLevel = False
        solarViz.storageUnitStateInMsg = scenario.solarReader
        solarViz.this.disown()
        
        # List of generic storage elements
        gsList = [[batteryPanel, solarViz]]
        
        # Create visualization binary directory if it doesn't exist
        vizfiles_dir = os.path.join(VIZ_DIR, "_VizFiles")
        if not os.path.exists(vizfiles_dir):
            try:
                os.makedirs(vizfiles_dir, exist_ok=True)
            except:
                print(f"Warning: Could not create directory {vizfiles_dir}")
        
        # Enable visualization
        binary_filename = "battery_fault_viz" if saveBinary else None
        if saveBinary:
            binary_path = os.path.join(vizfiles_dir, binary_filename)
        else:
            binary_path = None
            
        viz = vizSupport.enableUnityVisualization(
            scenario,
            scenario.get_DynModel().taskName,
            scenario.get_DynModel().scObject,
            rwEffectorList=scenario.get_DynModel().rwStateEffector,
            liveStream=not saveBinary,
            saveFile=binary_path,
            genericStorageList=gsList
        )
        
        # Set up instrument GUI to show battery panel
        vizSupport.setInstrumentGuiSetting(viz, 
                                        spacecraftName=scenario.get_DynModel().scObject.ModelTag,
                                        showGenericStoragePanel=True)

        # Add targets
        for target in scenario.targets:
            lat = target["lat"]
            lon = target["lon"]
            color = target.get("color", "red")
            alt = 0.0
            radius = 6371000.0 + alt
            lat_rad = lat * macros.D2R
            lon_rad = lon * macros.D2R
            x = radius * np.cos(lat_rad) * np.cos(lon_rad)
            y = radius * np.cos(lat_rad) * np.sin(lon_rad)
            z = radius * np.sin(lat_rad)
            location_position = [x, y, z]

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
            setMode=1,  # Standard camera mode (attached to body)
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="RW Camera",
            pointingVector_B=[0, 0, 0],  # Look at spacecraft center
            position_B=scenario.cameraLocation  # Camera position in body frame
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the battery fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Battery Fault Scenario =====")
    print(f"Save Binary: {saveBinary}")
    scenario = BatteryFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'battery_fault_viz.bin'")
        print("You can now open this file in Vizard for visualization.")

    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Battery Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)