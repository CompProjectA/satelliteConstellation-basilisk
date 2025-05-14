<<<<<<< HEAD
#!/usr/bin/env python
"""
friction_fault.py

A Basilisk scenario that simulates spacecraft dynamics with friction faults
in the reaction wheels and properly saves binary files for Vizard visualization.
"""
import inspect
import os
import sys
import numpy as np

from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Set paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import plotting.BSK_Plotting as BSK_plt
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class FrictionFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating friction faults in reaction wheels.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    """
    def __init__(self):
        super(FrictionFaultScenario, self).__init__()
        self.name = 'FrictionFaultScenario'
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

        self.oneTimeRWFaultFlag = 1
        self.repeatRWFaultFlag = 1
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.fault_magnitude = 0.0005
        self.fault_wheel_number = 3
        self.DynModels = self.get_DynModel()
        self.DynModels.RWFaultLog = []

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

    def log_outputs(self):
        """Configure message logging for analysis"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

    def pull_outputs(self, showPlots):
        """Process and plot simulation outputs"""
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = [
            np.delete(self.rwLogs[i].u_f, 0, 0) for i in range(numRW)
        ]

        # Get attitude logs
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)

        BSK_plt.clear_all_plots()
        timeData = np.delete(self.rwSpeedRec.times(), 0, 0) * macros.NANO2MIN
        
        # Generate plots
        BSK_plt.plot_attitude_error(timeData, sigma_BR, id=1)
        BSK_plt.plot_rate_error(timeData, omega_BR_B, id=2)
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW, id=3)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.DynModels.RWFaultLog, id=4)

        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["attitudeErrorNorm", "rateError", "RWSpeeds", "RWFriction"]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        return figureList
        
    def plot_target_visibility(self, timeData):
        """Plot target visibility based on spacecraft orbit"""
        return BSK_plt.plot_target_visibility(timeData, self.DynModels.scObject.hub.r_CN_N, self.targets)

def runScenario(scenario, saveBinary=True):
    """Run the friction fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        [f"self.get_DynModel().AddRWFault('friction',{scenario.fault_magnitude},{scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeRWFaultFlag=0"]
    )

    scenario.createNewEvent(
        "addRepeatedRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.repeatRWFaultFlag==1"],
        ["self.get_DynModel().PeriodicRWFault(360,'friction',0.1,1, self.TotalSim.CurrentNanos)", 
         "self.setEventActivity('addRepeatedRWFault',True)"]
    )

    viz = None
    if vizSupport.vizFound:
        viz = vizSupport.enableUnityVisualization(
            scenario,
            scenario.get_DynModel().taskName,
            scenario.get_DynModel().scObject,
            rwEffectorList=scenario.get_DynModel().rwStateEffector,
            liveStream=not saveBinary,
            saveFile="friction_fault_viz" if saveBinary else None
        )

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

        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="RW Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the friction fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Friction Fault Scenario =====")
    print(f"Save Binary: {saveBinary}")
    scenario = FrictionFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'friction_fault_viz.bin'")
        print("You can now open this file in Vizard for visualization.")

    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Friction Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)
=======
# Import necessary standard libraries
import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import utilities from Basilisk
from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Determine file path and append relevant directories to system path
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')

# Import simulation modules
from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

# Define the scenario class for RW fault analysis
class scenario_AddRWFault(BSKSim, BSKScenario):
    def __init__(self):
        super(scenario_AddRWFault, self).__init__()
        self.name = 'scenario_AddRWFault'  # Scenario name for identification
        self.msgRecList = {}  # Placeholder for any future message subscriptions
        self.sNavTransName = "sNavTransMsg"  # Navigation translation message name
        self.attGuidName = "attGuidMsg"  # Attitude guidance message name
        self.cameraLocation = [0.0, 3.0, 0.0]  # Camera position for visualization

        # Set the dynamics and flight software models
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        # Set up initial conditions and log desired outputs
        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection configuration
        self.oneTimeRWFaultFlag = 1  # Enable one-time fault
        self.repeatRWFaultFlag = 1  # Enable repeated fault
        self.oneTimeFaultTime = macros.min2nano(10.)  # Inject one-time fault at 10 minutes
        self.get_DynModel().RWFaultLog = []  # Initialize fault log

    def configure_initial_conditions(self):
        # Set classical orbital elements
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0  # Semi-major axis in meters
        oe.e = 0.01  # Eccentricity
        oe.i = 33.3 * macros.D2R  # Inclination (converted to radians)
        oe.Omega = 48.2 * macros.D2R  # RAAN
        oe.omega = 347.8 * macros.D2R  # Argument of perigee
        oe.f = 85.3 * macros.D2R  # True anomaly

        # Get dynamics model and compute position/velocity from orbital elements
        DynModel = self.get_DynModel()
        mu = DynModel.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)

        # Initialize spacecraft position and velocity
        DynModel.scObject.hub.r_CN_NInit = rN
        DynModel.scObject.hub.v_CN_NInit = vN

        # Initialize spacecraft attitude and angular velocity
        DynModel.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModel.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        # Set up message recording for RW speeds and friction torques
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # Record RW speed message
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        # Record RW output torque messages
        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def pull_outputs(self, showPlots):
        # Extract recorded RW data (speed and friction)
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = []
        for i in range(numRW):
            RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))

        # Estimate RW temperatures based on speed and friction
        RW_temperatures = self.calculate_temperatures(RW_speeds, RW_friction)

        # Plotting section
        BSK_plt.clear_all_plots()
        timeData = np.delete(self.rwSpeedRec.times(), 0, 0) * macros.NANO2MIN

        # Plot RW speeds and friction
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.get_DynModel().RWFaultLog)

        # Plot temperatures
        self.plot_rw_temperatures(timeData, RW_temperatures, numRW)

        # Return or show/save figures
        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["RWSpeeds", "RWFriction", "RWTemperatures"]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        return figureList

    def calculate_temperatures(self, rw_speeds, rw_friction):
        """Estimate temperatures based on power dissipated by RW friction."""
        numRW = len(rw_friction)
        num_samples = len(rw_speeds)
        temperatures = []

        T_ambient = 20.0  # Ambient temperature in Celsius

        for rw in range(numRW):
            temp = np.zeros(num_samples)
            temp[0] = T_ambient

            for i in range(1, num_samples):
                # Convert RPM to rad/s
                omega = rw_speeds[i, rw] * 2 * np.pi / 60 # power = torque × angular velocity

                # Compute power due to friction (P = torque × angular speed)
                P_friction = abs(rw_friction[rw][i] * omega) if i < len(rw_friction[rw]) else 0

                # Estimate temperature increase and cooling
                temp_rise = P_friction * 0.1  # Arbitrary scaling factor for heating
                cooling = (temp[i-1] - T_ambient) * 0.05  # Simplified cooling

                temp[i] = temp[i-1] + temp_rise - cooling

                # Clamp temperature between ambient and 100°C
                temp[i] = max(T_ambient, min(temp[i], 100.0))

            temperatures.append(temp)

        return temperatures

    def plot_rw_temperatures(self, timeData, RW_temperatures, numRW):
        """Generate plot of RW temperatures over time."""
        plt.figure()
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('Reaction Wheel Temperatures')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=22, color='orange', linestyle='--', alpha=0.7, label='Warning')
        plt.axhline(y=23, color='red', linestyle='--', alpha=0.7, label='Critical')
        plt.tight_layout()

    def runScenario(self):
        # Define simulation duration
        simulationTime = macros.min2nano(30.)
        self.modeRequest = "hillPoint"  # Request a specific control mode

        # If Unity viz is available, set up visualization
        if vizSupport.vizFound:
            viz = vizSupport.enableUnityVisualization(
                self,
                self.get_DynModel().taskName,
                self.get_DynModel().scObject,
                rwEffectorList=self.get_DynModel().rwStateEffector,
                liveStream=True,
                saveFile="friction_fault"
            )

            # Set up default camera
            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=self.get_DynModel().scObject.ModelTag,
                fieldOfView=30 * macros.D2R,
                displayName="RW Camera",
                pointingVector_B=[0, 0, 0],
                position_B=self.cameraLocation
            )

        # Run the simulation
        self.InitializeSimulation()
        self.ConfigureStopTime(simulationTime)
        self.ExecuteSimulation()

# Entry point for running the scenario as a script
def run(showPlots=True):
    scenario = scenario_AddRWFault()
    scenario.runScenario()  # Run the simulation
    scenario.pull_outputs(showPlots)  # Plot or save outputs

# Run the script if executed directly
if __name__ == "__main__":
    run(True)
>>>>>>> c93b1ec010aeecda176f654e7e5855c36bfc3ee0
