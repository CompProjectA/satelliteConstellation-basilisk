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
class FrictionFaultScenario(BSKSim, BSKScenario):
    def __init__(self):
        super(FrictionFaultScenario, self).__init__()
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

        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

        # Record RW output torque messages
        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def pull_outputs(self, showPlots):

        # FSW process outputs, remove first data point as it is before FSW is called
        attErrRec = self.msgRecList[self.attGuidName]

        # B refers to the body frame (attached to the spacecraft).
        # R refers to the reference frame (desired orientation).
        # So, sigma_BR gives the orientation difference between the body and reference frames.

        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)

        # Extract recorded RW data (speed and friction)
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = []
        for i in range(numRW):
            RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))

        # Estimate RW temperatures based on speed and friction
        self.no_cooling, self.with_cooling = self.calculate_temperatures(RW_speeds, RW_friction)

        # Plotting section
        BSK_plt.clear_all_plots()
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

        # Plot RW speeds and friction
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.get_DynModel().RWFaultLog)
        BSK_plt.plot_attitude_error(timeData, sigma_BR)

        # Plot temperatures
        self.plot_rw_temperature(timeData, self.no_cooling, numRW)
        self.plot_rw_C_temperature(timeData, self.with_cooling, numRW)

        # Return or show/save figures
        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["RWSpeeds", "RWFriction","RWTemperatures(c)","RWTemperatures","attitudeErrorNorm" ]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        return figureList

    def calculate_temperatures(self, rw_speeds, rw_friction):
        """Estimate temperatures with and without cooling based on RW friction."""
        numRW = len(rw_friction)
        num_samples = len(rw_speeds)
        no_cooling = []
        with_cooling = []

        T_ambient = 10.0  # Ambient temperature in Celsius

        for rw in range(numRW):
            temp_nc = np.zeros(num_samples)
            temp_c = np.zeros(num_samples)
            temp_nc[0] = T_ambient
            temp_c[0] = T_ambient

            for i in range(1, num_samples):
                omega = rw_speeds[i, rw] * 2 * np.pi / 60
                P_friction = abs(rw_friction[rw][i] * omega) if i < len(rw_friction[rw]) else 0
                temp_rise = P_friction * 0.2
                cooling = 0.005 * (temp_c[i-1] - T_ambient)  # Arbitrary cooling term

                temp_nc[i] = temp_nc[i-1] + temp_rise
                temp_c[i] = temp_c[i-1] + (temp_rise) - cooling

            no_cooling.append(temp_nc)
            with_cooling.append(temp_c)

        return no_cooling, with_cooling




    def plot_rw_temperature(self, timeData, RW_temperatures, numRW):
        """Generate plot of RW temperatures over time."""
        plt.figure()
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('Reaction Wheel Temperatures(Without Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning')
        plt.tight_layout()

    def plot_rw_C_temperature(self, timeData, RW_temperatures, numRW):
        """Generate plot of RW temperatures over time."""
        plt.figure()
        colors = ['blue', 'green', 'red', 'cyan']

        for idx in range(numRW):
            plt.plot(timeData, RW_temperatures[idx], color=colors[idx],
                     label=f'RW {idx+1}', linewidth=2)

        plt.xlabel('Time [min]')
        plt.ylabel('Temperature [°C]')
        plt.title('Reaction Wheel Temperatures(With Cooling)')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Draw warning/critical temperature lines
        plt.axhline(y=30, color='orange', linestyle='--', alpha=0.7, label='Warning')
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
