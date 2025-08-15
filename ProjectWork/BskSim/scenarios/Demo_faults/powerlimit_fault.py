import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from Basilisk.utilities import orbitalMotion, macros, vizSupport


# Get the file path information
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')

from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

###############################################################################
# SCENARIO DEFINITION
###############################################################################
class scenario_PowerLimitFault(BSKSim, BSKScenario):
    def __init__(self):
        super(scenario_PowerLimitFault, self).__init__()
        self.name = 'scenario_PowerLimitFault'

        # Declare message recorder names
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection parameters
        self.oneTimeRWFaultFlag = 1
        self.oneTimeFaultTime = macros.min2nano(15.)  # Inject fault at 30 minutes

        DynModels = self.get_DynModel()
        self.DynModels.RWFaultLog = []

        # Define Earth as the observation target
        self.targetPosition_N = np.array([0.0, 0.0, 0.0])  # Earth's center in inertial coordinates


        # Enable Vizard visualization
        if vizSupport.vizFound:
            viz = vizSupport.enableUnityVisualization(
                self,
                self.DynModels.taskName,
                self.DynModels.scObject,
                liveStream=True,
                saveFile="power_fault"
            )
            viz.settings.orbitLinesOn = 1
            viz.settings.showSpacecraftLabels = 1

            # Attach a forward-pointing camera
            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=self.DynModels.scObject.ModelTag,
                displayName="ScienceCam",
                fieldOfView=10 * macros.D2R,
                pointingVector_B=[0.0, 0.0, 0.0],  # Initial orientation
                position_B=[0.0, 1.5, 0.0]
            )
            
            self.viz = viz  # Store for later camera updates

    def configure_initial_conditions(self):
        # Set classical orbital elements
        oe = orbitalMotion.ClassicElements()
        oe.a = 1e7  # m
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModels = self.get_DynModel()
        mu = DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        DynModels.scObject.hub.r_CN_NInit = rN
        DynModels.scObject.hub.v_CN_NInit = vN

        # Set initial attitude and angular velocities
        DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # Reaction wheel speeds
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        # FSW outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])
        self.msgRecList[self.sNavTransName] = DynModel.simpleNavObject.transOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.sNavTransName])

        return

###############################################################################
# RUNNING THE SCENARIO WITH FAULT INJECTION
###############################################################################
def runScenario(scenario, powerLimit, faultIndex):
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    # Inject reaction wheel power limit fault
    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos >= self.oneTimeFaultTime and self.oneTimeRWFaultFlag == 1"],
        [f"self.DynModels.AddRWFault('powerLimit', {powerLimit}, {faultIndex}, self.TotalSim.CurrentNanos)",
         "self.oneTimeRWFaultFlag = 0"]
    )
    
    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    # Save visualization output file *after* simulation completes
    scenario.viz.settings.recordFile = "./_VizFiles/power_fault_UnityViz.bin"


###############################################################################
# PLOTTING CAMERA TRACKING PERFORMANCE
###############################################################################
def pull_outputs(scenario, showPlots, powerLimit, faultIndex):
    attErrRec = scenario.msgRecList[scenario.attGuidName]
    sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)  # Attitude error quaternion
    timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

    num_RW = 4
    RW_speeds = np.delete(scenario.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
    RW_speeds_norm = RW_speeds / np.max(RW_speeds, axis=0)  # Normalize speeds

    fault_time_min = scenario.oneTimeFaultTime * macros.NANO2MIN

    # Clear previous plots
    BSK_plt.clear_all_plots()

    # Plot Attitude Error
    plt.figure(figsize=(10, 5))
    plt.plot(timeData, np.linalg.norm(sigma_BR, axis=1), label="Attitude Error Norm")
    plt.axvline(fault_time_min, linestyle="--", color="red", label="Fault Injection")
    plt.xlabel("Time (min)")
    plt.ylabel("Attitude Error Norm")
    plt.title(f"Attitude Error (PL={powerLimit}W, FI={faultIndex})")
    plt.legend()
    plt.grid(True)

    # Plot RW Speed Changes
    plt.figure(figsize=(10, 5))
    for i in range(num_RW):
        plt.plot(timeData, RW_speeds_norm[:, i], label=f"RW {i+1} Speed (Norm)")
    plt.axvline(fault_time_min, linestyle="--", color="red", label="Fault Injection")
    plt.xlabel("Time (min)")
    plt.ylabel("Normalized RW Speed")
    plt.title(f"Reaction Wheel Speeds (PL={powerLimit}W, FI={faultIndex})")
    plt.legend()
    plt.grid(True)

    # Show plots
    if showPlots:
        plt.show()

###############################################################################
# MAIN EXECUTION
###############################################################################
if __name__ == "__main__":
    scenario = scenario_PowerLimitFault()
    runScenario(scenario, powerLimit=0.5, faultIndex=2)
    pull_outputs(scenario, showPlots=True, powerLimit=0.5, faultIndex=2)
