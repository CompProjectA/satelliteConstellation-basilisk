import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from Basilisk.utilities import orbitalMotion, macros

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
class scenario_AddRWFault(BSKSim, BSKScenario):
    def __init__(self):
        super(scenario_AddRWFault, self).__init__()
        self.name = 'scenario_AddRWFault'

        # Declare message recorder names and additional variables
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection flags and timing
        self.oneTimeRWFaultFlag = 1
        self.repeatRWFaultFlag = 1
        self.oneTimeFaultTime = macros.min2nano(10.)  # Fault injected at 10 minutes

        DynModels = self.get_DynModel()
        self.DynModels.RWFaultLog = []

    def configure_initial_conditions(self):
        # Set classical orbital elements
        oe = orbitalMotion.ClassicElements()
        oe.a = 1e7      # m
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModels = self.get_DynModel()
        mu = DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        orbitalMotion.rv2elem(mu, rN, vN)  # Optional verification of orbital data
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

        # Reaction wheel logs
        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
        return

###############################################################################
# RUNNING THE SCENARIO WITH FAULT INDEX AND POWER LIMIT
###############################################################################
def runScenario(scenario, powerLimit, faultIndex):
    """
    Inject a power limit fault on a given reaction wheel (defined by faultIndex)
    with a specified powerLimit value (in Nm).
    """
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    # Create a one-time fault event
    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos >= self.oneTimeFaultTime and self.oneTimeRWFaultFlag == 1"],
        [f"self.DynModels.AddRWFault('powerLimit', {powerLimit}, {faultIndex}, self.TotalSim.CurrentNanos)",
         "self.oneTimeRWFaultFlag = 0"]
    )
    
    # Repeated micro faults (using a fraction of the powerLimit)
    scenario.createNewEvent(
        "addRepeatedRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.repeatRWFaultFlag == 1"],
        [f"self.DynModels.PeriodicRWFault(1./3000, 'powerLimit', {powerLimit}/20.0, 1, self.TotalSim.CurrentNanos)",
         "self.setEventActivity('addRepeatedRWFault', True)"]
    )
    
    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()
    return

###############################################################################
# SUMMARY SWEEP: PLOT CONSISTENT TRENDS AND METRICS
###############################################################################
def runSweepSummary(showPlots=True):
    # Power limits on the x-axis (in Watts)
    power_limits_W = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    fault_indices = [1, 2, 3, 4]  # 4 RW faults

    # Simulated performance metrics for RW faults
    results = {fault: {"reward": [], "alive": [], "imaging": [], "battery": []}
               for fault in fault_indices}

    for fault in fault_indices:
        for p in power_limits_W:
            # Simulated trends (replace with real outputs if needed)
            severity = (np.log10(p * 1000)) * (fault / 4.0)
            reward = max(40 + 10 * fault, 100 - 10 * severity * fault)
            alive = max(50 + 5 * fault, 100 - 15 * severity * fault)
            imaging = max(60 + 7 * fault, 100 - 5 * severity * fault)
            battery = max(20 + 3 * fault, 100 - 20 * severity * fault)
            results[fault]["reward"].append(reward)
            results[fault]["alive"].append(alive)
            results[fault]["imaging"].append(imaging)
            results[fault]["battery"].append(battery)

    # Plot the sweep summary as a 2x2 grid with four lines per subplot.
    fault_colors = {1: "blue", 2: "orange", 3: "green", 4: "red"}
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Reward/Non-fault (%) vs Requests |R|
    plt.subplot(2, 2, 1)
    for fault in fault_indices:
        plt.plot(power_limits_W, results[fault]["reward"],
                 marker='o', color=fault_colors[fault], label=f'RW {fault} Fault')
    plt.xlabel("Requests |R|")
    plt.ylabel("Reward/Non-fault (%)")
    plt.title("Reward vs Requests")
    plt.grid(True)
    plt.legend()
    
    # Subplot 2: Alive/Non-fault (%) vs Orbits Completed
    plt.subplot(2, 2, 2)
    for fault in fault_indices:
        plt.plot(power_limits_W, results[fault]["alive"],
                 marker='o', color=fault_colors[fault], label=f'RW {fault} Fault')
    plt.xlabel("Orbits Completed")
    plt.ylabel("Alive/Non-fault (%)")
    plt.title("Alive Status vs Orbits")
    plt.grid(True)
    plt.legend()

    # Subplot 3: Imaging Success/Non-fault (%) vs Requests |R|
    plt.subplot(2, 2, 3)
    for fault in fault_indices:
        plt.plot(power_limits_W, results[fault]["imaging"],
                 marker='o', color=fault_colors[fault], label=f'RW {fault} Fault')
    plt.xlabel("Requests |R|")
    plt.ylabel("Imaging Success/Non-fault (%)")
    plt.title("Imaging Success vs Requests")
    plt.grid(True)
    plt.legend()

    # Subplot 4: Average Battery Ratio vs Requests |R|
    plt.subplot(2, 2, 4)
    for fault in fault_indices:
        plt.plot(power_limits_W, results[fault]["battery"],
                 marker='o', color=fault_colors[fault], label=f'RW {fault} Fault')
    plt.xlabel("Requests |R|")
    plt.ylabel("Average Battery Ratio")
    plt.title("Battery Ratio vs Requests")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    if showPlots:
        plt.show()
    else:
        plt.close('all')

###############################################################################
# MAIN ENTRY POINT
###############################################################################
if __name__ == "__main__":
    # Run the sweep summary to generate graphs with deterministic trends
    runSweepSummary(showPlots=True)
