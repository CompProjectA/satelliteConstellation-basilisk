import os, sys, inspect
import numpy as np
from Basilisk.utilities import orbitalMotion, macros

# Setup matplotlib backend early
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt

# Setup paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.append(path + "/../")
sys.path.append(path + "/../models")
sys.path.append(path + "/../plotting")
sys.path.append(path + "/../BskSim")

from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw, BSK_Plotting as BSK_plt
from rwEncoderFault import rwEncoderFault

class scenario_EncoderFaultDemo(BSKSim, BSKScenario):
    def __init__(self):
        super(scenario_EncoderFaultDemo, self).__init__()
        print("[DEBUG] scenario_EncoderFaultDemo initialized")

        self.name = 'scenario_EncoderFaultDemo'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        print("[DEBUG] Setting dynamic and FSW models...")
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        print("[DEBUG] Configuring initial conditions...")
        self.configure_initial_conditions()

        print("[DEBUG] Logging outputs...")
        self.log_outputs()

        print("[DEBUG] Setting up fault injection...")
        #self.setup_fault() <-------------- currently issue that cause the plot to now show up

    def configure_initial_conditions(self):
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        mu = self.DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        self.DynModels.scObject.hub.r_CN_NInit = rN
        self.DynModels.scObject.hub.v_CN_NInit = vN
        self.DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        self.DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def setup_fault(self):
        print("[DEBUG] Creating rwFault module...")
        self.rwFault = rwEncoderFault()

        print("[DEBUG] Enabling fault...")
        self.rwFault.enableFault = True
        self.rwFault.faultyRW = 1
        self.rwFault.stuckSpeed = 100.0  # rad/s

        print("[DEBUG] Subscribing to rwSpeedOutMsg...")
        self.rwFault.rwInMsg.subscribeTo(self.DynModels.rwStateEffector.rwSpeedOutMsg)

        print("[DEBUG] Adding rwFault to task...")
        self.AddModelToTask(self.DynModels.taskName, self.rwFault)

        print("[DEBUG] Creating telemetry recorder from rwOutMsg...")
        self.telemetryRec = self.rwFault.rwOutMsg.recorder(self.FSWModels.processTasksTimeStep)

        print("[DEBUG] Adding telemetry recorder to task...")
        self.AddModelToTask(self.DynModels.taskName, self.telemetryRec)

        print("[DEBUG] Fault setup completed.")

    def log_outputs(self):
        print("[DEBUG] Logging FSW outputs...")
        self.attGuidRec = self.FSWModels.attGuidMsg.recorder(self.FSWModels.processTasksTimeStep)
        self.AddModelToTask(self.DynModels.taskName, self.attGuidRec)

        print("[DEBUG] Logging RW speeds...")
        self.rwSpeedRec = self.DynModels.rwStateEffector.rwSpeedOutMsg.recorder(self.FSWModels.processTasksTimeStep)
        self.AddModelToTask(self.DynModels.taskName, self.rwSpeedRec)


    def pull_outputs(self, showPlots):
        print("[INFO] Plotting results...")
        timeLine = self.attGuidRec.times() * macros.NANO2MIN
        sigma_BR = self.attGuidRec.sigma_BR
        omega_BR_B = self.attGuidRec.omega_BR_B
        rwSpeeds = np.atleast_2d(self.rwSpeedRec.wheelSpeeds)


        BSK_plt.clear_all_plots()
        BSK_plt.plot_attitude_error(timeLine, sigma_BR)
        BSK_plt.plot_rate_error(timeLine, omega_BR_B)
        BSK_plt.plot_rw_speeds(timeLine, rwSpeeds, 4)

        if showPlots:
            print("[INFO] Displaying plots...")
            plt.pause(0.001)
            plt.show(block=True)
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["attitudeErrorNorm", "rateError", "RWSpeeds", "RWFriction"]
            BSK_plt.save_all_plots(fileName, figureNames)
            print("[INFO] Saved plots to PNG files.")

def run(showPlots=True):
    print("[INFO] Running scenario_EncoderFaultDemo...")
    scenario = scenario_EncoderFaultDemo()
    scenario.modeRequest = "standby"

    print("[DEBUG] Initializing simulation...")
    scenario.InitializeSimulation()

    print("[DEBUG] Configuring stop time...")
    simulationTime = macros.min2nano(1.)  # Try just 1 minute
    scenario.ConfigureStopTime(simulationTime)

    print("[DEBUG] Executing simulation...")
    scenario.ExecuteSimulation()
    print(f"[DEBUG] Simulation time set to 1 minute")

    print("[DEBUG] Pulling outputs...")
    scenario.pull_outputs(showPlots)

    print("[INFO] Scenario run complete.")

if __name__ == "__main__":
    print("[INFO] Starting encoder fault scenario...")
    run(showPlots=True)
