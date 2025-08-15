#!/usr/bin/env python
"""
rw_encoder_fault.py

A Basilisk scenario that simulates spacecraft dynamics with encoder faults
in the reaction wheels and properly saves binary files for Vizard visualization.
"""
import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


from Basilisk.utilities import (orbitalMotion, macros, vizSupport)
from faults.rw_fault import RWFaultScenario

# Set paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')
sys.path.append(path + '/../BskSim')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import BSK_Plotting as BSK_plt
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

# ... [keep all imports and initial setup as-is]

class scenario_EncoderFault(RWFaultScenario):
    def __init__(self):
        super(scenario_EncoderFault, self).__init__()
        self.name = 'scenario_EncoderFault'

        # Declare message recorder names
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        # Fault injection parameters
        self.encoderFaultFlag = 1
        self.encoderFaultTime = macros.min2nano(10.)  # Inject fault at 10 minutes

        DynModels = self.get_DynModel()
        self.DynModels.EncoderFaultLog = []

        self.targetPosition_N = np.array([0.0, 0.0, 0.0])  # Earth's center

        if vizSupport.vizFound:
            viz = vizSupport.enableUnityVisualization(
                self,
                self.DynModels.taskName,
                self.DynModels.scObject,
                liveStream=True,
                saveFile="encoder_fault"
            )
            viz.settings.orbitLinesOn = 1
            viz.settings.showSpacecraftLabels = 1

            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=self.DynModels.scObject.ModelTag,
                displayName="ScienceCam",
                fieldOfView=100 * macros.D2R,
                pointingVector_B=[0.0, 0.0, 0.0],
                position_B=[0.0, 1.5, 0.0]
            )
            self.viz = viz

    def configure_initial_conditions(self):
        oe = orbitalMotion.ClassicElements()
        oe.a = 1e7
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

        DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])
        self.msgRecList[self.sNavTransName] = DynModel.simpleNavObject.transOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.sNavTransName])

        return


def runScenario(scenario, faultIndex):
    simulationTime = macros.min2nano(60.)
    scenario.modeRequest = "hillPoint"

    print(f"Fault Injection Set for {scenario.encoderFaultTime * macros.NANO2MIN} min, Current Sim Time: {scenario.TotalSim.CurrentNanos * macros.NANO2MIN} min")


    # Inject encoder fault on specified reaction wheel
    scenario.createNewEvent(
        "injectEncoderFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos >= self.encoderFaultTime and self.encoderFaultFlag == 1"], 
        [f"self.FSWModels.inject_rw_encoder_fault({faultIndex}, self.TotalSim.CurrentNanos)",  
        "self.encoderFaultFlag = 0"]
    )


    
    print(f"Injecting fault at time: {scenario.encoderFaultTime * macros.NANO2MIN} minutes")





    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    # Save viz file
    scenario.viz.settings.recordFile = "./_VizFiles/encoder_fault_UnityViz.bin"



def pull_outputs(scenario, showPlots, faultIndex):
    attErrRec = scenario.msgRecList[scenario.attGuidName]
    sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
    timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN

    num_RW = 4
    RW_speeds = np.delete(scenario.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
    RW_speeds_norm = RW_speeds / np.max(RW_speeds, axis=0)

    print("Initial RW Speeds:", RW_speeds[0])  # First set of RW speeds
    print("Final RW Speeds:", RW_speeds[-1])  # Last set of RW speeds
    print("Initial Attitude Error:", np.linalg.norm(sigma_BR[0]))  # First value
    print("Final Attitude Error:", np.linalg.norm(sigma_BR[-1]))  # Last value

    fault_time_min = scenario.encoderFaultTime * macros.NANO2MIN

    BSK_plt.clear_all_plots()

    plt.figure(figsize=(8, 2))
    plt.plot(timeData, np.linalg.norm(sigma_BR, axis=1), label="Attitude Error Norm")
    plt.axvline(fault_time_min, linestyle="--", color="red", label="Encoder Fault")
    plt.xlabel("Time (min)")
    plt.ylabel("Attitude Error Norm")
    plt.title(f"Attitude Error (Encoder Fault on RW {faultIndex})")
    plt.legend()
    plt.grid(True)

    # New Plot: Difference in RW Speeds Post-Fault
    plt.figure(figsize=(8, 2))
    diff_start_idx = np.argmax(timeData >= fault_time_min)
    rw_diff = RW_speeds[diff_start_idx:] - RW_speeds[diff_start_idx]

    for i in range(num_RW):
        plt.plot(timeData[diff_start_idx:], rw_diff[:, i], label=f"RW {i+1} Δ Speed")

    plt.xlabel("Time (min)")
    plt.ylabel("Δ RW Speed (rad/s)")
    plt.title(f"Change in RW Speeds After Encoder Fault on RW {faultIndex}")
    plt.legend()
    plt.grid(True)

    plt.figure(figsize=(8, 2))
    for i in range(num_RW):
        plt.plot(timeData, RW_speeds[:, i], label=f"RW {i+1} Speed (rad/s)")
    plt.axvline(fault_time_min, linestyle="--", color="red", label="Encoder Fault")
    plt.xlabel("Time (min)")
    plt.ylabel("RW Speed (rad/s)")
    plt.title(f"Reaction Wheel Speeds (Encoder Fault on RW {faultIndex})")
    plt.legend()
    plt.grid(True)

    if showPlots:
        plt.show()


if __name__ == "__main__":
    scenario = scenario_EncoderFault()
    runScenario(scenario, faultIndex=2)
    pull_outputs(scenario, showPlots=True, faultIndex=2)
