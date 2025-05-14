#!/usr/bin/env python
"""
rw_fault.py

A simplified Basilisk scenario that simulates spacecraft dynamics with reaction wheel faults
and properly saves binary files for Vizard visualization.
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

# Import BSK modules - modify these paths as needed for your environment
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import BSK_Plotting as BSK_plt
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    print("Make sure the BSK modules are in your Python path.")
    sys.exit(1)

class RWFaultScenario(BSKSim, BSKScenario):
    def __init__(self):
        super(RWFaultScenario, self).__init__()
        self.name = 'RWFaultScenario'
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
        self.get_DynModel().RWFaultLog = []

    def configure_initial_conditions(self):
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
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

    def pull_outputs(self, showPlots):
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        RW_friction = [
            np.delete(self.rwLogs[i].u_f, 0, 0) for i in range(numRW)
        ]

        BSK_plt.clear_all_plots()
        timeData = np.delete(self.rwSpeedRec.times(), 0, 0) * macros.NANO2MIN
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, numRW)
        BSK_plt.plot_rw_friction(timeData, RW_friction, numRW, self.get_DynModel().RWFaultLog)

        figureList = {}
        if showPlots:
            BSK_plt.show_all_plots()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["RWSpeeds", "RWFriction"]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)

        return figureList

def runScenario(scenario, saveBinary=True):
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    scenario.createNewEvent(
        "addOneTimeRWFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        ["self.get_DynModel().AddRWFault('friction',0.0005,3, self.TotalSim.CurrentNanos)", 
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
            saveFile="rw_fault_viz" if saveBinary else None
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
    print("\n===== Running Improved RW Fault Scenario =====")
    print(f"Save Binary: {saveBinary}")
    scenario = RWFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'rw_fault_viz.bin'")
        print("You can now open this file in Vizard for visualization.")

    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Improved RW Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)
