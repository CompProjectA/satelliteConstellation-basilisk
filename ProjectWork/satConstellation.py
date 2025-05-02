import os
import numpy as np
from Basilisk.simulation import spacecraft, spacecraftLocation, simSynch
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion, simIncludeGravBody, vizSupport)
from Basilisk import __path__
bskPath = __path__[0]

def run(show_plots=True):
    print("Basilisk Path:", bskPath)
    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    satellites = []
    num_satellites = 4
    for i in range(num_satellites):
        sat = spacecraft.Spacecraft()
        sat.ModelTag = f"Satellite{i+1}"
        satellites.append(sat)
        scSim.AddModelToTask(simTaskName, sat)
    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    earth.isCentralBody = True
    mu = earth.mu
    for sat in satellites:
        gravFactory.addBodiesTo(sat)
    oe = orbitalMotion.ClassicElements()
    oe.a = 10000000.0
    oe.e = 0.0
    oe.i = 0.0
    oe.Omega = 0.0
    oe.omega = 0.0
    for i, sat in enumerate(satellites):
        oe.f = i * (360.0 / num_satellites) * macros.D2R
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        sat.hub.r_CN_NInit = rN
        sat.hub.v_CN_NInit = vN
    scLocationModules = []
    for i, sat in enumerate(satellites):
        scLocation = spacecraftLocation.SpacecraftLocation()
        scLocation.ModelTag = f"CommCheck{i+1}"
        scLocation.primaryScStateInMsg.subscribeTo(sat.scStateOutMsg)
        for j, other_sat in enumerate(satellites):
            if i != j:
                scLocation.addSpacecraftToModel(other_sat.scStateOutMsg)
        scLocation.rEquator = earth.radEquator
        scLocation.maximumRange = 1000000.0
        scSim.AddModelToTask(simTaskName, scLocation)
        scLocationModules.append(scLocation)
    if vizSupport.vizFound:
        print("Vizard Visualization Found. Enabling live streaming...")
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)
        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, satellites, liveStream=True)
        for i, scLocation in enumerate(scLocationModules):
            vizSupport.addLocation(viz, stationName=f"Satellite{i+1}", parentBodyName=f"Satellite{i+1}",
                                   r_GP_P=[0, 0, 0], gHat_P=[0, 1, 0], fieldOfView=np.pi/4, range=scLocation.maximumRange)
    else:
        print("Vizard Visualization Module Not Found. Check Basilisk installation.")
    simulationTime = macros.min2nano(10.0)
    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()

if __name__ == "__main__":
    run(show_plots=True)
