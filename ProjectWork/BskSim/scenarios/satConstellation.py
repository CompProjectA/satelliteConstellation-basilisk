import os
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.simulation import spacecraft, spacecraftLocation, simSynch
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion, simIncludeGravBody, unitTestSupport, vizSupport)
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

    num_satellites = 4
    satellites = []
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
    rLEO = 7000. * 1000
    oe.a = rLEO
    oe.e = 0.0001
    oe.i = 33.3 * macros.D2R
    oe.Omega = 48.2 * macros.D2R
    oe.omega = 347.8 * macros.D2R

    for i, sat in enumerate(satellites):
        oe.f = i * (360.0 / num_satellites) * macros.D2R
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        sat.hub.r_CN_NInit = rN
        sat.hub.v_CN_NInit = vN

    # Set up communication check modules
    scLocationModules = []
    accessRecorders = []

    for i, sat in enumerate(satellites):
        scLocation = spacecraftLocation.SpacecraftLocation()
        scLocation.ModelTag = f"CommCheck{i+1}"
        scLocation.primaryScStateInMsg.subscribeTo(sat.scStateOutMsg)

        # Add other satellites as communication targets
        for j, other_sat in enumerate(satellites):
            if i != j:
                scLocation.addSpacecraftToModel(other_sat.scStateOutMsg)

        
        scLocation.rEquator = earth.radEquator
        scLocation.rPolar = earth.radEquator * 0.98
        scLocation.aHat_B = [1.0, 0.0, 0.0] 
        scLocation.theta = np.radians(180.0)  
        scLocation.maximumRange = 2e7  # Larger range
        scSim.AddModelToTask(simTaskName, scLocation)
        scLocationModules.append(scLocation)

        # Record access results
        recorder = scLocation.scStateInMsgs[0].recorder()
        accessRecorders.append(recorder)
        scSim.AddModelToTask(simTaskName, recorder)

    # Enable visualization if available
    if vizSupport.vizFound:
        print("Vizard Visualization Found. Enabling live streaming...")
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)

        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, satellites, liveStream=True)
        for i, scLocation in enumerate(scLocationModules):
            vizSupport.addLocation(
                viz,
                stationName=f"Satellite{i+1}",
                parentBodyName=f"Satellite{i+1}",
                r_GP_P=[0, 0, 0],
                gHat_P=[0, 1, 0],
                fieldOfView=np.pi / 4,
                range=scLocation.maximumRange
            )
    else:
        print("Vizard Visualization Module Not Found. Check Basilisk installation.")

    simulationTime = macros.min2nano(10.0)
    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()

    # Plotting access for each satellite
    if show_plots:
        time_vec = accessRecorders[0].times() * macros.NANO2MIN
        plt.figure(figsize=(10, 6))

        for i in range(num_satellites):
            for j in range(num_satellites - 1):
                access = accessRecorders[i].r_BN_N
                plt.plot(time_vec, access, label=f'Sat{i+1} to Sat{j+1 if j < i else j+2}')

        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Satellites')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    run(show_plots=True)