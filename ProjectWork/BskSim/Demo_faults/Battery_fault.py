#
#  ISC License
#
#  Copyright (c) 2016, Autonomous Vehicle Systems Lab, University of Colorado at Boulder
#
#  Permission to use, copy, modify, and/or distribute this software for any
#  purpose with or without fee is hereby granted, provided that the above
#  copyright notice and this permission notice appear in all copies.
#
#  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
#  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
#  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

"""
Overview
--------

This script duplicates the basic orbit simulation in the scenario :ref:`scenarioBasicOrbit`.
The difference is that this version allows for the Basilisk simulation data to be live streamed to the
:ref:`vizard` visualization program.

The script is found in the folder ``basilisk/examples`` and executed by using::

    python3 scenarioBasicOrbitStream.py

To enable live data streaming, the ``enableUnityVisualization()`` method is provided with ``liveStream``
argument using::

    vizSupport.enableUnityVisualization(scSim, simTaskName, scObject
                                        , liveStream=True)

When starting Basilisk simulation it prints now to the terminal that it is trying to connect to Vizard::

    Waiting for Vizard at tcp://localhost:5556

Copy ``tcp://localhost:5556`` and open the Vizard application.  Enter this address in the connection field and select
"Direct Communication" mode as well as "Live Streaming".  After this the Basilisk simulation resumes and
will live stream the data to Vizard.

.. figure:: /_images/static/vizard-ImgStream.jpg
   :align: center
   :scale: 50 %

   Vizard Direct Communication Panel Illustration


To avoid the simulation running too quickly, this tutorial example script includes the ``clock_sync`` module that
enables a 50x realtime mode using::

    clockSync = clock_synch.ClockSynch()
    clockSync.accelFactor = 50.0
    scSim.AddModelToTask(simTaskName, clockSync)

This way a 10s simulation time step will take 0.2 seconds with the 50x speed up factor.

"""


#
# Basilisk Scenario Script and Integrated Test
#
# Purpose:  Integrated test of the spacecraft() and gravity modules.  Illustrates
#           a 3-DOV spacecraft on a range of orbit types with live Vizard data streaming.
# Author:   Hanspeter Schaub
# Creation Date:  Sept. 29, 2019
#



import os
import sys
import inspect
import matplotlib.pyplot as plt
import numpy as np
# The path to the location of Basilisk
# Used to get the location of supporting data.
from Basilisk import __path__


filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

sys.path.append(path + '/../')
sys.path.append(path + '/../modelsMultiSat')
sys.path.append(path + '/../plottingMultiSat')



bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])

# import simulation related support
from Basilisk.simulation import spacecraft
# general support file with common unit test functions
# import general simulation support files
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion,
                                simIncludeGravBody, unitTestSupport, vizSupport)

from Basilisk.simulation import simSynch

from Basilisk.simulation import simpleBattery 
from Basilisk.architecture import messaging
from Basilisk.simulation import simpleSolarPanel
from Basilisk.simulation import eclipse
from Basilisk.simulation import simplePowerSink








def run(show_plots, liveStream, timeStep, orbitCase, useSphericalHarmonics, planetCase):
    

    # Create simulation variable names
    simTaskName = "simTask"
    simProcessName = "simProcess"

   

    #  Create a sim module as an empty container
    scSim = SimulationBaseClass.SimBaseClass()

    #
    #  create the simulation process
    #
    dynProcess = scSim.CreateNewProcess(simProcessName)

    # create the dynamics task and specify the integration update time
    simulationTimeStep = macros.sec2nano(timeStep)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

    #
    #   setup the simulation tasks/objects
    #

    # initialize spacecraft object and set properties
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = "bskSat"

    cameraLocation = [0.0, 2.0, 0.0]

    targets = [ {"name": "Tokyo",          "lat":  35.68,   "lon": 139.77,   "color": "green"},
                {"name": "Sri Lanka",      "lat":   6.9271, "lon":  79.8612, "color": "orange"},
                {"name": "Central Africa", "lat":   1.5333, "lon":  17.6667, "color": "purple"}
    ]


    
    # Create and configure the battery
    battery = simpleBattery.SimpleBattery()
    battery.ModelTag = "satBattery"
    battery.storageCapacity   = 100.0 
    battery.storedCharge_Init = 50.0 
    scSim.AddModelToTask(simTaskName, battery)
 
  

    # Power Consumption 
    # Primary power sink (10 W)
    # Create & add the power sink 
    powerSink = simplePowerSink.SimplePowerSink()
    powerSink.ModelTag    = "powerSink"
    powerSink.nodePowerOut = -0.01       # sink 10 W
    scSim.AddModelToTask(simTaskName, powerSink)
    # hook the sink into the battery
    battery.addPowerNodeToModel(powerSink.nodePowerOutMsg)

    # 2) Secondary (fault-injection) load
    powerSinkFault = simplePowerSink.SimplePowerSink()
    powerSinkFault.ModelTag    = "powerSinkFault"
    powerSinkFault.nodePowerOut = 0.0     # start off
    scSim.AddModelToTask(simTaskName, powerSinkFault)
    battery.addPowerNodeToModel(powerSinkFault.nodePowerOutMsg)
   


################################################################################################
  
    # expose sink so your event can see it
    scSim.powerSink = powerSink

    # expose for events
    scSim.powerSinkFault = powerSinkFault

    # 1) record the battery status
    battRec = battery.batPowerOutMsg.recorder(simulationTimeStep)
    # expose to scSim
    scSim.battRec = battRec
    scSim.AddModelToTask(simTaskName, battRec)

    faultTime = macros.min2nano(5.0)

    # 2) compute the 20% threshold
    safeThresh = battery.storageCapacity * 0.2

    # 3) create the safe‐mode event, only firing once there's at least one sample:
    cond = (
        f"(len(self.battRec.storageLevel) > 0) and "
        f"(self.battRec.storageLevel[-1] <= {safeThresh})"
    )
    scSim.createNewEvent(
        "batterySafeMode",          # event name
        simulationTimeStep,         # check interval
        True,                       # enabled immediately
        [cond],                     # guarded condition
        [
            "print('*** ENTERING SAFE MODE: SOC ≤ 20% ***')",
            "self.powerSink.nodePowerOut = 0.0",
            "self.powerSinkFault.nodePowerOut = 0.0",
            "self.setEventActivity('batterySafeMode', False)"
        ]
    )


    scSim.createNewEvent(
        "powerSinkFaultEvent",        
        simulationTimeStep,               # how often to check
        True,                   
        [f"self.TotalSim.CurrentNanos >= {faultTime}"],   # condition
        [
        # start drawing 10 W 
        "self.powerSinkFault.nodePowerOut = -0.05",
        # disable this event so it only fires once
        "self.setEventActivity('powerSinkFaultEvent', False)"
        ]
    )


#################################################################################################

    
    

    

    



    # Solar Panel
    solarPanel = simpleSolarPanel.SimpleSolarPanel()
    solarPanel.ModelTag = "solarPanel"
    solarPanel.setPanelParameters([-1.0, -10.0, -1.0], 0.00001, 0.0000001)
    solarPanel.stateInMsg.subscribeTo(scObject.scStateOutMsg)
    scSim.AddModelToTask(simTaskName, solarPanel)
    battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)  
    


    
    rawSun = np.array([-1.0, -10.0, -1.0])
    sunDir = (rawSun / np.linalg.norm(rawSun)).tolist()
    sunMsgData = messaging.SpicePlanetStateMsgPayload()
    sunMsgData.PositionVector = sunDir
    sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
    solarPanel.sunInMsg.subscribeTo(sunMsg)
    
  
    # 1) make & add the recorder
    battRec = battery.batPowerOutMsg.recorder(simulationTimeStep)
    scSim.AddModelToTask(simTaskName, battRec)

    # 2) list its public attrs
    print(">>> battRec attrs:", [a for a in dir(battRec) if not a.startswith("_")])

        
 

    # Visualization Setup 
    gsList = []

    
    # Battery visualization
    batteryPanel = vizSupport.vizInterface.GenericStorage()
    batteryPanel.label = "Battery (%)"
    batteryPanel.units = "%"
    batteryPanel.minValue = 0
    batteryPanel.maxValue = 100

    batteryPanel.useStorageLevel = True
    batteryInMsg = messaging.PowerStorageStatusMsgReader()
    batteryInMsg.subscribeTo(battery.batPowerOutMsg)
    batteryPanel.batteryStateInMsg = batteryInMsg
    batteryPanel.this.disown()

    batteryPanel.thresholds = vizSupport.vizInterface.IntVector([20, 50, 80])

    batteryPanel.color = vizSupport.vizInterface.IntVector(
        vizSupport.toRGBA255("red") +
        vizSupport.toRGBA255("orange") +
        vizSupport.toRGBA255("yellow") +
        vizSupport.toRGBA255("green")
    )

    batteryPanel.thresholds = vizSupport.vizInterface.IntVector([0])
    batteryPanel.color      = vizSupport.vizInterface.IntVector(
        vizSupport.toRGBA255("red") +   # SOC ≤ 20% → red
        vizSupport.toRGBA255("green")   # SOC > 20% → green
    )
    



    solarViz = vizSupport.vizInterface.GenericStorage()
    solarViz.label           = "Solar Power"
    solarViz.units           = "W"
    solarViz.minValue        = 0.0
    solarViz.maxValue        = 20.0    # set to a bit above your expected peak

    solarViz.useStorageLevel = False  
    solarReader = messaging.PowerNodeUsageMsgReader()  # Change to match the correct message type
    solarReader.subscribeTo(solarPanel.nodePowerOutMsg)
    solarViz.powerNodeUsageInMsg = solarReader
    solarViz.this.disown()
    scSim.solarReader = solarReader


    sinkReader = messaging.PowerNodeUsageMsgReader()
    sinkReader.subscribeTo(powerSink.nodePowerOutMsg)
    scSim.sinkReader = sinkReader   

    sinkViz = vizSupport.vizInterface.GenericStorage()
    sinkViz.label           = "Power Sink"
    sinkViz.units           = "W"
    sinkViz.minValue        = -0.15
    sinkViz.maxValue        =  0.0
    sinkViz.useStorageLevel = True
    sinkViz.powerNodeUsageInMsg = sinkReader
    sinkViz.this.disown()
        


    scSim.solarReader = solarReader
    scSim.sinkReader  = sinkReader

    solarViz.powerNodeUsageInMsg = solarReader
    sinkViz .powerNodeUsageInMsg = sinkReader

    # now show all three panels on the same row
    gsList.append([batteryPanel])
   
    

   
  
 
    ##################################################### 


    # add spacecraft object to the simulation process
    scSim.AddModelToTask(simTaskName, scObject)

    # setup Gravity Body
    gravFactory = simIncludeGravBody.gravBodyFactory()
    if planetCase == 'Mars':
        planet = gravFactory.createMarsBarycenter()
        planet.isCentralBody = True           # ensure this is the central gravitational body
        if useSphericalHarmonics:
            planet.useSphericalHarmonicsGravityModel(bskPath + '/supportData/LocalGravData/GGM2BData.txt', 100)

    else:  # Earth
        planet = gravFactory.createEarth()
        planet.isCentralBody = True          # ensure this is the central gravitational body
        if useSphericalHarmonics:
            planet.useSphericalHarmonicsGravityModel(bskPath + '/supportData/LocalGravData/GGM03S-J2-only.txt', 2)
    mu = planet.mu

    # attach gravity model to spacecraft
    gravFactory.addBodiesTo(scObject)

    #
    #   setup orbit and simulation time
    #
    # setup the orbit using classical orbit elements
    oe = orbitalMotion.ClassicElements()
    rLEO = 7000. * 1000      # meters
    rGEO = 42000. * 1000     # meters
    if orbitCase == 'GEO':
        oe.a = rGEO
        oe.e = 0.00001
        oe.i = 0.0 * macros.D2R
    elif orbitCase == 'GTO':
        oe.a = (rLEO + rGEO) / 2.0
        oe.e = 1.0 - rLEO / oe.a
        oe.i = 0.0 * macros.D2R
    else:                   # LEO case, default case 0
        oe.a = rLEO
        oe.e = 0.0001
        oe.i = 33.3 * macros.D2R
    oe.Omega = 48.2 * macros.D2R
    oe.omega = 347.8 * macros.D2R
    oe.f = 85.3 * macros.D2R
    rN, vN = orbitalMotion.elem2rv(mu, oe)
    oe = orbitalMotion.rv2elem(mu, rN, vN)      # this stores consistent initial orbit elements
    # with circular or equatorial orbit, some angles are arbitrary

    #
    #   initialize Spacecraft States with the initialization variables
    #
    scObject.hub.r_CN_NInit = rN  # m   - r_BN_N
    scObject.hub.v_CN_NInit = vN  # m/s - v_BN_N

    # set the simulation time
    n = np.sqrt(mu / oe.a / oe.a / oe.a)
    P = 2. * np.pi / n
    
    if useSphericalHarmonics:
        simulationTime = macros.sec2nano(2. * P)
    else:
        simulationTime = macros.sec2nano(2 * P)

    #
    #   Setup data logging before the simulation is initialized
    #
    if useSphericalHarmonics:

        numDataPoints = 400
    else:
        numDataPoints = 100
    samplingTime = unitTestSupport.samplingTime(simulationTime, simulationTimeStep, numDataPoints)
    dataLog = scObject.scStateOutMsg.recorder(samplingTime)
    scSim.AddModelToTask(simTaskName, dataLog)

   
    # record the battery’s storage-level and net-power messages every time step
    batteryLog = battery.batPowerOutMsg.recorder(simulationTimeStep)
    scSim.AddModelToTask(simTaskName, batteryLog)

    

    if liveStream:
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)

        # if this scenario is to interface with the BSK Viz, uncomment the following line
        viz=vizSupport.enableUnityVisualization(scSim, simTaskName, scObject
                                            , liveStream=True
                                            , genericStorageList=gsList
                                           
                                            ) 


    
        vizSupport.setInstrumentGuiSetting(viz, 
                                            spacecraftName=scObject.ModelTag,
                                            showGenericStoragePanel=True)
        
        
        
        

        bodyName = planetCase.lower()   # 'earth' or 'mars'
        earthRadius = 6371000.0

        for tgt in targets:
            lat_r = tgt["lat"] * macros.D2R
            lon_r = tgt["lon"] * macros.D2R
            x = earthRadius * np.cos(lat_r) * np.cos(lon_r)
            y = earthRadius * np.cos(lat_r) * np.sin(lon_r)
            z = earthRadius * np.sin(lat_r)
            location_position = [x, y, z]

            vizSupport.addLocation(
                viz,
                stationName   = tgt["name"],
                parentBodyName= bodyName,
                r_GP_P        = location_position,
                color         = tgt["color"],
                range         = 2000e3 
            )
        

        # mount one camera on your sat 
        vizSupport.createStandardCamera(
            viz,
            setMode         = 1,                # 1 → spacecraft-attached mode
            spacecraftName  = scObject.ModelTag,
            displayName     = "BatteryCam",
            fieldOfView     = 30 * macros.D2R,
            pointingVector_B= [0.0, 0.0, 0.0],  # look forward in +X
            position_B      = cameraLocation
           
        )
        
        
      



    #
    #   initialize Simulation:  This function clears the simulation log, and runs the self_init()
    #   cross_init() and reset() routines on each module.
    #   If the routine InitializeSimulationAndDiscover() is run instead of InitializeSimulation(),
    #   then the all messages are auto-discovered that are shared across different BSK threads.
    #
    # Force message initialization


    scSim.InitializeSimulation()

    #
    #   configure a simulation stop time and execute the simulation run
    #
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()

    # debug: print the raw panel output



    #
    #   retrieve the logged data
    #
    posData = dataLog.r_BN_N
    velData = dataLog.v_BN_N

    np.set_printoptions(precision=16)

    #
    #   plot the results
    #
    # draw the inertial position vector components
    plt.close("all")  # clears out plots from earlier test runs
    plt.figure(1)
    fig = plt.gcf()
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False, style='plain')
    for idx in range(3):
        plt.plot(dataLog.times() * macros.NANO2SEC / P, posData[:, idx] / 1000.,
                 color=unitTestSupport.getLineColor(idx, 3),
                 label='$r_{BN,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [orbits]')
    plt.ylabel('Inertial Position [km]')
    figureList = {}
    pltName = fileName + "1" + orbitCase + str(int(useSphericalHarmonics))+ planetCase
    figureList[pltName] = plt.figure(1)

    if useSphericalHarmonics is False:
        # draw orbit in perifocal frame
        b = oe.a * np.sqrt(1 - oe.e * oe.e)
        p = oe.a * (1 - oe.e * oe.e)
        plt.figure(2, figsize=np.array((1.0, b / oe.a)) * 4.75, dpi=100)
        plt.axis(np.array([-oe.rApoap, oe.rPeriap, -b, b]) / 1000 * 1.25)
        # draw the planet
        fig = plt.gcf()
        ax = fig.gca()
        if planetCase == 'Mars':
            planetColor = '#884400'
        else:
            planetColor = '#008800'
        planetRadius = planet.radEquator / 1000
        ax.add_artist(plt.Circle((0, 0), planetRadius, color=planetColor))
        # draw the actual orbit
        rData = []
        fData = []
        for idx in range(0, len(posData)):
            oeData = orbitalMotion.rv2elem(mu, posData[idx], velData[idx])
            rData.append(oeData.rmag)
            fData.append(oeData.f + oeData.omega - oe.omega)
        plt.plot(rData * np.cos(fData) / 1000, rData * np.sin(fData) / 1000, color='#aa0000', linewidth=3.0
                 )
        # draw the full osculating orbit from the initial conditions
        fData = np.linspace(0, 2 * np.pi, 100)
        rData = []
        for idx in range(0, len(fData)):
            rData.append(p / (1 + oe.e * np.cos(fData[idx])))
        plt.plot(rData * np.cos(fData) / 1000, rData * np.sin(fData) / 1000, '--', color='#555555'
                 )
        plt.xlabel('$i_e$ Cord. [km]')
        plt.ylabel('$i_p$ Cord. [km]')
        plt.grid()

    else:
        plt.figure(2)
        fig = plt.gcf()
        ax = fig.gca()
        ax.ticklabel_format(useOffset=False, style='plain')
        smaData = []
        for idx in range(0, len(posData)):
            oeData = orbitalMotion.rv2elem(mu, posData[idx], velData[idx])
            smaData.append(oeData.a / 1000.)
        plt.plot(posData[:, 0] * macros.NANO2SEC / P, smaData, color='#aa0000',
                 )
        plt.xlabel('Time [orbits]')
        plt.ylabel('SMA [km]')

    pltName = fileName + "2" + orbitCase + str(int(useSphericalHarmonics)) + planetCase
    figureList[pltName] = plt.figure(2)

    # Retrieve the battery log you set up earlier
    timeData = batteryLog.times() * macros.NANO2SEC        
    storageData = batteryLog.storageLevel                  
    
    # Plot storage level
    plt.figure()
    plt.plot(timeData, storageData, label='Battery Stored Charge')
    
    # Mark the fault injection moment
    faultTime = macros.min2nano(60.0) * macros.NANO2SEC     
    plt.axvline(x=faultTime, color='r', linestyle='--', label='Fault Injected (–50 W sink)')
    
    plt.xlabel('Time [s]')
    plt.ylabel('Stored Charge [Wh]')
    plt.title('Battery Storage Level with Fault Injection')
    plt.legend()
    plt.grid()

    if show_plots:
        plt.show()

    # close the plots being saved off to avoid over-writing old and new figures
    plt.close("all")

    return figureList

#
# This statement below ensures that the unit test scrip can be run as a
# stand-along python script
#
if __name__ == "__main__":
    run(
        True,        # show_plots
        True,        # liveStream
        1.0,         # time step (s)
        'LEO',       # orbit Case (LEO, GTO, GEO)
        False,       # useSphericalHarmonics
        'Earth'      # planetCase (Earth, Mars)
    )

