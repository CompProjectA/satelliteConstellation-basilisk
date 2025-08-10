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
    scSim = SimulationBaseClass.SimBaseClass()

    #
    #  create the simulation process
    #
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(timeStep)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
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
   
    #for camara
    batteryReader = messaging.PowerStorageStatusMsgReader()
    batteryReader.subscribeTo(battery.batPowerOutMsg)      # listen for battery status
    scSim.batteryReader = batteryReader  
    
    print(dir(battery))
    powerSink = simplePowerSink.SimplePowerSink()
    powerSink.ModelTag    = "powerSink"
    powerSink.nodePowerOut = -0.01       # sink 10 W
    scSim.AddModelToTask(simTaskName, powerSink)
    # hook the sink into the battery
    battery.addPowerNodeToModel(powerSink.nodePowerOutMsg)


    scSim.powerSink = powerSink
    # compute the 1 min 
    faultTime = macros.min2nano(60.0)

    scSim.createNewEvent(
        "powerSinkFault",        
        simulationTimeStep,               # how often to check
        True,                   
        [f"self.TotalSim.CurrentNanos >= {faultTime}"],   # condition
        [
        # start drawing 10 W 
        "self.powerSink.nodePowerOut = -0.05",
        # disable this event so it only fires once
        "self.setEventActivity('powerSinkFault', False)"
        ]
    )
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
    scSim.AddModelToTask(simTaskName, battery)
    
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

    

    solarViz = vizSupport.vizInterface.GenericStorage()
    solarViz.label           = "Solar Power"
    solarViz.units           = "W"
    solarViz.minValue        = 0.0
    solarViz.maxValue        = 20.0    # set to a bit above your expected peak
    solarViz.useStorageLevel = False  # raw watts

    
    solarReader = messaging.PowerNodeUsageMsgReader()
    solarReader.subscribeTo(solarPanel.nodePowerOutMsg)

    
    solarViz.storageUnitStateInMsg = solarReader
    solarViz.this.disown()


   
    

    gsList.append([batteryPanel])
 


    # add spacecraft object to the simulation process
    scSim.AddModelToTask(simTaskName, scObject)

    # setup Gravity Body
    gravFactory = simIncludeGravBody.gravBodyFactory()
   
    planet = gravFactory.createEarth()
    planet.isCentralBody = True          # ensure this is the central gravitational body
    mu = planet.mu

    # attach gravity model to spacecraft
    gravFactory.addBodiesTo(scObject)
    oe = orbitalMotion.ClassicElements()
    rLEO = 7000. * 1000      # meters
    rGEO = 42000. * 1000     # meters
                     # LEO case, default case 0
    oe.a = (6371e3 + 1700e3)
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
    
   
    simulationTime = macros.sec2nano(3 * P)

    #
    #   Setup data logging before the simulation is initialized
    #
  
    numDataPoints = 100
    samplingTime = unitTestSupport.samplingTime(simulationTime, simulationTimeStep, numDataPoints)
    dataLog = scObject.scStateOutMsg.recorder(samplingTime)
    scSim.AddModelToTask(simTaskName, dataLog)

    # Battery state logger
    batteryLog = battery.batPowerOutMsg.recorder(samplingTime)
    scSim.AddModelToTask(simTaskName, batteryLog)

    

    if liveStream:
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)

        # if this scenario is to interface with the BSK Viz, uncomment the following line
        viz=vizSupport.enableUnityVisualization(scSim, simTaskName, scObject
                                            , liveStream=True
                                            , genericStorageList=gsList
                                            , saveFile=fileName
                                            ) 


    
        vizSupport.setInstrumentGuiSetting(viz, 
                                            spacecraftName=scObject.ModelTag,
                                            showGenericStoragePanel=True)
        

        

        for tgt in targets:
            lat = tgt["lat"]
            lon = tgt["lon"]
            color = tgt.get("color", "red")
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
                stationName=tgt["name"],
                parentBodyName="earth",
                r_GP_P=location_position,
                color=color
            )

        # mount one camera on your sat 
        vizSupport.createStandardCamera(
            viz,
            setMode         = 1,                # 1 → spacecraft-attached mode
            spacecraftName  = scObject.ModelTag,
            displayName     = "BatteryCam",
            fieldOfView     = 30 * macros.D2R,
            pointingVector_B= [0.0, 0.0, 0.0],  # look forward in +X
            position_B      = cameraLocation,
            
           
        )
        
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

