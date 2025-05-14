import os, inspect
import numpy as np
from Basilisk import __path__
from Basilisk.simulation import spacecraft, extForceTorque, simpleNav, ephemerisConverter
from Basilisk.utilities import SimulationBaseClass, macros, simIncludeGravBody, unitTestSupport, orbitalMotion
from Basilisk.architecture import messaging
from Basilisk.fswAlgorithms import mrpFeedback, attTrackingError, locationPointing
from Basilisk.utilities import vizSupport
from Basilisk.simulation import vizInterface

def run():
    # Setup
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path = os.path.dirname(os.path.abspath(filename))
    bskPath = __path__[0]

    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    scSim.SetProgressBar(True)
    dynProcess = scSim.CreateNewProcess(simProcessName)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, macros.sec2nano(1.0)))

    # Spacecraft
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = "EarthImager"

    # Gravity using SPICE
    gravFactory = simIncludeGravBody.gravBodyFactory()
    gravFactory.createBodies(['sun', 'earth', 'moon'])
    earth = gravFactory.gravBodies['earth']
    earth.isCentralBody = True
    spice = gravFactory.createSpiceInterface(bskPath + "/supportData/EphemerisData/",
                                             "2021 MAY 04 07:47:00.0",
                                             epochInMsg=True)
    spice.zeroBase = 'Earth'
    gravFactory.addBodiesTo(scObject)

    mu = earth.mu

    # Orbit
    oe = orbitalMotion.ClassicElements()
    oe.a = 6778000.0
    oe.e = 0.0001
    oe.i = 51.6 * macros.D2R
    oe.Omega = 48.2 * macros.D2R
    oe.omega = 347.8 * macros.D2R
    oe.f = 85.3 * macros.D2R
    rN, vN = orbitalMotion.elem2rv(mu, oe)
    scObject.hub.r_CN_NInit = rN
    scObject.hub.v_CN_NInit = vN
    scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
    scObject.hub.omega_BN_BInit = [[0.0], [0.0], [0.0]]
    scObject.hub.mHub = 750.0
    I = [900., 0., 0., 0., 800., 0., 0., 0., 600.]
    scObject.hub.IHubPntBc_B = unitTestSupport.np2EigenMatrix3d(I)

    # Dynamics modules
    extFT = extForceTorque.ExtForceTorque()
    scObject.addDynamicEffector(extFT)
    sNav = simpleNav.SimpleNav()
    sNav.scStateInMsg.subscribeTo(scObject.scStateOutMsg)

    scSim.AddModelToTask(simTaskName, scObject)
    scSim.AddModelToTask(simTaskName, extFT)
    scSim.AddModelToTask(simTaskName, sNav)
    scSim.AddModelToTask(simTaskName, spice)

    # Earth ephemeris converter
    ephem = ephemerisConverter.EphemerisConverter()
    ephem.ModelTag = "earthEphem"
    ephem.addSpiceInputMsg(spice.planetStateOutMsgs[0])  # earth
    scSim.AddModelToTask(simTaskName, ephem)

    # Location pointing (New York)
    targetLat, targetLon, targetAlt = 40.7128, -74.0060, 0.0
    locPoint = locationPointing.locationPointing()
    locPoint.ModelTag = "locationPointing"
    locPoint.useLLA = 1
    locPoint.r_LP_P_Init = [targetLat * macros.D2R, targetLon * macros.D2R, targetAlt]
    locPoint.planetRadius = 6378136.3
    locPoint.pHat_B = [0, 1, 0]
    locPoint.gsHat_B = [0, 0, 1]
    locPoint.scTransInMsg.subscribeTo(sNav.transOutMsg)
    locPoint.scAttInMsg.subscribeTo(sNav.attOutMsg)
    locPoint.celBodyInMsg.subscribeTo(ephem.ephemOutMsgs[0])
    scSim.AddModelToTask(simTaskName, locPoint)

    # Attitude tracking and control
    attErr = attTrackingError.attTrackingError()
    attErr.attRefInMsg.subscribeTo(locPoint.attRefOutMsg)
    attErr.attNavInMsg.subscribeTo(sNav.attOutMsg)
    scSim.AddModelToTask(simTaskName, attErr)

    vehConfig = messaging.VehicleConfigMsgPayload()
    vehConfig.ISCPntB_B = I
    vcMsg = messaging.VehicleConfigMsg().write(vehConfig)

    mrpCtrl = mrpFeedback.mrpFeedback()
    mrpCtrl.ModelTag = "mrpFeedback"
    mrpCtrl.guidInMsg.subscribeTo(attErr.attGuidOutMsg)
    mrpCtrl.vehConfigInMsg.subscribeTo(vcMsg)
    mrpCtrl.Ki = -1.0
    mrpCtrl.P = 35.0
    mrpCtrl.K = 0.8
    scSim.AddModelToTask(simTaskName, mrpCtrl)
    extFT.cmdTorqueInMsg.subscribeTo(mrpCtrl.cmdTorqueOutMsg)

    # Visualization
    if vizSupport.vizFound:
        camHUD = vizInterface.Transceiver()
        camHUD.r_SB_B = [0.0, 1.0, 0.0]
        camHUD.normalVector = [0.0, 1.0, 0.0]
        camHUD.fieldOfView = 5 * macros.D2R
        camHUD.label = "Camera"
        camHUD.transceiverState = 2
        camHUD.color = vizInterface.IntVector([0, 255, 0, 255])

        viz = vizSupport.enableUnityVisualization(scSim, simTaskName, scObject,
                                                  transceiverList=[[camHUD]],
                                                  liveStream=True)
        
        viz.epochInMsg.subscribeTo(gravFactory.epochMsg)
        viz.settings.orbitLinesOn = 1
        viz.settings.showMissionTime = 1
        vizSupport.createStandardCamera(
                                        viz=viz,
                                        setMode=1,
                                        spacecraftName=scObject.ModelTag,
                                        fieldOfView=5 * macros.D2R,
                                        displayName="EarthImagerCam",
                                        pointingVector_B=[0, 1, 0],
                                        position_B=[0, 1, 0]
                                    )

        try:
            viz.addEarthLocation("New York", np.radians(targetLat), np.radians(targetLon), targetAlt, colorIndex=2)
        except AttributeError:
            print("Note: Earth location marker not supported in your Basilisk version.")

    # Run simulation
    scSim.InitializeSimulation()
    simTime = macros.sec2nano(4 * 60)
    scSim.ConfigureStopTime(simTime)
    print("🛰️  Running simulation for 4 minutes...")
    scSim.ExecuteSimulation()
    print("✅ Done.")

if __name__ == "__main__":
    run()
