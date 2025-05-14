"""
scenario_SatCameraTarget.py

Main Basilisk scenario: Simulates a spacecraft using a camera to track
and point at various ground targets. Integrates Vizard visualization,
target switching, and camera attitude control.
"""

# Get current file path
import inspect
import os
import sys
import numpy as np

from Basilisk import __path__  # Added for proper Basilisk path handling
from Basilisk.utilities import orbitalMotion, macros, unitTestSupport, vizSupport
from Basilisk.architecture import messaging
from Basilisk.simulation import spacecraft, extForceTorque, simpleNav
from Basilisk.fswAlgorithms import mrpFeedback, attTrackingError, locationPointing

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
bskPath = __path__[0]  # Get Basilisk path for SPICE data

# Import master classes: simulation base class and scenario base class
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')
from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

# Import the target_plotting module for visualization
from target_plotting import plot_target_visibility

def target_ground_location(currentTime, targetLat, targetLon, targetAlt, navObject, lastPrintTime):
    """
    Compute attitude pointing commands to point at a ground target.
    
    Args:
        currentTime: Current simulation time (nanoseconds)
        targetLat: Target latitude in degrees
        targetLon: Target longitude in degrees  
        targetAlt: Target altitude in meters
        navObject: Navigation object providing spacecraft state
        lastPrintTime: Time of last status print
        
    Returns:
        sigma_BR: MRP attitude command
        lastPrintTime: Updated last print time
    """
    # Convert target coordinates to radians
    lat = targetLat * macros.D2R
    lon = targetLon * macros.D2R
    
    # Get spacecraft position and velocity
    r_BN_N = navObject.r_BN_N
    v_BN_N = navObject.v_BN_N
    
    # Earth parameters
    r_Earth = 6378136.3  # Earth radius in meters
    
    # Compute ECEF coordinates of the target
    r_TP_P = np.array([
        [(r_Earth + targetAlt) * np.cos(lat) * np.cos(lon)],
        [(r_Earth + targetAlt) * np.cos(lat) * np.sin(lon)],
        [(r_Earth + targetAlt) * np.sin(lat)]
    ])
    
    # Compute Earth rotation at current time
    timeMin = currentTime * macros.NANO2MIN
    omega_E = 2 * np.pi / (23.93 * 60.0)  # Earth rotation rate in rad/min
    theta = omega_E * timeMin
    
    # Rotation matrix from ECEF to ECI
    DCM_NP = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
    
    # Target position in ECI
    r_TN_N = DCM_NP @ r_TP_P
    
    # Vector from spacecraft to target
    r_BT_N = r_TN_N - r_BN_N
    
    # Normalize
    r_BT_N_unit = r_BT_N / np.linalg.norm(r_BT_N)
    
    # Define camera boresight in body frame
    r_camera_B = np.array([[0], [1], [0]])  # Camera points along body Y-axis
    
    # Compute required rotation from body to ECI frame
    DCM_BN = np.eye(3)  # Start with identity
    DCM_BN[:, 1] = r_BT_N_unit.flatten()  # Second column is target direction
    
    # Orthogonalize the frame
    z_temp = np.cross(DCM_BN[:, 1], r_BN_N.flatten())
    z_temp = z_temp / np.linalg.norm(z_temp)
    x_temp = np.cross(z_temp, DCM_BN[:, 1])
    x_temp = x_temp / np.linalg.norm(x_temp)
    
    DCM_BN[:, 0] = x_temp
    DCM_BN[:, 2] = z_temp
    
    # Convert to MRP
    sigma_BN = unitTestSupport.dcm_to_mrp(DCM_BN)
    sigma_BR = np.array([[sigma_BN[0][0]], [sigma_BN[1][0]], [sigma_BN[2][0]]])
    
    # Print status every 5 minutes
    if timeMin - lastPrintTime >= 5:
        # Check if target is in view (rough check)
        dot_product = np.dot(r_BT_N_unit.flatten(), -r_BN_N.flatten() / np.linalg.norm(r_BN_N))
        visible = dot_product <= 0  # Negative dot product means target is on visible side of Earth
        
        print(f"Time: {timeMin:.1f} min | Target: {targetLat:.1f}°N, {targetLon:.1f}°E | "
              f"Distance: {np.linalg.norm(r_BT_N)/1000:.1f} km | Visible: {'Yes' if visible else 'No'}")
        lastPrintTime = timeMin
    
    return sigma_BR, lastPrintTime

def setup_vizard_markers(viz, targets):
    """Create basic markers using most compatible method"""
    if viz is None:
        return
        
    try:
        # Earth radius
        r_Earth = 6371000.0
        
        for i, target in enumerate(targets):
            lat = target["lat"] * macros.D2R
            lon = target["lon"] * macros.D2R
            
            # Convert lat/lon to Earth-centered Cartesian coordinates
            r_Target = [
                r_Earth * np.cos(lat) * np.cos(lon),
                r_Earth * np.cos(lat) * np.sin(lon),
                r_Earth * np.sin(lat)
            ]
            
            # Try the simplest marker creation method that works in most versions
            try:
                viz.addEarthLocation(target["name"], lat, lon, 0.0, colorIndex=i+1)
                print(f"Added marker for {target['name']}")
            except AttributeError:
                try:
                    viz.createPointMarker(target["name"], r_Target, colorIndex=i+1)
                    print(f"Added marker for {target['name']} using point marker")
                except:
                    print(f"Could not add marker for {target['name']} - incompatible Vizard version")
    except Exception as e:
        print(f"Warning: Could not add markers: {e}")

class scenario_SatCameraTarget(BSKSim, BSKScenario):
    def __init__(self, use_vizard=True):
        super().__init__()
        self.name = 'scenario_SatCameraTarget'
        self.use_vizard = use_vizard
        
        # Initialize all event flags to prevent attribute errors
        self.batteryDropFlag = False
        self.powerCutFlag = False
        self.thrusterFaultFlag = False
        self.sensorNoiseFlag = False
        self.commDelayFlag = False
        self.commOutageFlag = False
        self.attitudeErrorFlag = False
        self.reactionWheelFailureFlag = False
        self.orbitPerturbationFlag = False
        self.payloadActivationFlag = False
        self.dataCollectionFlag = False
        self.downlinkFlag = False
        self.eclipseFlag = False
        self.targetInViewFlag = False
        self.cssFaultFlag = False
        self.sunPointingErrorFlag = False
        self.starTrackerFaultFlag = False
        self.gpsOutageFlag = False
        self.payloadFaultFlag = False
        self.memoryCorruptionFlag = False
        self.propulsionFailureFlag = False
        self.orbitalDecayFlag = False
        self.radiationDamageFlag = False
        self.timelineEventFlag = False
        
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.targets = [
            {"name": "New York", "lat": 40.7128, "lon": -74.0060},
            {"name": "Los Angeles", "lat": 34.0522, "lon": -118.2437},
            {"name": "London", "lat": 51.5074, "lon": -0.1278}
        ]
        self.currentTargetIndex = 0
        self.lastTargetSwitchTime = 0
        self._last_target_print = -1

        # Set up dynamics and FSW models
        self.configure_initial_conditions()
        self.viz = None
        
        # Fix: Create the attGuidMsg before configuring dynamics
        # Create and initialize the attitude guidance message
        self.attGuidMsg = messaging.AttGuidMsg()
        self.attGuidMsg.sigma_BR = [[0.0], [0.0], [0.0]]
        self.attGuidMsg.omega_BR_B = [[0.0], [0.0], [0.0]]
        self.attGuidMsg.domega_BR_B = [[0.0], [0.0], [0.0]]
            
        # Configure additional dynamics and log outputs AFTER viz setup
        self.configure_dynamics()
        
        # Set up visualization if requested
        if use_vizard and vizSupport.vizFound:
            self.setup_visualization()
            
        self.log_outputs()

    def configure_dynamics(self):
        """Additional dynamics configuration - for more robust simulation"""
        # Create attitude tracking error module
        self.attErr = attTrackingError.attTrackingError()
        self.attErr.ModelTag = "attTrackingError"
        
        # Connect navigation output to attitude error module
        self.attErr.attNavInMsg.subscribeTo(self.DynModels.simpleNavObject.attOutMsg)
        
        # Add attitude control
        self.mrpCtrl = mrpFeedback.mrpFeedback()
        self.mrpCtrl.ModelTag = "mrpFeedback"
        self.mrpCtrl.Ki = -1.0
        self.mrpCtrl.P = 35.0
        self.mrpCtrl.K = 0.8
        
        # Create external torque effector
        self.extFT = extForceTorque.ExtForceTorque()
        self.extFT.ModelTag = "extForceTorque"
        
        # Fix: First, ensure we have a properly crafted attitude message
        # that can be connected to the attitude error module
        self.attRefMsg = messaging.AttRefMsg()
        self.attRefMsg.write(messaging.AttRefMsgPayload())
        
        # Connect the attitude error to the messages
        self.attErr.attRefInMsg.subscribeTo(self.attRefMsg)
        self.mrpCtrl.guidInMsg.subscribeTo(self.attErr.attGuidOutMsg)
        
        # Connect torque commands to external force/torque object
        self.extFT.cmdTorqueInMsg.subscribeTo(self.mrpCtrl.cmdTorqueOutMsg)
        self.DynModels.scObject.addDynamicEffector(self.extFT)
        
        # Set up vehicle configuration message
        vehConfig = messaging.VehicleConfigMsgPayload()
        I = [900., 0., 0., 0., 800., 0., 0., 0., 600.]
        vehConfig.ISCPntB_B = I
        self.vcMsg = messaging.VehicleConfigMsg().write(vehConfig)
        self.mrpCtrl.vehConfigInMsg.subscribeTo(self.vcMsg)
        
        # Add models to task
        self.AddModelToTask(self.DynModels.taskName, self.attErr)
        self.AddModelToTask(self.DynModels.taskName, self.mrpCtrl)
        self.AddModelToTask(self.DynModels.taskName, self.extFT)

    def setup_visualization(self):
        """Setup Vizard visualization with proper initialization"""
        try:
            # Ensure the viz directory exists
            viz_dir = os.path.join(path, "_VizFiles")
            if not os.path.exists(viz_dir):
                os.makedirs(viz_dir)

            # Enable visualization with basic settings
            self.viz = vizSupport.enableUnityVisualization(
                self, 
                self.DynModels.taskName, 
                self.DynModels.scObject,
                liveStream=True, 
                saveFile=os.path.join(viz_dir, "vizardOutput_UnityViz.bin")
            )
            
            if self.viz is not None:
                print("Successfully connected to Vizard")
                # Configure visualization settings
                self.viz.settings.orbitLinesOn = 1
                self.viz.settings.showMissionTime = 1
                
                # Create camera model for visualization
                camHUD = vizSupport.vizInterface.Transceiver()
                camHUD.r_SB_B = [0.0, 1.0, 0.0]  # Camera position in body frame
                camHUD.normalVector = [0.0, 1.0, 0.0]  # Camera pointing direction
                camHUD.fieldOfView = 5 * macros.D2R  # 5 degree FOV
                camHUD.label = "Camera"
                camHUD.transceiverState = 2  # Active state
                camHUD.color = vizSupport.vizInterface.IntVector([0, 255, 0, 255])  # Green
                
                # Add camera to visualization
                self.viz.addActorToScene(camHUD)
                
                # Add ground target markers
                setup_vizard_markers(self.viz, self.targets)
                
                # Create a standard camera view for visualization
                vizSupport.createStandardCamera(
                    viz=self.viz,
                    setMode=1,  # Spacecraft centered view
                    spacecraftName=self.DynModels.scObject.ModelTag,
                    fieldOfView=30 * macros.D2R,  # 30 degree FOV for viewing
                    displayName="SpacecraftView",
                    pointingVector_B=[0, 1, 0],
                    position_B=[0, -20, 0]
                )
            else:
                print("Failed to connect to Vizard - running without visualization")
                
        except Exception as e:
            print(f"Failed to connect to Vizard: {e}")
            print("Running simulation without visualization...")
            self.viz = None

    def configure_initial_conditions(self):
        """Set up the initial orbital elements and spacecraft state"""
        # Set up orbital elements
        oe = orbitalMotion.ClassicElements()
        oe.a = 7000000.0  # Semi-major axis [m] - slightly higher orbit for better visibility
        oe.e = 0.001      # Low eccentricity
        oe.i = 51.6 * macros.D2R  # Inclination (ISS-like) [rad]
        oe.Omega = 48.2 * macros.D2R  # Right ascension of ascending node [rad]
        oe.omega = 347.8 * macros.D2R  # Argument of periapsis [rad]
        oe.f = 85.3 * macros.D2R  # True anomaly [rad]

        # Convert to position/velocity
        mu = self.DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)

        # Set spacecraft initial state
        self.DynModels.scObject.hub.r_CN_NInit = rN
        self.DynModels.scObject.hub.v_CN_NInit = vN
        self.DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        self.DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]
        
        # Set spacecraft physical properties
        self.DynModels.scObject.hub.mHub = 750.0  # Spacecraft mass [kg]
        I = [900., 0., 0., 0., 800., 0., 0., 0., 600.]  # Inertia tensor [kg*m^2]
        self.DynModels.scObject.hub.IHubPntBc_B = unitTestSupport.np2EigenMatrix3d(I)

    def log_outputs(self):
        """Set up message logging for analysis"""
        samplingTime = self.get_FswModel().processTasksTimeStep
        self.sNavTransName = "sNavTransMsg"
        self.attRefName = "attRefMsg"

        self.msgRecList = {
            self.attRefName: self.attRefMsg.recorder(samplingTime),
            self.sNavTransName: self.get_DynModel().simpleNavObject.transOutMsg.recorder(samplingTime)
        }

        for rec in self.msgRecList.values():
            self.AddModelToTask(self.DynModels.taskName, rec)

    def process_tasks(self, time):
        """Process tasks at each simulation step"""
        # Check if it's time to switch targets
        time_min = time * macros.NANO2MIN
        if time_min - self.lastTargetSwitchTime >= 5.0:
            self.currentTargetIndex = (self.currentTargetIndex + 1) % len(self.targets)
            self.lastTargetSwitchTime = time_min
            print(f"Switched to: {self.targets[self.currentTargetIndex]['name']}")

        # Get current target
        target = self.targets[self.currentTargetIndex]
        
        # Compute pointing to target
        sigma_BR, self._last_target_print = target_ground_location(
            time, target["lat"], target["lon"], 0.0,
            self.DynModels.simpleNavObject, self._last_target_print
        )

        # Update attitude reference message with new pointing command
        attData = messaging.AttRefMsgPayload()
        attData.sigma_RN = sigma_BR  # Set the reference attitude
        attData.omega_RN_N = [[0.0], [0.0], [0.0]]
        attData.domega_RN_N = [[0.0], [0.0], [0.0]]
        
        # Write to the attitude reference message
        self.attRefMsg.write(attData)

        # Call parent process_tasks
        super().process_tasks(time)

    def pull_outputs(self, showPlots):
        """Process simulation outputs and generate plots"""
        # Get message recorders
        attRec = self.msgRecList[self.attRefName]
        posRec = self.msgRecList[self.sNavTransName]

        # Process data
        timeData = np.delete(attRec.times(), 0) * macros.NANO2MIN
        sigma_BR = np.delete(attRec.sigma_RN, 0, 0)  # Fixed to use sigma_RN from AttRefMsg
        r_BN_N = np.delete(posRec.r_BN_N, 0, 0)

        # Generate plots
        BSK_plt.clear_all_plots()
        BSK_plt.plot_attitude_error(timeData, sigma_BR, id=1)

        # Generate the target visibility plot
        vis_fig = plot_target_visibility(timeData, r_BN_N, self.targets)

        # Show or save plots
        if showPlots:
            BSK_plt.show_all_plots()
            vis_fig.show()
        else:
            fname = os.path.splitext(os.path.basename(__file__))[0]
            figureList = BSK_plt.save_all_plots(fname, ["attitudeErrorNorm"])
            vis_fig.savefig(fname + "_targetVisibility.png")
            figureList[fname + "_targetVisibility"] = vis_fig
            return figureList 