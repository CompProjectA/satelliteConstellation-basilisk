"""
scenario_SatCameraTarget.py

This script demonstrates a spacecraft targeting ground locations with a camera.
It includes visualization, attitude control, target visibility analysis, 
and reaction wheel fault modeling.
"""

# Get current file path
import inspect
import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

from Basilisk.utilities import orbitalMotion, macros, unitTestSupport, vizSupport
from Basilisk.architecture import messaging
from Basilisk.simulation import spacecraft  # Import here for possible marker creation

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

# Import master classes: simulation base class and scenario base class
sys.path.append(path + '/../')
sys.path.append(path + '/../models')
sys.path.append(path + '/../plotting')
from BSK_masters import BSKSim, BSKScenario
import BSK_Dynamics, BSK_Fsw
import BSK_Plotting as BSK_plt

# Import custom camera utilities if available (otherwise use local functions)
try:
    from camera_utils import setup_camera, create_ground_target_markers, target_ground_location, update_camera_fov
    USING_CAMERA_UTILS = True
    print("Using imported camera_utils module")
except ImportError:
    USING_CAMERA_UTILS = False
    print("Using local camera functions")

# Create your own scenario child class
class scenario_SatCameraTarget(BSKSim, BSKScenario):
    def __init__(self, use_vizard=True):
        super(scenario_SatCameraTarget, self).__init__()
        self.name = 'scenario_SatCameraTarget'
        self.use_vizard = use_vizard

        # Declare message recorders
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        # Set dynamics and FSW models
        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        # Define target locations (lat, lon in degrees)
        self.targets = [
            {"name": "Target1", "lat": 40.7128, "lon": -74.0060},  # New York
            {"name": "Target2", "lat": 34.0522, "lon": -118.2437}, # Los Angeles
            {"name": "Target3", "lat": 51.5074, "lon": -0.1278}    # London
        ]
        self.currentTargetIndex = 0
        self.lastTargetSwitchTime = 0
        self._last_target_print = 0

        # Configure initial spacecraft state
        self.configure_initial_conditions()

        # Log relevant output messages
        self.log_outputs()

        # Initialize all GUI fault flags
        self.oneTimeRWFaultFlag = 1
        self.repeatRWFaultFlag = 1
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.batteryDropFlag = 0          
        self.powerCutFlag = 0              
        self.thrusterFaultFlag = 0        
        self.sensorNoiseFlag = 0          
        self.cssFaultFlag = 0             
        self.cssFaultIndex = -1           

        # Initialize fault logs
        self.DynModels.RWFaultLog = []

        # Vizard camera support
        self.viz = None
        self.camera = None
        if use_vizard and vizSupport.vizFound:
            try:
                print("Attempting to connect to Vizard...")
                self.viz = vizSupport.enableUnityVisualization(
                    self,
                    self.DynModels.taskName,
                    self.DynModels.scObject,
                    liveStream=True
                )
                
                # Configure visualization settings
                if hasattr(self.viz, 'settings'):
                    self.viz.settings.orbitLinesOn = 1
                    self.viz.settings.showSpacecraftLabels = 1
                    # These additional settings may or may not be available depending on Vizard version
                    if hasattr(self.viz.settings, 'enableEarthGrid'):
                        self.viz.settings.enableEarthGrid = 1
                    if hasattr(self.viz.settings, 'showOrbitPlane'):
                        self.viz.settings.showOrbitPlane = 1

                # Create and configure camera
                self.camera = self.setup_camera()
                
                # Create target markers in visualization
                self.createGroundTargetMarkers()
                
                print("Successfully connected to Vizard!")
            except Exception as e:
                print(f"WARNING: Failed to connect to Vizard: {e}")
                print("Continuing simulation without visualization...")

    def createGroundTargetMarkers(self):
        """Create visual markers for all ground targets using vizInterface"""
        if not hasattr(self, 'viz') or self.viz is None:
            return
                
        try:
            # Use a safer approach to create target markers with just basic vizSupport functionality
            for i, target in enumerate(self.targets):
                lat = target["lat"]
                lon = target["lon"]
                name = target["name"]
                
                # Create a model that represents the target - simple sphere
                try:
                    # Try to use any available methods in vizSupport for Earth locations
                    if hasattr(vizSupport, 'createTargetMarker'):
                        vizSupport.createTargetMarker(
                            self.viz, 
                            name,
                            lat * macros.D2R, 
                            lon * macros.D2R, 
                            0.0,
                            color=[1.0, 0.3, 0.3],
                            scale=150000.0
                        )
                        print(f"Created marker for {name} using vizSupport.createTargetMarker")
                        continue
                    
                    # If the viz has a method to add a shape to Earth, try that
                    if hasattr(self.viz, 'addShape'):
                        # Convert lat/lon to cartesian coordinates
                        lat_rad = lat * macros.D2R
                        lon_rad = lon * macros.D2R
                        r_Earth = 6371000.0  # Earth radius in meters
                        
                        # Create position vector in Earth-fixed frame
                        pos = [
                            r_Earth * np.cos(lat_rad) * np.cos(lon_rad),
                            r_Earth * np.cos(lat_rad) * np.sin(lon_rad),
                            r_Earth * np.sin(lat_rad)
                        ]
                        
                        self.viz.addShape(
                            f"Target_{name}", 
                            "Earth", 
                            "sphere.obj", 
                            pos, 
                            [0.0, 0.0, 0.0], 
                            150000.0, 
                            [1.0, 0.0, 0.0, 0.8]
                        )
                        print(f"Created marker for {name} using viz.addShape")
                        continue
                        
                    # Otherwise simply log the target info - we'll have to see it in text form
                    print(f"No available method to create visual marker for {name}")
                    
                except Exception as e:
                    print(f"Failed to create marker for {name}: {e}")
            
        except Exception as e:
            print(f"Target marker creation failed: {e}")
            
        # Always log the targets for reference
        print("\nGround Targets:")
        for target in self.targets:
            print(f"  {target['name']}: Lat {target['lat']:.4f}°, Lon {target['lon']:.4f}°")

    def configure_initial_conditions(self):
        # Configure Dynamics initial conditions
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0  # meters
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModels = self.get_DynModel()
        mu = DynModels.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        orbitalMotion.rv2elem(mu, rN, vN)
        DynModels.scObject.hub.r_CN_NInit = rN  # m   - r_CN_N
        DynModels.scObject.hub.v_CN_NInit = vN  # m/s - v_CN_N
        DynModels.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]  # sigma_BN_B
        DynModels.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]  # rad/s - omega_BN_B

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

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])

        return
    
    def target_ground_location(self, time, lat, lon, alt=0.0):
        """
        Create attitude command to point the camera at a specific ground location
        
        Args:
            time: Current simulation time (ns)
            lat: Latitude in degrees
            lon: Longitude in degrees
            alt: Altitude above sea level in meters (default: 0.0)
        
        Returns:
            sigma_BR: MRP attitude error
        """
        DynModel = self.get_DynModel()
        
        # Get current spacecraft position from navigation
        r_BN_N = DynModel.simpleNavObject.transOutMsg.read().r_BN_N
        
        # Convert lat/lon to ECEF coordinates (Earth fixed frame)
        lat_rad = lat * macros.D2R
        lon_rad = lon * macros.D2R
        r_Earth = 6371000.0  # Earth radius in meters
        
        # Target position in Earth-fixed frame
        r_Target_E = np.array([
            [(r_Earth + alt) * np.cos(lat_rad) * np.cos(lon_rad)], 
            [(r_Earth + alt) * np.cos(lat_rad) * np.sin(lon_rad)], 
            [(r_Earth + alt) * np.sin(lat_rad)]
        ])
        
        # Convert Earth-fixed to inertial
        time_hours = time * macros.NANO2HOUR
        earth_rotation = time_hours * 15.0 * macros.D2R  # Earth rotates ~15 degrees per hour
        
        c_rot = np.cos(earth_rotation)
        s_rot = np.sin(earth_rotation)
        
        R_N_E = np.array([
            [c_rot, -s_rot, 0],
            [s_rot, c_rot, 0],
            [0, 0, 1]
        ])
        
        r_Target_N = R_N_E @ r_Target_E
        
        # Vector from spacecraft to target in inertial frame
        r_BT_N = r_Target_N - np.array([[r_BN_N[0]], [r_BN_N[1]], [r_BN_N[2]]])
        r_BT_N_unit = r_BT_N / np.linalg.norm(r_BT_N)
        
        # Define desired camera pointing in spacecraft body frame
        # Assuming camera points along +Y body axis
        r_Camera_B = np.array([[0], [1], [0]])
        
        # Current rotation from inertial to body
        dcm_BN = DynModel.simpleNavObject.attOutMsg.read().dcm_BN
        
        # Convert target vector to body frame
        r_BT_B = dcm_BN @ r_BT_N_unit
        
        # Compute rotation from current to desired attitude
        v1 = r_Camera_B / np.linalg.norm(r_Camera_B)
        v2 = r_BT_B / np.linalg.norm(r_BT_B)
        
        # Cross product and dot product
        cross = np.cross(v1.T, v2.T).T
        dot = np.dot(v1.T, v2.T)
        
        # Handle parallel vectors case
        if np.isclose(np.linalg.norm(cross), 0):
            if dot > 0:  # Same direction
                return [[0], [0], [0]]  # No rotation needed
            else:  # Opposite direction - need 180° rotation, choose arbitrary axis
                return [[1], [0], [0]]  # 180° around x-axis
        
        # Convert axis-angle to MRP (Modified Rodrigues Parameters)
        axis = cross / np.linalg.norm(cross)
        angle = np.arccos(np.clip(dot, -1.0, 1.0))
        
        # MRP = tan(angle/4) * axis
        tan_angle_4 = np.tan(angle/4)
        sigma_BR = tan_angle_4 * axis
        
        # Print targeting info every 30 seconds (in sim time)
        static_time = int(time_hours * 120)  # Create a time bucket for every 30 sec
        if not hasattr(self, '_last_target_print') or self._last_target_print != static_time:
            self._last_target_print = static_time
            print(f"Targeting: {lat:.2f}°N, {lon:.2f}°E | Att Error: {np.linalg.norm(sigma_BR):.6f}")
        
        return sigma_BR

    def process_tasks(self, time):
        """Override to implement target switching and pointing"""
        # Switch targets every 5 minutes
        time_min = time * macros.NANO2MIN
        if time_min - self.lastTargetSwitchTime >= 5.0:
            self.currentTargetIndex = (self.currentTargetIndex + 1) % len(self.targets)
            self.lastTargetSwitchTime = time_min
            print(f"Switching to target: {self.targets[self.currentTargetIndex]['name']} at {time_min:.1f} minutes")
        
        # Get current target
        target = self.targets[self.currentTargetIndex]
        
        # Calculate pointing to current target
        sigma_BR = self.target_ground_location(time, target["lat"], target["lon"])
        
        # Set the attitude guidance command
        attGuidMsg = messaging.AttGuidMsg()
        attGuidMsg.sigma_BR = sigma_BR
        attGuidMsg.omega_BR_B = [[0.0], [0.0], [0.0]]  # No rate command
        attGuidMsg.domega_BR_B = [[0.0], [0.0], [0.0]]  # No acceleration command
        
        # Write to the FSW attitude guidance message
        self.get_FswModel().attGuidMsg.write(attGuidMsg)

        # Update camera field of view based on distance to target
        self.update_camera_fov(time)
        
        # Call the base class process_tasks
        super(scenario_SatCameraTarget, self).process_tasks(time)
    
    def update_camera_fov(self, time):
        """Update camera field of view based on distance to target"""
        if hasattr(self, 'camera') and hasattr(self, 'viz') and self.viz is not None and self.camera is not None:
            try:
                DynModel = self.get_DynModel()
                r_BN_N = DynModel.simpleNavObject.transOutMsg.read().r_BN_N
                sc_alt = np.linalg.norm(r_BN_N) - 6371000.0  # Approximate altitude above Earth
                # Adjust FOV based on altitude - narrower at higher altitudes
                fov = min(20.0, max(3.0, 10.0 * (500000.0 / sc_alt))) * macros.D2R
                self.camera.fieldOfView = fov
            except Exception as e:
                # Just log the error but don't halt execution
                print(f"Camera FOV update error: {e}")
    
    def plot_target_visibility(self, timeData):
        """
        Plot target visibility information
        
        Args:
            timeData: Array of simulation times in minutes
            
        Returns:
            fig: The matplotlib figure object that can be shown or saved
        """
        import matplotlib.pyplot as plt
        
        # Get spacecraft position data
        posData = self.msgRecList[self.sNavTransName]
        r_BN_N = np.array(posData.r_BN_N)
        
        # Calculate visibility for each target over time
        visibilityData = []
        targetNames = []
        
        for target in self.targets:
            targetNames.append(target["name"])
            visibility = []
            
            for i, t in enumerate(timeData):
                # Convert time to hours for Earth rotation calculation
                time_hours = t / 60.0  # minutes to hours
                earth_rotation = time_hours * 15.0 * macros.D2R  # Earth rotates ~15 degrees per hour
                
                # Create rotation matrix
                c_rot = np.cos(earth_rotation)
                s_rot = np.sin(earth_rotation)
                R_N_E = np.array([
                    [c_rot, -s_rot, 0],
                    [s_rot, c_rot, 0],
                    [0, 0, 1]
                ])
                
                # Convert lat/lon to ECEF coordinates
                lat_rad = target["lat"] * macros.D2R
                lon_rad = target["lon"] * macros.D2R
                r_Earth = 6371000.0  # Earth radius in meters
                
                # Target position in Earth-fixed frame
                r_Target_E = np.array([
                    (r_Earth) * np.cos(lat_rad) * np.cos(lon_rad),
                    (r_Earth) * np.cos(lat_rad) * np.sin(lon_rad),
                    (r_Earth) * np.sin(lat_rad)
                ])
                
                # Convert to inertial frame
                r_Target_N = R_N_E @ r_Target_E
                
                # Check if target is visible (above horizon)
                sc_pos = r_BN_N[i]
                sc_to_target = r_Target_N - sc_pos
                sc_to_center = -sc_pos
                
                # If angle between these vectors is less than 90 degrees, target is visible
                cos_angle = np.dot(sc_to_target, sc_to_center) / (np.linalg.norm(sc_to_target) * np.linalg.norm(sc_to_center))
                
                # Set a value that represents visibility (1 for visible, 0 for not visible)
                if cos_angle < 0:  # Angle > 90 degrees, target is visible
                    visibility.append(1)
                else:
                    visibility.append(0)
            
            visibilityData.append(visibility)
        
        # Create the plot
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        
        for i, vis in enumerate(visibilityData):
            ax.plot(timeData, [i+1]*len(timeData), 'k-', alpha=0.2)  # baseline
            visible_segments = []
            start_idx = None
            
            # Find continuous visible segments
            for j, v in enumerate(vis):
                if v == 1 and start_idx is None:
                    start_idx = j
                elif v == 0 and start_idx is not None:
                    visible_segments.append((timeData[start_idx], timeData[j-1], i+1))
                    start_idx = None
            
            # Handle case where visibility extends to the end
            if start_idx is not None:
                visible_segments.append((timeData[start_idx], timeData[-1], i+1))
            
            # Plot visible segments
            for start, end, level in visible_segments:
                ax.plot([start, end], [level, level], 'g-', linewidth=3)
        
        ax.set_yticks(range(1, len(targetNames)+1))
        ax.set_yticklabels(targetNames)
        ax.set_xlabel('Time [min]')
        ax.set_title('Target Visibility Timeline')
        ax.grid(True)
        
        # Return the figure instead of trying to store it in BSK_plt.figureList
        return fig

    def pull_outputs(self, showPlots):
        """Process and plot simulation outputs"""
        # Import matplotlib here to ensure it's available
        import matplotlib.pyplot as plt
        
        # FSW process outputs, remove first data point as it is before FSW is called
        attErrRec = self.msgRecList[self.attGuidName]

        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
        
        num_RW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
        RW_friction = []
        for i in range(num_RW):
            RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))

        # Plot results
        BSK_plt.clear_all_plots()
        timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
        BSK_plt.plot_attitude_error(timeData, sigma_BR, id=1)
        BSK_plt.plot_rate_error(timeData, omega_BR_B, id=2)
        BSK_plt.plot_rw_speeds(timeData, RW_speeds, num_RW, id=3)
        BSK_plt.plot_rw_friction(timeData, RW_friction, num_RW, self.DynModels.RWFaultLog, id=4)
        
        # Create the target visibility figure (using an explicit figure number)
        vis_fig = self.plot_target_visibility(timeData)
        
        # Make sure we have a local figureList dict
        figureList = {}
        
        if showPlots:
            BSK_plt.show_all_plots()
            plt.figure(vis_fig.number)  # Ensure visibility figure is shown
            plt.show()
        else:
            fileName = os.path.basename(os.path.splitext(__file__)[0])
            figureNames = ["attitudeErrorNorm", "rateError", "RWSpeeds", "RWFriction"]
            figureList = BSK_plt.save_all_plots(fileName, figureNames)
            
            # Manually add the visibility figure
            pltName = fileName + "_targetVisibility"
            figureList[pltName] = vis_fig
            vis_fig.savefig(pltName + ".png")

        return figureList

    def setup_camera(self):
        """Configure a camera for visualization in Vizard"""
        if not hasattr(self, 'viz') or self.viz is None:
            print("Vizard not available - skipping camera setup")
            return None
            
        try:
            # Create a camera with better positioning and narrower FOV for spacecraft imaging
            camera = vizSupport.createStandardCamera(
                self.viz,
                setMode=1,  # Camera mode (1 = spacecraft-fixed)
                spacecraftName=self.DynModels.scObject.ModelTag,
                displayName="ImagingCamera",
                fieldOfView=5 * macros.D2R,  # 5 degrees FOV - narrow like a real imaging camera
                pointingVector_B=[0.0, 1.0, 0.0],  # Camera points along Y-axis
                position_B=[0.0, 1.5, 0.0]  # Position on spacecraft
            )
            
            # Additional camera configuration
            if hasattr(camera, 'viewScaleFactor'):
                camera.viewScaleFactor = 1.0  # Adjust based on your needs
                
            # Some Vizard versions support these functions
            if hasattr(camera, 'setFrustumColors'):
                camera.setFrustumColors([0.2, 0.8, 0.2, 0.3], [1.0, 1.0, 0.0, 0.3])
                
            # Display camera FOV if supported
            if hasattr(camera, 'showFOV'):
                camera.showFOV = True
                
            print("Camera configured successfully")
            return camera
            
        except Exception as e:
            print(f"Error configuring camera: {e}")
            return None
    
    def setup_rw_fail_safe(self):
        """
        Set up reaction wheel fault-tolerance and handling
        This may help with the memory leak issues
        """
        # Add proper cleanup code for RW config messages
        try:
            # Clean up the RW configs by using proper Python SWIG memory management
            import gc
            
            DynModel = self.get_DynModel()
            if hasattr(DynModel, 'rwStateEffector'):
                # Force garbage collection to help with memory leaks
                for i in range(4):
                    if hasattr(DynModel.rwStateEffector, f'rwParamMsg{i}'):
                        setattr(DynModel.rwStateEffector, f'rwParamMsg{i}', None)
                
                gc.collect()
        except Exception as e:
            print(f"Error setting up RW fail safe: {e}")
    
def runScenario(scenario):
    """method to initialize and execute the scenario"""

    simulationTime = macros.min2nano(30.)
    
    # Set the mode request attribute required by the BSK event system
    scenario.modeRequest = "hillPoint"
    
    # Initialize the simulation
    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    
    # Setup RW fault tolerance to address memory leaks
    scenario.setup_rw_fail_safe()
    
    scenario.ExecuteSimulation()

    return


def run(showPlots=True, use_vizard=True, targets=None):
    """
        The scenarios can be run with the following parameters:

        Args:
            showPlots (bool): Determines if the script should display plots
            use_vizard (bool): Whether to attempt to connect to Vizard for visualization
            targets (list): Optional list of target dictionaries with "name", "lat", "lon" keys
    """
    scenario = scenario_SatCameraTarget(use_vizard=use_vizard)
    
    # Override default targets if provided
    if targets:
        scenario.targets = targets
    
    runScenario(scenario)
    figureList = scenario.pull_outputs(showPlots)
    return figureList

if __name__ == "__main__":
    # Define your custom targets here (optional)
    my_targets = [
        {"name": "Washington DC", "lat": 38.9072, "lon": -77.0369},
        {"name": "Beijing", "lat": 39.9042, "lon": 116.4074},
        {"name": "Sydney", "lat": -33.8688, "lon": 151.2093}
    ]
    
    # Run with visualization and custom targets
    run(True, use_vizard=True, targets=my_targets)