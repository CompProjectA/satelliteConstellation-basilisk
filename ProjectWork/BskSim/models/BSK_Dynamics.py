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

from Basilisk.utilities import macros


import numpy as np
from Basilisk import __path__
from Basilisk.simulation import ephemerisConverter
from Basilisk.simulation import (spacecraft, extForceTorque, simpleNav,
                                 reactionWheelStateEffector, coarseSunSensor, eclipse)
from Basilisk.simulation import thrusterDynamicEffector
from Basilisk.utilities import RigidBodyKinematics as rbk
from Basilisk.utilities import macros as mc
from Basilisk.utilities import simIncludeRW, simIncludeGravBody
from Basilisk.utilities import simIncludeThruster
from Basilisk.utilities import unitTestSupport as sp
from Basilisk.utilities import macros 

bskPath = __path__[0]


class BSKDynamicModels():
    """
    General bskSim simulation class that sets up the spacecraft simulation configuration.

    """
    def __init__(self, SimBase, dynRate):
        # define empty class variables
        self.sun = None
        self.earth = None
        self.moon = None
        self.epochMsg = None
        self.RW1 = None
        self.RW2 = None
        self.RW3 = None
        self.RW4 = None

        # Define process name, task name and task time-step
        self.processName = SimBase.DynamicsProcessName
        self.taskName = "DynamicsTask"
        self.processTasksTimeStep = mc.sec2nano(dynRate)

        # Create task
        SimBase.dynProc.addTask(SimBase.CreateNewTask(self.taskName, self.processTasksTimeStep))

        # Instantiate Dyn modules as objects
        self.scObject = spacecraft.Spacecraft()
        self.gravFactory = simIncludeGravBody.gravBodyFactory()
        self.rwFactory = simIncludeRW.rwFactory()
        self.extForceTorqueObject = extForceTorque.ExtForceTorque()
        self.simpleNavObject = simpleNav.SimpleNav()
        self.eclipseObject = eclipse.Eclipse()
        self.CSSConstellationObject = coarseSunSensor.CSSConstellation()
        self.rwStateEffector = reactionWheelStateEffector.ReactionWheelStateEffector()
        self.thrustersDynamicEffector = thrusterDynamicEffector.ThrusterDynamicEffector()
        self.EarthEphemObject = ephemerisConverter.EphemerisConverter()

        # Initialize all modules and write init one-time messages
        self.InitAllDynObjects()

        # Assign initialized modules to tasks
        SimBase.AddModelToTask(self.taskName, self.scObject, 201)
        SimBase.AddModelToTask(self.taskName, self.simpleNavObject, 109)
        SimBase.AddModelToTask(self.taskName, self.gravFactory.spiceObject, 200)
        SimBase.AddModelToTask(self.taskName, self.EarthEphemObject, 199)
        SimBase.AddModelToTask(self.taskName, self.CSSConstellationObject, 108)
        SimBase.AddModelToTask(self.taskName, self.eclipseObject, 204)
        SimBase.AddModelToTask(self.taskName, self.rwStateEffector, 301)
        SimBase.AddModelToTask(self.taskName, self.extForceTorqueObject, 300)
        

        SimBase.createNewEvent("addOneTimeRWFault", self.processTasksTimeStep, True,
            ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
            ["self.DynModels.AddRWFault('friction',0.05,2, self.TotalSim.CurrentNanos)", "self.oneTimeRWFaultFlag=0"])

        
        SimBase.createNewEvent("addRepeatedRWFault", self.processTasksTimeStep, True,
            ["self.repeatRWFaultFlag==1"],
            ["self.DynModels.PeriodicRWFault(1./3000,'friction',0.005,3, self.TotalSim.CurrentNanos)", "self.setEventActivity('addRepeatedRWFault',True)"])

        self.RWFaultLog = []
    # ------------------------------------------------------------------------------------------- #
    # These are module-initialization methods

    def SetSpacecraftHub(self):
        """
        Specify the spacecraft hub parameters.
        """
        self.scObject.ModelTag = "bskSat"
        # -- Crate a new variable for the sim sc inertia I_sc. Note: this is currently accessed from FSWClass
        self.I_sc = [900., 0., 0.,
                     0., 800., 0.,
                     0., 0., 600.]
        self.scObject.hub.mHub = 750.0  # kg - spacecraft mass
        self.scObject.hub.r_BcB_B = [[0.0], [0.0], [0.0]]  # m - position vector of body-fixed point B relative to CM
        self.scObject.hub.IHubPntBc_B = sp.np2EigenMatrix3d(self.I_sc)

    def SetGravityBodies(self):
        """
        Specify what gravitational bodies to include in the simulation
        """
        timeInitString = "2012 MAY 1 00:28:30.0"
        gravBodies = self.gravFactory.createBodies(['sun', 'earth', 'moon'])
        gravBodies['earth'].isCentralBody = True
        self.sun = 0
        self.earth = 1
        self.moon = 2

        self.gravFactory.addBodiesTo(self.scObject)
        self.gravFactory.createSpiceInterface(bskPath + '/supportData/EphemerisData/',
                                              timeInitString,
                                              epochInMsg=True)
        self.epochMsg = self.gravFactory.epochMsg

        self.gravFactory.spiceObject.zeroBase = 'Earth'

        self.EarthEphemObject.addSpiceInputMsg(self.gravFactory.spiceObject.planetStateOutMsgs[self.earth])

    def SetEclipseObject(self):
        """
        Specify what celestial object is causing an eclipse message.
        """
        self.eclipseObject.ModelTag = "eclipseObject"
        self.eclipseObject.sunInMsg.subscribeTo(self.gravFactory.spiceObject.planetStateOutMsgs[self.sun])
        # add all celestial objects in spiceObjects except for the sun (0th object)
        for c in range(1, len(self.gravFactory.spiceObject.planetStateOutMsgs)):
            self.eclipseObject.addPlanetToModel(self.gravFactory.spiceObject.planetStateOutMsgs[c])
        self.eclipseObject.addSpacecraftToModel(self.scObject.scStateOutMsg)

    def SetExternalForceTorqueObject(self):
        """Set the external force and torque object."""
        self.extForceTorqueObject.ModelTag = "externalDisturbance"
        self.scObject.addDynamicEffector(self.extForceTorqueObject)

    def SetSimpleNavObject(self):
        """Set the navigation sensor object."""
        self.simpleNavObject.ModelTag = "SimpleNavigation"
        self.simpleNavObject.scStateInMsg.subscribeTo(self.scObject.scStateOutMsg)


    def SetReactionWheelDynEffector(self):
        """Set the 4 reaction wheel devices."""
        # specify RW momentum capacity
        maxRWMomentum = 50.  # Nms
        
        # Store default Coulomb friction value
        self.default_coulomb_friction = 0.0005  # N⋅m

        # Define orthogonal RW pyramid
        # -- Pointing directions
        rwElAngle = np.array([40.0, 40.0, 40.0, 40.0])*mc.D2R
        rwAzimuthAngle = np.array([45.0, 135.0, 225.0, 315.0])*mc.D2R
        rwPosVector = [[0.8, 0.8, 1.79070],
                       [0.8, -0.8, 1.79070],
                       [-0.8, -0.8, 1.79070],
                       [-0.8, 0.8, 1.79070]
                       ]

        gsHat = (rbk.Mi(-rwAzimuthAngle[0], 3).dot(rbk.Mi(rwElAngle[0], 2))).dot(np.array([1, 0, 0]))
        self.RW1 = self.rwFactory.create('Honeywell_HR16',
                                         gsHat,
                                         maxMomentum=maxRWMomentum,
                                         rWB_B=rwPosVector[0],
                                         fCoulomb=self.default_coulomb_friction)  # Set default friction
        
        gsHat = (rbk.Mi(-rwAzimuthAngle[1], 3).dot(rbk.Mi(rwElAngle[1], 2))).dot(np.array([1, 0, 0]))
        self.RW2 = self.rwFactory.create('Honeywell_HR16',
                                         gsHat,
                                         maxMomentum=maxRWMomentum,
                                         rWB_B=rwPosVector[1],
                                         fCoulomb=self.default_coulomb_friction)  # Set default friction

        gsHat = (rbk.Mi(-rwAzimuthAngle[2], 3).dot(rbk.Mi(rwElAngle[2], 2))).dot(np.array([1, 0, 0]))
        self.RW3 = self.rwFactory.create('Honeywell_HR16',
                                         gsHat,
                                         maxMomentum=maxRWMomentum,
                                         rWB_B=rwPosVector[2],
                                         fCoulomb=self.default_coulomb_friction)  # Set default friction
            
        gsHat = (rbk.Mi(-rwAzimuthAngle[3], 3).dot(rbk.Mi(rwElAngle[3], 2))).dot(np.array([1, 0, 0]))
        self.RW4 = self.rwFactory.create('Honeywell_HR16',
                                         gsHat,
                                         maxMomentum=maxRWMomentum,
                                         rWB_B=rwPosVector[3],
                                         fCoulomb=self.default_coulomb_friction)  # Set default friction

        self.rwFactory.addToSpacecraft("RWA", self.rwStateEffector, self.scObject)
        
        print(f"Reaction wheels configured with default Coulomb friction: {self.default_coulomb_friction} N⋅m")


# 3. Add a method to get current friction values (useful for diagnostics):

    def get_rw_friction_values(self):
        """Get current friction values for all reaction wheels"""
        friction_values = {}
        if hasattr(self, 'RW1') and self.RW1:
            friction_values['RW1'] = self.RW1.fCoulomb
        if hasattr(self, 'RW2') and self.RW2:
            friction_values['RW2'] = self.RW2.fCoulomb
        if hasattr(self, 'RW3') and self.RW3:
            friction_values['RW3'] = self.RW3.fCoulomb
        if hasattr(self, 'RW4') and self.RW4:
            friction_values['RW4'] = self.RW4.fCoulomb
        return friction_values


    def SetThrusterStateEffector(self):
        """Set the 8 ACS thrusters."""
        # Make a fresh TH factory instance, this is critical to run multiple times
        thFactory = simIncludeThruster.thrusterFactory()

        # 8 thrusters are modeled that act in pairs to provide the desired torque
        thPos = [
            [825.5/1000.0, 880.3/1000.0, 1765.3/1000.0],
            [825.5/1000.0, 880.3/1000.0, 260.4/1000.0],
            [880.3/1000.0, 825.5/1000.0, 1765.3/1000.0],
            [880.3/1000.0, 825.5/1000.0, 260.4/1000.0],
            [-825.5/1000.0, -880.3/1000.0, 1765.3/1000.0],
            [-825.5/1000.0, -880.3/1000.0, 260.4/1000.0],
            [-880.3/1000.0, -825.5/1000.0, 1765.3/1000.0],
            [-880.3/1000.0, -825.5/1000.0, 260.4/1000.0]
                 ]
        thDir = [
            [0.0, -1.0, 0.0],
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0]
        ]
        for pos_B, dir_B in zip(thPos, thDir):
            thFactory.create(
                'MOOG_Monarc_1'
                , pos_B
                , dir_B
            )
        # create thruster object container and tie to spacecraft object
        thFactory.addToSpacecraft("ACS Thrusters",
                                  self.thrustersDynamicEffector,
                                  self.scObject)

    def SetCSSConstellation(self):
        """Set the 8 CSS sensors"""
        self.CSSConstellationObject.ModelTag = "cssConstellation"

        def setupCSS(cssDevice):
            cssDevice.fov = 80. * mc.D2R         # half-angle field of view value
            cssDevice.scaleFactor = 2.0
            cssDevice.sunInMsg.subscribeTo(self.gravFactory.spiceObject.planetStateOutMsgs[self.sun])
            cssDevice.stateInMsg.subscribeTo(self.scObject.scStateOutMsg)
            cssDevice.sunEclipseInMsg.subscribeTo(self.eclipseObject.eclipseOutMsgs[0])
            cssDevice.this.disown()

        # setup CSS sensor normal vectors in body frame components
        nHat_B_List = [
            [0.0, 0.707107, 0.707107],
            [0.707107, 0., 0.707107],
            [0.0, -0.707107, 0.707107],
            [-0.707107, 0., 0.707107],
            [0.0, -0.965926, -0.258819],
            [-0.707107, -0.353553, -0.612372],
            [0., 0.258819, -0.965926],
            [0.707107, -0.353553, -0.612372]
        ]
        numCSS = len(nHat_B_List)

        # store all
        cssList = []
        for nHat_B, i in zip(nHat_B_List, list(range(1,numCSS+1))):
            CSS = coarseSunSensor.CoarseSunSensor()
            setupCSS(CSS)
            CSS.ModelTag = "CSS" + str(i)
            CSS.nHat_B = np.array(nHat_B)
            cssList.append(CSS)

        # assign the list of CSS devices to the CSS array class
        self.CSSConstellationObject.sensorList = coarseSunSensor.CSSVector(cssList)

    # Method for adding reaction wheel faults

    def AddRWFault(self, faultType, fault, faultRW, currentTime):
        """
        Adds a friction fault to the reaction wheel with proper GUI parameter handling.
        
        Parameters:
        faultType (str): Type of fault ('friction', 'power_limit', etc.)
        fault (float): Fault magnitude value
        faultRW (int): Reaction wheel index (0-3, displayed as RW 1-4)
        currentTime (int): Current simulation time in nanoseconds
        
        Returns:
        bool: True if fault was successfully applied
        """
        # Ensure RWFaultLog exists
        if not hasattr(self, 'RWFaultLog'):
            self.RWFaultLog = []
            
        # Log the fault with proper time conversion
        fault_time_min = currentTime * macros.NANO2MIN
        self.RWFaultLog.append([faultType, fault, faultRW, fault_time_min])
        
        # Validate RW index
        if faultRW < 0 or faultRW > 3:
            print(f"ERROR: Invalid RW index {faultRW}. Must be 0-3.")
            return False
        
        # Get the reaction wheel list
        rw_list = [self.RW1, self.RW2, self.RW3, self.RW4]
        
        if faultType == "friction":
            # Apply friction fault to the specified wheel
            if rw_list[faultRW] is not None:
                # Get current friction value
                current_friction = rw_list[faultRW].fCoulomb
                
                # Apply the fault (add to existing friction)
                rw_list[faultRW].fCoulomb = current_friction + fault
                
                # Log the fault application with user-friendly numbering
                print(f"FRICTION FAULT APPLIED:")
                print(f"  - Reaction Wheel: RW{faultRW + 1}")
                print(f"  - Previous friction: {current_friction:.6f} N⋅m")
                print(f"  - Fault magnitude: {fault:.6f} N⋅m")
                print(f"  - New total friction: {rw_list[faultRW].fCoulomb:.6f} N⋅m")
                print(f"  - Time: {fault_time_min:.2f} minutes")
                
                return True
            else:
                print(f"ERROR: RW{faultRW + 1} not available")
                return False
                
        elif faultType == "power_limit":
            # Power limit fault implementation
            print(f"POWER LIMIT FAULT APPLIED:")
            print(f"  - Reaction Wheel: RW{faultRW + 1}")
            print(f"  - Power limit: {fault} W")
            print(f"  - Time: {fault_time_min:.2f} minutes")
            
            # Store power limit for use in control logic
            if not hasattr(self, 'power_limits'):
                self.power_limits = {}
            self.power_limits[faultRW] = fault
            
            return True
            
        elif faultType == "encoder":
            # Encoder fault implementation
            print(f"ENCODER FAULT APPLIED:")
            print(f"  - Reaction Wheel: RW{faultRW + 1}")
            print(f"  - Error magnitude: {fault}%")
            print(f"  - Time: {fault_time_min:.2f} minutes")
            
            # Store encoder error for use in measurement
            if not hasattr(self, 'encoder_errors'):
                self.encoder_errors = {}
            self.encoder_errors[faultRW] = fault / 100.0  # Convert percentage to fraction
            
            return True
            
        elif faultType == "battery":
            # Battery fault is handled separately in the main simulation
            print(f"BATTERY FAULT NOTED (handled by battery simulation):")
            print(f"  - Additional drain: {fault} W")
            print(f"  - Time: {fault_time_min:.2f} minutes")
            return True
            
        else:
            print(f"ERROR: Unknown fault type '{faultType}'")
            return False


    def get_rw_friction_values(self):
        """
        Get current friction values for all reaction wheels with RW 1-4 labels
        
        Returns:
        dict: Dictionary with RW1-4 labels and current friction values
        """
        friction_values = {}
        
        if hasattr(self, 'RW1') and self.RW1:
            friction_values['RW1'] = self.RW1.fCoulomb
        if hasattr(self, 'RW2') and self.RW2:
            friction_values['RW2'] = self.RW2.fCoulomb
        if hasattr(self, 'RW3') and self.RW3:
            friction_values['RW3'] = self.RW3.fCoulomb
        if hasattr(self, 'RW4') and self.RW4:
            friction_values['RW4'] = self.RW4.fCoulomb
            
        return friction_values


    def PeriodicRWFault(self, probability, faultType, fault, faultRW, currentTime):
        """
        Adds a fault periodically based on probability.
        
        Parameters:
        probability (float): Chance of fault occurring per update (0.0 to 1.0)
        faultType (str): Type of fault
        fault (float): Fault magnitude
        faultRW (int): Reaction wheel index (0-3)
        currentTime (int): Current time in nanoseconds
        """
        if np.random.uniform() < probability:
            print(f"PERIODIC FAULT TRIGGERED (probability {probability:.2%})")
            self.AddRWFault(faultType, fault, faultRW, currentTime)
            return True
        return False
        
        
    

    
    # Global call to initialize every module
    def InitAllDynObjects(self):
        """
        Initialize all the dynamics objects.
        """
        self.SetSpacecraftHub()
        self.SetGravityBodies()
        self.SetExternalForceTorqueObject()
        self.SetSimpleNavObject()
        self.SetEclipseObject()
        self.SetCSSConstellation()

        self.SetReactionWheelDynEffector()
        self.SetThrusterStateEffector()

