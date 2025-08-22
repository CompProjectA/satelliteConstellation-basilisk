import os
from message_data import MessageData
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.simulation import spacecraft, spacecraftLocation, simSynch,extForceTorque,simpleNav
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion, simIncludeGravBody, unitTestSupport, vizSupport)
from Basilisk import __path__
from Basilisk.architecture import messaging
from Basilisk.simulation import simplePowerSink, simpleSolarPanel
from Basilisk.utilities import SimulationBaseClass
from satellite_definition import LeadingSatellite, ChildSatellite

bskPath = __path__[0]

# Define a simple message payload for strings
from Basilisk.architecture import messaging
from Basilisk.architecture import sysModel

def thereIsOne(value, lst):
    for item in lst:
        if item == value:
            return True
    return False
#return the first index that matches value
def getFirstMatch(value, lst):
    for i,item in enumerate(lst):
        if item == value:
            return i
    return None
def checkNext10Seconds(index, list):
    print("Checking next 10 seconds for leader...")
    for i in range (10):
        if list[index + i + 1]!=1:
            return False
    return True
def run(show_plots=True):
    print("Basilisk Path:", bskPath)

    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    simulationTime = 1200
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    num_satellites = 16
    leaders = []
    satellites = []

    cameraLocation = [0.0, 2.0, 0.0]
    #target area
    targets = [ {"name": "Tokyo",          "lat":  35.68,   "lon": 139.77,   "color": "green"},
                {"name": "Sri Lanka",      "lat":   6.9271, "lon":  79.8612, "color": "orange"},
                {"name": "Central Africa", "lat":   1.5333, "lon":  17.6667, "color": "purple"}
    ]
    #Define variables and parameters for battery fault
    numDataPoints = 100
    solarPanel = None
    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    earth.isCentralBody = True
    mu = earth.mu
    earth_radius = earth.radEquator

     # Ensure the sun is not the central body for the orbits of the satellites
  
    rLEO_base = 7000. * 1000  # 7000 km
   
    # Define the angular speed based on the base velocity and position
    # angular_speed=vN_base/ np.linalg.norm(rN_base)
    # print(rN_base)
    # Formation offsets relative to the leader - TIGHT FORMATION
    for i in range (num_satellites):
        if(i>=0 and i<4):
            oe = orbitalMotion.ClassicElements()
            oe.a = 7000e3   
            oe.e = 0.0  # Circular
            oe.i = (55+i*1) * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = 45 * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            if i == 0:
                sat = LeadingSatellite(i, rN, vN,gravFactory, f"Satellite{i+1}")
                leaders.append(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm(earth_radius, [0.2,-0.4, 0.2], 1000e3)
            else:
                sat = ChildSatellite(i, rN, vN,gravFactory, leaders[0])
                leaders[0].add_child(sat)
                satellites.append(sat.sc)
               
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm()
        if(i>=4 and i<8):
            oe = orbitalMotion.ClassicElements()
            oe.a = 7000e3   
            oe.e = 0.0  # Circular
            oe.i = (55+(i-4)*1) * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = 225 * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            if i == 4:
                sat = LeadingSatellite(i, rN, -vN,gravFactory, f"Satellite{i+1}")
                leaders.append(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm(earth_radius, [-0.2,0.4, -0.2], 1000e3)
            else:
                sat = ChildSatellite(i, rN, -vN,gravFactory, leaders[1])
                leaders[1].add_child(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm()
        if(i>=8 and i<12):
            oe = orbitalMotion.ClassicElements()
            oe.a = 7000e3   
            oe.e = 0.0  # Circular
            oe.i = 56 * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = (135+(i-8)*2) * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            if i == 8:
                sat = LeadingSatellite(i, rN, vN,gravFactory, f"Satellite{i+1}")
                leaders.append(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm(earth_radius, [-0.2,-0.4, -0.4], 1000e3)
            else:
                sat = ChildSatellite(i, rN, vN,gravFactory, leaders[2])
                leaders[2].add_child(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm()
        if(i>=12 and i<16):
            oe = orbitalMotion.ClassicElements()
            oe.a = 7000e3   
            oe.e = 0.0  # Circular
            oe.i = 56 * macros.D2R
            oe.Omega = 0
            oe.omega = 0    
            oe.f = (315+(i-8)*2) * macros.D2R  # Same angular position
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            print(rN)
            if i == 12:
                sat = LeadingSatellite(i, rN, -vN,gravFactory, f"Satellite{i+1}")
                leaders.append(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm(earth_radius, [0.2,0.1, 0.1], 1000e3)
            else:
                sat = ChildSatellite(i, rN, -vN,gravFactory, leaders[3])
                leaders[3].add_child(sat)
                satellites.append(sat.sc)
                scSim.AddModelToTask(simTaskName, sat.sc)
                sat.setup_comm()
    
    # Set up communication modules
    scLocationModules = []
    scLocationLeadersModules = []
    accessRecorders = []
    accessLeadersRecorders = []
    maximum_range = 2000e3  # Increased maximum range to 5000 km
    print(f"Maximum Range: {maximum_range/1000} km")
   # Communicatio between leaders
    for i, leader in enumerate(leaders):
        scLocation = spacecraftLocation.SpacecraftLocation()
        scLocation.ModelTag = f"CommCheck between leaders from {i+1}"
        scLocation.rEquator = earth_radius
        scLocation.rPolar = earth_radius * 0.98
        scLocation.aHat_B = [-0.2,-0.4, -0.4]  # Pointing in the +Y direction
        scLocation.theta = np.radians(30.0)  # Adjust as needed to be narrower
        scLocation.maximumRange = maximum_range
        scLocation.primaryScStateInMsg.subscribeTo(leader.sc.scStateOutMsg)
        for j, leaderFollowed in enumerate(leaders):
            if j != i:
               scLocation.addSpacecraftToModel(leaderFollowed.sc.scStateOutMsg)
        scLocationLeadersModules.append(scLocation)
        scSim.AddModelToTask(simTaskName, scLocation)
        for j in range(3):
            recorder_name = f"AccessRecorder{j+1}"
            recorder = scLocation.accessOutMsgs[j].recorder()
            recorder.ModelTag = recorder_name
            accessLeadersRecorders.append(recorder)
            scSim.AddModelToTask(simTaskName, recorder)

    #Recorders and modules initialization for parent-child
    for i, leader in enumerate(leaders):
        scLocationModules.append(leader.comm_module) 
        scSim.AddModelToTask(simTaskName, leader.comm_module)
        for j in range(3):
                recorder_name = f"AccessRecorder{j+1}"
                recorder = leader.comm_module.accessOutMsgs[j].recorder()
                recorder.ModelTag = recorder_name
                accessRecorders.append(recorder)
                scSim.AddModelToTask(simTaskName, recorder)

   

    # Battery state logger
 
    # Visualization
    if vizSupport.vizFound:
        print("Vizard Visualization Found. Enabling live streaming...")
        clockSync = simSynch.ClockSynch()
        clockSync.accelFactor = 50.0
        scSim.AddModelToTask(simTaskName, clockSync)
        
        viz = vizSupport.enableUnityVisualization(scSim,
                                                 simTaskName,
                                                 satellites,
            
                                                 liveStream=True,

                                                 )

     
        
        for i, scLocation in enumerate(leaders):
            vizSupport.addLocation(
                viz,
                stationName=f"Satellite{i*4 + 1}",
                parentBodyName=f"Satellite{i*4 + 1}",
                r_GP_P=[0, 0, 0],
                gHat_P=scLocation.aHat_B,
                fieldOfView=np.pi / 4,
                range=scLocation.comm_module.maximumRange
            )
        for i, scLocation in enumerate(scLocationLeadersModules):
            vizSupport.addLocation(
                viz,
                stationName=f"Satellite{i*4 + 1}",
                parentBodyName=f"Satellite{i*4 + 1}",
                r_GP_P=[0, 0, 0],
                gHat_P=[-0.2,-0.4, -0.4],
                fieldOfView=np.pi / 4,
                color="red",
                range=scLocation.maximumRange
            )

    else:
        print("Vizard Visualization Module Not Found. Check Basilisk installation.")
   
     # Simulate for 1200 seconds (20 minutes)
    stepSize = 120  # seconds
    checkWindow = 10  # seconds
    messageLeadSent = False
    currentTime = 0
    scSim.InitializeSimulation()
    #checking communication access for 10 seconds before sending the actual message 
    while currentTime < simulationTime:
        nextStopTime = macros.sec2nano(currentTime + stepSize)
        scSim.ConfigureStopTime(nextStopTime)
        scSim.ExecuteSimulation()

        # Get the time array and find indices for the last 10 seconds
        times = accessRecorders[0].times() * macros.NANO2SEC
        hasAccesstoSat2 = accessRecorders[0].hasAccess
        hasAccesstoSat3 = accessRecorders[1].hasAccess
        hasAccesstoSat4 = accessRecorders[2].hasAccess
        hasAccesstoLeader3to2 = accessLeadersRecorders[7].hasAccess
        # Find indices where time is in (currentTime + 110, currentTime + 120]
        timeStart = currentTime + stepSize - checkWindow  # e.g., 110
        timeEnd = currentTime + stepSize  # e.g., 120
        indices = (times > timeStart) & (times <= timeEnd)
        leadCheckindices = (times>0)&(times<=currentTime+stepSize)
        #Check and send hello messsage between parent and child, with condition of all children accessible before sending to one
        print("Checking 10 secs access to the child satellite...")
        if ((hasAccesstoSat2[indices] == 1).all() and
           (hasAccesstoSat3[indices] == 1).all() and
           (hasAccesstoSat4[indices] == 1).all()):
            print("check successful")
            leader = leaders[0]
            message = MessageData(f"Hello from{leader.model_tag}", (currentTime+stepSize)/60, leader.children[0])
            #Send a message from the leader1 to first child
            leader.sendMessage(f"Hello from{leader.model_tag}", (currentTime+stepSize)/60, leader.children[0] )
            if(message in leader.messageOutHistory):
                print("Sending message to the first child satellite successful!")
            else:
                print("Warning: Can not send the message")
                print(message in leader.messageOutHistory)
        #send message only once
        print (thereIsOne(1,hasAccesstoLeader3to2[leadCheckindices]))
        if((messageLeadSent==False)and(thereIsOne(1,hasAccesstoLeader3to2[leadCheckindices]))):
            #in 10 next seconds
            if((getFirstMatch(1,hasAccesstoLeader3to2[leadCheckindices])<(currentTime+stepSize-10))and
               checkNext10Seconds(getFirstMatch(1,hasAccesstoLeader3to2[leadCheckindices]), hasAccesstoLeader3to2[leadCheckindices])):
                print("check successful")
                leadSender = leaders[2]
                leadReceiver = leaders[1]
                timeSent = (getFirstMatch(1,hasAccesstoLeader3to2[leadCheckindices])+10)/60
                message = MessageData(f"Hello from{leadSender.model_tag}", timeSent, leadReceiver)
                #Send a message from the leader1 to second child
                leadSender.sendMessageToLead(f"Hello from{leadSender.model_tag}", timeSent, leadReceiver )
                if(message in leadSender.messageOutHistory):
                    print("Sending message from lead 3( sat 9) to lead 2 ( sat 2) successful!")
                    messageLeadSent = True
                else:
                    print("Warning: Can not send the message")
                    print(message in leadSender.messageOutHistory)
        currentTime += stepSize
    #show history in sat 1 to sat 2
    for i, message in enumerate(leaders[0].messageOutHistory):
        print(f"Message {i+1} to {message.objectActive.model_tag} from Satellite 1 at time {message.timeSent:.2f} minutes: {message.message_content}")
    # show message history in sat 5 (2nd leaders) ( message sent only once)
    for i, message in enumerate(leaders[2].messageOutHistory):
        print(f"Message {i+1} to {message.objectActive.model_tag} from Satellite 9 at time {message.timeSent:.2f} minutes: {message.message_content}")
    # Plotting access
    if show_plots:
        time_vec = accessLeadersRecorders[0].times() * macros.NANO2MIN
        plt.figure(figsize=(10, 6))

        
        for i in range(1,4 ):
                    # Plot access from leader to others and from others to leader
                    access_data = accessRecorders[i-1].hasAccess
                    label = f'Sat1 to Satellite{i+1}'
                    plt.plot(time_vec, access_data, label=label)
        for message in leaders[0].messageOutHistory:
            plt.annotate(f'Msg → Sat2', xy=(message.timeSent, 1), xytext=(message.timeSent, 1.01),
                     arrowprops=dict(arrowstyle='->', color='black'), fontsize=8, rotation=45)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader (Sat1) and Children (Sat2-4)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(10, 6))

        
        for i in range(4,7):
                    # Plot access from leader to others and from others to leader
                    access_data = accessRecorders[i-1].hasAccess
                    label = f'Sat5 to Satellite{i+1}th'
                    plt.plot(time_vec, access_data, label=label)
        # for message in leaders[0].messageOutHistory:
        #     plt.annotate(f'Msg → Sat2', xy=(message.timeSent, 1), xytext=(message.timeSent, 1.01),
        #              arrowprops=dict(arrowstyle='->', color='black'), fontsize=8, rotation=45)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader (Sat5) and Children (Sat6-8)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.figure(figsize=(10, 6))

        
        for i in range(7,10 ):
                    # Plot access from leader to others and from others to leader
                    access_data = accessRecorders[i-1].hasAccess
                    label = f'Sat9 to SatLeader{i+1}th'
                    plt.plot(time_vec, access_data, label=label)
        # for message in leaders[0].messageOutHistory:
        #     plt.annotate(f'Msg → Sat2', xy=(message.timeSent, 1), xytext=(message.timeSent, 1.01),
        #              arrowprops=dict(arrowstyle='->', color='black'), fontsize=8, rotation=45)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader (Sat9) and Children (Sat10-12)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.figure(figsize=(10, 6))

        
        for i in range(10,13):
                    # Plot access from leader to others and from others to leader
                    access_data = accessRecorders[i-1].hasAccess
                    label = f'Sat13 to SatLeader{i+1}th'
                    plt.plot(time_vec, access_data, label=label)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader (Sat13) and Children (Sat14-16)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.figure(figsize=(10, 6))
        for i in range(7,10):
                    # Plot access from leader to others and from others to leader
                    access_data = accessLeadersRecorders[i-1].hasAccess
                    label = f'Sat9 to SatLeader{i-6}th'
                    plt.plot(time_vec, access_data, label=label)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Has Access (1 = Yes, 0 = No)')
        plt.title('Access Between Leader3 to other leads')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
       
       

if __name__ == "__main__":
    run(show_plots=True)

