
import numpy as np
from Basilisk.simulation import spacecraft, spacecraftLocation
from Basilisk.simulation import simpleBattery 
from Basilisk.architecture import messaging
from message_data import MessageData

class ChildSatellite:
    def initializeBattery(self):
        self.battery = simpleBattery.SimpleBattery()
        self.battery.ModelTag = f"satBattery {self.index}"
        self.battery.storageCapacity   = 100.0 
        self.battery.storedCharge_Init = 50.0 
        self.batteryReader = messaging.PowerStorageStatusMsgReader()
        self.batteryReader.subscribeTo(self.battery.batPowerOutMsg)

    def __init__(self, index, rN, vN,gravFactory, leading_sat):
        self.index = index
        self.battery = None
        self.batteryReader = None
        self.initializeBattery()
        self.model_tag = f"Satellite{index+1}"
        self.messageInHistory: list[MessageData] = []
        self.messageOutHistory: list[MessageData] = []
        self.sc = spacecraft.Spacecraft()
        gravFactory.addBodiesTo(self.sc)
        self.sc.ModelTag = self.model_tag
        self.sc.hub.r_CN_NInit=rN
        self.sc.hub.v_CN_NInit=vN
        self.leader = leading_sat
        self.comm_module = leading_sat.comm_module
    
    def setup_comm(self):
        self.comm_module.addSpacecraftToModel(self.sc.scStateOutMsg)
    
    def sendMessage(self, message,timeSent):
        self.writeOut(MessageData(message, timeSent, self.leader))
        self.leader.writeIn(MessageData(message, timeSent, self))
        
    def writeIn(self, message:MessageData):
        self.messageInHistory.append(message)
    def writeOut(self, message:MessageData):
        self.messageOutHistory.append(message)

class LeadingSatellite:
    def __init__(self, index, rN, vN,gravFactory, model_tag="LeadingSat"):
        self.index = index
        self.model_tag = model_tag
        self.messageInHistory: list[MessageData] = []
        self.messageOutHistory: list[MessageData] = []
        self.sc = spacecraft.Spacecraft()
        gravFactory.addBodiesTo(self.sc)
        self.sc.ModelTag = self.model_tag
        self.sc.hub.r_CN_NInit=rN
        self.sc.hub.v_CN_NInit=vN
        self.aHat_B=[0.0,0.0,0.0]
        self.comm_module = spacecraftLocation.SpacecraftLocation()
        self.comm_module.ModelTag = f"CommCheck_{model_tag}"
        self.children = []

    def setup_comm(self, earth_radius,aHat_B, maximum_range):
        self.aHat_B=aHat_B
        self.comm_module.primaryScStateInMsg.subscribeTo(self.sc.scStateOutMsg)
        self.comm_module.rEquator = earth_radius
        self.comm_module.rPolar = earth_radius * 0.98
        self.comm_module.aHat_B = aHat_B
        self.comm_module.theta = np.radians(30.0)
        self.comm_module.maximumRange = maximum_range
        
    def add_child(self, child_sat):
        self.children.append(child_sat)
    def sendMessageToLead(self, message, timeSent, lead):
        self.writeOut(MessageData(message, timeSent, lead))
        lead.writeIn(MessageData(message, timeSent, self))
    def sendMessage(self, message,timeSent, child_sat):
        if(child_sat in self.children):
            self.writeOut(MessageData(message, timeSent, child_sat))
            child_sat.writeIn(MessageData(message, timeSent, self))
    
    def writeOut(self, message:MessageData):
        self.messageOutHistory.append(message)

    def writeIn(self, message:MessageData):
        self.messageInHistory.append(message)