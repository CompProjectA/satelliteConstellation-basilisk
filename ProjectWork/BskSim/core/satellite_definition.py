# satellite_definition.py
import numpy as np
from Basilisk.simulation import spacecraft, spacecraftLocation, simpleBattery
from Basilisk.architecture import messaging
from message_data import MessageData

class ChildSatellite:
    """Child satellite that shares the leader's comm module and can exchange messages with the leader."""

    def __init__(self, index, rN, vN, gravFactory, leading_sat):
        self.index = index
        self.model_tag = f"Satellite{index+1}"
        self.sc = spacecraft.Spacecraft()
        gravFactory.addBodiesTo(self.sc)
        self.sc.ModelTag = self.model_tag
        self.sc.hub.r_CN_NInit = rN
        self.sc.hub.v_CN_NInit = vN

        # message logs
        self.messageInHistory = []
        self.messageOutHistory = []

        # comm
        self.leader = leading_sat
        self.comm_module = leading_sat.comm_module  # share leader's location model

        # simple battery (kept lightweight; you can tune thresholds in GUI later)
        self.battery = simpleBattery.SimpleBattery()
        self.battery.ModelTag = f"satBattery_{self.index}"
        self.battery.storageCapacity = 100.0
        self.battery.storedCharge_Init = 50.0
        self.batteryReader = messaging.PowerStorageStatusMsgReader()
        self.batteryReader.subscribeTo(self.battery.batPowerOutMsg)


    def setup_comm(self):
        # index of the next accessOutMsg before adding this child
        idx = len(self.comm_module.accessOutMsgs)
        self.comm_module.addSpacecraftToModel(self.sc.scStateOutMsg)
        self.access_idx = idx

    # --- messaging API ---
    def sendMessage(self, message, timeSent):
        """Send message from child -> leader"""
        md = MessageData(message, timeSent, self.leader)
        self.writeOut(md)
        self.leader.writeIn(MessageData(message, timeSent, self))

    def writeIn(self, message: MessageData):
        self.messageInHistory.append(message)

    def writeOut(self, message: MessageData):
        self.messageOutHistory.append(message)


class LeadingSatellite:
    """Leader satellite with its own comm module and a list of children."""

    def __init__(self, index, rN, vN, gravFactory, model_tag="LeadingSat"):
        self.index = index
        self.model_tag = model_tag
        self.sc = spacecraft.Spacecraft()
        gravFactory.addBodiesTo(self.sc)
        self.sc.ModelTag = self.model_tag
        self.sc.hub.r_CN_NInit = rN
        self.sc.hub.v_CN_NInit = vN

        # message logs
        self.messageInHistory = []
        self.messageOutHistory = []

        # comm module owned by leader
        self.comm_module = spacecraftLocation.SpacecraftLocation()
        self.comm_module.ModelTag = f"Comm_{model_tag}"
        self.children = []

    def setup_comm(self, earth_radius, aHat_B, maximum_range):
        self.comm_module.primaryScStateInMsg.subscribeTo(self.sc.scStateOutMsg)
        self.comm_module.rEquator = earth_radius
        self.comm_module.rPolar = earth_radius * 0.98
        self.comm_module.aHat_B = aHat_B
        self.comm_module.theta = np.radians(30.0)  # ~60° total FOV
        self.comm_module.maximumRange = maximum_range

    def add_child(self, child_sat):
        self.children.append(child_sat)

    # --- messaging API ---
    def sendMessageToLead(self, message, timeSent, lead):
        """Leader -> Leader"""
        md = MessageData(message, timeSent, lead)
        self.writeOut(md)
        lead.writeIn(MessageData(message, timeSent, self))

    def sendMessage(self, message, timeSent, child_sat):
        """Leader -> Child"""
        if child_sat in self.children:
            md = MessageData(message, timeSent, child_sat)
            self.writeOut(md)
            child_sat.writeIn(MessageData(message, timeSent, self))

    def writeOut(self, message: MessageData):
        self.messageOutHistory.append(message)

    def writeIn(self, message: MessageData):
        self.messageInHistory.append(message)
