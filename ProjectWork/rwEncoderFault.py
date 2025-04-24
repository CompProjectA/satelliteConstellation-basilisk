from Basilisk.architecture import messaging
from Basilisk.architecture import sysModel

class rwEncoderFault(sysModel.SysModel):
    def __init__(self):
        super(rwEncoderFault, self).__init__()
        self.ModelTag = "rwEncoderFault"
        self.rwInMsg = messaging.RWSpeedMsgReader()
        self.rwOutMsg = messaging.RWSpeedMsg()
        self.enableFault = False
        self.stuckSpeed = 0.0
        self.faultyRW = 0

    def Reset(self, currentTimeNanos):
        print("[DEBUG] rwEncoderFault Reset called")

    def UpdateState(self, currentTimeNanos):
        data = self.rwInMsg.read()
        newData = messaging.RWSpeedMsgPayload()
        newData = data
        if self.enableFault:
            print(f"[DEBUG] Fault injecting stuck speed on RW {self.faultyRW}")
            newData.wheelSpeeds[self.faultyRW] = self.stuckSpeed
        self.rwOutMsg.write(newData, currentTimeNanos)