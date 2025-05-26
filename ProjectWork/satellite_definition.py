# satellite_network.py

from Basilisk.simulation import spacecraft, spacecraftLocation
from Basilisk.utilities import orbitalMotion, macros

class LeadingSatellite:
    def __init__(self, index, mu, model_tag="LeadingSat"):
        self.index = index
        self.model_tag = model_tag
        self.sc = spacecraft.Spacecraft()
        self.sc.ModelTag = self.model_tag
        self.comm_module = spacecraftLocation.SpacecraftLocation()
        self.comm_module.ModelTag = f"CommCheck_{model_tag}"
        self.mu = mu
        self.children = []

    def set_orbit(self, oe):
        oe.f = self.index * 0.0  
        rN, vN = orbitalMotion.elem2rv(self.mu, oe)
        self.sc.hub.r_CN_NInit = rN
        self.sc.hub.v_CN_NInit = vN

    def setup_comm(self, earth_radius):
        self.comm_module.primaryScStateInMsg.subscribeTo(self.sc.scStateOutMsg)
        self.comm_module.rEquator = earth_radius
        self.comm_module.rPolar = earth_radius * 0.98
        self.comm_module.aHat_B = [1.0, 0.0, 0.0]
        self.comm_module.theta = macros.D2R * 180.0
        self.comm_module.maximumRange = 2e7

    def add_child(self, child_sat):
        self.children.append(child_sat)
        self.comm_module.addSpacecraftToModel(child_sat.sc.scStateOutMsg)


class ChildSatellite:
    def __init__(self, index, mu):
        self.index = index
        self.model_tag = f"ChildSat{index+1}"
        self.sc = spacecraft.Spacecraft()
        self.sc.ModelTag = self.model_tag
        self.comm_module = spacecraftLocation.SpacecraftLocation()
        self.comm_module.ModelTag = f"CommCheck_{self.model_tag}"
        self.mu = mu

    def set_orbit(self, oe, phase_deg):
        oe.f = phase_deg * macros.D2R
        rN, vN = orbitalMotion.elem2rv(self.mu, oe)
        self.sc.hub.r_CN_NInit = rN
        self.sc.hub.v_CN_NInit = vN

    def setup_comm(self, earth_radius):
        self.comm_module.primaryScStateInMsg.subscribeTo(self.sc.scStateOutMsg)
        self.comm_module.rEquator = earth_radius
        self.comm_module.rPolar = earth_radius * 0.98
        self.comm_module.aHat_B = [1.0, 0.0, 0.0]
        self.comm_module.theta = macros.D2R * 180.0
        self.comm_module.maximumRange = 2e7
