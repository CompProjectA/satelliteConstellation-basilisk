#!/usr/bin/env python
"""
cluster_integration.py
Cluster/formation integration built around LeadingSatellite and ChildSatellite.
- Creates clusters from config
- Supports individual S/C
- Intra-cluster comm (leader <-> children)
- Optional inter-cluster comm (leader <-> leader)
- Simple fault event hook (log at time)

Drop this file in the same folder as satellite_definition.py and message_data.py
"""

import numpy as np
from Basilisk.simulation import spacecraft, spacecraftLocation
from Basilisk.utilities import macros, orbitalMotion

# --- local modules (same folder) ---
try:
    from satellite_definition import LeadingSatellite, ChildSatellite
    from message_data import MessageData
    SATELLITE_CLASSES_AVAILABLE = True
except ImportError:
    SATELLITE_CLASSES_AVAILABLE = False
    print("Warning: Could not import satellite_definition / message_data")

class ClusterManager:
    """Manages satellite clusters with leader-child hierarchy using proper classes"""

    def __init__(self):
        self.clusters = {}
        self.leaders = []
        self.all_satellites = []           # list[spacecraft.Spacecraft]
        self.communication_modules = []    # list[spacecraftLocation.SpacecraftLocation]
        self.access_recorders = []

    # ---- core creators ----
    def create_cluster_from_config(self, cluster_config, satellites_config, scSim, simTaskName, gravFactory, mu):
        """
        cluster_config: {'name': str, 'formation': 'Tight'|'Loose'|...}
        satellites_config: list of sat configs for the cluster
        """
        cluster_name = cluster_config['name']
        cluster_sats = []
        leader = None
        children = []

        # ensure leader first
        satellites_config.sort(key=lambda x: 0 if x.get('role') == 'leader' else 1)

        sat_index = len(self.all_satellites)

        for sat_cfg in satellites_config:
            # classic elements -> r,v
            oe = orbitalMotion.ClassicElements()
            oe.a = sat_cfg['orbit']['a'] * 1000.0
            oe.e = sat_cfg['orbit']['e']
            oe.i = sat_cfg['orbit']['i'] * np.pi / 180.0
            oe.Omega = sat_cfg['orbit']['Omega'] * np.pi / 180.0
            oe.omega = sat_cfg['orbit']['omega'] * np.pi / 180.0
            oe.f = sat_cfg['orbit']['f'] * np.pi / 180.0
            rN, vN = orbitalMotion.elem2rv(mu, oe)

            # optional position offset (formation)
            if 'position_offset' in sat_cfg:
                rN = rN + np.array(sat_cfg['position_offset'], dtype=float) * 1000.0

            if SATELLITE_CLASSES_AVAILABLE and sat_cfg.get('role') == 'leader':
                # leader
                leader = LeadingSatellite(
                    index=sat_index, rN=rN, vN=vN, gravFactory=gravFactory, model_tag=sat_cfg['name']
                )
                earth_radius = 6371000.0
                comm_range = sat_cfg['communication']['range'] * 1000.0
                aHat_B = sat_cfg['communication']['aHat_B']
                leader.setup_comm(earth_radius, aHat_B, comm_range)

                self.leaders.append(leader)
                cluster_sats.append(leader)
                self.all_satellites.append(leader.sc)
                scSim.AddModelToTask(simTaskName, leader.sc)
                scSim.AddModelToTask(simTaskName, leader.comm_module)

                if sat_cfg.get('fault', {}).get('enabled'):
                    self._apply_fault_to_satellite(leader.sc, sat_cfg['fault'], scSim)

            elif SATELLITE_CLASSES_AVAILABLE:
                # child
                if leader is None:
                    print(f"Error: No leader found for child '{sat_cfg['name']}' in cluster '{cluster_name}'")
                    continue

                child = ChildSatellite(index=sat_index, rN=rN, vN=vN, gravFactory=gravFactory, leading_sat=leader)
                child.setup_comm()

                leader.add_child(child)
                children.append(child)
                cluster_sats.append(child)
                self.all_satellites.append(child.sc)

                scSim.AddModelToTask(simTaskName, child.sc)

                if sat_cfg.get('fault', {}).get('enabled'):
                    self._apply_fault_to_satellite(child.sc, sat_cfg['fault'], scSim)
            else:
                # fallback: plain spacecraft if classes missing
                sc = spacecraft.Spacecraft()
                sc.ModelTag = sat_cfg['name']
                gravFactory.addBodiesTo(sc)
                sc.hub.r_CN_NInit = rN
                sc.hub.v_CN_NInit = vN
                scSim.AddModelToTask(simTaskName, sc)
                cluster_sats.append(sc)
                self.all_satellites.append(sc)

            sat_index += 1

        self.clusters[cluster_name] = {
            'config': cluster_config,
            'satellites': cluster_sats,
            'leader': leader,
            'children': children
        }
        return cluster_sats

    def create_individual_satellite(self, sat_cfg, scSim, simTaskName, gravFactory, mu):
        """Create a non-cluster spacecraft."""
        oe = orbitalMotion.ClassicElements()
        oe.a = sat_cfg['orbit']['a'] * 1000.0
        oe.e = sat_cfg['orbit']['e']
        oe.i = sat_cfg['orbit']['i'] * np.pi / 180.0
        oe.Omega = sat_cfg['orbit']['Omega'] * np.pi / 180.0
        oe.omega = sat_cfg['orbit']['omega'] * np.pi / 180.0
        oe.f = sat_cfg['orbit']['f'] * np.pi / 180.0
        rN, vN = orbitalMotion.elem2rv(mu, oe)

        sc = spacecraft.Spacecraft()
        sc.ModelTag = sat_cfg['name']
        gravFactory.addBodiesTo(sc)
        sc.hub.r_CN_NInit = rN
        sc.hub.v_CN_NInit = vN
        sc.hub.sigma_BNInit = [[0.01], [0.02], [-0.01]]
        sc.hub.omega_BN_BInit = [[0.0001], [-0.0002], [0.0001]]

        scSim.AddModelToTask(simTaskName, sc)
        self.all_satellites.append(sc)

        if sat_cfg.get('fault', {}).get('enabled'):
            self._apply_fault_to_satellite(sc, sat_cfg['fault'], scSim)

        return sc

    # ---- comms between cluster leaders ----
    def setup_inter_cluster_communication(self, scSim, simTaskName, earth_radius=6371000.0):
        """Create wide-FOV comm location models on each leader and connect to other leaders."""
        if len(self.leaders) < 2:
            return

        for i, leader in enumerate(self.leaders):
            scLoc = spacecraftLocation.SpacecraftLocation()
            scLoc.ModelTag = f"InterCluster_{leader.model_tag}"
            scLoc.rEquator = earth_radius
            scLoc.rPolar = earth_radius * 0.98
            scLoc.aHat_B = [-0.2, -0.4, -0.4]
            scLoc.theta = np.radians(45.0)
            scLoc.maximumRange = 5_000e3

            scLoc.primaryScStateInMsg.subscribeTo(leader.sc.scStateOutMsg)
            # add the *other* leaders as targets
            for j, other in enumerate(self.leaders):
                if i == j:
                    continue
                scLoc.addSpacecraftToModel(other.sc.scStateOutMsg)

            self.communication_modules.append(scLoc)
            scSim.AddModelToTask(simTaskName, scLoc)

            # record all access outputs exposed by this model
            for k in range(len(scLoc.accessOutMsgs)):
                rec = scLoc.accessOutMsgs[k].recorder()
                rec.ModelTag = f"AccessRec_{leader.model_tag}_{k}"
                self.access_recorders.append(rec)
                scSim.AddModelToTask(simTaskName, rec)

    # ---- utilities ----
    def _apply_fault_to_satellite(self, sc: spacecraft.Spacecraft, fault_cfg, scSim):
        """Simple event that logs a fault trigger at the configured minute timestamp."""
        sc.faultConfig = dict(fault_cfg)
        sc.faultInjected = False

        fault_time_nano = macros.min2nano(fault_cfg.get('time', 0.0))
        event_name = f"fault_{sc.ModelTag}"

        # Just log + disable the event. You can hook real fault logic here.
        scSim.createNewEvent(
            event_name,
            macros.sec2nano(1.0),
            True,
            [f"self.TotalSim.CurrentNanos >= {fault_time_nano}"],
            [
                (
                    "print("
                    f"'FAULT: {fault_cfg.get('type','unknown')} for {sc.ModelTag} at', "
                    "self.TotalSim.CurrentNanos * 1.6666666666666667e-11, 'minutes')"
                ),
                f"self.setEventActivity('{event_name}', False)"
            ]
        )

    def get_all_spacecraft_objects(self):
        """Return a flat list of spacecraft objects (clusters + individuals)."""
        s_list = []
        # clusters
        for cdata in self.clusters.values():
            for sat in cdata['satellites']:
                if hasattr(sat, 'sc'):
                    s_list.append(sat.sc)
                else:
                    s_list.append(sat)
        # uniques
        for sc in self.all_satellites:
            if sc not in s_list:
                s_list.append(sc)
        return s_list

    def send_message_in_cluster(self, cluster_name, message_content, time_min, from_leader=True, to_child_index=0):
        """Send a message within a cluster using the classes' messaging API."""
        cluster = self.clusters.get(cluster_name)
        if not cluster or not cluster.get('leader'):
            print(f"Cluster '{cluster_name}' not found or has no leader.")
            return False

        leader = cluster['leader']
        children = cluster.get('children', [])
        if from_leader:
            if 0 <= to_child_index < len(children):
                leader.sendMessage(message_content, time_min, children[to_child_index])
                return True
        else:
            if 0 <= to_child_index < len(children):
                children[to_child_index].sendMessage(message_content, time_min)
                return True
        return False

    def send_inter_cluster_message(self, from_cluster, to_cluster, message_content, time_min):
        """Send a message between leaders of two clusters."""
        a = self.clusters.get(from_cluster, {}).get('leader')
        b = self.clusters.get(to_cluster, {}).get('leader')
        if a and b:
            a.sendMessageToLead(message_content, time_min, b)
            return True
        return False

    def get_communication_status(self):
        status = {
            'clusters': len(self.clusters),
            'leaders': len(self.leaders),
            'total_satellites': len(self.all_satellites),
            'comm_modules': len(self.communication_modules),
            'access_recorders': len(self.access_recorders),
        }
        for name, data in self.clusters.items():
            leader = data['leader']
            status[f'{name}_children'] = len(getattr(leader, 'children', [])) if leader else 0
        return status

    def get_message_history(self):
        """Gather sent/received logs from leaders and children."""
        hist = {'sent': [], 'received': []}

        for leader in self.leaders:
            for msg in getattr(leader, 'messageOutHistory', []):
                hist['sent'].append({
                    'from': leader.model_tag,
                    'to': getattr(msg.objectActive, 'model_tag', 'Unknown'),
                    'content': msg.message_content,
                    'time': msg.timeSent
                })
            for msg in getattr(leader, 'messageInHistory', []):
                hist['received'].append({
                    'to': leader.model_tag,
                    'from': getattr(msg.objectActive, 'model_tag', 'Unknown'),
                    'content': msg.message_content,
                    'time': msg.timeSent
                })

        for cdata in self.clusters.values():
            for child in cdata.get('children', []):
                for msg in getattr(child, 'messageOutHistory', []):
                    hist['sent'].append({
                        'from': child.model_tag,
                        'to': getattr(msg.objectActive, 'model_tag', 'Unknown'),
                        'content': msg.message_content,
                        'time': msg.timeSent
                    })
                for msg in getattr(child, 'messageInHistory', []):
                    hist['received'].append({
                        'to': child.model_tag,
                        'from': getattr(msg.objectActive, 'model_tag', 'Unknown'),
                        'content': msg.message_content,
                        'time': msg.timeSent
                    })
        return hist


def integrate_clusters_with_simulation(config, cluster_manager, scSim, simTaskName, gravFactory, mu):
    """
    Build clusters and individuals from a SimulationConfig-like object with .spacecraft_list
    Each spacecraft entry should contain:
      - name, orbit{a,e,i,Omega,omega,f}, type ('cluster_member' or other), role ('leader'|'child'), cluster, communication{range,aHat_B}, fault{enabled,...}
    """
    all_spacecraft = []

    # group cluster members vs individuals
    grouped = {}
    singles = []
    for sat_cfg in getattr(config, 'spacecraft_list', []):
        if sat_cfg.get('type') == 'cluster_member':
            grouped.setdefault(sat_cfg.get('cluster'), []).append(sat_cfg)
        else:
            singles.append(sat_cfg)

    # create clusters
    for cname, sats in grouped.items():
        cluster_cfg = {'name': cname, 'formation': sats[0].get('formation', 'Tight') if sats else 'Tight'}
        cluster_manager.create_cluster_from_config(cluster_cfg, sats, scSim, simTaskName, gravFactory, mu)

    # individuals
    for sat_cfg in singles:
        sc = cluster_manager.create_individual_satellite(sat_cfg, scSim, simTaskName, gravFactory, mu)
        all_spacecraft.append(sc)

    # inter-cluster comm if needed
    if len(cluster_manager.clusters) > 1:
        cluster_manager.setup_inter_cluster_communication(scSim, simTaskName)

    # final list for viz etc.
    return cluster_manager.get_all_spacecraft_objects()
