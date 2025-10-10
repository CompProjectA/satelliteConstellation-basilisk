#!/usr/bin/env python
"""
formation_geometry.py
Proper formation flying implementation for satellite clusters.
This replaces the scattered orbital positioning with tight formations.
"""

import numpy as np
from Basilisk.utilities import orbitalMotion, macros
from typing import List
from dataclasses import dataclass

class FormationGeometry:
    """Calculate proper formation positions for satellite clusters"""
    
    @staticmethod
    def calculate_formation_positions(num_sats, formation_type, leader_orbit, separation_km=10.0, mu=3.986e14):
        """
        Calculate relative positions for satellites in a formation.
        
        Parameters:
        - num_sats: Number of satellites in the formation
        - formation_type: Type of formation ("Circle", "Line", "Diamond", "Triangle", "Leader-Follower")
        - leader_orbit: Dict with orbital elements {a, e, i, Omega, omega, f} in SI units
        - separation_km: Separation distance between satellites in km
        - mu: Gravitational parameter
        
        Returns:
        - List of dicts with {rN, vN} for each satellite
        """
        positions_velocities = []
        
        # Convert separation to meters
        separation_m = separation_km * 1000.0
        
        # Leader orbital elements
        oe_leader = orbitalMotion.ClassicElements()
        oe_leader.a = leader_orbit['a']
        oe_leader.e = leader_orbit['e']
        oe_leader.i = leader_orbit['i']
        oe_leader.Omega = leader_orbit['Omega']
        oe_leader.omega = leader_orbit['omega']
        oe_leader.f = leader_orbit['f']
        
        # Get leader position/velocity
        rN_leader, vN_leader = orbitalMotion.elem2rv(mu, oe_leader)
        
        # Calculate orbital period and mean motion
        T = 2 * np.pi * np.sqrt(oe_leader.a**3 / mu)
        n = 2 * np.pi / T  # Mean motion (rad/s)
        
        # Normalize formation type
        formation = formation_type.lower().strip()
        
        if formation in ["circle", "ring"]:
            # Circle/Ring formation in the orbital plane
            positions_velocities = FormationGeometry._calculate_circle_formation(
                num_sats, oe_leader, separation_m, mu
            )
            
        elif formation in ["line", "column", "along-track"]:
            # Line formation along the velocity vector
            positions_velocities = FormationGeometry._calculate_line_formation(
                num_sats, oe_leader, separation_m, mu
            )
            
        elif formation in ["diamond", "box", "square"]:
            # Diamond/Box formation in the orbital plane
            positions_velocities = FormationGeometry._calculate_diamond_formation(
                num_sats, oe_leader, separation_m, mu
            )
            
        elif formation in ["triangle"]:
            # Triangle formation
            positions_velocities = FormationGeometry._calculate_triangle_formation(
                num_sats, oe_leader, separation_m, mu
            )
            
        elif formation in ["leader-follower", "train", "trail"]:
            # Leader-Follower formation (along-track separation)
            positions_velocities = FormationGeometry._calculate_leader_follower_formation(
                num_sats, oe_leader, separation_m, mu
            )
            
        else:
            # Default to circle formation
            print(f"Unknown formation type '{formation_type}', defaulting to circle")
            positions_velocities = FormationGeometry._calculate_circle_formation(
                num_sats, oe_leader, separation_m, mu
            )
        
        return positions_velocities
    
    @staticmethod
    def _calculate_circle_formation(num_sats, oe_leader, separation_m, mu):
        """Circle formation in the orbital plane"""
        positions = []
        
        # Leader at center
        rN_leader, vN_leader = orbitalMotion.elem2rv(mu, oe_leader)
        positions.append({'rN': rN_leader, 'vN': vN_leader})
        
        if num_sats > 1:
            # Place other satellites in a circle around the leader
            # Use Hill's equations for relative motion
            radius_deg = separation_m / (oe_leader.a * np.pi / 180)  # Convert to degrees of true anomaly
            
            for i in range(1, num_sats):
                angle = 2 * np.pi * (i - 1) / (num_sats - 1)
                
                # Small variations in true anomaly and inclination
                oe_sat = orbitalMotion.ClassicElements()
                oe_sat.a = oe_leader.a
                oe_sat.e = oe_leader.e
                oe_sat.i = oe_leader.i + (separation_m / oe_leader.a) * np.sin(angle) * 0.001
                oe_sat.Omega = oe_leader.Omega
                oe_sat.omega = oe_leader.omega
                oe_sat.f = oe_leader.f + radius_deg * np.cos(angle) * macros.D2R
                
                rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
                positions.append({'rN': rN, 'vN': vN})
        
        return positions
    
    @staticmethod
    def _calculate_line_formation(num_sats, oe_leader, separation_m, mu):
        """Line formation along the velocity vector"""
        positions = []
        
        # Calculate spacing in true anomaly
        # For along-track separation, we use small differences in true anomaly
        spacing_deg = (separation_m / oe_leader.a) * (180 / np.pi)
        
        # Center the formation around the leader position
        start_offset = -(num_sats - 1) * spacing_deg / 2
        
        for i in range(num_sats):
            oe_sat = orbitalMotion.ClassicElements()
            oe_sat.a = oe_leader.a
            oe_sat.e = oe_leader.e
            oe_sat.i = oe_leader.i
            oe_sat.Omega = oe_leader.Omega
            oe_sat.omega = oe_leader.omega
            oe_sat.f = oe_leader.f + (start_offset + i * spacing_deg) * macros.D2R
            
            rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
            positions.append({'rN': rN, 'vN': vN})
        
        return positions
    
    @staticmethod
    def _calculate_diamond_formation(num_sats, oe_leader, separation_m, mu):
        """Diamond/Box formation in the orbital plane"""
        positions = []
        
        # Leader at center
        rN_leader, vN_leader = orbitalMotion.elem2rv(mu, oe_leader)
        positions.append({'rN': rN_leader, 'vN': vN_leader})
        
        if num_sats > 1:
            # Diamond pattern offsets
            spacing_deg = (separation_m / oe_leader.a) * (180 / np.pi)
            
            # Define diamond positions (along-track and cross-track)
            diamond_offsets = [
                (spacing_deg, 0),      # Front
                (0, spacing_deg*0.5),  # Right
                (-spacing_deg, 0),     # Back
                (0, -spacing_deg*0.5), # Left
            ]
            
            for i in range(1, min(num_sats, 5)):
                along_track, cross_track = diamond_offsets[i-1]
                
                oe_sat = orbitalMotion.ClassicElements()
                oe_sat.a = oe_leader.a
                oe_sat.e = oe_leader.e
                oe_sat.i = oe_leader.i + cross_track * 0.001 * macros.D2R
                oe_sat.Omega = oe_leader.Omega + cross_track * 0.001 * macros.D2R
                oe_sat.omega = oe_leader.omega
                oe_sat.f = oe_leader.f + along_track * macros.D2R
                
                rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
                positions.append({'rN': rN, 'vN': vN})
            
            # If more than 5 satellites, add them in a second layer
            if num_sats > 5:
                for i in range(5, num_sats):
                    angle = 2 * np.pi * (i - 5) / max(1, num_sats - 5)
                    
                    oe_sat = orbitalMotion.ClassicElements()
                    oe_sat.a = oe_leader.a
                    oe_sat.e = oe_leader.e
                    oe_sat.i = oe_leader.i + spacing_deg * 0.002 * np.sin(angle) * macros.D2R
                    oe_sat.Omega = oe_leader.Omega
                    oe_sat.omega = oe_leader.omega
                    oe_sat.f = oe_leader.f + spacing_deg * 2 * np.cos(angle) * macros.D2R
                    
                    rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
                    positions.append({'rN': rN, 'vN': vN})
        
        return positions
    
    @staticmethod
    def _calculate_triangle_formation(num_sats, oe_leader, separation_m, mu):
        """Triangle formation"""
        positions = []
        
        spacing_deg = (separation_m / oe_leader.a) * (180 / np.pi)
        
        # Triangle vertices
        if num_sats >= 1:
            # Leader at front vertex
            rN_leader, vN_leader = orbitalMotion.elem2rv(mu, oe_leader)
            positions.append({'rN': rN_leader, 'vN': vN_leader})
        
        if num_sats >= 2:
            # Left vertex
            oe_sat = orbitalMotion.ClassicElements()
            oe_sat.a = oe_leader.a
            oe_sat.e = oe_leader.e
            oe_sat.i = oe_leader.i + spacing_deg * 0.0005 * macros.D2R
            oe_sat.Omega = oe_leader.Omega - spacing_deg * 0.0005 * macros.D2R
            oe_sat.omega = oe_leader.omega
            oe_sat.f = oe_leader.f - spacing_deg * macros.D2R
            
            rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
            positions.append({'rN': rN, 'vN': vN})
        
        if num_sats >= 3:
            # Right vertex
            oe_sat = orbitalMotion.ClassicElements()
            oe_sat.a = oe_leader.a
            oe_sat.e = oe_leader.e
            oe_sat.i = oe_leader.i - spacing_deg * 0.0005 * macros.D2R
            oe_sat.Omega = oe_leader.Omega + spacing_deg * 0.0005 * macros.D2R
            oe_sat.omega = oe_leader.omega
            oe_sat.f = oe_leader.f - spacing_deg * macros.D2R
            
            rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
            positions.append({'rN': rN, 'vN': vN})
        
        # Additional satellites go in the center
        if num_sats > 3:
            for i in range(3, num_sats):
                oe_sat = orbitalMotion.ClassicElements()
                oe_sat.a = oe_leader.a
                oe_sat.e = oe_leader.e
                oe_sat.i = oe_leader.i
                oe_sat.Omega = oe_leader.Omega
                oe_sat.omega = oe_leader.omega
                oe_sat.f = oe_leader.f - spacing_deg * 0.5 * macros.D2R
                
                rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
                positions.append({'rN': rN, 'vN': vN})
        
        return positions
    
    @staticmethod
    def _calculate_leader_follower_formation(num_sats, oe_leader, separation_m, mu):
        """Leader-Follower formation (train configuration)"""
        positions = []
        
        # Calculate spacing in true anomaly for along-track separation
        spacing_deg = (separation_m / oe_leader.a) * (180 / np.pi)
        
        for i in range(num_sats):
            oe_sat = orbitalMotion.ClassicElements()
            oe_sat.a = oe_leader.a
            oe_sat.e = oe_leader.e
            oe_sat.i = oe_leader.i
            oe_sat.Omega = oe_leader.Omega
            oe_sat.omega = oe_leader.omega
            # Each satellite trails the previous one
            oe_sat.f = oe_leader.f - i * spacing_deg * macros.D2R
            
            rN, vN = orbitalMotion.elem2rv(mu, oe_sat)
            positions.append({'rN': rN, 'vN': vN})
        
        return positions


# Example usage for creating 4 clusters with different formations
def create_constellation_with_formations():
    """
    Example function showing how to create a constellation with 4 clusters,
    each using a different formation type.
    """
    from Basilisk.utilities import simIncludeGravBody
    from satellite_definition import LeadingSatellite, ChildSatellite
    
    # Earth parameters
    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    mu = earth.mu
    
    clusters = []
    satellites = []
    
    # Define 4 clusters with different formations
    cluster_configs = [
        {
            "name": "Cluster1",
            "num_sats": 4,
            "formation": "Diamond",
            "leader_orbit": {
                'a': 7000e3,  # 7000 km altitude
                'e': 0.001,
                'i': 55 * macros.D2R,
                'Omega': 0 * macros.D2R,
                'omega': 0 * macros.D2R,
                'f': 45 * macros.D2R
            },
            "separation": 10.0  # 10 km separation
        },
        {
            "name": "Cluster2", 
            "num_sats": 4,
            "formation": "Line",
            "leader_orbit": {
                'a': 7000e3,
                'e': 0.001,
                'i': 55 * macros.D2R,
                'Omega': 0 * macros.D2R,
                'omega': 0 * macros.D2R,
                'f': 225 * macros.D2R
            },
            "separation": 15.0  # 15 km separation
        },
        {
            "name": "Cluster3",
            "num_sats": 4,
            "formation": "Triangle",
            "leader_orbit": {
                'a': 7000e3,
                'e': 0.001,
                'i': 56 * macros.D2R,
                'Omega': 0 * macros.D2R,
                'omega': 0 * macros.D2R,
                'f': 135 * macros.D2R
            },
            "separation": 12.0
        },
        {
            "name": "Cluster4",
            "num_sats": 4,
            "formation": "Circle",
            "leader_orbit": {
                'a': 7000e3,
                'e': 0.001,
                'i': 56 * macros.D2R,
                'Omega': 0 * macros.D2R,
                'omega': 0 * macros.D2R,
                'f': 315 * macros.D2R
            },
            "separation": 8.0
        }
    ]
    
    # Create each cluster
    for cluster_config in cluster_configs:
        # Calculate formation positions
        positions = FormationGeometry.calculate_formation_positions(
            num_sats=cluster_config["num_sats"],
            formation_type=cluster_config["formation"],
            leader_orbit=cluster_config["leader_orbit"],
            separation_km=cluster_config["separation"],
            mu=mu
        )
        
        # Create satellites for this cluster
        cluster_sats = []
        leader = None
        
        for i, pos_vel in enumerate(positions):
            if i == 0:
                # Create leader
                sat = LeadingSatellite(
                    index=len(satellites),
                    rN=pos_vel['rN'],
                    vN=pos_vel['vN'],
                    gravFactory=gravFactory,
                    model_tag=f"{cluster_config['name']}_Leader"
                )
                leader = sat
                satellites.append(sat)
                cluster_sats.append(sat)
            else:
                # Create child
                sat = ChildSatellite(
                    index=len(satellites),
                    rN=pos_vel['rN'],
                    vN=pos_vel['vN'],
                    gravFactory=gravFactory,
                    leading_sat=leader
                )
                leader.add_child(sat)
                satellites.append(sat)
                cluster_sats.append(sat)
        
        clusters.append({
            "config": cluster_config,
            "satellites": cluster_sats,
            "leader": leader
        })
    
    return clusters, satellites


if __name__ == "__main__":
    # Test the formation geometry calculations
    test_orbit = {
        'a': 7000e3,
        'e': 0.001,
        'i': 55 * macros.D2R,
        'Omega': 0,
        'omega': 0,
        'f': 0
    }
    
    formations = ["Circle", "Line", "Diamond", "Triangle", "Leader-Follower"]
    
    for formation in formations:
        print(f"\n{formation} Formation (4 satellites, 10 km separation):")
        positions = FormationGeometry.calculate_formation_positions(
            num_sats=4,
            formation_type=formation,
            leader_orbit=test_orbit,
            separation_km=10.0
        )
        
        for i, pos_vel in enumerate(positions):
            r_mag = np.linalg.norm(pos_vel['rN'])
            v_mag = np.linalg.norm(pos_vel['vN'])
            alt_km = (r_mag - 6371e3) / 1000
            print(f"  Sat {i+1}: Alt={alt_km:.1f} km, V={v_mag/1000:.3f} km/s")


@dataclass
class CartesianMember:
    name: str
    role: str     # 'leader' or 'child'
    offset_km: list  # [dx, dy, dz] from leader

class CartesianFormations:
    @staticmethod
    def ring(n:int, sep_km:float) -> list[CartesianMember]:
        members = [CartesianMember("Leader", "leader", [0.0, 0.0, 0.0])]
        if n <= 1: return members
        R = max(2.0, sep_km)
        for k in range(1, n):
            th = 2.0 * np.pi * (k-1) / max(1, (n-1))
            members.append(CartesianMember(f"Sat{k+1}", "child",
                                           [R*np.cos(th), R*np.sin(th), 0.0]))
        return members

    @staticmethod
    def column(n:int, sep_km:float) -> list[CartesianMember]:
        members = [CartesianMember("Leader", "leader", [0.0, 0.0, 0.0])]
        d = max(2.0, sep_km*0.5)
        for k in range(1, n):
            members.append(CartesianMember(f"Sat{k+1}", "child", [k*d, 0.2*k, 0.0]))
        return members

    @staticmethod
    def box(n:int, sep_km:float) -> list[CartesianMember]:
        members = [CartesianMember("Leader", "leader", [0.0, 0.0, 0.0])]
        if n == 1: return members
        side = max(3.0, sep_km)
        ring1 = [( side, 0, 0),(0, side, 0),(-side, 0, 0),(0,-side, 0)]
        ring2 = [( 2*side,  2*side, 0),(-2*side,  2*side, 0),
                 (-2*side, -2*side, 0),( 2*side, -2*side, 0)]
        coords = ring1 + ring2
        for i,(dx,dy,dz) in enumerate(coords[:max(0,n-1)], start=2):
            members.append(CartesianMember(f"Sat{i}", "child", [dx,dy,dz]))
        r3 = 1.5*side
        k = len(members)
        while k < n:
            ang = 2.0*np.pi*(k-1)/max(1,n-1)
            members.append(CartesianMember(f"Sat{k}", "child",
                                           [r3*np.cos(ang), r3*np.sin(ang), 0.0]))
            k += 1
        return members

    @staticmethod
    def train(n:int, sep_km:float) -> list[CartesianMember]:
        members = [CartesianMember("Leader", "leader", [0.0, 0.0, 0.0])]
        d = max(8.0, sep_km)
        for k in range(1, n):
            members.append(CartesianMember(f"Sat{k+1}", "child", [k*d, 0.0, 0.0]))
        return members
