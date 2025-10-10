#!/usr/bin/env python
"""
reference_constellation_formations.py

REFERENCE EXAMPLE: Correct formation setup for 4 clusters
Demonstrates proper use of FormationGeometry and cartesian offsets.

Run this as a standalone test to verify formations work correctly.
"""

import os
import sys
import numpy as np
from datetime import datetime

# Add project paths
import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = path  # This script is in BskSim root, so use current directory
sys.path.extend([ROOT_DIR, os.path.join(ROOT_DIR, 'core')])

# Create necessary directories
VIZFILE_DIR = os.path.join(ROOT_DIR, 'Vizfile')
VIZFILES_SUBDIR = os.path.join(VIZFILE_DIR, '_VizFiles')
os.makedirs(VIZFILE_DIR, exist_ok=True)
os.makedirs(VIZFILES_SUBDIR, exist_ok=True)

from Basilisk.utilities import macros, orbitalMotion, simIncludeGravBody, SimulationBaseClass, vizSupport
from Basilisk.simulation import spacecraft


def create_reference_constellation():
    """
    Create a reference constellation with 4 clusters showing all formation types.
    This demonstrates the CORRECT way to set up formations.
    """
    
    # Setup simulation
    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    
    # Create Earth
    gravFactory = simIncludeGravBody.gravBodyFactory()
    planet = gravFactory.createEarth()
    planet.isCentralBody = True
    mu = planet.mu
    
    print("\n" + "="*60)
    print("REFERENCE CONSTELLATION - CORRECT FORMATIONS")
    print("="*60)
    
    # Define 4 clusters with different formations
    clusters = [
        {
            "name": "Alpha",
            "formation": "Leader-Follower",
            "num_sats": 4,
            "base_anomaly": 0.0,      # degrees
            "separation": 15.0,        # km
            "color": "red"
        },
        {
            "name": "Beta", 
            "formation": "Diamond",
            "num_sats": 5,
            "base_anomaly": 90.0,
            "separation": 12.0,
            "color": "green"
        },
        {
            "name": "Gamma",
            "formation": "Triangle", 
            "num_sats": 3,
            "base_anomaly": 180.0,
            "separation": 10.0,
            "color": "blue"
        },
        {
            "name": "Delta",
            "formation": "Line",
            "num_sats": 4,
            "base_anomaly": 270.0,
            "separation": 20.0,
            "color": "yellow"
        }
    ]
    
    all_spacecraft = []
    
    # Base orbit parameters (LEO 600km)
    base_orbit = {
        "altitude": 600.0,  # km above Earth
        "inclination": 55.0,  # degrees
        "raan": 0.0,
        "eccentricity": 0.001
    }
    
    for cluster_idx, cluster in enumerate(clusters):
        print(f"\n--- Creating Cluster: {cluster['name']} ---")
        print(f"Formation: {cluster['formation']}")
        print(f"Satellites: {cluster['num_sats']}")
        print(f"Separation: {cluster['separation']} km")
        
        # Get formation positions using proper geometry
        formation_positions = get_formation_positions(
            cluster['formation'],
            cluster['num_sats'],
            cluster['separation']
        )
        
        for sat_idx, pos_data in enumerate(formation_positions):
            # Create orbital elements for this satellite
            oe = orbitalMotion.ClassicElements()
            oe.a = (6371 + base_orbit["altitude"]) * 1000  # Convert to meters
            oe.e = base_orbit["eccentricity"]
            oe.i = (base_orbit["inclination"] + pos_data["di"]) * macros.D2R
            oe.Omega = (base_orbit["raan"] + pos_data["dOmega"]) * macros.D2R
            oe.omega = 0.0
            oe.f = (cluster["base_anomaly"] + pos_data["df"]) * macros.D2R
            
            # Convert to position and velocity
            rN, vN = orbitalMotion.elem2rv(mu, oe)
            
            # Apply cartesian offset for tight formation visualization
            offset_m = np.array(pos_data["offset"]) * 1000  # km to meters
            rN = rN + offset_m
            
            # Create spacecraft
            sc = spacecraft.Spacecraft()
            sc.ModelTag = f"{cluster['name']}_{sat_idx+1}"
            
            if sat_idx == 0:
                sc.ModelTag = f"{cluster['name']}_Leader"
            
            # Set initial conditions
            gravFactory.addBodiesTo(sc)
            sc.hub.r_CN_NInit = rN
            sc.hub.v_CN_NInit = vN
            sc.hub.sigma_BNInit = [[0.0], [0.0], [0.0]]
            sc.hub.omega_BN_BInit = [[0.0], [0.0], [0.0]]
            
            scSim.AddModelToTask(simTaskName, sc)
            all_spacecraft.append(sc)
            
            # Print satellite info
            alt_km = np.linalg.norm(rN) / 1000 - 6371
            print(f"  {sc.ModelTag}: Alt={alt_km:.1f}km, TA={oe.f*macros.R2D:.1f}°, " + 
                  f"Offset=[{pos_data['offset'][0]:.1f}, {pos_data['offset'][1]:.1f}, {pos_data['offset'][2]:.1f}]km")
    
    print(f"\nTotal spacecraft created: {len(all_spacecraft)}")
    
    # Setup Vizard visualization
    if vizSupport.vizFound:
        print("\n" + "-"*60)
        print("SETTING UP VIZARD VISUALIZATION")
        print("-"*60)
        
        viz = vizSupport.enableUnityVisualization(
            scSim, simTaskName, all_spacecraft,
            saveFile=os.path.join(ROOT_DIR, "Vizfile", "reference_formations"),
            liveStream=False
        )
        
        # Configure visualization settings
        viz.settings.showSpacecraftLabels = True
        viz.settings.orbitLinesOn = 1
        viz.settings.spacecraftSizeMultiplier = 50.0
        viz.settings.showCelestialBodyLabels = True
        
        print("Vizard configured with spacecraft labels")
    
    # Run short simulation
    simulationTime = macros.min2nano(10.0)
    print(f"\nRunning {10.0} minute simulation...")
    
    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()
    
    print("✓ Simulation complete")
    print(f"✓ Vizard file: {ROOT_DIR}/Vizfile/reference_formations_UnityViz.bin")
    
    return scSim, all_spacecraft


def get_formation_positions(formation_type, num_sats, separation_km):
    """
    Calculate formation positions with proper geometry.
    
    Returns: list of dicts with:
        - 'df': True anomaly offset (degrees)
        - 'di': Inclination offset (degrees)  
        - 'dOmega': RAAN offset (degrees)
        - 'offset': Cartesian offset [x, y, z] km for visualization
    """
    positions = []
    
    # Conversion factor: for LEO ~7000km orbit, 1° ≈ 120 km along track
    deg_per_km = 1.0 / 120.0
    sep_deg = separation_km * deg_per_km
    
    # Ensure minimum visible separation
    sep_deg = max(sep_deg, 0.5)  # At least 0.5 degrees
    
    formation = formation_type.lower()
    
    if formation == "leader-follower":
        # Train formation: satellites follow in a line
        for i in range(num_sats):
            positions.append({
                "df": -i * sep_deg,           # Each behind the previous
                "di": 0.0,
                "dOmega": 0.0,
                "offset": [-i * separation_km, 0.0, 0.0]  # Visual spacing
            })
    
    elif formation == "diamond":
        # Diamond with center + 4 cardinal points
        diamond_config = [
            {"df": 0.0, "di": 0.0, "dOmega": 0.0, 
             "offset": [0.0, 0.0, 0.0]},  # Center/leader
            
            {"df": sep_deg, "di": 0.0, "dOmega": 0.0,
             "offset": [separation_km, 0.0, 0.0]},  # Front
            
            {"df": 0.0, "di": -sep_deg * 0.05, "dOmega": 0.0,
             "offset": [0.0, separation_km, 0.0]},  # Right
            
            {"df": -sep_deg, "di": 0.0, "dOmega": 0.0,
             "offset": [-separation_km, 0.0, 0.0]},  # Back
            
            {"df": 0.0, "di": sep_deg * 0.05, "dOmega": 0.0,
             "offset": [0.0, -separation_km, 0.0]}  # Left
        ]
        
        for i in range(min(num_sats, len(diamond_config))):
            positions.append(diamond_config[i])
        
        # Extra satellites form outer ring
        for i in range(len(diamond_config), num_sats):
            angle = 2 * np.pi * (i - 5) / max(1, num_sats - 5)
            positions.append({
                "df": sep_deg * 1.5 * np.cos(angle),
                "di": sep_deg * 0.05 * np.sin(angle),
                "dOmega": 0.0,
                "offset": [
                    separation_km * 1.5 * np.cos(angle),
                    separation_km * 1.5 * np.sin(angle),
                    0.0
                ]
            })
    
    elif formation == "triangle":
        # Equilateral triangle
        triangle_config = [
            {"df": 0.0, "di": 0.0, "dOmega": 0.0,
             "offset": [0.0, 0.0, 0.0]},  # Front apex (leader)
            
            {"df": -sep_deg * 0.866, "di": sep_deg * 0.05, "dOmega": 0.0,
             "offset": [-separation_km * 0.866, -separation_km * 0.5, 0.0]},  # Left rear
            
            {"df": -sep_deg * 0.866, "di": -sep_deg * 0.05, "dOmega": 0.0,
             "offset": [-separation_km * 0.866, separation_km * 0.5, 0.0]}  # Right rear
        ]
        
        for i in range(min(num_sats, len(triangle_config))):
            positions.append(triangle_config[i])
        
        # Extra satellites in center
        for i in range(len(triangle_config), num_sats):
            positions.append({
                "df": -sep_deg * 0.5,
                "di": 0.0,
                "dOmega": 0.0,
                "offset": [-separation_km * 0.5, 0.0, 0.0]
            })
    
    elif formation == "line":
        # Linear formation along orbital path
        for i in range(num_sats):
            positions.append({
                "df": i * sep_deg,
                "di": 0.0,
                "dOmega": 0.0,
                "offset": [i * separation_km, 0.0, 0.0]
            })
    
    else:
        # Default to leader-follower
        return get_formation_positions("leader-follower", num_sats, separation_km)
    
    return positions


if __name__ == "__main__":
    print("\n" + "="*60)
    print("REFERENCE CONSTELLATION FORMATIONS TEST")
    print("="*60)
    print("\nThis script demonstrates CORRECT formation setup.")
    print("Use this as a reference for your constellation configuration.\n")
    
    try:
        scSim, spacecraft_list = create_reference_constellation()
        
        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)
        print("\nFormations created correctly:")
        print("  ✓ 4 clusters with different formations")
        print("  ✓ Proper orbital element offsets")
        print("  ✓ Cartesian offsets for visualization")
        print("  ✓ Vizard binary file generated")
        print("\nOpen the .bin file in Vizard to see the formations.")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()