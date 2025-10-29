#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
spacecraft_simulation.py

COMPLETE MERGED VERSION - Core simulation functionality for the Spacecraft Constellation Fault Simulator.

This file merges:
- Your original enhanced engine (fault loader, REAL ML, plot APIs new/old, ClusterManager integration)
- Claude's “Enhanced 2.0” visualization helpers (coverage cones & verification, color-matched target connections)
- Compatibility helpers for both *lightweight* and *rich* cluster formats

Key features (all preserved; nothing removed):
- Cluster integration (ClusterManager) when available
- 4-cluster validation + leader checks (Claude feature)
- Enhanced Vizard styling per cluster (color/size/labels) (Claude feature)
- Enhanced target coverage cones + satellite-target lines (Claude feature)
- Optional enhanced plot generation via ConstellationPlotter (Claude feature)
- Real fault simulation for EACH satellite (fault loader)
- Battery visualization panels (if battery faults enabled)
- Full plotting coverage (supports both old/new plots APIs)
- Detailed summary + binary verification + cameras/targets
- REAL ML detection via fault_detection_router (supports LSTM + Isolation Forest)
"""

import os
import sys
import inspect
import numpy as np
from datetime import datetime
import traceback

import matplotlib
matplotlib.use('Agg')  # Headless plotting OK (GUI embeds set backend in their own modules)
import matplotlib.pyplot as plt
import logging

# --- Compact logging / spam control ---
# Quiet common noisy libs
logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# Simulation logger with duplicate filter
sim_logger = logging.getLogger('simulation')
sim_logger.setLevel(logging.INFO)

class DuplicateFilter(logging.Filter):
    """Filter out duplicate log messages (only every 10th repeat)."""
    def __init__(self):
        super().__init__()
        self.last = None
        self.count = 0
    def filter(self, record):
        key = (record.module, record.levelno, record.msg)
        if key == self.last:
            self.count += 1
            return (self.count % 10) == 0
        self.last = key
        self.count = 0
        return True

sim_logger.addFilter(DuplicateFilter())


# ----- Path setup -----
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
FAULTS_DIR = os.path.join(ROOT_DIR, 'faults')
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plots')
LOGS_DIR = os.path.join(ROOT_DIR, 'logs')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, FAULTS_DIR, MODELS_DIR, PLOTTING_DIR])

# ----- Basilisk modules -----
try:
    from Basilisk.utilities import macros, vizSupport
    from Basilisk.utilities import orbitalMotion
    from Basilisk.utilities import simIncludeGravBody
    from Basilisk.simulation import spacecraft
    from Basilisk.utilities import SimulationBaseClass
    from Basilisk import __path__ as bsk_path
    BASILISK_AVAILABLE = True
    BSK_INSTALL_PATH = bsk_path[0]
except Exception as e:
    print(f"ERROR: Could not import Basilisk modules: {e}")
    print("Make sure Basilisk is properly installed.")
    BASILISK_AVAILABLE = False
    BSK_INSTALL_PATH = None

# ----- REAL ML fault detection integration -----
# Real ML fault detection integration with ROUTER
try:
    from fault_detection_router import run_fault_detection_on_scenario
    REAL_ML_AVAILABLE = True
    print("Fault detection router available (supports LSTM + Isolation Forest)")
except ImportError as e:
    REAL_ML_AVAILABLE = False
    print(f"Fault detection router not available: {e}")

# ----- Cluster integration (existing feature) -----
try:
    from cluster_integration import ClusterManager, integrate_clusters_with_simulation
    CLUSTER_INTEGRATION_AVAILABLE = True
    print("Cluster integration module loaded successfully")
except Exception as e:
    print(f"Warning: Could not import cluster integration: {e}")
    CLUSTER_INTEGRATION_AVAILABLE = False

# ----- Plotting imports - support BOTH old and new APIs -----
PLOTS_AVAILABLE = False
PLOTS_API = {"new": False, "old": False}

# Try new plots API
try:
    from plots import (
        generate_constellation_overview_plots,     # new
        generate_cluster_communication_plots,      # new
        generate_inter_satellite_distance_plots,   # new
        generate_fault_plots                       # both
    )
    PLOTS_AVAILABLE = True
    PLOTS_API["new"] = True
    print("Plotting module (new API) loaded successfully")
except Exception as e:
    print(f"Warning: Could not import new plotting API: {e}")

# Try old plots API additions
try:
    from plots import (
        generate_constellation_plots,              # old
        create_fault_config_for_real_simulation    # old (CRITICAL for real sim)
    )
    PLOTS_AVAILABLE = True
    PLOTS_API["old"] = True
    print("Plotting module (old API) loaded successfully")
except Exception as e:
    print(f"Note: Old plotting API not fully available: {e}")

# Optional enhanced plotter (Claude-style extras)
CONSTELLATION_PLOTTER_AVAILABLE = False
try:
    from plots import ConstellationPlotter
    CONSTELLATION_PLOTTER_AVAILABLE = True
    print("ConstellationPlotter available")
except Exception:
    pass

# ---------------------------
# Fault loader / module checks
# ---------------------------
FAULT_IMPORT_STATUS = {}
AVAILABLE_FAULTS = []
FAILED_FAULTS = []
FAULT_LOADER_AVAILABLE = False


def check_fault_modules(verbose=False):
    """Check and report status of all fault modules (quiet by default)."""
    global FAULT_IMPORT_STATUS, AVAILABLE_FAULTS, FAILED_FAULTS, FAULT_LOADER_AVAILABLE

    FAULT_IMPORT_STATUS, AVAILABLE_FAULTS, FAILED_FAULTS = {}, [], []

    if verbose:
        print("\n" + "="*60)
        print("FAULT MODULE IMPORT STATUS")
        print("="*60)

    try:
        from fault_loader import (
            get_fault_scenario_class, create_scenario, run_scenario,
            apply_fault_to_spacecraft, extract_fault_data_from_scenario,
            get_available_fault_types
        )
        FAULT_LOADER_AVAILABLE = True
        if verbose:
            print("Fault loader imported successfully")
        try:
            AVAILABLE_FAULTS = get_available_fault_types()
            if verbose:
                print(f"Available fault types: {AVAILABLE_FAULTS}")
        except Exception as e:
            if verbose:
                print(f"Could not get available fault types: {e}")
    except ImportError as e:
        FAULT_LOADER_AVAILABLE = False
        if verbose:
            print(f"Fault loader import failed: {e}")

    # Check individual modules silently unless verbose
    fault_modules_to_check = [
        ("friction_fault", "FrictionFaultScenario"),
        ("powerlimit_fault", "PowerLimitFaultScenario"),
        ("encoder_fault", "EncoderFaultScenario"),
        ("battery_fault", "BatteryFaultScenario")
    ]
    for module_name, class_name in fault_modules_to_check:
        try:
            module = __import__(f"faults.{module_name}", fromlist=[class_name])
            ok = hasattr(module, class_name)
            FAULT_IMPORT_STATUS[module_name] = ok
            if not ok:
                FAILED_FAULTS.append(module_name)
        except Exception:
            FAULT_IMPORT_STATUS[module_name] = False
            FAILED_FAULTS.append(module_name)

    if verbose:
        print("="*60 + "\n")

    return FAULT_LOADER_AVAILABLE, AVAILABLE_FAULTS, FAILED_FAULTS



def check_plots_module(verbose=True):
    """Verify plots module availability."""
    global PLOTS_AVAILABLE

    if PLOTS_API["new"] or PLOTS_API["old"]:
        if verbose:
            apis = []
            if PLOTS_API["new"]:
                apis.append("new")
            if PLOTS_API["old"]:
                apis.append("old")
            print(f"Plots module is available (APIs: {', '.join(apis)})")
        PLOTS_AVAILABLE = True
    else:
        if verbose:
            print("Plots module NOT available")
        PLOTS_AVAILABLE = False

    return PLOTS_AVAILABLE


def inject_fault_into_spacecraft(sc_object, fault_config, current_time_nano):
    """
    Inject/log fault into a spacecraft based on configuration.
    Returns True if injection conditions were met at current_time_nano.
    """
    if not fault_config.get("enabled", False):
        return False

    fault_type = fault_config["type"]
    fault_magnitude = fault_config["magnitude"]
    fault_wheel = fault_config["wheel"]
    fault_time_nano = macros.min2nano(fault_config["time"])

    if current_time_nano < fault_time_nano:
        return False  # Not time yet

    if fault_type not in AVAILABLE_FAULTS and FAULT_LOADER_AVAILABLE:
        print(f"ERROR: Fault type '{fault_type}' is not available")
        print(f"Available types: {AVAILABLE_FAULTS}")
        return False

    if FAULT_LOADER_AVAILABLE:
        try:
            if not hasattr(sc_object, 'fault_log'):
                sc_object.fault_log = []
            sc_object.fault_log.append({
                'type': fault_type,
                'magnitude': fault_magnitude,
                'wheel': fault_wheel,
                'time': current_time_nano
            })
            print(f"Fault logged: {fault_type} on wheel {fault_wheel} at {current_time_nano*macros.NANO2MIN:.2f} min")
            return True
        except Exception as e:
            print(f"Could not inject/log fault: {e}")
            return False
    else:
        print(f"Fault loader not available for {fault_type}")
        return False


def calculate_target_visibility(spacecraft_pos, target_pos, earth_radius=6371000):
    """Simple LOS check to a target (kept for downstream use)."""
    try:
        sc_pos = np.array(spacecraft_pos)
        tgt_pos = np.array(target_pos)
        sc_to_target = tgt_pos - sc_pos
        if np.linalg.norm(sc_to_target) == 0:
            return False
        sc_altitude = np.linalg.norm(sc_pos) - earth_radius
        if sc_altitude > 200000:  # 200 km
            earth_center = np.array([0, 0, 0])
            sc_to_earth = earth_center - sc_pos
            cos_angle = np.dot(sc_to_target, sc_to_earth) / (np.linalg.norm(sc_to_target) * np.linalg.norm(sc_to_earth))
            elevation_angle = np.pi/2 - np.arccos(np.clip(cos_angle, -1, 1))
            return elevation_angle > np.radians(10)
        return False
    except Exception:
        return False


# ---------------------------------------------------------------------
#  Target & Config classes (names/shape kept for GUI import consistency)
# ---------------------------------------------------------------------
class TargetDefinition:
    """Class to hold target location definitions."""
    def __init__(self, name, latitude, longitude, color="red"):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.color = color
        self.assigned_to = []

    def to_dict(self):
        return {
            "name": self.name,
            "lat": self.latitude,
            "lon": self.longitude,
            "color": self.color,
            "assigned_to": self.assigned_to
        }


class SimulationConfig:
    """Configuration class for simulation parameters."""
    def __init__(self):
        # Default simulation time
        self.simulation_time = 30.0  # minutes

        # Spacecraft list for constellation support
        self.spacecraft_list = []

        # Default fault magnitudes for proper scaling
        self.DEFAULT_FAULT_MAGNITUDES = {
            "friction": 0.0005,    # N·m
            "power_limit": 0.5,    # W
            "encoder": 20.0,       # %
            "battery": 50.0        # W
        }

        # Legacy single-sat fault parameters (back-compat if used elsewhere)
        self.fault_type = "friction"
        self.fault_magnitude = 0.0005
        self.fault_wheel_number = 3
        self.fault_time = 10.0  # minutes

        # Periodic fault parameters
        self.enable_periodic_fault = False
        self.periodic_fault_interval = 360  # seconds
        self.periodic_fault_magnitude = 0.1
        self.periodic_fault_wheel = 1

        # Output settings
        self.binary_filename = "spacecraft_viz"
        self.show_plots = True
        self.save_plots = True
        self.save_binary = True

        # Camera configuration
        self.camera_position = [0.0, 0.0, 10.0]
        self.camera_fov = 80.0
        self.active_camera_name = None

        # Targets (defaults; GUI may overwrite)
        self.targets = [
            TargetDefinition("Melbourne", -37.8136, 144.9631, "red"),
            TargetDefinition("New York", 40.71, -74.00, "blue"),
            TargetDefinition("Tokyo", 35.68, 139.77, "green"),
            TargetDefinition("London", 51.51, -0.13, "yellow")
        ]
        # Optional: self.clusters (either lightweight or rich) can be attached by GUI

    def validate(self):
        """Validate configuration parameters."""
        if not BASILISK_AVAILABLE:
            raise ValueError("Basilisk is not available. Please install Basilisk to run simulations.")

        if self.simulation_time <= 0:
            raise ValueError("Simulation time must be positive (in minutes)")

        if self.spacecraft_list:
            for sc in self.spacecraft_list:
                # Orbit keys must exist and be in km/deg as expected by the GUI
                if sc["orbit"]["a"] <= 6371:
                    raise ValueError(f"Spacecraft {sc['name']}: Semi-major axis must be greater than Earth radius (6371 km)")
                if sc["orbit"]["e"] < 0 or sc["orbit"]["e"] >= 1.0:
                    raise ValueError(f"Spacecraft {sc['name']}: Eccentricity must be between 0 and 1")

                # Fault parameters if enabled
                if sc["fault"]["enabled"]:
                    if sc["fault"]["time"] >= self.simulation_time:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time must be less than simulation time")
                    if sc["fault"]["time"] < 0:
                        raise ValueError(f"Spacecraft {sc['name']}: Fault time cannot be negative")

                    if FAULT_LOADER_AVAILABLE and sc["fault"]["type"] not in AVAILABLE_FAULTS:
                        print(f"Warning: Spacecraft {sc['name']}: Fault type '{sc['fault']['type']}' may not be available")

                    if sc["fault"]["wheel"] not in range(4):
                        raise ValueError(f"Spacecraft {sc['name']}: Fault wheel number must be between 0 and 3")


# ---------------------------
# Helpers for cluster plotting (existing)
# ---------------------------
def _name_to_index_map(sc_objects):
    """Return dict mapping ModelTag -> index for spacecraft list."""
    mapping = {}
    for idx, sc in enumerate(sc_objects):
        tag = getattr(sc, "ModelTag", None)
        if tag:
            mapping[tag] = idx
    return mapping


def _cluster_data_from_manager(cluster_manager, sc_objects):
    """
    Build {cluster_name: {'leader': idx, 'children': [idx,...]}} from ClusterManager.
    """
    cluster_data = {}
    if not cluster_manager:
        return cluster_data

    name_index = _name_to_index_map(sc_objects)
    for cname, cinfo in cluster_manager.clusters.items():
        leader_obj = cinfo.get('leader')
        child_objs = cinfo.get('children', [])
        leader_idx = None
        child_indices = []

        if leader_obj is not None:
            ltag = getattr(leader_obj, 'model_tag', None)
            leader_idx = name_index.get(ltag, None)

        for ch in child_objs:
            ctag = getattr(ch, 'model_tag', None)
            idx = name_index.get(ctag, None)
            if idx is not None:
                child_indices.append(idx)

        if leader_idx is not None and child_indices:
            cluster_data[cname] = {'leader': leader_idx, 'children': child_indices}

    return cluster_data


def _cluster_data_from_config(config, sc_objects):
    """
    Build {cluster_name: {'leader': idx, 'children': [idx,...]}} from config alone.
    """
    cluster_data = {}
    if not hasattr(config, "spacecraft_list"):
        return cluster_data

    name_index = _name_to_index_map(sc_objects)

    for sc in config.spacecraft_list:
        cname = sc.get('cluster')
        if cname:
            cluster_data.setdefault(cname, {'leader': None, 'children': []})

    for sc in config.spacecraft_list:
        cname = sc.get('cluster')
        if not cname:
            continue
        idx = name_index.get(sc.get('name'))
        if idx is None:
            continue
        if sc.get('role') == 'leader':
            cluster_data[cname]['leader'] = idx
        elif sc.get('role') == 'child':
            cluster_data[cname]['children'].append(idx)

    cluster_data = {
        k: v for k, v in cluster_data.items()
        if v['leader'] is not None and len(v['children']) > 0
    }
    return cluster_data


# ---------------------------
# Light→Rich cluster unification (ours) + Claude-compatible builders
# ---------------------------
ALLOWED_FORMATIONS = {"Line", "Triangle", "Diamond", "Leader-Follower"}


def _build_rich_from_lightweight(lightweight, sc_objects, sc_list):
    """Convert lightweight cluster map to 'rich' form with names and indices."""
    rich = {}
    idx_to_name = {i: getattr(sc, "ModelTag", f"SC{i}") for i, sc in enumerate(sc_objects)}

    # Best-effort metadata per sat from spacecraft_list
    meta = {sc["name"]: {
        "formation": sc.get("formation"),
        "orbit": sc.get("orbit_name", "LEO 600km"),
        "separation": sc.get("separation", 10.0)}
        for sc in sc_list
    }

    for cname, data in (lightweight or {}).items():
        leader_idx = data.get("leader")
        child_idx = data.get("children", [])
        sats, leader_entry, children = [], None, []
        if leader_idx is not None:
            lname = idx_to_name.get(leader_idx, f"SC{leader_idx}")
            m = meta.get(lname, {})
            leader_entry = {"name": lname, "role": "leader", "index": leader_idx,
                            "formation": m.get("formation"), "orbit": m.get("orbit"), "separation": m.get("separation")}
            sats.append(leader_entry)
        for ci in child_idx:
            n = idx_to_name.get(ci, f"SC{ci}")
            m = meta.get(n, {})
            entry = {"name": n, "role": "child", "index": ci,
                     "formation": m.get("formation"), "orbit": m.get("orbit"), "separation": m.get("separation")}
            sats.append(entry)
            children.append(entry)
        if sats:
            formation = (leader_entry or {}).get("formation") or "Line"
            rich[cname] = {
                "satellites": sats,
                "leader": leader_entry or sats[0],
                "children": children,
                "formation": formation,
                "orbit": (leader_entry or {}).get("orbit") or "LEO 600km",
                "separation": (leader_entry or {}).get("separation") or 10.0
            }
    return rich


def _normalize_rich_clusters(config, sc_objects):
    """
    Returns a unified 'rich' cluster dict if clusters provided (either already rich
    or converted from lightweight); empty dict if none.
    """
    clusters = getattr(config, "clusters", {}) or {}
    # Already 'rich'?
    if clusters and all(isinstance(v, dict) and "satellites" in v for v in clusters.values()):
        return clusters
    # Lightweight?
    if clusters and all(isinstance(v, dict) and "leader" in v for v in clusters.values()):
        return _build_rich_from_lightweight(clusters, sc_objects, getattr(config, "spacecraft_list", []))
    return {}


def validate_cluster_configuration_rich(rich_clusters):
    """Enforce max 4 clusters + exactly one leader + optional formation check."""
    if not rich_clusters:
        return
    if len(rich_clusters) > 4:
        raise ValueError(f"ERROR: Maximum 4 clusters allowed per project requirements. Got {len(rich_clusters)}")
    for cname, c in rich_clusters.items():
        sats = c.get("satellites", [])
        leader_count = sum(1 for s in sats if s.get("role") == "leader")
        if leader_count != 1:
            raise ValueError(f"ERROR: Cluster '{cname}' must have exactly 1 leader, found {leader_count}")
        formation = c.get("formation")
        if formation and formation not in ALLOWED_FORMATIONS:
            raise ValueError(f"ERROR: Invalid formation '{formation}' for cluster '{cname}'")


def _print_cluster_summary(rich_clusters):
    if not rich_clusters:
        return
    print("\nCLUSTER CONFIGURATION:")
    print("-" * 40)
    parts = []
    for name, data in rich_clusters.items():
        leader = (data.get("leader") or {}).get("name", "?")
        children = len(data.get("children", []))
        formation = data.get("formation", "-")
        parts.append(f"{name}[Leader={leader}, {children} children, {formation}]")
        print(f"  {name}:")
        print(f"    Formation: {formation}")
        print(f"    Satellites: {len(data.get('satellites', []))}")
        print(f"    Leader: {leader}")
        print(f"    Separation: {data.get('separation', 10.0)} km")
    print(f"\nClusters({len(rich_clusters)}): {', '.join(parts)}")
    print("="*60 + "\n")


# ---------- Claude's lightweight builders (exposed for compatibility) ----------
def build_rich_clusters_for_sim(config, sc_index_map):
    """Build rich cluster dictionary for simulation (Claude’s helper retained)."""
    rich = {}

    clusters_src = getattr(config, 'clusters', {})

    if isinstance(clusters_src, dict):
        items = clusters_src.items()
    elif isinstance(clusters_src, list):
        items = [(c.get('name', f"Cluster{i}"), c) for i, c in enumerate(clusters_src)]
    else:
        return {}

    for cname, cluster_data in items:
        if isinstance(cluster_data, dict) and 'leader' in cluster_data:
            leader_idx = cluster_data['leader']
            child_indices = cluster_data.get('children', [])

            satellites = []
            leader_entry = None

            for idx in [leader_idx] + child_indices:
                if idx < len(config.spacecraft_list):
                    sc_cfg = config.spacecraft_list[idx]
                    entry = {
                        'name': sc_cfg['name'],
                        'role': 'leader' if idx == leader_idx else 'child',
                        'index': sc_index_map.get(sc_cfg['name'], idx),
                        'formation': cluster_data.get('formation', sc_cfg.get('formation', 'Line')),
                        'orbit': sc_cfg.get('orbit_name', 'Unknown'),
                        'separation': cluster_data.get('separation', sc_cfg.get('separation', 10.0))
                    }
                    satellites.append(entry)
                    if entry['role'] == 'leader':
                        leader_entry = entry

            rich[cname] = {
                'satellites': satellites,
                'leader': leader_entry,
                'children': [s for s in satellites if s['role'] == 'child'],
                'formation': cluster_data.get('formation', 'Line'),
                'orbit': cluster_data.get('orbit', 'Unknown'),
                'separation': cluster_data.get('separation', 10.0)
            }

    return rich


def validate_cluster_configuration(rich_clusters):
    """Validate cluster configuration (Claude’s version retained)."""
    if not rich_clusters:
        return

    if len(rich_clusters) > 4:
        raise ValueError(f"Maximum 4 clusters allowed. Got {len(rich_clusters)}")

    for cname, c in rich_clusters.items():
        sats = c.get("satellites", [])
        leader_count = sum(1 for s in sats if s.get("role") == "leader")
        if leader_count != 1:
            raise ValueError(f"Cluster '{cname}' must have exactly 1 leader, found {leader_count}")


# ---------------------------
# Enhanced Visualization (merged)
# ---------------------------
def enhance_vizard_settings(viz, rich_clusters):
    """
    Apply per-cluster visualization with:
    - Larger spacecraft (leaders bigger)
    - Formation lines between cluster members
    - Distinct cluster colors
    - Enhanced target visibility
    """
    if not viz or not rich_clusters:
        return viz

    cluster_colors = {
        0: [1.0, 0.1, 0.1, 1.0],   # Bright Red
        1: [0.1, 0.9, 0.9, 1.0],   # Bright Cyan
        2: [0.2, 1.0, 0.2, 1.0],   # Bright Green
        3: [1.0, 0.8, 0.0, 1.0]    # Bright Yellow/Orange
    }

    try:
        cidx = 0
        for cname, cdata in rich_clusters.items():
            color = cluster_colors.get(cidx, [0.8, 0.8, 0.8, 1.0])

            leader_idx = None
            child_indices = []

            for sat in cdata.get('satellites', []):
                idx = sat.get('index', 0)
                role = sat.get('role', 'child')

                try:
                    viz.scData[idx].scColor = color
                except Exception:
                    pass

                try:
                    if role == 'leader':
                        viz.scData[idx].spacecraftSize = 400.0
                        leader_idx = idx
                        viz.scData[idx].modelDictionaryKey = f"{cname}_L"
                    else:
                        viz.scData[idx].spacecraftSize = 200.0
                        child_indices.append(idx)
                        sat_name = sat.get('name', f"SC{idx}")
                        viz.scData[idx].modelDictionaryKey = f"{cname}_{sat_name[-4:]}"
                except Exception as e:
                    print(f"Could not set size for SC {idx}: {e}")

                try:
                    viz.scData[idx].showOrbitLine = True
                    viz.scData[idx].orbitColor = color
                    viz.scData[idx].orbitLineWidth = 6.0
                    viz.scData[idx].orbitLineAlpha = 0.9
                except Exception:
                    pass

            cidx += 1

        # Global settings
        try:
            viz.settings.showSpacecraftLabels = 1
            viz.settings.spacecraftLabelSize = 16
            viz.settings.spacecraftLabelOffset = [0, 0, 50]
            viz.settings.orbitLinesOn = 1
            viz.settings.trueTrajectoryLinesOn = 1

            viz.settings.showFormationLines = 1
            viz.settings.formationLineWidth = 5.0
            viz.settings.formationLineAlpha = 0.85
            viz.settings.formationLineColor = [1.0, 1.0, 0.0, 0.85]

            viz.settings.spacecraftCSOn = 1
            viz.settings.spacecraftCSsize = 100.0

            viz.settings.showLocationLabels = 1
            viz.settings.locationLabelSize = 14
            viz.settings.showLocationCones = 1
            viz.settings.locationConeColor = [0.3, 0.8, 0.3, 0.5]
            viz.settings.showLocationCommLines = 1
            viz.settings.locationCommLineWidth = 3.0

            viz.settings.ambientLightOn = 1
            viz.settings.ambientBrightness = 0.3

            viz.settings.planetCSon = 1
            viz.settings.viewCelestialBodiesAsPoints = 0
            viz.settings.showOrbitGrid = 0

            print("\nEnhanced Vizard settings applied:")
            print("  - Spacecraft sizes: Leaders=400, Children=200")
            print("  - Orbit lines: width=6.0, alpha=0.9")
            print("  - Formation lines: enabled, width=5.0")
            print("  - Coverage cones & target lines: enabled")

        except Exception as e:
            print(f"Could not set global viz settings: {e}")

    except Exception as e:
        print(f"Vizard enhance error (non-fatal): {e}")

    return viz


# ----- Enhanced target/coverage utilities from Claude (added; not replacing) -----
def calculate_ground_target_position(lat_deg, lon_deg, altitude_m=100000.0):
    """Calculate target position with specified altitude (meters)."""
    earth_radius = 6371000.0
    lat_rad = lat_deg * macros.D2R
    lon_rad = lon_deg * macros.D2R
    radius = earth_radius + altitude_m

    x = radius * np.cos(lat_rad) * np.cos(lon_rad)
    y = radius * np.cos(lat_rad) * np.sin(lon_rad)
    z = radius * np.sin(lat_rad)

    return [x, y, z]


def _unit(v):
    """Return unit vector for list-like v, fall back to [0,0,1] if zero-length."""
    arr = np.array(v, dtype=float)
    n = np.linalg.norm(arr)
    if n < 1e-12:
        return [0.0, 0.0, 1.0]
    return (arr / n).tolist()


def add_enhanced_targets_to_vizard(viz, config, sc_objects, rich_clusters, verbose=True):
    """
    Add targets with complete coverage visualization:
    - Large markers, coverage cones
    - Color matching to cluster/satellite
    - Connection lines enabled
    """
    if not viz:
        return {"success": False, "error": "No viz object"}

    cluster_colors = {
        0: [1.0, 0.2, 0.2, 1.0],
        1: [0.2, 1.0, 0.8, 1.0],
        2: [0.3, 1.0, 0.3, 1.0],
        3: [1.0, 0.9, 0.1, 1.0]
    }

    # Build satellite name → index/color map
    sat_name_to_idx = {}
    sat_idx_to_color = {}

    for idx, sc in enumerate(sc_objects):
        sat_name = getattr(sc, "ModelTag", f"SC{idx}")
        sat_name_to_idx[sat_name] = idx

        if rich_clusters:
            for cidx, (cname, cdata) in enumerate(rich_clusters.items()):
                for sat_info in cdata.get('satellites', []):
                    if sat_info.get('name') == sat_name:
                        sat_idx_to_color[idx] = cluster_colors.get(cidx, [0.8, 0.8, 0.8, 1.0])
                        break

    targets_added = 0
    coverage_report = {
        "total_targets": len(config.targets),
        "visible_targets": 0,
        "assigned_targets": 0,
        "unassigned_targets": 0,
        "target_details": []
    }

    if verbose:
        print("\n" + "="*70)
        print("TARGET COVERAGE VISUALIZATION")
        print("="*70)

    for target in config.targets:
        assigned_to = target.assigned_to if hasattr(target, 'assigned_to') else []

        if not assigned_to:
            coverage_report["unassigned_targets"] += 1
            if verbose:
                print(f"\nTarget '{target.name}': UNASSIGNED (not added to Vizard)")
            continue

        coverage_report["assigned_targets"] += 1

        try:
            target_pos = calculate_ground_target_position(
                target.latitude,
                target.longitude,
                altitude_m=100000.0
            )

            primary_sat_name = assigned_to[0]
            primary_sat_idx = sat_name_to_idx.get(primary_sat_name, 0)

            if primary_sat_idx in sat_idx_to_color:
                target_color_rgb = sat_idx_to_color[primary_sat_idx]
                # Convert [0..1] RGBA to 0..255 RGBA for Vizard color if needed
                try:
                    target_viz_color = vizSupport.toRGBA255(target_color_rgb)
                except Exception:
                    target_viz_color = target.color if hasattr(target, 'color') else 'red'
            else:
                target_viz_color = target.color if hasattr(target, 'color') else 'red'

            coverage_range = 4000000.0  # 4000 km

            # IMPORTANT FIX: use gHat_P instead of unsupported normalVector
            vizSupport.addLocation(
                viz,
                stationName=target.name,
                parentBodyName="earth",
                r_GP_P=target_pos,
                gHat_P=_unit(target_pos),            # <— boresight = surface normal
                fieldOfView=70.0 * macros.D2R,
                color=target_viz_color,
                range=coverage_range
            )

            targets_added += 1
            coverage_report["visible_targets"] += 1

            if verbose:
                sat_list = ", ".join(assigned_to[:2])
                if len(assigned_to) > 2:
                    sat_list += f" +{len(assigned_to)-2}"
                print(f"\nTarget '{target.name}':")
                print(f"  Pos: ({target.latitude:.2f}°, {target.longitude:.2f}°)")
                print(f"  Assigned to: {sat_list}")
                print(f"  Coverage: 4000 km, 70° FOV")

        except Exception as e:
            if verbose:
                print(f"\nCould not add target {target.name}: {e}")

    if verbose:
        print("\n" + "="*70)
        print("TARGET COVERAGE SUMMARY")
        print("="*70)
        print(f"Total:     {coverage_report['total_targets']}")
        print(f"Visible:   {coverage_report['visible_targets']}")
        print(f"Assigned:  {coverage_report['assigned_targets']}")
        print(f"Unassigned:{coverage_report['unassigned_targets']}")
        print("="*70 + "\n")

    coverage_report["success"] = True
    coverage_report["targets_added"] = targets_added

    return coverage_report


def enable_satellite_to_target_lines(viz, verbose=True):
    """Enable visual connection lines from satellites to targets."""
    if not viz:
        return
    try:
        viz.settings.showLocationCommLines = 1
        viz.settings.locationCommLineWidth = 4.0
        viz.settings.locationCommLineAlpha = 0.9
        viz.settings.showLocationCones = 1
        viz.settings.locationMarkerSize = 25.0

        if verbose:
            print("Satellite-to-target connection lines enabled; coverage cones shown; target markers enlarged.")
    except Exception as e:
        if verbose:
            print(f"Could not fully enable target connections: {e}")


def create_optimal_camera(viz, sc_objects, config):
    """
    Create optimally positioned camera for formation viewing.
    """
    if not viz or not sc_objects:
        return False

    camera_created = False

    if config.active_camera_name:
        for sc in sc_objects:
            if getattr(sc, "ModelTag", None) == config.active_camera_name:
                try:
                    vizSupport.createStandardCamera(
                        viz,
                        setMode=1,
                        spacecraftName=sc.ModelTag,
                        fieldOfView=config.camera_fov * macros.D2R,
                        displayName=f"{sc.ModelTag} View",
                        pointingVector_B=[0, 0, -1],
                        position_B=config.camera_position
                    )
                    camera_created = True
                    print(f"Camera created for {sc.ModelTag} at {config.camera_position}")
                    break
                except Exception as e:
                    print(f"Could not create user camera: {e}")

    if not camera_created and sc_objects:
        try:
            vizSupport.createStandardCamera(
                viz,
                setMode=1,
                spacecraftName=sc_objects[0].ModelTag,
                fieldOfView=70.0 * macros.D2R,
                displayName="Constellation Overview",
                pointingVector_B=[0, 0, -1],
                position_B=[0.0, 0.0, 200.0]
            )
            camera_created = True
            print("Camera: Constellation Overview (wide-angle)")
        except Exception as e:
            print(f"Could not create overview camera: {e}")

    if viz and not camera_created:
        try:
            vizSupport.createStandardCamera(
                viz,
                setMode=0,
                displayName="Earth View",
                fieldOfView=80.0 * macros.D2R
            )
            print("Camera: Earth-centered fallback")
            camera_created = True
        except Exception:
            pass

    return camera_created


def setup_enhanced_visualization(scSim, simTaskName, sc_objects, config, rich_clusters, viz_filename):
    """
    Complete visualization setup with all enhancements (kept for compatibility with Claude’s entrypoints).
    """
    viz = None

    if not vizSupport.vizFound:
        print("Vizard not available")
        return None

    try:
        print("\n" + "="*60)
        print("SETTING UP ENHANCED VIZARD VISUALIZATION")
        print("="*60)

        viz = vizSupport.enableUnityVisualization(
            scSim, simTaskName, sc_objects,
            saveFile=viz_filename,
            liveStream=False
        )

        print(f"Vizard initialized: {len(sc_objects)} spacecraft")

        # Apply enhanced styling
        viz = enhance_vizard_settings(viz, rich_clusters)

        # Add targets with coverage
        _ = add_enhanced_targets_to_vizard(
            viz, config, sc_objects, rich_clusters, verbose=True
        )

        # Enable connection lines
        enable_satellite_to_target_lines(viz, verbose=True)

        # Create camera
        create_optimal_camera(viz, sc_objects, config)

        print("="*60 + "\n")
        return viz

    except Exception as e:
        print(f"Visualization setup failed: {e}")
        traceback.print_exc()
        return None


# ---------------------------------------------------------------------
# Main Simulation Function (your existing flow + enhanced target handling)
# ---------------------------------------------------------------------
def run_custom_simulation(config, fault_detection_tab=None):
    """
    Run a customized simulation based on the configuration.

    Args:
        config: SimulationConfig object
        fault_detection_tab: Optional FaultDetectionTab for ML model selection

    Returns:
        tuple: (scenario, viz, figureList, output_dir, ml_results)
    """
    # Check modules (quiet)
    print("\nChecking module availability...")
    check_fault_modules(verbose=False)
    check_plots_module(verbose=False)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(LOGS_DIR, f"sim_results_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    # Validate configuration
    config.validate()
    print("Configuration validation passed")

    if not config.spacecraft_list:
        print("ERROR: No spacecraft configured in simulation")
        return None, None, {}, output_dir

    print("\n" + "="*60)
    print("SPACECRAFT CONSTELLATION SIMULATION")
    print("="*60)
    print(f"Simulation Duration: {config.simulation_time} minutes")
    print(f"Number of Spacecraft: {len(config.spacecraft_list)}")

    # Clusters?
    has_clusters = any(sc.get('cluster') for sc in config.spacecraft_list) or \
                   any(sc.get('type') == 'cluster_member' for sc in config.spacecraft_list)

    # Time setup
    simulationTime = macros.min2nano(config.simulation_time)
    print("\nTime Configuration:")
    print(f"  Input: {config.simulation_time} minutes")
    print(f"  Simulation: {simulationTime} nanoseconds")
    print(f"  Verification: {simulationTime/1e9/60:.1f} minutes")

    # Create simulation
    simTaskName = "simTask"
    simProcessName = "simProcess"
    scSim = SimulationBaseClass.SimBaseClass()
    dynProcess = scSim.CreateNewProcess(simProcessName)
    simulationTimeStep = macros.sec2nano(1.0)  # 1 second step
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    print("Simulation time step: 1 second")

    # Gravity
    gravFactory = simIncludeGravBody.gravBodyFactory()
    planet = gravFactory.createEarth()
    planet.isCentralBody = True
    mu = planet.mu
    print(f"Earth gravity parameter: {mu/1e14:.3f} x 10^14 m^3/s^2")

    # Create spacecraft
    sc_objects = []
    cluster_manager = None
    fault_spacecraft_count = 0

    if has_clusters and CLUSTER_INTEGRATION_AVAILABLE:
        print("\n" + "-"*50)
        print("USING CLUSTER INTEGRATION")
        print("-"*50)
        cluster_manager = ClusterManager()
        sc_objects = integrate_clusters_with_simulation(
            config, cluster_manager, scSim, simTaskName, gravFactory, mu
        )
        print(f"Created {len(sc_objects)} spacecraft objects using ClusterManager")
        try:
            comm_status = cluster_manager.get_communication_status()
            print(f"Communication setup: {comm_status}")
        except Exception:
            pass
    else:
        print("\n" + "-"*50)
        print("CREATING INDIVIDUAL SPACECRAFT")
        print("-"*50)

        for i, sc_config in enumerate(config.spacecraft_list):
            sc = spacecraft.Spacecraft()
            sc.ModelTag = sc_config["name"]

            # Orbit (km/deg -> m/rad)
            oe = orbitalMotion.ClassicElements()
            orbit = sc_config["orbit"]
            oe.a = orbit["a"] * 1000
            oe.e = orbit["e"]
            oe.i = orbit["i"] * macros.D2R
            oe.Omega = orbit["Omega"] * macros.D2R
            oe.omega = orbit["omega"] * macros.D2R
            oe.f = orbit["f"] * macros.D2R

            rN, vN = orbitalMotion.elem2rv(mu, oe)
            sc.hub.r_CN_NInit = rN
            sc.hub.v_CN_NInit = vN
            sc.hub.sigma_BNInit = [[0.01], [0.02], [-0.01]]
            sc.hub.omega_BN_BInit = [[0.0001], [-0.0002], [0.0001]]

            gravFactory.addBodiesTo(sc)
            scSim.AddModelToTask(simTaskName, sc)
            sc_objects.append(sc)

            # Fault setup (no per-sat chatter)
            if sc_config["fault"]["enabled"]:
                ft = sc_config["fault"]["type"]
                fm = sc_config["fault"]["magnitude"]
                if ft in config.DEFAULT_FAULT_MAGNITUDES and abs(fm - 0.0005) < 1e-12:
                    sc_config["fault"]["magnitude"] = config.DEFAULT_FAULT_MAGNITUDES[ft]
                sc.faultConfig = sc_config["fault"].copy()
                sc.faultInjected = False
                fault_spacecraft_count += 1

        print(f"Created {len(sc_objects)} spacecraft"
              + (f" | faults configured on {fault_spacecraft_count}" if fault_spacecraft_count else ""))

    # Normalize rich clusters & validate
    rich_clusters = _normalize_rich_clusters(config, sc_objects)
    if rich_clusters:
        validate_cluster_configuration_rich(rich_clusters)
        _print_cluster_summary(rich_clusters)

    # ---------------- Visualization ----------------
    viz = None
    if config.save_binary and BASILISK_AVAILABLE and vizSupport.vizFound:
        print("\n" + "-"*50)
        print("SETTING UP VIZARD VISUALIZATION")
        print("-"*50)

        viz_files_dir = os.path.join(VIZ_DIR, "_VizFiles")
        os.makedirs(viz_files_dir, exist_ok=True)

        binary_filename = os.path.basename(config.binary_filename)
        binary_full_path = os.path.join(VIZ_DIR, binary_filename)
        print(f"Binary output (base): {binary_filename}_UnityViz.bin")

        # Battery visualization if battery faults exist
        has_battery_fault = any(sat["fault"]["enabled"] and sat["fault"]["type"] == "battery"
                                for sat in config.spacecraft_list)
        gsList = [[] for _ in range(len(config.spacecraft_list))]
        battery_components = {}

        if has_battery_fault:
            try:
                from Basilisk.simulation import simpleBattery, simpleSolarPanel, simplePowerSink
                from Basilisk.architecture import messaging

                print("Battery fault detected - setting up battery visualization panels...")

                for i, sc_config in enumerate(config.spacecraft_list):
                    if sc_config["fault"]["enabled"] and sc_config["fault"]["type"] == "battery":
                        sc = sc_objects[i]

                        # Battery
                        battery = simpleBattery.SimpleBattery()
                        battery.ModelTag = f"{sc.ModelTag}_battery"
                        battery.storageCapacity = 100.0  # Wh
                        battery.storedCharge_Init = 50.0  # Wh
                        scSim.AddModelToTask(simTaskName, battery)

                        # Solar panel
                        solarPanel = simpleSolarPanel.SimpleSolarPanel()
                        solarPanel.ModelTag = f"{sc.ModelTag}_solarPanel"
                        solarPanel.setPanelParameters([-1.0, 0.0, 0.0], 0.05, 0.3)
                        solarPanel.stateInMsg.subscribeTo(sc.scStateOutMsg)
                        scSim.AddModelToTask(simTaskName, solarPanel)
                        battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)

                        # Sun message
                        sunMsgData = messaging.SpicePlanetStateMsgPayload()
                        sunMsgData.PositionVector = [-1.0, 0.0, 0.0]
                        sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
                        solarPanel.sunInMsg.subscribeTo(sunMsg)

                        # Power sinks
                        powerSink = simplePowerSink.SimplePowerSink()
                        powerSink.ModelTag = f"{sc.ModelTag}_powerSink"
                        powerSink.nodePowerOut = -0.01  # 10 W nominal
                        scSim.AddModelToTask(simTaskName, powerSink)
                        battery.addPowerNodeToModel(powerSink.nodePowerOutMsg)

                        powerSinkFault = simplePowerSink.SimplePowerSink()
                        powerSinkFault.ModelTag = f"{sc.ModelTag}_powerSinkFault"
                        powerSinkFault.nodePowerOut = 0.0  # Initially off
                        scSim.AddModelToTask(simTaskName, powerSinkFault)
                        battery.addPowerNodeToModel(powerSinkFault.nodePowerOutMsg)

                        battery_components[sc.ModelTag] = {
                            'battery': battery,
                            'solarPanel': solarPanel,
                            'powerSink': powerSink,
                            'powerSinkFault': powerSinkFault
                        }

                        # Generic storage panels
                        batteryPanel = vizSupport.vizInterface.GenericStorage()
                        batteryPanel.label = f"{sc.ModelTag} Battery (%)"
                        batteryPanel.units = "%"
                        batteryPanel.minValue = 0
                        batteryPanel.maxValue = 100
                        batteryPanel.useStorageLevel = True

                        batteryReader = messaging.PowerStorageStatusMsgReader()
                        batteryReader.subscribeTo(battery.batPowerOutMsg)
                        batteryPanel.batteryStateInMsg = batteryReader
                        batteryPanel.this.disown()

                        solarViz = vizSupport.vizInterface.GenericStorage()
                        solarViz.label = f"{sc.ModelTag} Solar Power"
                        solarViz.units = "W"
                        solarViz.minValue = 0.0
                        solarViz.maxValue = 60.0
                        solarViz.useStorageLevel = False

                        solarReader = messaging.PowerNodeUsageMsgReader()
                        solarReader.subscribeTo(solarPanel.nodePowerOutMsg)
                        solarViz.powerNodeUsageInMsg = solarReader
                        solarViz.this.disown()

                        gsList[i] = [batteryPanel, solarViz]

                        # Fault activation event
                        fault_time_nano = macros.min2nano(sc_config["fault"]["time"])
                        scSim.__dict__[f'powerSinkFault_{i}'] = powerSinkFault

                        scSim.createNewEvent(
                            f"batteryFault_{sc.ModelTag}",
                            simulationTimeStep,
                            True,
                            [f"self.TotalSim.CurrentNanos >= {fault_time_nano}"],
                            [
                                f"self.powerSinkFault_{i}.nodePowerOut = -0.05",
                                f"print('Battery fault activated for {sc.ModelTag}')",
                                f"self.setEventActivity('batteryFault_{sc.ModelTag}', False)"
                            ]
                        )

                        print(f"Set up battery visualization for {sc.ModelTag}")

            except Exception as e:
                print(f"Battery visualization setup failed (continuing without panels): {e}")

        # Enable viz (with/without generic storage panels)
        try:
            has_generic_storage = any(len(panels) > 0 for panels in gsList)
            if has_generic_storage:
                viz = vizSupport.enableUnityVisualization(
                    scSim,
                    simTaskName,
                    sc_objects,
                    saveFile=binary_full_path,
                    liveStream=False,
                    genericStorageList=gsList
                )
                # Turn on instrument GUI for those sats
                for i, sc in enumerate(sc_objects):
                    if sc.ModelTag in battery_components:
                        vizSupport.setInstrumentGuiSetting(
                            viz,
                            spacecraftName=sc.ModelTag,
                            showGenericStoragePanel=True
                        )
            else:
                viz = vizSupport.enableUnityVisualization(
                    scSim,
                    simTaskName,
                    sc_objects,
                    saveFile=binary_full_path,
                    liveStream=False
                )

            print(f"Vizard enabled for {len(sc_objects)} spacecraft")

            # Baseline viz settings
            if hasattr(viz, 'settings'):
                viz.settings.showSpacecraftLabels = True
                viz.settings.orbitLinesOn = 1
                viz.settings.spacecraftSizeMultiplier = 50.0
                viz.settings.spacecraftCSOn = True
                viz.settings.showCelestialBodyLabels = True
                viz.settings.planetCSon = 1
                viz.settings.viewCelestialBodiesAsPoints = 0
                viz.settings.showLocationCommLines = 1
                viz.settings.showLocationCones = 1
                viz.settings.showLocationLabels = 1
                print("Vizard settings configured (size multiplier: 50.0)")

            # Enhanced cluster styling
            viz = enhance_vizard_settings(viz, rich_clusters)

            # Prefer enhanced target addition with coverage/lines:
            try:
                _ = add_enhanced_targets_to_vizard(viz, config, sc_objects, rich_clusters, verbose=True)
                enable_satellite_to_target_lines(viz, verbose=True)
            except Exception as e:
                print(f"Enhanced target injection failed; falling back to simple add: {e}")
                # Fallback simple markers
                targets_added = 0
                for target in config.targets:
                    if target.assigned_to:
                        try:
                            lat, lon, color = target.latitude, target.longitude, target.color
                            alt = 50000.0  # 50 km
                            radius = 6371000.0 + alt
                            lat_rad = lat * macros.D2R
                            lon_rad = lon * macros.D2R
                            x = radius * np.cos(lat_rad) * np.cos(lon_rad)
                            y = radius * np.cos(lat_rad) * np.sin(lon_rad)
                            z = radius * np.sin(lat_rad)

                            vizSupport.addLocation(
                                viz,
                                stationName=target.name,
                                parentBodyName="earth",
                                r_GP_P=[x, y, z],
                                color=color,
                                range=2000000.0
                            )
                            targets_added += 1
                            print(f"Added target: {target.name}")
                        except Exception as e2:
                            print(f"Could not add target {target.name}: {e2}")
                print(f"Added {targets_added} targets to visualization")

            # Cameras
            _ = create_optimal_camera(viz, sc_objects, config)

        except Exception as e:
            print(f"Failed to set up visualization: {e}")
            viz = None
    else:
        if BASILISK_AVAILABLE and not vizSupport.vizFound:
            print("Vizard visualization module not found")
        if not config.save_binary:
            print("Binary file generation is disabled")

    # ---------------- RUN SIM ----------------
    print("\n" + "-"*50)
    print("RUNNING SIMULATION")
    print("-"*50)

    sim_logger.info("Initializing simulation...")
    scSim.InitializeSimulation()

    print(f"Setting stop time to {simulationTime} ns...")
    scSim.ConfigureStopTime(simulationTime)

    print(f"Starting simulation for {config.simulation_time} minutes...")
    start_time = datetime.now()
    scSim.ExecuteSimulation()
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()

    print("Simulation completed!")
    print(f"Wall-clock time: {duration:.2f} s | Simulated: {config.simulation_time:.1f} min | Speed: {(config.simulation_time*60)/max(duration,1):.1f}x")

    current_sim_time = scSim.TotalSim.CurrentNanos
    actual_sim_minutes = current_sim_time / 1e9 / 60
    print(f"Actual simulation time: {actual_sim_minutes:.2f} minutes")
    if abs(actual_sim_minutes - config.simulation_time) > 0.1:
        print("WARNING: Simulation time mismatch!")

    # ---------------- PLOTS (condensed logging) ----------------
    figureList = {}
    if (config.show_plots or config.save_plots) and PLOTS_AVAILABLE:
        print("\n" + "-"*50)
        print("GENERATING PLOTS")
        print("-"*50)

        time_data = np.linspace(0, config.simulation_time, max(100, int(config.simulation_time * 2)))
        plot_counts = {'constellation': 0, 'cluster': 0, 'distance': 0, 'fault': 0}

        # Constellation
        try:
            constellation_plots = generate_constellation_overview_plots(sc_objects, time_data, mu)
            figureList.update(constellation_plots)
            plot_counts['constellation'] = len(constellation_plots)
        except Exception as e:
            print(f"Constellation plots failed: {e}")

        # Cluster
        try:
            has_clusters = any(sc.get('cluster') for sc in config.spacecraft_list) or \
                           any(sc.get('type') == 'cluster_member' for sc in config.spacecraft_list)
            if has_clusters:
                cluster_data = _cluster_data_from_manager(cluster_manager, sc_objects) if cluster_manager else {}
                if not cluster_data:
                    cluster_data = _cluster_data_from_config(config, sc_objects)
                if cluster_data:
                    cluster_plots = generate_cluster_communication_plots(cluster_data, sc_objects, time_data)
                    figureList.update(cluster_plots)
                    plot_counts['cluster'] = len(cluster_plots)
        except Exception as e:
            print(f"Cluster plots failed: {e}")

        # Distances
        try:
            distance_plots = generate_inter_satellite_distance_plots(sc_objects, time_data, mu)
            figureList.update(distance_plots)
            plot_counts['distance'] = len(distance_plots)
        except Exception as e:
            print(f"Distance plots failed: {e}")

        # Faults (REAL sim per sat with faults)
        try:
            if FAULT_LOADER_AVAILABLE:
                fault_sats = [sc for sc in config.spacecraft_list if sc["fault"]["enabled"]]
                if fault_sats:
                    print(f"Generating fault plots for {len(fault_sats)} spacecraft...")
                for sc_config in fault_sats:
                    ft = sc_config["fault"]["type"]
                    if ft not in AVAILABLE_FAULTS:
                        continue
                    params = {
                        'fault_magnitude': sc_config["fault"]["magnitude"],
                        'fault_wheel': sc_config["fault"]["wheel"],
                        'fault_time_min': sc_config["fault"]["time"],
                        'simulation_time_min': config.simulation_time
                    }
                    real_fault_cfg = create_fault_config_for_real_simulation(ft, params)
                    fps = generate_fault_plots(
                        fault_type=ft,
                        fault_data=real_fault_cfg,
                        time_data=time_data,
                        fault_time_min=sc_config["fault"]["time"],
                        spacecraft_name=sc_config["name"]
                    )
                    if fps:
                        figureList.update(fps)
                        plot_counts['fault'] += len(fps)
        except Exception as e:
            print(f"Fault plots failed: {e}")

        print("\nPlot generation complete:")
        print(f"  Constellation: {plot_counts['constellation']}")
        print(f"  Cluster:       {plot_counts['cluster']}")
        print(f"  Distance:      {plot_counts['distance']}")
        print(f"  Fault:         {plot_counts['fault']}")
        print(f"  Total:         {sum(plot_counts.values())} plots")

        # Save Figures
        if config.save_plots and figureList:
            print(f"\nSaving plots to disk...")
            os.makedirs(PLOTTING_DIR, exist_ok=True)
            ts = datetime.now().strftime("%Y%m%d%H%M%S")

            saved_count = 0
            for name, fig in list(figureList.items()):
                try:
                    clean_name = name.replace(" ", "_").replace("/", "_")
                    plot_filename = f"{ts}_{clean_name}.png"
                    plot_path = os.path.join(PLOTTING_DIR, plot_filename)
                    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                    saved_count += 1
                    plt.close(fig)
                    del figureList[name]
                except Exception as e:
                    print(f"  ✗ Error saving plot {name}: {e}")

            print(f"Saved {saved_count} plots to {PLOTTING_DIR}")
        else:
            # Close figures if not saving, to free memory
            for _, fig in figureList.items():
                try:
                    plt.close(fig)
                except Exception:
                    pass
            figureList.clear()

    # Scenario object (as your GUI expects)
    class ConstellationScenario:
        def __init__(self, scSim, sc_objects, config, cluster_manager=None):
            self.TotalSim = scSim
            self.sc_objects = sc_objects
            self.targets = [t.to_dict() for t in config.targets]
            self.fault_type = getattr(config, 'fault_type', None)
            self.actual_sim_time = actual_sim_minutes
            self.cluster_manager = cluster_manager

    scenario = ConstellationScenario(scSim, sc_objects, config, cluster_manager)

    # ---------------- REAL ML Detection ----------------
    print("\n" + "="*60)
    print("FAULT DETECTION (LSTM OR ISOLATION FOREST)")
    print("="*60)

    ml_results = None
    if REAL_ML_AVAILABLE and fault_detection_tab:
        try:
            print("Running fault detection via router...")
            print(f"Selected model type: {fault_detection_tab.model_type_var.get()}")

            # Use the router to select between LSTM and Isolation Forest
            ml_results = run_fault_detection_on_scenario(
                scenario=scenario,
                fault_detection_tab=fault_detection_tab,
                scenario_config=config,
                output_dir=output_dir
            )

            if ml_results:
                print("FAULT DETECTION COMPLETED!")
                summary = ml_results.get('summary', {})
                print(f"   Spacecraft: {summary.get('total_spacecraft', 0)}")
                print(f"   Detections: {summary.get('total_detections', 0)}")
                if summary.get('detection_times'):
                    print(f"   First Detection: {min(summary.get('detection_times', [])):.1f} min")

        except Exception as e:
            print(f"Fault detection error: {e}")
            import traceback
            traceback.print_exc()
    elif not fault_detection_tab:
        print("WARNING: fault_detection_tab not provided - skipping fault detection")
    else:
        print("Fault detection router not available")

    # ---------------- Summary file ----------------
    try:
        summary_path = os.path.join(output_dir, "simulation_summary.txt")
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write("SPACECRAFT CONSTELLATION SIMULATION SUMMARY\n")
            f.write("="*50 + "\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("MODULE STATUS:\n")
            f.write(f"  Basilisk: {'Available' if BASILISK_AVAILABLE else 'Not available'}\n")
            f.write(f"  Cluster Integration: {'Available' if CLUSTER_INTEGRATION_AVAILABLE else 'Not available'}\n")
            f.write(f"  Plotting: {'Available' if PLOTS_AVAILABLE else 'Not available'}")
            if PLOTS_AVAILABLE:
                apis = []
                if PLOTS_API["new"]:
                    apis.append("new")
                if PLOTS_API["old"]:
                    apis.append("old")
                f.write(f" (APIs: {', '.join(apis)})")
            f.write("\n")
            f.write(f"  Fault Loader: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}\n")
            f.write(f"  Available Faults: {AVAILABLE_FAULTS}\n")
            if FAILED_FAULTS:
                f.write(f"  Failed Faults: {FAILED_FAULTS}\n")
            f.write("\n")

            f.write("SIMULATION CONFIGURATION:\n")
            f.write(f"  Duration: {config.simulation_time} minutes\n")
            f.write(f"  Time Step: 1 second\n")
            f.write(f"  Spacecraft: {len(config.spacecraft_list)}\n")
            rich_clusters = _normalize_rich_clusters(config, [])  # for count only
            if rich_clusters:
                f.write(f"  Clusters: {len(rich_clusters)}\n")
            f.write("\n")

            f.write("SPACECRAFT DETAILS:\n")
            for i, sc in enumerate(config.spacecraft_list):
                f.write(f"\nSpacecraft {i+1}: {sc['name']}\n")
                f.write(f"  Type: {sc.get('type','individual')}\n")
                if sc.get('cluster'):
                    f.write(f"  Cluster: {sc['cluster']}\n")
                    f.write(f"  Role: {sc.get('role','')}\n")
                f.write(f"  Orbit:\n")
                f.write(f"    - Semi-major axis: {sc['orbit']['a']} km\n")
                f.write(f"    - Altitude: {sc['orbit']['a'] - 6371:.1f} km\n")
                f.write(f"    - Eccentricity: {sc['orbit']['e']}\n")
                f.write(f"    - Inclination: {sc['orbit']['i']}°\n")
                f.write(f"    - RAAN: {sc['orbit'].get('Omega', 'NA')}°\n")

                if sc["fault"]["enabled"]:
                    f.write("  Fault: ENABLED\n")
                    f.write(f"    - Type: {sc['fault']['type']}\n")
                    f.write(f"    - Magnitude: {sc['fault']['magnitude']}\n")
                    f.write(f"    - Wheel: {sc['fault']['wheel']}\n")
                    f.write(f"    - Injection time: {sc['fault']['time']} minutes\n")
                else:
                    f.write("  Fault: DISABLED\n")

            f.write("\nOUTPUT:\n")
            f.write(f"  Results Directory: {output_dir}\n")
            f.write(f"  Plots Generated: {len(figureList) if figureList else 0}\n")

            if config.save_binary:
                vizard_paths = [
                    os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
                    os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
                ]

                for vp in vizard_paths:
                    if os.path.exists(vp):
                        size_mb = os.path.getsize(vp)/(1024*1024)
                        f.write(f"  Binary File: {vp} ({size_mb:.2f} MB)\n")
                        break

        print(f"Summary saved: {summary_path}")

    except Exception as e:
        print(f"Could not save simulation summary: {e}")

    # ---------------- Binary verification ----------------
    if config.save_binary:
        print("\n" + "-"*50)
        print("BINARY FILE VERIFICATION")
        print("-"*50)

        vizard_paths = [
            os.path.join(VIZ_DIR, f"{config.binary_filename}_UnityViz.bin"),
            os.path.join(VIZ_DIR, "_VizFiles", f"{config.binary_filename}_UnityViz.bin")
        ]

        binary_found = False
        for vp in vizard_paths:
            if os.path.exists(vp):
                binary_found = True
                size_mb = os.path.getsize(vp)/(1024*1024)
                print(f"Binary file created: {os.path.basename(vp)}")
                print(f"Location: {vp}")
                print(f"Size: {size_mb:.2f} MB")
                print(f"Duration: {config.simulation_time} minutes")
                print(f"Spacecraft: {len(sc_objects)}")
                break

        if not binary_found:
            print("Binary file not found in expected locations")

    print("\n" + "="*60)
    print("SIMULATION COMPLETED SUCCESSFULLY")
    print("="*60)
    if rich_clusters:
        print(f"Clusters: {len(rich_clusters)}")
    print(f"Total Satellites: {len(sc_objects)}")
    print(f"Simulation Time: {config.simulation_time} minutes")
    print(f"Real Fault Simulations: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}")

    return scenario, viz, figureList, output_dir



# ============= MODULE TEST =============
if __name__ == "__main__":
    print("\nTesting spacecraft_simulation module...")
    check_fault_modules(verbose=True)
    check_plots_module(verbose=True)

    print("\nModule status:")
    print(f"  Basilisk: {'Available' if BASILISK_AVAILABLE else 'Not available'}")
    print(f"  Cluster Integration: {'Available' if CLUSTER_INTEGRATION_AVAILABLE else 'Not available'}")
    print(f"  Plotting: {'Available' if PLOTS_AVAILABLE else 'Not available'}")
    if PLOTS_AVAILABLE:
        apis = []
        if PLOTS_API["new"]:
            apis.append("new")
        if PLOTS_API["old"]:
            apis.append("old")
        print(f"    APIs: {', '.join(apis)}")
    print(f"  Fault Loader: {'Available' if FAULT_LOADER_AVAILABLE else 'Not available'}")

    if BASILISK_AVAILABLE:
        config = SimulationConfig()
        config.simulation_time = 5.0
        config.spacecraft_list = [
            {
                "name": "TestSat1",
                "type": "cluster_member",
                "cluster": "ALPHA",
                "role": "leader",
                "orbit": {"a": 6971, "e": 0.01, "i": 55.0, "Omega": 0.0, "omega": 0.0, "f": 0.0},
                "orbit_name": "Test Orbit",
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 2.0,
                    "enabled": True,
                    "periodic": {"enabled": False, "interval": 360, "magnitude": 0.1, "wheel": 3}
                },
                "camera": {"position": [0.0, 0.0, 15.0], "fov": 80.0, "enabled": True},
                "communication": {"range": 2000.0, "fov": 30.0, "aHat_B": [0.0, 0.0, -1.0]},
                "targets": []
            },
            {
                "name": "TestSat2",
                "type": "cluster_member",
                "cluster": "ALPHA",
                "role": "child",
                "orbit": {"a": 6971, "e": 0.01, "i": 55.0, "Omega": 0.0, "omega": 0.0, "f": 10.0},
                "orbit_name": "Test Orbit",
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 2,
                    "time": 3.0,
                    "enabled": False,
                    "periodic": {"enabled": False, "interval": 360, "magnitude": 0.1, "wheel": 3}
                },
                "camera": {"position": [0.0, 0.0, 15.0], "fov": 80.0, "enabled": False},
                "communication": {"range": 2000.0, "fov": 30.0, "aHat_B": [0.0, 0.0, -1.0]},
                "targets": []
            }
        ]
        # Lightweight cluster example (GUI path)
        config.clusters = {"ALPHA": {"leader": 0, "children": [1]}}

        try:
            config.validate()
            print("Configuration validation passed")
            print("\nRunning test simulation...")
            scenario, viz, figureList, output_dir = run_custom_simulation(config)
            if scenario:
                print("\nTest simulation completed successfully!")
                print(f"Results saved to: {output_dir}")
            else:
                print("\nTest simulation failed")
        except Exception as e:
            print(f"Test failed: {e}")
            traceback.print_exc()
    else:
        print("\nCannot run test - Basilisk not available")

    print("\nspacecraft_simulation.py module test complete!")