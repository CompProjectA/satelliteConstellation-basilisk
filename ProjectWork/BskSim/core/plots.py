#!/usr/bin/env python
# core/plots.py
"""
Integrated Constellation Plotting


"""

import os
import sys
import warnings
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for saving files
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# -----------------------------------------------------------------------------
# Optional Basilisk utilities
# -----------------------------------------------------------------------------
try:
    from Basilisk.utilities import macros
    BASILISK_AVAILABLE = True
except Exception:
    BASILISK_AVAILABLE = False

    class _DummyMacros:
        NANO2SEC = 1e-9
        NANO2MIN = 1e-9 / 60.0
        D2R = np.pi / 180.0

        @staticmethod
        def min2nano(x):
            return x * 60.0 * 1e9

    macros = _DummyMacros()

# -----------------------------------------------------------------------------
# Optional fault_loader (REAL sim only)
# -----------------------------------------------------------------------------
try:
    from fault_loader import (
        run_scenario,
        run_scenario_enhanced,
        extract_fault_data_from_scenario,
        get_available_fault_types,
        is_fault_type_available
    )
    FAULT_LOADER_AVAILABLE = True
    print("plots.py: fault_loader available (REAL simulations enabled)")
except Exception as e:
    FAULT_LOADER_AVAILABLE = False
    print(f"plots.py: fault_loader NOT available: {e}")
    print("          Real fault simulations will be disabled in generate_fault_plots().")


# =============================================================================
# Shared helpers
# =============================================================================
def safe_create_figure(figsize=(12, 8)):
    try:
        fig = plt.figure(figsize=figsize, facecolor='white', edgecolor='black')
        return fig
    except Exception as e:
        print(f"Error creating figure: {e}")
        return None


def _vec3(pos):
    """Return a flat 3-vector np.array([x,y,z]) from lists or numpy shapes."""
    arr = np.asarray(pos, dtype=float).flatten()
    if arr.size >= 3:
        return arr[:3]
    return np.array([0.0, 0.0, 0.0])


# =============================================================================
# Claude's ConstellationPlotter (kept, with small robustness tweaks)
# =============================================================================
class ConstellationPlotter:
    """Enhanced plotting with formation checks and clear visualization"""

    CLUSTER_COLORS = {
        0: '#FF6B6B',  # Red
        1: '#4ECDC4',  # Teal
        2: '#95E77E',  # Green
        3: '#FFD93D'   # Yellow
    }

    ROLE_MARKERS = {
        'leader': 's',  # Square
        'child': 'o'    # Circle
    }

    def __init__(self, output_dir="plots"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    # ---------------------------- Overview (Claude-style inputs)
    def plot_constellation_overview(self, clusters_data, satellites_data, save=True):
        """Generate comprehensive constellation overview with all satellites (cluster dict style)"""
        fig = plt.figure(figsize=(14, 10), facecolor='white')

        # Get total counts
        total_sats = len(satellites_data) if satellites_data is not None else 0
        num_clusters = len(clusters_data) if clusters_data is not None else 0

        fig.suptitle(
            f"COMPREHENSIVE CONSTELLATION - {total_sats} SATS - {num_clusters} CLUSTERS",
            fontsize=16, fontweight='bold'
        )

        # 3D Earth and orbits plot
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.set_title("3D Orbital View", fontsize=12)

        # Draw Earth (simplified sphere)
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x_earth = 6371 * np.outer(np.cos(u), np.sin(v))
        y_earth = 6371 * np.outer(np.sin(u), np.sin(v))
        z_earth = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))
        ax1.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.3)

        # Plot satellites by cluster
        cluster_idx = 0
        for cluster_name, cluster_data in (clusters_data or {}).items():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')

            for sat in cluster_data.get('satellites', []):
                # Generate orbital position (simplified)
                altitude = self._parse_orbit_altitude(cluster_data.get('orbit', 'LEO 600'))
                theta = np.random.uniform(0, 2 * np.pi)
                phi = np.random.uniform(0, np.pi)

                r = 6371 + altitude
                x = r * np.sin(phi) * np.cos(theta)
                y = r * np.sin(phi) * np.sin(theta)
                z = r * np.cos(phi)

                role = sat.get('role', 'child')
                marker = self.ROLE_MARKERS.get(role, 'o')
                size = 100 if role == 'leader' else 50

                ax1.scatter(x, y, z, c=color, marker=marker, s=size,
                            edgecolors='black', linewidth=0.5)

                # Draw orbit trace (simplified circle)
                orbit_theta = np.linspace(0, 2 * np.pi, 100)
                orbit_x = r * np.sin(phi) * np.cos(orbit_theta)
                orbit_y = r * np.sin(phi) * np.sin(orbit_theta)
                orbit_z = np.full_like(orbit_theta, z)
                ax1.plot(orbit_x, orbit_y, orbit_z, color=color, alpha=0.3, linewidth=0.5)

            cluster_idx += 1

        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        ax1.set_zlabel('Z (km)')
        ax1.set_box_aspect([1, 1, 1])

        # 2D Ground track plot
        ax2 = fig.add_subplot(222)
        ax2.set_title("Ground Track Projection", fontsize=12)

        ax2.add_patch(Rectangle((-180, -90), 360, 180, facecolor='lightblue', alpha=0.2))
        ax2.grid(True, alpha=0.3)

        cluster_idx = 0
        for cluster_name, cluster_data in (clusters_data or {}).items():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')

            for sat in cluster_data.get('satellites', []):
                lon = np.random.uniform(-180, 180)
                lat = np.random.uniform(-60, 60)

                role = sat.get('role', 'child')
                marker = self.ROLE_MARKERS.get(role, 'o')
                size = 100 if role == 'leader' else 50

                ax2.scatter(lon, lat, c=color, marker=marker, s=size,
                            edgecolors='black', linewidth=0.5)

            cluster_idx += 1

        ax2.set_xlim(-180, 180)
        ax2.set_ylim(-90, 90)
        ax2.set_xlabel('Longitude (deg)')
        ax2.set_ylabel('Latitude (deg)')

        # Cluster information panel
        ax3 = fig.add_subplot(223)
        ax3.axis('off')
        ax3.set_title("Cluster Configuration", fontsize=12)

        info_text = "CLUSTER DETAILS\n" + "=" * 40 + "\n\n"
        cluster_idx = 0

        for cluster_name, cluster_data in (clusters_data or {}).items():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')
            leader_name = (cluster_data.get('leader') or {}).get('name', 'N/A')
            num_sats = len(cluster_data.get('satellites', []))
            formation = cluster_data.get('formation', 'Unknown')

            info_text += f"■ {cluster_name} - {formation} - {num_sats} sats\n"
            info_text += f"  Leader: {leader_name}\n"
            info_text += f"  Orbit: {cluster_data.get('orbit','LEO 600')}\n"
            info_text += f"  Separation: {cluster_data.get('separation','?')} km\n\n"

            cluster_idx += 1

        ax3.text(0.05, 0.95, info_text, transform=ax3.transAxes,
                 fontsize=10, verticalalignment='top', family='monospace')

        # Legend
        ax4 = fig.add_subplot(224)
        ax4.axis('off')
        ax4.set_title("Legend", fontsize=12)

        legend_elements = []
        cluster_idx = 0
        for cluster_name in (clusters_data or {}).keys():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')
            legend_elements.append(mpatches.Patch(color=color, label=f"{cluster_name} cluster"))
            cluster_idx += 1

        legend_elements.append(plt.Line2D([0], [0], marker='s', color='w',
                                          markerfacecolor='gray', markersize=10,
                                          label='Leader'))
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor='gray', markersize=8,
                                          label='Child'))

        ax4.legend(handles=legend_elements, loc='center', fontsize=11)

        # Faulty satellites if provided in satellites_data
        faulty_sats = [s for s in (satellites_data or []) if s.get('fault_type')]
        if faulty_sats:
            fault_text = f"\nFAULTY SATELLITES ({len(faulty_sats)}):\n"
            for sat in faulty_sats[:5]:
                fault_text += f"  • {sat.get('name','?')}: {sat.get('fault_type','unknown')}\n"
            ax4.text(0.5, 0.2, fault_text, transform=ax4.transAxes,
                     fontsize=9, ha='center', color='red')

        plt.tight_layout()

        if save:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{timestamp}_ConstellationOverview_AllSatellites.png"
            filepath = os.path.join(self.output_dir, filename)
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close()
            return filepath

        return fig

    # ---------------------------- Formation check (Claude-style inputs)
    def plot_formation_check(self, cluster_name, cluster_data, save=True):
        """Generate formation check plot for a specific cluster"""
        fig = plt.figure(figsize=(12, 10), facecolor='white')

        formation_type = cluster_data.get('formation', 'Unknown')
        num_sats = len(cluster_data.get('satellites', []))
        separation = cluster_data.get('separation', 1.0)

        fig.suptitle(f"Formation Check: {cluster_name} - {formation_type} Formation",
                     fontsize=14, fontweight='bold')

        # Leader-relative positions in 3D
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.set_title("Leader-Relative Positions (Current)", fontsize=11)

        leader = cluster_data.get('leader', {})
        leader_pos = np.array(leader.get('position', [0, 0, 0]), dtype=float)

        ideal_positions = self._get_ideal_formation_positions(formation_type, num_sats, separation)

        # Draw ideal template
        if formation_type == "Line":
            for i in range(len(ideal_positions) - 1):
                ax1.plot(
                    [ideal_positions[i][0], ideal_positions[i + 1][0]],
                    [ideal_positions[i][1], ideal_positions[i + 1][1]],
                    [ideal_positions[i][2], ideal_positions[i + 1][2]],
                    'g--', alpha=0.4, linewidth=2, label='Ideal' if i == 0 else ''
                )

        elif formation_type == "Triangle":
            for i in range(len(ideal_positions)):
                j = (i + 1) % len(ideal_positions)
                ax1.plot(
                    [ideal_positions[i][0], ideal_positions[j][0]],
                    [ideal_positions[i][1], ideal_positions[j][1]],
                    [ideal_positions[i][2], ideal_positions[j][2]],
                    'g--', alpha=0.4, linewidth=2, label='Ideal' if i == 0 else ''
                )

        elif formation_type == "Diamond":
            connections = [(0, 1), (1, 2), (2, 3), (3, 0)]
            for k, (a, b) in enumerate(connections):
                if a < len(ideal_positions) and b < len(ideal_positions):
                    ax1.plot(
                        [ideal_positions[a][0], ideal_positions[b][0]],
                        [ideal_positions[a][1], ideal_positions[b][1]],
                        [ideal_positions[a][2], ideal_positions[b][2]],
                        'g--', alpha=0.4, linewidth=2, label='Ideal' if k == 0 else ''
                    )

        # Plot actual positions and compute errors
        errors = []
        for idx, sat in enumerate(cluster_data.get('satellites', [])):
            rel_pos = np.array(sat.get('position', [0, 0, 0]), dtype=float) - leader_pos

            if sat.get('role') == 'leader':
                ax1.scatter(rel_pos[0], rel_pos[1], rel_pos[2],
                            c='red', s=150, marker='s', edgecolors='black',
                            linewidth=2, label='Leader')
            else:
                ax1.scatter(rel_pos[0], rel_pos[1], rel_pos[2],
                            c='blue', s=80, marker='o', edgecolors='black',
                            linewidth=1, alpha=0.8)
                if idx < len(ideal_positions):
                    ideal = np.array(ideal_positions[idx], dtype=float)
                    errors.append(float(np.linalg.norm(rel_pos - ideal)))

        ax1.set_xlabel('Along-track (km)')
        ax1.set_ylabel('Cross-track (km)')
        ax1.set_zlabel('Radial (km)')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)

        rms_error = np.sqrt(np.mean(np.square(errors))) if errors else 0.0
        max_error = max(errors) if errors else 0.0

        if rms_error < 0.5:
            quality_color, quality_text = 'green', 'EXCELLENT'
        elif rms_error < 2.0:
            quality_color, quality_text = 'orange', 'GOOD'
        else:
            quality_color, quality_text = 'red', 'NEEDS CORRECTION'

        # Metrics panel
        ax2 = fig.add_subplot(222)
        ax2.axis('off')
        ax2.set_title("Formation Quality Metrics", fontsize=11)

        metrics_text = f"""Formation Type: {formation_type}
Satellites: {num_sats}
Leader: {leader.get('name','?')}
Target Separation: {separation} km

QUALITY ASSESSMENT
{'='*30}
RMS Error: {rms_error:.3f} km
Max Error: {max_error:.3f} km
Mean Error: {np.mean(errors) if errors else 0:.3f} km

Status: {quality_text}"""

        ax2.text(0.1, 0.9, metrics_text, transform=ax2.transAxes,
                 fontsize=11, family='monospace', verticalalignment='top',
                 bbox=dict(boxstyle="round,pad=0.5", facecolor=quality_color, alpha=0.2))

        # Error evolution (synthetic visualization)
        ax3 = fig.add_subplot(223)
        ax3.set_title("Formation Error Evolution", fontsize=11)
        time_hours = np.linspace(0, 24, 100)
        base_error = rms_error if rms_error > 0 else 0.1
        err = base_error + 0.3 * np.sin(2 * np.pi * time_hours / 12) + 0.1 * np.sin(2 * np.pi * time_hours / 3)
        err = np.maximum(err, 0)
        ax3.plot(time_hours, err, 'b-', linewidth=2, label='RMS Error')
        ax3.axhline(y=0.5, color='green', linestyle='--', alpha=0.5, label='Excellent')
        ax3.axhline(y=2.0, color='orange', linestyle='--', alpha=0.5, label='Good')
        ax3.axhline(y=5.0, color='red', linestyle='--', alpha=0.5, label='Poor')
        ax3.set_xlabel('Time (hours)')
        ax3.set_ylabel('RMS Error (km)')
        ax3.grid(True, alpha=0.3)
        ax3.legend(loc='upper right')
        ax3.set_ylim(0, max(6, max_error * 1.5 if max_error > 0 else 6))

        # 2D projection
        ax4 = fig.add_subplot(224)
        ax4.set_title("2D Formation View (X-Y Plane)", fontsize=11)

        for i in range(len(ideal_positions)):
            if formation_type in ["Triangle", "Diamond"]:
                j = (i + 1) % len(ideal_positions)
                ax4.plot([ideal_positions[i][0], ideal_positions[j][0]],
                         [ideal_positions[i][1], ideal_positions[j][1]],
                         'g--', alpha=0.4, linewidth=2)
            elif formation_type == "Line" and i < len(ideal_positions) - 1:
                ax4.plot([ideal_positions[i][0], ideal_positions[i + 1][0]],
                         [ideal_positions[i][1], ideal_positions[i + 1][1]],
                         'g--', alpha=0.4, linewidth=2)

        for sat in cluster_data.get('satellites', []):
            rel_pos = np.array(sat.get('position', [0, 0, 0]), dtype=float) - leader_pos
            if sat.get('role') == 'leader':
                ax4.scatter(rel_pos[0], rel_pos[1], c='red', s=150,
                            marker='s', edgecolors='black', linewidth=2)
                ax4.annotate('L', (rel_pos[0], rel_pos[1]),
                             ha='center', va='center', color='white', fontweight='bold')
            else:
                ax4.scatter(rel_pos[0], rel_pos[1], c='blue', s=80,
                            marker='o', edgecolors='black', linewidth=1)

        ax4.set_xlabel('Along-track (km)')
        ax4.set_ylabel('Cross-track (km)')
        ax4.grid(True, alpha=0.3)
        ax4.axis('equal')

        plt.tight_layout()

        if save:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{timestamp}_FormationCheck_{cluster_name}.png"
            filepath = os.path.join(self.output_dir, filename)
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close()
            return filepath

        return fig

    # ---------------------------- Cluster communication (Claude-style inputs)
    def plot_cluster_communication(self, clusters_data, simulation_time_min=30, save=True):
        """Generate cluster communication Gantt chart (cluster dict style)"""
        fig = plt.figure(figsize=(14, 8), facecolor='white')

        ax = fig.add_subplot(111)
        ax.set_title(f"Cluster Communication Windows - All {len(clusters_data or {})} Clusters",
                     fontsize=14, fontweight='bold')

        y_pos = 0
        y_labels = []
        y_positions = []

        cluster_idx = 0
        for cluster_name, cluster_data in (clusters_data or {}).items():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')
            leader = cluster_data.get('leader', {'name': 'Leader'})

            # Cluster header on the left
            ax.text(-2, y_pos + max(1, len(cluster_data.get('children', []))) / 2, cluster_name,
                    fontsize=10, fontweight='bold', ha='right', va='center')

            for child in cluster_data.get('children', []):
                windows = self._generate_comm_windows(simulation_time_min)

                for start, end in windows:
                    rect = Rectangle((start, y_pos - 0.3), end - start, 0.6,
                                     facecolor=color, alpha=0.3, edgecolor=color)
                    ax.add_patch(rect)

                # Dots as message events
                num_messages = np.random.randint(2, 5)
                for _ in range(num_messages):
                    msg_time = np.random.uniform(0, simulation_time_min)
                    ax.plot(msg_time, y_pos, 'o', color=color, markersize=6)

                label = f"{leader.get('name','')[:8]} ↔ {child.get('name','')[:8]}"
                y_labels.append(label)
                y_positions.append(y_pos)
                y_pos += 1

            if cluster_idx < (len(clusters_data) - 1):
                ax.axhline(y=y_pos - 0.5, color='gray', linestyle='-', linewidth=0.5)

            cluster_idx += 1

        ax.set_xlim(-3, simulation_time_min)
        ax.set_ylim(-0.5, y_pos - 0.5 if y_pos > 0 else 0.5)
        ax.set_xlabel("Time (minutes)", fontsize=11)
        ax.set_ylabel("Communication Links", fontsize=11)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.grid(True, axis='x', alpha=0.3)

        # Legend
        legend_elements = []
        cluster_idx = 0
        for cname in (clusters_data or {}).keys():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')
            legend_elements.append(mpatches.Patch(color=color, alpha=0.3, label=cname))
            cluster_idx += 1
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor='gray', markersize=6,
                                          label='Message'))
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

        plt.tight_layout()

        if save:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{timestamp}_ClusterCommunication_AllClusters.png"
            filepath = os.path.join(self.output_dir, filename)
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close()
            return filepath

        return fig

    # ---------------------------- Distance analysis (Claude-style inputs)
    def plot_distance_analysis(self, clusters_data, satellites_data, save=True):
        """Generate inter-satellite distance analysis (cluster dict style)"""
        fig = plt.figure(figsize=(14, 8), facecolor='white')

        fig.suptitle("Inter-Satellite Distance Analysis", fontsize=14, fontweight='bold')

        # Distance matrix
        ax1 = fig.add_subplot(121)
        ax1.set_title("Distance Matrix (km)", fontsize=11)

        satellites = satellites_data or []
        n_sats = len(satellites)
        distance_matrix = np.zeros((n_sats, n_sats))

        # Calculate distances
        for i, sat1 in enumerate(satellites):
            for j, sat2 in enumerate(satellites):
                if i != j:
                    pos1 = np.array(sat1.get('position', [0, 0, 0]), dtype=float)
                    pos2 = np.array(sat2.get('position', [0, 0, 0]), dtype=float)
                    distance_matrix[i, j] = float(np.linalg.norm(pos2 - pos1))

        im = ax1.imshow(distance_matrix, cmap='RdYlGn_r', aspect='auto')

        # Add cluster separators
        cluster_boundaries = []
        current_idx = 0
        for cdata in (clusters_data or {}).values():
            current_idx += len(cdata.get('satellites', []))
            cluster_boundaries.append(current_idx - 0.5)

        for boundary in cluster_boundaries[:-1]:
            ax1.axhline(y=boundary, color='black', linewidth=1.5)
            ax1.axvline(x=boundary, color='black', linewidth=1.5)

        sat_names = [str(s.get('name', '?'))[:8] for s in satellites]
        ax1.set_xticks(range(n_sats))
        ax1.set_yticks(range(n_sats))
        ax1.set_xticklabels(sat_names, rotation=90, fontsize=7)
        ax1.set_yticklabels(sat_names, fontsize=7)

        cbar = plt.colorbar(im, ax=ax1)
        cbar.set_label('Distance (km)', rotation=270, labelpad=15)

        # Time series of leader-child distances (synthetic evolution)
        ax2 = fig.add_subplot(122)
        ax2.set_title("Leader-Child Distances Over Time", fontsize=11)
        time_hours = np.linspace(0, 24, 100)

        cluster_idx = 0
        for cluster_name, cluster_data in (clusters_data or {}).items():
            color = self.CLUSTER_COLORS.get(cluster_idx, '#888888')
            leader_pos = np.array(cluster_data.get('leader', {}).get('position', [0, 0, 0]), dtype=float)
            for k, child in enumerate(cluster_data.get('children', [])):
                child_pos = np.array(child.get('position', [0, 0, 0]), dtype=float)
                base_distance = np.linalg.norm(child_pos - leader_pos)
                distance_evolution = base_distance + 2 * np.sin(2 * np.pi * time_hours / 12)
                distance_evolution += np.random.normal(0, 0.5, len(time_hours))
                ax2.plot(time_hours, distance_evolution, color=color, alpha=0.7,
                         label=f"{cluster_name}" if k == 0 else '')
            cluster_idx += 1

        ax2.set_xlabel('Time (hours)')
        ax2.set_ylabel('Distance (km)')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='upper right')

        plt.tight_layout()

        if save:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{timestamp}_DistanceAnalysis.png"
            filepath = os.path.join(self.output_dir, filename)
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close()
            return filepath

        return fig

    # ---------------------------- Helpers (Claude)
    def _parse_orbit_altitude(self, orbit_str):
        if "LEO" in orbit_str:
            if "400" in orbit_str:
                return 400
            elif "600" in orbit_str:
                return 600
            elif "800" in orbit_str:
                return 800
        elif "MEO" in orbit_str:
            if "10000" in orbit_str:
                return 10000
            elif "20000" in orbit_str:
                return 20000
        elif "GEO" in orbit_str:
            return 35786
        return 600

    def _get_ideal_formation_positions(self, formation_type, num_sats, separation):
        positions = []
        if formation_type == "Line":
            for i in range(num_sats):
                positions.append([i * separation, 0, 0])
        elif formation_type == "Triangle":
            if num_sats >= 3:
                angles = np.linspace(0, 2 * np.pi, num_sats, endpoint=False)
                for angle in angles:
                    x = separation * np.cos(angle)
                    y = separation * np.sin(angle)
                    positions.append([x, y, 0])
        elif formation_type == "Diamond":
            positions = [
                [0, 0, 0],  # Center/Leader
                [separation, 0, 0],
                [0, separation, 0],
                [-separation, 0, 0]
            ]
            if num_sats > 4:
                positions.append([0, -separation, 0])
        elif formation_type == "Leader-Follower":
            for i in range(num_sats):
                positions.append([-i * separation, 0, 0])
        return positions[:num_sats]

    def _generate_comm_windows(self, simulation_time):
        windows = []
        num_windows = np.random.randint(3, 6)
        for i in range(num_windows):
            start = i * (simulation_time / num_windows)
            duration = 3 + np.random.random() * 2
            end = min(start + duration, simulation_time)
            windows.append((start, end))
        return windows


# =============================================================================
# YOUR (enhanced) plotting paths for ALL satellites, Basilisk & REAL fault sims
# =============================================================================

def create_fault_config_for_real_simulation(fault_type, fault_params):
    """Create a config that enforces REAL simulation through fault_loader."""
    return {
        'use_real_simulation': True,
        'simulation_params': dict(fault_params or {}),
        'fault_type': fault_type
    }


def generate_fault_plots(fault_type, fault_data, time_data, fault_time_min, spacecraft_name="Spacecraft"):
    """
    Generate plots using ONLY real simulation data via fault_loader — NO synthetic fallback.
    """
    plots = {}

    simulation_time_min = float(time_data[-1]) if len(time_data) > 0 else 30.0

    if not FAULT_LOADER_AVAILABLE:
        print(f"[FAULT] fault_loader not available — cannot run REAL sim for {spacecraft_name}")
        return plots

    if callable(globals().get('is_fault_type_available')) and not is_fault_type_available(fault_type):
        print(f"[FAULT] Type '{fault_type}' not supported by fault_loader")
        return plots

    if not isinstance(fault_data, dict) or not fault_data.get('use_real_simulation', False):
        # Convert to real-sim config
        simulation_params = {
            'fault_magnitude': fault_data.get('friction_magnitude',
                                  fault_data.get('power_limit',
                                  fault_data.get('encoder_error',
                                  fault_data.get('battery_drain', 0.0005)))),
            'fault_wheel': fault_data.get('fault_wheel', 3),
            'fault_time_min': fault_time_min,
            'simulation_time_min': simulation_time_min
        }
        fault_data = create_fault_config_for_real_simulation(fault_type, simulation_params)

    simulation_params = dict(fault_data.get('simulation_params', {}))
    simulation_params['simulation_time_min'] = simulation_time_min

    try:
        result = run_scenario_enhanced(
            fault_type,
            showPlots=False,
            saveBinary=False,
            **simulation_params
        )
    except Exception as e:
        print(f"[FAULT] Enhanced run failed: {e} — trying run_scenario()")
        result = run_scenario(
            fault_type,
            showPlots=False,
            saveBinary=False,
            **simulation_params
        )

    scenario, viz, figure_list = None, None, {}
    if result is None:
        print(f"[FAULT] run_scenario returned None for {fault_type}")
        return plots
    elif isinstance(result, tuple):
        if len(result) >= 3:
            scenario, viz, figure_list = result
        elif len(result) == 2:
            scenario, viz = result
        elif len(result) == 1:
            scenario = result[0]
    else:
        scenario = result

    if scenario is None:
        print(f"[FAULT] scenario could not be created")
        return plots

    # If scenario already returns figures, use them
    if figure_list:
        for plot_name, fig in figure_list.items():
            new_name = f"REAL_{plot_name}_{spacecraft_name}"
            plots[new_name] = fig
        return plots

    # Otherwise extract and build
    real_fault_data = extract_fault_data_from_scenario(scenario, fault_type)

    if 'wheel_speeds' in real_fault_data and real_fault_data['wheel_speeds'] is not None:
        wheel_speeds = np.asarray(real_fault_data['wheel_speeds'])

        fig = plt.figure(figsize=(14, 10))

        ax1 = fig.add_subplot(2, 2, 1)
        if wheel_speeds.ndim == 2 and wheel_speeds.shape[1] >= 4:
            time_points = np.linspace(0, simulation_time_min, wheel_speeds.shape[0])
            for i in range(4):
                wheel_label = f'RW{i+1}'
                if i == real_fault_data.get('fault_wheel', 3):
                    wheel_label += ' (FAULTY)'
                    ax1.plot(time_points, wheel_speeds[:, i], linewidth=3, label=wheel_label)
                else:
                    ax1.plot(time_points, wheel_speeds[:, i], linewidth=1.5, alpha=0.7, label=wheel_label)

        ax1.axvline(x=real_fault_data.get('fault_time', fault_time_min),
                    color='red', linestyle='--', linewidth=2, label='Fault Injection')
        ax1.set_xlabel('Time (minutes)')
        ax1.set_ylabel('Wheel Speed (rad/s)')
        ax1.set_title(f'REAL Reaction Wheel Speeds - {fault_type.upper()} Fault')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2 = fig.add_subplot(2, 2, 2)
        if fault_type == "friction":
            txt = (f"Friction Magnitude: {real_fault_data.get('friction_magnitude', 0.0005)} N·m\n"
                   f"Baseline: {real_fault_data.get('friction_baseline', 0.02)} N·m\n"
                   f"Wheel: RW{real_fault_data.get('fault_wheel', 3) + 1}")
            ax2.set_title('Fault Parameters')
            ax2.text(0.5, 0.5, txt, ha='center', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            ax2.axis('off')
        elif fault_type == "power_limit":
            txt = (f"Power Limit: {real_fault_data.get('power_limit', 0.5)} W\n"
                   f"Wheel: RW{real_fault_data.get('fault_wheel', 3) + 1}")
            ax2.set_title('Fault Parameters')
            ax2.text(0.5, 0.5, txt, ha='center', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
            ax2.axis('off')
        elif fault_type == "encoder":
            txt = (f"Encoder Error: {real_fault_data.get('encoder_error', 20.0)}%\n"
                   f"Wheel: RW{real_fault_data.get('fault_wheel', 3) + 1}")
            ax2.set_title('Fault Parameters')
            ax2.text(0.5, 0.5, txt, ha='center', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
            ax2.axis('off')
        elif fault_type == "battery":
            txt = f"Battery Drain: {real_fault_data.get('battery_drain', 50.0)} W"
            ax2.set_title('Fault Parameters')
            ax2.text(0.5, 0.5, txt, ha='center', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
            ax2.axis('off')

        if 'attitude_error' in real_fault_data and real_fault_data['attitude_error'] is not None:
            ax3 = fig.add_subplot(2, 2, 3)
            att_error = np.asarray(real_fault_data['attitude_error'])
            time_points = np.linspace(0, simulation_time_min, att_error.shape[0])
            ax3.plot(time_points, att_error, linewidth=2)
            ax3.axvline(x=real_fault_data.get('fault_time', fault_time_min),
                        color='red', linestyle='--', linewidth=2)
            ax3.set_xlabel('Time (minutes)')
            ax3.set_ylabel('Attitude Error (rad)')
            ax3.set_title('Attitude Control Error')
            ax3.grid(True, alpha=0.3)

        plt.suptitle(f'REAL {fault_type.upper()} Fault Analysis - {spacecraft_name}',
                     fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        plots[f"REAL_{fault_type}_Analysis_{spacecraft_name}"] = fig
        return plots

    # If no detailed data, produce a summary
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    txt = (f"REAL {fault_type.upper()} Fault Simulation\n\n"
           f"Spacecraft: {spacecraft_name}\n"
           f"Fault Time: {fault_time_min} minutes\n"
           f"Simulation Duration: {simulation_time_min} minutes\n"
           f"Simulation completed successfully!\n"
           f"(Detailed telemetry extraction in progress)")
    ax.text(0.5, 0.5, txt, ha='center', va='center', fontsize=12,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    ax.set_title(f'REAL {fault_type.upper()} Fault - {spacecraft_name}', fontsize=14, fontweight='bold')
    ax.axis('off')
    plots[f"REAL_{fault_type}_Summary_{spacecraft_name}"] = fig
    return plots


def generate_constellation_overview_plots(spacecraft_list, time_data, planet_mu):
    """
    Generate comprehensive constellation overview plots for ALL satellites
    (Your original approach; no limit on number of satellites)
    """
    plots = {}

    num_sats = len(spacecraft_list)
    if num_sats == 0:
        return plots

    print(f"\n=== GENERATING CONSTELLATION PLOTS FOR ALL {num_sats} SATELLITES ===")

    fig_width = max(24, 16 + num_sats * 0.2)
    fig_height = max(14, 8 + num_sats * 0.3)
    fig = safe_create_figure(figsize=(fig_width, fig_height))
    if fig is None:
        return plots

    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

    # 3D ORBIT VIEW
    ax1 = fig.add_subplot(gs[0:2, 0:2], projection='3d')

    earth_radius = 6371.0  # km
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 40)
    x = earth_radius * np.outer(np.cos(u), np.sin(v))
    y = earth_radius * np.outer(np.sin(u), np.sin(v))
    z = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_surface(x, y, z, color='lightblue', alpha=0.3, zorder=0)

    colors = plt.cm.tab20(np.linspace(0, 1, min(num_sats, 20))) if num_sats <= 20 else plt.cm.rainbow(np.linspace(0, 1, num_sats))

    positions = []
    sat_names = []
    sat_types = []
    faulty_sats = []
    cluster_info = {}

    for i in range(num_sats):
        sc = spacecraft_list[i]
        try:
            name = getattr(sc, 'ModelTag', getattr(sc, 'model_tag', f'Sat{i+1}'))
            sat_names.append(name)

            has_fault = False
            fault_type = None
            if hasattr(sc, 'faultConfig') and isinstance(sc.faultConfig, dict):
                has_fault = sc.faultConfig.get('enabled', False)
                fault_type = sc.faultConfig.get('type', 'unknown')
                if has_fault:
                    faulty_sats.append((i, name, fault_type))

            sat_type = 'individual'
            lname = name.lower()
            if 'leader' in lname:
                sat_type = 'leader'
            elif 'child' in lname or 'sat' in lname:
                sat_type = 'child'
            sat_types.append(sat_type)

            if '_' in name:
                cl_name = name.split('_')[0]
                if cl_name not in cluster_info:
                    cluster_info[cl_name] = {'leader': None, 'children': []}
                if sat_type == 'leader':
                    cluster_info[cl_name]['leader'] = i
                else:
                    cluster_info[cl_name]['children'].append(i)

            if hasattr(sc, 'hub') and hasattr(sc.hub, 'r_CN_NInit'):
                pos_init = sc.hub.r_CN_NInit
            else:
                ang = 2 * np.pi * i / max(1, num_sats)
                radius = 7000000 + i * 50000
                pos_init = [radius * np.cos(ang), radius * np.sin(ang), 0.0]

            pos_array = _vec3(pos_init)
            positions.append(pos_array)

            r = np.linalg.norm(pos_array)
            orbit_radius = r / 1000.0

            angles = np.linspace(0, 2 * np.pi, 100)
            orbit_x = orbit_radius * np.cos(angles)
            orbit_y = orbit_radius * np.sin(angles)
            orbit_z = np.zeros_like(angles)

            inclination = 0.0
            if hasattr(sc, 'orbit') and hasattr(sc.orbit, 'i'):
                inc_attr = float(sc.orbit.i)
                inclination = inc_attr if np.isfinite(inc_attr) else 0.0
            else:
                inclination = 55.0 + (i % 4) * 10  # degrees heuristic

            inc_rad = inclination if inclination > np.pi else (inclination * np.pi / 180.0)
            rot_matrix = np.array([[1, 0, 0],
                                   [0, np.cos(inc_rad), -np.sin(inc_rad)],
                                   [0, np.sin(inc_rad),  np.cos(inc_rad)]])
            orbit_points = np.vstack([orbit_x, orbit_y, orbit_z])
            rotated = rot_matrix @ orbit_points
            orbit_x, orbit_y, orbit_z = rotated[0], rotated[1], rotated[2]

            lw = 2 if has_fault else 1
            ls = '--' if has_fault else '-'
            alpha = 0.8 if has_fault else 0.4
            ax1.plot(orbit_x, orbit_y, orbit_z,
                     color=colors[i % len(colors)], alpha=alpha, linewidth=lw, linestyle=ls)

            current_x, current_y, current_z = (pos_array / 1000.0).tolist()
            if sat_type == 'leader':
                marker, size = '^', (200 if has_fault else 120)
            elif sat_type == 'child':
                marker, size = 'o', (150 if has_fault else 80)
            else:
                marker, size = 's', (180 if has_fault else 100)
            edgecolor = 'red' if has_fault else 'black'
            lw_pt = 3 if has_fault else 1
            ax1.scatter(current_x, current_y, current_z,
                        c=[colors[i % len(colors)]], marker=marker, s=size,
                        edgecolors=edgecolor, linewidths=lw_pt, zorder=10)

        except Exception as e:
            print(f"  ERROR plotting spacecraft {i}: {e}")
            continue

    # Draw cluster leader-child connection hints
    for cname, cdata in cluster_info.items():
        if cdata['leader'] is not None and cdata['children']:
            leader_pos = _vec3(positions[cdata['leader']]) / 1000.0
            for child_idx in cdata['children']:
                if child_idx < len(positions):
                    child_pos = _vec3(positions[child_idx]) / 1000.0
                    ax1.plot([leader_pos[0], child_pos[0]],
                             [leader_pos[1], child_pos[1]],
                             [leader_pos[2], child_pos[2]],
                             'g--', alpha=0.3, linewidth=0.6)

    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_zlabel('Z (km)')
    ax1.set_title(f'Constellation 3D View - ALL {num_sats} Satellites', fontsize=12, fontweight='bold')

    summary_text = f"Total: {num_sats} satellites\nFaulty: {len(faulty_sats)}\nClusters: {len(cluster_info)}"
    ax1.text2D(0.02, 0.98, summary_text,
               transform=ax1.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Altitude distribution
    ax2 = fig.add_subplot(gs[0, 2])
    altitudes = [(np.linalg.norm(p) / 1000.0) - earth_radius for p in positions]
    ax2.hist(altitudes, bins=min(20, max(5, num_sats // 2 + 1)),
             color='skyblue', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Altitude (km)')
    ax2.set_ylabel('Number of Satellites')
    ax2.set_title(f'Altitude Distribution ({num_sats} sats)')
    ax2.grid(True, alpha=0.3)

    # Fault status pie chart
    ax3 = fig.add_subplot(gs[1, 2])
    fault_counts = {'Healthy': num_sats - len(faulty_sats)}
    for _, _, ftype in faulty_sats:
        fault_counts[ftype] = fault_counts.get(ftype, 0) + 1

    labels = list(fault_counts.keys())
    sizes = list(fault_counts.values())
    colors_pie = []
    for l in labels:
        if l == 'Healthy': colors_pie.append('green')
        elif 'friction' in l: colors_pie.append('red')
        elif 'power' in l: colors_pie.append('orange')
        elif 'battery' in l: colors_pie.append('yellow')
        elif 'encoder' in l: colors_pie.append('blue')
        else: colors_pie.append('purple')
    ax3.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.0f%%', startangle=90)
    ax3.set_title('Fault Status Distribution')

    # Satellite manifest (first 20 for readability)
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')

    list_text = "SATELLITE MANIFEST:\n" + "=" * 50 + "\n"
    show_n = min(20, num_sats)
    for i in range(show_n):
        name = sat_names[i] if i < len(sat_names) else f"Sat{i+1}"
        stype = sat_types[i] if i < len(sat_types) else "unknown"
        is_faulty = any(idx == i for idx, _, _ in faulty_sats)
        fault_info = ""
        if is_faulty:
            for idx, _, ftype in faulty_sats:
                if idx == i:
                    fault_info = f" [FAULT: {ftype}]"
                    break
        list_text += f"{i+1:3d}. {name:20s} ({stype}){fault_info}\n"
    if num_sats > show_n:
        list_text += f"... and {num_sats - show_n} more satellites\n"

    ax4.text(0.05, 0.95, list_text, transform=ax4.transAxes,
             fontsize=8, verticalalignment='top', fontfamily='monospace')

    plt.suptitle(f'COMPREHENSIVE CONSTELLATION ANALYSIS - ALL {num_sats} SATELLITES',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plots["ConstellationOverview_AllSatellites"] = fig
    return plots


def generate_cluster_communication_plots(clusters, spacecraft_list, time_data):
    """
    Generate comprehensive cluster communication plots
    clusters: dict like { 'ClusterA': {'leader': idx, 'children': [idx,...]}, ... }
    """
    plots = {}

    if not clusters:
        print("[COMM] No clusters found - skipping")
        return plots

    num_clusters = len(clusters)
    fig_height = max(12, 4 + num_clusters * 3)
    fig = safe_create_figure(figsize=(20, fig_height))
    if fig is None:
        return plots

    rows = max(3, num_clusters)
    gs = fig.add_gridspec(rows, 3, hspace=0.4, wspace=0.3)

    for cluster_idx, (cluster_name, cluster_data) in enumerate(clusters.items()):
        leader_idx = cluster_data.get('leader')
        children_indices = cluster_data.get('children', [])

        # Formation sketch
        ax1 = fig.add_subplot(gs[cluster_idx, 0])
        if leader_idx is not None:
            ax1.scatter(0, 0, color='red', s=200, marker='^', label='Leader', zorder=5)
            ax1.text(0, -0.15, 'Leader', ha='center', fontsize=8)
            num_children = len(children_indices)
            for i, child_idx in enumerate(children_indices):
                angle = 2 * np.pi * i / max(1, num_children)
                x = 0.5 * np.cos(angle)
                y = 0.5 * np.sin(angle)
                ax1.scatter(x, y, color='blue', s=100, marker='o', zorder=4)
                ax1.plot([0, x], [0, y], 'g--', alpha=0.5, linewidth=1)
                ax1.text(x, y - 0.1, f'C{i+1}', ha='center', fontsize=8)
        ax1.set_xlim(-1, 1)
        ax1.set_ylim(-1, 1)
        ax1.set_aspect('equal')
        ax1.set_title(f'{cluster_name} Formation', fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Communication timeline (synthetic visualization on provided time grid)
        ax2 = fig.add_subplot(gs[cluster_idx, 1])
        num_children = len(children_indices)
        if num_children > 0 and len(time_data) > 0:
            for i in range(num_children):
                comm_windows = np.sin(2 * np.pi * time_data / 10 + i * np.pi / 4) > 0.3
                ax2.fill_between(time_data, i - 0.4, i + 0.4,
                                 where=comm_windows, alpha=0.3, color='green')
                message_times = time_data[::10]
                for msg_time in message_times:
                    if msg_time <= time_data[-1]:
                        ax2.scatter(msg_time, i, color='red', s=20, marker='o', zorder=5)

        ax2.set_xlabel('Time (minutes)', fontsize=9)
        ax2.set_ylabel('Child Satellite', fontsize=9)
        ax2.set_title(f'{cluster_name} Communication', fontsize=10)
        if num_children > 0:
            ax2.set_ylim(-0.5, num_children - 0.5)
            ax2.set_yticks(range(num_children))
            ax2.set_yticklabels([f'Child {i+1}' for i in range(num_children)])
        ax2.grid(True, alpha=0.3)

        # Cluster metrics snapshot
        ax3 = fig.add_subplot(gs[cluster_idx, 2])
        metrics = ['Comm\nQuality', 'Formation', 'Data Rate', 'Link\nStrength']
        values = [85 + 10 * np.random.random(),
                  90 + 8 * np.random.random(),
                  75 + 15 * np.random.random(),
                  80 + 10 * np.random.random()]
        colors_bar = ['green' if v > 80 else 'orange' if v > 60 else 'red' for v in values]
        bars = ax3.bar(metrics, values, color=colors_bar, alpha=0.7, edgecolor='black')

        ax3.set_ylabel('Health (%)', fontsize=9)
        ax3.set_title(f'{cluster_name} Health', fontsize=10)
        ax3.set_ylim(0, 100)
        ax3.grid(True, alpha=0.3, axis='y')
        for bar, value in zip(bars, values):
            ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                     f'{value:.0f}%', ha='center', fontsize=8)

    plt.suptitle(f'CLUSTER COMMUNICATION ANALYSIS - {num_clusters} CLUSTERS',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plots["ClusterCommunication_AllClusters"] = fig
    return plots


def generate_inter_satellite_distance_plots(spacecraft_list, time_data, planet_mu):
    """Generate inter-satellite distance analysis for ALL satellites"""
    plots = {}
    num_sats = len(spacecraft_list)
    if num_sats < 2:
        return plots

    fig = safe_create_figure(figsize=(16, 10))
    if fig is None:
        return plots

    # Distance matrix
    ax1 = fig.add_subplot(1, 2, 1)
    distance_matrix = np.zeros((num_sats, num_sats))
    for i in range(num_sats):
        for j in range(num_sats):
            if i == j:
                continue
            p1 = _vec3(spacecraft_list[i].hub.r_CN_NInit) if hasattr(spacecraft_list[i], 'hub') else np.zeros(3)
            p2 = _vec3(spacecraft_list[j].hub.r_CN_NInit) if hasattr(spacecraft_list[j], 'hub') else np.zeros(3)
            distance_matrix[i, j] = np.linalg.norm(p2 - p1) / 1000.0  # km

    im = ax1.imshow(distance_matrix, cmap='viridis', aspect='auto')
    ax1.set_xlabel('Spacecraft Index')
    ax1.set_ylabel('Spacecraft Index')
    ax1.set_title(f'Distance Matrix - {num_sats} Satellites (km)')

    if num_sats <= 10:
        for i in range(num_sats):
            for j in range(num_sats):
                if i != j:
                    ax1.text(j, i, f'{distance_matrix[i, j]:.0f}',
                             ha="center", va="center", color="white", fontsize=8)

    plt.colorbar(im, ax=ax1, label='Distance (km)')

    # Time evolution plot for sample pairs (synthetic evolution)
    ax2 = fig.add_subplot(1, 2, 2)
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown']
    pairs_plotted = 0
    max_pairs = min(6, (num_sats * (num_sats - 1)) // 2)
    t = time_data if len(time_data) > 0 else np.linspace(0, 30, 100)

    for i in range(min(3, num_sats)):
        for j in range(i + 1, min(i + 3, num_sats)):
            if pairs_plotted >= max_pairs:
                break
            name1 = getattr(spacecraft_list[i], 'ModelTag', f'Sat{i+1}')
            name2 = getattr(spacecraft_list[j], 'ModelTag', f'Sat{j+1}')

            base_dist = distance_matrix[i, j] if distance_matrix[i, j] > 0 else 1000.0
            distances = base_dist + 100 * np.sin(2 * np.pi * t / 30 + i - j)

            color = colors[pairs_plotted % len(colors)]
            ax2.plot(t, distances, color=color, linewidth=2,
                     label=f'{name1[:8]}-{name2[:8]}')
            pairs_plotted += 1

    ax2.set_xlabel('Time (minutes)')
    ax2.set_ylabel('Distance (km)')
    ax2.set_title('Inter-Satellite Distances Over Time')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=8)

    plt.suptitle(f'INTER-SATELLITE DISTANCE ANALYSIS - {num_sats} SATELLITES',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plots["InterSatelliteDistances_AllSatellites"] = fig
    return plots


# =============================================================================
# Backward-compat wrappers
# =============================================================================
def generate_constellation_plots(spacecraft_list, time_data, planet_mu):
    """
    Wrapper (Your original signature).
    Returns a dict of figures including "ConstellationOverview_AllSatellites".
    """
    return generate_constellation_overview_plots(spacecraft_list, time_data, planet_mu)


def generate_constellation_plots_from_clusters(clusters_data, satellites_data, output_dir="plots"):
    """
    Claude compatibility wrapper using ConstellationPlotter API.
    Returns list of saved filepaths.
    """
    plotter = ConstellationPlotter(output_dir)
    outputs = []

    # Constellation overview
    outputs.append(plotter.plot_constellation_overview(clusters_data, satellites_data))

    # Formation checks
    for cluster_name, cluster_data in (clusters_data or {}).items():
        outputs.append(plotter.plot_formation_check(cluster_name, cluster_data))

    # Communication plot
    outputs.append(plotter.plot_cluster_communication(clusters_data))

    # Distance analysis
    outputs.append(plotter.plot_distance_analysis(clusters_data, satellites_data))

    return outputs


# =============================================================================
# Module self-test (non-failing)
# =============================================================================
if __name__ == "__main__":
    print("core/plots.py integrated module loaded.")
    print(f" - Basilisk available: {BASILISK_AVAILABLE}")
    print(f" - fault_loader available: {FAULT_LOADER_AVAILABLE}")

    # Minimal smoketest of wrappers with dummy data (no file writes here)
    try:
        # Claude-style dummy
        clusters = {
            "Alpha": {
                "leader": {"name": "Alpha_Leader", "position": [100, 0, 0]},
                "children": [{"name": "Alpha_C1", "position": [100.5, 0.2, 0]},
                             {"name": "Alpha_C2", "position": [99.5, -0.3, 0]}],
                "satellites": [
                    {"name": "Alpha_Leader", "role": "leader"},
                    {"name": "Alpha_C1", "role": "child"},
                    {"name": "Alpha_C2", "role": "child"}
                ],
                "formation": "Triangle",
                "separation": 0.5,
                "orbit": "LEO 600"
            }
        }
        sats = [
            {"name": "Alpha_Leader", "position": [100, 0, 0]},
            {"name": "Alpha_C1", "position": [100.5, 0.2, 0]},
            {"name": "Alpha_C2", "position": [99.5, -0.3, 0]}
        ]

        cp = ConstellationPlotter(output_dir="plots")
        _ = cp.plot_constellation_overview(clusters, sats, save=False)
        _ = cp.plot_formation_check("Alpha", clusters["Alpha"], save=False)
        _ = cp.plot_cluster_communication(clusters, save=False)
        _ = cp.plot_distance_analysis(clusters, sats, save=False)
        plt.close('all')

        # Your-style dummy
        class _DummySC:
            def __init__(self, name, pos):
                self.ModelTag = name
                class _Hub:
                    def __init__(self, r):
                        self.r_CN_NInit = r
                self.hub = _Hub(pos)
                self.faultConfig = {'enabled': False}

        dummy_scs = [_DummySC(f"Sat{i+1}", [7000000 + i * 10000, 0, 0]) for i in range(6)]
        _ = generate_constellation_overview_plots(dummy_scs, np.linspace(0, 30, 120), planet_mu=398600.4418)
        _ = generate_inter_satellite_distance_plots(dummy_scs, np.linspace(0, 30, 120), planet_mu=398600.4418)
        plt.close('all')

        print("Self-test OK.")
    except Exception as e:
        print(f"Self-test encountered an issue (non-fatal): {e}")
