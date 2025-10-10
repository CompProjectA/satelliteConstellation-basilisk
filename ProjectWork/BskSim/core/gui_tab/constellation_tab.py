#!/usr/bin/env python
"""
core/gui_tab/constellation_tab.py


"""

import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use("Agg")  # GUI is embedded; figures saved offscreen
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed by mpl for 3D proj)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Optional helper formations (kept from your code)
try:
    from formation_geometry import FormationGeometry, CartesianFormations  # type: ignore
except Exception:  # Safe fallback if module not present
    FormationGeometry = None
    class _Dummy:
        def __getattr__(self, *_args, **_kwargs):
            raise AttributeError("formation_geometry not available")
    CartesianFormations = _Dummy()


# --- BaseTab fallback (so this file can be imported directly) -----------------
try:
    from .base_tab import BaseTab  # type: ignore
except Exception:  # pragma: no cover
    class BaseTab:
        def __init__(self, parent_app, parent_frame):
            self.parent_app = parent_app
            self.parent_frame = parent_frame

        def add_log(self, message):
            if hasattr(self.parent_app, "add_log"):
                self.parent_app.add_log(message)


class ConstellationTab(BaseTab):
    """Constellation Management (clusters, orbits, individuals, comms)"""

    # Project constraints (Claude)
    ALLOWED_FORMATIONS = ["Leader-Follower", "Line", "Triangle", "Diamond"]
    MAX_CLUSTERS = 4

    # Stable palette (shared with plots/Vizard via config.cluster_colors)
    PALETTE = ["#2ecc71", "#e74c3c", "#9b59b6", "#f39c12"]

    def __init__(self, parent_app, parent_frame):
        super().__init__(parent_app, parent_frame)

        # References to app-wide data
        self.satellites = parent_app.satellites
        self.current_satellite_index = parent_app.current_satellite_index

        # Cluster & orbit state
        self.clusters = []                # list of dicts: {name, leader, children, satellites, formation, ...}
        self.cluster_configurations = []  # reserved (compat)
        self.orbit_configurations = [
            {
                "name": "LEO 600km",
                "altitude": 600,
                "inclination": 55.0,
                "raan": 0.0,
                "satellites": [],
                "description": "Low Earth Orbit for cluster deployment",
            },
            {
                "name": "LEO 700km",
                "altitude": 700,
                "inclination": 56.0,
                "raan": 0.0,
                "satellites": [],
                "description": "Standard LEO for constellations",
            },
        ]

        # Communication defaults
        self.comm_range = tk.DoubleVar(value=2000.0)  # km
        self.comm_fov = tk.DoubleVar(value=30.0)      # deg

        # Formation defaults
        self.formation_separation = tk.DoubleVar(value=10.0)  # km

        # Build UI
        self.create_tab_ui()

    # --------------------------------------------------------------------- UI
    def create_tab_ui(self):
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.constellation_notebook = ttk.Notebook(main_container)
        self.constellation_notebook.pack(fill=tk.BOTH, expand=True)

        cluster_frame = ttk.Frame(self.constellation_notebook)
        orbit_frame = ttk.Frame(self.constellation_notebook)
        individual_frame = ttk.Frame(self.constellation_notebook)
        comm_frame = ttk.Frame(self.constellation_notebook)

        self.constellation_notebook.add(cluster_frame, text="Cluster Management")
        self.constellation_notebook.add(orbit_frame, text="Orbit Management")
        self.constellation_notebook.add(individual_frame, text="Individual Satellites")
        self.constellation_notebook.add(comm_frame, text="Communication")

        self._create_cluster_tab(cluster_frame)
        self._create_orbit_management_tab(orbit_frame)
        self._create_individual_tab(individual_frame)
        self._create_comm_tab(comm_frame)

    # ------------------------------------------------------------ Cluster Tab
    def _create_cluster_tab(self, parent):
        # Create
        create_frame = ttk.LabelFrame(parent, text="Create New Cluster", padding=10)
        create_frame.pack(fill=tk.X, padx=10, pady=10)

        basic = ttk.Frame(create_frame)
        basic.pack(fill=tk.X, pady=5)

        ttk.Label(basic, text="Cluster Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.cluster_name_var = tk.StringVar(value="")
        ttk.Entry(basic, textvariable=self.cluster_name_var, width=20).grid(row=0, column=1, padx=5)

        ttk.Label(basic, text="Number of Satellites:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.sats_per_cluster_var = tk.IntVar(value=4)
        sats_combo = ttk.Combobox(
            basic, textvariable=self.sats_per_cluster_var,
            values=[2, 3, 4, 5, 6], width=10, state="readonly"
        )
        sats_combo.grid(row=0, column=3, padx=5)
        sats_combo.current(2)  # default 4

        formation_row = ttk.Frame(create_frame)
        formation_row.pack(fill=tk.X, pady=5)

        ttk.Label(formation_row, text="Formation Type:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.formation_var = tk.StringVar(value=self.ALLOWED_FORMATIONS[0])
        formation_combo = ttk.Combobox(
            formation_row, textvariable=self.formation_var,
            values=self.ALLOWED_FORMATIONS, width=18, state="readonly"
        )
        formation_combo.grid(row=0, column=1, padx=5)

        ttk.Label(formation_row, text="Separation (km):").grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Entry(formation_row, textvariable=self.formation_separation, width=10).grid(row=0, column=3, padx=5)

        self.cartesian_mode = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            formation_row, text="Cartesian (no-orbit) mode", variable=self.cartesian_mode
        ).grid(row=0, column=4, padx=10)

        # Orbit selection (single orbit per cluster)
        orbit_frame = ttk.LabelFrame(create_frame, text="Cluster Orbit", padding=5)
        orbit_frame.pack(fill=tk.X, pady=5)

        ttk.Label(orbit_frame, text="Select Orbit:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.primary_orbit_var = tk.StringVar()
        orbit_combo = ttk.Combobox(
            orbit_frame, textvariable=self.primary_orbit_var,
            values=[o["name"] for o in self.orbit_configurations],
            width=20, state="readonly",
        )
        orbit_combo.grid(row=0, column=1, padx=5)
        if self.orbit_configurations:
            self.primary_orbit_var.set(self.orbit_configurations[0]["name"])
            orbit_combo.current(0)

        ttk.Label(
            orbit_frame,
            text="All satellites share this orbit; small offsets create the formation.",
            foreground="gray",
        ).grid(row=1, column=0, columnspan=2, pady=3, sticky=tk.W)

        # Buttons
        btns = ttk.Frame(create_frame)
        btns.pack(fill=tk.X, pady=10)
        ttk.Button(btns, text="Create Cluster", command=self.create_cluster, style="Run.TButton").pack(side=tk.LEFT, padx=5)
        ttk.Button(btns, text="Clear Form", command=self.clear_cluster_form).pack(side=tk.LEFT, padx=5)

        # List
        list_frame = ttk.LabelFrame(parent, text="Existing Clusters", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = ("Cluster", "Leader", "Children", "Formation", "Status")
        self.cluster_tree = ttk.Treeview(list_frame, columns=columns, show="tree headings", height=10)
        self.cluster_tree.heading("#0", text="")
        for col in columns:
            self.cluster_tree.heading(col, text=col)
        self.cluster_tree.column("#0", width=30)
        self.cluster_tree.column("Cluster", width=130)
        self.cluster_tree.column("Leader", width=150)
        self.cluster_tree.column("Children", width=90)
        self.cluster_tree.column("Formation", width=120)
        self.cluster_tree.column("Status", width=100)
        self.cluster_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        sb = ttk.Scrollbar(list_frame, orient="vertical", command=self.cluster_tree.yview)
        self.cluster_tree.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        # Actions
        actions = ttk.Frame(list_frame)
        actions.pack(fill=tk.X, pady=5)
        ttk.Button(actions, text="View Details", command=self.view_cluster_details).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Modify Cluster", command=self.modify_cluster).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Delete Cluster", command=self.delete_cluster).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Test Communication", command=self.test_cluster_communication).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="View Formation", command=self.view_formation).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Export Config", command=self.export_config).pack(side=tk.LEFT, padx=5)

        self.update_cluster_tree()

    def _update_orbit_assignments(self, *_args):
        """No-op (compat) — clusters use a single selected orbit."""
        return

    # ------------------------------------------------------------ Create/Modify
    def create_cluster(self):
        # Cap at 4 clusters (Claude requirement)
        if len(self.clusters) >= self.MAX_CLUSTERS:
            messagebox.showwarning("Cluster Limit", f"Only {self.MAX_CLUSTERS} clusters are allowed.")
            return

        name = self.cluster_name_var.get().strip()
        if not name:
            messagebox.showerror("Input Error", "Please enter a cluster name")
            return
        if any(c["name"] == name for c in self.clusters):
            messagebox.showerror("Duplicate Name", f"Cluster '{name}' already exists")
            return

        n = int(self.sats_per_cluster_var.get())
        formation = self.formation_var.get()
        separation = float(self.formation_separation.get())

        if not self.orbit_configurations:
            messagebox.showwarning("No Orbits", "Create or import at least one orbit first.")
            return

        # Base orbit for the entire cluster
        primary_orbit = next((o for o in self.orbit_configurations if o["name"] == self.primary_orbit_var.get()),
                             self.orbit_configurations[0])

        cluster_base_angle = (len(self.clusters) * 90.0) % 360.0  # spread clusters around the Earth

        base_orbit = {
            "a": primary_orbit["altitude"] + 6371,
            "e": 0.001,
            "i": primary_orbit["inclination"],
            "Omega": primary_orbit.get("raan", 0.0),
            "omega": 0.0,
            "f": cluster_base_angle,
        }

        # Small orbital offsets to shape formation
        offsets = self._get_formation_offsets(n, formation, separation)

        # Optional cartesian offsets (visual spacing used by plots/Vizard)
        cart = None
        try:
            # Map formation to helper shapes
            f = formation.lower()
            if "leader-follower" in f:
                cart = CartesianFormations.train(n, separation)
            elif "line" in f:
                cart = CartesianFormations.column(n, separation)
            elif "triangle" in f or "diamond" in f:
                cart = CartesianFormations.box(n, separation)
            else:
                cart = CartesianFormations.box(n, separation)
        except Exception:
            cart = None  # fall back to storing just numeric offsets we compute below

        # Build cluster record
        cluster = {
            "name": name,
            "type": "cluster",
            "formation": formation,
            "leader": None,
            "children": [],
            "satellites": [],
            "separation": separation,
            "cartesian_mode": bool(self.cartesian_mode.get()),
            "color_idx": len(self.clusters) % len(self.PALETTE),
            "color_hex": self.PALETTE[len(self.clusters) % len(self.PALETTE)],
        }

        for i in range(n):
            is_leader = (i == 0)
            sat_name = f"{name}_{'Leader' if is_leader else f'Sat{i+1}'}"

            sat_orbit = {
                "a": base_orbit["a"],
                "e": base_orbit["e"],
                "i": base_orbit["i"] + offsets[i]["di"],
                "Omega": base_orbit["Omega"] + offsets[i]["dOmega"],
                "omega": base_orbit["omega"],
                "f": base_orbit["f"] + offsets[i]["df"],
            }

            # Communication pointing defaults
            aHat = [0.2, -0.4, 0.2] if is_leader else [0.0, 0.0, -1.0]

            pos_off = [0.0, 0.0, 0.0]
            if cart is not None:
                try:
                    pos_off = cart[i].offset_km  # provided by helper
                except Exception:
                    pos_off = [0.0, 0.0, 0.0]

            sat = {
                "name": sat_name,
                "type": "cluster_member",
                "cluster": name,
                "role": "leader" if is_leader else "child",
                "cartesian_mode": bool(self.cartesian_mode.get()),
                "position_offset_km": pos_off,
                "formation": formation,
                "orbit": sat_orbit,
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 15.0,
                    "enabled": False,
                },
                "camera": {
                    "position": [10.0, 10.0, 10.0],
                    "fov": 90.0,
                    "enabled": is_leader,
                },
                "communication": {
                    "range": float(self.comm_range.get()),
                    "fov": float(self.comm_fov.get()),
                    "aHat_B": aHat,
                },
                "targets": [],
                "orbit_name": primary_orbit["name"],
            }

            self.satellites.append(sat)
            cluster["satellites"].append(sat_name)
            if is_leader:
                cluster["leader"] = sat_name
            else:
                cluster["children"].append(sat_name)

            primary_orbit["satellites"].append(sat_name)

        # Book-keeping
        self.clusters.append(cluster)
        cluster["_members"] = [next(s for s in self.satellites if s["name"] == nm) for nm in cluster["satellites"]]

        # UI updates
        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        if hasattr(self.parent_app, "update_satellite_dropdowns"):
            self.parent_app.update_satellite_dropdowns()
        if hasattr(self.parent_app, "update_status_counts"):
            self.parent_app.update_status_counts()

        # Log
        self.add_log(f"Created cluster '{name}' with {n} satellites in {formation} formation")
        self.add_log(f"  Formation separation: {separation} km")
        self.add_log(f"  All satellites in orbit: {primary_orbit['name']}")
        self.add_log(f"  Leader: {cluster['leader']}")
        self.add_log(f"  Children: {', '.join(cluster['children'])}")

        # Warn if limit reached
        if len(self.clusters) >= self.MAX_CLUSTERS:
            messagebox.showinfo(
                "Cluster Limit Reached",
                f"Maximum {self.MAX_CLUSTERS} clusters created (project requirement)."
            )

    def _get_formation_offsets(self, num_sats, formation, separation_km):
        """
        Small orbital deltas producing tight formations.
        Returns: list of dicts {'df','di','dOmega'} (degrees)
        """
        offsets = []
        deg_per_km = 1.0 / 116.0  # approx along-track conversion
        sep_deg = separation_km * deg_per_km

        f = (formation or "").lower()

        if "line" in f:
            for i in range(num_sats):
                offsets.append({"df": i * sep_deg, "di": 0.0, "dOmega": 0.0})

        elif "triangle" in f:
            if num_sats >= 1:
                offsets.append({"df": 0.0, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 2:
                offsets.append({"df": -sep_deg, "di": +sep_deg * 0.005, "dOmega": 0.0})
            if num_sats >= 3:
                offsets.append({"df": -sep_deg, "di": -sep_deg * 0.005, "dOmega": 0.0})
            for i in range(3, num_sats):
                offsets.append({"df": -sep_deg * 0.5, "di": 0.0, "dOmega": 0.0})

        elif "diamond" in f:
            if num_sats >= 1:
                offsets.append({"df": 0.0, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 2:
                offsets.append({"df": +sep_deg, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 3:
                offsets.append({"df": -sep_deg, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 4:
                offsets.append({"df": 0.0, "di": +sep_deg * 0.01, "dOmega": 0.0})
            for i in range(4, num_sats):
                offsets.append({"df": (i - 3) * sep_deg * 0.5, "di": 0.0, "dOmega": 0.0})

        else:  # Leader-Follower (train)
            for i in range(num_sats):
                offsets.append({"df": -i * sep_deg, "di": 0.0, "dOmega": 0.0})

        return offsets

    # ------------------------------------------------------------ Cluster ops
    def view_cluster_details(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to view")
            return
        cluster_name = self.cluster_tree.item(sel[0])["values"][0]
        cluster = next((c for c in self.clusters if c["name"] == cluster_name), None)
        if not cluster:
            return

        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Cluster Details: {cluster_name}")
        w.geometry("600x520")
        frame = ttk.Frame(w, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        txt = tk.Text(frame, wrap=tk.WORD, font=("Consolas", 10))
        sb = ttk.Scrollbar(frame, command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        details = [
            f"CLUSTER: {cluster['name']}",
            "=" * 60,
            f"Formation: {cluster['formation']}",
            f"Leader: {cluster['leader']}",
            f"Children: {', '.join(cluster['children']) or '-'}",
            f"Total Satellites: {len(cluster['satellites'])}",
            f"Separation: {cluster.get('separation', 10.0)} km",
            "",
            "SATELLITES",
            "-" * 60,
        ]
        for nm in cluster["satellites"]:
            sat = next((s for s in self.satellites if s["name"] == nm), None)
            if not sat:
                continue
            details += [
                f"{sat['name']}:",
                f"  Role: {sat['role']}",
                f"  Orbit: {sat.get('orbit_name','')}",
                f"  Altitude: {sat['orbit']['a'] - 6371:.0f} km",
                f"  Inc: {sat['orbit']['i']:.1f}°, RAAN: {sat['orbit']['Omega']:.1f}°, TA: {sat['orbit']['f']:.1f}°",
                f"  Camera: {'Enabled' if sat['camera']['enabled'] else 'Disabled'}",
                f"  Fault: {'Enabled' if sat['fault']['enabled'] else 'Disabled'}",
                ""
            ]
        details += [
            "=" * 60, "COMMUNICATION",
            f"Range: {self.comm_range.get():.0f} km | FOV: {self.comm_fov.get():.1f}°",
        ]

        txt.insert(tk.END, "\n".join(details))
        txt.configure(state="disabled")
        ttk.Button(w, text="Close", command=w.destroy).pack(pady=8)

    def modify_cluster(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to modify")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        cluster = next((c for c in self.clusters if c["name"] == name), None)
        if not cluster:
            return

        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Modify Cluster: {name}")
        w.geometry("380x240")
        f = ttk.Frame(w, padding=20); f.pack(fill=tk.BOTH, expand=True)

        ttk.Label(f, text="Formation Type:").grid(row=0, column=0, sticky=tk.W, pady=5)
        form_var = tk.StringVar(value=cluster["formation"])
        ttk.Combobox(f, textvariable=form_var, values=self.ALLOWED_FORMATIONS, width=20, state="readonly").grid(row=0, column=1, pady=5)

        ttk.Label(f, text="Separation (km):").grid(row=1, column=0, sticky=tk.W, pady=5)
        sep_var = tk.DoubleVar(value=cluster.get("separation", 10.0))
        ttk.Entry(f, textvariable=sep_var, width=20).grid(row=1, column=1, pady=5)

        def apply():
            cluster["formation"] = form_var.get()
            cluster["separation"] = float(sep_var.get())
            self._update_cluster_formation(cluster)
            self.update_cluster_tree()
            self.add_log(f"Modified cluster: {name}")
            w.destroy()

        row = ttk.Frame(f); row.grid(row=2, column=0, columnspan=2, pady=18)
        ttk.Button(row, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Cancel", command=w.destroy).pack(side=tk.LEFT, padx=5)

    def delete_cluster(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to delete")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        if not messagebox.askyesno("Confirm Delete", f"Delete cluster '{name}' and all its satellites?"):
            return

        cluster = next((c for c in self.clusters if c["name"] == name), None)
        if not cluster:
            return

        for sat_name in list(cluster["satellites"]):
            self.satellites[:] = [s for s in self.satellites if s["name"] != sat_name]
            for orbit in self.orbit_configurations:
                if sat_name in orbit["satellites"]:
                    orbit["satellites"].remove(sat_name)

        self.clusters.remove(cluster)

        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        if hasattr(self.parent_app, "update_satellite_dropdowns"):
            self.parent_app.update_satellite_dropdowns()
        if hasattr(self.parent_app, "update_status_counts"):
            self.parent_app.update_status_counts()

        self.add_log(f"Deleted cluster: {name}")

    def test_cluster_communication(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to test")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        cluster = next((c for c in self.clusters if c["name"] == name), None)
        if not cluster:
            return

        msgs = [f"Testing communication in cluster: {name}",
                f"Leader: {cluster['leader']} broadcasting to children..."]
        for child in cluster["children"]:
            msgs.append(f"  → {child}: Message received")
        msgs.append("\nChildren reporting to leader:")
        for child in cluster["children"]:
            msgs.append(f"  ← {child}: Status OK")

        messagebox.showinfo("Communication Test", "\n".join(msgs))
        self.add_log(f"Communication test completed for cluster: {name}")

    def _add_satellite_to_cluster(self, cluster_name):
        cluster = next((c for c in self.clusters if c["name"] == cluster_name), None)
        if not cluster:
            return

        choices = [o["name"] for o in self.orbit_configurations]
        orbit_name = simpledialog.askstring("Select Orbit", f"Enter orbit name for new satellite\nAvailable: {', '.join(choices)}")
        if not orbit_name or orbit_name not in choices:
            return
        orbit = next(o for o in self.orbit_configurations if o["name"] == orbit_name)

        child_num = len(cluster["children"]) + 1
        sat_name = f"{cluster_name}_Child{child_num}"

        satellite = {
            "name": sat_name,
            "type": "cluster_member",
            "cluster": cluster_name,
            "role": "child",
            "formation": cluster["formation"],
            "orbit": {
                "a": orbit["altitude"] + 6371,
                "e": 0.01,
                "i": orbit["inclination"],
                "Omega": orbit["raan"],
                "omega": 0.0,
                "f": child_num * 10.0,
            },
            "fault": {"type": "friction", "magnitude": 0.0005, "wheel": 3, "time": 15.0, "enabled": False},
            "camera": {"position": [0.0, 0.0, 15.0], "fov": 80.0, "enabled": False},
            "communication": {"range": float(self.comm_range.get()), "fov": float(self.comm_fov.get()), "aHat_B": [0.0, 0.0, -1.0]},
            "targets": [],
            "orbit_name": orbit_name,
            "position_offset_km": [0.0, 0.0, 0.0],
        }

        self.satellites.append(satellite)
        cluster["satellites"].append(sat_name)
        cluster["children"].append(sat_name)
        orbit["satellites"].append(sat_name)

        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        self.add_log(f"Added {sat_name} to cluster {cluster_name}")

    def _update_cluster_formation(self, cluster):
        n = len(cluster["satellites"])
        formation = cluster["formation"]
        separation = cluster.get("separation", 10.0)

        # compute new leader-relative offsets (visual) and TA adjustments (kept minimal)
        positions = self._calculate_cluster_positions(n, formation, leader_anomaly=0.0, phase_offset=0.0, separation=separation)

        for idx, sat_name in enumerate(cluster["satellites"]):
            sat = next((s for s in self.satellites if s["name"] == sat_name), None)
            if not sat:
                continue
            sat["orbit"]["f"] += positions[idx]["anomaly"]  # small tweak
            sat["position_offset_km"] = positions[idx].get("offset", [0.0, 0.0, 0.0])

    def _calculate_cluster_positions(self, num_sats, formation, leader_anomaly, phase_offset, separation):
        # Normalize to 4 formations
        f = (formation or "").strip().lower()
        if f == "line":
            key = "column"
        elif f == "triangle":
            key = "box"  # close visual
        elif f == "diamond":
            key = "box"
        else:
            key = "train"  # Leader-Follower

        positions = []

        def leader_entry():
            return {"anomaly": leader_anomaly, "offset": [0.0, 0.0, 0.0]}

        if key == "column":  # line
            dA = max(2.0, separation * 0.5)
            positions.append(leader_entry())
            for k in range(1, num_sats):
                positions.append({"anomaly": leader_anomaly + k * dA, "offset": [0.0, k * 0.2, 0.0]})

        elif key == "box":  # diamond/triangle
            side = max(3.0, separation)
            positions.append(leader_entry())
            ring1 = [( side,  0.0), (0.0,  side), (-side, 0.0), (0.0, -side)]
            ring2 = [( 2*side,  2*side), (-2*side,  2*side), (-2*side, -2*side), ( 2*side, -2*side)]
            coords = ring1 + ring2
            for (dx, dy) in coords[:max(0, num_sats - 1)]:
                positions.append({"anomaly": leader_anomaly, "offset": [dx, dy, 0.0]})
            while len(positions) < num_sats:
                angle = 2 * np.pi * (len(positions) - 1) / max(1, num_sats - 1)
                positions.append({"anomaly": leader_anomaly, "offset": [1.5 * side * np.cos(angle), 1.5 * side * np.sin(angle), 0.0]})

        else:  # train (Leader-Follower)
            dA = max(8.0, separation)
            positions.append(leader_entry())
            for k in range(1, num_sats):
                positions.append({"anomaly": leader_anomaly + k * dA, "offset": [0.0, 0.0, 0.0]})

        return positions

    # ------------------------------------------------------------ Formation view
    def view_formation(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showwarning("No Selection", "Please select a cluster to view")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        self.show_formation_check_plot(name)

    def show_formation_check_plot(self, cluster_name):
        cluster = next((c for c in self.clusters if c["name"] == cluster_name), None)
        if not cluster:
            return

        # Gather positions: use saved 'position_offset_km' (leader-relative)
        positions = []
        for nm in cluster["satellites"]:
            sat = next((s for s in self.satellites if s["name"] == nm), None)
            pos = sat.get("position_offset_km", [0.0, 0.0, 0.0]) if sat else [0.0, 0.0, 0.0]
            positions.append(pos)

        # Compute simple RMS from ideal (0 for leader; others target at separation distance)
        rms_error = self.calculate_formation_error(cluster, positions)

        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Formation Check: {cluster_name}")
        w.geometry("900x700")

        fig = Figure(figsize=(10, 8), facecolor="white")

        # 3D leader-relative scatter
        ax1 = fig.add_subplot(221, projection="3d")
        ax1.set_title(f"{cluster_name} — {cluster['formation']} (Leader-Relative)")
        # Ideal template lines (thin)
        self.plot_ideal_formation(ax1, cluster["formation"], len(cluster["satellites"]), cluster.get("separation", 10.0))

        for i, nm in enumerate(cluster["satellites"]):
            x, y, z = positions[i]
            if nm == cluster["leader"]:
                ax1.scatter(x, y, z, c="red", s=100, marker="s", label="Leader")
            else:
                ax1.scatter(x, y, z, c=cluster["color_hex"], s=50, marker="o")
        ax1.set_xlabel("Along-track (km)")
        ax1.set_ylabel("Cross-track (km)")
        ax1.set_zlabel("Radial (km)")
        ax1.legend()

        # Status panel
        ax2 = fig.add_subplot(222); ax2.axis("off")
        if rms_error < 0.5:
            status_color, status_text = "green", "EXCELLENT"
        elif rms_error < 2.0:
            status_color, status_text = "yellow", "GOOD"
        else:
            status_color, status_text = "red", "NEEDS CORRECTION"
        metrics = f"""Formation Quality
Formation: {cluster['formation']}
Satellites: {len(cluster['satellites'])}
Leader: {cluster['leader']}
Separation: {cluster.get('separation', 10.0)} km

RMS Error: {rms_error:.3f} km
Status: {status_text}"""
        ax2.text(0.1, 0.5, metrics, fontsize=11, family="monospace",
                 bbox=dict(boxstyle="round,pad=0.6", facecolor=status_color, alpha=0.25))

        # Error vs time (synthetic)
        ax3 = fig.add_subplot(223)
        t = np.linspace(0, 24, 100)
        e = rms_error + 0.2 * np.sin(2 * np.pi * t / 12)
        ax3.plot(t, e, linewidth=2)
        ax3.axhline(0.5, color="g", linestyle="--", label="Excellent")
        ax3.axhline(2.0, color="y", linestyle="--", label="Good")
        ax3.set_xlabel("Time (hours)")
        ax3.set_ylabel("RMS Error (km)")
        ax3.set_title("Formation Error Over Time")
        ax3.grid(True, alpha=0.3); ax3.legend()

        # 2D intra-cluster comm sketch
        ax4 = fig.add_subplot(224)
        ax4.set_title("Intra-Cluster Communication (schematic)")
        leader_pos = positions[0] if positions else [0, 0, 0]
        for i, nm in enumerate(cluster["children"]):
            idx = cluster["satellites"].index(nm)
            child_pos = positions[idx]
            ax4.plot([leader_pos[0], child_pos[0]], [leader_pos[1], child_pos[1]], color=cluster["color_hex"], alpha=0.6)
        ax4.scatter(leader_pos[0], leader_pos[1], c="red", s=120, marker="s", label="Leader")
        for i, nm in enumerate(cluster["children"]):
            idx = cluster["satellites"].index(nm)
            x, y, _ = positions[idx]
            ax4.scatter(x, y, c=cluster["color_hex"], s=70)
        ax4.set_xlabel("Along-track (km)"); ax4.set_ylabel("Cross-track (km)")
        ax4.grid(True, alpha=0.3); ax4.legend()

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=w)
        canvas.draw(); canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Save
        bar = ttk.Frame(w); bar.pack(fill=tk.X, pady=6)
        def save_plot():
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            fname = f"plots/{ts}_FormationCheck_{cluster_name}.png"
            fig.savefig(fname, dpi=150, bbox_inches="tight")
            messagebox.showinfo("Saved", f"Plot saved to {fname}")
        ttk.Button(bar, text="Save Plot", command=save_plot).pack()

    def plot_ideal_formation(self, ax, formation_type, num_sats, separation):
        f = (formation_type or "").lower()
        if "line" in f:
            for i in range(num_sats):
                ax.plot([i * separation, i * separation], [0, 0], [-0.1, 0.1], "g--", alpha=0.3, linewidth=1)
        elif "triangle" in f:
            ang = np.linspace(0, 2*np.pi, 3, endpoint=False)
            x = separation * np.cos(ang); y = separation * np.sin(ang)
            for i in range(3):
                j = (i + 1) % 3
                ax.plot([x[i], x[j]], [y[i], y[j]], [0, 0], "g--", alpha=0.3)
        elif "diamond" in f:
            pts = [[ separation, 0, 0],[0, separation, 0],[-separation, 0, 0],[0,-separation, 0]]
            for i in range(4):
                j = (i + 1) % 4
                ax.plot([pts[i][0], pts[j][0]],[pts[i][1], pts[j][1]],[0,0],"g--", alpha=0.3)
        else:  # Leader-Follower
            for i in range(num_sats):
                ax.plot([-i * separation, -i * separation], [0, 0], [-0.1, 0.1], "g--", alpha=0.3, linewidth=1)

    def calculate_formation_error(self, cluster, positions):
        # Simple RMS distance from ideal ring/line center (leader at origin)
        if not positions:
            return 0.0
        sep = float(cluster.get("separation", 10.0))
        f = (cluster.get("formation","") or "").lower()
        errs = []
        for i, p in enumerate(positions):
            x, y, z = p
            if i == 0:
                errs.append(np.linalg.norm([x, y, z]))  # leader ~ 0
                continue
            if "line" in f or "leader-follower" in f:
                target = np.array([i*sep, 0.0, 0.0])
            elif "diamond" in f:
                targets = [np.array([ sep,0,0]), np.array([0, sep,0]), np.array([-sep,0,0]), np.array([0,-sep,0])]
                target = targets[(i-1) % len(targets)]
            elif "triangle" in f:
                ang = 2*np.pi*((i-1) % 3)/3
                target = np.array([sep*np.cos(ang), sep*np.sin(ang), 0.0])
            else:
                target = np.array([0.0, 0.0, 0.0])
            errs.append(np.linalg.norm(np.array([x,y,z]) - target))
        return float(np.sqrt(np.mean(np.square(errs))))

    def export_config(self):
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"cluster_config_{ts}.txt"
        with open(fname, "w", encoding="utf-8") as f:
            f.write("Cluster Configuration Export\n")
            f.write(f"Generated: {datetime.now()}\n")
            f.write(f"Total Clusters: {len(self.clusters)}/{self.MAX_CLUSTERS}\n\n")
            for c in self.clusters:
                f.write(f"Cluster: {c['name']}\n")
                f.write(f"  Formation: {c['formation']}\n")
                f.write(f"  Satellites: {len(c['satellites'])}\n")
                f.write(f"  Leader: {c['leader']}\n")
                f.write(f"  Children: {len(c['children'])}\n")
                f.write(f"  Separation: {c.get('separation', 10.0)} km\n")
                f.write(f"  Color: {c.get('color_hex','')}\n\n")
        messagebox.showinfo("Export Complete", f"Configuration exported to {fname}")

    # ------------------------------------------------------------ Orbit Tab
    def _create_orbit_management_tab(self, parent):
        create_frame = ttk.LabelFrame(parent, text="Create New Orbit", padding=10)
        create_frame.pack(fill=tk.X, padx=10, pady=10)

        grid = ttk.Frame(create_frame); grid.pack(fill=tk.X, pady=5)
        ttk.Label(grid, text="Orbit Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.new_orbit_name_var = tk.StringVar()
        ttk.Entry(grid, textvariable=self.new_orbit_name_var, width=20).grid(row=0, column=1, padx=5)

        ttk.Label(grid, text="Altitude (km):").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.new_orbit_alt_var = tk.DoubleVar(value=600.0)
        ttk.Entry(grid, textvariable=self.new_orbit_alt_var, width=15).grid(row=0, column=3, padx=5)

        ttk.Label(grid, text="Inclination (deg):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.new_orbit_inc_var = tk.DoubleVar(value=55.0)
        ttk.Entry(grid, textvariable=self.new_orbit_inc_var, width=15).grid(row=1, column=1, padx=5)

        ttk.Label(grid, text="RAAN (deg):").grid(row=1, column=2, sticky=tk.W, padx=5)
        self.new_orbit_raan_var = tk.DoubleVar(value=0.0)
        ttk.Entry(grid, textvariable=self.new_orbit_raan_var, width=15).grid(row=1, column=3, padx=5)

        ttk.Button(create_frame, text="Create Orbit", command=self.create_new_orbit).pack(pady=10)

        orbits_frame = ttk.LabelFrame(parent, text="Orbit Library", padding=10)
        orbits_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = ("Name", "Altitude", "Inclination", "RAAN", "Satellites")
        self.orbit_tree = ttk.Treeview(orbits_frame, columns=columns, show="headings", height=8)
        for col in columns:
            self.orbit_tree.heading(col, text=col)
        self.orbit_tree.column("Name", width=140)
        self.orbit_tree.column("Altitude", width=100)
        self.orbit_tree.column("Inclination", width=100)
        self.orbit_tree.column("RAAN", width=100)
        self.orbit_tree.column("Satellites", width=100)
        self.orbit_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        sb = ttk.Scrollbar(orbits_frame, orient="vertical", command=self.orbit_tree.yview)
        self.orbit_tree.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        row = ttk.Frame(orbits_frame); row.pack(fill=tk.X, pady=5)
        ttk.Button(row, text="Modify Orbit", command=self.modify_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Delete Orbit", command=self.delete_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="View Satellites", command=self.view_orbit_satellites).pack(side=tk.LEFT, padx=5)

        self.update_orbit_tree()

    def create_new_orbit(self):
        name = self.new_orbit_name_var.get().strip()
        if not name:
            messagebox.showwarning("Missing Name", "Please provide an orbit name")
            return
        if name in [o["name"] for o in self.orbit_configurations]:
            messagebox.showwarning("Duplicate Name", f"Orbit '{name}' already exists")
            return

        new_orbit = {
            "name": name,
            "altitude": float(self.new_orbit_alt_var.get()),
            "inclination": float(self.new_orbit_inc_var.get()),
            "raan": float(self.new_orbit_raan_var.get()),
            "satellites": [],
            "description": f"Custom orbit at {self.new_orbit_alt_var.get()} km",
        }
        self.orbit_configurations.append(new_orbit)
        self.update_orbit_tree()
        self.add_log(f"Created new orbit: {name}")

        # Reset form
        self.new_orbit_name_var.set("")
        self.new_orbit_alt_var.set(600.0)
        self.new_orbit_inc_var.set(55.0)
        self.new_orbit_raan_var.set(0.0)

    def modify_orbit(self):
        sel = self.orbit_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select an orbit to modify")
            return
        name = self.orbit_tree.item(sel[0])["values"][0]
        orbit = next((o for o in self.orbit_configurations if o["name"] == name), None)
        if not orbit:
            return

        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Modify Orbit: {name}")
        w.geometry("360x260")
        f = ttk.Frame(w, padding=20); f.pack(fill=tk.BOTH, expand=True)

        ttk.Label(f, text="Altitude (km):").grid(row=0, column=0, sticky=tk.W, pady=5)
        alt = tk.DoubleVar(value=orbit["altitude"]); ttk.Entry(f, textvariable=alt, width=20).grid(row=0, column=1, pady=5)
        ttk.Label(f, text="Inclination (deg):").grid(row=1, column=0, sticky=tk.W, pady=5)
        inc = tk.DoubleVar(value=orbit["inclination"]); ttk.Entry(f, textvariable=inc, width=20).grid(row=1, column=1, pady=5)
        ttk.Label(f, text="RAAN (deg):").grid(row=2, column=0, sticky=tk.W, pady=5)
        raan = tk.DoubleVar(value=orbit["raan"]); ttk.Entry(f, textvariable=raan, width=20).grid(row=2, column=1, pady=5)

        def apply():
            orbit["altitude"] = float(alt.get())
            orbit["inclination"] = float(inc.get())
            orbit["raan"] = float(raan.get())
            # push changes to satellites on this orbit
            for nm in orbit["satellites"]:
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    sat["orbit"]["a"] = orbit["altitude"] + 6371
                    sat["orbit"]["i"] = orbit["inclination"]
                    sat["orbit"]["Omega"] = orbit["raan"]
            self.update_orbit_tree()
            self.add_log(f"Modified orbit: {name}")
            w.destroy()

        row = ttk.Frame(f); row.grid(row=3, column=0, columnspan=2, pady=18)
        ttk.Button(row, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Cancel", command=w.destroy).pack(side=tk.LEFT, padx=5)

    def delete_orbit(self):
        sel = self.orbit_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select an orbit to delete")
            return
        name = self.orbit_tree.item(sel[0])["values"][0]
        orbit = next((o for o in self.orbit_configurations if o["name"] == name), None)
        if not orbit:
            return
        if orbit["satellites"]:
            messagebox.showwarning("Cannot Delete", f"Cannot delete orbit with {len(orbit['satellites'])} satellites")
            return
        if not messagebox.askyesno("Confirm Delete", f"Delete orbit '{name}'?"):
            return
        self.orbit_configurations.remove(orbit)
        self.update_orbit_tree()
        self.add_log(f"Deleted orbit: {name}")

    def view_orbit_satellites(self):
        sel = self.orbit_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select an orbit")
            return
        name = self.orbit_tree.item(sel[0])["values"][0]
        orbit = next((o for o in self.orbit_configurations if o["name"] == name), None)
        if not orbit:
            return
        if orbit["satellites"]:
            txt = "\n".join([f"• {s}" for s in orbit["satellites"]])
            messagebox.showinfo(f"Satellites in {name}", f"Satellites ({len(orbit['satellites'])}):\n\n{txt}")
        else:
            messagebox.showinfo(f"Satellites in {name}", "No satellites in this orbit")

    # ------------------------------------------------------------ Individuals
    def _create_individual_tab(self, parent):
        add = ttk.LabelFrame(parent, text="Add Individual Satellite", padding=10)
        add.pack(fill=tk.X, padx=10, pady=10)

        row = ttk.Frame(add); row.pack(fill=tk.X, pady=5)
        ttk.Label(row, text="Satellite Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.individual_sat_name_var = tk.StringVar()
        ttk.Entry(row, textvariable=self.individual_sat_name_var, width=20).grid(row=0, column=1, padx=5)

        ttk.Label(row, text="Orbit:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.individual_orbit_var = tk.StringVar()
        orbit_combo = ttk.Combobox(row, textvariable=self.individual_orbit_var, width=20, state="readonly")
        orbit_combo.grid(row=0, column=3, padx=5)
        orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        if orbit_combo["values"]:
            orbit_combo.current(0)

        ttk.Label(row, text="True Anomaly (deg):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.individual_anomaly_var = tk.DoubleVar(value=0.0)
        ttk.Entry(row, textvariable=self.individual_anomaly_var, width=20).grid(row=1, column=1, padx=5)

        ttk.Button(add, text="Add Satellite", command=self.add_individual_satellite).pack(pady=10)

        list_frame = ttk.LabelFrame(parent, text="All Satellites", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        container = ttk.Frame(list_frame); container.pack(fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(container); sb.pack(side=tk.RIGHT, fill=tk.Y)

        self.satellite_listbox = tk.Listbox(
            container, selectmode=tk.SINGLE, yscrollcommand=sb.set, font=("Segoe UI", 10)
        )
        self.satellite_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.config(command=self.satellite_listbox.yview)

        actions = ttk.Frame(list_frame); actions.pack(fill=tk.X, pady=5)
        ttk.Button(actions, text="Remove Satellite", command=self.remove_satellite).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Assign to Orbit", command=self.assign_satellite_to_orbit).pack(side=tk.LEFT, padx=5)

        self.update_satellite_listbox()

    def add_individual_satellite(self):
        name = self.individual_sat_name_var.get().strip()
        orb_name = self.individual_orbit_var.get()
        anomaly = float(self.individual_anomaly_var.get())

        if not name:
            messagebox.showwarning("Missing Name", "Please provide a satellite name")
            return
        if name in [s["name"] for s in self.satellites]:
            messagebox.showwarning("Duplicate Name", f"Satellite '{name}' already exists")
            return

        orbit = next((o for o in self.orbit_configurations if o["name"] == orb_name), None)
        if not orbit:
            messagebox.showerror("Invalid Orbit", "Selected orbit not found")
            return

        sat = {
            "name": name,
            "type": "individual",
            "cluster": None,
            "role": "independent",
            "orbit": {
                "a": orbit["altitude"] + 6371, "e": 0.01, "i": orbit["inclination"],
                "Omega": orbit["raan"], "omega": 0.0, "f": anomaly
            },
            "fault": {"type": "friction", "magnitude": 0.0005, "wheel": 3, "time": 15.0, "enabled": False},
            "camera": {"position": [0.0, 0.0, 15.0], "fov": 80.0, "enabled": True},
            "communication": {"range": float(self.comm_range.get()), "fov": float(self.comm_fov.get()), "aHat_B": [0.0, 0.0, -1.0]},
            "targets": [],
            "orbit_name": orb_name,
        }
        self.satellites.append(sat)
        orbit["satellites"].append(name)

        self.update_satellite_listbox()
        self.update_orbit_tree()
        if hasattr(self.parent_app, "update_satellite_dropdowns"):
            self.parent_app.update_satellite_dropdowns()
        if hasattr(self.parent_app, "update_status_counts"):
            self.parent_app.update_status_counts()

        self.individual_sat_name_var.set("")
        self.individual_anomaly_var.set(0.0)
        self.add_log(f"Added individual satellite: {name}")

    def remove_satellite(self):
        sel = self.satellite_listbox.curselection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a satellite to remove")
            return
        item_text = self.satellite_listbox.get(sel[0])

        if "[INDIVIDUAL]" not in item_text:
            messagebox.showwarning("Cannot Remove", "Cluster satellites must be removed via Cluster Management")
            return

        sat_name = item_text.split("]")[1].strip().split("(")[0].strip()
        if not messagebox.askyesno("Confirm Remove", f"Remove satellite '{sat_name}'?"):
            return

        sat = next((s for s in self.satellites if s["name"] == sat_name), None)
        if sat:
            self.satellites.remove(sat)
            orbit_name = sat.get("orbit_name")
            if orbit_name:
                orbit = next((o for o in self.orbit_configurations if o["name"] == orbit_name), None)
                if orbit and sat_name in orbit["satellites"]:
                    orbit["satellites"].remove(sat_name)

            self.update_satellite_listbox()
            self.update_orbit_tree()
            if hasattr(self.parent_app, "update_satellite_dropdowns"):
                self.parent_app.update_satellite_dropdowns()
            if hasattr(self.parent_app, "update_status_counts"):
                self.parent_app.update_status_counts()

            self.add_log(f"Removed satellite: {sat_name}")

    def assign_satellite_to_orbit(self):
        sel = self.satellite_listbox.curselection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a satellite")
            return
        item_text = self.satellite_listbox.get(sel[0])
        if "[INDIVIDUAL]" not in item_text:
            messagebox.showinfo("Cannot Reassign", "Only individual satellites can be reassigned")
            return

        sat_name = item_text.split("]")[1].strip().split("(")[0].strip()
        sat = next((s for s in self.satellites if s["name"] == sat_name), None)
        if not sat:
            return

        choices = [o["name"] for o in self.orbit_configurations]
        new_name = simpledialog.askstring("Select Orbit", f"Enter new orbit for {sat_name}\nAvailable: {', '.join(choices)}")
        if not new_name or new_name not in choices:
            return

        old = sat.get("orbit_name")
        new_orbit = next(o for o in self.orbit_configurations if o["name"] == new_name)
        if old:
            old_orbit = next((o for o in self.orbit_configurations if o["name"] == old), None)
            if old_orbit and sat_name in old_orbit["satellites"]:
                old_orbit["satellites"].remove(sat_name)

        sat["orbit_name"] = new_name
        sat["orbit"]["a"] = new_orbit["altitude"] + 6371
        sat["orbit"]["i"] = new_orbit["inclination"]
        sat["orbit"]["Omega"] = new_orbit["raan"]
        new_orbit["satellites"].append(sat_name)

        self.update_satellite_listbox()
        self.update_orbit_tree()
        self.add_log(f"Reassigned {sat_name} to orbit: {new_name}")

    # ------------------------------------------------------------ Comm Tab
    def _create_comm_tab(self, parent):
        settings = ttk.LabelFrame(parent, text="Communication Settings", padding=10)
        settings.pack(fill=tk.X, padx=10, pady=10)

        row = ttk.Frame(settings); row.pack(fill=tk.X, pady=5)
        ttk.Label(row, text="Range (km):").grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Entry(row, textvariable=self.comm_range, width=15).grid(row=0, column=1, padx=5)
        ttk.Label(row, text="FOV (deg):").grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Entry(row, textvariable=self.comm_fov, width=15).grid(row=0, column=3, padx=5)
        ttk.Button(row, text="Refresh", command=lambda: [self.update_communication_plot(), self.comm_canvas.draw()]).grid(row=0, column=4, padx=5)
        ttk.Button(row, text="Help", command=self.show_communication_help).grid(row=0, column=5, padx=5)
        ttk.Label(row, text="(Bars = comm possible, Dots = messages)", font=("Arial", 8), foreground="gray").grid(row=1, column=0, columnspan=6, sticky=tk.W)

        plot_frame = ttk.LabelFrame(parent, text="Communication Windows", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.comm_figure = Figure(figsize=(10, 6), dpi=80)
        self.comm_ax = self.comm_figure.add_subplot(111)
        self.comm_canvas = FigureCanvasTkAgg(self.comm_figure, plot_frame)
        self.comm_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        ttk.Button(plot_frame, text="Update Communication Plot", command=self.update_communication_plot).pack(pady=8)

        self.update_communication_plot()

    def update_communication_plot(self):
        self.comm_ax.clear()

        if not self.clusters:
            self.comm_ax.text(0.5, 0.5, "No clusters configured\nCreate clusters in the Cluster Management tab",
                              ha="center", va="center", fontsize=12)
            self.comm_ax.set_xlim(0, 1); self.comm_ax.set_ylim(0, 1)
            self.comm_canvas.draw(); return

        t = np.linspace(0, 30, 100)
        y = 0
        title_text = f"Communication Windows for {len(self.clusters)} Clusters"
        self.comm_ax.text(15, -1.5, "Colored bars = communication windows | Dots = messages",
                          ha="center", fontsize=9, style="italic", color="gray")

        legend_elems = [Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", markersize=8, label="Message sent"),
                        Patch(facecolor="gray", alpha=0.3, label="Comm window")]

        for idx, c in enumerate(self.clusters):
            color = self.PALETTE[c.get("color_idx", idx % len(self.PALETTE))]
            leader = c.get("leader", "Leader")
            children = c.get("children", [])

            # cluster band
            self.comm_ax.axhspan(y - 0.5, y + len(children) + 0.5, color=color, alpha=0.06, zorder=0)
            self.comm_ax.text(-3, y + len(children) / 2, f"[{c['name'].upper()}]",
                              fontsize=11, fontweight="bold", color=color, va="center", ha="right",
                              bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=color, linewidth=2))

            # legend entry per cluster
            legend_elems.append(Patch(facecolor=color, edgecolor=color, alpha=0.3, label=f"{c['name']} cluster"))

            for j, child in enumerate(children):
                phase = j * np.pi / 4
                comm_window = np.sin(2 * np.pi * t / 10 + phase) > 0.3

                in_window = False; start = 0.0
                for k in range(len(t)):
                    if comm_window[k] and not in_window:
                        start = t[k]; in_window = True
                    elif not comm_window[k] and in_window:
                        self.comm_ax.barh(y, t[k] - start, left=start, height=0.6,
                                          color=color, alpha=0.3, edgecolor=color, linewidth=1)
                        in_window = False
                if in_window:
                    self.comm_ax.barh(y, t[-1] - start, left=start, height=0.6,
                                      color=color, alpha=0.3, edgecolor=color, linewidth=1)

                # message dots
                for mt in t[::12]:
                    k = np.argmin(np.abs(t - mt))
                    if comm_window[k]:
                        self.comm_ax.scatter(mt, y, color=color, s=80, marker="o", zorder=5, edgecolors="white", linewidths=1)

                self.comm_ax.text(-2, y, f"{leader[:12]} → {child[:12]}", fontsize=9, va="center", color="black")
                y += 1

            y += 0.8  # gap between clusters

        self.comm_ax.set_xlabel("Time (minutes)", fontsize=11, fontweight="bold")
        self.comm_ax.set_ylabel("Communication Links", fontsize=11, fontweight="bold")
        self.comm_ax.set_title(title_text, fontsize=13, fontweight="bold", pad=12)
        self.comm_ax.set_xlim(-3.5, 31); self.comm_ax.set_ylim(-2, max(1, y))
        self.comm_ax.grid(True, alpha=0.2, axis="x", linestyle="--")
        for tmark in [0, 5, 10, 15, 20, 25, 30]:
            self.comm_ax.axvline(x=tmark, color="gray", linewidth=0.5, alpha=0.3, linestyle=":")

        # dedupe legend labels
        seen = set(); final = []
        for h in legend_elems:
            lbl = h.get_label()
            if lbl not in seen:
                final.append(h); seen.add(lbl)

        self.comm_ax.legend(handles=final, loc="upper left", fontsize=9, framealpha=0.9, ncol=2)
        self.comm_canvas.draw()

    def show_communication_help(self):
        messagebox.showinfo(
            "Communication Plot Help",
            """WHAT YOU SEE
• Each row is a leader→child link
• Links are grouped by cluster color
• Bars = time windows when link is available
• Dots = actual messages (only inside windows)

WHY WINDOWS COME/GO
• Relative motion causes satellites to move in/out of range
• FOV and Earth occlusion also gate communication""",
        )

    # ------------------------------------------------------------ Update lists
    def update_cluster_tree(self):
        for item in self.cluster_tree.get_children():
            self.cluster_tree.delete(item)
        for c in self.clusters:
            status = "Active" if c.get("leader") else "Inactive"
            node = self.cluster_tree.insert("", "end", values=(c["name"], c.get("leader","None"),
                                                               len(c.get("children", [])), c.get("formation",""), status))
            # children rows with satellite details
            for nm in c.get("satellites", []):
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    self.cluster_tree.insert(node, "end", values=(f"  {sat['name']}",
                                                                  sat.get("role", ""), sat.get("orbit_name", ""),
                                                                  f"{sat['orbit']['a']-6371:.0f}km",
                                                                  "Fault" if sat["fault"]["enabled"] else "OK"))

    def update_orbit_tree(self):
        for item in self.orbit_tree.get_children():
            self.orbit_tree.delete(item)
        for o in self.orbit_configurations:
            self.orbit_tree.insert("", "end", values=(o["name"], f"{o['altitude']:.0f} km",
                                                      f"{o['inclination']:.1f}°", f"{o['raan']:.1f}°",
                                                      len(o["satellites"])))

    def update_satellite_listbox(self):
        self.satellite_listbox.delete(0, tk.END)
        # clusters
        for c in self.clusters:
            self.satellite_listbox.insert(tk.END, f"[CLUSTER] {c['name']}")
            for nm in c["satellites"]:
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    role = "L" if sat["role"] == "leader" else "C"
                    orbit = sat.get("orbit_name", "Unknown")
                    alt = sat["orbit"]["a"] - 6371
                    fault = "F" if sat["fault"]["enabled"] else ""
                    self.satellite_listbox.insert(tk.END, f"  [{role}] {sat['name']} ({orbit}, {alt:.0f}km) {fault}")
        # individuals
        for sat in self.satellites:
            if sat["type"] == "individual":
                orbit = sat.get("orbit_name", "Unknown")
                alt = sat["orbit"]["a"] - 6371
                fault = "F" if sat["fault"]["enabled"] else ""
                self.satellite_listbox.insert(tk.END, f"[INDIVIDUAL] {sat['name']} ({orbit}, {alt:.0f}km) {fault}")

    def clear_cluster_form(self):
        self.cluster_name_var.set("")
        self.sats_per_cluster_var.set(4)
        self.formation_var.set(self.ALLOWED_FORMATIONS[0])
        self.formation_separation.set(10.0)

    # ------------------------------------------------------------ Debug helper
    def debug_cluster_status(self):
        print("\n=== CLUSTER DEBUG INFO ===")
        print(f"Number of clusters: {len(self.clusters)}")
        for c in self.clusters:
            print(f"  Cluster '{c['name']}':")
            print(f"    - Leader: {c.get('leader', 'None')}")
            print(f"    - Children: {c.get('children', [])}")
            print(f"    - Formation: {c.get('formation', 'Unknown')}")
        print(f"Total satellites: {len(self.satellites)}")
        print("=========================\n")
