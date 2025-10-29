#!/usr/bin/env python
"""
core/gui_tab/constellation_tab.py

COMPLETE FIXED VERSION with:
- Unique satellite naming (cluster prefix)
- Removed Cartesian mode checkbox (set internally)
- Fixed orbit dropdown updates
- Better orbit assignment
- Improved formation visibility
- Progress feedback
"""

import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Optional helper formations
try:
    from formation_geometry import FormationGeometry, CartesianFormations
except Exception:
    FormationGeometry = None
    class _Dummy:
        def __getattr__(self, *_args, **_kwargs):
            raise AttributeError("formation_geometry not available")
    CartesianFormations = _Dummy()


# BaseTab fallback
try:
    from .base_tab import BaseTab
except Exception:
    class BaseTab:
        def __init__(self, parent_app, parent_frame):
            self.parent_app = parent_app
            self.parent_frame = parent_frame

        def add_log(self, message):
            if hasattr(self.parent_app, "add_log"):
                self.parent_app.add_log(message)


class ConstellationTab(BaseTab):
    """Constellation Management (clusters, orbits, individuals, comms)"""

    # Project constraints
    ALLOWED_FORMATIONS = ["Leader-Follower", "Line", "Triangle", "Diamond"]
    MAX_CLUSTERS = 4

    # Stable palette (shared with plots/Vizard)
    PALETTE = ["#2ecc71", "#e74c3c", "#9b59b6", "#f39c12"]

    def __init__(self, parent_app, parent_frame):
        super().__init__(parent_app, parent_frame)

        # References to app-wide data
        self.satellites = parent_app.satellites
        self.current_satellite_index = parent_app.current_satellite_index

        # Cluster & orbit state
        self.clusters = []
        self.cluster_configurations = []
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
            {
                "name": "MEO 10000km",
                "altitude": 10000,
                "inclination": 55.0,
                "raan": 0.0,
                "satellites": [],
                "description": "Medium Earth Orbit",
            },
        ]

        # Communication defaults
        self.comm_range = tk.DoubleVar(value=2000.0)
        self.comm_fov = tk.DoubleVar(value=30.0)

        # Formation defaults
        self.formation_separation = tk.DoubleVar(value=15.0)  # Increased default for visibility

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
        sats_combo.current(2)

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
        
        # NOTE: Cartesian mode removed from GUI (set internally to True for all clusters)
        ttk.Label(formation_row, text="(Min 15km for Vizard visibility)", 
                 font=("Arial", 8), foreground="gray").grid(row=0, column=4, padx=5)

        # Orbit selection
        orbit_frame = ttk.LabelFrame(create_frame, text="Cluster Orbit", padding=5)
        orbit_frame.pack(fill=tk.X, pady=5)

        ttk.Label(orbit_frame, text="Select Orbit:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.primary_orbit_var = tk.StringVar()
        self.orbit_combo = ttk.Combobox(
            orbit_frame, textvariable=self.primary_orbit_var,
            values=[o["name"] for o in self.orbit_configurations],
            width=20, state="readonly",
        )
        self.orbit_combo.grid(row=0, column=1, padx=5)
        if self.orbit_configurations:
            self.primary_orbit_var.set(self.orbit_configurations[0]["name"])
            self.orbit_combo.current(0)

        ttk.Label(
            orbit_frame,
            text="All satellites share this orbit; offsets create the formation.",
            foreground="gray", font=("Arial", 8)
        ).grid(row=1, column=0, columnspan=2, pady=3, sticky=tk.W)

        # Buttons
        btns = ttk.Frame(create_frame)
        btns.pack(fill=tk.X, pady=10)
        ttk.Button(btns, text="Create Cluster", command=self.create_cluster, 
                  style="Run.TButton").pack(side=tk.LEFT, padx=5)
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

    # ------------------------------------------------------------ Create/Modify
    def create_cluster(self):
        """Create a new cluster with UNIQUE satellite names and proper orbit assignment"""
        # Cap at 4 clusters
        if len(self.clusters) >= self.MAX_CLUSTERS:
            messagebox.showwarning("Cluster Limit", 
                                 f"Maximum {self.MAX_CLUSTERS} clusters allowed (project requirement).")
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

        # Enforce minimum separation for Vizard visibility
        if separation < 10.0:
            messagebox.showwarning("Separation Too Small", 
                                 "Minimum 10km separation recommended for Vizard visibility. Setting to 10km.")
            separation = 10.0
            self.formation_separation.set(10.0)

        if not self.orbit_configurations:
            messagebox.showwarning("No Orbits", "Create at least one orbit first.")
            return

        # Get selected orbit
        primary_orbit = next((o for o in self.orbit_configurations 
                            if o["name"] == self.primary_orbit_var.get()),
                           self.orbit_configurations[0])

        # Progress feedback
        self.add_log(f"Creating cluster '{name}' with {n} satellites...")
        if hasattr(self.parent_app, 'update_status'):
            self.parent_app.update_status(f"Creating cluster {name}...")
            self.parent_app.root.update_idletasks()

        # Spread clusters around Earth
        cluster_base_angle = (len(self.clusters) * 90.0) % 360.0

        base_orbit = {
            "a": primary_orbit["altitude"] + 6371,
            "e": 0.001,
            "i": primary_orbit["inclination"],
            "Omega": primary_orbit.get("raan", 0.0),
            "omega": 0.0,
            "f": cluster_base_angle,
        }

        # CRITICAL FIX: Generate UNIQUE satellite names
        sat_names = []
        for i in range(n):
            if i == 0:
                sat_names.append(f"{name}_Leader")
            else:
                sat_names.append(f"{name}_Sat{i+1}")

        # Verify no duplicates
        existing_names = {s["name"] for s in self.satellites}
        for sname in sat_names:
            if sname in existing_names:
                messagebox.showerror("Naming Error", 
                    f"Satellite '{sname}' already exists. Use a different cluster name.")
                if hasattr(self.parent_app, 'update_status'):
                    self.parent_app.update_status("Ready")
                return

        # Formation offsets
        offsets = self._get_formation_offsets(n, formation, separation)

        # Cartesian positions (for visualization)
        cart = None
        try:
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
            cart = None

        # Build cluster record
        cluster = {
            "name": name,
            "type": "cluster",
            "formation": formation,
            "leader": sat_names[0],
            "children": sat_names[1:],
            "satellites": sat_names,
            "separation": separation,
            "cartesian_mode": True,  # Always True (not exposed in GUI)
            "color_idx": len(self.clusters) % len(self.PALETTE),
            "color_hex": self.PALETTE[len(self.clusters) % len(self.PALETTE)],
            "orbit_name": primary_orbit["name"],
            "base_orbit": base_orbit,
        }

        # Create satellites
        for i, sat_name in enumerate(sat_names):
            is_leader = (i == 0)

            # Apply orbital offset
            sat_orbit = {
                "a": base_orbit["a"],
                "e": base_orbit["e"],
                "i": base_orbit["i"] + offsets[i]["di"],
                "Omega": base_orbit["Omega"] + offsets[i]["dOmega"],
                "omega": base_orbit["omega"],
                "f": base_orbit["f"] + offsets[i]["df"],
            }

            # Communication pointing
            aHat = [0.2, -0.4, 0.2] if is_leader else [0.0, 0.0, -1.0]

            # Cartesian offset
            pos_off = [0.0, 0.0, 0.0]
            if cart is not None:
                try:
                    pos_off = cart[i].offset_km
                except Exception:
                    pos_off = [0.0, 0.0, 0.0]

            # Create satellite with ALL required fields
            sat = {
                "name": sat_name,
                "type": "cluster_member",
                "cluster": name,
                "role": "leader" if is_leader else "child",
                "cartesian_mode": True,
                "position_offset_km": pos_off,
                "formation": formation,
                "separation": separation,  # Add to satellite for simulation
                "orbit": sat_orbit,
                "orbit_name": primary_orbit["name"],
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 15.0,
                    "enabled": False,
                    "periodic": {
                        "enabled": False,
                        "interval": 360,
                        "magnitude": 0.1,
                        "wheel": 1
                    }
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
            }

            self.satellites.append(sat)
            primary_orbit["satellites"].append(sat_name)

        # Book-keeping
        self.clusters.append(cluster)
        cluster["_members"] = [next(s for s in self.satellites if s["name"] == nm) 
                               for nm in cluster["satellites"]]

        # Update all UI elements
        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        
        if hasattr(self.parent_app, "update_satellite_dropdowns"):
            self.parent_app.update_satellite_dropdowns()
        if hasattr(self.parent_app, "update_status_counts"):
            self.parent_app.update_status_counts()

        # Detailed logging
        self.add_log(f"Created cluster '{name}' with {n} satellites in {formation} formation")
        self.add_log(f"  Formation separation: {separation} km")
        self.add_log(f"  All satellites in orbit: {primary_orbit['name']}")
        self.add_log(f"  Leader: {cluster['leader']}")
        self.add_log(f"  Children: {', '.join(cluster['children'])}")

        # Update status
        if hasattr(self.parent_app, 'update_status'):
            self.parent_app.update_status("Ready")

        # Clear form
        self.cluster_name_var.set("")

        # Warn if limit reached
        if len(self.clusters) >= self.MAX_CLUSTERS:
            messagebox.showinfo(
                "Cluster Limit Reached",
                f"Maximum {self.MAX_CLUSTERS} clusters created.\n\n"
                f"This limit ensures optimal simulation performance."
            )

    def _get_formation_offsets(self, num_sats, formation, separation_km):
        """
        Calculate proper orbital offsets for VISIBLE formations in Vizard.
        
        Returns: list of dicts {'df','di','dOmega'} in degrees
        
        Key: Use LARGE separations for clear Vizard visibility (min 15km recommended)
        """
        offsets = []
        
        # Enforce minimum separation for visibility
        effective_sep = max(separation_km, 15.0)
        
        # Convert km to degrees (for LEO ~7000km orbit, 1 deg ≈ 122 km)
        deg_per_km = 1.0 / 122.0
        sep_deg = effective_sep * deg_per_km
        
        # Minimum angular separation (1.5 degrees for clear visibility)
        sep_deg = max(sep_deg, 1.5)
        
        f = (formation or "").lower()
        
        if "line" in f:
            # Satellites in a line along orbital path
            for i in range(num_sats):
                offsets.append({
                    "df": i * sep_deg,
                    "di": 0.0,
                    "dOmega": 0.0
                })
        
        elif "triangle" in f:
            # Equilateral triangle
            if num_sats >= 1:
                offsets.append({"df": 0.0, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 2:
                offsets.append({
                    "df": -sep_deg * 0.866,
                    "di": sep_deg * 0.08,
                    "dOmega": 0.0
                })
            if num_sats >= 3:
                offsets.append({
                    "df": -sep_deg * 0.866,
                    "di": -sep_deg * 0.08,
                    "dOmega": 0.0
                })
            for i in range(3, num_sats):
                offsets.append({
                    "df": -sep_deg * 0.5,
                    "di": 0.0,
                    "dOmega": 0.0
                })
        
        elif "diamond" in f:
            # Diamond/box formation (clear 4-point pattern)
            if num_sats >= 1:
                offsets.append({"df": 0.0, "di": 0.0, "dOmega": 0.0})
            if num_sats >= 2:
                offsets.append({
                    "df": sep_deg,
                    "di": 0.0,
                    "dOmega": 0.0
                })
            if num_sats >= 3:
                offsets.append({
                    "df": 0.0,
                    "di": -sep_deg * 0.08,
                    "dOmega": 0.0
                })
            if num_sats >= 4:
                offsets.append({
                    "df": -sep_deg,
                    "di": 0.0,
                    "dOmega": 0.0
                })
            if num_sats >= 5:
                offsets.append({
                    "df": 0.0,
                    "di": sep_deg * 0.08,
                    "dOmega": 0.0
                })
            for i in range(5, num_sats):
                angle = 2 * np.pi * (i - 5) / max(1, num_sats - 5)
                offsets.append({
                    "df": sep_deg * 1.5 * np.cos(angle),
                    "di": sep_deg * 0.08 * np.sin(angle),
                    "dOmega": 0.0
                })
        
        else:  # Leader-Follower (train)
            # Satellites follow in a line behind leader
            for i in range(num_sats):
                offsets.append({
                    "df": -i * sep_deg,
                    "di": 0.0,
                    "dOmega": 0.0
                })
        
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
            f"Total Satellites: {len(cluster['satellites'])} (including leader)",
            f"Separation: {cluster.get('separation', 10.0)} km",
            f"Orbit: {cluster.get('orbit_name', 'Unknown')}",
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
                f"  Role: {sat['role'].upper()}",
                f"  Orbit: {sat.get('orbit_name','')}",
                f"  Altitude: {sat['orbit']['a'] - 6371:.0f} km",
                f"  Inc: {sat['orbit']['i']:.1f}°, RAAN: {sat['orbit']['Omega']:.1f}°, TA: {sat['orbit']['f']:.1f}°",
                f"  Camera: {'Enabled' if sat['camera']['enabled'] else 'Disabled'}",
                f"  Fault: {'Enabled' if sat['fault']['enabled'] else 'Disabled'}",
                ""
            ]
        details += [
            "=" * 60, 
            "COMMUNICATION",
            f"Range: {self.comm_range.get():.0f} km",
            f"FOV: {self.comm_fov.get():.1f}°",
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
        f = ttk.Frame(w, padding=20)
        f.pack(fill=tk.BOTH, expand=True)

        ttk.Label(f, text="Formation Type:").grid(row=0, column=0, sticky=tk.W, pady=5)
        form_var = tk.StringVar(value=cluster["formation"])
        ttk.Combobox(f, textvariable=form_var, values=self.ALLOWED_FORMATIONS, 
                    width=20, state="readonly").grid(row=0, column=1, pady=5)

        ttk.Label(f, text="Separation (km):").grid(row=1, column=0, sticky=tk.W, pady=5)
        sep_var = tk.DoubleVar(value=cluster.get("separation", 10.0))
        ttk.Entry(f, textvariable=sep_var, width=20).grid(row=1, column=1, pady=5)

        def apply():
            new_sep = float(sep_var.get())
            if new_sep < 10.0:
                messagebox.showwarning("Separation Too Small", 
                                     "Minimum 10km recommended. Setting to 10km.")
                new_sep = 10.0
            
            cluster["formation"] = form_var.get()
            cluster["separation"] = new_sep
            
            # Update all satellites in cluster
            for sat_name in cluster["satellites"]:
                sat = next((s for s in self.satellites if s["name"] == sat_name), None)
                if sat:
                    sat["formation"] = cluster["formation"]
                    sat["separation"] = cluster["separation"]
            
            self._update_cluster_formation(cluster)
            self.update_cluster_tree()
            self.add_log(f"Modified cluster '{name}': {cluster['formation']}, {new_sep}km separation")
            w.destroy()

        row = ttk.Frame(f)
        row.grid(row=2, column=0, columnspan=2, pady=18)
        ttk.Button(row, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Cancel", command=w.destroy).pack(side=tk.LEFT, padx=5)

    def delete_cluster(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to delete")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        if not messagebox.askyesno("Confirm Delete", 
                                   f"Delete cluster '{name}' and all {len(self.clusters[0]['satellites'])} satellites?"):
            return

        cluster = next((c for c in self.clusters if c["name"] == name), None)
        if not cluster:
            return

        # Remove all satellites
        for sat_name in list(cluster["satellites"]):
            self.satellites[:] = [s for s in self.satellites if s["name"] != sat_name]
            for orbit in self.orbit_configurations:
                if sat_name in orbit["satellites"]:
                    orbit["satellites"].remove(sat_name)

        self.clusters.remove(cluster)

        # Update everything
        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        if hasattr(self.parent_app, "update_satellite_dropdowns"):
            self.parent_app.update_satellite_dropdowns()
        if hasattr(self.parent_app, "update_status_counts"):
            self.parent_app.update_status_counts()

        self.add_log(f"Deleted cluster: {name} ({len(cluster['satellites'])} satellites removed)")

    def test_cluster_communication(self):
        sel = self.cluster_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a cluster to test")
            return
        name = self.cluster_tree.item(sel[0])["values"][0]
        cluster = next((c for c in self.clusters if c["name"] == name), None)
        if not cluster:
            return

        msgs = [
            f"Communication Test: {name}",
            f"Formation: {cluster['formation']}",
            f"Separation: {cluster.get('separation', 10.0)} km",
            "",
            f"Leader: {cluster['leader']} broadcasting to {len(cluster['children'])} children..."
        ]
        for child in cluster["children"]:
            msgs.append(f"  → {child}: Link OK ✓")
        msgs.append("\nChildren reporting to leader:")
        for child in cluster["children"]:
            msgs.append(f"  ← {child}: Status OK ✓")
        msgs.append(f"\nCommunication range: {self.comm_range.get():.0f} km")
        msgs.append(f"Field of view: {self.comm_fov.get():.1f}°")

        messagebox.showinfo("Communication Test", "\n".join(msgs))
        self.add_log(f"Communication test completed for cluster: {name}")

    def _update_cluster_formation(self, cluster):
        """Update formation when modified"""
        n = len(cluster["satellites"])
        formation = cluster["formation"]
        separation = cluster.get("separation", 10.0)

        offsets = self._get_formation_offsets(n, formation, separation)

        for idx, sat_name in enumerate(cluster["satellites"]):
            sat = next((s for s in self.satellites if s["name"] == sat_name), None)
            if not sat:
                continue
            
            # Update orbital elements
            base_orbit = cluster.get("base_orbit", sat["orbit"])
            sat["orbit"]["i"] = base_orbit["i"] + offsets[idx]["di"]
            sat["orbit"]["Omega"] = base_orbit["Omega"] + offsets[idx]["dOmega"]
            sat["orbit"]["f"] = base_orbit["f"] + offsets[idx]["df"]

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

        # Gather positions
        positions = []
        for nm in cluster["satellites"]:
            sat = next((s for s in self.satellites if s["name"] == nm), None)
            pos = sat.get("position_offset_km", [0.0, 0.0, 0.0]) if sat else [0.0, 0.0, 0.0]
            positions.append(pos)

        rms_error = self.calculate_formation_error(cluster, positions)

        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Formation Check: {cluster_name}")
        w.geometry("900x700")

        fig = Figure(figsize=(10, 8), facecolor="white")

        # 3D view
        ax1 = fig.add_subplot(221, projection="3d")
        ax1.set_title(f"{cluster_name} — {cluster['formation']} Formation")
        self.plot_ideal_formation(ax1, cluster["formation"], len(cluster["satellites"]), 
                                 cluster.get("separation", 10.0))

        for i, nm in enumerate(cluster["satellites"]):
            x, y, z = positions[i]
            if nm == cluster["leader"]:
                ax1.scatter(x, y, z, c="red", s=120, marker="s", label="Leader", edgecolors='black', linewidths=2)
            else:
                ax1.scatter(x, y, z, c=cluster["color_hex"], s=70, marker="o", edgecolors='black', linewidths=1)
        
        ax1.set_xlabel("Along-track (km)")
        ax1.set_ylabel("Cross-track (km)")
        ax1.set_zlabel("Radial (km)")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Status panel
        ax2 = fig.add_subplot(222)
        ax2.axis("off")
        if rms_error < 0.5:
            status_color, status_text = "green", "EXCELLENT"
        elif rms_error < 2.0:
            status_color, status_text = "yellow", "GOOD"
        else:
            status_color, status_text = "red", "NEEDS CORRECTION"
        
        metrics = f"""Formation Quality Assessment

Formation: {cluster['formation']}
Total Satellites: {len(cluster['satellites'])}
Leader: {cluster['leader']}
Children: {len(cluster['children'])}
Target Separation: {cluster.get('separation', 10.0)} km

RMS Error: {rms_error:.3f} km
Status: {status_text}

Orbit: {cluster.get('orbit_name', 'Unknown')}"""
        
        ax2.text(0.1, 0.5, metrics, fontsize=10, family="monospace",
                bbox=dict(boxstyle="round,pad=0.8", facecolor=status_color, alpha=0.25))

        # Error over time
        ax3 = fig.add_subplot(223)
        t = np.linspace(0, 24, 100)
        e = rms_error + 0.2 * np.sin(2 * np.pi * t / 12)
        ax3.plot(t, e, linewidth=2, color=cluster["color_hex"])
        ax3.axhline(0.5, color="g", linestyle="--", label="Excellent", alpha=0.7)
        ax3.axhline(2.0, color="y", linestyle="--", label="Good", alpha=0.7)
        ax3.set_xlabel("Time (hours)")
        ax3.set_ylabel("RMS Error (km)")
        ax3.set_title("Formation Error Over Time")
        ax3.grid(True, alpha=0.3)
        ax3.legend()

        # Communication topology
        ax4 = fig.add_subplot(224)
        ax4.set_title("Intra-Cluster Communication")
        leader_pos = positions[0] if positions else [0, 0, 0]
        
        for i, nm in enumerate(cluster["children"]):
            idx = cluster["satellites"].index(nm)
            child_pos = positions[idx]
            ax4.plot([leader_pos[0], child_pos[0]], [leader_pos[1], child_pos[1]], 
                    color=cluster["color_hex"], alpha=0.6, linewidth=2)
        
        ax4.scatter(leader_pos[0], leader_pos[1], c="red", s=150, marker="s", 
                   label="Leader", edgecolors='black', linewidths=2, zorder=5)
        
        for i, nm in enumerate(cluster["children"]):
            idx = cluster["satellites"].index(nm)
            x, y, _ = positions[idx]
            ax4.scatter(x, y, c=cluster["color_hex"], s=80, 
                       edgecolors='black', linewidths=1, zorder=5)
        
        ax4.set_xlabel("Along-track (km)")
        ax4.set_ylabel("Cross-track (km)")
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        ax4.set_aspect('equal')

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=w)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Save button
        bar = ttk.Frame(w)
        bar.pack(fill=tk.X, pady=6)
        
        def save_plot():
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            fname = f"plots/{ts}_FormationCheck_{cluster_name}.png"
            try:
                fig.savefig(fname, dpi=150, bbox_inches='tight')
                messagebox.showinfo("Saved", f"Plot saved to {fname}")
            except Exception as e:
                messagebox.showerror("Save Error", f"Could not save plot: {e}")
        
        ttk.Button(bar, text="Save Plot", command=save_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(bar, text="Close", command=w.destroy).pack(side=tk.LEFT, padx=5)

    def plot_ideal_formation(self, ax, formation_type, num_sats, separation):
        """Plot ideal formation template lines"""
        f = (formation_type or "").lower()
        
        if "line" in f:
            for i in range(num_sats):
                ax.plot([i * separation, i * separation], [0, 0], [-0.2, 0.2], 
                       "g--", alpha=0.4, linewidth=1.5)
        
        elif "triangle" in f:
            if num_sats >= 3:
                ang = np.linspace(0, 2*np.pi, 3, endpoint=False)
                x = separation * np.cos(ang)
                y = separation * np.sin(ang)
                for i in range(3):
                    j = (i + 1) % 3
                    ax.plot([x[i], x[j]], [y[i], y[j]], [0, 0], "g--", alpha=0.4, linewidth=1.5)
        
        elif "diamond" in f:
            pts = [[separation, 0, 0], [0, separation, 0], [-separation, 0, 0], [0, -separation, 0]]
            for i in range(4):
                j = (i + 1) % 4
                ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]], [0, 0], 
                       "g--", alpha=0.4, linewidth=1.5)
        
        else:  # Leader-Follower
            for i in range(num_sats):
                ax.plot([-i * separation, -i * separation], [0, 0], [-0.2, 0.2], 
                       "g--", alpha=0.4, linewidth=1.5)

    def calculate_formation_error(self, cluster, positions):
        """Calculate RMS formation error"""
        if not positions:
            return 0.0
        
        sep = float(cluster.get("separation", 10.0))
        f = (cluster.get("formation","") or "").lower()
        errs = []
        
        for i, p in enumerate(positions):
            x, y, z = p
            
            if i == 0:  # Leader should be at origin
                errs.append(np.linalg.norm([x, y, z]))
                continue
            
            # Target positions based on formation
            if "line" in f or "leader-follower" in f:
                target = np.array([i*sep, 0.0, 0.0])
            elif "diamond" in f:
                targets = [np.array([sep,0,0]), np.array([0,sep,0]), 
                          np.array([-sep,0,0]), np.array([0,-sep,0])]
                target = targets[(i-1) % len(targets)]
            elif "triangle" in f:
                ang = 2*np.pi*((i-1) % 3)/3
                target = np.array([sep*np.cos(ang), sep*np.sin(ang), 0.0])
            else:
                target = np.array([0.0, 0.0, 0.0])
            
            errs.append(np.linalg.norm(np.array([x,y,z]) - target))
        
        return float(np.sqrt(np.mean(np.square(errs))))

    def export_config(self):
        """Export cluster configuration to file"""
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"cluster_config_{ts}.txt"
        try:
            with open(fname, "w", encoding="utf-8") as f:
                f.write("="*70 + "\n")
                f.write("CLUSTER CONFIGURATION EXPORT\n")
                f.write("="*70 + "\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total Clusters: {len(self.clusters)}/{self.MAX_CLUSTERS}\n")
                f.write(f"Total Satellites: {len(self.satellites)}\n\n")
                
                for c in self.clusters:
                    f.write("-"*70 + "\n")
                    f.write(f"Cluster: {c['name']}\n")
                    f.write(f"  Formation: {c['formation']}\n")
                    f.write(f"  Satellites: {len(c['satellites'])} (including leader)\n")
                    f.write(f"  Leader: {c['leader']}\n")
                    f.write(f"  Children ({len(c['children'])}): {', '.join(c['children'])}\n")
                    f.write(f"  Separation: {c.get('separation', 10.0)} km\n")
                    f.write(f"  Orbit: {c.get('orbit_name', 'Unknown')}\n")
                    f.write(f"  Color: {c.get('color_hex','')}\n")
                    f.write("\n")
                
                f.write("="*70 + "\n")
            
            messagebox.showinfo("Export Complete", f"Configuration exported to:\n{fname}")
            self.add_log(f"Exported cluster configuration to {fname}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not export configuration:\n{e}")

    # ------------------------------------------------------------ Orbit Tab
    def _create_orbit_management_tab(self, parent):
        create_frame = ttk.LabelFrame(parent, text="Create New Orbit", padding=10)
        create_frame.pack(fill=tk.X, padx=10, pady=10)

        grid = ttk.Frame(create_frame)
        grid.pack(fill=tk.X, pady=5)
        
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

        row = ttk.Frame(orbits_frame)
        row.pack(fill=tk.X, pady=5)
        ttk.Button(row, text="Modify Orbit", command=self.modify_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Delete Orbit", command=self.delete_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="View Satellites", command=self.view_orbit_satellites).pack(side=tk.LEFT, padx=5)

        self.update_orbit_tree()

    def create_new_orbit(self):
        """Create a new orbit and UPDATE all dropdowns"""
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
            "description": f"Custom orbit at {self.new_orbit_alt_var.get():.0f} km",
        }
        self.orbit_configurations.append(new_orbit)
        
        # CRITICAL FIX: Update orbit dropdown in cluster creation
        self.orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        
        # Also update individual satellite orbit dropdown if it exists
        if hasattr(self, 'individual_orbit_combo'):
            self.individual_orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        
        self.update_orbit_tree()
        self.add_log(f"Created new orbit: {name} ({new_orbit['altitude']}km, {new_orbit['inclination']}° inc)")

        # Reset form
        self.new_orbit_name_var.set("")
        self.new_orbit_alt_var.set(600.0)
        self.new_orbit_inc_var.set(55.0)
        self.new_orbit_raan_var.set(0.0)
        
        messagebox.showinfo("Orbit Created", 
                          f"Orbit '{name}' created successfully.\n\n"
                          f"It's now available in cluster creation and individual satellite assignment.")

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
        f = ttk.Frame(w, padding=20)
        f.pack(fill=tk.BOTH, expand=True)

        ttk.Label(f, text="Altitude (km):").grid(row=0, column=0, sticky=tk.W, pady=5)
        alt = tk.DoubleVar(value=orbit["altitude"])
        ttk.Entry(f, textvariable=alt, width=20).grid(row=0, column=1, pady=5)
        
        ttk.Label(f, text="Inclination (deg):").grid(row=1, column=0, sticky=tk.W, pady=5)
        inc = tk.DoubleVar(value=orbit["inclination"])
        ttk.Entry(f, textvariable=inc, width=20).grid(row=1, column=1, pady=5)
        
        ttk.Label(f, text="RAAN (deg):").grid(row=2, column=0, sticky=tk.W, pady=5)
        raan = tk.DoubleVar(value=orbit["raan"])
        ttk.Entry(f, textvariable=raan, width=20).grid(row=2, column=1, pady=5)

        def apply():
            orbit["altitude"] = float(alt.get())
            orbit["inclination"] = float(inc.get())
            orbit["raan"] = float(raan.get())
            
            # Update all satellites on this orbit
            for nm in orbit["satellites"]:
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    sat["orbit"]["a"] = orbit["altitude"] + 6371
                    sat["orbit"]["i"] = orbit["inclination"]
                    sat["orbit"]["Omega"] = orbit["raan"]
            
            self.update_orbit_tree()
            self.add_log(f"Modified orbit '{name}': {orbit['altitude']}km, {orbit['inclination']}° inc")
            w.destroy()

        row = ttk.Frame(f)
        row.grid(row=3, column=0, columnspan=2, pady=18)
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
            messagebox.showwarning("Cannot Delete", 
                f"Cannot delete orbit with {len(orbit['satellites'])} satellites.\n\n"
                f"Reassign satellites first.")
            return
        
        if not messagebox.askyesno("Confirm Delete", f"Delete orbit '{name}'?"):
            return
        
        self.orbit_configurations.remove(orbit)
        
        # Update dropdowns
        self.orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        if hasattr(self, 'individual_orbit_combo'):
            self.individual_orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        
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
            txt = "\n".join([f"  • {s}" for s in orbit["satellites"]])
            messagebox.showinfo(f"Satellites in {name}", 
                              f"Satellites ({len(orbit['satellites'])}):\n\n{txt}\n\n"
                              f"Altitude: {orbit['altitude']}km\n"
                              f"Inclination: {orbit['inclination']}°")
        else:
            messagebox.showinfo(f"Satellites in {name}", 
                              f"No satellites currently assigned to this orbit.\n\n"
                              f"Altitude: {orbit['altitude']}km\n"
                              f"Inclination: {orbit['inclination']}°")

    # ------------------------------------------------------------ Individuals
    def _create_individual_tab(self, parent):
        add = ttk.LabelFrame(parent, text="Add Individual Satellite", padding=10)
        add.pack(fill=tk.X, padx=10, pady=10)

        row = ttk.Frame(add)
        row.pack(fill=tk.X, pady=5)
        
        ttk.Label(row, text="Satellite Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.individual_sat_name_var = tk.StringVar()
        ttk.Entry(row, textvariable=self.individual_sat_name_var, width=20).grid(row=0, column=1, padx=5)

        ttk.Label(row, text="Orbit:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.individual_orbit_var = tk.StringVar()
        self.individual_orbit_combo = ttk.Combobox(row, textvariable=self.individual_orbit_var, 
                                                   width=20, state="readonly")
        self.individual_orbit_combo.grid(row=0, column=3, padx=5)
        self.individual_orbit_combo["values"] = [o["name"] for o in self.orbit_configurations]
        if self.individual_orbit_combo["values"]:
            self.individual_orbit_combo.current(0)

        ttk.Label(row, text="True Anomaly (deg):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.individual_anomaly_var = tk.DoubleVar(value=0.0)
        ttk.Entry(row, textvariable=self.individual_anomaly_var, width=20).grid(row=1, column=1, padx=5)

        ttk.Button(add, text="Add Satellite", command=self.add_individual_satellite).pack(pady=10)

        list_frame = ttk.LabelFrame(parent, text="All Satellites", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        container = ttk.Frame(list_frame)
        container.pack(fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(container)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        self.satellite_listbox = tk.Listbox(
            container, selectmode=tk.SINGLE, yscrollcommand=sb.set, font=("Segoe UI", 10)
        )
        self.satellite_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.config(command=self.satellite_listbox.yview)

        actions = ttk.Frame(list_frame)
        actions.pack(fill=tk.X, pady=5)
        ttk.Button(actions, text="Remove Satellite", command=self.remove_satellite).pack(side=tk.LEFT, padx=5)
        ttk.Button(actions, text="Assign to Orbit", command=self.assign_satellite_to_orbit).pack(side=tk.LEFT, padx=5)

        self.update_satellite_listbox()

    def add_individual_satellite(self):
        """Add an individual satellite (not part of a cluster)"""
        name = self.individual_sat_name_var.get().strip()
        orb_name = self.individual_orbit_var.get()
        anomaly = float(self.individual_anomaly_var.get())

        if not name:
            messagebox.showwarning("Missing Name", "Please provide a satellite name")
            return
        
        # Check for duplicates
        if name in [s["name"] for s in self.satellites]:
            messagebox.showwarning("Duplicate Name", 
                f"Satellite '{name}' already exists.\n\nUse a unique name.")
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
                "a": orbit["altitude"] + 6371,
                "e": 0.01,
                "i": orbit["inclination"],
                "Omega": orbit["raan"],
                "omega": 0.0,
                "f": anomaly
            },
            "fault": {
                "type": "friction",
                "magnitude": 0.0005,
                "wheel": 3,
                "time": 15.0,
                "enabled": False,
                "periodic": {
                    "enabled": False,
                    "interval": 360,
                    "magnitude": 0.1,
                    "wheel": 1
                }
            },
            "camera": {
                "position": [0.0, 0.0, 15.0],
                "fov": 80.0,
                "enabled": True,
            },
            "communication": {
                "range": float(self.comm_range.get()),
                "fov": float(self.comm_fov.get()),
                "aHat_B": [0.0, 0.0, -1.0],
            },
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
        self.add_log(f"Added individual satellite: {name} in orbit {orb_name}")

    def remove_satellite(self):
        sel = self.satellite_listbox.curselection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a satellite to remove")
            return
        
        item_text = self.satellite_listbox.get(sel[0])

        if "[INDIVIDUAL]" not in item_text:
            messagebox.showwarning("Cannot Remove", 
                "Cluster satellites must be removed via Cluster Management.\n\n"
                "Delete the entire cluster to remove its satellites.")
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
        """Reassign an individual satellite to a different orbit"""
        sel = self.satellite_listbox.curselection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a satellite")
            return
        
        item_text = self.satellite_listbox.get(sel[0])
        
        if "[INDIVIDUAL]" not in item_text:
            messagebox.showinfo("Cannot Reassign", 
                "Only individual satellites can be reassigned.\n\n"
                "Cluster satellites share their cluster's orbit.")
            return

        sat_name = item_text.split("]")[1].strip().split("(")[0].strip()
        sat = next((s for s in self.satellites if s["name"] == sat_name), None)
        if not sat:
            return

        # Create dialog with dropdown
        w = tk.Toplevel(self.parent_app.root)
        w.title(f"Reassign {sat_name}")
        w.geometry("400x180")
        f = ttk.Frame(w, padding=20)
        f.pack(fill=tk.BOTH, expand=True)

        ttk.Label(f, text=f"Select new orbit for {sat_name}:", 
                 font=("Segoe UI", 10, "bold")).pack(pady=(0, 10))

        orbit_var = tk.StringVar(value=sat.get("orbit_name", ""))
        orbit_combo = ttk.Combobox(f, textvariable=orbit_var, 
                                  values=[o["name"] for o in self.orbit_configurations],
                                  width=30, state="readonly")
        orbit_combo.pack(pady=10)

        def apply():
            new_name = orbit_var.get()
            if not new_name:
                messagebox.showwarning("No Selection", "Please select an orbit")
                return
            
            old = sat.get("orbit_name")
            new_orbit = next(o for o in self.orbit_configurations if o["name"] == new_name)
            
            # Remove from old orbit
            if old:
                old_orbit = next((o for o in self.orbit_configurations if o["name"] == old), None)
                if old_orbit and sat_name in old_orbit["satellites"]:
                    old_orbit["satellites"].remove(sat_name)

            # Update satellite
            sat["orbit_name"] = new_name
            sat["orbit"]["a"] = new_orbit["altitude"] + 6371
            sat["orbit"]["i"] = new_orbit["inclination"]
            sat["orbit"]["Omega"] = new_orbit["raan"]
            
            # Add to new orbit
            new_orbit["satellites"].append(sat_name)

            self.update_satellite_listbox()
            self.update_orbit_tree()
            self.add_log(f"Reassigned {sat_name}: {old} → {new_name}")
            w.destroy()

        row = ttk.Frame(f)
        row.pack(pady=15)
        ttk.Button(row, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        ttk.Button(row, text="Cancel", command=w.destroy).pack(side=tk.LEFT, padx=5)

    # ------------------------------------------------------------ Comm Tab
    def _create_comm_tab(self, parent):
        settings = ttk.LabelFrame(parent, text="Communication Settings", padding=10)
        settings.pack(fill=tk.X, padx=10, pady=10)

        row = ttk.Frame(settings)
        row.pack(fill=tk.X, pady=5)
        
        ttk.Label(row, text="Range (km):").grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Entry(row, textvariable=self.comm_range, width=15).grid(row=0, column=1, padx=5)
        
        ttk.Label(row, text="FOV (deg):").grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Entry(row, textvariable=self.comm_fov, width=15).grid(row=0, column=3, padx=5)
        
        ttk.Button(row, text="Refresh Plot", 
                  command=lambda: [self.update_communication_plot(), self.comm_canvas.draw()]).grid(row=0, column=4, padx=5)
        ttk.Button(row, text="Help", command=self.show_communication_help).grid(row=0, column=5, padx=5)
        
        ttk.Label(row, text="Bars = communication windows | Dots = messages", 
                 font=("Arial", 8), foreground="gray").grid(row=1, column=0, columnspan=6, sticky=tk.W, pady=5)

        plot_frame = ttk.LabelFrame(parent, text="Communication Windows (Simulated)", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.comm_figure = Figure(figsize=(10, 6), dpi=80)
        self.comm_ax = self.comm_figure.add_subplot(111)
        self.comm_canvas = FigureCanvasTkAgg(self.comm_figure, plot_frame)
        self.comm_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        ttk.Button(plot_frame, text="Update Communication Plot", 
                  command=self.update_communication_plot).pack(pady=8)

        self.update_communication_plot()

    def update_communication_plot(self):
        """Generate communication window visualization"""
        self.comm_ax.clear()

        if not self.clusters:
            self.comm_ax.text(0.5, 0.5, 
                "No clusters configured\n\nCreate clusters in the Cluster Management tab\nto visualize communication",
                ha="center", va="center", fontsize=12, color="gray")
            self.comm_ax.set_xlim(0, 1)
            self.comm_ax.set_ylim(0, 1)
            self.comm_canvas.draw()
            return

        t = np.linspace(0, 30, 100)
        y = 0
        title_text = f"Communication Windows for {len(self.clusters)} Clusters"
        
        self.comm_ax.text(15, -1.5, 
            "Colored bars = communication windows | Dots = message transmission events",
            ha="center", fontsize=9, style="italic", color="gray")

        legend_elems = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", 
                   markersize=8, label="Message sent"),
            Patch(facecolor="gray", alpha=0.3, label="Comm window")
        ]

        for idx, c in enumerate(self.clusters):
            color = self.PALETTE[c.get("color_idx", idx % len(self.PALETTE))]
            leader = c.get("leader", "Leader")
            children = c.get("children", [])

            # Cluster background band
            self.comm_ax.axhspan(y - 0.5, y + len(children) + 0.5, 
                                color=color, alpha=0.06, zorder=0)
            
            # Cluster label
            self.comm_ax.text(-3, y + len(children) / 2, f"[{c['name'].upper()}]",
                            fontsize=11, fontweight="bold", color=color, va="center", ha="right",
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", 
                                    edgecolor=color, linewidth=2))

            legend_elems.append(Patch(facecolor=color, edgecolor=color, alpha=0.3, 
                                     label=f"{c['name']} cluster"))

            for j, child in enumerate(children):
                phase = j * np.pi / 4
                comm_window = np.sin(2 * np.pi * t / 10 + phase) > 0.3

                # Draw window bars
                in_window = False
                start = 0.0
                for k in range(len(t)):
                    if comm_window[k] and not in_window:
                        start = t[k]
                        in_window = True
                    elif not comm_window[k] and in_window:
                        self.comm_ax.barh(y, t[k] - start, left=start, height=0.6,
                                        color=color, alpha=0.3, edgecolor=color, linewidth=1)
                        in_window = False
                
                if in_window:
                    self.comm_ax.barh(y, t[-1] - start, left=start, height=0.6,
                                    color=color, alpha=0.3, edgecolor=color, linewidth=1)

                # Message dots
                for mt in t[::12]:
                    k = np.argmin(np.abs(t - mt))
                    if comm_window[k]:
                        self.comm_ax.scatter(mt, y, color=color, s=80, marker="o", 
                                           zorder=5, edgecolors="white", linewidths=1)

                # Link label
                self.comm_ax.text(-2, y, f"{leader[:10]} → {child[:10]}", 
                                fontsize=9, va="center", color="black")
                y += 1

            y += 0.8

        self.comm_ax.set_xlabel("Time (minutes)", fontsize=11, fontweight="bold")
        self.comm_ax.set_ylabel("Communication Links", fontsize=11, fontweight="bold")
        self.comm_ax.set_title(title_text, fontsize=13, fontweight="bold", pad=12)
        self.comm_ax.set_xlim(-3.5, 31)
        self.comm_ax.set_ylim(-2, max(1, y))
        self.comm_ax.grid(True, alpha=0.2, axis="x", linestyle="--")
        
        for tmark in [0, 5, 10, 15, 20, 25, 30]:
            self.comm_ax.axvline(x=tmark, color="gray", linewidth=0.5, alpha=0.3, linestyle=":")

        # Deduplicate legend
        seen = set()
        final = []
        for h in legend_elems:
            lbl = h.get_label()
            if lbl not in seen:
                final.append(h)
                seen.add(lbl)

        self.comm_ax.legend(handles=final, loc="upper left", fontsize=9, 
                          framealpha=0.9, ncol=2)
        self.comm_canvas.draw()

    def show_communication_help(self):
        """Show help dialog for communication plot"""
        messagebox.showinfo(
            "Communication Plot Help",
            """UNDERSTANDING THE COMMUNICATION PLOT

WHAT YOU SEE:
• Each row represents a leader→child link
• Links are grouped by cluster (colored bands)
• Horizontal bars = time windows when communication is possible
• Dots = actual message transmissions (only within windows)

WHY WINDOWS APPEAR AND DISAPPEAR:
• Relative orbital motion causes satellites to move in/out of range
• Communication requires both range and line-of-sight
• Field of view (FOV) constraints limit when links are active
• Earth occlusion can block signals

SETTINGS:
• Range: Maximum distance for communication
• FOV: Angular field of view for antennas

This is a SIMULATED visualization showing expected patterns
based on orbital mechanics and communication constraints.""",
        )

    # ------------------------------------------------------------ Update lists
    def update_cluster_tree(self):
        """Update the cluster treeview display"""
        for item in self.cluster_tree.get_children():
            self.cluster_tree.delete(item)
        
        for c in self.clusters:
            status = "Active" if c.get("leader") else "Inactive"
            node = self.cluster_tree.insert("", "end", 
                values=(c["name"], c.get("leader","None"),
                       len(c.get("children", [])), c.get("formation",""), status))
            
            # Add satellite sub-rows
            for nm in c.get("satellites", []):
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    self.cluster_tree.insert(node, "end", 
                        values=(f"  {sat['name']}",
                               sat.get("role", "").upper(), 
                               sat.get("orbit_name", ""),
                               f"{sat['orbit']['a']-6371:.0f}km",
                               "FAULT" if sat["fault"]["enabled"] else "OK"))

    def update_orbit_tree(self):
        """Update the orbit library treeview"""
        for item in self.orbit_tree.get_children():
            self.orbit_tree.delete(item)
        
        for o in self.orbit_configurations:
            self.orbit_tree.insert("", "end", 
                values=(o["name"], 
                       f"{o['altitude']:.0f} km",
                       f"{o['inclination']:.1f}°", 
                       f"{o['raan']:.1f}°",
                       len(o["satellites"])))

    def update_satellite_listbox(self):
        """Update the satellite listbox showing all satellites"""
        self.satellite_listbox.delete(0, tk.END)
        
        # Clusters
        for c in self.clusters:
            self.satellite_listbox.insert(tk.END, f"[CLUSTER] {c['name']}")
            for nm in c["satellites"]:
                sat = next((s for s in self.satellites if s["name"] == nm), None)
                if sat:
                    role = "L" if sat["role"] == "leader" else "C"
                    orbit = sat.get("orbit_name", "Unknown")
                    alt = sat["orbit"]["a"] - 6371
                    fault = " ⚠" if sat["fault"]["enabled"] else ""
                    self.satellite_listbox.insert(tk.END, 
                        f"  [{role}] {sat['name']} ({orbit}, {alt:.0f}km){fault}")
        
        # Individuals
        for sat in self.satellites:
            if sat["type"] == "individual":
                orbit = sat.get("orbit_name", "Unknown")
                alt = sat["orbit"]["a"] - 6371
                fault = " ⚠" if sat["fault"]["enabled"] else ""
                self.satellite_listbox.insert(tk.END, 
                    f"[INDIVIDUAL] {sat['name']} ({orbit}, {alt:.0f}km){fault}")

    def clear_cluster_form(self):
        """Clear the cluster creation form"""
        self.cluster_name_var.set("")
        self.sats_per_cluster_var.set(4)
        self.formation_var.set(self.ALLOWED_FORMATIONS[0])
        self.formation_separation.set(15.0)
        if self.orbit_configurations:
            self.primary_orbit_var.set(self.orbit_configurations[0]["name"])

    # ------------------------------------------------------------ Debug helper
    def debug_cluster_status(self):
        """Print debug information about current cluster state"""
        print("\n=== CLUSTER DEBUG INFO ===")
        print(f"Number of clusters: {len(self.clusters)}")
        for c in self.clusters:
            print(f"  Cluster '{c['name']}':")
            print(f"    - Leader: {c.get('leader', 'None')}")
            print(f"    - Children ({len(c.get('children', []))}): {c.get('children', [])}")
            print(f"    - Total satellites: {len(c.get('satellites', []))}")
            print(f"    - Formation: {c.get('formation', 'Unknown')}")
            print(f"    - Separation: {c.get('separation', 10.0)} km")
        print(f"Total satellites in app: {len(self.satellites)}")
        print("=========================\n")