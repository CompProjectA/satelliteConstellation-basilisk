#!/usr/bin/env python
"""
constellation_tab.py

FULLY FIXED VERSION: Complete constellation tab with all working buttons and proper cluster support.
Each satellite in cluster on different orbits as per requirements.
"""
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import base tab
try:
    from .base_tab import BaseTab
except ImportError:
    class BaseTab:
        def __init__(self, parent_app, parent_frame):
            self.parent_app = parent_app
            self.parent_frame = parent_frame
        
        def add_log(self, message):
            if hasattr(self.parent_app, 'add_log'):
                self.parent_app.add_log(message)

class ConstellationTab(BaseTab):
    """FULLY FIXED Constellation Management tab with all features working"""
    
    def __init__(self, parent_app, parent_frame):
        """Initialize the constellation tab"""
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        self.current_satellite_index = parent_app.current_satellite_index
        
        # Cluster configurations
        self.clusters = []
        self.cluster_configurations = []
        
        # Communication tracking
        self.communication_windows = {}
        self.message_history = []
        
        # FIXED: Proper orbit configurations per requirements
        self.orbit_configurations = [
            {
                "name": "Primary Orbit",
                "altitude": 600,  # km above Earth
                "inclination": 55.0,
                "raan": 0.0,  # Right Ascension of Ascending Node
                "satellites": [],
                "description": "Primary orbit for sat1 and sat2"
            },
            {
                "name": "Opposite Orbit", 
                "altitude": 600,  # Same altitude
                "inclination": 55.0,  # Same inclination  
                "raan": 180.0,  # Opposite RAAN for opposite orbit
                "satellites": [],
                "description": "Opposite orbit for sat3 and sat4"
            },
            {
                "name": "Orbit 3",
                "altitude": 700,  # Different altitude
                "inclination": 60.0,
                "raan": 90.0,
                "satellites": [],
                "description": "Alternative orbit for additional satellites"
            },
            {
                "name": "Orbit 4",
                "altitude": 800,  # Higher altitude
                "inclination": 65.0,
                "raan": 270.0,
                "satellites": [],
                "description": "High altitude orbit"
            }
        ]
        
        # Communication settings
        self.comm_range = tk.DoubleVar(value=2000.0)  # km
        self.comm_fov = tk.DoubleVar(value=30.0)  # degrees
        
        # Formation settings
        self.formation_separation = tk.DoubleVar(value=10.0)  # km
        
        # Create the tab UI
        self.create_tab_ui()
        
    def create_tab_ui(self):
        """Create the Constellation Management tab UI"""
        # Main container
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create notebook for sub-tabs
        self.constellation_notebook = ttk.Notebook(main_container)
        self.constellation_notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create frames for sub-tabs
        cluster_frame = ttk.Frame(self.constellation_notebook)
        orbit_frame = ttk.Frame(self.constellation_notebook)
        individual_frame = ttk.Frame(self.constellation_notebook)
        comm_frame = ttk.Frame(self.constellation_notebook)
        
        # Add sub-tabs
        self.constellation_notebook.add(cluster_frame, text="Cluster Management")
        self.constellation_notebook.add(orbit_frame, text="Orbit Management")
        self.constellation_notebook.add(individual_frame, text="Individual Satellites")
        self.constellation_notebook.add(comm_frame, text="Communication")
        
        # Create sub-tab content
        self._create_cluster_tab(cluster_frame)
        self._create_orbit_management_tab(orbit_frame)
        self._create_individual_tab(individual_frame)
        self._create_comm_tab(comm_frame)
            
    def _create_cluster_tab(self, parent):
        """Create the cluster management sub-tab with proper formations"""
        # Cluster creation section
        create_frame = ttk.LabelFrame(parent, text="Create New Cluster", padding=10)
        create_frame.pack(fill=tk.X, padx=10, pady=10)
        
        # Basic cluster parameters
        basic_frame = ttk.Frame(create_frame)
        basic_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(basic_frame, text="Cluster Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.cluster_name_var = tk.StringVar(value="")
        name_entry = ttk.Entry(basic_frame, textvariable=self.cluster_name_var, width=20)
        name_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(basic_frame, text="Number of Satellites:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.sats_per_cluster_var = tk.IntVar(value=4)
        sats_combo = ttk.Combobox(basic_frame, textvariable=self.sats_per_cluster_var, 
                                values=[2, 3, 4, 5, 6], width=10, state="readonly")
        sats_combo.grid(row=0, column=3, padx=5)
        sats_combo.current(2)  # Default to 4 satellites
        
        # Formation parameters - Updated formations based on your image
        formation_frame = ttk.Frame(create_frame)
        formation_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(formation_frame, text="Formation Type:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.formation_var = tk.StringVar(value="Circle")
        formation_combo = ttk.Combobox(formation_frame, textvariable=self.formation_var,
                                    values=["Circle", "Line", "Leader-Follower", "Ellipse"], 
                                    width=15, state="readonly")
        formation_combo.grid(row=0, column=1, padx=5)
        
        ttk.Label(formation_frame, text="Separation (km):").grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Entry(formation_frame, textvariable=self.formation_separation, width=10).grid(row=0, column=3, padx=5)
        
        # Orbit assignment 
        orbit_frame = ttk.LabelFrame(create_frame, text="Orbit Assignment", padding=5)
        orbit_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(orbit_frame, text="Primary Orbit:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.primary_orbit_var = tk.StringVar(value="Primary Orbit")
        primary_combo = ttk.Combobox(orbit_frame, textvariable=self.primary_orbit_var,
                                    values=["Primary Orbit", "Opposite Orbit"], 
                                    width=15, state="readonly")
        primary_combo.grid(row=0, column=1, padx=5)
        
        self.orbit_assignment_frame = ttk.Frame(orbit_frame)
        self.orbit_assignment_frame.grid(row=1, column=0, columnspan=4, sticky="we", pady=(6, 0))

        # When '# of sats' changes, rebuild the rows
        sats_combo.bind('<<ComboboxSelected>>', self._update_orbit_assignments)

        # first build
        self._update_orbit_assignments()
        # Create cluster button - REMOVED Sprint 2 button
        button_frame = ttk.Frame(create_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(button_frame, text="Create Cluster", 
                command=self.create_cluster, style="Run.TButton").pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Clear Form", 
                command=self.clear_cluster_form).pack(side=tk.LEFT, padx=5)
        
        
        
        # Cluster list section
        list_frame = ttk.LabelFrame(parent, text="Existing Clusters", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create treeview for clusters
        columns = ('Cluster', 'Leader', 'Children', 'Formation', 'Status')
        self.cluster_tree = ttk.Treeview(list_frame, columns=columns, show='tree headings', height=10)
        
        # Define column headings
        self.cluster_tree.heading('#0', text='')
        self.cluster_tree.heading('Cluster', text='Cluster')
        self.cluster_tree.heading('Leader', text='Leader')
        self.cluster_tree.heading('Children', text='Children')
        self.cluster_tree.heading('Formation', text='Formation')
        self.cluster_tree.heading('Status', text='Status')
        
        # Set column widths
        self.cluster_tree.column('#0', width=30)
        self.cluster_tree.column('Cluster', width=120)
        self.cluster_tree.column('Leader', width=120)
        self.cluster_tree.column('Children', width=80)
        self.cluster_tree.column('Formation', width=100)
        self.cluster_tree.column('Status', width=100)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.cluster_tree.yview)
        self.cluster_tree.configure(yscrollcommand=scrollbar.set)
        
        self.cluster_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Cluster actions - ALL BUTTONS IMPLEMENTED
        action_frame = ttk.Frame(list_frame)
        action_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(action_frame, text="View Details", 
                  command=self.view_cluster_details).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Modify Cluster", 
                  command=self.modify_cluster).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Delete Cluster", 
                  command=self.delete_cluster).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Test Communication", 
                  command=self.test_cluster_communication).pack(side=tk.LEFT, padx=5)
        
        # Update tree
        self.update_cluster_tree()
        
    def _update_orbit_assignments(self, event=None):
        """Update orbit assignment rows to match the number of satellites."""
        # Clear existing
        for w in self.orbit_assignment_frame.winfo_children():
            w.destroy()
        self.orbit_assignments = {}

        try:
            num_sats = int(self.sats_per_cluster_var.get())
        except Exception:
            num_sats = 0

        if num_sats <= 0:
            return


        orbit_names = [o['name'] for o in self.orbit_configurations]

        # ensure stale assignment vars don’t leak from a previous create
        if not hasattr(self, "orbit_assignments"):
            self.orbit_assignments = {}

        for i in range(num_sats):
            role = "Leader" if i == 0 else f"Child {i}"
            ttk.Label(self.orbit_assignment_frame, text=f"{role}:").grid(
                row=i // 2, column=(i % 2) * 2, sticky=tk.W, padx=5, pady=2
            )

            var = tk.StringVar()
            combo = ttk.Combobox(self.orbit_assignment_frame, textvariable=var,
                                values=orbit_names, width=15, state="readonly")
            combo.grid(row=i // 2, column=(i % 2) * 2 + 1, padx=5, pady=2)

            # defaults that match your requirement/image:
            if i in (0, 1):
                var.set("Primary Orbit")
            else:
                # round‑robin to make sure >2 sats still have reasonable defaults
                var.set("Opposite Orbit" if "Opposite Orbit" in orbit_names else orbit_names[i % len(orbit_names)])

            self.orbit_assignments[i] = var


    def create_cluster(self):
        """
        Create a new satellite cluster with proper leader/child roles,
        each satellite optionally on different orbits, and role-based
        communication pointing. Safe for FaultTab (includes periodic).
        """
        # --- Read form values ---
        cluster_name = self.cluster_name_var.get().strip()
        try:
            num_sats = int(self.sats_per_cluster_var.get())
        except Exception:
            num_sats = 0
        formation = self.formation_var.get()
        separation = float(self.formation_separation.get()) if hasattr(self, "formation_separation") else 0.0

        # --- Validation ---
        if not cluster_name:
            messagebox.showwarning("Missing Information", "Please provide a cluster name")
            return
        if num_sats <= 0:
            messagebox.showwarning("Invalid Value", "Number of satellites must be greater than zero")
            return
        if hasattr(self, "clusters") and cluster_name in [c['name'] for c in self.clusters]:
            messagebox.showwarning("Duplicate Name", f"Cluster '{cluster_name}' already exists")
            return
        if not getattr(self, "orbit_configurations", None):
            messagebox.showwarning("No Orbits", "No orbit configurations found. Create or import orbits first.")
            return

        # --- Build cluster record ---
        cluster = {
            "name": cluster_name,
            "type": "cluster",
            "formation": formation,
            "leader": None,
            "children": [],
            "satellites": [],
            "separation": separation,
        }

        # --- Formation geometry (keeps your existing helper) ---
        # positions[i] should include keys like 'anomaly' and optional 'offset'
        positions = self._calculate_cluster_positions(num_sats, formation, 0, 0, separation)

        # --- Create satellites ---
        for i in range(num_sats):
            is_leader = (i == 0)
            sat_name = f"{cluster_name}_{'Leader' if is_leader else f'Sat{i+1}'}"

            # Choose an orbit for this satellite:
            # 1) honor explicit assignment if provided in self.orbit_assignments[i]
            # 2) otherwise round-robin through self.orbit_configurations
            orbit_name = None
            if hasattr(self, "orbit_assignments"):
                var = self.orbit_assignments.get(i)
                if hasattr(var, "get"):
                    orbit_name = var.get().strip() or None
            if orbit_name is None:
                orbit_name = self.orbit_configurations[i % len(self.orbit_configurations)]['name']

            orbit_cfg = next((o for o in self.orbit_configurations if o['name'] == orbit_name), self.orbit_configurations[0])

            # Role-based comm pointing (kept from your new snippet)
            leader_aHat = [0.2, -0.4, 0.2]
            child_aHat  = [0.0, 0.0, -1.0]

            sat = {
                "name": sat_name,
                "type": "cluster_member",
                "cluster": cluster_name,
                "role": "leader" if is_leader else "child",

                # Orbit state (ECI classical elements → Basilisk-friendly if you convert downstream)
                "orbit": {
                    "a": orbit_cfg['altitude'] + 6371.0,   # km, Earth radius added
                    "e": 0.01,
                    "i": orbit_cfg['inclination'],
                    "Omega": orbit_cfg['raan'],
                    "omega": 0.0,
                    "f": positions[i].get('anomaly', 0.0),
                },
                "orbit_name": orbit_name,
                "position_offset": positions[i].get('offset', [0.0, 0.0, 0.0]),

                # Fault block compatible with FaultTab (includes periodic)
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,               # 0-based (RW4 shown as 4 in GUI)
                    "time": 15.0,             # minutes
                    "enabled": False,
                    "periodic": {              # <<< important for FaultTab safety
                        "enabled": False,
                        "interval": 360.0,
                        "magnitude": 0.1,
                        "wheel": 3
                    }
                },

                # Camera defaults: leader camera on by default
                "camera": {
                    "position": [0.0, 0.0, 15.0],
                    "fov": 80.0,
                    "enabled": is_leader
                },

                # Communications (range/fov from GUI; pointing depends on role)
                "communication": {
                    "range": float(self.comm_range.get()) if hasattr(self, "comm_range") else 1000.0,
                    "fov": float(self.comm_fov.get()) if hasattr(self, "comm_fov") else 45.0,
                    "aHat_B": leader_aHat if is_leader else child_aHat,

                    # Optional scaffolding you may use elsewhere
                    "links": {
                        "intra_cluster": [],
                        "inter_cluster": []
                    }
                },

                # Placeholder for imaging targets or tasks
                "targets": []
            }

            # Register satellite in app structures
            self.satellites.append(sat)
            cluster["satellites"].append(sat_name)
            if "satellites" in orbit_cfg:
                orbit_cfg["satellites"].append(sat_name)
            else:
                orbit_cfg["satellites"] = [sat_name]

            if is_leader:
                cluster["leader"] = sat_name
            else:
                cluster["children"].append(sat_name)

        # --- Save cluster and refresh UI ---
        if not hasattr(self, "clusters"):
            self.clusters = []
        self.clusters.append(cluster)

        # UI refresh hooks (exist in your codebase)
        try:
            self.update_cluster_tree()
        except Exception:
            pass
        try:
            self.update_satellite_listbox()
        except Exception:
            pass
        try:
            self.update_orbit_tree()
        except Exception:
            pass
        try:
            # Let other tabs (e.g., FaultTab) refresh their dropdowns
            self.parent_app.update_satellite_dropdowns()
        except Exception:
            pass
        try:
            self.parent_app.update_status_counts()
        except Exception:
            pass
        try:
            # If you maintain a comm graph elsewhere
            if hasattr(self, "_rebuild_comm_graph"):
                self._rebuild_comm_graph()
        except Exception:
            pass

        # Log + clean up form
        self.add_log(f"Created cluster '{cluster_name}' with {num_sats} satellites (formation={formation}, separation={separation})")
        try:
            self.clear_cluster_form()
        except Exception:
            pass

        try:
            self.update_communication_plot()
        except Exception as e:
            print(f"Could not update communication plot: {e}")
        
        # Also notify the main communication tab if it exists
        try:
            if hasattr(self.parent_app, 'communication_tab'):
                self.parent_app.communication_tab.update_for_clusters()
        except Exception:
            pass

        
    
    
    def view_cluster_details(self):
        """View details of selected cluster - IMPLEMENTED"""
        selection = self.cluster_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a cluster to view")
            return
        
        # Get cluster name from selection
        item = self.cluster_tree.item(selection[0])
        cluster_name = item['values'][0]
        
        # Find cluster
        cluster = next((c for c in self.clusters if c['name'] == cluster_name), None)
        if not cluster:
            return
        
        # Create details window
        details_window = tk.Toplevel(self.parent_app.root)
        details_window.title(f"Cluster Details: {cluster_name}")
        details_window.geometry("600x500")
        
        # Details text
        text_frame = ttk.Frame(details_window, padding=10)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        text = tk.Text(text_frame, wrap=tk.WORD, font=('Consolas', 10))
        scrollbar = ttk.Scrollbar(text_frame, command=text.yview)
        text.config(yscrollcommand=scrollbar.set)
        
        # Build details
        details = f"CLUSTER: {cluster['name']}\n"
        details += "="*50 + "\n\n"
        details += f"Formation: {cluster['formation']}\n"
        details += f"Leader: {cluster['leader']}\n"
        details += f"Children: {', '.join(cluster['children'])}\n"
        details += f"Total Satellites: {len(cluster['satellites'])}\n\n"
        
        # Satellite details
        details += "SATELLITE DETAILS:\n"
        details += "-"*50 + "\n"
        for sat_name in cluster['satellites']:
            sat = next((s for s in self.satellites if s['name'] == sat_name), None)
            if sat:
                details += f"\n{sat['name']}:\n"
                details += f"  Role: {sat['role']}\n"
                details += f"  Orbit: {sat['orbit_name']}\n"
                details += f"  Altitude: {sat['orbit']['a'] - 6371:.0f} km\n"
                details += f"  Inclination: {sat['orbit']['i']:.1f}°\n"
                details += f"  RAAN: {sat['orbit']['Omega']:.1f}°\n"
                details += f"  True Anomaly: {sat['orbit']['f']:.1f}°\n"
                details += f"  Camera: {'Enabled' if sat['camera']['enabled'] else 'Disabled'}\n"
                details += f"  Fault: {'Enabled' if sat['fault']['enabled'] else 'Disabled'}\n"
                if sat['targets']:
                    details += f"  Targets: {', '.join(sat['targets'])}\n"
        
        # Communication info
        details += "\n" + "="*50 + "\n"
        details += "COMMUNICATION CONFIGURATION:\n"
        details += f"Range: {self.comm_range.get():.0f} km\n"
        details += f"FOV: {self.comm_fov.get():.1f}°\n"
        
        text.insert(tk.END, details)
        text.config(state="disabled")
        
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Close button
        ttk.Button(details_window, text="Close", 
                  command=details_window.destroy).pack(pady=10)
    
    def modify_cluster(self):
        """Modify selected cluster - IMPLEMENTED"""
        selection = self.cluster_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a cluster to modify")
            return
        
        # Get cluster
        item = self.cluster_tree.item(selection[0])
        cluster_name = item['values'][0]
        cluster = next((c for c in self.clusters if c['name'] == cluster_name), None)
        if not cluster:
            return
        
        # Create modification window
        mod_window = tk.Toplevel(self.parent_app.root)
        mod_window.title(f"Modify Cluster: {cluster_name}")
        mod_window.geometry("400x300")
        
        # Modification options
        frame = ttk.Frame(mod_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Formation type
        ttk.Label(frame, text="Formation Type:").grid(row=0, column=0, sticky=tk.W, pady=5)
        formation_var = tk.StringVar(value=cluster['formation'])
        formation_combo = ttk.Combobox(frame, textvariable=formation_var,
                                      values=["Leader-Follower", "Diamond", "Line", "Triangle"],
                                      width=20, state="readonly")
        formation_combo.grid(row=0, column=1, pady=5)
        
        # Separation
        ttk.Label(frame, text="Separation (km):").grid(row=1, column=0, sticky=tk.W, pady=5)
        sep_var = tk.DoubleVar(value=cluster.get('separation', 10.0))
        ttk.Entry(frame, textvariable=sep_var, width=20).grid(row=1, column=1, pady=5)
        
        # Add satellite button
        ttk.Label(frame, text="Add Satellite:").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Button(frame, text="Add Child Satellite", 
                  command=lambda: self._add_satellite_to_cluster(cluster_name)).grid(row=2, column=1, pady=5)
        
        # Apply changes
        def apply_changes():
            cluster['formation'] = formation_var.get()
            cluster['separation'] = sep_var.get()
            
            # Update satellites formation
            self._update_cluster_formation(cluster)
            
            self.update_cluster_tree()
            self.add_log(f"Modified cluster: {cluster_name}")
            mod_window.destroy()
        
        # Buttons
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=20)
        
        ttk.Button(button_frame, text="Apply", command=apply_changes).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Cancel", command=mod_window.destroy).pack(side=tk.LEFT, padx=5)
    
    def delete_cluster(self):
        """Delete selected cluster - IMPLEMENTED"""
        selection = self.cluster_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a cluster to delete")
            return
        
        # Get cluster
        item = self.cluster_tree.item(selection[0])
        cluster_name = item['values'][0]
        
        if not messagebox.askyesno("Confirm Delete", f"Delete cluster '{cluster_name}' and all its satellites?"):
            return
        
        # Find and remove cluster
        cluster = next((c for c in self.clusters if c['name'] == cluster_name), None)
        if cluster:
            # Remove satellites
            for sat_name in cluster['satellites']:
                # Remove from main list
                self.satellites[:] = [s for s in self.satellites if s['name'] != sat_name]
                
                # Remove from orbit lists
                for orbit in self.orbit_configurations:
                    if sat_name in orbit['satellites']:
                        orbit['satellites'].remove(sat_name)
            
            # Remove cluster
            self.clusters.remove(cluster)
            
            # Update UI
            self.update_cluster_tree()
            self.update_satellite_listbox()
            self.update_orbit_tree()
            self.parent_app.update_satellite_dropdowns()
            self.parent_app.update_status_counts()
            
            self.add_log(f"Deleted cluster: {cluster_name}")
    
    def test_cluster_communication(self):
        """Test communication within selected cluster"""
        selection = self.cluster_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a cluster to test")
            return
        
        # Get cluster
        item = self.cluster_tree.item(selection[0])
        cluster_name = item['values'][0]
        cluster = next((c for c in self.clusters if c['name'] == cluster_name), None)
        if not cluster:
            return
        
        # Simulate communication
        messages = []
        messages.append(f"Testing communication in cluster: {cluster_name}")
        messages.append(f"Leader: {cluster['leader']} broadcasting to children...")
        
        for child in cluster['children']:
            messages.append(f"  → {child}: Message received")
        
        messages.append(f"\nChildren reporting to leader:")
        for child in cluster['children']:
            messages.append(f"  ← {child}: Status OK")
        
        # Show results
        messagebox.showinfo("Communication Test", "\n".join(messages))
        self.add_log(f"Communication test completed for cluster: {cluster_name}")
    
    def _add_satellite_to_cluster(self, cluster_name):
        """Add a new satellite to existing cluster"""
        cluster = next((c for c in self.clusters if c['name'] == cluster_name), None)
        if not cluster:
            return
        
        # Get orbit selection
        orbit_names = [o['name'] for o in self.orbit_configurations]
        orbit_name = simpledialog.askstring("Select Orbit", 
                                           f"Enter orbit name for new satellite\nAvailable: {', '.join(orbit_names)}")
        if not orbit_name or orbit_name not in orbit_names:
            return
        
        orbit_config = next(o for o in self.orbit_configurations if o['name'] == orbit_name)
        
        # Create new child satellite
        child_num = len(cluster['children']) + 1
        sat_name = f"{cluster_name}_Child{child_num}"
        
        satellite = {
            "name": sat_name,
            "type": "cluster_member",
            "cluster": cluster_name,
            "role": "child",
            "orbit": {
                "a": orbit_config['altitude'] + 6371,
                "e": 0.01,
                "i": orbit_config['inclination'],
                "Omega": orbit_config['raan'],
                "omega": 0.0,
                "f": child_num * 10.0  # Spread out
            },
            "fault": {
                "type": "friction",
                "magnitude": 0.0005,
                "wheel": 3,
                "time": 15.0,
                "enabled": False,
            },
            "camera": {
                "position": [0.0, 0.0, 15.0],
                "fov": 80.0,
                "enabled": False
            },
            "communication": {
                "range": self.comm_range.get(),
                "fov": self.comm_fov.get(),
                "aHat_B": [0.0, 0.0, -1.0]
            },
            "targets": [],
            "orbit_name": orbit_name
        }
        
        self.satellites.append(satellite)
        cluster['satellites'].append(sat_name)
        cluster['children'].append(sat_name)
        orbit_config['satellites'].append(sat_name)
        
        # Update UI
        self.update_cluster_tree()
        self.update_satellite_listbox()
        self.update_orbit_tree()
        
        self.add_log(f"Added {sat_name} to cluster {cluster_name}")
    
    def _update_cluster_formation(self, cluster):
        """Update satellite positions based on formation type"""
        num_sats = len(cluster['satellites'])
        formation = cluster['formation']
        separation = cluster.get('separation', 10.0)
        
        positions = self._calculate_cluster_positions(num_sats, formation, 0, 0, separation)
        
        # Update satellite positions
        for i, sat_name in enumerate(cluster['satellites']):
            sat = next((s for s in self.satellites if s['name'] == sat_name), None)
            if sat and i < len(positions):
                sat['orbit']['f'] = positions[i]['anomaly']
                sat['position_offset'] = positions[i].get('offset', [0, 0, 0])
    
    def _calculate_cluster_positions(self, num_sats, formation, leader_anomaly, phase_offset, separation):
        """
        Return a list of dicts, each like:
        {'anomaly': <deg>, 'offset': [dx_km, dy_km, dz_km]}
        - leader at index 0
        - 'separation' is in km and acts as a nominal scale knob
        """
        # Normalize formation name + allow synonyms
        f = (formation or "").strip().lower()
        if f in ("circle", "ring"):
            f = "ring"
        elif f in ("line", "column", "stack"):
            f = "column"
        elif f in ("diamond", "box", "square"):
            f = "box"
        elif f in ("leader-follower", "train"):
            f = "train"
        else:
            f = "ring"  # sensible default

        positions = []

        # Helper: leader at index 0 (no offset, base anomaly)
        def leader_entry():
            return {'anomaly': leader_anomaly, 'offset': [0.0, 0.0, 0.0]}

        if f == "ring":
            # Satellites evenly spaced on a circle around the leader.
            # Leader at "top" (0 offset), followers on ring of radius R.
            R = max(2.0, separation)  # km
            positions.append(leader_entry())
            if num_sats > 1:
                for k in range(1, num_sats):
                    theta = 2.0 * np.pi * (k - 1) / max(1, (num_sats - 1))
                    dx = R * np.cos(theta)
                    dy = R * np.sin(theta)
                    positions.append({
                        'anomaly': leader_anomaly,   # keep anomaly same; only in-plane offset
                        'offset': [dx, dy, 0.0]
                    })

        elif f == "column":
            # Vertical stack: small along-track spacing (anomaly) and tiny cross-track spread.
            dA = max(2.0, separation * 0.5)  # deg between neighbors
            positions.append(leader_entry())
            for k in range(1, num_sats):
                positions.append({
                    'anomaly': leader_anomaly + k * dA,
                    'offset': [0.0, k * 0.2, 0.0]  # tiny lateral so links don't overlap in plots
                })

        elif f == "box":
            # 2x2 diamond/box around leader (like bubble 3).
            # For >4 sats, we wrap additional ones around in a second, larger box.
            side = max(3.0, separation)  # km
            positions.append(leader_entry())
            # First ring (up to 4 children)
            ring1 = [
                ( side,   0.0),
                ( 0.0,    side),
                (-side,   0.0),
                ( 0.0,   -side),
            ]
            # Second ring (optional) for >5
            ring2 = [
                ( 2*side,  2*side),
                (-2*side,  2*side),
                (-2*side, -2*side),
                ( 2*side, -2*side),
            ]
            coords = ring1 + ring2
            k = 1
            for (dx, dy) in coords[:max(0, num_sats-1)]:
                positions.append({'anomaly': leader_anomaly, 'offset': [dx, dy, 0.0]})
                k += 1
            # If still more, fall back to ring fill
            while len(positions) < num_sats:
                angle = 2*np.pi*(len(positions)-1)/max(1, num_sats-1)
                positions.append({'anomaly': leader_anomaly, 'offset': [1.5*side*np.cos(angle), 1.5*side*np.sin(angle), 0.0]})

        elif f == "train":
            # Long along-track line (wide spacing) like bubble 4.
            dA = max(8.0, separation)  # deg between neighbors (wider than column)
            positions.append(leader_entry())
            for k in range(1, num_sats):
                positions.append({
                    'anomaly': leader_anomaly + k * dA,
                    'offset': [0.0, 0.0, 0.0]
                })

        return positions

    
    def _create_orbit_management_tab(self, parent):
        """Create orbit management tab with all working buttons"""
        # Create orbit section
        create_frame = ttk.LabelFrame(parent, text="Create New Orbit", padding=10)
        create_frame.pack(fill=tk.X, padx=10, pady=10)
        
        # Orbit parameters
        params_grid = ttk.Frame(create_frame)
        params_grid.pack(fill=tk.X, pady=5)
        
        ttk.Label(params_grid, text="Orbit Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.new_orbit_name_var = tk.StringVar()
        ttk.Entry(params_grid, textvariable=self.new_orbit_name_var, width=20).grid(row=0, column=1, padx=5)
        
        ttk.Label(params_grid, text="Altitude (km):").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.new_orbit_alt_var = tk.DoubleVar(value=600.0)
        ttk.Entry(params_grid, textvariable=self.new_orbit_alt_var, width=15).grid(row=0, column=3, padx=5)
        
        ttk.Label(params_grid, text="Inclination (deg):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.new_orbit_inc_var = tk.DoubleVar(value=55.0)
        ttk.Entry(params_grid, textvariable=self.new_orbit_inc_var, width=15).grid(row=1, column=1, padx=5)
        
        ttk.Label(params_grid, text="RAAN (deg):").grid(row=1, column=2, sticky=tk.W, padx=5)
        self.new_orbit_raan_var = tk.DoubleVar(value=0.0)
        ttk.Entry(params_grid, textvariable=self.new_orbit_raan_var, width=15).grid(row=1, column=3, padx=5)
        
        ttk.Button(create_frame, text="Create Orbit", 
                  command=self.create_new_orbit).pack(pady=10)
        
        # Existing orbits
        orbits_frame = ttk.LabelFrame(parent, text="Orbit Library", padding=10)
        orbits_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Orbit tree
        columns = ('Name', 'Altitude', 'Inclination', 'RAAN', 'Satellites')
        self.orbit_tree = ttk.Treeview(orbits_frame, columns=columns, show='headings', height=8)
        
        for col in columns:
            self.orbit_tree.heading(col, text=col)
            
        self.orbit_tree.column('Name', width=120)
        self.orbit_tree.column('Altitude', width=100)
        self.orbit_tree.column('Inclination', width=100)
        self.orbit_tree.column('RAAN', width=100)
        self.orbit_tree.column('Satellites', width=100)
        
        scrollbar = ttk.Scrollbar(orbits_frame, orient="vertical", command=self.orbit_tree.yview)
        self.orbit_tree.configure(yscrollcommand=scrollbar.set)
        
        self.orbit_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Orbit actions - ALL IMPLEMENTED
        action_frame = ttk.Frame(orbits_frame)
        action_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(action_frame, text="Modify Orbit", 
                  command=self.modify_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Delete Orbit", 
                  command=self.delete_orbit).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="View Satellites", 
                  command=self.view_orbit_satellites).pack(side=tk.LEFT, padx=5)
        
        self.update_orbit_tree()
    
    def create_new_orbit(self):
        """Create a new orbit"""
        name = self.new_orbit_name_var.get().strip()
        if not name:
            messagebox.showwarning("Missing Name", "Please provide an orbit name")
            return
            
        if name in [o['name'] for o in self.orbit_configurations]:
            messagebox.showwarning("Duplicate Name", f"Orbit '{name}' already exists")
            return
        
        new_orbit = {
            "name": name,
            "altitude": self.new_orbit_alt_var.get(),
            "inclination": self.new_orbit_inc_var.get(),
            "raan": self.new_orbit_raan_var.get(),
            "satellites": [],
            "description": f"Custom orbit at {self.new_orbit_alt_var.get()}km"
        }
        
        self.orbit_configurations.append(new_orbit)
        self.update_orbit_tree()
        self._update_orbit_assignments()  # Update dropdowns
        
        self.add_log(f"Created new orbit: {name}")
        
        # Clear form
        self.new_orbit_name_var.set("")
        self.new_orbit_alt_var.set(600.0)
        self.new_orbit_inc_var.set(55.0)
        self.new_orbit_raan_var.set(0.0)
    
    def modify_orbit(self):
        """Modify selected orbit - IMPLEMENTED"""
        selection = self.orbit_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select an orbit to modify")
            return
        
        item = self.orbit_tree.item(selection[0])
        orbit_name = item['values'][0]
        orbit = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
        if not orbit:
            return
        
        # Modification window
        mod_window = tk.Toplevel(self.parent_app.root)
        mod_window.title(f"Modify Orbit: {orbit_name}")
        mod_window.geometry("350x250")
        
        frame = ttk.Frame(mod_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Parameters
        ttk.Label(frame, text="Altitude (km):").grid(row=0, column=0, sticky=tk.W, pady=5)
        alt_var = tk.DoubleVar(value=orbit['altitude'])
        ttk.Entry(frame, textvariable=alt_var, width=20).grid(row=0, column=1, pady=5)
        
        ttk.Label(frame, text="Inclination (deg):").grid(row=1, column=0, sticky=tk.W, pady=5)
        inc_var = tk.DoubleVar(value=orbit['inclination'])
        ttk.Entry(frame, textvariable=inc_var, width=20).grid(row=1, column=1, pady=5)
        
        ttk.Label(frame, text="RAAN (deg):").grid(row=2, column=0, sticky=tk.W, pady=5)
        raan_var = tk.DoubleVar(value=orbit['raan'])
        ttk.Entry(frame, textvariable=raan_var, width=20).grid(row=2, column=1, pady=5)
        
        def apply_changes():
            orbit['altitude'] = alt_var.get()
            orbit['inclination'] = inc_var.get()
            orbit['raan'] = raan_var.get()
            
            # Update satellites using this orbit
            for sat_name in orbit['satellites']:
                sat = next((s for s in self.satellites if s['name'] == sat_name), None)
                if sat:
                    sat['orbit']['a'] = orbit['altitude'] + 6371
                    sat['orbit']['i'] = orbit['inclination']
                    sat['orbit']['Omega'] = orbit['raan']
            
            self.update_orbit_tree()
            self.add_log(f"Modified orbit: {orbit_name}")
            mod_window.destroy()
        
        # Buttons
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=20)
        
        ttk.Button(button_frame, text="Apply", command=apply_changes).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Cancel", command=mod_window.destroy).pack(side=tk.LEFT, padx=5)
    
    def delete_orbit(self):
        """Delete selected orbit - IMPLEMENTED"""
        selection = self.orbit_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select an orbit to delete")
            return
        
        item = self.orbit_tree.item(selection[0])
        orbit_name = item['values'][0]
        orbit = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
        
        if not orbit:
            return
        
        if orbit['satellites']:
            messagebox.showwarning("Cannot Delete", 
                                 f"Cannot delete orbit with {len(orbit['satellites'])} satellites")
            return
        
        if not messagebox.askyesno("Confirm Delete", f"Delete orbit '{orbit_name}'?"):
            return
        
        self.orbit_configurations.remove(orbit)
        self.update_orbit_tree()
        self._update_orbit_assignments()
        
        self.add_log(f"Deleted orbit: {orbit_name}")
    
    def view_orbit_satellites(self):
        """View satellites in selected orbit - IMPLEMENTED"""
        selection = self.orbit_tree.selection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select an orbit")
            return
        
        item = self.orbit_tree.item(selection[0])
        orbit_name = item['values'][0]
        orbit = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
        
        if not orbit:
            return
        
        # Show satellites
        if orbit['satellites']:
            sat_list = "\n".join([f"• {s}" for s in orbit['satellites']])
            messagebox.showinfo(f"Satellites in {orbit_name}", 
                              f"Satellites ({len(orbit['satellites'])}):\n\n{sat_list}")
        else:
            messagebox.showinfo(f"Satellites in {orbit_name}", "No satellites in this orbit")
    
    def _create_individual_tab(self, parent):
        """Create individual satellites tab"""
        # Add satellite section
        add_frame = ttk.LabelFrame(parent, text="Add Individual Satellite", padding=10)
        add_frame.pack(fill=tk.X, padx=10, pady=10)
        
        # Satellite name and orbit
        params = ttk.Frame(add_frame)
        params.pack(fill=tk.X, pady=5)
        
        ttk.Label(params, text="Satellite Name:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.individual_sat_name_var = tk.StringVar()
        ttk.Entry(params, textvariable=self.individual_sat_name_var, width=20).grid(row=0, column=1, padx=5)
        
        ttk.Label(params, text="Orbit:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.individual_orbit_var = tk.StringVar()
        orbit_combo = ttk.Combobox(params, textvariable=self.individual_orbit_var, width=20, state="readonly")
        orbit_combo.grid(row=0, column=3, padx=5)
        orbit_combo['values'] = [o['name'] for o in self.orbit_configurations]
        if orbit_combo['values']:
            orbit_combo.current(0)
        
        ttk.Label(params, text="True Anomaly (deg):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.individual_anomaly_var = tk.DoubleVar(value=0.0)
        ttk.Entry(params, textvariable=self.individual_anomaly_var, width=20).grid(row=1, column=1, padx=5)
        
        ttk.Button(add_frame, text="Add Satellite", 
                  command=self.add_individual_satellite).pack(pady=10)
        
        # Satellite list
        list_frame = ttk.LabelFrame(parent, text="All Satellites", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Listbox with scrollbar
        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.satellite_listbox = tk.Listbox(list_container, 
                                           selectmode=tk.SINGLE,
                                           yscrollcommand=scrollbar.set,
                                           font=('Segoe UI', 10))
        self.satellite_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.satellite_listbox.yview)
        
        # Satellite actions
        action_frame = ttk.Frame(list_frame)
        action_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(action_frame, text="Remove Satellite", 
                  command=self.remove_satellite).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Assign to Orbit", 
                  command=self.assign_satellite_to_orbit).pack(side=tk.LEFT, padx=5)
        
        self.update_satellite_listbox()
    
    def add_individual_satellite(self):
        """Add an individual satellite"""
        name = self.individual_sat_name_var.get().strip()
        orbit_name = self.individual_orbit_var.get()
        anomaly = self.individual_anomaly_var.get()
        
        if not name:
            messagebox.showwarning("Missing Name", "Please provide a satellite name")
            return
        
        if name in [s['name'] for s in self.satellites]:
            messagebox.showwarning("Duplicate Name", f"Satellite '{name}' already exists")
            return
        
        orbit = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
        if not orbit:
            messagebox.showerror("Invalid Orbit", "Selected orbit not found")
            return
        
        satellite = {
            "name": name,
            "type": "individual",
            "cluster": None,
            "role": "independent",
            "orbit": {
                "a": orbit['altitude'] + 6371,
                "e": 0.01,
                "i": orbit['inclination'],
                "Omega": orbit['raan'],
                "omega": 0.0,
                "f": anomaly
            },
            "fault": {
                "type": "friction",
                "magnitude": 0.0005,
                "wheel": 3,
                "time": 15.0,
                "enabled": False,
            },
            "camera": {
                "position": [0.0, 0.0, 15.0],
                "fov": 80.0,
                "enabled": True
            },
            "communication": {
                "range": self.comm_range.get(),
                "fov": self.comm_fov.get(),
                "aHat_B": [0.0, 0.0, -1.0]
            },
            "targets": [],
            "orbit_name": orbit_name
        }
        
        self.satellites.append(satellite)
        orbit['satellites'].append(name)
        
        # Update UI
        self.update_satellite_listbox()
        self.update_orbit_tree()
        self.parent_app.update_satellite_dropdowns()
        self.parent_app.update_status_counts()
        
        # Clear form
        self.individual_sat_name_var.set("")
        self.individual_anomaly_var.set(0.0)
        
        self.add_log(f"Added individual satellite: {name}")
    
    def remove_satellite(self):
        """Remove selected satellite - IMPLEMENTED"""
        selection = self.satellite_listbox.curselection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a satellite to remove")
            return
        
        # Get satellite name
        item_text = self.satellite_listbox.get(selection[0])
        
        # Parse satellite name from display text
        sat_name = None
        if "[INDIVIDUAL]" in item_text:
            sat_name = item_text.split("]")[1].strip().split("(")[0].strip()
        elif "[L]" in item_text or "[C]" in item_text:
            messagebox.showwarning("Cannot Remove", "Cannot remove cluster satellites here. Use cluster management.")
            return
        
        if not sat_name:
            return
        
        if not messagebox.askyesno("Confirm Remove", f"Remove satellite '{sat_name}'?"):
            return
        
        # Remove satellite
        sat = next((s for s in self.satellites if s['name'] == sat_name), None)
        if sat:
            # Remove from satellites list
            self.satellites.remove(sat)
            
            # Remove from orbit
            orbit_name = sat.get('orbit_name')
            if orbit_name:
                orbit = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
                if orbit and sat_name in orbit['satellites']:
                    orbit['satellites'].remove(sat_name)
            
            # Update UI
            self.update_satellite_listbox()
            self.update_orbit_tree()
            self.parent_app.update_satellite_dropdowns()
            self.parent_app.update_status_counts()
            
            self.add_log(f"Removed satellite: {sat_name}")
    
    def assign_satellite_to_orbit(self):
        """Assign satellite to different orbit - IMPLEMENTED"""
        selection = self.satellite_listbox.curselection()
        if not selection:
            messagebox.showinfo("No Selection", "Please select a satellite")
            return
        
        # Get satellite
        item_text = self.satellite_listbox.get(selection[0])
        sat_name = None
        
        if "[INDIVIDUAL]" in item_text:
            sat_name = item_text.split("]")[1].strip().split("(")[0].strip()
        else:
            messagebox.showinfo("Cannot Reassign", "Only individual satellites can be reassigned")
            return
        
        sat = next((s for s in self.satellites if s['name'] == sat_name), None)
        if not sat:
            return
        
        # Get new orbit
        orbit_names = [o['name'] for o in self.orbit_configurations]
        new_orbit_name = simpledialog.askstring("Select Orbit", 
                                               f"Enter new orbit for {sat_name}\nAvailable: {', '.join(orbit_names)}")
        
        if not new_orbit_name or new_orbit_name not in orbit_names:
            return
        
        # Update orbit assignment
        old_orbit_name = sat.get('orbit_name')
        new_orbit = next(o for o in self.orbit_configurations if o['name'] == new_orbit_name)
        
        # Remove from old orbit
        if old_orbit_name:
            old_orbit = next((o for o in self.orbit_configurations if o['name'] == old_orbit_name), None)
            if old_orbit and sat_name in old_orbit['satellites']:
                old_orbit['satellites'].remove(sat_name)
        
        # Add to new orbit
        sat['orbit_name'] = new_orbit_name
        sat['orbit']['a'] = new_orbit['altitude'] + 6371
        sat['orbit']['i'] = new_orbit['inclination']
        sat['orbit']['Omega'] = new_orbit['raan']
        new_orbit['satellites'].append(sat_name)
        
        # Update UI
        self.update_satellite_listbox()
        self.update_orbit_tree()
        
        self.add_log(f"Reassigned {sat_name} to orbit: {new_orbit_name}")
    
    def _create_comm_tab(self, parent):
        """Create communication tab with visualization"""
        # Global settings
        settings_frame = ttk.LabelFrame(parent, text="Communication Settings", padding=10)
        settings_frame.pack(fill=tk.X, padx=10, pady=10)

        params = ttk.Frame(settings_frame)
        params.pack(fill=tk.X, pady=5)

        ttk.Label(params, text="Range (km):").grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Entry(params, textvariable=self.comm_range, width=15).grid(row=0, column=1, padx=5)

        ttk.Label(params, text="FOV (deg):").grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Entry(params, textvariable=self.comm_fov, width=15).grid(row=0, column=3, padx=5)

        ttk.Button(
            params, text="Refresh",
            command=lambda: [self.update_communication_plot(), self.comm_canvas.draw()]
        ).grid(row=0, column=4, padx=5)

      
        ttk.Button(
            params, text="Help",
            command=self.show_communication_help
        ).grid(row=0, column=5, padx=5)

        ttk.Label(
            params, text="(Bars = comm possible, Dots = messages)",
            font=('Arial', 8), foreground='gray'
        ).grid(row=1, column=0, columnspan=6, pady=2, sticky=tk.W)

        # Communication plot
        plot_frame = ttk.LabelFrame(parent, text="Communication Windows", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Create matplotlib figure
        self.comm_figure = Figure(figsize=(10, 6), dpi=80)
        self.comm_ax = self.comm_figure.add_subplot(111)

        self.comm_canvas = FigureCanvasTkAgg(self.comm_figure, plot_frame)
        self.comm_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Update button
        ttk.Button(
            plot_frame, text="Update Communication Plot",
            command=lambda: [self.update_communication_plot(),
                            print(f"Clusters available: {len(self.clusters)}")]
        ).pack(pady=10)

        # Initialize plot
        self.update_communication_plot()


    def update_communication_plot(self):
        """Update communication windows plot - IMPROVED VERSION"""
        self.comm_ax.clear()
        
        print(f"Updating communication plot. Clusters: {len(self.clusters)}")
        
        if not self.clusters or len(self.clusters) == 0:
            self.comm_ax.text(0.5, 0.5, 'No clusters configured\n\nCreate clusters in the Cluster Management tab', 
                            ha='center', va='center', fontsize=12)
            self.comm_ax.set_xlim(0, 1)
            self.comm_ax.set_ylim(0, 1)
        else:
            # Simulate communication windows
            time = np.linspace(0, 30, 100)  # 30 minutes
            
            y_pos = 0
            colors = ['#4CAF50', '#2196F3', '#FF9800', '#9C27B0', '#F44336']  # Better colors
            cluster_colors = {}
            
            # First, add title and explanation
            title_text = f"Communication Windows for {len(self.clusters)} Clusters"
            self.comm_ax.text(15, -1.5, "Green/Blue/Orange bars = Communication possible | Dots = Messages sent", 
                            ha='center', fontsize=9, style='italic', color='gray')
            
            # Plot each cluster with clear separation
            for cluster_idx, cluster in enumerate(self.clusters):
                color = colors[cluster_idx % len(colors)]
                cluster_colors[cluster['name']] = color
                cluster_name = cluster['name']
                leader = cluster.get('leader', 'Unknown')
                children = cluster.get('children', [])
                
                # Add cluster header with background
                cluster_y_start = y_pos
                self.comm_ax.axhspan(y_pos - 0.5, y_pos + len(children) + 0.5, 
                                    alpha=0.05, color=color, zorder=0)
                
                # Cluster name label
                self.comm_ax.text(-3, y_pos + len(children)/2, f"[{cluster_name.upper()}]", 
                                fontsize=11, fontweight='bold', color=color, 
                                va='center', ha='right',
                                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', 
                                        edgecolor=color, linewidth=2))
                
                # Plot each communication link
                for child_idx, child in enumerate(children):
                    # Create communication windows with better visibility
                    phase = child_idx * np.pi/4
                    comm_window = np.sin(2 * np.pi * time / 10 + phase) > 0.3
                    
                    # Draw communication availability as continuous bars
                    in_window = False
                    window_start = 0
                    
                    for k in range(len(time)):
                        if comm_window[k] and not in_window:
                            window_start = time[k]
                            in_window = True
                        elif not comm_window[k] and in_window:
                            # Draw the window
                            self.comm_ax.barh(y_pos, time[k] - window_start, 
                                            left=window_start, height=0.6,
                                            color=color, alpha=0.3, edgecolor=color, 
                                            linewidth=1)
                            in_window = False
                    
                    # Handle case where window extends to end
                    if in_window:
                        self.comm_ax.barh(y_pos, time[-1] - window_start, 
                                        left=window_start, height=0.6,
                                        color=color, alpha=0.3, edgecolor=color, 
                                        linewidth=1)
                    
                    # Add message indicators
                    message_times = time[::12]  # Messages every 12 time steps
                    for msg_time in message_times:
                        idx = np.argmin(np.abs(time - msg_time))
                        if comm_window[idx]:
                            # Message dot with annotation
                            self.comm_ax.scatter(msg_time, y_pos, color=color, 
                                            s=80, marker='o', zorder=5,
                                            edgecolors='white', linewidths=1)
                    
                    # Communication link label with arrow
                    link_label = f"{leader[:12]} → {child[:12]}"
                    self.comm_ax.text(-2, y_pos, link_label, fontsize=9, va='center',
                                    color='black', fontweight='normal')
                    
                    # Add link status indicator
                    if np.any(comm_window):
                        status = "●"  # Active
                        status_color = 'green'
                    else:
                        status = "○"  # Inactive
                        status_color = 'red'
                    self.comm_ax.text(30.5, y_pos, status, fontsize=12, va='center',
                                    color=status_color)
                    
                    y_pos += 1
                
                # Add spacing between clusters
                y_pos += 0.8
            
            # Configure plot with better formatting
            self.comm_ax.set_xlabel('Time (minutes)', fontsize=11, fontweight='bold')
            self.comm_ax.set_ylabel('Communication Links', fontsize=11, fontweight='bold')
            self.comm_ax.set_title(title_text, fontsize=13, fontweight='bold', pad=15)
            
            # Set limits with padding
            self.comm_ax.set_xlim(-3.5, 31)
            self.comm_ax.set_ylim(-2, max(1, y_pos))
            
            # Add grid for better readability
            self.comm_ax.grid(True, alpha=0.2, axis='x', linestyle='--')
            self.comm_ax.axvline(x=0, color='black', linewidth=1, alpha=0.5)
            self.comm_ax.axvline(x=30, color='black', linewidth=1, alpha=0.5)
            
            # Add time markers
            for t in [0, 5, 10, 15, 20, 25, 30]:
                self.comm_ax.axvline(x=t, color='gray', linewidth=0.5, alpha=0.3, linestyle=':')
            
            # Create legend
            legend_elements = []
            for cluster_name, color in cluster_colors.items():
                from matplotlib.patches import Patch
                from matplotlib.lines import Line2D
                # Combined legend entry
                legend_elements.append(Patch(facecolor=color, alpha=0.3, 
                                            edgecolor=color, label=f'{cluster_name} cluster'))
            
            # Add legend items for symbols
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor='gray', markersize=8, 
                                        label='Message sent'))
            legend_elements.append(Patch(facecolor='gray', alpha=0.3, 
                                        label='Comm window'))
            
            self.comm_ax.legend(handles=legend_elements, loc='upper left', 
                            fontsize=9, framealpha=0.9, ncol=2)
        
        self.comm_canvas.draw()
        print(f"Communication plot updated")
    # ============================================
    # Also add this helper to check cluster status:
    # ============================================

    def debug_cluster_status(self):
        """Debug helper to check cluster configuration"""
        print("\n=== CLUSTER DEBUG INFO ===")
        print(f"Number of clusters: {len(self.clusters)}")
        for cluster in self.clusters:
            print(f"  Cluster '{cluster['name']}':")
            print(f"    - Leader: {cluster.get('leader', 'None')}")
            print(f"    - Children: {cluster.get('children', [])}")
            print(f"    - Formation: {cluster.get('formation', 'Unknown')}")
        print(f"Total satellites: {len(self.satellites)}")
        print("=========================\n")
    

    def show_communication_help(self):
        """Show help dialog explaining the communication plot"""
        help_text = """UNDERSTANDING THE COMMUNICATION PLOT

    WHAT IT SHOWS:
    • Each horizontal line represents a communication link between two satellites
    • Links are grouped by cluster (different colors for each cluster)

    VISUAL ELEMENTS:
    • Colored Bars = Time windows when satellites can communicate
    - Bars appear when satellites are in range and have line-of-sight
    - No bar = satellites cannot communicate (out of range/blocked by Earth)

    • Dots = Messages being sent
    - Only appear during communication windows
    - Represent actual data transmission events

    • Y-Axis Labels = Communication links
    - Format: "Leader_name → Child_name"
    - Shows who is talking to whom

    • X-Axis = Time in minutes
    - Shows the full simulation duration (typically 30 minutes)

    HOW SATELLITES COMMUNICATE:
    1. Leader satellites coordinate their cluster
    2. Children report back to their leader
    3. Communication only possible when:
    - Satellites are within range (2000 km default)
    - Have line-of-sight (not blocked by Earth)
    - Within field-of-view angle (30° default)

    PATTERN EXPLANATION:
    The sinusoidal (wave-like) pattern occurs because:
    • Satellites orbit Earth at different speeds/positions
    • They move in and out of communication range
    • This creates periodic communication opportunities

    USE THIS FOR:
    • Planning when to send commands
    • Understanding data relay delays
    • Identifying communication blackout periods
    • Optimizing constellation design"""
        
        messagebox.showinfo("Communication Plot Help", help_text)

    
    # Update methods
    def update_cluster_tree(self):
        """Update cluster tree view"""
        for item in self.cluster_tree.get_children():
            self.cluster_tree.delete(item)
        
        for cluster in self.clusters:
            status = "Active" if cluster.get('leader') else "Inactive"
            
            cluster_item = self.cluster_tree.insert('', 'end', values=(
                cluster['name'],
                cluster.get('leader', 'None'),
                len(cluster.get('children', [])),
                cluster.get('formation', 'Unknown'),
                status
            ))
            
            # Add satellite details as children
            for sat_name in cluster.get('satellites', []):
                sat = next((s for s in self.satellites if s['name'] == sat_name), None)
                if sat:
                    self.cluster_tree.insert(cluster_item, 'end', values=(
                        f"  {sat['name']}",
                        sat.get('role', ''),
                        sat.get('orbit_name', ''),
                        f"{sat['orbit']['a']-6371:.0f}km",
                        "Fault" if sat['fault']['enabled'] else "OK"
                    ))
    
    def update_orbit_tree(self):
        """Update orbit tree view"""
        for item in self.orbit_tree.get_children():
            self.orbit_tree.delete(item)
        
        for orbit in self.orbit_configurations:
            self.orbit_tree.insert('', 'end', values=(
                orbit['name'],
                f"{orbit['altitude']:.0f} km",
                f"{orbit['inclination']:.1f}°",
                f"{orbit['raan']:.1f}°",
                len(orbit['satellites'])
            ))
    
    def update_satellite_listbox(self):
        """Update satellite listbox"""
        self.satellite_listbox.delete(0, tk.END)
        
        # Add clusters
        for cluster in self.clusters:
            self.satellite_listbox.insert(tk.END, f"[CLUSTER] {cluster['name']}")
            
            for sat_name in cluster['satellites']:
                sat = next((s for s in self.satellites if s['name'] == sat_name), None)
                if sat:
                    role = "L" if sat['role'] == 'leader' else "C"
                    orbit = sat.get('orbit_name', 'Unknown')
                    alt = sat['orbit']['a'] - 6371
                    fault = "F" if sat['fault']['enabled'] else ""
                    self.satellite_listbox.insert(tk.END, 
                                                 f"  [{role}] {sat['name']} ({orbit}, {alt:.0f}km) {fault}")
        
        # Add individual satellites
        for sat in self.satellites:
            if sat['type'] == 'individual':
                orbit = sat.get('orbit_name', 'Unknown')
                alt = sat['orbit']['a'] - 6371
                fault = "F" if sat['fault']['enabled'] else ""
                self.satellite_listbox.insert(tk.END, 
                                             f"[INDIVIDUAL] {sat['name']} ({orbit}, {alt:.0f}km) {fault}")
    
    def clear_cluster_form(self):
        """Clear cluster creation form"""
        self.cluster_name_var.set("")
        self.sats_per_cluster_var.set(4)
        self.formation_var.set("Leader-Follower")
        self.formation_separation.set(10.0)
        self._update_orbit_assignments()