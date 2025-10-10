#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
target_tab.py

Implements the Target Management tab with improved target-satellite connections.

Notes:
- This file now also exposes helpers to export targets into the exact
  structures expected by `spacecraft_simulation.py` without requiring a separate
  targets.py module. No existing code here has been removed; only safe additions.
"""
import tkinter as tk
from tkinter import ttk, messagebox, colorchooser
from .base_tab import BaseTab
import numpy as np

# ---- Optional Basilisk bits (only used for math constants if present) ----
try:
    from Basilisk.utilities import macros
    _D2R = macros.D2R
except Exception:
    _D2R = np.pi / 180.0


class TargetTab(BaseTab):
    """Target Management tab with improved target visibility"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the target tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        self.targets = parent_app.targets
        
        # Create the tab UI
        self.create_tab_ui()
        
        # Auto-assign targets to satellites for better Vizard visibility
        self.auto_assign_targets_on_startup()

    # ---------------------------------------------------------------------
    # NEW: helpers for other modules (no removals of your code)
    # ---------------------------------------------------------------------
    def build_sim_targets(self):
        """
        Convert the UI target dicts into TargetDefinition objects that
        `spacecraft_simulation.py` expects (no circular import at import-time).

        Returns
        -------
        list[TargetDefinition]
        """
        try:
            # Lazy import to avoid circulars during app startup
            from spacecraft_simulation import TargetDefinition
        except Exception:
            # If import fails, provide a clear error early
            raise ImportError(
                "Could not import TargetDefinition from spacecraft_simulation. "
                "Ensure spacecraft_simulation.py is on PYTHONPATH and importable."
            )
        out = []
        for t in self.targets:
            td = TargetDefinition(
                name=t.get("name", "Target"),
                latitude=float(t.get("lat", 0.0)),
                longitude=float(t.get("lon", 0.0)),
                color=t.get("color", "#FF0000"),
                priority=int(t.get("priority", 1))
            )
            # Preserve existing assignments
            td.assigned_to = list(t.get("assigned_to", []))
            out.append(td)
        return out

    def export_targets_for_config(self):
        """
        Return targets as a list of plain dicts (useful if another part of your
        app wants to write config files or skip the class objects).

        Keys: name, lat, lon, color, priority, assigned_to
        """
        out = []
        for t in self.targets:
            out.append({
                "name": t.get("name", "Target"),
                "lat": float(t.get("lat", 0.0)),
                "lon": float(t.get("lon", 0.0)),
                "color": t.get("color", "#FF0000"),
                "priority": int(t.get("priority", 1)),
                "assigned_to": list(t.get("assigned_to", []))
            })
        return out

    def coverage_summary(self):
        """
        Summarize coverage status for quick logging/GUI.

        Returns
        -------
        dict with total, assigned, visible_estimate, unassigned
        """
        total = len(self.targets)
        assigned = sum(1 for t in self.targets if t.get("assigned_to"))
        visible = 0
        for t in self.targets:
            st = self.calculate_target_coverage_status(t)
            if st.get("can_be_seen"):
                visible += 1
        return {
            "total_targets": total,
            "assigned_targets": assigned,
            "visible_estimate": visible,
            "unassigned_targets": total - assigned
        }

    @staticmethod
    def calculate_ground_target_position(lat_deg, lon_deg, altitude_m=100_000.0):
        """
        (Utility) Convert lat/lon to simple ECEF-like position for any later use.
        Matches what our simulation helpers do (mean radius + altitude).
        """
        earth_radius = 6_371_000.0
        lat = float(lat_deg) * _D2R
        lon = float(lon_deg) * _D2R
        r = earth_radius + float(altitude_m)
        x = r * np.cos(lat) * np.cos(lon)
        y = r * np.cos(lat) * np.sin(lon)
        z = r * np.sin(lat)
        return [float(x), float(y), float(z)]
    # ---------------------------------------------------------------------

    def auto_assign_targets_on_startup(self):
        """Auto-assign targets to satellites on startup for immediate Vizard visibility"""
        try:
            if self.satellites and self.targets:
                # Assign each target to a satellite if not already assigned
                for i, target in enumerate(self.targets):
                    if not target.get("assigned_to", []):
                        # Assign to satellite by cycling through available satellites
                        sat_index = i % len(self.satellites)
                        satellite_name = self.satellites[sat_index]["name"]
                        target["assigned_to"] = [satellite_name]
                        self.parent_app.add_log(f"Auto-assigned {target['name']} to {satellite_name}")
                
                # Update the UI to reflect auto-assignments
                self.update_target_assignments()
                self.update_target_listbox()
                
                # Update constellation tab to show assignments
                try:
                    self.parent_app.constellation_tab.update_satellite_listbox()
                except:
                    pass
                    
        except Exception as e:
            print(f"Could not auto-assign targets: {e}")
        
    def create_tab_ui(self):
        """Create the Target Management tab UI with improved connections"""
        # Main container
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Add info frame at top for better user guidance
        info_frame = ttk.Frame(main_container)
        info_frame.pack(fill=tk.X, pady=(0, 10))
        
        info_text = "Targets are ground locations tracked by satellites. Only assigned targets will be visible in Vizard."
        ttk.Label(info_frame, text=info_text, style="Info.TLabel", wraplength=1000).pack(anchor=tk.W)
        
        # Auto-assign button for user convenience
        auto_frame = ttk.Frame(info_frame)
        auto_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(auto_frame, text="Auto-Assign All Targets", 
                  command=self.auto_assign_targets,
                  style="Run.TButton").pack(side=tk.LEFT, padx=5)
        
        ttk.Button(auto_frame, text="Clear All Assignments", 
                  command=self.clear_all_assignments).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(auto_frame, text="Check Coverage", 
                  command=self.check_target_coverage).pack(side=tk.LEFT, padx=5)
        
        # Top frame for target list and details
        top_frame = ttk.Frame(main_container)
        top_frame.pack(fill=tk.BOTH, expand=True)
        
        # Split top frame into left (list) and right (details + assignment)
        left_frame = ttk.Frame(top_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=(0, 5))
        
        right_frame = ttk.Frame(top_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))
        
        # Target list on the left (fixed width)
        list_frame = ttk.LabelFrame(left_frame, text="Target List", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True)
        list_frame.configure(width=350)  # Fixed width for consistent layout
        
        # Buttons for target management
        btn_frame = ttk.Frame(list_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        
        add_btn = ttk.Button(btn_frame, text="Add Target", 
                             command=self.add_new_target)
        add_btn.pack(side=tk.LEFT, padx=5)
        
        remove_btn = ttk.Button(btn_frame, text="Remove Target", 
                                command=self.remove_target)
        remove_btn.pack(side=tk.LEFT, padx=5)
        
        # Generate random targets button
        generate_btn = ttk.Button(btn_frame, text="Generate Random", 
                                 command=self.generate_random_targets)
        generate_btn.pack(side=tk.RIGHT, padx=5)
        
        # Target listbox with scrollbar
        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.target_listbox = tk.Listbox(list_container, 
                                        selectmode=tk.SINGLE,
                                        yscrollcommand=scrollbar.set,
                                        font=('Segoe UI', 10),
                                        exportselection=False)
        self.target_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.target_listbox.yview)
        
        # Populate the listbox
        self.update_target_listbox()
        
        # Bind selection event
        self.target_listbox.bind('<<ListboxSelect>>', self.on_target_selected)
        
        # Right side - split into details and assignment (top) and map (bottom)
        details_assignment_frame = ttk.Frame(right_frame)
        details_assignment_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Target details on the top-left of right side
        details_frame = ttk.LabelFrame(details_assignment_frame, text="Target Details", padding=10)
        details_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
        
        # Name frame
        name_frame = ttk.Frame(details_frame)
        name_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(name_frame, text="Name:").pack(side=tk.LEFT)
        self.target_name_var = tk.StringVar()
        name_entry = ttk.Entry(name_frame, textvariable=self.target_name_var)
        name_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        # Location parameters
        loc_frame = ttk.Frame(details_frame)
        loc_frame.pack(fill=tk.X, pady=5)
        
        # Latitude
        lat_frame = ttk.Frame(loc_frame)
        lat_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(lat_frame, text="Latitude:").pack(side=tk.LEFT)
        self.target_lat_var = tk.DoubleVar()
        lat_entry = ttk.Entry(lat_frame, textvariable=self.target_lat_var, width=10)
        lat_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(lat_frame, text="(-90 to 90)", style="Info.TLabel").pack(side=tk.LEFT)
        
        # Longitude
        lon_frame = ttk.Frame(loc_frame)
        lon_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(lon_frame, text="Longitude:").pack(side=tk.LEFT)
        self.target_lon_var = tk.DoubleVar()
        lon_entry = ttk.Entry(lon_frame, textvariable=self.target_lon_var, width=10)
        lon_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(lon_frame, text="(-180 to 180)", style="Info.TLabel").pack(side=tk.LEFT)
        
        # Priority
        priority_frame = ttk.Frame(loc_frame)
        priority_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(priority_frame, text="Priority:").pack(side=tk.LEFT)
        self.target_priority_var = tk.IntVar(value=1)
        priority_spinbox = ttk.Spinbox(priority_frame, from_=1, to=5, 
                                      textvariable=self.target_priority_var, width=5)
        priority_spinbox.pack(side=tk.LEFT, padx=5)
        ttk.Label(priority_frame, text="(1-5, higher = more important)", 
                style="Info.TLabel").pack(side=tk.LEFT)
        
        # Color selection
        color_frame = ttk.Frame(details_frame)
        color_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(color_frame, text="Color:").pack(side=tk.LEFT)
        
        self.target_color_var = tk.StringVar(value="#FF0000")
        self.color_preview = tk.Canvas(color_frame, width=30, height=20, bg=self.target_color_var.get())
        self.color_preview.pack(side=tk.LEFT, padx=5)
        
        color_btn = ttk.Button(color_frame, text="Select Color", 
                              command=self.choose_target_color)
        color_btn.pack(side=tk.LEFT, padx=5)
        
        # Update button
        update_btn = ttk.Button(details_frame, text="Update Target", 
                               command=self.update_current_target)
        update_btn.pack(pady=10)
        
        # Target assignment frame on the top-right of right side
        assign_frame = ttk.LabelFrame(details_assignment_frame, text="Target Assignment", padding=10)
        assign_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))
        
        # Assignment status
        status_frame = ttk.Frame(assign_frame)
        status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(status_frame, text="Assignment Status:", font=('Segoe UI', 10, 'bold')).pack(anchor=tk.W)
        self.assignment_status_label = ttk.Label(status_frame, text="No target selected", style="Info.TLabel")
        self.assignment_status_label.pack(anchor=tk.W, pady=2)
        
        # Coverage status
        self.coverage_status_label = ttk.Label(status_frame, text="", style="Info.TLabel")
        self.coverage_status_label.pack(anchor=tk.W, pady=2)
        
        # Satellite selection
        sat_frame = ttk.Frame(assign_frame)
        sat_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(sat_frame, text="Assign to Satellite:").pack(side=tk.LEFT)
        
        self.assign_satellite_var = tk.StringVar()
        self.assign_satellite_combo = ttk.Combobox(sat_frame, textvariable=self.assign_satellite_var, state="readonly")
        self.assign_satellite_combo.pack(side=tk.LEFT, padx=5)
        self.update_satellite_dropdown()
        
        # Assignment buttons
        assign_btn_frame = ttk.Frame(assign_frame)
        assign_btn_frame.pack(fill=tk.X, pady=5)
        
        assign_btn = ttk.Button(assign_btn_frame, text="Assign Target", 
                               command=self.assign_target)
        assign_btn.pack(side=tk.LEFT, padx=5)
        
        unassign_btn = ttk.Button(assign_btn_frame, text="Unassign Target", 
                                 command=self.unassign_target)
        unassign_btn.pack(side=tk.LEFT, padx=5)
        
        # Current assignments display
        assignments_frame = ttk.LabelFrame(assign_frame, text="Current Assignments", padding=5)
        assignments_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.assignments_text = tk.Text(assignments_frame, height=4, width=25, state="disabled")
        self.assignments_text.pack(fill=tk.BOTH, expand=True)
        
        # Target map at the bottom with proper size
        map_frame = ttk.LabelFrame(right_frame, text="Target Coverage Map", padding=10)
        map_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        
        # Simple world map using matplotlib with good size
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        
        # Good figure size for visibility
        self.map_figure = Figure(figsize=(10, 6), dpi=100)
        self.map_ax = self.map_figure.add_subplot(111)
        
        # Create canvas and toolbar for map interaction
        self.map_canvas = FigureCanvasTkAgg(self.map_figure, map_frame)
        self.map_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add navigation toolbar for zoom/pan functionality
        toolbar_frame = ttk.Frame(map_frame)
        toolbar_frame.pack(fill=tk.X, pady=(5, 0))
        
        self.map_toolbar = NavigationToolbar2Tk(self.map_canvas, toolbar_frame)
        self.map_toolbar.update()
        
        # Draw the initial map
        self.update_map()
        
        # Select the first target if available
        if self.targets:
            self.target_listbox.selection_set(0)
            self.on_target_selected(None)
            
    def update_satellite_dropdown(self):
        """Update the satellite dropdown with current satellites"""
        values = [sat['name'] for sat in self.satellites]
        self.assign_satellite_combo['values'] = values
        if values:
            # If current selection not valid, reset to first
            current = self.assign_satellite_var.get()
            if current not in values:
                self.assign_satellite_combo.current(0)
        else:
            self.assign_satellite_var.set("")
            
    def update_target_listbox(self):
        """Update the target listbox with current targets and coverage info"""
        self.target_listbox.delete(0, tk.END)
        for target in self.targets:
            # Add assignment count and coverage status to display
            assignment_count = len(target.get("assigned_to", []))
            coverage_status = self.calculate_target_coverage_status(target)
            
            display_name = f"{target['name']} "
            if assignment_count > 0:
                display_name += f"[{assignment_count} assigned] "
                if coverage_status["can_be_seen"]:
                    display_name += "VISIBLE"
                else:
                    display_name += "ASSIGNED"
            else:
                display_name += "[unassigned]"
            
            self.target_listbox.insert(tk.END, display_name)
            
    def calculate_target_coverage_status(self, target):
        """Calculate if target can be seen by assigned satellites"""
        if not target.get("assigned_to", []):
            return {"can_be_seen": False, "reason": "No satellites assigned"}
        
        # Simple coverage calculation based on satellite altitudes
        assigned_satellites = []
        for sat_name in target["assigned_to"]:
            for sat in self.satellites:
                if sat["name"] == sat_name:
                    assigned_satellites.append(sat)
                    break
        
        if not assigned_satellites:
            return {"can_be_seen": False, "reason": "Assigned satellites not found"}
        
        # Calculate if any satellite can see the target
        # For simplicity, assume satellites above 200km can see targets
        can_see = False
        best_altitude = 0
        
        for sat in assigned_satellites:
            altitude = sat["orbit"]["a"] - 6371  # Altitude in km
            if altitude > 200:  # Minimum altitude for good coverage
                can_see = True
                best_altitude = max(best_altitude, altitude)
        
        if can_see:
            return {
                "can_be_seen": True, 
                "reason": f"Visible from {best_altitude:.0f}km altitude"
            }
        else:
            return {
                "can_be_seen": False, 
                "reason": f"Satellites too low (need >200km altitude)"
            }
            
    def on_target_selected(self, event):
        """Handle target selection event"""
        selection = self.target_listbox.curselection()
        if selection:
            index = selection[0]
            self.load_target_details(index)
            self.update_target_assignments()
            
    def load_target_details(self, index):
        """Load target details into UI fields"""
        if 0 <= index < len(self.targets):
            target = self.targets[index]
            self.target_name_var.set(target["name"])
            self.target_lat_var.set(target["lat"])
            self.target_lon_var.set(target["lon"])
            self.target_color_var.set(target["color"])
            self.color_preview.config(bg=target["color"])
            
            # Set priority if it exists, otherwise default to 1
            self.target_priority_var.set(target.get("priority", 1))
            
            # Update assignment status with coverage info
            assigned_to = target.get("assigned_to", [])
            coverage_status = self.calculate_target_coverage_status(target)
            
            if assigned_to:
                status_text = f"Assigned to: {', '.join(assigned_to)}"
                if coverage_status["can_be_seen"]:
                    status_text += " - VISIBLE IN VIZARD"
                    self.assignment_status_label.config(text=status_text, foreground="green")
                else:
                    status_text += " - ASSIGNED BUT MAY NOT BE VISIBLE"
                    self.assignment_status_label.config(text=status_text, foreground="orange")
                
                # Show coverage details
                self.coverage_status_label.config(text=coverage_status["reason"], foreground="blue")
            else:
                self.assignment_status_label.config(text="Not assigned - WILL NOT BE VISIBLE IN VIZARD", foreground="red")
                self.coverage_status_label.config(text="", foreground="black")
            
    def update_target_assignments(self):
        """Update the target assignments display"""
        selection = self.target_listbox.curselection()
        if selection:
            index = selection[0]
            if index < len(self.targets):
                target = self.targets[index]
                
                self.assignments_text.config(state="normal")
                self.assignments_text.delete(1.0, tk.END)
                
                assigned_to = target.get("assigned_to", [])
                if assigned_to:
                    for sat_name in assigned_to:
                        # Find satellite to show altitude info
                        for sat in self.satellites:
                            if sat["name"] == sat_name:
                                altitude = sat["orbit"]["a"] - 6371
                                self.assignments_text.insert(tk.END, f"• {sat_name} ({altitude:.0f}km)\n")
                                break
                        else:
                            self.assignments_text.insert(tk.END, f"• {sat_name}\n")
                    
                    coverage_status = self.calculate_target_coverage_status(target)
                    if coverage_status["can_be_seen"]:
                        self.assignments_text.insert(tk.END, f"\nStatus: VISIBLE in Vizard\n")
                        self.assignments_text.insert(tk.END, f"Reason: {coverage_status['reason']}")
                    else:
                        self.assignments_text.insert(tk.END, f"\nStatus: MAY NOT BE VISIBLE\n")
                        self.assignments_text.insert(tk.END, f"Reason: {coverage_status['reason']}")
                else:
                    self.assignments_text.insert(tk.END, "No assignments\n\nStatus: NOT VISIBLE in Vizard")
                    
                self.assignments_text.config(state="disabled")
        
        # Update the map
        self.update_map()
        
    def check_target_coverage(self):
        """Check coverage for all targets and show results"""
        if not self.targets:
            messagebox.showinfo("No Targets", "No targets defined to check coverage.")
            return
            
        # Create coverage report
        coverage_window = tk.Toplevel(self.parent_app.root)
        coverage_window.title("Target Coverage Report")
        coverage_window.geometry("600x400")
        coverage_window.transient(self.parent_app.root)
        
        # Create scrolled text for report
        from tkinter import scrolledtext
        report_text = scrolledtext.ScrolledText(coverage_window, wrap=tk.WORD, height=20)
        report_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Generate report
        report = "TARGET COVERAGE ANALYSIS REPORT\n"
        report += "=" * 50 + "\n\n"
        
        visible_count = 0
        assigned_count = 0
        unassigned_count = 0
        
        for target in self.targets:
            report += f"Target: {target['name']}\n"
            report += f"Location: {target['lat']:.2f}°, {target['lon']:.2f}°\n"
            
            assigned_to = target.get("assigned_to", [])
            if assigned_to:
                assigned_count += 1
                report += f"Assigned to: {', '.join(assigned_to)}\n"
                
                # Get coverage status
                coverage_status = self.calculate_target_coverage_status(target)
                if coverage_status["can_be_seen"]:
                    visible_count += 1
                    report += f"Visibility: VISIBLE in Vizard\n"
                    report += f"Coverage: {coverage_status['reason']}\n"
                else:
                    report += f"Visibility: ASSIGNED but may not be visible\n"
                    report += f"Issue: {coverage_status['reason']}\n"
            else:
                unassigned_count += 1
                report += f"Assignment: UNASSIGNED\n"
                report += f"Visibility: NOT VISIBLE in Vizard\n"
            
            report += "\n" + "-" * 30 + "\n\n"
        
        # Summary
        report += "SUMMARY:\n"
        report += f"Total Targets: {len(self.targets)}\n"
        report += f"Assigned Targets: {assigned_count}\n"
        report += f"Visible Targets: {visible_count}\n"
        report += f"Unassigned Targets: {unassigned_count}\n\n"
        
        if visible_count < len(self.targets):
            report += "RECOMMENDATIONS:\n"
            report += "• Assign unassigned targets to satellites\n"
            report += "• Increase satellite altitudes above 200km for better coverage\n"
            report += "• Use 'Auto-Assign All Targets' button for quick assignment\n"
        
        report_text.insert(tk.END, report)
        report_text.config(state="disabled")
        
        # Close button
        ttk.Button(coverage_window, text="Close", 
                  command=coverage_window.destroy).pack(pady=10)
            
    def add_new_target(self):
        """Add a new target"""
        # Default values for new target
        new_target = {
            "name": f"Target{len(self.targets) + 1}",
            "lat": 0.0,
            "lon": 0.0,
            "color": "#FF0000",
            "priority": 1,
            "assigned_to": []
        }
        
        self.targets.append(new_target)
        
        # Update UI
        self.update_target_listbox()
        self.target_listbox.selection_clear(0, tk.END)
        self.target_listbox.selection_set(len(self.targets) - 1)
        self.on_target_selected(None)
        
        self.parent_app.add_log(f"Added new target: {new_target['name']}")
    
    def generate_random_targets(self):
        """Generate random target locations around the world"""
        # Ask how many targets to generate
        from tkinter import simpledialog
        num_targets = simpledialog.askinteger("Generate Targets", 
                                            "How many targets to generate?", 
                                            initialvalue=5, minvalue=1, maxvalue=20)
        if not num_targets:
            return
            
        # Clear existing targets if requested
        if self.targets:
            clear = messagebox.askyesno("Clear Existing Targets", 
                                      "Do you want to clear existing targets?")
            if clear:
                self.targets.clear()
        
        # Define some interesting locations
        predefined_cities = [
            {"name": "New York", "lat": 40.71, "lon": -74.00},
            {"name": "London", "lat": 51.51, "lon": -0.13},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77},
            {"name": "Sydney", "lat": -33.87, "lon": 151.21},
            {"name": "Rio de Janeiro", "lat": -22.91, "lon": -43.17},
            {"name": "Cape Town", "lat": -33.92, "lon": 18.42},
            {"name": "Moscow", "lat": 55.75, "lon": 37.62},
            {"name": "Beijing", "lat": 39.91, "lon": 116.40},
            {"name": "Mumbai", "lat": 19.07, "lon": 72.88},
            {"name": "Los Angeles", "lat": 34.05, "lon": -118.24},
            {"name": "Dubai", "lat": 25.20, "lon": 55.27},
            {"name": "Paris", "lat": 48.85, "lon": 2.35},
            {"name": "Cairo", "lat": 30.04, "lon": 31.23},
            {"name": "Singapore", "lat": 1.35, "lon": 103.82},
            {"name": "Melbourne", "lat": -37.81, "lon": 144.96},
            {"name": "Santiago", "lat": -33.46, "lon": -70.65},
            {"name": "Berlin", "lat": 52.52, "lon": 13.40},
            {"name": "Toronto", "lat": 43.65, "lon": -79.38},
            {"name": "Mexico City", "lat": 19.43, "lon": -99.13},
            {"name": "Istanbul", "lat": 41.01, "lon": 28.97}
        ]
        
        # Define colors to cycle through
        colors = ['#FF0000', '#0000FF', '#00FF00', '#FFFF00', '#FF00FF', 
                 '#00FFFF', '#FFA500', '#800080', '#008000', '#A52A2A']
        
        # Generate random targets
        import random
        for i in range(num_targets):
            # Use predefined cities if available
            if i < len(predefined_cities):
                city = predefined_cities[i]
                name = city["name"]
                lat = city["lat"]
                lon = city["lon"]
            else:
                # Generate random position
                lat = random.uniform(-80, 80)  # Avoid extreme polar regions
                lon = random.uniform(-180, 180)
                name = f"Target{len(self.targets) + 1}"
            
            # Use cycling colors
            color = colors[i % len(colors)]
            
            # Random priority
            priority = random.randint(1, 5)
            
            # Create new target
            new_target = {
                "name": name,
                "lat": lat,
                "lon": lon,
                "color": color,
                "priority": priority,
                "assigned_to": []
            }
            
            self.targets.append(new_target)
        
        # Auto-assign the new targets
        self.auto_assign_targets()
        
        # Update UI
        self.update_target_listbox()
        if self.targets:
            self.target_listbox.selection_clear(0, tk.END)
            self.target_listbox.selection_set(0)
            self.on_target_selected(None)
        
        self.parent_app.add_log(f"Generated {num_targets} random targets and auto-assigned them")
        
    def remove_target(self):
        """Remove the selected target"""
        selection = self.target_listbox.curselection()
        if selection:
            index = selection[0]
            target_name = self.targets[index]["name"]
            
            # Ask for confirmation
            confirm = messagebox.askyesno("Confirm Removal", 
                                         f"Are you sure you want to remove {target_name}?")
            if confirm:
                # Remove target
                self.targets.pop(index)
                
                # Update UI
                self.update_target_listbox()
                if self.targets:
                    self.target_listbox.selection_set(0)
                    self.on_target_selected(None)
                
                # Update constellation tab to reflect target changes
                try:
                    self.parent_app.constellation_tab.update_satellite_listbox()
                except:
                    pass
                
                self.parent_app.add_log(f"Removed target: {target_name}")
            
    def update_current_target(self):
        """Update the current target with UI values"""
        selection = self.target_listbox.curselection()
        if selection:
            index = selection[0]
            
            # Update with new values
            self.targets[index]["name"] = self.target_name_var.get()
            self.targets[index]["lat"] = self.target_lat_var.get()
            self.targets[index]["lon"] = self.target_lon_var.get()
            self.targets[index]["color"] = self.target_color_var.get()
            self.targets[index]["priority"] = self.target_priority_var.get()
            
            # Update UI
            self.update_target_listbox()
            self.target_listbox.selection_set(index)
            self.update_map()
            
            # Update constellation tab to reflect target changes
            try:
                self.parent_app.constellation_tab.update_satellite_listbox()
            except:
                pass
            
            self.parent_app.add_log(f"Updated target: {self.targets[index]['name']}")
            
    def choose_target_color(self):
        """Open color chooser dialog for target color"""
        color = colorchooser.askcolor(self.target_color_var.get())[1]
        if color:
            self.target_color_var.set(color)
            self.color_preview.config(bg=color)
            
    def assign_target(self):
        """Assign the selected target to a satellite"""
        target_selection = self.target_listbox.curselection()
        sat_name = self.assign_satellite_var.get()
        
        if target_selection and sat_name:
            target_index = target_selection[0]
            if target_index < len(self.targets):
                target = self.targets[target_index]
                
                # Initialize assigned_to if it doesn't exist
                if "assigned_to" not in target:
                    target["assigned_to"] = []
                
                # Check if already assigned
                if sat_name not in target["assigned_to"]:
                    target["assigned_to"].append(sat_name)
                    self.update_target_assignments()
                    self.update_target_listbox()
                    
                    # Keep selection
                    self.target_listbox.selection_set(target_index)
                    self.on_target_selected(None)
                    
                    # Update constellation tab to reflect target changes
                    try:
                        self.parent_app.constellation_tab.update_satellite_listbox()
                    except:
                        pass
                    
                    # Check coverage
                    coverage_status = self.calculate_target_coverage_status(target)
                    if coverage_status["can_be_seen"]:
                        self.parent_app.add_log(f"Assigned {target['name']} to {sat_name} - VISIBLE IN VIZARD")
                    else:
                        self.parent_app.add_log(f"Assigned {target['name']} to {sat_name} - {coverage_status['reason']}")
                else:
                    messagebox.showinfo("Already Assigned", 
                                      f"{target['name']} is already assigned to {sat_name}")
                
    def unassign_target(self):
        """Unassign the selected target from a satellite"""
        target_selection = self.target_listbox.curselection()
        sat_name = self.assign_satellite_var.get()
        
        if target_selection and sat_name:
            target_index = target_selection[0]
            if target_index < len(self.targets):
                target = self.targets[target_index]
                
                # Initialize assigned_to if it doesn't exist
                if "assigned_to" not in target:
                    target["assigned_to"] = []
                
                # Check if assigned
                if sat_name in target["assigned_to"]:
                    target["assigned_to"].remove(sat_name)
                    self.update_target_assignments()
                    self.update_target_listbox()
                    
                    # Keep selection
                    self.target_listbox.selection_set(target_index)
                    self.on_target_selected(None)
                    
                    # Update constellation tab to reflect target changes
                    try:
                        self.parent_app.constellation_tab.update_satellite_listbox()
                    except:
                        pass
                    
                    self.parent_app.add_log(f"Unassigned {target['name']} from {sat_name}")
                else:
                    messagebox.showinfo("Not Assigned", 
                                      f"{target['name']} is not assigned to {sat_name}")
    
    def auto_assign_targets(self):
        """Automatically assign targets to satellites for optimal coverage"""
        if not self.satellites or not self.targets:
            messagebox.showinfo("Cannot Assign", 
                              "Need at least one satellite and one target to auto-assign.")
            return
            
        # Get satellite names
        sat_names = [sat["name"] for sat in self.satellites]
        
        # Sort targets by priority (higher priority first)
        sorted_targets = sorted(self.targets, 
                               key=lambda t: t.get("priority", 1), 
                               reverse=True)
        
        # Simple assignment strategy: distribute targets evenly among satellites
        assignments_per_sat = {name: 0 for name in sat_names}
        
        for target in sorted_targets:
            # Find satellite with fewest assignments
            min_assignments = min(assignments_per_sat.values()) if assignments_per_sat else 0
            candidates = [name for name, count in assignments_per_sat.items() 
                        if count == min_assignments] if assignments_per_sat else []
            
            # Assign to first candidate
            if candidates:
                assigned_sat = candidates[0]
                target["assigned_to"] = [assigned_sat]  # Replace existing assignments
                assignments_per_sat[assigned_sat] += 1
            else:
                # Fallback (no satellites?): keep existing assignment or none
                target.setdefault("assigned_to", [])
        
        # Update UI
        self.update_target_assignments()
        self.update_target_listbox()
        
        # Update constellation tab to reflect target changes
        try:
            self.parent_app.constellation_tab.update_satellite_listbox()
        except:
            pass
        
        # Check coverage for assigned targets
        visible_count = 0
        for target in self.targets:
            coverage_status = self.calculate_target_coverage_status(target)
            if coverage_status["can_be_seen"]:
                visible_count += 1
        
        self.parent_app.add_log(f"Auto-assigned {len(self.targets)} targets to {len(self.satellites)} satellites")
        self.parent_app.add_log(f"{visible_count} targets will be VISIBLE in Vizard")
        
    def clear_all_assignments(self):
        """Clear all target assignments"""
        confirm = messagebox.askyesno("Confirm Clear", 
                                    "Clear all target assignments? Targets will not be visible in Vizard.")
        if confirm:
            for target in self.targets:
                target["assigned_to"] = []
            
            # Update UI
            self.update_target_assignments()
            self.update_target_listbox()
            
            # Update constellation tab
            try:
                self.parent_app.constellation_tab.update_satellite_listbox()
            except:
                pass
            
            self.parent_app.add_log("Cleared all target assignments - targets will not be visible in Vizard")
        
    def update_map(self):
        """Update the target map display with coverage information"""
        # Clear previous plot
        self.map_ax.clear()
        
        # Draw simple world map background
        self.map_ax.set_xlim(-180, 180)  
        self.map_ax.set_ylim(-90, 90)
        self.map_ax.grid(True, alpha=0.3)
        self.map_ax.set_title('Target Coverage Map (Assigned targets visible in Vizard)', fontsize=14, fontweight='bold')
        self.map_ax.set_xlabel('Longitude', fontsize=12)
        self.map_ax.set_ylabel('Latitude', fontsize=12)
        
        # Try to add simple continent outlines if available
        try:
            # Simple approximation of continent outlines
            continents = [
                # North America
                [(-168, 66), (-125, 72), (-91, 83), (-60, 74), (-54, 52), 
                 (-60, 30), (-85, 25), (-120, 32), (-130, 54), (-168, 66)],
                # South America
                [(-85, 15), (-35, 12), (-35, -35), (-75, -55), (-85, 15)],
                # Europe and Africa
                [(-10, 55), (40, 60), (50, 40), (40, 10), (50, -30), 
                 (18, -35), (-10, 5), (-10, 55)],
                # Asia and Australia
                [(50, 40), (100, 60), (140, 50), (140, 20), (100, 10), 
                 (140, -10), (155, -20), (130, -40), (110, -20), (75, -20), (50, 40)]
            ]
            
            # Draw each continent
            for continent in continents:
                xs, ys = zip(*continent)
                self.map_ax.plot(xs, ys, 'k-', linewidth=1.0, alpha=0.7)
                
        except Exception:
            # If we can't draw continents, just note it
            pass
        
        # Plot targets with coverage information
        for target in self.targets:
            # Extract lat/lon
            lat = target["lat"]
            lon = target["lon"]
            
            # Get color
            color = target["color"]
            
            # Determine marker size based on priority
            base_size = 100
            size = base_size + (target.get("priority", 1) * 30)
            
            # Determine marker style and color based on assignment and coverage
            assigned_to = target.get("assigned_to", [])
            coverage_status = self.calculate_target_coverage_status(target)
            
            if assigned_to:
                if coverage_status["can_be_seen"]:
                    marker = '*'  # Star for visible targets
                    alpha = 1.0
                    edge_color = 'green'
                    line_width = 3
                    status_text = "VISIBLE"
                else:
                    marker = 'o'  # Circle for assigned but potentially not visible
                    alpha = 0.8
                    edge_color = 'orange'
                    line_width = 2
                    status_text = "ASSIGNED"
            else:
                marker = 'X'  # X for unassigned targets
                alpha = 0.6
                edge_color = 'red'
                line_width = 2
                status_text = "NOT VISIBLE"
                
            # Draw marker
            self.map_ax.scatter(lon, lat, color=color, s=size, 
                               marker=marker, zorder=5, 
                               alpha=alpha,
                               edgecolor=edge_color, linewidth=line_width)
            
            # Add label with assignment and coverage status
            label_offset = 4  # Degrees offset for label
            if assigned_to:
                if coverage_status["can_be_seen"]:
                    label_text = f"{target['name']}\n{status_text}"
                    bbox_color = 'lightgreen'
                else:
                    label_text = f"{target['name']}\n{status_text}"
                    bbox_color = 'lightyellow'
            else:
                label_text = f"{target['name']}\n{status_text}"
                bbox_color = 'lightcoral'
                
            self.map_ax.text(lon, lat + label_offset, label_text, 
                            ha='center', va='bottom', fontsize=9, fontweight='bold',
                            bbox=dict(facecolor=bbox_color, alpha=0.8, boxstyle='round,pad=0.3',
                                     edgecolor=color, linewidth=1))
        
        # Add legend with coverage information
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='*', color='w', markerfacecolor='green', 
                   markersize=12, label='Visible Targets (VISIBLE in Vizard)', 
                   markeredgecolor='green', markeredgewidth=3),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', 
                   markersize=10, label='Assigned Targets (check altitude)', 
                   markeredgecolor='orange', markeredgewidth=2),
            Line2D([0], [0], marker='X', color='w', markerfacecolor='red', 
                   markersize=10, label='Unassigned Targets (NOT VISIBLE)', 
                   markeredgecolor='red', markeredgewidth=2)
        ]
        self.map_ax.legend(handles=legend_elements, loc='upper right')
        
        # Update canvas
        try:
            self.map_figure.tight_layout()
        except Exception:
            pass
        self.map_canvas.draw()
