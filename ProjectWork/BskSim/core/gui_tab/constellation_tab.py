#!/usr/bin/env python
"""
constellation_tab.py

Implements the Constellation Management tab for the Spacecraft Constellation Fault Simulator.
FIXED: Simplified satellite names, better orbital altitudes for visibility, improved display
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from .base_tab import BaseTab

class ConstellationTab(BaseTab):
    """Constellation Management tab implementation"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the constellation tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        self.current_satellite_index = parent_app.current_satellite_index
        
        # Create the tab UI
        self.create_tab_ui()
        
    def create_tab_ui(self):
        """Create the Constellation Management tab UI"""
        # Split into two frames
        left_frame = ttk.Frame(self.parent_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
        
        right_frame = ttk.Frame(self.parent_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))
        
        # Satellite list on the left
        list_frame = ttk.LabelFrame(left_frame, text="Satellite List", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True)
        
        # Buttons for satellite management
        btn_frame = ttk.Frame(list_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        
        add_btn = ttk.Button(btn_frame, text="Add Satellite", 
                             command=self.add_new_satellite)
        add_btn.pack(side=tk.LEFT, padx=5)
        
        # Add multi-satellite option
        multi_frame = ttk.Frame(btn_frame)
        multi_frame.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(multi_frame, text="Add Multiple:").pack(side=tk.LEFT)
        self.num_satellites_var = tk.IntVar(value=4)
        num_satellites_spin = ttk.Spinbox(multi_frame, from_=2, to=8, 
                                          textvariable=self.num_satellites_var, width=5)
        num_satellites_spin.pack(side=tk.LEFT, padx=5)
        
        add_multi_btn = ttk.Button(multi_frame, text="Add Constellation", 
                                  command=self.add_multiple_satellites)
        add_multi_btn.pack(side=tk.LEFT, padx=5)
        
        remove_btn = ttk.Button(btn_frame, text="Remove Satellite", 
                                command=self.remove_satellite)
        remove_btn.pack(side=tk.RIGHT, padx=5)
        
        # Display options
        display_frame = ttk.Frame(list_frame)
        display_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(display_frame, text="Display:").pack(side=tk.LEFT)
        self.show_details_var = tk.BooleanVar(value=False)
        details_check = ttk.Checkbutton(display_frame, text="Show Details", 
                                      variable=self.show_details_var,
                                      command=self.update_satellite_listbox)
        details_check.pack(side=tk.LEFT, padx=5)
        
        # Satellite listbox with scrollbar
        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.satellite_listbox = tk.Listbox(list_container, 
                                           selectmode=tk.SINGLE,
                                           yscrollcommand=scrollbar.set,
                                           font=('Segoe UI', 10),
                                           exportselection=False)
        self.satellite_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.satellite_listbox.yview)
        
        # Populate the listbox
        self.update_satellite_listbox()
        
        # Bind selection event
        self.satellite_listbox.bind('<<ListboxSelect>>', self.on_satellite_selected)
        
        # Satellite details on the right
        details_frame = ttk.LabelFrame(right_frame, text="Satellite Details", padding=10)
        details_frame.pack(fill=tk.BOTH, expand=True)
        
        # Name frame
        name_frame = ttk.Frame(details_frame)
        name_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(name_frame, text="Name:").pack(side=tk.LEFT)
        self.sat_name_var = tk.StringVar()
        name_entry = ttk.Entry(name_frame, textvariable=self.sat_name_var)
        name_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        # Orbit parameters
        orbit_frame = ttk.LabelFrame(details_frame, text="Orbit Parameters", padding=10)
        orbit_frame.pack(fill=tk.X, pady=5)
        
        # FIXED: Better orbital parameters for improved visibility
        orbit_params = [
            ("Semi-major axis (km):", "a", 8000.0),  # FIXED: Higher altitude for better visibility
            ("Eccentricity:", "e", 0.05),            # Small eccentricity for stable orbit
            ("Inclination (deg):", "i", 55.0),       # Good for ground coverage
            ("RAAN (deg):", "Omega", 45.0),          # Right ascension of ascending node
            ("Arg. of Periapsis (deg):", "omega", 30.0),  # Argument of periapsis
            ("True Anomaly (deg):", "f", 0.0)        # True anomaly - user specified
        ]
        
        self.orbit_vars = {}
        for i, (label_text, param_name, default_val) in enumerate(orbit_params):
            param_frame = ttk.Frame(orbit_frame)
            param_frame.pack(fill=tk.X, pady=2)
            
            ttk.Label(param_frame, text=label_text).pack(side=tk.LEFT)
            
            self.orbit_vars[param_name] = tk.DoubleVar(value=default_val)
            entry = ttk.Entry(param_frame, textvariable=self.orbit_vars[param_name], width=10)
            entry.pack(side=tk.RIGHT)
        
        # Orbital period calculator and presets
        period_frame = ttk.LabelFrame(orbit_frame, text="Orbital Period Control", padding=5)
        period_frame.pack(fill=tk.X, pady=5)
        
        # Current period display
        period_info_frame = ttk.Frame(period_frame)
        period_info_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(period_info_frame, text="Current Period:").pack(side=tk.LEFT)
        self.period_label = ttk.Label(period_info_frame, text="-- minutes", font=('Segoe UI', 9, 'italic'))
        self.period_label.pack(side=tk.LEFT, padx=5)
        
        # Period presets - FIXED: Better visibility presets
        preset_frame = ttk.Frame(period_frame)
        preset_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(preset_frame, text="Quick Presets:").pack(side=tk.LEFT)
        
        preset_5min_btn = ttk.Button(preset_frame, text="5 min (Fast)", 
                                    command=lambda: self.set_orbital_period(5.0))
        preset_5min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_15min_btn = ttk.Button(preset_frame, text="15 min (Good)", 
                                     command=lambda: self.set_orbital_period(15.0))
        preset_15min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_30min_btn = ttk.Button(preset_frame, text="30 min (Best)", 
                                     command=lambda: self.set_orbital_period(30.0))
        preset_30min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_90min_btn = ttk.Button(preset_frame, text="90 min (ISS)", 
                                     command=lambda: self.set_orbital_period(90.0))
        preset_90min_btn.pack(side=tk.LEFT, padx=2)
        
        # Custom period setter
        custom_frame = ttk.Frame(period_frame)
        custom_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(custom_frame, text="Custom Period (min):").pack(side=tk.LEFT)
        self.custom_period_var = tk.DoubleVar(value=30.0)  # FIXED: Default to 30 minutes
        custom_entry = ttk.Entry(custom_frame, textvariable=self.custom_period_var, width=8)
        custom_entry.pack(side=tk.LEFT, padx=5)
        
        custom_btn = ttk.Button(custom_frame, text="Apply", 
                               command=lambda: self.set_orbital_period(self.custom_period_var.get()))
        custom_btn.pack(side=tk.LEFT, padx=2)
        
        # Calculate period button
        calc_btn = ttk.Button(preset_frame, text="Calculate", 
                             command=self.calculate_current_period)
        calc_btn.pack(side=tk.RIGHT, padx=5)
        
        # Add helpful notes - FIXED: Better guidance
        notes_frame = ttk.Frame(orbit_frame)
        notes_frame.pack(fill=tk.X, pady=5)
        note_text = ("FIXED: Default altitude is now 1629km for better visibility.\n"
                    "30-minute orbits provide good target viewing opportunities.\n"
                    "Higher altitudes = slower orbits, better for camera work.")
        ttk.Label(notes_frame, text=note_text, style="Info.TLabel", wraplength=400).pack(anchor=tk.W)
        
        # Update button
        update_btn = ttk.Button(details_frame, text="Update Satellite", 
                               command=self.update_current_satellite)
        update_btn.pack(pady=10)
        
        # Bind orbit parameter changes to update period display
        for var in self.orbit_vars.values():
            var.trace('w', self.on_orbit_param_changed)
        
        # Select the first satellite
        if self.satellites:
            self.satellite_listbox.selection_set(0)
            self.on_satellite_selected(None)
    
    def calculate_orbital_period(self, semi_major_axis_km):
        """
        Calculate orbital period for given semi-major axis
        
        Parameters:
        semi_major_axis_km (float): Semi-major axis in kilometers
        
        Returns:
        float: Orbital period in minutes
        """
        # Earth's gravitational parameter (m³/s²)
        mu_earth = 3.986004418e14
        
        # Convert semi-major axis to meters
        a_m = semi_major_axis_km * 1000
        
        # Calculate period using Kepler's third law: T = 2π√(a³/μ)
        period_sec = 2 * np.pi * np.sqrt((a_m**3) / mu_earth)
        period_min = period_sec / 60.0
        
        return period_min
    
    def calculate_semi_major_axis_for_period(self, period_minutes):
        """
        Calculate semi-major axis for desired orbital period
        
        Parameters:
        period_minutes (float): Desired orbital period in minutes
        
        Returns:
        float: Semi-major axis in kilometers
        """
        # Earth's gravitational parameter (m³/s²)
        mu_earth = 3.986004418e14
        
        # Convert period to seconds
        period_sec = period_minutes * 60.0
        
        # Calculate semi-major axis using Kepler's third law: a = ∛(μT²/4π²)
        a_m = ((mu_earth * (period_sec**2)) / (4 * (np.pi**2)))**(1/3)
        a_km = a_m / 1000.0
        
        return a_km
    
    def set_orbital_period(self, period_minutes):
        """Set orbital period by adjusting semi-major axis"""
        try:
            # Calculate required semi-major axis
            required_a = self.calculate_semi_major_axis_for_period(period_minutes)
            
            # Check if the altitude is reasonable (above Earth's surface)
            altitude = required_a - 6371.0  # Earth radius
            if altitude < 200:  # Less than 200 km altitude
                messagebox.showwarning("Low Altitude Warning", 
                                     f"Warning: This orbit has altitude of {altitude:.1f} km.\n"
                                     f"Very low orbits may be unrealistic due to atmospheric drag.")
            
            # Update the semi-major axis
            self.orbit_vars["a"].set(required_a)
            
            # Update period display
            self.calculate_current_period()
            
            self.add_log(f"Set orbital period to {period_minutes:.1f} minutes (altitude: {altitude:.1f} km)")
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not set orbital period: {e}")
    
    def calculate_current_period(self):
        """Calculate and display the current orbital period"""
        try:
            a = self.orbit_vars["a"].get()
            period = self.calculate_orbital_period(a)
            altitude = a - 6371.0  # Earth radius
            
            self.period_label.config(text=f"{period:.1f} min (alt: {altitude:.0f} km)")
            
        except Exception as e:
            self.period_label.config(text="Error calculating")
    
    def on_orbit_param_changed(self, *args):
        """Called when orbit parameters change"""
        # Update period display when semi-major axis changes
        self.calculate_current_period()
            
    def update_satellite_listbox(self):
        """FIXED: Simplified satellite display names"""
        self.satellite_listbox.delete(0, tk.END)
        
        show_details = self.show_details_var.get()
        
        for i, sat in enumerate(self.satellites):
            # FIXED: Start with simple name
            display_name = sat["name"]
            
            if show_details:
                # Add status indicators only if details are enabled
                if sat["fault"]["enabled"]:
                    display_name += " [FAULT]"
                    
                if sat["camera"]["enabled"]:
                    display_name += " [CAMERA]"
                    
                # Add target count
                target_count = sum(1 for target in self.parent_app.targets 
                                  if sat["name"] in target["assigned_to"])
                if target_count > 0:
                    display_name += f" [{target_count}T]"
                
                # Add period info
                try:
                    period = self.calculate_orbital_period(sat["orbit"]["a"])
                    display_name += f" ({period:.0f}min)"
                except:
                    pass
                    
            self.satellite_listbox.insert(tk.END, display_name)
            
    def on_satellite_selected(self, event):
        """Handle satellite selection event"""
        selection = self.satellite_listbox.curselection()
        if selection:
            index = selection[0]
            self.current_satellite_index.set(index)
            self.load_satellite_details(index)
            
    def load_satellite_details(self, index):
        """Load satellite details into UI fields"""
        if 0 <= index < len(self.satellites):
            sat = self.satellites[index]
            self.sat_name_var.set(sat["name"])
            
            # Load orbit parameters
            for param, var in self.orbit_vars.items():
                var.set(sat["orbit"][param])
            
            # Update period display
            self.calculate_current_period()
                
    def add_new_satellite(self):
        """Add a new satellite to the constellation"""
        # Create a default name based on number of satellites
        name = f"Satellite{len(self.satellites) + 1}"
        
        # Calculate appropriate true anomaly to evenly space satellites
        true_anomaly = 0.0
        if self.satellites:
            # If satellites exist, calculate a true anomaly that spaces them evenly
            true_anomaly = len(self.satellites) * (360.0 / (len(self.satellites) + 1))
        
        # FIXED: Better orbital parameters for visibility
        new_satellite = {
            "name": name,
            "orbit": {
                "a": 8000.0,   # FIXED: Higher altitude for better visibility (1629 km altitude)
                "e": 0.05,     # Small eccentricity for stable orbit
                "i": 55.0,     # Inclination in degrees (good for coverage)
                "Omega": 45.0, # Right ascension of ascending node
                "omega": 30.0, # Argument of periapsis
                "f": true_anomaly  # True anomaly - evenly spaced
            },
            "fault": {
                "type": "friction",
                "magnitude": 0.0005,
                "wheel": 3,
                "time": 10.0,  # 10 minutes fault injection
                "enabled": False,
                "periodic": {
                    "enabled": False,
                    "interval": 360,  # 6 minutes
                    "magnitude": 0.1,
                    "wheel": 1
                }
            },
            "camera": {
                "position": [0.0, 0.0, 5.0],  # FIXED: Higher camera position for better view
                "fov": 70.0,  # Field of view in degrees
                "enabled": False  # Disable by default
            },
            "targets": []  # Assigned targets
        }
        
        self.satellites.append(new_satellite)
        
        # Update the UI
        self.update_satellite_listbox()
        self.satellite_listbox.selection_clear(0, tk.END)
        self.satellite_listbox.selection_set(len(self.satellites) - 1)
        self.on_satellite_selected(None)
        
        # Update all satellite dropdowns
        self.parent_app.update_satellite_dropdowns()
        
        self.parent_app.add_log(f"Added new satellite: {name}")
        return new_satellite
        
    def add_multiple_satellites(self):
        """Add multiple satellites to form a constellation"""
        num_satellites = self.num_satellites_var.get()
        
        # Check if the number is valid
        if num_satellites < 2 or num_satellites > 8:
            messagebox.showwarning("Invalid Input", 
                                  "Please enter a number between 2 and 8 satellites.")
            return
            
        # Clear existing satellites if user confirms
        confirm = messagebox.askyesno("Confirm Constellation Creation", 
                                    f"This will clear existing satellites and create a new constellation with {num_satellites} satellites. Continue?")
        if not confirm:
            return
            
        # Clear existing satellites
        self.satellites.clear()
        
        # FIXED: Better constellation with improved visibility
        for i in range(num_satellites):
            name = f"Satellite{i+1}"
            true_anomaly = i * (360.0 / num_satellites)  # Even spacing
            
            # Create satellite with better orbital parameters
            satellite = {
                "name": name,
                "orbit": {
                    "a": 8000.0,   # FIXED: Higher altitude for better visibility
                    "e": 0.05,     # Small eccentricity for stable orbit
                    "i": 55.0,     # Inclination in degrees
                    "Omega": 45.0, # Right ascension of ascending node
                    "omega": 30.0, # Argument of periapsis
                    "f": true_anomaly  # Evenly spaced true anomalies
                },
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 10.0,  # 10 minutes fault injection
                    "enabled": False,  # Don't enable faults by default
                    "periodic": {
                        "enabled": False,
                        "interval": 360,  # 6 minutes
                        "magnitude": 0.1,
                        "wheel": 1
                    }
                },
                "camera": {
                    "position": [0.0, 0.0, 5.0],  # FIXED: Better camera position
                    "fov": 70.0,  # Field of view in degrees
                    "enabled": True if i == 0 else False  # Enable camera only for first satellite
                },
                "targets": []  # Assigned targets
            }
            
            self.satellites.append(satellite)
                
        # Update UI
        self.update_satellite_listbox()
        self.satellite_listbox.selection_clear(0, tk.END)
        self.satellite_listbox.selection_set(0)
        self.on_satellite_selected(None)
        
        # Update all satellite dropdowns
        self.parent_app.update_satellite_dropdowns()
        
        # Update target assignments
        self.parent_app.update_target_assignments()
        
        self.parent_app.add_log(f"Created a new constellation with {num_satellites} satellites (improved visibility)")
        
    def remove_satellite(self):
        """Remove the selected satellite"""
        selection = self.satellite_listbox.curselection()
        if selection:
            index = selection[0]
            sat_name = self.satellites[index]["name"]
            
            # Ask for confirmation
            confirm = messagebox.askyesno("Confirm Removal", 
                                         f"Are you sure you want to remove {sat_name}?")
            if confirm:
                # Remove satellite
                self.satellites.pop(index)
                
                # Update UI
                self.update_satellite_listbox()
                if self.satellites:
                    self.satellite_listbox.selection_set(0)
                    self.on_satellite_selected(None)
                
                # Update dropdowns
                self.parent_app.update_satellite_dropdowns()
                
                # Update target assignments
                for target in self.parent_app.targets:
                    if sat_name in target["assigned_to"]:
                        target["assigned_to"].remove(sat_name)
                
                self.parent_app.update_target_assignments()
                self.parent_app.add_log(f"Removed satellite: {sat_name}")
            
    def update_current_satellite(self):
        """Update the current satellite with UI values"""
        selection = self.satellite_listbox.curselection()
        if selection:
            index = selection[0]
            old_name = self.satellites[index]["name"]
            
            # Update with new values
            self.satellites[index]["name"] = self.sat_name_var.get()
            
            # Update orbit parameters
            for param, var in self.orbit_vars.items():
                self.satellites[index]["orbit"][param] = var.get()
                
            # Update UI
            self.update_satellite_listbox()
            self.satellite_listbox.selection_set(index)
            
            # Update target assignments if name changed
            if old_name != self.satellites[index]["name"]:
                for target in self.parent_app.targets:
                    if old_name in target["assigned_to"]:
                        target["assigned_to"].remove(old_name)
                        target["assigned_to"].append(self.satellites[index]["name"])
            
            # Update dropdowns
            self.parent_app.update_satellite_dropdowns()
            
            self.parent_app.add_log(f"Updated satellite: {self.satellites[index]['name']}")