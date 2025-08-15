#!/usr/bin/env python
"""
constellation_tab.py

Implements the Constellation Management tab with multiple orbit support.
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np

# Import base tab - handle import error gracefully
try:
    from .base_tab import BaseTab
except ImportError:
    # Fallback BaseTab if import fails
    class BaseTab:
        def __init__(self, parent_app, parent_frame):
            self.parent_app = parent_app
            self.parent_frame = parent_frame
        
        def add_help_button(self, parent, title, topic=None, command=None):
            pass
        
        def add_log(self, message):
            if hasattr(self.parent_app, 'add_log'):
                self.parent_app.add_log(message)

class ConstellationTab(BaseTab):
    """Constellation Management tab with multiple orbit support"""
    
    def __init__(self, parent_app, parent_frame):
        """Initialize the constellation tab"""
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        self.current_satellite_index = parent_app.current_satellite_index
        
        # Orbit configurations - support multiple orbits
        self.orbit_configurations = [
            {
                "name": "Default Orbit",
                "altitude": 600,  # km above Earth
                "inclination": 53.0,
                "satellites": []
            },
            {
                "name": "MEO Navigation", 
                "altitude": 1200,  # km above Earth
                "inclination": 55.0,
                "satellites": []
            },
            {
                "name": "High Coverage",
                "altitude": 2000,  # km above Earth
                "inclination": 98.0,  # Near polar
                "satellites": []
            }
        ]

        
        # Create the tab UI
        self.create_tab_ui()
        
    def create_tab_ui(self):
        """Create the Constellation Management tab UI with orbit support"""
        # Main container
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Orbit Management Section
        orbit_frame = ttk.LabelFrame(main_container, text="Orbit Management", padding=10)
        orbit_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Orbit selection and creation
        orbit_controls = ttk.Frame(orbit_frame)
        orbit_controls.pack(fill=tk.X, pady=5)
        
        ttk.Label(orbit_controls, text="Select Orbit:").pack(side=tk.LEFT)
        
        self.current_orbit_var = tk.StringVar()
        self.orbit_combo = ttk.Combobox(orbit_controls, textvariable=self.current_orbit_var)
        self.orbit_combo.pack(side=tk.LEFT, padx=5)
        self.update_orbit_combo()
        
        # Orbit buttons
        ttk.Button(orbit_controls, text="New Orbit", 
                  command=self.create_new_orbit).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(orbit_controls, text="Delete Orbit", 
                  command=self.delete_orbit).pack(side=tk.LEFT, padx=5)
        
        # Quick constellation buttons
        quick_frame = ttk.Frame(orbit_frame)
        quick_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(quick_frame, text="Quick Setup:").pack(side=tk.LEFT)
        
        ttk.Button(quick_frame, text="4-Sat Default", 
                command=lambda: self.create_constellation(4, "Default Orbit")).pack(side=tk.LEFT, padx=2)

        ttk.Button(quick_frame, text="6-Sat MEO", 
                command=lambda: self.create_constellation(6, "MEO Navigation")).pack(side=tk.LEFT, padx=2)

        ttk.Button(quick_frame, text="2-Sat High", 
                command=lambda: self.create_constellation(2, "High Coverage")).pack(side=tk.LEFT, padx=2)

        
        # Split main area
        split_frame = ttk.Frame(main_container)
        split_frame.pack(fill=tk.BOTH, expand=True)
        
        left_frame = ttk.Frame(split_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
        
        right_frame = ttk.Frame(split_frame)
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
        
        remove_btn = ttk.Button(btn_frame, text="Remove Satellite", 
                                command=self.remove_satellite)
        remove_btn.pack(side=tk.LEFT, padx=5)
        
        # Satellite count display
        self.sat_count_label = ttk.Label(btn_frame, text="Total: 0 satellites")
        self.sat_count_label.pack(side=tk.RIGHT, padx=5)
        
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
        
        # Orbit assignment
        orbit_assign_frame = ttk.Frame(details_frame)
        orbit_assign_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(orbit_assign_frame, text="Assigned Orbit:").pack(side=tk.LEFT)
        self.sat_orbit_var = tk.StringVar()
        orbit_assign_combo = ttk.Combobox(orbit_assign_frame, textvariable=self.sat_orbit_var)
        orbit_assign_combo.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        orbit_assign_combo['values'] = [orbit['name'] for orbit in self.orbit_configurations]
        
        # Orbit parameters (read-only display)
        orbit_info_frame = ttk.LabelFrame(details_frame, text="Orbit Parameters", padding=10)
        orbit_info_frame.pack(fill=tk.X, pady=5)
        
        self.orbit_info_text = tk.Text(orbit_info_frame, height=6, state="disabled", font=('Consolas', 9))
        self.orbit_info_text.pack(fill=tk.X)
        
        # Position in orbit
        position_frame = ttk.Frame(details_frame)
        position_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(position_frame, text="Position in Orbit (degrees):").pack(side=tk.LEFT)
        self.true_anomaly_var = tk.DoubleVar(value=0.0)
        position_entry = ttk.Entry(position_frame, textvariable=self.true_anomaly_var, width=10)
        position_entry.pack(side=tk.LEFT, padx=5)
        
        # Coverage info
        coverage_frame = ttk.LabelFrame(details_frame, text="Coverage Information", padding=5)
        coverage_frame.pack(fill=tk.X, pady=5)
        
        self.coverage_info_label = ttk.Label(coverage_frame, text="Select satellite to view coverage", 
                                           style="Info.TLabel")
        self.coverage_info_label.pack(fill=tk.X)
        
        # Update button
        update_btn = ttk.Button(details_frame, text="Update Satellite", 
                               command=self.update_current_satellite)
        update_btn.pack(pady=10)
        
        # Select the first satellite if available
        if self.satellites:
            self.satellite_listbox.selection_set(0)
            self.on_satellite_selected(None)
    
    def update_orbit_combo(self):
        """Update orbit combination dropdown"""
        orbit_names = [orbit['name'] for orbit in self.orbit_configurations]
        self.orbit_combo['values'] = orbit_names
        if orbit_names:
            self.orbit_combo.current(0)
    
    def create_new_orbit(self):
        """Create a new orbit configuration"""
        # Simple dialog for orbit parameters
        dialog = tk.Toplevel(self.parent_app.root)
        dialog.title("Create New Orbit")
        dialog.geometry("500x400")
        dialog.transient(self.parent_app.root)
        dialog.grab_set()
        
        frame = ttk.Frame(dialog, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Orbit name
        ttk.Label(frame, text="Orbit Name:").pack(anchor=tk.W)
        name_var = tk.StringVar(value=f"Custom Orbit {len(self.orbit_configurations)+1}")
        ttk.Entry(frame, textvariable=name_var).pack(fill=tk.X, pady=5)
        
        # Altitude
        ttk.Label(frame, text="Altitude (km):").pack(anchor=tk.W, pady=(10,0))
        alt_var = tk.DoubleVar(value=1000)
        ttk.Entry(frame, textvariable=alt_var).pack(fill=tk.X, pady=5)
        
        # Inclination
        ttk.Label(frame, text="Inclination (degrees):").pack(anchor=tk.W, pady=(10,0))
        inc_var = tk.DoubleVar(value=60.0)
        ttk.Entry(frame, textvariable=inc_var).pack(fill=tk.X, pady=5)

        
        # Number of satellites
        ttk.Label(frame, text="Number of Satellites:").pack(anchor=tk.W, pady=(10,0))
        num_sats_var = tk.IntVar(value=4)
        ttk.Spinbox(frame, from_=1, to=20, textvariable=num_sats_var).pack(fill=tk.X, pady=5)


        
        def create_orbit():
            new_orbit = {
                "name": name_var.get(),
                "altitude": alt_var.get(),
                "inclination": inc_var.get(),
                "satellites": []
            }
            self.orbit_configurations.append(new_orbit)
            self.update_orbit_combo()

                # Select the newly created orbit
            self.current_orbit_var.set(new_orbit["name"])
            self.orbit_combo.set(new_orbit["name"])
        
            dialog.destroy()
            self.add_log(f"Created new orbit: {new_orbit['name']} at {new_orbit['altitude']}km")
        
            num_sats = num_sats_var.get()
            self.create_constellation(num_sats, new_orbit["name"])
            self.update_satellite_listbox()


        ttk.Button(frame, text="Create Orbit", command=create_orbit).pack(pady=20)
        ttk.Button(frame, text="Cancel", command=dialog.destroy).pack()
    
    def delete_orbit(self):
        """Delete selected orbit configuration"""
        current_orbit = self.current_orbit_var.get()
        if not current_orbit:
            return
            
        # Check if orbit has satellites
        orbit_config = next((o for o in self.orbit_configurations if o['name'] == current_orbit), None)
        if orbit_config and orbit_config['satellites']:
            messagebox.showwarning("Cannot Delete", 
                                 f"Orbit '{current_orbit}' contains {len(orbit_config['satellites'])} satellites. "
                                 f"Remove satellites first.")
            return
        
        confirm = messagebox.askyesno("Delete Orbit", f"Delete orbit '{current_orbit}'?")
        if confirm:
            self.orbit_configurations = [o for o in self.orbit_configurations if o['name'] != current_orbit]
            self.update_orbit_combo()
            
            self.add_log(f"Deleted orbit: {current_orbit}")


    def create_constellation(self, num_satellites, orbit_name):
        """Create a constellation of satellites in specified orbit"""
        # Find the orbit configuration
        orbit_config = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
        if not orbit_config:
            messagebox.showerror("Error", f"Orbit '{orbit_name}' not found")
            return
        
        confirm = messagebox.askyesno("Create Constellation", 
                                    f"Create {num_satellites} satellites in {orbit_name} orbit?\n"
                                    f"Altitude: {orbit_config['altitude']}km")
        if not confirm:
            return
        
        # Calculate orbital parameters
        altitude_km = orbit_config['altitude']
        semi_major_axis = 6371 + altitude_km  # Earth radius + altitude
        inclination = orbit_config['inclination']
        
        # Fixed RAAN for each orbit type to ensure separation between different orbits
        orbit_raan_map = {
            "Default Orbit": 0.0,
            "MEO Navigation": 60.0,
            "High Coverage": 120.0
        }
        raan_degrees = orbit_raan_map.get(orbit_name, 0.0)

        # Create satellites
        for i in range(num_satellites):
            name = f"{orbit_name.replace(' ', '')}_Sat{i+1}"
            true_anomaly = i * (360.0 / num_satellites)  # Even spacing
            
            # Create satellite with proper orbital parameters
            satellite = {
                "name": name,
                "orbit": {
                    "a": semi_major_axis,  # Semi-major axis in km
                    "e": 0.01,  # Small eccentricity for stable orbit
                    "i": inclination,  # Inclination from orbit config
                    "Omega": raan_degrees,  # Consistent RAAN for this orbit type
                    "omega": 0.0,  # Argument of periapsis
                    "f": true_anomaly  # True anomaly - evenly spaced
                },
                "fault": {
                    "type": "friction",
                    "magnitude": 0.0005,
                    "wheel": 3,
                    "time": 15.0,  # 15 minutes fault injection
                    "enabled": False,
                    "periodic": {
                        "enabled": False,
                        "interval": 360,
                        "magnitude": 0.1,
                        "wheel": 1
                    }
                },
                "camera": {
                    "position": [0.0, 0.0, 15.0],  # Higher camera position for better target view
                    "fov": 80.0,  # Wider FOV
                    "enabled": True if i == 0 else False  # Enable camera only for first satellite
                },
                "targets": [],
                "orbit_name": orbit_name  # Track which orbit this satellite belongs to
            }
            
            self.satellites.append(satellite)
            orbit_config['satellites'].append(name)
        
        # Update UI
        self.update_satellite_listbox()
        self.parent_app.update_satellite_dropdowns()
        
        self.add_log(f"Created {num_satellites}-satellite constellation in {orbit_name} orbit")
        self.add_log(f"Satellites at {altitude_km}km altitude, RAAN: {raan_degrees}°")
        
    def calculate_orbital_period(self, semi_major_axis_km):
        """Calculate orbital period for given semi-major axis"""
        mu_earth = 3.986004418e14  # Earth's gravitational parameter (m³/s²)
        a_m = semi_major_axis_km * 1000  # Convert to meters
        period_sec = 2 * np.pi * np.sqrt((a_m**3) / mu_earth)
        return period_sec / 60.0  # Return in minutes
    
    def calculate_coverage_info(self, satellite):
        """Calculate coverage information for satellite"""
        try:
            altitude = satellite["orbit"]["a"] - 6371  # Altitude above Earth
            
            # Calculate coverage radius (simplified)
            earth_radius = 6371  # km
            satellite_height = altitude
            
            # Maximum range for communication/observation (simplified calculation)
            horizon_distance = np.sqrt(2 * earth_radius * satellite_height + satellite_height**2)
            
            # Coverage area (approximate)
            coverage_area = np.pi * horizon_distance**2
            
            # Orbital period
            period = self.calculate_orbital_period(satellite["orbit"]["a"])
            
            return {
                "altitude": altitude,
                "horizon_distance": horizon_distance,
                "coverage_area": coverage_area,
                "period": period
            }
        except:
            return None
            
    def update_satellite_listbox(self):
        """Update the satellite listbox with orbit information"""
        self.satellite_listbox.delete(0, tk.END)
        
        orbit_counts = {}
        
        for i, sat in enumerate(self.satellites):
            orbit_name = sat.get("orbit_name", "Unknown")
            
            # Count satellites per orbit
            if orbit_name not in orbit_counts:
                orbit_counts[orbit_name] = 0
            orbit_counts[orbit_name] += 1
            
            # Display format: "Satellite Name [Orbit] (altitude)"
            altitude = sat["orbit"]["a"] - 6371
            display_name = f"{sat['name']} [{orbit_name}] ({altitude:.0f}km)"
            
            self.satellite_listbox.insert(tk.END, display_name)
        
        # Update satellite count
        total_sats = len(self.satellites)
        orbit_summary = ", ".join([f"{count} in {orbit}" for orbit, count in orbit_counts.items()])
        self.sat_count_label.config(text=f"Total: {total_sats} satellites ({orbit_summary})")
            
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
            self.sat_orbit_var.set(sat.get("orbit_name", "Unknown"))
            self.true_anomaly_var.set(sat["orbit"]["f"])
            
            # Update orbit info display
            self.orbit_info_text.config(state="normal")
            self.orbit_info_text.delete(1.0, tk.END)
            
            orbit_info = f"Semi-major axis: {sat['orbit']['a']:.1f} km\n"
            orbit_info += f"Altitude: {sat['orbit']['a'] - 6371:.1f} km\n"
            orbit_info += f"Eccentricity: {sat['orbit']['e']:.3f}\n"
            orbit_info += f"Inclination: {sat['orbit']['i']:.1f}°\n"
            orbit_info += f"True Anomaly: {sat['orbit']['f']:.1f}°\n"
            
            period = self.calculate_orbital_period(sat['orbit']['a'])
            orbit_info += f"Orbital Period: {period:.1f} minutes"
            
            self.orbit_info_text.insert(tk.END, orbit_info)
            self.orbit_info_text.config(state="disabled")
            
            # Update coverage info
            coverage = self.calculate_coverage_info(sat)
            if coverage:
                coverage_text = f"Altitude: {coverage['altitude']:.0f}km, "
                coverage_text += f"Horizon: {coverage['horizon_distance']:.0f}km, "
                coverage_text += f"Period: {coverage['period']:.1f}min"
                self.coverage_info_label.config(text=coverage_text)
            else:
                self.coverage_info_label.config(text="Coverage calculation error")
                


    def add_new_satellite(self):
        """Add a new satellite to selected orbit"""
        current_orbit = self.current_orbit_var.get()
        if not current_orbit:
            messagebox.showwarning("No Orbit Selected", "Please select an orbit first")
            return
        
        # Find orbit configuration
        orbit_config = next((o for o in self.orbit_configurations if o['name'] == current_orbit), None)
        if not orbit_config:
            messagebox.showerror("Error", f"Orbit '{current_orbit}' not found")
            return
        
        # Create satellite name
        sat_count_in_orbit = len(orbit_config['satellites'])
        name = f"{current_orbit.replace(' ', '')}_Sat{sat_count_in_orbit + 1}"
        
        # Calculate position in orbit
        if orbit_config['satellites']:
            # Space evenly among existing satellites
            total_sats = len(orbit_config['satellites']) + 1
            true_anomaly = sat_count_in_orbit * (360.0 / total_sats)
        else:
            true_anomaly = 0.0
        
        # Calculate orbital parameters from orbit config
        altitude_km = orbit_config['altitude']
        semi_major_axis = 6371 + altitude_km
        inclination = orbit_config['inclination']
        
        # Fixed RAAN for each orbit type
        orbit_raan_map = {
            "Default Orbit": 0.0,
            "MEO Navigation": 60.0,
            "High Coverage": 120.0
        }
        raan_degrees = orbit_raan_map.get(current_orbit, 0.0)
        
        # Create new satellite
        new_satellite = {
            "name": name,
            "orbit": {
                "a": semi_major_axis,
                "e": 0.01,
                "i": inclination,
                "Omega": raan_degrees,  # Consistent RAAN for orbit type
                "omega": 0.0,
                "f": true_anomaly
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
                "position": [0.0, 0.0, 15.0],  # Higher position for target viewing
                "fov": 80.0,
                "enabled": False
            },
            "targets": [],
            "orbit_name": current_orbit
        }
        
        self.satellites.append(new_satellite)
        orbit_config['satellites'].append(name)
        
        # Update UI
        self.update_satellite_listbox()
        self.satellite_listbox.selection_clear(0, tk.END)
        self.satellite_listbox.selection_set(len(self.satellites) - 1)
        self.on_satellite_selected(None)
        
        self.parent_app.update_satellite_dropdowns()
        self.add_log(f"Added satellite: {name} to {current_orbit} orbit (RAAN: {raan_degrees}°)")
        
    def remove_satellite(self):
        """Remove the selected satellite"""
        selection = self.satellite_listbox.curselection()
        if selection:
            index = selection[0]
            sat = self.satellites[index]
            sat_name = sat["name"]
            orbit_name = sat.get("orbit_name", "Unknown")
            
            # Ask for confirmation
            confirm = messagebox.askyesno("Confirm Removal", 
                                         f"Remove {sat_name} from {orbit_name} orbit?")
            if confirm:
                # Remove from orbit configuration
                orbit_config = next((o for o in self.orbit_configurations if o['name'] == orbit_name), None)
                if orbit_config and sat_name in orbit_config['satellites']:
                    orbit_config['satellites'].remove(sat_name)
                
                # Remove satellite
                self.satellites.pop(index)
                
                # Update UI
                self.update_satellite_listbox()
                if self.satellites:
                    self.satellite_listbox.selection_set(0)
                    self.on_satellite_selected(None)
                
                # Update dropdowns and target assignments
                self.parent_app.update_satellite_dropdowns()
                
                for target in self.parent_app.targets:
                    if sat_name in target["assigned_to"]:
                        target["assigned_to"].remove(sat_name)
                
                self.parent_app.update_target_assignments()
                self.add_log(f"Removed satellite: {sat_name} from {orbit_name}")
            
    def update_current_satellite(self):
        """Update the current satellite with UI values"""
        selection = self.satellite_listbox.curselection()
        if selection:
            index = selection[0]
            old_name = self.satellites[index]["name"]
            new_orbit = self.sat_orbit_var.get()
            
            # Update name and orbit assignment
            self.satellites[index]["name"] = self.sat_name_var.get()
            self.satellites[index]["orbit"]["f"] = self.true_anomaly_var.get()
            
            # Handle orbit change
            old_orbit = self.satellites[index].get("orbit_name", "Unknown") 
            if new_orbit != old_orbit:
                # Remove from old orbit
                old_orbit_config = next((o for o in self.orbit_configurations if o['name'] == old_orbit), None)
                if old_orbit_config and old_name in old_orbit_config['satellites']:
                    old_orbit_config['satellites'].remove(old_name)
                
                # Add to new orbit and update parameters
                new_orbit_config = next((o for o in self.orbit_configurations if o['name'] == new_orbit), None)
                if new_orbit_config:
                    new_orbit_config['satellites'].append(self.satellites[index]["name"])
                    
                    # Update orbital parameters to match new orbit
                    altitude_km = new_orbit_config['altitude']
                    semi_major_axis = 6371 + altitude_km
                    inclination = new_orbit_config['inclination']
                    
                    self.satellites[index]["orbit"]["a"] = semi_major_axis
                    self.satellites[index]["orbit"]["i"] = inclination
                    self.satellites[index]["orbit_name"] = new_orbit
                    
                    self.add_log(f"Moved {self.satellites[index]['name']} to {new_orbit} orbit")
            
            # Update target assignments if name changed
            if old_name != self.satellites[index]["name"]:
                for target in self.parent_app.targets:
                    if old_name in target["assigned_to"]:
                        target["assigned_to"].remove(old_name)
                        target["assigned_to"].append(self.satellites[index]["name"])
            
            # Update UI
            self.update_satellite_listbox()
            self.satellite_listbox.selection_set(index)
            self.load_satellite_details(index)
            
            self.parent_app.update_satellite_dropdowns()
            self.add_log(f"Updated satellite: {self.satellites[index]['name']}")

    def add_log(self, message):
        """Add a message to the log"""
        if hasattr(self.parent_app, 'add_log'):
            self.parent_app.add_log(message)
        else:
            print(f"Log: {message}")