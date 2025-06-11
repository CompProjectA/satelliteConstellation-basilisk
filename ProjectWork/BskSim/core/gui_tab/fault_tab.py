#!/usr/bin/env python
"""
fault_tab.py

Implements the Fault Configuration tab for the Spacecraft Constellation Fault Simulator.
Corrected to show wheel numbers as 1-4 consistently throughout the interface.
"""
import tkinter as tk
from tkinter import ttk, messagebox
from .base_tab import BaseTab

class FaultTab(BaseTab):
    """Fault Configuration tab implementation"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the fault tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        
        # Create the tab UI
        self.create_tab_ui()
        
        # Bind fault type change after UI is created
        self.fault_type_var.trace('w', self.on_fault_type_change)

    def get_active_satellite_index(self):
        """
        Get the index of the currently selected satellite in the fault tab
        
        Returns:
        int: Index of the active satellite, or None if no satellite is selected
        """
        sat_name = self.fault_satellite_var.get()
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                return i
        return None
    
    def set_active_satellite(self, index):
        """
        Set the active satellite by index
        
        Parameters:
        index (int): Index of the satellite to set as active
        """
        if 0 <= index < len(self.satellites):
            self.fault_satellite_combo.current(index)
            self.load_fault_config(index)
            
    def on_fault_type_change(self, *args):
        """Handle fault type change to update default magnitudes"""
        fault_type = self.fault_type_var.get()
        
        # Default magnitudes per fault type
        DEFAULT_FAULT_MAGNITUDES = {
            "friction": 0.0005,    # N⋅m - constant Coulomb friction
            "power_limit": 0.5,    # W - power limitation
            "encoder": 20.0,       # % - encoder error percentage
            "battery": 50.0        # W - additional power drain
        }
        
        # If the current magnitude is still the old default, update it
        current_magnitude = self.fault_mag.get()
        if abs(current_magnitude - 0.0005) < 1e-6 or fault_type != "friction":
            if fault_type in DEFAULT_FAULT_MAGNITUDES:
                new_magnitude = DEFAULT_FAULT_MAGNITUDES[fault_type]
                self.fault_mag.set(new_magnitude)
                
                # Log the change if not loading from existing config
                if hasattr(self, '_loading_config') and not self._loading_config:
                    self.parent_app.add_log(f"Updated {fault_type} fault magnitude to {new_magnitude}")
                    
    def create_tab_ui(self):
        """Create the Fault Configuration tab UI"""
        # Create a notebook for sub-tabs
        self.fault_notebook = ttk.Notebook(self.parent_frame)
        self.fault_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create frames for sub-tabs
        config_frame = ttk.Frame(self.fault_notebook)
        summary_frame = ttk.Frame(self.fault_notebook)
        
        # Add sub-tabs
        self.fault_notebook.add(config_frame, text="Fault Configuration")
        self.fault_notebook.add(summary_frame, text="Fault Summary")
        
        # Create fault configuration content
        self._create_fault_config_tab(config_frame)
        
        # Create fault summary content
        self._create_fault_summary_tab(summary_frame)
        
        # Load fault configuration for the first satellite
        if self.satellites:
            self.load_fault_config(0)
            self.update_fault_summary()

    def _create_fault_config_tab(self, parent):
        """Create the fault configuration sub-tab content"""
        # Satellite selection at the top
        select_frame = ttk.Frame(parent)
        select_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(select_frame, text="Select Satellite:").pack(side=tk.LEFT)
        
        self.fault_satellite_var = tk.StringVar()
        self.fault_satellite_combo = ttk.Combobox(select_frame, textvariable=self.fault_satellite_var)
        self.fault_satellite_combo.pack(side=tk.LEFT, padx=5)
        self.update_satellite_dropdown()
        
        # Bind selection change event
        self.fault_satellite_combo.bind('<<ComboboxSelected>>', self.on_fault_satellite_changed)
        
        # Current fault status
        status_frame = ttk.Frame(select_frame)
        status_frame.pack(side=tk.RIGHT, padx=5)
        
        self.fault_status_label = ttk.Label(status_frame, text="No fault enabled")
        self.fault_status_label.pack(side=tk.RIGHT)
        
        # Main fault configuration
        fault_config_frame = ttk.LabelFrame(parent, text="Fault Configuration", padding=10)
        fault_config_frame.pack(fill=tk.BOTH, expand=True)
        
        # Enable fault checkbox
        self.fault_enabled_var = tk.BooleanVar(value=False)
        enable_check = ttk.Checkbutton(fault_config_frame, text="Enable Fault", 
                                    variable=self.fault_enabled_var,
                                    command=self.update_fault_config)
        enable_check.pack(anchor=tk.W, pady=5)
        
        # Fault type selection
        type_frame = ttk.Frame(fault_config_frame)
        type_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(type_frame, text="Fault Type:").pack(side=tk.LEFT)
        
        self.fault_type_var = tk.StringVar(value="friction")
        fault_types = [
            ("Friction", "friction"),
            ("Power Limit", "power_limit"),
            ("Encoder", "encoder"),
            ("Battery", "battery")
        ]
        
        fault_type_frame = ttk.Frame(type_frame)
        fault_type_frame.pack(side=tk.LEFT, padx=5)
        
        for text, value in fault_types:
            ttk.Radiobutton(fault_type_frame, text=text, value=value, 
                        variable=self.fault_type_var,
                        command=self.update_fault_config).pack(side=tk.LEFT, padx=5)
        
        # Fault description 
        self.fault_description_frame = ttk.LabelFrame(fault_config_frame, text="Fault Description", padding=10)
        self.fault_description_frame.pack(fill=tk.X, pady=5)
        
        self.fault_description_label = ttk.Label(self.fault_description_frame, 
                                            text="Select a fault type to see its description",
                                            wraplength=500, justify=tk.LEFT)
        self.fault_description_label.pack(fill=tk.X, pady=5)
        
        # Fault parameters
        self.params_frame = ttk.LabelFrame(fault_config_frame, text="Fault Parameters", padding=10)
        self.params_frame.pack(fill=tk.X, pady=5)
        
        # Initialize variables for fault parameters
        self.fault_mag = tk.DoubleVar(value=0.0005)
        self.fault_time = tk.DoubleVar(value=10.0)
        self.fault_wheel = tk.IntVar(value=4)  # Default to RW4 (displayed as 4, stored as 3)
        
        # Create parameter widgets (updated based on fault type)
        self.create_parameter_widgets()
        
        # Periodic fault section
        periodic_frame = ttk.LabelFrame(fault_config_frame, text="Periodic Fault", padding=10)
        periodic_frame.pack(fill=tk.X, pady=5)
        
        # Enable periodic fault
        self.periodic_enabled_var = tk.BooleanVar(value=False)
        periodic_check = ttk.Checkbutton(periodic_frame, text="Enable Periodic Fault", 
                                        variable=self.periodic_enabled_var,
                                        command=self.toggle_periodic)
        periodic_check.pack(anchor=tk.W, pady=5)
        
        # Periodic parameters
        periodic_params_frame = ttk.Frame(periodic_frame)
        periodic_params_frame.pack(fill=tk.X)
        
        # Interval
        interval_frame = ttk.Frame(periodic_params_frame)
        interval_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(interval_frame, text="Interval (sec):").pack(side=tk.LEFT)
        
        self.periodic_interval_var = tk.DoubleVar(value=360.0)
        self.periodic_interval_entry = ttk.Entry(interval_frame, textvariable=self.periodic_interval_var, width=10, state="disabled")
        self.periodic_interval_entry.pack(side=tk.LEFT, padx=5)
        
        # Periodic magnitude
        p_magnitude_frame = ttk.Frame(periodic_params_frame)
        p_magnitude_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(p_magnitude_frame, text="Magnitude:").pack(side=tk.LEFT)
        
        self.periodic_magnitude_var = tk.DoubleVar(value=0.1)
        self.periodic_magnitude_entry = ttk.Entry(p_magnitude_frame, textvariable=self.periodic_magnitude_var, width=10, state="disabled")
        self.periodic_magnitude_entry.pack(side=tk.LEFT, padx=5)
        
        # Periodic wheel
        p_wheel_frame = ttk.Frame(periodic_params_frame)
        p_wheel_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(p_wheel_frame, text="Wheel Number:").pack(side=tk.LEFT)
        
        self.periodic_wheel_var = tk.IntVar(value=1)
        self.periodic_wheel_spinbox = ttk.Spinbox(p_wheel_frame, from_=1, to=4, textvariable=self.periodic_wheel_var, width=5, state="disabled")
        self.periodic_wheel_spinbox.pack(side=tk.LEFT, padx=5)
        
        # Button frame
        button_frame = ttk.Frame(fault_config_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        # Update button
        update_btn = ttk.Button(button_frame, text="Apply Fault Configuration", 
                            command=self.apply_fault_config)
        update_btn.pack(side=tk.LEFT, padx=5)
        
        # Apply to all satellites button
        apply_all_btn = ttk.Button(button_frame, text="Apply to All Satellites", 
                                command=self.apply_to_all_satellites)
        apply_all_btn.pack(side=tk.RIGHT, padx=5)

    def _create_fault_summary_tab(self, parent):
        """Create the fault summary sub-tab content"""
        # Main frame
        main_frame = ttk.Frame(parent)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Title and description
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(title_frame, text="Fault Summary Overview", 
                font=('Segoe UI', 12, 'bold')).pack(side=tk.LEFT)
        
        ttk.Label(title_frame, text="View all configured faults across the constellation", 
                style="Info.TLabel").pack(side=tk.LEFT, padx=(20, 0))
        
        # Summary controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(control_frame, text="Refresh Summary", 
                command=self.update_fault_summary).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(control_frame, text="Clear All Faults", 
                command=self.clear_all_faults).pack(side=tk.LEFT, padx=5)
        
        # Statistics
        stats_frame = ttk.LabelFrame(control_frame, text="Statistics", padding=5)
        stats_frame.pack(side=tk.RIGHT, padx=10)
        
        self.fault_stats_label = ttk.Label(stats_frame, text="Total: 0 satellites | Faults: 0 enabled")
        self.fault_stats_label.pack()
        
        # Create treeview for fault summary
        tree_frame = ttk.LabelFrame(main_frame, text="Satellite Fault Status", padding=10)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        columns = ('Satellite', 'Fault Status', 'Type', 'Wheel', 'Time (min)', 'Magnitude', 'Periodic')
        self.fault_summary_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', height=15)
        
        # Define column headings
        self.fault_summary_tree.heading('Satellite', text='Satellite')
        self.fault_summary_tree.heading('Fault Status', text='Status')
        self.fault_summary_tree.heading('Type', text='Fault Type')
        self.fault_summary_tree.heading('Wheel', text='Wheel')
        self.fault_summary_tree.heading('Time (min)', text='Time')
        self.fault_summary_tree.heading('Magnitude', text='Magnitude')
        self.fault_summary_tree.heading('Periodic', text='Periodic')
        
        # Set column widths
        self.fault_summary_tree.column('Satellite', width=150)
        self.fault_summary_tree.column('Fault Status', width=100)
        self.fault_summary_tree.column('Type', width=120)
        self.fault_summary_tree.column('Wheel', width=80)
        self.fault_summary_tree.column('Time (min)', width=100)
        self.fault_summary_tree.column('Magnitude', width=120)
        self.fault_summary_tree.column('Periodic', width=100)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.fault_summary_tree.yview)
        self.fault_summary_tree.configure(yscrollcommand=scrollbar.set)
        
        self.fault_summary_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def clear_all_faults(self):
        """Clear all faults from all satellites"""
        if not messagebox.askyesno("Clear All Faults", 
                                  "Are you sure you want to disable all faults on all satellites?"):
            return
        
        # Clear faults on all satellites
        for satellite in self.satellites:
            fault = satellite["fault"]
            fault["enabled"] = False
            fault["periodic"]["enabled"] = False
        
        # Update UI
        self.load_fault_config(self.fault_satellite_combo.current())
        self.update_fault_summary()
        
        # Update constellation tab
        try:
            self.parent_app.constellation_tab.update_satellite_listbox()
        except:
            pass
        
        self.parent_app.add_log("Cleared all faults from all satellites")
        messagebox.showinfo("Success", "All faults have been disabled")
            
    def update_satellite_dropdown(self):
        """Update the satellite dropdown with current satellites"""
        self.fault_satellite_combo['values'] = [sat['name'] for sat in self.satellites]
        if self.satellites:
            self.fault_satellite_combo.current(0)
        
            
    def on_fault_satellite_changed(self, event):
        """Handle fault satellite selection change"""
        sat_name = self.fault_satellite_var.get()
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                self.load_fault_config(i)
                break
                
    def load_fault_config(self, index):
        """Load fault configuration for the specified satellite"""
        if 0 <= index < len(self.satellites):
            # Set flag to prevent magnitude change logging when loading
            self._loading_config = True
            
            fault = self.satellites[index]["fault"]
            
            self.fault_enabled_var.set(fault["enabled"])
            self.fault_type_var.set(fault["type"])
            self.fault_mag.set(fault["magnitude"])
            self.fault_wheel.set(fault["wheel"] + 1)  # Convert 0-3 to 1-4 for display
            self.fault_time.set(fault["time"])
            
            periodic = fault["periodic"]
            self.periodic_enabled_var.set(periodic["enabled"])
            self.periodic_interval_var.set(periodic["interval"])
            self.periodic_magnitude_var.set(periodic["magnitude"])
            self.periodic_wheel_var.set(periodic["wheel"] + 1)  # Convert 0-3 to 1-4 for display
            
            self.update_fault_config()
            self.update_fault_status()
            
            # Clear loading flag
            self._loading_config = False
            
    def update_fault_config(self):
        """Update the fault configuration UI based on selections"""
        # Clear previous parameter widgets
        for widget in self.params_frame.winfo_children():
            widget.destroy()
            
        fault_type = self.fault_type_var.get()
        
        # Update the fault description
        self.update_fault_description()
        
        # Create parameter widgets based on fault type
        if fault_type == "friction":
            self.create_friction_params()
        elif fault_type == "power_limit":
            self.create_power_limit_params()
        elif fault_type == "encoder":
            self.create_encoder_params()
        elif fault_type == "battery":
            self.create_battery_params()
        
        # Update periodic fault widgets state
        self.toggle_periodic()
        
        # Update fault status
        self.update_fault_status()
        self.update_fault_summary()
        
    def update_fault_description(self):
        """Update the fault description based on selected fault type"""
        fault_type = self.fault_type_var.get()
        
        descriptions = {
            "friction": "Friction fault increases the constant Coulomb friction in the reaction wheel, simulating mechanical issues like bearing damage. The default Coulomb friction is 0.0005 N⋅m. Higher values cause slower wheel operation and increased temperatures.",
            "power_limit": "Power limit fault restricts the maximum electrical power available to the reaction wheel, simulating power system limitations. Lower magnitude values create more severe effects by limiting available power.",
            "encoder": "Encoder fault causes measurement errors in the reaction wheel speed feedback, leading to attitude control errors. This can result in oscillations and instability. Error is specified as a percentage.",
            "battery": "Battery fault simulates increased power consumption or battery degradation. When battery state of charge falls below 20%, the spacecraft enters safe mode to preserve remaining charge. Magnitude is additional power drain in Watts."
        }
        
        description = descriptions.get(fault_type, "No description available")
        self.fault_description_label.config(text=description)
        
    def update_fault_status(self):
        """Update the fault status indicator"""
        sat_name = self.fault_satellite_var.get()
        
        if self.fault_enabled_var.get():
            fault_type = self.fault_type_var.get().replace('_', ' ').title()
            wheel_num = self.fault_wheel.get()
            self.fault_status_label.config(text=f"Fault enabled: {fault_type} on RW{wheel_num}")
        else:
            self.fault_status_label.config(text="No fault enabled")
        
    def create_parameter_widgets(self):
        """Create initial parameter widgets based on default fault type"""
        self.create_friction_params()
        
    def create_friction_params(self):
        """Create friction fault parameter widgets"""
        # Magnitude
        mag_frame = ttk.Frame(self.params_frame)
        mag_frame.pack(fill=tk.X, pady=2)
        
        self.fault_mag_label = ttk.Label(mag_frame, text="Magnitude (N⋅m):")
        self.fault_mag_label.pack(side=tk.LEFT)
        ttk.Entry(mag_frame, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(mag_frame, text="(additional Coulomb friction torque)", style="Info.TLabel").pack(side=tk.LEFT, padx=5)
        
        # Wheel
        wheel_frame = ttk.Frame(self.params_frame)
        wheel_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(wheel_frame, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(wheel_frame, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)
        
        # Time
        time_frame = ttk.Frame(self.params_frame)
        time_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(time_frame, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)
        
    def create_power_limit_params(self):
        """Create power limit fault parameter widgets"""
        # Power limit
        limit_frame = ttk.Frame(self.params_frame)
        limit_frame.pack(fill=tk.X, pady=2)
        
        self.fault_mag_label = ttk.Label(limit_frame, text="Power Limit (W):")
        self.fault_mag_label.pack(side=tk.LEFT)
        ttk.Entry(limit_frame, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(limit_frame, text="(maximum power in Watts)", style="Info.TLabel").pack(side=tk.LEFT, padx=5)
        
        # Wheel
        wheel_frame = ttk.Frame(self.params_frame)
        wheel_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(wheel_frame, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(wheel_frame, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)
        
        # Time
        time_frame = ttk.Frame(self.params_frame)
        time_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(time_frame, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)
        
    def create_encoder_params(self):
        """Create encoder fault parameter widgets"""
        # Add magnitude field for encoder
        mag_frame = ttk.Frame(self.params_frame)
        mag_frame.pack(fill=tk.X, pady=2)
        
        self.fault_mag_label = ttk.Label(mag_frame, text="Error Magnitude (%):")
        self.fault_mag_label.pack(side=tk.LEFT)
        ttk.Entry(mag_frame, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(mag_frame, text="(encoder error percentage)", style="Info.TLabel").pack(side=tk.LEFT, padx=5)
        
        # Wheel
        wheel_frame = ttk.Frame(self.params_frame)
        wheel_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(wheel_frame, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(wheel_frame, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)
        
        # Time
        time_frame = ttk.Frame(self.params_frame)
        time_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(time_frame, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)
        
        # Info
        ttk.Label(self.params_frame, text="Encoder faults cause measurement errors in wheel speed feedback.", 
                 style="Info.TLabel").pack(pady=5)
        
    def create_battery_params(self):
        """Create battery fault parameter widgets"""
        # Power drain
        drain_frame = ttk.Frame(self.params_frame)
        drain_frame.pack(fill=tk.X, pady=2)
        
        self.fault_mag_label = ttk.Label(drain_frame, text="Power Drain (W):")
        self.fault_mag_label.pack(side=tk.LEFT)
        ttk.Entry(drain_frame, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(drain_frame, text="(additional power drain in Watts)", style="Info.TLabel").pack(side=tk.LEFT, padx=5)
        
        # Time
        time_frame = ttk.Frame(self.params_frame)
        time_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(time_frame, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)
        
        # Info
        ttk.Label(self.params_frame, text="Safe mode activates when battery charge falls below 20%.", 
                 style="Info.TLabel").pack(pady=5)
        
    def toggle_periodic(self):
        """Enable/disable periodic fault entry fields"""
        state = "normal" if self.periodic_enabled_var.get() else "disabled"
        self.periodic_interval_entry.config(state=state)
        self.periodic_magnitude_entry.config(state=state)
        self.periodic_wheel_spinbox.config(state=state)

    def apply_to_all_satellites(self):
        """Apply the current fault configuration to all satellites"""
        # Ask for confirmation
        if not messagebox.askyesno("Confirm", 
                                "Apply the current fault configuration to all satellites?"):
            return
            
        # Get current configuration
        fault_enabled = self.fault_enabled_var.get()
        fault_type = self.fault_type_var.get()
        fault_mag = self.fault_mag.get()
        fault_wheel = self.fault_wheel.get() - 1  # Convert back to 0-based for internal storage
        fault_time = self.fault_time.get()
        
        periodic_enabled = self.periodic_enabled_var.get()
        periodic_interval = self.periodic_interval_var.get()
        periodic_magnitude = self.periodic_magnitude_var.get()
        periodic_wheel = self.periodic_wheel_var.get() - 1  # Convert back to 0-based
        
        # Apply to all satellites
        for satellite in self.satellites:
            fault = satellite["fault"]
            
            fault["enabled"] = fault_enabled
            fault["type"] = fault_type
            fault["magnitude"] = fault_mag
            fault["wheel"] = fault_wheel
            fault["time"] = fault_time
            
            periodic = fault["periodic"]
            periodic["enabled"] = periodic_enabled
            periodic["interval"] = periodic_interval
            periodic["magnitude"] = periodic_magnitude
            periodic["wheel"] = periodic_wheel
        
        # Update constellation tab to show fault indicators
        try:
            self.parent_app.constellation_tab.update_satellite_listbox()
        except:
            pass
            
        self.parent_app.add_log(f"Applied fault configuration to all satellites")
        messagebox.showinfo("Success", "Fault configuration applied to all satellites")
        self.update_fault_summary()
        
    def apply_fault_config(self):
        """Apply the fault configuration to the selected satellite"""
        sat_name = self.fault_satellite_var.get()
        sat_index = -1
        
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                sat_index = i
                break
                
        if sat_index >= 0:
            fault = self.satellites[sat_index]["fault"]
            
            # Get current values from GUI
            fault_enabled = self.fault_enabled_var.get()
            fault_type = self.fault_type_var.get()
            fault_magnitude = self.fault_mag.get()
            fault_wheel = self.fault_wheel.get() - 1  # Convert from 1-4 to 0-3 for internal storage
            fault_time = self.fault_time.get()
            
            # Default magnitudes per fault type
            DEFAULT_FAULT_MAGNITUDES = {
                "friction": 0.0005,    # N⋅m - default Coulomb friction
                "power_limit": 0.5,    # W - realistic power limit
                "encoder": 20.0,       # % - noticeable encoder error
                "battery": 50.0        # W - significant battery drain
            }
            
            # Scale fault magnitude if it's still the default 0.0005
            if abs(fault_magnitude - 0.0005) < 1e-6 and fault_type in DEFAULT_FAULT_MAGNITUDES:
                scaled_magnitude = DEFAULT_FAULT_MAGNITUDES[fault_type]
                if fault_type != "friction":  # Don't log for friction since 0.0005 is correct
                    self.parent_app.add_log(f"Scaling {fault_type} fault magnitude from {fault_magnitude} to {scaled_magnitude}")
                fault_magnitude = scaled_magnitude
                # Update the GUI display to show the new value
                self.fault_mag.set(fault_magnitude)
            
            # Apply fault configuration
            fault["enabled"] = fault_enabled
            fault["type"] = fault_type
            fault["magnitude"] = fault_magnitude
            fault["wheel"] = fault_wheel
            fault["time"] = fault_time
            
            # Apply periodic fault configuration
            periodic = fault["periodic"]
            periodic["enabled"] = self.periodic_enabled_var.get()
            periodic["interval"] = self.periodic_interval_var.get()
            periodic_magnitude = self.periodic_magnitude_var.get()
            periodic_wheel = self.periodic_wheel_var.get() - 1  # Convert from 1-4 to 0-3
            
            # Also scale periodic fault magnitude if needed
            if abs(periodic_magnitude - 0.1) < 1e-6 and fault_type in DEFAULT_FAULT_MAGNITUDES:
                # Use 20% of the main fault magnitude as default for periodic
                scaled_periodic = DEFAULT_FAULT_MAGNITUDES[fault_type] * 0.2
                periodic_magnitude = scaled_periodic
                self.periodic_magnitude_var.set(periodic_magnitude)
                if periodic["enabled"]:
                    self.parent_app.add_log(f"Scaled periodic {fault_type} magnitude to {scaled_periodic}")
            
            periodic["magnitude"] = periodic_magnitude
            periodic["wheel"] = periodic_wheel
            
            # Update fault status
            self.update_fault_status()
            
            # Update constellation tab to show fault indicator
            try:
                self.parent_app.constellation_tab.update_satellite_listbox()
            except:
                pass
            
            # Log the configuration with actual values (display wheel as 1-based)
            if fault_enabled:
                self.parent_app.add_log(
                    f"Applied {fault_type} fault to {sat_name}: "
                    f"magnitude={fault_magnitude}, wheel=RW{fault_wheel+1}, time={fault_time}min"
                )
            else:
                self.parent_app.add_log(f"Disabled fault for {sat_name}")
                
            self.update_fault_summary()

    def update_fault_summary(self):
        """Update the fault summary table showing all satellites and their faults"""
        # Check if the summary tree exists yet
        if not hasattr(self, 'fault_summary_tree'):
            return
            
        # Clear existing items
        for item in self.fault_summary_tree.get_children():
            self.fault_summary_tree.delete(item)
        
        # Count statistics
        total_satellites = len(self.satellites)
        enabled_faults = 0
        
        # Add each satellite's fault configuration
        for sat in self.satellites:
            fault = sat["fault"]
            
            if fault["enabled"]:
                enabled_faults += 1
            
            # Determine status
            if fault["enabled"]:
                status = "ENABLED"
                status_color = 'red'
            else:
                status = "Disabled"
                status_color = 'gray'
            
            # Format magnitude based on fault type with proper units
            magnitude_str = str(fault["magnitude"])
            if fault["type"] == "friction":
                magnitude_str = f"{fault['magnitude']} N⋅m"
            elif fault["type"] == "power_limit":
                magnitude_str = f"{fault['magnitude']} W"
            elif fault["type"] == "battery":
                magnitude_str = f"{fault['magnitude']} W"
            elif fault["type"] == "encoder":
                magnitude_str = f"{fault['magnitude']}%"
            
            # Periodic status with RW 1-4 display
            periodic_str = "No"
            if fault["periodic"]["enabled"]:
                periodic_wheel_display = fault['periodic']['wheel'] + 1  # Convert 0-3 to 1-4
                periodic_str = f"Yes (RW{periodic_wheel_display}, {fault['periodic']['interval']}s)"
            
            # Display wheel number as 1-4 instead of 0-3
            wheel_display = f"RW{fault['wheel'] + 1}"
            
            # Insert row with all values properly formatted
            item = self.fault_summary_tree.insert('', 'end', values=(
                sat["name"],
                status,
                fault["type"].replace('_', ' ').title(),
                wheel_display,
                f"{fault['time']} min",
                magnitude_str,
                periodic_str
            ))
            
            # Color-code the status
            if fault["enabled"]:
                self.fault_summary_tree.item(item, tags=('enabled',))
            else:
                self.fault_summary_tree.item(item, tags=('disabled',))
        
        # Configure tags for coloring
        self.fault_summary_tree.tag_configure('enabled', foreground='darkred', font=('Segoe UI', 10, 'bold'))
        self.fault_summary_tree.tag_configure('disabled', foreground='gray')
        
        # Update statistics label
        if hasattr(self, 'fault_stats_label'):
            self.fault_stats_label.config(text=f"Total: {total_satellites} satellites | Faults: {enabled_faults} enabled")
            
        # Log summary update
        self.parent_app.add_log(f"Fault summary updated: {enabled_faults}/{total_satellites} faults enabled")