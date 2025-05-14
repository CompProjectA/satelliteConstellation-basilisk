#!/usr/bin/env python
"""
spacecraft_fault_gui.py

A modern Tkinter GUI for simulating and visualizing spacecraft faults,
including reaction wheel and battery faults, with 3D visualization capabilities.
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import sys
import os
import json
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
from datetime import datetime

# Fix path resolution to work with new project structure
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CORE_DIR = os.path.join(ROOT_DIR, 'core')
PLOTS_DIR = os.path.join(ROOT_DIR, 'plots')
ASSETS_DIR = os.path.join(ROOT_DIR, 'assets')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, CORE_DIR])

# Import CLI core functionality
from core.spacecraft_fault_cli import SimulationConfig, TargetDefinition, run_custom_simulation
# Import styling module
from style import setup_style


# Import help content from separate module (if available)
try:
    from core.rw_fault_help import (
        get_help_content, 
        get_general_help, 
        get_fault_description, 
        get_parameter_description
    )
except ImportError:
    # Fallback if help module is missing
    def get_help_content(topic): 
        return f"Help content for {topic} not available. Please ensure rw_fault_help.py is in the core directory."
    def get_general_help(section): 
        return f"General help for {section} not available. Please ensure rw_fault_help.py is in the core directory."
    def get_fault_description(fault_type): 
        descriptions = {
            "Friction": "Friction fault adds additional friction to the reaction wheel, simulating mechanical issues like bearing damage.",
            "Power Limit": "Power limit fault restricts the maximum power available to the reaction wheel, simulating power system failures.",
            "Encoder": "Encoder fault causes measurement errors in the reaction wheel speed feedback, leading to control issues.",
            "Battery": "Battery fault simulates increased power consumption or battery degradation, which can lead to power system failure."
        }
        return descriptions.get(fault_type, f"No description available for {fault_type}.")
    def get_parameter_description(fault_type): 
        descriptions = {
            "friction": "Friction fault adds a constant friction torque to the reaction wheel, requiring more motor torque to overcome.",
            "power_limit": "Power limit fault restricts how much electrical power is available to the reaction wheel motor.",
            "encoder": "Encoder fault creates errors in the speed measurements fed back to the control system.",
            "battery": "Battery fault increases power consumption, simulating damaged or degraded power system components."
        }
        return descriptions.get(fault_type, "")

class FaultSimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Spacecraft Fault Simulator")
        self.root.geometry("1000x750")
        self.root.minsize(900, 700)  # Set minimum window size
        
        # Try to load application icon
        try:
            if os.path.exists(os.path.join(ASSETS_DIR, "satellite_icon.ico")):
                self.root.iconbitmap(os.path.join(ASSETS_DIR, "satellite_icon.ico"))
        except Exception:
            pass  # Skip if icon loading fails
        
        self.config = SimulationConfig()
        
        # Apply styling to the application
        self.style = setup_style(self.root)
        
        # Create the GUI widgets
        self.create_widgets()

    def create_widgets(self):
        """Create the main GUI layout and widgets"""
        # Create main frame with padding
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Add title with icon
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Try to load icon if available
        try:
            logo_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "satellite_logo.png")).resize((32, 32)))
            logo_label = ttk.Label(title_frame, image=logo_img, background='#f5f5f5')
            logo_label.image = logo_img  # Keep a reference
            logo_label.pack(side=tk.LEFT, padx=(0, 10))
        except Exception:
            pass  # Skip icon if not available
                
        ttk.Label(title_frame, text="Spacecraft Fault Simulator", 
                style='Title.TLabel').pack(side=tk.LEFT)

        # IMPORTANT: Action buttons at the TOP of the window where they'll always be visible
        action_frame = ttk.Frame(main_frame, padding=(0, 5, 0, 15))
        action_frame.pack(fill=tk.X)
        
        # Run simulation button with distinct styling - made larger and more prominent
        run_btn = ttk.Button(action_frame, text="RUN SIMULATION", command=self.run_simulation, style='Run.TButton')
        run_btn.pack(side=tk.LEFT, padx=(0, 10), pady=5, ipadx=10, ipady=5)
        
        # Other buttons
        ttk.Button(action_frame, text="Export Config", command=self.export_config).pack(side=tk.LEFT, padx=5, pady=5)
        ttk.Button(action_frame, text="Help", command=self.show_general_help).pack(side=tk.LEFT, padx=5, pady=5)

        # Create notebook for tabs
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Create each tab
        self.sim_tab = ttk.Frame(self.notebook)
        self.fault_tab = ttk.Frame(self.notebook)
        self.viz_tab = ttk.Frame(self.notebook)

        self.notebook.add(self.sim_tab, text="Simulation Settings")
        self.notebook.add(self.fault_tab, text="Fault Configuration")
        self.notebook.add(self.viz_tab, text="Visualization")

        self.create_simulation_tab()
        self.create_fault_tab()
        self.create_viz_tab()

        # Status bar at the bottom
        status_frame = ttk.Frame(main_frame, style='Status.TFrame')
        status_frame.pack(fill=tk.X, pady=(10, 0))

        self.status_label = ttk.Label(status_frame, text="Ready", style='Status.TLabel')
        self.status_label.pack(side=tk.LEFT, padx=10, fill=tk.X, expand=True)
        
        # Add timestamp to status bar
        self.time_label = ttk.Label(status_frame, text=datetime.now().strftime("%Y-%m-%d %H:%M"), 
                                style='Status.TLabel')
        self.time_label.pack(side=tk.RIGHT, padx=10)
        
        # Update time every minute
        self.update_time()

    def update_time(self):
        """Update the time display in the status bar"""
        self.time_label.config(text=datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.root.after(60000, self.update_time)  # Update every minute

    def create_simulation_tab(self):
        """Create the Simulation Settings tab"""
        f = self.sim_tab
        
        # Add frame with styling
        sim_frame = ttk.LabelFrame(f, text="Simulation Parameters")
        sim_frame.pack(fill=tk.X, padx=5, pady=5)

        # Add info button
        self.add_help_button(sim_frame, "Simulation Parameters", 
                           "Simulation Parameters")

        # Simulation time
        time_frame = ttk.Frame(sim_frame)
        time_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(time_frame, text="Simulation Time (min):").pack(side=tk.LEFT, padx=(0, 10))
        self.sim_time = tk.DoubleVar(value=self.config.simulation_time)
        ttk.Entry(time_frame, textvariable=self.sim_time, width=10).pack(side=tk.LEFT)
        ttk.Label(time_frame, text="(Recommended: 30-60 minutes for full orbit)", 
                 style='Info.TLabel').pack(side=tk.LEFT, padx=10)

        # Output section
        out_frame = ttk.LabelFrame(f, text="Output Options")
        out_frame.pack(fill=tk.X, padx=5, pady=10)
        
        self.add_help_button(out_frame, "Output Options", 
                           "Output Options")

        # Checkboxes and filename
        option_frame = ttk.Frame(out_frame)
        option_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.show_plots = tk.BooleanVar(value=self.config.show_plots)
        self.save_binary = tk.BooleanVar(value=self.config.save_binary)
        self.save_plots = tk.BooleanVar(value=True)  # New option to save plots
        self.binary_filename = tk.StringVar(value=self.config.binary_filename)
        
        ttk.Checkbutton(option_frame, text="Show Plots", variable=self.show_plots).pack(side=tk.LEFT, padx=(0, 20))
        ttk.Checkbutton(option_frame, text="Save Plots", variable=self.save_plots).pack(side=tk.LEFT, padx=(0, 20))
        ttk.Checkbutton(option_frame, text="Save Binary for Vizard", variable=self.save_binary).pack(side=tk.LEFT, padx=(0, 20))
        
        # Filename section
        file_frame = ttk.Frame(out_frame)
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(file_frame, text="Binary Filename:").pack(side=tk.LEFT, padx=(0, 10))
        ttk.Entry(file_frame, textvariable=self.binary_filename, width=25).pack(side=tk.LEFT)
        ttk.Label(file_frame, text="(Without extension)", style='Info.TLabel').pack(side=tk.LEFT, padx=10)
        
        # Plot directory selection
        plot_frame = ttk.Frame(out_frame)
        plot_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(plot_frame, text="Plot Directory:").pack(side=tk.LEFT, padx=(0, 10))
        self.plot_dir = tk.StringVar(value=PLOTS_DIR)
        plot_entry = ttk.Entry(plot_frame, textvariable=self.plot_dir, width=40)
        plot_entry.pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(plot_frame, text="Browse...", command=self.browse_plot_dir).pack(side=tk.LEFT)
        
        # Add simulation information box
        info_frame = ttk.LabelFrame(f, text="About This Simulator")
        info_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=10)
        
        info_text = """
The Spacecraft Fault Simulator lets you model and analyze the effects of 
different fault types on spacecraft attitude control and power systems:

- Run realistic simulations of spacecraft dynamics in Earth orbit
- Inject various fault types into spacecraft hardware
- Visualize the effects on spacecraft attitude and performance
- Analyze the resulting data through comprehensive plots

Each fault type represents a different failure mode that can occur in real spacecraft:
- Friction faults simulate mechanical issues like bearing problems in reaction wheels
- Power limit faults simulate electrical power restrictions to the reaction wheels
- Encoder faults simulate sensor feedback errors in reaction wheels
- Battery faults simulate power system degradation or increased consumption
        
The simulator uses the Basilisk astrodynamics framework for high-fidelity modeling.
Results can be viewed in plots or visualized in 3D using Vizard visualization.
        """
        
        info_label = ttk.Label(info_frame, text=info_text, wraplength=950, justify="left")
        info_label.pack(padx=10, pady=10, fill=tk.BOTH)

    def browse_plot_dir(self):
        """Open directory browser to select plot output directory"""
        dir_path = filedialog.askdirectory(
            initialdir=self.plot_dir.get(),
            title="Select Directory for Plot Output"
        )
        if dir_path:
            self.plot_dir.set(dir_path)

    def create_fault_tab(self):
        """Create the Fault Configuration tab"""
        f = self.fault_tab
        
        # Fault selection dropdown
        fault_type_frame = ttk.LabelFrame(f, text="Fault Type")
        fault_type_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Add help button
        self.add_help_button(fault_type_frame, "Fault Types", 
                           "Fault Types")
        
        select_frame = ttk.Frame(fault_type_frame)
        select_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(select_frame, text="Select Fault Type:").pack(side=tk.LEFT, padx=(0, 10))
        
        # Fault types with descriptions
        self.fault_types = ["Friction", "Power Limit", "Encoder", "Battery"]
        self.fault_type = tk.StringVar(value=self.fault_types[0])
        fault_dropdown = ttk.Combobox(select_frame, textvariable=self.fault_type, values=self.fault_types, state="readonly", width=15)
        fault_dropdown.pack(side=tk.LEFT, padx=(0, 10))
        fault_dropdown.bind("<<ComboboxSelected>>", self.on_fault_type_change)
        
        # Add type description label that updates with selection
        self.fault_description = ttk.Label(fault_type_frame, 
                                       text=get_fault_description("Friction"),
                                       wraplength=950)
        self.fault_description.pack(padx=5, pady=5)
        
        # Frame for fault parameters (will be updated based on selected fault type)
        self.fault_param_frame = ttk.LabelFrame(f, text="Fault Parameters")
        self.fault_param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Add help button for parameters
        self.param_help_button = self.add_help_button(self.fault_param_frame, "Friction Fault", 
                                                "friction")
        
        # Initialize variables for fault parameters
        self.fault_mag = tk.DoubleVar(value=self.config.fault_magnitude)
        self.fault_time = tk.DoubleVar(value=self.config.fault_time)
        self.fault_wheel = tk.IntVar(value=self.config.fault_wheel_number)
        
        # Power limit specific variables
        self.power_limit = tk.DoubleVar(value=0.5)  # Default power limit value
        
        # Battery specific variables
        self.battery_drain = tk.DoubleVar(value=0.05)  # Default battery drain in kW (50W)
        
        # Periodic fault parameters
        self.periodic_frame = ttk.LabelFrame(f, text="Periodic Fault (Optional)")
        self.periodic_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Add help button for periodic faults
        self.add_help_button(self.periodic_frame, "Periodic Faults", "periodic")
        
        # Use pack for all periodic fault widgets to avoid grid/pack conflict
        periodic_control_frame = ttk.Frame(self.periodic_frame)
        periodic_control_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.enable_periodic = tk.BooleanVar(value=self.config.enable_periodic_fault)
        ttk.Checkbutton(periodic_control_frame, text="Enable Periodic Fault", 
                      variable=self.enable_periodic, command=self.toggle_periodic).pack(side=tk.LEFT)
        
        # Add tooltip/hint about periodic faults
        ttk.Label(periodic_control_frame, text="Adds cycling faults that repeat at regular intervals", 
                 style='Info.TLabel').pack(side=tk.LEFT, padx=20)
        
        # Parameters frame using pack layout
        periodic_params_frame = ttk.Frame(self.periodic_frame)
        periodic_params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.periodic_interval = tk.DoubleVar(value=self.config.periodic_fault_interval)
        self.periodic_mag = tk.DoubleVar(value=self.config.periodic_fault_magnitude)
        self.periodic_wheel = tk.IntVar(value=self.config.periodic_fault_wheel)
        
        # First row of parameters
        interval_frame = ttk.Frame(periodic_params_frame)
        interval_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(interval_frame, text="Interval (sec):").pack(side=tk.LEFT, padx=(0, 5))
        self.periodic_interval_entry = ttk.Entry(interval_frame, textvariable=self.periodic_interval, width=10, state="disabled")
        self.periodic_interval_entry.pack(side=tk.LEFT)
        
        # Second row of parameters
        mag_frame = ttk.Frame(periodic_params_frame)
        mag_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(mag_frame, text="Magnitude:").pack(side=tk.LEFT, padx=(0, 5))
        self.periodic_mag_entry = ttk.Entry(mag_frame, textvariable=self.periodic_mag, width=10, state="disabled")
        self.periodic_mag_entry.pack(side=tk.LEFT)
        
        # Third row of parameters
        wheel_frame = ttk.Frame(periodic_params_frame)
        wheel_frame.pack(side=tk.LEFT)
        ttk.Label(wheel_frame, text="Wheel (0-3):").pack(side=tk.LEFT, padx=(0, 5))
        self.periodic_wheel_spinbox = ttk.Spinbox(wheel_frame, from_=0, to=3, textvariable=self.periodic_wheel, width=4, state="disabled")
        self.periodic_wheel_spinbox.pack(side=tk.LEFT)
        
        # Show friction fault parameters by default
        self.show_friction_params()
        
        # Add fault response visualization frame (shows what to expect)
        viz_frame = ttk.LabelFrame(f, text="Expected Fault Response")
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Try to load and display sample fault response graph
        try:
            response_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "friction_response.png")).resize((950, 250)))
            self.response_label = ttk.Label(viz_frame, image=response_img)
            self.response_label.image = response_img  # Keep a reference
            self.response_label.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        except:
            # Fallback to text description if image not available
            self.response_label = ttk.Label(viz_frame, 
                                         text="Friction faults typically cause increased power consumption and reduced wheel speed response.",
                                         wraplength=950)
            self.response_label.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

    def on_fault_type_change(self, event):
        """Handle fault type selection change"""
        # Clear previous content
        for widget in self.fault_param_frame.winfo_children():
            if not isinstance(widget, ttk.Button):  # Keep the help button
                widget.destroy()
            
        # Show parameters based on fault type
        fault_type = self.fault_type.get()
        
        # Update help button content and fault description
        if fault_type == "Friction":
            self.param_help_button.configure(command=lambda: self.show_help_content("Friction Fault", "friction"))
            self.fault_description.configure(text=get_fault_description("Friction"))
            self.show_friction_params()
            
            # Update response visualization
            try:
                response_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "friction_response.png")).resize((950, 250)))
                self.response_label.configure(image=response_img)
                self.response_label.image = response_img
            except:
                self.response_label.configure(text="Friction faults typically cause increased power consumption and reduced wheel speed response.")
                
        elif fault_type == "Power Limit":
            self.param_help_button.configure(command=lambda: self.show_help_content("Power Limit Fault", "power_limit"))
            self.fault_description.configure(text=get_fault_description("Power Limit"))
            self.show_power_limit_params()
            
            # Update response visualization
            try:
                response_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "power_response.png")).resize((950, 250)))
                self.response_label.configure(image=response_img)
                self.response_label.image = response_img
            except:
                self.response_label.configure(text="Power limit faults cause torque saturation when the wheel requires more power than the limit allows.")
                
        elif fault_type == "Encoder":
            self.param_help_button.configure(command=lambda: self.show_help_content("Encoder Fault", "encoder"))
            self.fault_description.configure(text=get_fault_description("Encoder"))
            self.show_encoder_params()
            
            # Update response visualization
            try:
                response_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "encoder_response.png")).resize((950, 250)))
                self.response_label.configure(image=response_img)
                self.response_label.image = response_img
            except:
                self.response_label.configure(text="Encoder faults cause discrepancies between commanded and actual wheel speeds, leading to control oscillations.")

        elif fault_type == "Battery":
            self.param_help_button.configure(command=lambda: self.show_help_content("Battery Fault", "battery"))
            self.fault_description.configure(text=get_fault_description("Battery"))
            self.show_battery_params()
            
            # Update response visualization
            try:
                response_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "battery_response.png")).resize((950, 250)))
                self.response_label.configure(image=response_img)
                self.response_label.image = response_img
            except:
                self.response_label.configure(text="Battery faults cause increased power consumption, leading to faster battery depletion and potential power system failure.")

    def show_friction_params(self):
        """Display friction fault parameters"""
        # Use pack layout consistently for all parameter frames
        param_frame = ttk.Frame(self.fault_param_frame)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Row 1: Magnitude, Wheel, and Time
        row1 = ttk.Frame(param_frame)
        row1.pack(fill=tk.X, pady=5)
        
        # Magnitude
        mag_frame = ttk.Frame(row1)
        mag_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(mag_frame, text="Magnitude:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(mag_frame, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT)
        ttk.Label(mag_frame, text="(0.0001 - 0.001 typical)").pack(side=tk.LEFT, padx=5)
        
        # Wheel
        wheel_frame = ttk.Frame(row1)
        wheel_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(wheel_frame, text="Wheel (0-3):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Spinbox(wheel_frame, from_=0, to=3, textvariable=self.fault_wheel, width=4).pack(side=tk.LEFT)
        
        # Time
        time_frame = ttk.Frame(row1)
        time_frame.pack(side=tk.LEFT)
        ttk.Label(time_frame, text="Time (min):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT)
        
        # Add description with details
        ttk.Label(self.fault_param_frame, text=get_parameter_description("friction"), 
                 wraplength=950, style='Info.TLabel').pack(padx=5, pady=5)

    def show_power_limit_params(self):
        """Display power limit fault parameters"""
        # Use pack layout consistently for all parameter frames
        param_frame = ttk.Frame(self.fault_param_frame)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Row 1: Power Limit, Wheel, and Time
        row1 = ttk.Frame(param_frame)
        row1.pack(fill=tk.X, pady=5)
        
        # Power Limit
        limit_frame = ttk.Frame(row1)
        limit_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(limit_frame, text="Power Limit (W):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(limit_frame, textvariable=self.power_limit, width=10).pack(side=tk.LEFT)
        ttk.Label(limit_frame, text="(0.2 - 1.0 typical)").pack(side=tk.LEFT, padx=5)
        
        # Wheel
        wheel_frame = ttk.Frame(row1)
        wheel_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(wheel_frame, text="Wheel (0-3):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Spinbox(wheel_frame, from_=0, to=3, textvariable=self.fault_wheel, width=4).pack(side=tk.LEFT)
        
        # Time
        time_frame = ttk.Frame(row1)
        time_frame.pack(side=tk.LEFT)
        ttk.Label(time_frame, text="Time (min):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT)
        
        # Add description with details
        ttk.Label(self.fault_param_frame, text=get_parameter_description("power_limit"), 
                 wraplength=950, style='Info.TLabel').pack(padx=5, pady=5)

    def show_encoder_params(self):
        """Display encoder fault parameters"""
        # Use pack layout consistently for all parameter frames
        param_frame = ttk.Frame(self.fault_param_frame)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Row 1: Wheel and Time
        row1 = ttk.Frame(param_frame)
        row1.pack(fill=tk.X, pady=5)
        
        # Wheel
        wheel_frame = ttk.Frame(row1)
        wheel_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(wheel_frame, text="Wheel (0-3):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Spinbox(wheel_frame, from_=0, to=3, textvariable=self.fault_wheel, width=4).pack(side=tk.LEFT)
        
        # Time
        time_frame = ttk.Frame(row1)
        time_frame.pack(side=tk.LEFT)
        ttk.Label(time_frame, text="Time (min):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT)
        
        # Add description with details
        ttk.Label(self.fault_param_frame, text=get_parameter_description("encoder"), 
                 wraplength=950, style='Info.TLabel').pack(padx=5, pady=5)

    def show_battery_params(self):
        """Display battery fault parameters"""
        # Use pack layout consistently for all parameter frames
        param_frame = ttk.Frame(self.fault_param_frame)
        param_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Row 1: Power Drain and Time
        row1 = ttk.Frame(param_frame)
        row1.pack(fill=tk.X, pady=5)
        
        # Power Drain
        drain_frame = ttk.Frame(row1)
        drain_frame.pack(side=tk.LEFT, padx=(0, 20))
        ttk.Label(drain_frame, text="Power Drain (kW):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(drain_frame, textvariable=self.battery_drain, width=10).pack(side=tk.LEFT)
        ttk.Label(drain_frame, text="(0.01 - 0.1 typical)").pack(side=tk.LEFT, padx=5)
        
        # Time
        time_frame = ttk.Frame(row1)
        time_frame.pack(side=tk.LEFT)
        ttk.Label(time_frame, text="Time (min):").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(time_frame, textvariable=self.fault_time, width=10).pack(side=tk.LEFT)
        
        # Add description with details
        ttk.Label(self.fault_param_frame, text=get_parameter_description("battery"), 
                 wraplength=950, style='Info.TLabel').pack(padx=5, pady=5)

    def toggle_periodic(self):
        """Enable/disable periodic fault entry fields"""
        state = "normal" if self.enable_periodic.get() else "disabled"
        self.periodic_interval_entry.config(state=state)
        self.periodic_mag_entry.config(state=state)
        self.periodic_wheel_spinbox.config(state=state)

    def create_viz_tab(self):
        """Create the Visualization tab"""
        f = self.viz_tab

        # Target Management Section
        target_frame = ttk.LabelFrame(f, text="Target Management")
        target_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Add help button
        self.add_help_button(target_frame, "Target Management", "targets")

        # Add controls to add a new target
        control_frame = ttk.Frame(target_frame)
        control_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.new_target_name = tk.StringVar()
        self.new_target_lat = tk.DoubleVar()
        self.new_target_lon = tk.DoubleVar()
        self.new_target_color = tk.StringVar(value="red")
        
        # Row 1: Name, Lat, Lon inputs
        input_frame = ttk.Frame(control_frame)
        input_frame.pack(fill=tk.X, pady=5)
        
        # Name
        name_frame = ttk.Frame(input_frame)
        name_frame.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Label(name_frame, text="Name:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(name_frame, textvariable=self.new_target_name, width=15).pack(side=tk.LEFT)
        
        # Latitude
        lat_frame = ttk.Frame(input_frame)
        lat_frame.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Label(lat_frame, text="Latitude:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(lat_frame, textvariable=self.new_target_lat, width=10).pack(side=tk.LEFT)
        
        # Longitude
        lon_frame = ttk.Frame(input_frame)
        lon_frame.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Label(lon_frame, text="Longitude:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(lon_frame, textvariable=self.new_target_lon, width=10).pack(side=tk.LEFT)
        
        # Color
        color_frame = ttk.Frame(input_frame)
        color_frame.pack(side=tk.LEFT)
        ttk.Label(color_frame, text="Color:").pack(side=tk.LEFT, padx=(0, 5))
        color_combo = ttk.Combobox(color_frame, textvariable=self.new_target_color, 
                               values=["red", "blue", "green", "yellow", "purple", "orange", "cyan", "magenta"], 
                               width=8, state="readonly")
        color_combo.pack(side=tk.LEFT)
        
        # Button row
        button_frame = ttk.Frame(target_frame)
        button_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(button_frame, text="Add Target", command=self.add_target).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Remove Selected", command=self.remove_target).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear All", command=self.clear_targets).pack(side=tk.LEFT, padx=5)
        
        # Add info text
        ttk.Label(target_frame, text="Targets are ground locations visible to the spacecraft during orbit. The simulation will show when each target is visible.",
                 style='Info.TLabel').pack(fill=tk.X, padx=5, pady=5)

        # Target List
        list_frame = ttk.Frame(target_frame)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Add a scrollbar for the treeview
        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.target_list = ttk.Treeview(list_frame, columns=('name', 'lat', 'lon', 'color'), show='headings', height=8)
        self.target_list.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        
        # Connect scrollbar to the treeview
        scrollbar.config(command=self.target_list.yview)
        self.target_list.config(yscrollcommand=scrollbar.set)
        
        # Configure column headings and widths
        for col, width in (('name', 150), ('lat', 100), ('lon', 100), ('color', 100)):
            self.target_list.heading(col, text=col.capitalize())
            self.target_list.column(col, width=width, anchor=tk.CENTER)
        
        # Camera settings
        camera_frame = ttk.LabelFrame(f, text="Camera Settings")
        camera_frame.pack(fill=tk.X, padx=5, pady=10)
        
        # Add help button
        self.add_help_button(camera_frame, "Camera Settings", "camera")
        
        camera_controls = ttk.Frame(camera_frame)
        camera_controls.pack(fill=tk.X, padx=5, pady=5)
        
        self.camera_x = tk.DoubleVar(value=self.config.camera_position[0])
        self.camera_y = tk.DoubleVar(value=self.config.camera_position[1])
        self.camera_z = tk.DoubleVar(value=self.config.camera_position[2])
        
        # Use pack layout for consistency
        pos_label = ttk.Label(camera_controls, text="Camera Position (Body Frame):")
        pos_label.pack(side=tk.LEFT, padx=(0, 10))
        
        # X coordinate
        x_frame = ttk.Frame(camera_controls)
        x_frame.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Label(x_frame, text="X:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(x_frame, textvariable=self.camera_x, width=8).pack(side=tk.LEFT)
        
        # Y coordinate
        y_frame = ttk.Frame(camera_controls)
        y_frame.pack(side=tk.LEFT, padx=(0, 10))
        ttk.Label(y_frame, text="Y:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(y_frame, textvariable=self.camera_y, width=8).pack(side=tk.LEFT)
        
        # Z coordinate
        z_frame = ttk.Frame(camera_controls)
        z_frame.pack(side=tk.LEFT)
        ttk.Label(z_frame, text="Z:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Entry(z_frame, textvariable=self.camera_z, width=8).pack(side=tk.LEFT)
        
        # Add preset buttons
        presets_frame = ttk.Frame(camera_frame)
        presets_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(presets_frame, text="Presets:").pack(side=tk.LEFT, padx=(0,5))
        ttk.Button(presets_frame, text="Side View", 
                  command=lambda: self.set_camera_preset(0.0, 2.0, 0.0)).pack(side=tk.LEFT, padx=5)
        ttk.Button(presets_frame, text="Top View", 
                  command=lambda: self.set_camera_preset(0.0, 0.0, 2.0)).pack(side=tk.LEFT, padx=5)
        ttk.Button(presets_frame, text="Front View", 
                  command=lambda: self.set_camera_preset(2.0, 0.0, 0.0)).pack(side=tk.LEFT, padx=5)
        
        # Add camera info text
        ttk.Label(camera_frame, 
                 text="The camera position determines the viewpoint in the 3D visualization. Adjust X, Y, and Z coordinates to change perspective.",
                 style='Info.TLabel').pack(fill=tk.X, padx=5, pady=5)
        
        # Visualization preview frame
        preview_frame = ttk.LabelFrame(f, text="Visualization Preview")
        preview_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Try to load and display visualization preview
        try:
            preview_img = ImageTk.PhotoImage(Image.open(os.path.join(ASSETS_DIR, "vizard_preview.png")).resize((950, 300)))
            preview_label = ttk.Label(preview_frame, image=preview_img)
            preview_label.image = preview_img  # Keep a reference
            preview_label.pack(fill=tk.BOTH, expand=True)
        except:
            # Fallback to text description if image not available
            ttk.Label(preview_frame, 
                     text="""
The Vizard visualization shows a detailed 3D model of the spacecraft with reaction wheels. 
It allows you to observe the spacecraft's attitude changes in response to faults.
After running the simulation, load the binary file in Vizard to view the animation.
                     """,
                     wraplength=950).pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        # Load initial targets
        for target in self.config.targets:
            self.target_list.insert('', tk.END, values=(target.name, target.latitude, target.longitude, target.color))

    def set_camera_preset(self, x, y, z):
        """Set camera position to a preset view"""
        self.camera_x.set(x)
        self.camera_y.set(y)
        self.camera_z.set(z)

    def add_target(self):
        """Add a new target to the list"""
        try:
            name = self.new_target_name.get()
            lat = self.new_target_lat.get()
            lon = self.new_target_lon.get()
            color = self.new_target_color.get()
            
            if not name:
                messagebox.showerror("Error", "Target name cannot be empty")
                return
                
            if lat < -90 or lat > 90:
                messagebox.showerror("Error", "Latitude must be between -90 and 90 degrees")
                return
                
            if lon < -180 or lon > 180:
                messagebox.showerror("Error", "Longitude must be between -180 and 180 degrees")
                return
                
            self.target_list.insert('', tk.END, values=(name, lat, lon, color))
            
            # Clear entries
            self.new_target_name.set("")
            self.new_target_lat.set(0.0)
            self.new_target_lon.set(0.0)
            
        except Exception as e:
            messagebox.showerror("Error Adding Target", str(e))

    def remove_target(self):
        """Remove selected targets from the list"""
        selection = self.target_list.selection()
        if selection:
            for item in selection:
                self.target_list.delete(item)
        else:
            messagebox.showinfo("Info", "Please select one or more targets to remove")
            
    def clear_targets(self):
        """Remove all targets from the list"""
        if messagebox.askyesno("Confirm", "Are you sure you want to remove all targets?"):
            for item in self.target_list.get_children():
                self.target_list.delete(item)

    def add_help_button(self, parent, title, topic=None):
        """Add a help button to a frame with a popup dialog for the given topic"""
        help_button = ttk.Button(parent, text="?", width=2, style='Help.TButton',
                               command=lambda: self.show_help_content(title, topic))
        help_button.pack(side=tk.RIGHT, padx=5, pady=5, anchor='ne')
        return help_button

    def show_help_content(self, title, topic=None):
        """Show help content in a popup dialog"""
        help_window = tk.Toplevel(self.root)
        help_window.title(f"Help: {title}")
        help_window.geometry("600x500")
        help_window.transient(self.root)  # Make window modal
        help_window.grab_set()
        
        # Create a frame for the help content
        frame = ttk.Frame(help_window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Add title
        ttk.Label(frame, text=title, style='Title.TLabel').pack(fill=tk.X, pady=(0, 10))
        
        # Create scrolled text widget for content
        help_text = scrolledtext.ScrolledText(frame, wrap=tk.WORD, width=70, height=20)
        help_text.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Insert the content
        content = get_help_content(topic) if topic else "No help content available for this topic."
        help_text.insert(tk.END, content)
        
        # Make it read-only
        help_text.config(state=tk.DISABLED)
        
        # Add close button
        ttk.Button(frame, text="Close", command=help_window.destroy).pack(pady=10)
        
    def show_general_help(self):
        """Show general help about the simulator"""
        help_window = tk.Toplevel(self.root)
        help_window.title("Spacecraft Fault Simulator Help")
        help_window.geometry("800x600")
        help_window.transient(self.root)
        help_window.grab_set()
        
        # Create notebook for help topics
        notebook = ttk.Notebook(help_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create tabs for different help sections
        overview_tab = ttk.Frame(notebook)
        simulation_tab = ttk.Frame(notebook)
        faults_tab = ttk.Frame(notebook)
        viz_tab = ttk.Frame(notebook)
        
        notebook.add(overview_tab, text="Overview")
        notebook.add(simulation_tab, text="Simulation")
        notebook.add(faults_tab, text="Fault Types")
        notebook.add(viz_tab, text="Visualization")
        
        # Overview tab content
        overview_frame = ttk.Frame(overview_tab, padding=10)
        overview_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(overview_frame, text="Spacecraft Fault Simulator", style='Title.TLabel').pack(pady=(0, 10))
        
        overview_text = scrolledtext.ScrolledText(overview_frame, wrap=tk.WORD)
        overview_text.pack(fill=tk.BOTH, expand=True)
        overview_text.insert(tk.END, get_general_help("overview"))
        overview_text.config(state=tk.DISABLED)
        
        # Simulation tab content
        sim_frame = ttk.Frame(simulation_tab, padding=10)
        sim_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(sim_frame, text="Simulation Settings", style='Title.TLabel').pack(pady=(0, 10))
        
        sim_text = scrolledtext.ScrolledText(sim_frame, wrap=tk.WORD)
        sim_text.pack(fill=tk.BOTH, expand=True)
        sim_text.insert(tk.END, get_general_help("simulation"))
        sim_text.config(state=tk.DISABLED)
        
        # Faults tab content
        faults_frame = ttk.Frame(faults_tab, padding=10)
        faults_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(faults_frame, text="Fault Types", style='Title.TLabel').pack(pady=(0, 10))
        
        faults_text = scrolledtext.ScrolledText(faults_frame, wrap=tk.WORD)
        faults_text.pack(fill=tk.BOTH, expand=True)
        faults_text.insert(tk.END, get_general_help("faults"))
        faults_text.config(state=tk.DISABLED)
        
        # Visualization tab content
        viz_frame = ttk.Frame(viz_tab, padding=10)
        viz_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(viz_frame, text="Visualization", style='Title.TLabel').pack(pady=(0, 10))
        
        viz_text = scrolledtext.ScrolledText(viz_frame, wrap=tk.WORD)
        viz_text.pack(fill=tk.BOTH, expand=True)
        viz_text.insert(tk.END, get_general_help("visualization"))
        viz_text.config(state=tk.DISABLED)
        
        # Add close button
        ttk.Button(help_window, text="Close", command=help_window.destroy).pack(pady=10)

 
    def run_simulation(self):
        """Run the simulation with the current settings"""
        try:
            # Validate fields
            if self.sim_time.get() <= 0:
                raise ValueError("Simulation time must be positive")
                
            if self.fault_time.get() >= self.sim_time.get():
                raise ValueError("Fault time must be less than simulation time")
                
            if self.save_binary.get() and not self.binary_filename.get():
                raise ValueError("Binary filename cannot be empty when Save Binary is enabled")
            
            # Update base config
            self.config.simulation_time = self.sim_time.get()
            self.config.show_plots = self.show_plots.get()
            self.config.save_binary = self.save_binary.get()
            self.config.binary_filename = self.binary_filename.get()
            
            # Set plot directory for saving plots
            plots_dir = self.plot_dir.get()
            if self.save_plots.get() and not os.path.exists(plots_dir):
                try:
                    os.makedirs(plots_dir, exist_ok=True)
                    print(f"Created plots directory: {plots_dir}")
                except Exception as e:
                    messagebox.showwarning("Warning", f"Could not create plots directory: {e}")
            
            # Update fault parameters based on selected fault type
            fault_type = self.fault_type.get()
            self.config.fault_type = fault_type.lower().replace(" ", "_")
            
            if fault_type == "Friction":
                if self.fault_mag.get() <= 0:
                    raise ValueError("Friction magnitude must be positive")
                self.config.fault_magnitude = self.fault_mag.get()
            elif fault_type == "Power Limit":
                if self.power_limit.get() <= 0:
                    raise ValueError("Power limit must be positive")
                self.config.fault_magnitude = self.power_limit.get()
            elif fault_type == "Battery":
                if self.battery_drain.get() <= 0:
                    raise ValueError("Battery power drain must be positive")
                self.config.fault_magnitude = self.battery_drain.get()
            # For encoder fault, we don't need a magnitude parameter
            
            self.config.fault_time = self.fault_time.get()
            self.config.fault_wheel_number = self.fault_wheel.get()
            
            # Update periodic fault settings
            self.config.enable_periodic_fault = self.enable_periodic.get()
            if self.enable_periodic.get():
                if self.periodic_interval.get() <= 0:
                    raise ValueError("Periodic fault interval must be positive")
                if self.periodic_mag.get() <= 0:
                    raise ValueError("Periodic fault magnitude must be positive")
                    
                self.config.periodic_fault_interval = self.periodic_interval.get()
                self.config.periodic_fault_magnitude = self.periodic_mag.get()
                self.config.periodic_fault_wheel = self.periodic_wheel.get()
            
            # Update camera position
            self.config.camera_position = [self.camera_x.get(), self.camera_y.get(), self.camera_z.get()]
            
            # Update targets
            self.config.targets = []
            for item in self.target_list.get_children():
                vals = self.target_list.item(item)['values']
                self.config.targets.append(TargetDefinition(vals[0], float(vals[1]), float(vals[2]), vals[3]))

            # Update status
            self.status_label.config(text="Running simulation...")
            self.root.update()

            # Run the simulation
            scenario, viz, figureList, output_dir = run_custom_simulation(self.config)

            # If save plots is enabled, save them to the specified directory
            saved_plots = []
            if self.save_plots.get() and figureList:
                # Make sure the plot directory exists
                if not os.path.exists(plots_dir):
                    os.makedirs(plots_dir, exist_ok=True)
                    
                for name, fig in figureList.items():
                    try:
                        # Create a more descriptive filename with timestamp and fault type
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        plot_filename = f"{name}_{self.config.fault_type}_{timestamp}.png"
                        plot_path = os.path.join(plots_dir, plot_filename)
                        
                        # Save the figure with high resolution
                        fig.savefig(plot_path, bbox_inches='tight', dpi=300)
                        saved_plots.append((name, plot_path))
                        print(f"Saved plot: {plot_path}")
                    except Exception as e:
                        print(f"Warning: Failed to save plot {name}: {e}")

            # Update status upon completion
            self.status_label.config(text="Simulation completed successfully")
            
            # Look for the binary file in multiple potential locations
            binary_found = False
            binary_path = ""
            if self.config.save_binary:
                binary_paths = [
                    os.path.join(VIZ_DIR, f"{self.config.binary_filename}_UnityViz.bin"),
                    os.path.join(VIZ_DIR, "_VizFiles", f"{self.config.binary_filename}_UnityViz.bin")
                ]
                
                for path in binary_paths:
                    if os.path.exists(path):
                        binary_found = True
                        binary_path = path
                        break
            
            # Show results summary
            summary = f"""
    Simulation completed successfully!

    Results saved to: {output_dir}

    Fault type: {fault_type}
    Fault injected at: {self.config.fault_time} minutes
    """
            if fault_type in ["Friction", "Power Limit", "Encoder"]:
                summary += f"Affected wheel: {self.config.fault_wheel_number}\n"

            if self.save_plots.get():
                if saved_plots:
                    summary += f"\nPlots saved to: {plots_dir}\n"
                    for i, (name, plot) in enumerate(saved_plots):
                        if i < 4:  # Only list the first 4 plots to avoid overly long message
                            summary += f" - {name}: {os.path.basename(plot)}\n"
                    if len(saved_plots) > 4:
                        summary += f" - ...and {len(saved_plots) - 4} more plots\n"
                else:
                    summary += "\nNo plots were generated or saved.\n"
                
            if self.config.save_binary:
                if binary_found:
                    summary += f"\nVisualization binary saved as: {binary_path}"
                    summary += "\nYou can open this file in Vizard for 3D visualization."
                else:
                    summary += "\nWarning: Visualization binary file was not created."
                    summary += "\nCheck if Vizard visualization is properly configured."
            
            # Show the message box
            messagebox.showinfo("Simulation Complete", summary)
            
            # If plots were created and should be shown, open the fault-specific ones
            if self.show_plots.get() and saved_plots:
                try:
                    import subprocess
                    import platform
                    
                    # Define fault-specific plot names to look for
                    fault_specific_plots = {
                        "Friction": ["RWFriction", "attitudeErrorNorm"],
                        "Power Limit": ["RWPower", "attitudeErrorNorm"],
                        "Encoder": ["RWEncoder", "attitudeErrorNorm"],
                        "Battery": ["BatteryStatus", "attitudeErrorNorm"]
                    }
                    
                    # Get the relevant plots for this fault type
                    plots_to_open = []
                    relevant_plot_types = fault_specific_plots.get(fault_type, ["attitudeErrorNorm"])
                    
                    # Find the matching plots
                    for plot_type in relevant_plot_types:
                        for name, path in saved_plots:
                            if plot_type in name:
                                plots_to_open.append(path)
                                break
                    
                    # If no specific plots found, use the first one
                    if not plots_to_open and saved_plots:
                        plots_to_open = [saved_plots[0][1]]
                    
                    # Open the selected plots with the default image viewer
                    for plot_path in plots_to_open:
                        if platform.system() == 'Windows':
                            os.startfile(plot_path)
                        elif platform.system() == 'Darwin':  # macOS
                            subprocess.call(['open', plot_path])
                        else:  # Linux
                            subprocess.call(['xdg-open', plot_path])
                        
                        print(f"Opened plot: {plot_path}")
                        
                except Exception as e:
                    print(f"Warning: Could not open plot: {e}")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text=f"Error: {str(e)}")



    def export_config(self):
        """Export the current configuration to a JSON file"""
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".json",
                                                     filetypes=[("JSON files", "*.json")],
                                                     title="Save Configuration")
            if file_path:
                # Update config object with current GUI values
                self.config.simulation_time = self.sim_time.get()
                self.config.show_plots = self.show_plots.get()
                self.config.save_binary = self.save_binary.get()
                self.config.binary_filename = self.binary_filename.get()
                self.config.fault_type = self.fault_type.get().lower().replace(" ", "_")
                
                # Get fault magnitude based on fault type
                if self.fault_type.get() == "Friction":
                    self.config.fault_magnitude = self.fault_mag.get()
                elif self.fault_type.get() == "Power Limit":
                    self.config.fault_magnitude = self.power_limit.get()
                elif self.fault_type.get() == "Battery":
                    self.config.fault_magnitude = self.battery_drain.get()
                
                self.config.fault_time = self.fault_time.get()
                self.config.fault_wheel_number = self.fault_wheel.get()
                self.config.enable_periodic_fault = self.enable_periodic.get()
                self.config.periodic_fault_interval = self.periodic_interval.get()
                self.config.periodic_fault_magnitude = self.periodic_mag.get()
                self.config.periodic_fault_wheel = self.periodic_wheel.get()
                self.config.camera_position = [self.camera_x.get(), self.camera_y.get(), self.camera_z.get()]
                
                # Create a dictionary representation of the configuration
                config_dict = {
                    "simulation_time": self.config.simulation_time,
                    "show_plots": self.config.show_plots,
                    "save_binary": self.config.save_binary,
                    "binary_filename": self.config.binary_filename,
                    "fault_type": self.config.fault_type,
                    "fault_magnitude": self.config.fault_magnitude,
                    "fault_time": self.config.fault_time,
                    "fault_wheel_number": self.config.fault_wheel_number,
                    "enable_periodic_fault": self.config.enable_periodic_fault,
                    "periodic_fault_interval": self.config.periodic_fault_interval,
                    "periodic_fault_magnitude": self.config.periodic_fault_magnitude,
                    "periodic_fault_wheel": self.config.periodic_fault_wheel,
                    "camera_position": self.config.camera_position,
                    "targets": [
                        {
                            "name": self.target_list.item(i)['values'][0],
                            "latitude": self.target_list.item(i)['values'][1],
                            "longitude": self.target_list.item(i)['values'][2],
                            "color": self.target_list.item(i)['values'][3]
                        } for i in self.target_list.get_children()
                    ]
                }
                
                # Save to file with nice formatting
                with open(file_path, 'w') as f:
                    json.dump(config_dict, f, indent=4)
                    
                self.status_label.config(text=f"Configuration saved to {os.path.basename(file_path)}")
                messagebox.showinfo("Export Success", f"Configuration saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))
            self.status_label.config(text="Export failed")

if __name__ == "__main__":
    root = tk.Tk()
    app = FaultSimulationGUI(root)
    root.mainloop()