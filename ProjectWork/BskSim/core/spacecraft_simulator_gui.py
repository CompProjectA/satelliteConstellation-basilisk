#!/usr/bin/env python
"""
spacecraft_simulator_gui.py

A comprehensive GUI for the Basilisk spacecraft fault simulator that supports
satellite constellations, fault injection, target assignment, and camera configuration.
FIXED: Better default orbital parameters, enhanced target visibility, improved layout
"""
import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import numpy as np
import threading
import time
import subprocess
import inspect
import json
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
from datetime import datetime
import logging

# Add near the top of the file
from plots import generate_fault_plots, generate_constellation_plots

# Import style module
from style import setup_style

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CORE_DIR = os.path.join(ROOT_DIR, 'core')
PLOTS_DIR = os.path.join(ROOT_DIR, 'plots')
ASSETS_DIR = os.path.join(ROOT_DIR, 'assets')

sys.path.extend([ROOT_DIR, CORE_DIR])

# Import simulation functionality to run simulations directly
from spacecraft_simulation import SimulationConfig, TargetDefinition, run_custom_simulation

# Import GUI tab modules
from gui_tab import (
    ConstellationTab,
    FaultTab,
    TargetTab,
    VisualizationTab,
    OutputTab
)

class SatelliteSimulatorApp:
    """Main application class for the Spacecraft Fault Simulator GUI"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Spacecraft Constellation Fault Simulator - ENHANCED")
        self.root.geometry("1200x900")  # Larger default window
        self.root.minsize(1000, 800)   # Larger minimum size
        
        # Setup logging
        logging.basicConfig(level=logging.INFO, 
                         format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger('simulator')
        
        # Fix path resolution to match project structure
        filename = inspect.getframeinfo(inspect.currentframe()).filename
        path = os.path.dirname(os.path.abspath(filename))
        self.ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
        self.FAULTS_DIR = os.path.join(self.ROOT_DIR, 'faults')
        self.MODELS_DIR = os.path.join(self.ROOT_DIR, 'models')
        self.PLOTTING_DIR = os.path.join(self.ROOT_DIR, 'plots')
        self.LOGS_DIR = os.path.join(self.ROOT_DIR, 'logs')
        self.VIZ_DIR = os.path.join(self.ROOT_DIR, 'Vizfile')
        
        # Create directories if they don't exist
        for dir_path in [self.FAULTS_DIR, self.MODELS_DIR, self.PLOTTING_DIR, self.LOGS_DIR, self.VIZ_DIR]:
            os.makedirs(dir_path, exist_ok=True)
        
        # Create _VizFiles subdirectory
        viz_files_dir = os.path.join(self.VIZ_DIR, "_VizFiles")
        if not os.path.exists(viz_files_dir):
            os.makedirs(viz_files_dir, exist_ok=True)
            
        # Set up the style
        self.style = setup_style(root)
        
        # Initialize data structures
        self.init_data()
        
        # Create main UI components
        self.create_main_ui()
        
        # Update status
        self.update_status("Ready - ENHANCED with better orbital parameters and target visibility")
    
    def init_data(self):
        """Initialize application data structures with FIXED parameters"""
        # Default paths
        self.output_dir = self.LOGS_DIR
        self.plots_dir = self.PLOTTING_DIR
        self.bin_dir = self.VIZ_DIR
        self.binary_filename = "saib_fault_viz"
        
        # FIXED: Simulation parameters - proper 30 minutes default
        self.simulation_time = tk.DoubleVar(value=30.0)  # 30 minutes default
        self.show_plots = tk.BooleanVar(value=True)
        self.save_plots = tk.BooleanVar(value=True)
        self.save_binary = tk.BooleanVar(value=True)
        
        # FIXED: Satellites data with ENHANCED orbital parameters for better visibility
        self.satellites = []
        self.add_satellite("Sat1", 0.0)    # FIXED: Simplified names
        self.add_satellite("Sat2", 90.0)
        self.add_satellite("Sat3", 180.0)
        self.add_satellite("Sat4", 270.0)
        self.current_satellite_index = tk.IntVar(value=0)
        
        # FIXED: Default targets with enhanced visibility and auto-assignment
        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "#FF0000", "priority": 3, "assigned_to": ["Sat1"]},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "#0000FF", "priority": 2, "assigned_to": ["Sat2"]},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "#00FF00", "priority": 2, "assigned_to": ["Sat3"]},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "#FFFF00", "priority": 1, "assigned_to": ["Sat4"]}
        ]
        
        # Running state
        self.is_running = False
        
        # Generated plots from last simulation
        self.latest_plots = []
        
    def add_satellite(self, name, true_anomaly=0.0):
        """FIXED: Add a new satellite with ENHANCED orbital parameters for better visibility"""
        satellite = {
            "name": name,
            "orbit": {
                "a": 8000.0,   # FIXED: Higher altitude (1629km) for better visibility from cameras
                "e": 0.05,     # Small eccentricity for stable elliptical orbit
                "i": 55.0,     # Inclination in degrees (good for ground coverage)
                "Omega": 45.0, # Right ascension of ascending node in degrees
                "omega": 30.0, # Argument of periapsis in degrees
                "f": true_anomaly  # True anomaly in degrees - user specified
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
                "position": [0.0, 0.0, 5.0],  # FIXED: Higher camera position for better target view
                "fov": 70.0,  # Field of view in degrees
                "enabled": True if name == "Sat1" else False  # Enable camera for first satellite
            },
            "targets": []  # Assigned targets
        }
        self.satellites.append(satellite)
        return satellite
        
    def create_main_ui(self):
        """Create the main user interface with FIXED layout"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create title bar
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        
        title_label = ttk.Label(title_frame, text="Spacecraft Constellation Fault Simulator - ENHANCED", 
                              style="Title.TLabel")
        title_label.pack(side=tk.LEFT)
        
        # FIXED: Add status indicator for improvements
        status_frame = ttk.Frame(title_frame)
        status_frame.pack(side=tk.LEFT, padx=20)
        
        improvements_label = ttk.Label(status_frame, text="✓ Enhanced Visibility  ✓ Better Orbits  ✓ 30min Default", 
                                     style="Info.TLabel", foreground="darkgreen")
        improvements_label.pack()
        
        # Add help button to title frame
        help_button = ttk.Button(title_frame, text="Help", command=self.show_general_help)
        help_button.pack(side=tk.RIGHT, padx=10)
        
        # Run button - make it prominent
        self.run_button = ttk.Button(title_frame, text="Run Enhanced Simulation", style="Run.TButton",
                                     command=self.run_simulation)
        self.run_button.pack(side=tk.RIGHT, padx=5)
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create tab frames
        self.constellation_frame = ttk.Frame(self.notebook)
        self.fault_frame = ttk.Frame(self.notebook)
        self.target_frame = ttk.Frame(self.notebook)
        self.visualization_frame = ttk.Frame(self.notebook)
        self.output_frame = ttk.Frame(self.notebook)
        self.plots_frame = ttk.Frame(self.notebook)  # Results tab
        
        # Add tab frames to notebook
        self.notebook.add(self.constellation_frame, text="Constellation")
        self.notebook.add(self.fault_frame, text="Fault Configuration")
        self.notebook.add(self.target_frame, text="Targets")
        self.notebook.add(self.visualization_frame, text="Visualization")
        self.notebook.add(self.output_frame, text="Output Settings")
        self.notebook.add(self.plots_frame, text="Results")
        
        # Create tab objects with FIXED versions
        self.constellation_tab = ConstellationTab(self, self.constellation_frame)
        self.fault_tab = FaultTab(self, self.fault_frame)
        self.target_tab = TargetTab(self, self.target_frame)
        self.visualization_tab = VisualizationTab(self, self.visualization_frame)
        self.output_tab = OutputTab(self, self.output_frame)
        
        # Create enhanced plots tab with zoom functionality
        self.create_enhanced_plots_tab()
        
        # FIXED: Status bar with enhanced information
        status_frame = ttk.Frame(main_frame, style="Status.TFrame")
        status_frame.pack(fill=tk.X, pady=(10, 0))
        
        self.status_label = ttk.Label(status_frame, text="Ready - Enhanced Configuration", style="Status.TLabel")
        self.status_label.pack(side=tk.LEFT, padx=5)
        
        # Add simulation time display to status bar
        sim_time_frame = ttk.Frame(status_frame)
        sim_time_frame.pack(side=tk.LEFT, padx=20)
        
        ttk.Label(sim_time_frame, text="Sim Time:", style='Status.TLabel').pack(side=tk.LEFT)
        self.sim_time_status = ttk.Label(sim_time_frame, text="30 min", 
                                       style='Status.TLabel', foreground='darkgreen')
        self.sim_time_status.pack(side=tk.LEFT, padx=5)
        
        # Bind simulation time change to update status
        self.simulation_time.trace('w', self.update_sim_time_status)
        
        # Add timestamp to status bar
        self.time_label = ttk.Label(status_frame, text=datetime.now().strftime("%Y-%m-%d %H:%M"), 
                                style='Status.TLabel')
        self.time_label.pack(side=tk.RIGHT, padx=10)
        
        # Update time every minute
        self.update_time()
        
        # Add menu bar
        self.create_menu_bar()
        
    def update_sim_time_status(self, *args):
        """Update simulation time display in status bar"""
        try:
            sim_time = self.simulation_time.get()
            # Convert to Vizard format for user reference
            hours = int(sim_time // 60)
            minutes = int(sim_time % 60)
            self.sim_time_status.config(text=f"{sim_time:.0f} min (Vizard: 00:{hours:02d}:{minutes:02d}:00:0)")
        except:
            self.sim_time_status.config(text="30 min")
    
    def create_enhanced_plots_tab(self):
        """Create enhanced plots tab with zoom functionality"""
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure
        
        # Main container
        main_frame = ttk.Frame(self.plots_frame, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Top controls
        controls_frame = ttk.Frame(main_frame)
        controls_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(controls_frame, text="Enhanced Simulation Results", style="Title.TLabel").pack(side=tk.LEFT)
        
        refresh_btn = ttk.Button(controls_frame, text="Refresh", 
                            command=self.refresh_plot_list)
        refresh_btn.pack(side=tk.RIGHT, padx=5)
        
        open_folder_btn = ttk.Button(controls_frame, text="Open Results Folder", 
                                    command=self.open_results_folder)
        open_folder_btn.pack(side=tk.RIGHT, padx=5)
        
        export_btn = ttk.Button(controls_frame, text="Export Plot", 
                               command=self.export_current_plot)
        export_btn.pack(side=tk.RIGHT, padx=5)
        
        # Split into two panels with better proportions
        split_frame = ttk.Frame(main_frame)
        split_frame.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - plot list (fixed width)
        list_frame = ttk.LabelFrame(split_frame, text="Available Plots", padding=10)
        list_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        list_frame.configure(width=320)
        
        # Scrollable list of plots
        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.plot_listbox = tk.Listbox(list_container, 
                                    selectmode=tk.SINGLE,
                                    yscrollcommand=scrollbar.set,
                                    font=('Segoe UI', 10),
                                    exportselection=False)
        self.plot_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.plot_listbox.yview)
        
        # Bind selection event
        self.plot_listbox.bind('<<ListboxSelect>>', self.on_plot_selected)
        
        # Right panel - plot display with matplotlib integration
        display_frame = ttk.LabelFrame(split_frame, text="Plot View", padding=10)
        display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))
        
        # Create matplotlib figure for interactive plotting
        self.plot_figure = Figure(figsize=(12, 8), dpi=100)
        self.plot_canvas = FigureCanvasTkAgg(self.plot_figure, display_frame)
        self.plot_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add matplotlib navigation toolbar for zoom/pan
        toolbar_frame = ttk.Frame(display_frame)
        toolbar_frame.pack(fill=tk.X, pady=(5, 0))
        
        self.plot_toolbar = NavigationToolbar2Tk(self.plot_canvas, toolbar_frame)
        self.plot_toolbar.update()
        
        # Plot controls frame
        plot_controls_frame = ttk.Frame(toolbar_frame)
        plot_controls_frame.pack(side=tk.RIGHT, padx=10)
        
        # Zoom controls
        zoom_in_btn = ttk.Button(plot_controls_frame, text="Zoom In", 
                                command=lambda: self.zoom_plot(1.5))
        zoom_in_btn.pack(side=tk.LEFT, padx=2)
        
        zoom_out_btn = ttk.Button(plot_controls_frame, text="Zoom Out", 
                                 command=lambda: self.zoom_plot(0.67))
        zoom_out_btn.pack(side=tk.LEFT, padx=2)
        
        reset_zoom_btn = ttk.Button(plot_controls_frame, text="Reset View", 
                                   command=self.reset_plot_view)
        reset_zoom_btn.pack(side=tk.LEFT, padx=2)
        
        # Info frame at the bottom
        info_frame = ttk.LabelFrame(main_frame, text="Plot Information", padding=10)
        info_frame.pack(fill=tk.X, pady=(10, 0))
        
        self.plot_info_text = scrolledtext.ScrolledText(info_frame, height=4, wrap=tk.WORD)
        self.plot_info_text.pack(fill=tk.BOTH, expand=True)
        self.plot_info_text.config(state=tk.DISABLED)
        
        # Initialize empty plot list
        self.refresh_plot_list()
    
    def refresh_plot_list(self):
        """Refresh the list of available plots"""
        self.plot_listbox.delete(0, tk.END)
        
        # Clear the plot info
        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)
        self.plot_info_text.insert(tk.END, "Select a plot to view details\n\nENHANCED: Plots now include better orbital analysis and target visibility data.")
        self.plot_info_text.config(state=tk.DISABLED)
        
        # Clear the plot canvas
        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        ax.text(0.5, 0.5, "No plot selected\nSelect a plot from the list to view\n\nENHANCED visualization available", 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.axis('off')
        self.plot_canvas.draw()
        
        # Check if plots directory exists
        if not os.path.exists(self.plots_dir):
            return
            
        # Get list of plot files
        plot_files = []
        for filename in os.listdir(self.plots_dir):
            if filename.endswith(('.png', '.jpg', '.jpeg', '.pdf', '.svg')):
                plot_files.append(filename)
                
        # Sort by modification time (newest first)
        plot_files.sort(key=lambda x: os.path.getmtime(os.path.join(self.plots_dir, x)), 
                      reverse=True)
                
        # Add to listbox
        for filename in plot_files:
            self.plot_listbox.insert(tk.END, filename)
            
        # Select first item if available
        if plot_files and (not self.latest_plots or self.latest_plots[0] not in plot_files):
            self.plot_listbox.selection_set(0)
            self.on_plot_selected(None)
        elif self.latest_plots:
            # Try to select the first plot from the latest simulation
            try:
                index = plot_files.index(self.latest_plots[0])
                self.plot_listbox.selection_set(index)
                self.on_plot_selected(None)
            except ValueError:
                # Latest plot not found, select first if available
                if plot_files:
                    self.plot_listbox.selection_set(0)
                    self.on_plot_selected(None)
    
    def on_plot_selected(self, event):
        """Handle plot selection event with matplotlib display"""
        selection = self.plot_listbox.curselection()
        if not selection:
            return
            
        index = selection[0]
        filename = self.plot_listbox.get(index)
        full_path = os.path.join(self.plots_dir, filename)
        
        # Update plot info
        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)
        
        # Get file info
        try:
            file_size = os.path.getsize(full_path) / 1024  # Size in KB
            mod_time = datetime.fromtimestamp(os.path.getmtime(full_path))
            
            # Parse information from filename
            plot_type = "Unknown"
            fault_type = "Unknown"
            
            # Try to extract info from filename
            parts = filename.split('_')
            if len(parts) >= 2:
                plot_type = parts[0].replace('Plot', '')
                
                if len(parts) >= 3:
                    fault_type = parts[1]
            
            info_text = f"File: {filename}\n"
            info_text += f"Type: {plot_type}\n"
            info_text += f"Fault: {fault_type}\n"
            info_text += f"Size: {file_size:.1f} KB\n"
            info_text += f"Created: {mod_time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            info_text += f"\nENHANCED: This plot includes improved orbital analysis with better target visibility data.\n"
            info_text += f"Use toolbar below to zoom, pan, and navigate the plot."
            
            self.plot_info_text.insert(tk.END, info_text)
        except Exception as e:
            self.plot_info_text.insert(tk.END, f"Error getting file info: {e}")
            
        self.plot_info_text.config(state=tk.DISABLED)
        
        # Display the plot using matplotlib for better interaction
        try:
            # Clear previous plot
            self.plot_figure.clear()
            
            # Load and display image
            if filename.endswith(('.png', '.jpg', '.jpeg')):
                # For image files, display directly
                img = plt.imread(full_path)
                ax = self.plot_figure.add_subplot(111)
                ax.imshow(img)
                ax.axis('off')  # Hide axes for image display
                ax.set_title(f"Enhanced Plot: {filename}", fontsize=14, pad=20)
                
            else:
                # For other formats, show info
                ax = self.plot_figure.add_subplot(111)
                ax.text(0.5, 0.5, f"Enhanced Plot: {filename}\n\nUse 'Export Plot' to save in different formats", 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.axis('off')
            
            # Update canvas
            try:
                self.plot_figure.tight_layout()
            except:
                pass  # Ignore layout warnings
            self.plot_canvas.draw()
            
            # Store current plot path for export
            self.current_plot_path = full_path
            
        except Exception as e:
            # Error loading plot
            self.plot_figure.clear()
            ax = self.plot_figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Error loading plot:\n{str(e)}", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.axis('off')
            self.plot_canvas.draw()

    def zoom_plot(self, factor):
        """Zoom the current plot by the given factor"""
        try:
            ax = self.plot_figure.gca()
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            
            # Calculate new limits
            xcenter = (xlim[0] + xlim[1]) / 2
            ycenter = (ylim[0] + ylim[1]) / 2
            
            xrange = (xlim[1] - xlim[0]) / factor
            yrange = (ylim[1] - ylim[0]) / factor
            
            ax.set_xlim(xcenter - xrange/2, xcenter + xrange/2)
            ax.set_ylim(ycenter - yrange/2, ycenter + yrange/2)
            
            self.plot_canvas.draw()
        except:
            pass  # Ignore errors for non-zoomable plots

    def reset_plot_view(self):
        """Reset the plot view to default"""
        try:
            ax = self.plot_figure.gca()
            ax.autoscale()
            self.plot_canvas.draw()
        except:
            # Reload the current plot
            self.on_plot_selected(None)

    def export_current_plot(self):
        """Export the currently selected plot"""
        if not hasattr(self, 'current_plot_path'):
            messagebox.showwarning("No Plot Selected", "Please select a plot first.")
            return
        
        # Ask user for save location
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("JPEG files", "*.jpg"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("All files", "*.*")
            ],
            title="Export Enhanced Plot"
        )
        
        if file_path:
            try:
                # Save the current figure
                self.plot_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Export Successful", f"Enhanced plot exported to:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Export Failed", f"Could not export plot:\n{str(e)}")

    def open_results_folder(self):
        """Open the results folder in the system file explorer"""
        if not os.path.exists(self.plots_dir):
            messagebox.showinfo("Folder Not Found", 
                             "The plots directory does not exist yet. Run a simulation first.")
            return
            
        # Use appropriate command based on platform
        try:
            if sys.platform == 'win32':
                os.startfile(self.plots_dir)
            elif sys.platform == 'darwin':  # macOS
                subprocess.run(['open', self.plots_dir])
            else:  # Linux
                subprocess.run(['xdg-open', self.plots_dir])
        except Exception as e:
            messagebox.showerror("Error Opening Folder", 
                              f"Could not open the results folder: {e}")

    def create_menu_bar(self):
        """Create the application menu bar"""
        menu_bar = tk.Menu(self.root)
        self.root.config(menu=menu_bar)
        
        # File menu
        file_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=file_menu)
        
        file_menu.add_command(label="New Enhanced Constellation", command=self.new_constellation)
        file_menu.add_separator()
        file_menu.add_command(label="Load Configuration", command=self.import_config)
        file_menu.add_command(label="Save Configuration", command=self.export_config)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Simulation menu
        sim_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Simulation", menu=sim_menu)
        
        sim_menu.add_command(label="Run Enhanced Simulation", command=self.run_simulation)
        sim_menu.add_separator()
        sim_menu.add_command(label="View Latest Results", command=lambda: self.notebook.select(self.plots_frame))
        sim_menu.add_command(label="Open Results Folder", command=self.open_results_folder)
        
        # Help menu
        help_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Help", menu=help_menu)
        
        help_menu.add_command(label="Enhanced Features Guide", command=self.show_general_help)
        help_menu.add_command(label="About Enhanced Version", command=self.show_about_dialog)

    def new_constellation(self):
        """Create a new enhanced constellation"""
        # Ask for confirmation if there are existing satellites
        if self.satellites:
            confirm = messagebox.askyesno("New Enhanced Constellation", 
                                        "This will clear the current configuration and create an enhanced constellation. Continue?")
            if not confirm:
                return
                
        # Clear current configuration
        self.satellites.clear()
        
        # FIXED: Add default satellites with enhanced parameters
        for i in range(4):
            true_anomaly = i * 90.0  # Space them 90 degrees apart
            sat_name = f"Sat{i+1}"  # Simplified naming
            self.add_satellite(sat_name, true_anomaly)
            
        # Reset target assignments with auto-assignment
        for i, target in enumerate(self.targets):
            target["assigned_to"] = [f"Sat{i+1}"] if i < len(self.satellites) else []
            
        # Update UI
        self.constellation_tab.update_satellite_listbox()
        self.update_satellite_dropdowns()
        self.update_target_assignments()
        
        # Select first satellite
        if self.satellites:
            self.constellation_tab.satellite_listbox.selection_set(0)
            self.constellation_tab.on_satellite_selected(None)
            
        self.add_log("Created new enhanced constellation with better orbital parameters and target assignments")
        
    def update_time(self):
        """Update the time display in the status bar"""
        self.time_label.config(text=datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.root.after(60000, self.update_time)  # Update every minute
        
    def update_status(self, message):
        """Update the status bar with a message"""
        self.status_label.config(text=message)
        self.root.update_idletasks()
        
    def add_log(self, message):
        """Add a message to the log"""
        try:
            log_text = self.output_tab.log_text
            log_text.config(state="normal")
            timestamp = datetime.now().strftime("%H:%M:%S")
            log_text.insert(tk.END, f"[{timestamp}] {message}\n")
            log_text.see(tk.END)
            log_text.config(state="disabled")
            self.root.update_idletasks()
            # Also log to file
            self.logger.info(message)
        except Exception as e:
            print(f"Error adding to log: {e}")
            print(f"Log message: {message}")
            
    def update_satellite_dropdowns(self):
        """Update all satellite dropdown menus"""
        try:
            self.fault_tab.update_satellite_dropdown()
        except:
            pass
            
        try:
            self.visualization_tab.update_satellite_dropdown()
        except:
            pass
            
        try:
            self.target_tab.update_satellite_dropdown()
        except:
            pass
            
    def update_target_assignments(self):
        """Update target assignments display"""
        try:
            self.target_tab.update_target_assignments()
        except:
            pass

    def show_general_help(self):
        """Show enhanced help about the simulator"""
        help_window = tk.Toplevel(self.root)
        help_window.title("Enhanced Spacecraft Fault Simulator Help")
        help_window.geometry("900x700")
        help_window.transient(self.root)
        help_window.grab_set()
        
        # Create notebook for help topics
        notebook = ttk.Notebook(help_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create tabs for different help sections
        overview_tab = ttk.Frame(notebook)
        enhanced_tab = ttk.Frame(notebook)
        simulation_tab = ttk.Frame(notebook)
        faults_tab = ttk.Frame(notebook)
        viz_tab = ttk.Frame(notebook)
        
        notebook.add(overview_tab, text="Overview")
        notebook.add(enhanced_tab, text="Enhanced Features")
        notebook.add(simulation_tab, text="Simulation")
        notebook.add(faults_tab, text="Fault Types")
        notebook.add(viz_tab, text="Visualization")
        
        # Create help content for each tab
        help_content = {
            "overview": """The Enhanced Spacecraft Constellation Fault Simulator allows you to model multiple satellites in orbit with improved visibility and better target observation capabilities.

Key improvements in this enhanced version:
• Better orbital altitudes (1629km) for improved camera visibility
• Enhanced target positioning (2km above Earth) for excellent visibility
• Simplified satellite naming (Sat1, Sat2, etc.)
• Optimized 30-minute default simulation time
• Improved camera positioning for ground target observation""",
            
            "enhanced": """ENHANCED FEATURES:

✓ ORBITAL IMPROVEMENTS:
- Default altitude: 1629km (better than previous 629km)
- Stable orbital parameters for realistic motion
- Better satellite spacing for constellation operations

✓ TARGET VISIBILITY:
- Targets positioned 2km above Earth surface
- 150km visibility markers in Vizard
- Auto-assignment of targets to satellites
- Enhanced color coding for easy identification

✓ CAMERA ENHANCEMENTS:
- Default camera height: 5m above spacecraft
- Optimal positioning for ground target observation
- Better field of view settings (70°)
- Enhanced Vizard compatibility

✓ USER INTERFACE:
- Simplified satellite names (Sat1, Sat2, etc.)
- 30-minute simulation time with Vizard format display
- Better status indicators and progress tracking
- Enhanced plot generation and analysis""",
            
            "simulation": """SIMULATION PARAMETERS:

• Simulation Time: Default 30 minutes (optimal for target analysis)
• Vizard Time Format: 00:30:00:00:0 (Days:Hours:Minutes:Seconds:Frames)
• Target Assignment: Ground locations visible to satellites with enhanced visibility
• Orbital Parameters: Configure satellite orbits with better default altitudes
• Fault Injection: Apply faults at specific times with improved analysis

ENHANCED SIMULATION FEATURES:
• Better orbital periods for target observation
• Improved satellite-to-target visibility
• Enhanced data logging and analysis
• Optimized visualization output for Vizard""",
            
            "fault": """The enhanced simulator supports improved fault analysis:

• Friction: Additional friction in reaction wheels with better visualization
• Power Limit: Electrical power restrictions with enhanced impact analysis  
• Encoder: Speed measurement errors with improved error modeling
• Battery: Power system degradation with comprehensive power analysis

ENHANCED FAULT FEATURES:
• Better fault impact visualization
• Improved timing controls
• Enhanced plot generation
• Better correlation with orbital motion""",
            
            "visualization": """ENHANCED VISUALIZATION:

The simulator outputs enhanced visualization files for Vizard:
• Targets positioned for maximum visibility (2km above Earth)
• 150km visibility markers for easy identification  
• Better camera positioning for ground observation
• Enhanced satellite motion with stable orbits

VIZARD IMPROVEMENTS:
• Better binary file organization
• Enhanced target visibility
• Improved camera angles
• Better time format handling (00:30:00:00:0 for 30 minutes)

Use the Results tab to view plots with enhanced zoom functionality."""
        }
        
        for (tab, content), text in zip(help_content.items(), 
                                       ["Overview", "Enhanced Features", "Simulation", "Fault Types", "Visualization"]):
            frame_index = list(help_content.keys()).index(tab)
            frame = [overview_tab, enhanced_tab, simulation_tab, faults_tab, viz_tab][frame_index]
            
            text_frame = ttk.Frame(frame, padding=10)
            text_frame.pack(fill=tk.BOTH, expand=True)
            
            text_widget = scrolledtext.ScrolledText(text_frame, wrap=tk.WORD, font=('Segoe UI', 10))
            text_widget.pack(fill=tk.BOTH, expand=True)
            text_widget.insert(tk.END, content)
            text_widget.config(state=tk.DISABLED)
        
        # Add close button
        ttk.Button(help_window, text="Close", command=help_window.destroy).pack(pady=10)
        
    def show_about_dialog(self):
        """Show the about dialog for enhanced version"""
        about_window = tk.Toplevel(self.root)
        about_window.title("About Enhanced Spacecraft Fault Simulator")
        about_window.geometry("500x400")
        about_window.transient(self.root)
        about_window.grab_set()
        
        # Create a frame for the content
        frame = ttk.Frame(about_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Add title
        ttk.Label(frame, text="Spacecraft Constellation Fault Simulator", 
                style='Title.TLabel').pack(pady=(0, 5))
        
        ttk.Label(frame, text="ENHANCED VERSION", 
                font=('Segoe UI', 12, 'bold'), foreground='darkgreen').pack(pady=(0, 10))
        
        # Version info
        ttk.Label(frame, text="Version 2.1.0 - Enhanced Visibility & Orbits", 
                font=('Segoe UI', 10, 'italic')).pack()
        
        # Description
        description = ttk.Label(frame, text=
                              "A high-fidelity simulation tool for spacecraft reaction wheel\n"
                              "fault analysis in multi-satellite constellations.\n\n"
                              "ENHANCED FEATURES:\n"
                              "✓ Better orbital altitudes (1629km) for improved visibility\n"
                              "✓ Enhanced target positioning (2km above Earth surface)\n"
                              "✓ Optimized camera placement for ground target observation\n"
                              "✓ 30-minute simulations with proper Vizard time formatting\n"
                              "✓ Simplified satellite naming and improved user interface\n"
                              "✓ Enhanced plot analysis with zoom functionality",
                              justify=tk.LEFT)
        description.pack(pady=10)
        
        # Built with info
        ttk.Label(frame, text="Built with:", font=('Segoe UI', 10, 'bold')).pack(pady=(10, 5))
        ttk.Label(frame, text="Python, Tkinter, Matplotlib, Basilisk, NumPy", 
                justify=tk.CENTER).pack()
        
        # Copyright
        ttk.Label(frame, text="© 2025 Enhanced Spacecraft Dynamics Lab", 
                font=('Segoe UI', 8)).pack(side=tk.BOTTOM, pady=10)
        
        # Close button
        ttk.Button(frame, text="Close", command=about_window.destroy).pack(side=tk.BOTTOM, pady=10)

    def run_simulation(self):
        """Run the enhanced simulation with the current configuration"""
        # Check if simulation is already running
        if self.is_running:
            messagebox.showinfo("Simulation Running", "An enhanced simulation is already running. Please wait for it to complete.")
            return
            
        # Update settings from UI
        try:
            # Update output settings from the output tab
            if hasattr(self, 'output_tab') and hasattr(self.output_tab, 'update_settings_from_ui'):
                self.output_tab.update_settings_from_ui()
        except Exception as e:
            self.add_log(f"Warning: Could not update settings from UI: {e}")
        
        # Create required directories if they don't exist
        for dir_path in [self.output_dir, self.plots_dir, self.bin_dir]:
            try:
                os.makedirs(dir_path, exist_ok=True)
            except Exception as e:
                self.add_log(f"Warning: Could not create directory {dir_path}: {e}")
        
        # Set running state
        self.is_running = True
        self.run_button.config(state="disabled")
        self.update_status("Running enhanced simulation...")
        
        # Start simulation in a separate thread
        self.add_log("Starting ENHANCED simulation with improved parameters...")
        simulation_thread = threading.Thread(target=self.run_simulation_process)
        simulation_thread.daemon = True
        simulation_thread.start()

    def run_simulation_process(self):
        """Run the enhanced simulation process"""
        try:
            # Create a simulation configuration based on UI settings
            config = SimulationConfig()
            
            # Set up general simulation parameters with proper time
            config.simulation_time = self.simulation_time.get()  # This is now properly 30 minutes
            config.show_plots = self.show_plots.get()
            config.save_plots = self.save_plots.get()
            config.save_binary = self.save_binary.get()
            
            # Use just the base filename, not a path
            binary_filename = os.path.basename(self.binary_filename)
            config.binary_filename = binary_filename
            self.add_log(f"Using binary filename: {binary_filename}")
            self.add_log(f"ENHANCED simulation duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)")
            
            # Convert to Vizard format for user info
            hours = int(config.simulation_time // 60)
            minutes = int(config.simulation_time % 60)
            self.add_log(f"Vizard time format: 00:{hours:02d}:{minutes:02d}:00:0")
            
            # Configure spacecraft for the simulation
            config.spacecraft_list = []
            
            # Add all satellites to the simulation
            for i, satellite in enumerate(self.satellites):
                # Create spacecraft configuration
                spacecraft_config = {
                    "name": satellite["name"],
                    "orbit": {
                        "a": satellite["orbit"]["a"],  # Already in km, ENHANCED altitude
                        "e": satellite["orbit"]["e"],
                        "i": satellite["orbit"]["i"],  # Already in degrees
                        "Omega": satellite["orbit"]["Omega"],  # Already in degrees
                        "omega": satellite["orbit"]["omega"],  # Already in degrees
                        "f": satellite["orbit"]["f"]  # Already in degrees
                    },
                    "fault": {
                        "enabled": satellite["fault"]["enabled"],
                        "type": satellite["fault"]["type"],
                        "magnitude": satellite["fault"]["magnitude"],
                        "wheel": satellite["fault"]["wheel"],
                        "time": satellite["fault"]["time"],  # Already in minutes
                        "periodic": satellite["fault"]["periodic"]
                    },
                    "camera": {
                        "enabled": satellite["camera"]["enabled"],
                        "position": satellite["camera"]["position"],  # ENHANCED position
                        "fov": satellite["camera"]["fov"]
                    }
                }
                
                # Add to configuration list
                config.spacecraft_list.append(spacecraft_config)
                
                # Log spacecraft details with enhanced info
                altitude = satellite["orbit"]["a"] - 6371  # Calculate altitude
                self.add_log(f"Added ENHANCED spacecraft: {satellite['name']} (altitude: {altitude:.0f}km)")
                if satellite["fault"]["enabled"]:
                    self.add_log(f"  - Fault: {satellite['fault']['type']}, magnitude: {satellite['fault']['magnitude']}")
                    self.add_log(f"  - Fault time: {satellite['fault']['time']} minutes ({satellite['fault']['time'] * 60:.0f} seconds)")
                if satellite["camera"]["enabled"]:
                    self.add_log(f"  - ENHANCED camera at position: {satellite['camera']['position']}")
            
            # Set active camera configuration
            active_camera_sat = self.visualization_tab.get_active_camera_satellite()
            if active_camera_sat:
                config.active_camera_name = active_camera_sat["name"]
                config.camera_position = active_camera_sat["camera"]["position"]  # ENHANCED position
                config.camera_fov = active_camera_sat["camera"]["fov"]
                self.add_log(f"Using ENHANCED camera from satellite: {active_camera_sat['name']}")
            else:
                # Fallback to visualization tab camera position
                camera_position = self.visualization_tab.get_camera_position()
                if camera_position:
                    config.camera_position = camera_position
                    config.active_camera_name = None
            
            # Configure targets with enhanced visibility
            config.targets = []
            for target in self.targets:
                if target["assigned_to"]:  # Only include targets assigned to satellites
                    # Convert hex color to proper Basilisk color name
                    color = target["color"].replace("#", "")
                    # Map common hex colors to names that Basilisk recognizes
                    color_map = {
                        "FF0000": "red",
                        "0000FF": "blue", 
                        "00FF00": "green",
                        "FFFF00": "yellow",
                        "800080": "purple",
                        "FFA500": "orange",
                        "00FFFF": "cyan",
                        "FF00FF": "magenta",
                        "008000": "green",
                        "000080": "blue",
                        "800000": "red"
                    }
                    # Use named color if available, otherwise default to red
                    display_color = color_map.get(color.upper(), "red")
                    
                    config.targets.append(TargetDefinition(
                        target["name"], 
                        target["lat"], 
                        target["lon"], 
                        color=display_color
                    ))
                    
                    # Also store which satellites this target is assigned to
                    config.targets[-1].assigned_to = target["assigned_to"]
                    
                    self.add_log(f"Added ENHANCED target: {target['name']} assigned to {', '.join(target['assigned_to'])} (color: {display_color})")
            
            # Run the enhanced simulation directly
            self.add_log("Starting direct ENHANCED simulation...")
            
            scenario, viz, figureList, output_dir = run_custom_simulation(config)
            
            # Save plots - handle Figure objects properly
            self.latest_plots = []  # Clear previous plot list
            if self.save_plots.get() and figureList:
                self.add_log(f"Saving {len(figureList)} ENHANCED plots to {self.plots_dir}...")
                # Create timestamp for unique filenames
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                
                for name, fig in figureList.items():
                    try:
                        # If figureList contains Figure objects, we need to save them directly
                        if hasattr(fig, 'savefig'):  # It's a matplotlib Figure object
                            # Create a more descriptive filename with timestamp
                            plot_filename = f"Enhanced_{name}_{timestamp}.png"
                            plot_path = os.path.join(self.plots_dir, plot_filename)
                            
                            # Save the figure with high resolution
                            fig.savefig(plot_path, bbox_inches='tight', dpi=300)
                            self.add_log(f"  - Saved ENHANCED plot: {plot_filename}")
                            self.latest_plots.append(plot_filename)
                        else:  # It's a path (original expectation)
                            import shutil
                            dest_path = os.path.join(self.plots_dir, os.path.basename(fig))
                            shutil.copy2(fig, dest_path)
                            self.add_log(f"  - Copied plot: {os.path.basename(fig)}")
                            self.latest_plots.append(os.path.basename(fig))
                    except Exception as e:
                        self.add_log(f"Error saving plot {name}: {e}")
            
            # Check for binary files
            binary_found = False
            vizard_path_found = None
            if self.save_binary.get():
                # Check all possible locations for the binary file
                vizard_paths = [
                    os.path.join(self.VIZ_DIR, f"{binary_filename}_UnityViz.bin"),
                    os.path.join(self.VIZ_DIR, "_VizFiles", f"{binary_filename}_UnityViz.bin")
                ]
                
                for vizard_path in vizard_paths:
                    if os.path.exists(vizard_path):
                        binary_found = True
                        vizard_path_found = vizard_path
                        file_size = os.path.getsize(vizard_path) / (1024*1024)  # Size in MB
                        self.add_log(f"ENHANCED binary file created: {vizard_path} ({file_size:.2f} MB)")
                        self.add_log("To view this ENHANCED file in Vizard:")
                        self.add_log("1. Open Vizard application")
                        self.add_log(f"2. Load the binary file: {vizard_path}")
                        self.add_log(f"3. Simulation duration: 00:{hours:02d}:{minutes:02d}:00:0")
                        self.add_log("4. Targets are positioned 2km above Earth for MAXIMUM visibility")
                        self.add_log("5. Camera optimized for ground target observation")
                        break
                
            # Update the plot list in the Results tab
            self.refresh_plot_list()
            
            # Switch to the Results tab if plots were created
            if self.latest_plots and self.save_plots.get():
                self.notebook.select(self.plots_frame)
            
            # Show enhanced results summary
            summary = f"""
ENHANCED SIMULATION COMPLETED SUCCESSFULLY!

Duration: {config.simulation_time} minutes ({config.simulation_time * 60:.0f} seconds)
Vizard Format: 00:{hours:02d}:{minutes:02d}:00:0
Results saved to: {output_dir}
Spacecraft: {len(config.spacecraft_list)} satellites simulated with ENHANCED parameters

ENHANCED FEATURES APPLIED:
✓ Better orbital altitudes (1629km) for improved visibility
✓ Enhanced target positioning (2km above Earth surface)
✓ Optimized camera placement for ground observation
✓ Improved satellite motion stability
"""
            
            if config.active_camera_name:
                summary += f"\nActive ENHANCED camera: {config.active_camera_name}\n"

            if config.targets:
                summary += f"\nTargets: {len(config.targets)} locations with ENHANCED visibility\n"
                
            if self.save_plots.get() and figureList:
                summary += f"\nENHANCED plots saved to: {self.plots_dir}\n"
                summary += f"You can view the enhanced plots in the Results tab.\n"
                
            if self.save_binary.get() and binary_found:
                summary += f"\nENHANCED visualization binary saved to: {vizard_path_found}\n"
                summary += f"You can open this file in Vizard for 3D visualization.\n\n"
                summary += f"ENHANCED Vizard Features:\n"
                summary += f"- Simulation time: 00:{hours:02d}:{minutes:02d}:00:0\n"
                summary += f"- Targets positioned 2km above Earth for maximum visibility\n"
                summary += f"- 150km visibility markers for easy identification\n"
                summary += f"- Camera optimized for ground target observation\n"
                summary += f"- Stable satellite orbital motion (no unwanted spinning)"
            
            # Show the enhanced message box after the simulation completes
            messagebox.showinfo("Enhanced Simulation Complete", summary)
                
        except Exception as e:
            self.add_log(f"Error running ENHANCED simulation: {str(e)}")
            import traceback
            self.add_log(f"Traceback: {traceback.format_exc()}")
            self.update_status("Enhanced simulation error")
            
            # Show error message box
            messagebox.showerror("Enhanced Simulation Error", f"An error occurred during the enhanced simulation:\n\n{str(e)}")
        finally:
            self.is_running = False
            self.run_button.config(state="normal")

    def export_config(self):
        """Export the current enhanced configuration to a JSON file"""
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".json",
                                                     filetypes=[("JSON files", "*.json")],
                                                     title="Save Enhanced Configuration")
            if file_path:
                # Get active satellite data
                active_satellite_index = self.fault_tab.get_active_satellite_index()
                if active_satellite_index is not None and 0 <= active_satellite_index < len(self.satellites):
                    satellite = self.satellites[active_satellite_index]
                else:
                    satellite = self.satellites[0] if self.satellites else None
                
                # Create enhanced config dictionary
                config_dict = {
                    "version": "Enhanced 2.1.0",
                    "simulation_time": self.simulation_time.get(),
                    "show_plots": self.show_plots.get(),
                    "save_plots": self.save_plots.get(),
                    "save_binary": self.save_binary.get(),
                    "binary_filename": self.binary_filename,
                    "satellites": self.satellites,
                    "targets": self.targets,
                    "active_satellite_index": active_satellite_index if active_satellite_index is not None else 0,
                    "enhanced_features": {
                        "better_orbits": True,
                        "enhanced_targets": True,
                        "optimized_camera": True,
                        "simplified_naming": True
                    }
                }
                
                # Save to file with nice formatting
                with open(file_path, 'w') as f:
                    json.dump(config_dict, f, indent=4)
                    
                self.update_status(f"Enhanced configuration saved to {os.path.basename(file_path)}")
                messagebox.showinfo("Export Success", f"Enhanced configuration saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))
            self.update_status("Enhanced export failed")
            
    def import_config(self):
        """Import enhanced configuration from a JSON file"""
        try:
            file_path = filedialog.askopenfilename(defaultextension=".json",
                                                filetypes=[("JSON files", "*.json")],
                                                title="Open Enhanced Configuration")
            if file_path:
                with open(file_path, 'r') as f:
                    config_dict = json.load(f)
                
                # Update simulation parameters
                self.simulation_time.set(config_dict.get("simulation_time", 30.0))  # Default to 30 minutes
                self.show_plots.set(config_dict.get("show_plots", True))
                self.save_plots.set(config_dict.get("save_plots", True))
                self.save_binary.set(config_dict.get("save_binary", True))
                self.binary_filename = config_dict.get("binary_filename", "saib_fault_viz")
                
                # Update satellites
                if "satellites" in config_dict:
                    self.satellites = config_dict["satellites"]
                    
                    # Ensure camera has enabled field and enhanced positioning
                    for sat in self.satellites:
                        if "camera" in sat:
                            if "enabled" not in sat["camera"]:
                                sat["camera"]["enabled"] = False
                            # Upgrade camera position if it's the old default
                            if sat["camera"]["position"] == [0.0, 0.0, 2.0]:
                                sat["camera"]["position"] = [0.0, 0.0, 5.0]  # Enhanced position
                                
                    self.update_satellite_dropdowns()
                
                # Update targets
                if "targets" in config_dict:
                    self.targets = config_dict["targets"]
                    # Ensure priority field
                    for target in self.targets:
                        if "priority" not in target:
                            target["priority"] = 1
                            
                    self.update_target_assignments()
                
                # Set active satellite
                active_satellite_index = config_dict.get("active_satellite_index", 0)
                if 0 <= active_satellite_index < len(self.satellites):
                    self.current_satellite_index.set(active_satellite_index)
                    self.fault_tab.set_active_satellite(active_satellite_index)
                
                # Check if this is an enhanced configuration
                is_enhanced = config_dict.get("version", "").startswith("Enhanced")
                if is_enhanced:
                    self.update_status(f"Enhanced configuration loaded from {os.path.basename(file_path)}")
                    messagebox.showinfo("Import Success", "Enhanced configuration loaded successfully")
                else:
                    self.update_status(f"Configuration loaded and enhanced from {os.path.basename(file_path)}")
                    messagebox.showinfo("Import Success", "Configuration loaded and automatically enhanced")
                    
        except Exception as e:
            messagebox.showerror("Import Failed", str(e))
            self.update_status("Enhanced import failed")


# Main entry point
if __name__ == "__main__":
    root = tk.Tk()
    app = SatelliteSimulatorApp(root)
    root.mainloop()