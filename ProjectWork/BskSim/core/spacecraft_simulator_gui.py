#!/usr/bin/env python
"""
spacecraft_simulator_gui.py

Main GUI application for the Basilisk Spacecraft Constellation Fault Simulator.
Provides a modular interface for configuring and running spacecraft simulations.
"""
import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import threading
import time
import subprocess
import json
from datetime import datetime
import logging
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend 
from PIL import ImageTk, Image

# Set up paths
import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
CORE_DIR = os.path.join(ROOT_DIR, 'core')
PLOTS_DIR = os.path.join(ROOT_DIR, 'plots')
ASSETS_DIR = os.path.join(ROOT_DIR, 'assets')

sys.path.extend([ROOT_DIR, CORE_DIR])

# Import simulation functionality
from spacecraft_simulation import SimulationConfig, TargetDefinition, run_custom_simulation

# Import GUI modules
from gui_tab import (
    ConstellationTab,
    FaultTab,
    TargetTab,
    VisualizationTab,
    OutputTab,
    ResultsTab
)

# Import help content
from rw_fault_help import get_general_help

# Import style
from style import setup_style


class SatelliteSimulatorApp:
    """Main application class for the Spacecraft Fault Simulator GUI"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Spacecraft Constellation Fault Simulator")
        self.root.geometry("1200x900")
        self.root.minsize(1000, 800)
        
        # Setup logging
        self._setup_logging()
        
        # Setup directories
        self._setup_directories()
        
        # Setup style
        self.style = setup_style(root)

        # Handle window close properly
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Initialize data
        self._initialize_data()
        
        # Create UI
        self._create_ui()
        
        # Update initial state
        self.update_status("Ready")
        self.update_status_counts()

    def _on_closing(self):
        """Handle window closing"""
        try:
            # Stop any running simulation
            self.is_running = False
            
            # Close any matplotlib figures
            import matplotlib.pyplot as plt
            plt.close('all')
        except:
            pass
        
        # Destroy the window
        self.root.quit()
            
    def _setup_logging(self):
        """Configure application logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('simulator')
        
    def _setup_directories(self):
        """Set up required directories"""
        self.ROOT_DIR = ROOT_DIR
        self.FAULTS_DIR = os.path.join(self.ROOT_DIR, 'faults')
        self.MODELS_DIR = os.path.join(self.ROOT_DIR, 'models')
        self.PLOTTING_DIR = os.path.join(self.ROOT_DIR, 'plots')
        self.LOGS_DIR = os.path.join(self.ROOT_DIR, 'logs')
        self.VIZ_DIR = os.path.join(self.ROOT_DIR, 'Vizfile')
        
        # Create directories
        for dir_path in [self.FAULTS_DIR, self.MODELS_DIR, self.PLOTTING_DIR, 
                        self.LOGS_DIR, self.VIZ_DIR]:
            os.makedirs(dir_path, exist_ok=True)
            
        # Create _VizFiles subdirectory
        viz_files_dir = os.path.join(self.VIZ_DIR, "_VizFiles")
        os.makedirs(viz_files_dir, exist_ok=True)
        
    def _initialize_data(self):
        """Initialize application data structures"""
        # Paths
        self.output_dir = self.LOGS_DIR
        self.plots_dir = self.PLOTTING_DIR
        self.bin_dir = self.VIZ_DIR
        self.binary_filename = "spacecraft_viz"
        
        # Simulation parameters
        self.simulation_time = tk.DoubleVar(value=30.0)
        self.show_plots = tk.BooleanVar(value=True)
        self.save_plots = tk.BooleanVar(value=True)
        self.save_binary = tk.BooleanVar(value=True)
        
        
        # Satellite data
        self.satellites = []
        self._initialize_default_satellites()
        self.current_satellite_index = tk.IntVar(value=0)
        
        # Target data
        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "#FF0000", 
             "priority": 3, "assigned_to": ["Satellite1"]},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "#0000FF", 
             "priority": 2, "assigned_to": ["Satellite2"]},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "#00FF00", 
             "priority": 2, "assigned_to": ["Satellite3"]},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "#FFFF00", 
             "priority": 1, "assigned_to": ["Satellite4"]}
        ]
        
        # State
        self.is_running = False
        self.latest_plots = []
        
    def _initialize_default_satellites(self):
        """Create default satellite constellation"""
        for i in range(4):
            self.add_satellite(f"Satellite{i+1}", i * 90.0, altitude=600)
            
    def _create_ui(self):
        """Create the main user interface"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title bar
        self._create_title_bar(main_frame)
        
        # Notebook for tabs
        self._create_notebook(main_frame)
        
        # Status bar
        self._create_status_bar(main_frame)
        
        # Menu bar
        self._create_menu_bar()
        
        # Update time periodically
        self.update_time()
        
    def _create_title_bar(self, parent):
        """Create the application title bar"""
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        title_frame.columnconfigure(1, weight=1)
        
        # Logo
        try:
            logo_path = os.path.join(CORE_DIR, "swinburne.png")
            if os.path.exists(logo_path):
                logo = Image.open(logo_path)
                logo = logo.resize((100, 50), Image.Resampling.LANCZOS)
                self.logo_photo = ImageTk.PhotoImage(logo)
                ttk.Label(title_frame, image=self.logo_photo).grid(
                    row=0, column=0, sticky="w", padx=(5, 10)
                )
        except Exception as e:
            self.logger.warning(f"Could not load logo: {e}")
            
        # Title
        ttk.Label(
            title_frame, 
            text="Spacecraft Constellation Fault Simulator",
            font=('Segoe UI', 16, 'bold'),
            anchor="center"
        ).grid(row=0, column=1)
        
        # Buttons
        btn_frame = ttk.Frame(title_frame)
        btn_frame.grid(row=0, column=2, sticky="e", padx=10)
        
        ttk.Button(btn_frame, text="Help", command=self.show_help).pack(
            side=tk.RIGHT, padx=5
        )
        
        self.run_button = ttk.Button(
            btn_frame, 
            text="Run Simulation",
            style="Run.TButton",
            command=self.run_simulation
        )
        self.run_button.pack(side=tk.RIGHT, padx=5)
        
    def _create_notebook(self, parent):
        """Create the tabbed interface"""
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create tab frames
        self.constellation_frame = ttk.Frame(self.notebook)
        self.fault_frame = ttk.Frame(self.notebook)
        self.target_frame = ttk.Frame(self.notebook)
        self.visualization_frame = ttk.Frame(self.notebook)
        self.output_frame = ttk.Frame(self.notebook)
        self.results_frame = ttk.Frame(self.notebook)
        
        # Add tabs
        self.notebook.add(self.constellation_frame, text="Constellation")
        self.notebook.add(self.fault_frame, text="Fault Configuration")
        self.notebook.add(self.target_frame, text="Targets")
        self.notebook.add(self.visualization_frame, text="Visualization")
        self.notebook.add(self.output_frame, text="Output Settings")
        self.notebook.add(self.results_frame, text="Results")
        
        # Create tab objects
        self.constellation_tab = ConstellationTab(self, self.constellation_frame)
        self.fault_tab = FaultTab(self, self.fault_frame)
        self.target_tab = TargetTab(self, self.target_frame)
        self.visualization_tab = VisualizationTab(self, self.visualization_frame)
        self.output_tab = OutputTab(self, self.output_frame)
        self.results_tab = ResultsTab(self, self.results_frame)
        
    def _create_status_bar(self, parent):
        """Create the status bar"""
        status_frame = ttk.Frame(parent, style="Status.TFrame")
        status_frame.pack(fill=tk.X, pady=(10, 0))
        
        # Status message
        self.status_label = ttk.Label(status_frame, text="Ready", style="Status.TLabel")
        self.status_label.pack(side=tk.LEFT, padx=5)
        
        # Satellite count
        self.sat_count_status = ttk.Label(
            status_frame, 
            text="0 satellites",
            style='Status.TLabel',
            foreground='blue'
        )
        self.sat_count_status.pack(side=tk.LEFT, padx=20)
        
        # Target count
        self.target_status = ttk.Label(
            status_frame,
            text="0/0 targets assigned",
            style='Status.TLabel',
            foreground='green'
        )
        self.target_status.pack(side=tk.LEFT, padx=20)
        
        # Simulation time
        sim_time_frame = ttk.Frame(status_frame)
        sim_time_frame.pack(side=tk.LEFT, padx=20)
        
        ttk.Label(sim_time_frame, text="Sim Time:", style='Status.TLabel').pack(side=tk.LEFT)
        self.sim_time_status = ttk.Label(
            sim_time_frame,
            text="30 min",
            style='Status.TLabel',
            foreground='darkgreen'
        )
        self.sim_time_status.pack(side=tk.LEFT, padx=5)
        
        # Bind simulation time change
        self.simulation_time.trace('w', self._update_sim_time_status)
        
        # Current time
        self.time_label = ttk.Label(
            status_frame,
            text=datetime.now().strftime("%Y-%m-%d %H:%M"),
            style='Status.TLabel'
        )
        self.time_label.pack(side=tk.RIGHT, padx=10)


    

    def _create_menu_bar(self):
        """Create the application menu bar"""
        menu_bar = tk.Menu(self.root)
        self.root.config(menu=menu_bar)
        
        # File menu
        file_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="New Constellation", command=self.new_constellation)
        file_menu.add_separator()
        file_menu.add_command(label="Load Configuration", command=self.import_config)
        file_menu.add_command(label="Save Configuration", command=self.export_config)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Simulation menu
        sim_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Simulation", menu=sim_menu)
        sim_menu.add_command(label="Run Simulation", command=self.run_simulation)
        sim_menu.add_separator()
        sim_menu.add_command(label="View Results", command=self._view_results)
        sim_menu.add_command(label="Open Results Folder", command=self.open_results_folder)
        
        # Help menu
        help_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_command(label="About", command=self.show_about)
        
    # Public methods
    def add_satellite(self, name, true_anomaly=0.0, altitude=600):
        """Add a new satellite to the constellation"""
        satellite = {
            "name": name,
            "orbit": {
                "a": 6371 + altitude,
                "e": 0.01,
                "i": 55.0,
                "Omega": 45.0,
                "omega": 30.0,
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
                "position": [0.0, 0.0, 15.0],
                "fov": 80.0,
                "enabled": name == "Satellite1"
            },
            "targets": [],
            "orbit_name": "Default Orbit"
        }
        self.satellites.append(satellite)
        return satellite
        
    def update_status(self, message):
        """Update the status bar message"""
        self.status_label.config(text=message)
        self.root.update_idletasks()
        
    def update_status_counts(self):
        """Update the status bar counts"""
        try:
            # Satellite count
            self.sat_count_status.config(text=f"{len(self.satellites)} satellites")
            
            # Target count
            assigned = len([t for t in self.targets if t.get('assigned_to', [])])
            total = len(self.targets)
            self.target_status.config(text=f"{assigned}/{total} targets assigned")
        except:
            pass
            
    def add_log(self, message):
        """Add a message to the application log"""
        try:
            self.output_tab.add_log_entry(message)  # Changed from add_log_message to add_log_entry
            self.logger.info(message)
        except Exception as e:
            print(f"Log: {message}")
            print(f"Error adding to log: {e}")
            
    def update_satellite_dropdowns(self):
        """Update all satellite dropdown menus"""
        try:
            self.fault_tab.update_satellite_dropdown()
            self.visualization_tab.update_satellite_dropdown()
            self.target_tab.update_satellite_dropdown()
        except:
            pass
            
    def update_target_assignments(self):
        """Update target assignment displays"""
        try:
            self.target_tab.update_target_assignments()
        except:
            pass
            
    def update_time(self):
        """Update the time display"""
        self.time_label.config(text=datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.root.after(60000, self.update_time)  # Update every minute
        
    # Simulation methods
    def run_simulation(self):
        """Run the simulation"""
        if self.is_running:
            messagebox.showinfo("Simulation Running", 
                              "A simulation is already running. Please wait.")
            return
            
        # Update settings from UI
        try:
            if hasattr(self.output_tab, 'update_settings_from_ui'):
                self.output_tab.update_settings_from_ui()
        except Exception as e:
            self.add_log(f"Warning: Could not update settings: {e}")
            
        # Set state
        self.is_running = True
        self.run_button.config(state="disabled")
        self.update_status("Running simulation...")
        
        # Start in thread
        self.add_log("Starting simulation...")
        thread = threading.Thread(target=self._run_simulation_process)
        thread.daemon = True
        thread.start()
            

    def _run_simulation_process(self):
        """Run the simulation process"""
        try:
            # Create configuration
            config = SimulationConfig()
            
            # General parameters
            config.simulation_time = self.simulation_time.get()
            config.show_plots = self.show_plots.get()
            config.save_plots = self.save_plots.get()
            config.save_binary = self.save_binary.get()
            config.binary_filename = os.path.basename(self.binary_filename)
            
            # Spacecraft configuration
            config.spacecraft_list = []
            for sat in self.satellites:
                config.spacecraft_list.append({
                    "name": sat["name"],
                    "orbit": sat["orbit"],
                    "fault": sat["fault"],
                    "camera": sat["camera"]
                })
                
            # Camera configuration
            active_camera = self.visualization_tab.get_active_camera_satellite()
            if active_camera:
                config.active_camera_name = active_camera["name"]
                config.camera_position = active_camera["camera"]["position"]
                config.camera_fov = active_camera["camera"]["fov"]
                
            # Target configuration
            config.targets = []
            for target in self.targets:
                if target["assigned_to"]:
                    t = TargetDefinition(
                        target["name"],
                        target["lat"],
                        target["lon"],
                        self._convert_color(target["color"])
                    )
                    t.assigned_to = target["assigned_to"]
                    config.targets.append(t)
                    
            # Log summary
            self.add_log(f"Simulation: {len(config.spacecraft_list)} satellites, "
                        f"{len(config.targets)} targets, {config.simulation_time} minutes")
            
            # Run simulation
            scenario, viz, figureList, output_dir = run_custom_simulation(config)
            
            # Handle results
            self._handle_simulation_results(figureList, output_dir, config)
            
            # Show completion message
            self._show_completion_message(config, output_dir)
            
        except Exception as e:
            self.add_log(f"Simulation error: {e}")
            self.update_status("Simulation error")
            messagebox.showerror("Simulation Error", f"Error: {str(e)}")
            
        finally:
            # Clean up matplotlib in thread-safe way
            try:
                import matplotlib.pyplot as plt
                plt.close('all')
            except:
                pass
                
            self.is_running = False
            self.run_button.config(state="normal")
            self.update_status("Ready")
            
    def _convert_color(self, hex_color):
        """Convert hex color to name for Basilisk"""
        color_map = {
            "#FF0000": "red",
            "#0000FF": "blue",
            "#00FF00": "green",
            "#FFFF00": "yellow",
            "#800080": "purple",
            "#FFA500": "orange",
            "#00FFFF": "cyan",
            "#FF00FF": "magenta"
        }
        return color_map.get(hex_color.upper(), "red")
        
    def _handle_simulation_results(self, figureList, output_dir, config):
        """Handle simulation results"""
        self.latest_plots = []
        
        if self.save_plots.get() and figureList:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.add_log(f"Saving {len(figureList)} plots...")
            
            # Import matplotlib here to ensure correct backend
            import matplotlib.pyplot as plt
            
            for name, fig in figureList.items():
                try:
                    if hasattr(fig, 'savefig'):
                        filename = f"{name}_{timestamp}.png"
                        path = os.path.join(self.plots_dir, filename)
                        fig.savefig(path, bbox_inches='tight', dpi=300)
                        self.latest_plots.append(filename)
                        self.add_log(f"Saved: {filename}")
                        
                        # Close figure immediately after saving
                        plt.close(fig)
                except Exception as e:
                    self.add_log(f"Error saving {name}: {e}")
                    
        # Refresh results tab
        self.results_tab.refresh_plot_list()
        
        # Switch to results if plots created
        if self.latest_plots and self.save_plots.get():
            self.notebook.select(self.results_frame)
            
    def _show_completion_message(self, config, output_dir):
        """Show simulation completion message"""
        message = f"Simulation completed!\n\n"
        message += f"Duration: {config.simulation_time} minutes\n"
        message += f"Satellites: {len(config.spacecraft_list)}\n"
        message += f"Targets: {len(config.targets)}\n"
        message += f"Results: {output_dir}\n"
        
        if self.latest_plots:
            message += f"\n{len(self.latest_plots)} plots saved"
            
        messagebox.showinfo("Simulation Complete", message)
        
    # File operations
    def new_constellation(self):
        """Create a new constellation"""
        if self.satellites and messagebox.askyesno("New Constellation",
                                                   "Clear current configuration?"):
            self.satellites.clear()
            self._initialize_default_satellites()
            
            # Reset targets
            for i, target in enumerate(self.targets):
                target["assigned_to"] = [f"Satellite{i+1}"] if i < 4 else []
                
            # Update UI
            self.constellation_tab.update_satellite_listbox()
            self.update_satellite_dropdowns()
            self.update_target_assignments()
            self.update_status_counts()
            
            self.add_log("Created new constellation")
            
    def import_config(self):
        """Import configuration from file"""
        file_path = filedialog.askopenfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json")],
            title="Open Configuration"
        )
        
        if not file_path:
            return
            
        try:
            with open(file_path, 'r') as f:
                config = json.load(f)
                
            # Apply configuration
            self._apply_configuration(config)
            
            self.update_status(f"Loaded: {os.path.basename(file_path)}")
            messagebox.showinfo("Import Success", "Configuration loaded successfully")
            
        except Exception as e:
            messagebox.showerror("Import Failed", str(e))
            self.update_status("Import failed")
            
    def export_config(self):
        """Export configuration to file"""
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json")],
            title="Save Configuration"
        )
        
        if not file_path:
            return
            
        try:
            config = self._create_configuration()
            
            with open(file_path, 'w') as f:
                json.dump(config, f, indent=4)
                
            self.update_status(f"Saved: {os.path.basename(file_path)}")
            messagebox.showinfo("Export Success", f"Configuration saved to {file_path}")
            
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))
            self.update_status("Export failed")
            
    def _create_configuration(self):
        """Create configuration dictionary"""
        return {
            "version": "1.0",
            "simulation_time": self.simulation_time.get(),
            "show_plots": self.show_plots.get(),
            "save_plots": self.save_plots.get(),
            "save_binary": self.save_binary.get(),
            "binary_filename": self.binary_filename,
            "satellites": self.satellites,
            "targets": self.targets,
            "active_satellite_index": self.current_satellite_index.get()
        }
        
    def _apply_configuration(self, config):
        """Apply loaded configuration"""
        # Simulation parameters
        self.simulation_time.set(config.get("simulation_time", 30.0))
        self.show_plots.set(config.get("show_plots", True))
        self.save_plots.set(config.get("save_plots", True))
        self.save_binary.set(config.get("save_binary", True))
        self.binary_filename = config.get("binary_filename", "spacecraft_viz")
        
        # Satellites
        if "satellites" in config:
            self.satellites = config["satellites"]
            self.update_satellite_dropdowns()
            
        # Targets
        if "targets" in config:
            self.targets = config["targets"]
            self.update_target_assignments()
            
        # Active satellite
        index = config.get("active_satellite_index", 0)
        if 0 <= index < len(self.satellites):
            self.current_satellite_index.set(index)
            self.fault_tab.set_active_satellite(index)
            
        self.update_status_counts()
        
    # UI operations
    def open_results_folder(self):
        """Open the results folder"""
        if not os.path.exists(self.plots_dir):
            messagebox.showinfo("Folder Not Found", "No results yet. Run a simulation first.")
            return
            
        try:
            if sys.platform == 'win32':
                os.startfile(self.plots_dir)
            elif sys.platform == 'darwin':
                subprocess.run(['open', self.plots_dir])
            else:
                subprocess.run(['xdg-open', self.plots_dir])
        except Exception as e:
            messagebox.showerror("Error", f"Could not open folder: {e}")
            
    def show_help(self):
        """Show help dialog"""
        help_window = tk.Toplevel(self.root)
        help_window.title("Spacecraft Fault Simulator Help")
        help_window.geometry("800x600")
        help_window.transient(self.root)
        help_window.grab_set()
        
        # Create notebook
        notebook = ttk.Notebook(help_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Add tabs
        topics = [
             ("Overview", get_general_help("overview")),
            ("Constellation", get_general_help("constellation")),
            ("Faults", get_general_help("fault")),
            ("Targets", get_general_help("target")),
            ("Visualization", get_general_help("visualization")),
            ("Output", get_general_help("output")),
            ("Simulation", get_general_help("simulation"))
        ]
        
        for tab_name, content in topics:
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=tab_name)

            text = tk.Text(frame, wrap=tk.WORD, font=('Segoe UI', 10))
            scrollbar = ttk.Scrollbar(frame, command=text.yview)
            text.config(yscrollcommand=scrollbar.set)

            text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

            text.insert(tk.END, content)
            text.config(state=tk.DISABLED)

            
        # Close button
        ttk.Button(help_window, text="Close", command=help_window.destroy).pack(pady=10)
        
    def show_about(self):
        """Show about dialog"""
        about_window = tk.Toplevel(self.root)
        about_window.title("About")
        about_window.geometry("500x400")
        about_window.transient(self.root)
        about_window.grab_set()
        
        frame = ttk.Frame(about_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        ttk.Label(
            frame,
            text="Spacecraft Constellation Fault Simulator",
            style='Title.TLabel'
        ).pack(pady=(0, 20))
        
        # Description
        description = """A high-fidelity simulation tool for spacecraft reaction wheel
fault analysis in multi-satellite constellations.

Features:
• Multi-orbit constellation modeling
• Target visibility analysis
• Realistic orbital mechanics
• Comprehensive fault injection
• 3D visualization with Vizard
• Advanced data analysis

Built with Python, Tkinter, Matplotlib, and Basilisk

© 2025 Spacecraft Dynamics Lab"""
        
        ttk.Label(frame, text=description, justify=tk.LEFT).pack()
        
        # Close button
        ttk.Button(frame, text="Close", command=about_window.destroy).pack(
            side=tk.BOTTOM, pady=20
        )
        
    # Delegate methods for Results tab
    def refresh_plot_list(self):
        """Delegate to ResultsTab"""
        if hasattr(self, 'results_tab'):
            self.results_tab.refresh_plot_list()
            
    
        
    # Private helper methods
    def _update_sim_time_status(self, *args):
        """Update simulation time in status bar"""
        try:
            sim_time = self.simulation_time.get()
            self.sim_time_status.config(text=f"{sim_time:.0f} min")
        except:
            self.sim_time_status.config(text="30 min")
            
    def _view_results(self):
        """Switch to results tab"""
        self.notebook.select(self.results_frame)


# Main entry point
if __name__ == "__main__":
    root = tk.Tk()
    app = SatelliteSimulatorApp(root)
    root.mainloop()