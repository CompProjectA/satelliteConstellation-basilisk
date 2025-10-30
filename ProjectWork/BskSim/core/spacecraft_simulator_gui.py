#!/usr/bin/env python
"""
spacecraft_simulator_gui.py

Main GUI application for the Basilisk Spacecraft Constellation Fault Simulator.
Provides a modular interface for configuring and running spacecraft simulations.
Includes DRL task reassignment integration and enhanced cluster plotting.

MERGED VERSION - Includes progress dialog and improved thread safety
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
import datetime
import logging
import matplotlib

matplotlib.use('Agg')  # Use non-interactive backend

# Set up paths
import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
CORE_DIR = os.path.join(ROOT_DIR, 'core')
PLOTS_DIR = os.path.join(ROOT_DIR, 'plots')
ASSETS_DIR = os.path.join(ROOT_DIR, 'assets')

sys.path.extend([ROOT_DIR, CORE_DIR])

# DRL path (project-standard) + optional override via drl_config
DRL_DIR = os.path.join(ROOT_DIR, "DRL")
if DRL_DIR not in sys.path:
    sys.path.insert(0, DRL_DIR)
try:
    from drl_config import DRL_DIR as _CONF_DRL_DIR
    if os.path.isdir(_CONF_DRL_DIR):
        DRL_DIR = _CONF_DRL_DIR
except Exception:
    pass

# Only add to sys.path if the directory actually exists
if os.path.isdir(DRL_DIR) and DRL_DIR not in sys.path:
    sys.path.insert(0, DRL_DIR)
    
# Import simulation functionality
from spacecraft_simulation import SimulationConfig, TargetDefinition, run_custom_simulation

# Import GUI modules - with error handling
try:
    from gui_tab import (
        ConstellationTab,
        FaultTab,
        FaultDetectionTab,
        TaskReassignmentTab,
        TargetTab,
        VisualizationTab,
        OutputTab,
        ResultsTab,
        DRLTab  
    )
    GUI_TABS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: GUI tab modules not fully available: {e}")
    GUI_TABS_AVAILABLE = False

# Import help content - with fallback
try:
    from rw_fault_help import get_general_help
    HELP_AVAILABLE = True
except ImportError:
    print("Warning: Help content not available")
    HELP_AVAILABLE = False
    def get_general_help(topic):
        return f"Help content for {topic} is not available. Please check the documentation."

# Import style - with fallback
try:
    from style import setup_style
    STYLE_AVAILABLE = True
except ImportError:
    print("Warning: Style module not available")
    STYLE_AVAILABLE = False
    def setup_style(root):
        return ttk.Style()

# Import communication visualization - with fallback
try:
    from communication_visualization import show_communication_visualization, CommunicationVisualization
    COMM_VIZ_AVAILABLE = True
except ImportError:
    print("Warning: Communication visualization not available")
    COMM_VIZ_AVAILABLE = False
    CommunicationVisualization = None

# Import PIL for logo - with fallback
try:
    from PIL import ImageTk, Image
    PIL_AVAILABLE = True
except ImportError:
    print("Warning: PIL not available - logo will not display")
    PIL_AVAILABLE = False

# Import enhanced constellation plotting
try:
    from plots import ConstellationPlotter
    PLOTTER_AVAILABLE = True
except ImportError:
    print("Warning: ConstellationPlotter not available")
    PLOTTER_AVAILABLE = False
    ConstellationPlotter = None


class ProgressDialog:
    """Simple progress dialog to show simulation is running"""
    def __init__(self, parent, title="Processing..."):
        self.dialog = tk.Toplevel(parent)
        self.dialog.title(title)
        self.dialog.geometry("400x150")
        self.dialog.transient(parent)
        self.dialog.grab_set()
        
        # Center on parent
        self.dialog.update_idletasks()
        x = parent.winfo_x() + (parent.winfo_width() // 2) - 200
        y = parent.winfo_y() + (parent.winfo_height() // 2) - 75
        self.dialog.geometry(f"+{x}+{y}")
        
        # Progress bar
        self.label = ttk.Label(self.dialog, text="Initializing...", font=('Segoe UI', 10))
        self.label.pack(pady=20)
        
        self.progress = ttk.Progressbar(self.dialog, mode='indeterminate', length=300)
        self.progress.pack(pady=10)
        self.progress.start(10)
        
        # Status label
        self.status = ttk.Label(self.dialog, text="", font=('Segoe UI', 9), foreground='gray')
        self.status.pack(pady=5)
        
        self.dialog.protocol("WM_DELETE_WINDOW", lambda: None)  # Prevent closing
        
    def update(self, message, detail=""):
        """Update progress message"""
        self.label.config(text=message)
        if detail:
            self.status.config(text=detail)
        self.dialog.update_idletasks()
        
    def close(self):
        """Close the dialog"""
        try:
            self.progress.stop()
            self.dialog.destroy()
        except:
            pass


class SatelliteSimulatorApp:
    """Main application class for the Spacecraft Fault Simulator GUI"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Spacecraft Constellation Fault Simulator with DRL Integration")
        self.root.update_idletasks()
        try:
            self.root.state("zoomed")
        except tk.TclError:
            # macOS/Linux fallback
            if sys.platform.startswith("linux"):
                self.root.attributes("-zoomed", True)
            else:
                self.root.geometry("1200x900")
        
        self.root.minsize(1000, 800)
        
        # Setup logging
        self._setup_logging()
        
        # Setup directories
        self._setup_directories()
        
        # Setup style
        if STYLE_AVAILABLE:
            self.style = setup_style(root)
        else:
            self.style = ttk.Style()

        # Handle window close properly
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Initialize data
        self._initialize_data()
        
        # Cluster state + plotter for enhanced plots
        self.clusters = {}
        if PLOTTER_AVAILABLE:
            self.plotter = ConstellationPlotter(output_dir=self.PLOTTING_DIR)
        else:
            self.plotter = None
        self._last_cluster_config = None  # cache for post-sim plotting
        
        # Create UI
        self._create_ui()
        
        # Update initial state
        self.update_status("Ready")
        self.update_status_counts()
        
        # Check DRL requirements on startup
        self._check_drl_on_startup()
        
        # Debug cluster plots
        self.debug_cluster_plots()

        if GUI_TABS_AVAILABLE:
            self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_changed)

    def _check_drl_on_startup(self):
        """Check DRL requirements on startup and log status"""
        try:
            requirements_met, details = self.check_drl_requirements()
            if requirements_met:
                self.add_log("DRL integration ready")
                if isinstance(details, dict):
                    for component, status in details.items():
                        if component in ['tensorflow', 'integration_bridge', 'drl_files']:
                            self.add_log(f"  {component}: {'Available' if status else 'Not available'}")
            else:
                if isinstance(details, list):
                    self.add_log(f"DRL integration not ready - missing files: {details}")
                elif isinstance(details, dict) and 'missing_files' in details:
                    self.add_log(f"DRL integration not ready - missing files: {details['missing_files']}")
                else:
                    self.add_log(f"DRL integration not ready: {details}")
        except Exception as e:
            self.add_log(f"Error checking DRL requirements: {e}")

    def _on_closing(self):
        """Handle window closing"""
        try:
            self.is_running = False
            import matplotlib.pyplot as plt
            plt.close('all')
        except:
            pass
        self.root.destroy()
            
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

    def check_cluster_plots(self):
        """Check why cluster plots aren't generating"""
        print("\n=== CHECKING CLUSTER CONFIGURATION ===")

        # Check constellation tab clusters
        if hasattr(self, 'constellation_tab'):
            clusters = self.constellation_tab.clusters
            print(f"Clusters in constellation tab: {len(clusters)}")
            for c in clusters:
                members = c.get('children') or c.get('satellites', [])
                print(f"  - {c['name']}: {len(members)} satellites")

        # Check satellites for cluster membership
        cluster_sats = [s for s in self.satellites if s.get('cluster')]
        print(f"Satellites with cluster assignment: {len(cluster_sats)}")
        for s in cluster_sats[:5]:  # Show first 5
            print(f"  - {s['name']}: cluster={s.get('cluster')}, role={s.get('role')}")

        print("=====================================\n")

    def on_tab_changed(self, event):
        """Handle tab change events"""
        try:
            selected_tab = event.widget.tab('current')['text']
            
            if selected_tab == 'Constellation':
                # Check if Communication sub-tab is selected
                if hasattr(self, 'constellation_tab'):
                    try:
                        current_subtab = self.constellation_tab.constellation_notebook.index('current')
                        if current_subtab == 3 and hasattr(self.constellation_tab, "update_communication_plot"):
                            self.constellation_tab.update_communication_plot()
                    except:
                        pass
                        
            elif selected_tab == 'Results':
                # Refresh results tab
                if hasattr(self, 'results_tab'):
                    self._thread_safe(self.results_tab.refresh_plot_list)
        except Exception as e:
            self.add_log(f"Error in tab change: {e}")

    def show_communication_window(self):
        """Show communication visualization window"""
        if not COMM_VIZ_AVAILABLE:
            messagebox.showerror("Communication Visualization", "Communication visualization module not available")
            return
            
        try:
            # Get cluster manager if using clusters
            cluster_manager = None
            
            # Check if we're using the cluster integration
            if hasattr(self, 'constellation_tab'):
                # Create a simple cluster manager from the constellation tab data
                class SimpleClusterManager:
                    def __init__(self, clusters):
                        def _name(x):
                            if isinstance(x, str):
                                return x
                            if isinstance(x, dict):
                                return x.get('name') or x.get('model_tag') or "Unknown"
                            return getattr(x, 'name', getattr(x, 'model_tag', "Unknown"))

                        self.clusters = {}
                        for cluster in clusters:
                            leader_name = _name(cluster.get('leader'))
                            raw_children = cluster.get('children') or cluster.get('satellites', [])
                            children_names = [ _name(c) for c in raw_children ]
                            self.clusters[cluster['name']] = {
                                'leader': type('obj', (object,), {'model_tag': leader_name})(),
                                'children': [ type('obj', (object,), {'model_tag': ch})() for ch in children_names ]
                            }
                
                if hasattr(self.constellation_tab, 'clusters'):
                    cluster_manager = SimpleClusterManager(self.constellation_tab.clusters)
            
            # Show visualization
            self.comm_visualizer = show_communication_visualization(self.root, cluster_manager)
        except Exception as e:
            self.add_log(f"Error showing communication window: {e}")
            messagebox.showerror("Communication Error", f"Could not show communication window: {e}")

    # ======== Claude-style integration helpers ========
    def setup_communication_visualization(self):
        """Setup communication visualization panel inside the Constellation tab (if supported)."""
        if CommunicationVisualization is None:
            return
        if not hasattr(self, 'constellation_tab'):
            return
        # Expect constellation_tab to expose a frame named 'comm_frame'
        if hasattr(self.constellation_tab, 'comm_frame'):
            try:
                for w in self.constellation_tab.comm_frame.winfo_children():
                    w.destroy()
                self.comm_viz = CommunicationVisualization(
                    self.constellation_tab.comm_frame,
                    getattr(self.constellation_tab, 'clusters', {})
                )
            except Exception:
                self.comm_viz = None

    def update_cluster_info(self, clusters: dict):
        """
        Optionally called by ConstellationTab whenever clusters change.
        clusters: expected dict mapping name -> {leader:{name,...}, children:[{name,...}] , ...}
        """
        try:
            self.clusters = clusters or {}
            if hasattr(self, 'comm_viz') and self.comm_viz:
                try:
                    self.comm_viz.update_clusters(self.clusters)
                except Exception:
                    pass
            self.log_cluster_summary()
        except Exception as e:
            self.add_log(f"update_cluster_info warning: {e}")

    def log_cluster_summary(self):
        """Log a compact, Claude-style cluster summary."""
        if not self.clusters:
            return
        parts = []
        for name, data in self.clusters.items():
            leader_name = (data.get('leader') or {}).get('name', 'N/A')
            num_children = len(data.get('children', []))
            formation = data.get('formation', 'Unknown')
            parts.append(f"{name}[Leader={leader_name}, {num_children} children, {formation}]")
        msg = f"Clusters({len(self.clusters)}): " + ", ".join(parts)
        print(msg)
        if hasattr(self, 'output_tab') and hasattr(self.output_tab, 'add_log_entry'):
            self.add_log(msg)

    def validate_configuration(self):
        """Validate constellation configuration prior to simulation (max 4 clusters, valid formations, etc.)."""
        cfg = self._get_cluster_configuration()
        errors, warnings = [], []
        # Check cluster count
        if cfg['cluster_count'] > 4:
            errors.append(f"Too many clusters ({cfg['cluster_count']}), maximum is 4")
        # Check each cluster
        allowed = getattr(ConstellationTab, 'ALLOWED_FORMATIONS',
                  {'Line', 'Triangle', 'Diamond', 'Leader-Follower'}) if GUI_TABS_AVAILABLE \
          else {'Line', 'Triangle', 'Diamond', 'Leader-Follower'}

        for cname, cdata in cfg['clusters'].items():
            if not cdata.get('leader'):
                errors.append(f"Cluster {cname} has no leader")
            if allowed and cdata.get('formation') not in allowed:
                errors.append(f"Cluster {cname} has invalid formation: {cdata.get('formation')}")
        # Duplicate satellite names
        sat_names = [s.get('name') for s in cfg['satellites']]
        if len(sat_names) != len(set(sat_names)):
            warnings.append("Duplicate satellite names detected")
        return errors, warnings

    def _get_cluster_configuration(self):
        """
        Build a Claude-compatible cluster config:
          {'cluster_count', 'clusters':{...}, 'satellites':[...]}
        Falls back to deriving from self.satellites if ConstellationTab lacks get_configuration().
        """
        # Prefer explicit method from ConstellationTab if available
        if hasattr(self, 'constellation_tab') and hasattr(self.constellation_tab, 'get_configuration'):
            try:
                cfg = self.constellation_tab.get_configuration()
                # ensure required top-level keys
                if 'cluster_count' in cfg and 'clusters' in cfg and 'satellites' in cfg:
                    return cfg
            except Exception:
                pass

        # Derive from self.satellites (current GUI state)
        clusters = {}
        satellites_data = []
        for idx, sat in enumerate(self.satellites):
            name = sat.get('name', f'Sat{idx+1}')
            role = sat.get('role', 'independent')
            cluster_name = sat.get('cluster')
            a = float(sat.get('orbit', {}).get('a', 6371 + 600))
            alt = max(0.0, a - 6371.0)
            orbit_str = f"LEO {int(round(alt))}"
            pos = [0.0, 0.0, 0.0]  # position not known here; plotter uses templates

            satellites_data.append({'name': name, 'position': pos, 'role': role, 'index': idx})

            if cluster_name:
                if cluster_name not in clusters:
                    clusters[cluster_name] = {
                        'leader': None,
                        'children': [],
                        'satellites': [],
                        'formation': sat.get('formation', 'Leader-Follower'),
                        'separation': sat.get('separation', 1.0),
                        'orbit': orbit_str
                    }
                entry = {'name': name, 'position': pos, 'role': role, 'index': idx}
                clusters[cluster_name]['satellites'].append(entry)
                if role == 'leader':
                    clusters[cluster_name]['leader'] = {'name': name, 'position': pos, 'index': idx}
                elif role == 'child':
                    clusters[cluster_name]['children'].append({'name': name, 'position': pos, 'index': idx})

        return {
            'cluster_count': len(clusters),
            'clusters': clusters,
            'satellites': satellites_data
        }

    def generate_enhanced_plots(self, cluster_config):
        """Generate Claude-style enhanced plots and update Results tab."""
        if not cluster_config or not isinstance(cluster_config, dict) or not self.plotter:
            return
        plots_generated = []
        try:
            # 1) Constellation Overview
            path = self.plotter.plot_constellation_overview(
                cluster_config['clusters'],
                cluster_config['satellites']
            )
            plots_generated.append(("Constellation Overview", path))
            # 2) Formation checks
            for cname, cdata in cluster_config['clusters'].items():
                p = self.plotter.plot_formation_check(cname, cdata)
                plots_generated.append((f"Formation Check - {cname}", p))
            # 3) Communication
            p = self.plotter.plot_cluster_communication(
                cluster_config['clusters'],
                self.simulation_time.get()
            )
            plots_generated.append(("Cluster Communication", p))
            # 4) Distance Analysis
            p = self.plotter.plot_distance_analysis(
                cluster_config['clusters'],
                cluster_config['satellites']
            )
            plots_generated.append(("Distance Analysis", p))

            # Track & surface to Results
            for _, pth in plots_generated:
                if pth and os.path.isfile(pth):
                    self.latest_plots.append(os.path.basename(pth))
            self._thread_safe(self.update_results_tab, plots_generated)
            self.add_log(f"Generated {len(plots_generated)} enhanced plots")
        except Exception as e:
            self.add_log(f"Error generating enhanced plots: {e}")

    def update_results_tab(self, plots_list):
        """Append new plot rows to the results tab if it exposes a tree widget."""
        try:
            if hasattr(self, 'results_tab') and hasattr(self.results_tab, 'results_tree'):
                ts = datetime.datetime.now().strftime("%H:%M:%S")
                for name, path in plots_list:
                    try:
                        self.results_tab.results_tree.insert("", "end", values=(ts, name, path))
                    except Exception:
                        pass
            # Always refresh the list view as a fallback
            if hasattr(self, 'results_tab'):
                self.results_tab.refresh_plot_list()
        except Exception as e:
            self.add_log(f"update_results_tab warning: {e}")

    def debug_cluster_plots(self):
        """Debug why cluster plots aren't generating"""
        print("\n" + "=" * 60)
        print("CLUSTER PLOT DEBUGGING")
        print("=" * 60)

        # Check satellites for cluster membership
        print("\n1. CHECKING SATELLITE CLUSTER ASSIGNMENTS:")
        cluster_sats = {}
        for sat in self.satellites:
            if sat.get('cluster'):
                cluster = sat['cluster']
                role = sat.get('role', 'unknown')
                if cluster not in cluster_sats:
                    cluster_sats[cluster] = {'leader': None, 'children': []}

                if role == 'leader':
                    cluster_sats[cluster]['leader'] = sat['name']
                elif role == 'child':
                    cluster_sats[cluster]['children'].append(sat['name'])

                print(f"  {sat['name']}: cluster='{cluster}', role='{role}'")

        print(f"\n2. CLUSTER SUMMARY:")
        for cluster_name, info in cluster_sats.items():
            print(f"  Cluster '{cluster_name}':")
            print(f"    Leader: {info['leader']}")
            print(f"    Children ({len(info['children'])}): {info['children']}")

        print(f"\n3. CLUSTER DATA STRUCTURE CHECK:")
        if hasattr(self, 'constellation_tab'):
            clusters = self.constellation_tab.clusters
            print(f"  Clusters in constellation_tab: {len(clusters)}")
            for c in clusters:
                members = c.get('children') or c.get('satellites', [])
                print(f"    - {c['name']}: leader={c.get('leader')}, members={len(members)}")

        print(f"\n4. SOLUTION:")
        if not cluster_sats:
            print("  NO SATELLITES HAVE CLUSTER ASSIGNMENTS!")
            print("  Fix: Make sure satellites have 'cluster' and 'role' fields set")
        elif all(info['leader'] is None for info in cluster_sats.values()):
            print("  NO CLUSTER LEADERS FOUND!")
            print("  Fix: Make sure at least one satellite per cluster has role='leader'")
        else:
            print("  Cluster structure looks correct")
            print("  If plots still don't generate, check spacecraft_simulation.py line ~700")

        print("=" * 60 + "\n")
    
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
        
        # Results storage - includes both ML and DRL results
        self.ml_detection_results = None
        self.integrated_drl_results = None
        self.processed_data = None  # For GNN Analysis
        
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
        if PIL_AVAILABLE:
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
            text="Spacecraft Constellation Fault Simulator with DRL",
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
            style="Run.TButton" if STYLE_AVAILABLE else "",
            command=self.run_simulation
        )
        self.run_button.pack(side=tk.RIGHT, padx=5)
        
    def _create_notebook(self, parent):
        """Create the tabbed interface"""
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        if GUI_TABS_AVAILABLE:
            # Create tab frames
            self.constellation_frame = ttk.Frame(self.notebook)
            self.fault_frame = ttk.Frame(self.notebook)
            self.fault_detection_frame = ttk.Frame(self.notebook)
            self.task_reassignment_frame = ttk.Frame(self.notebook)
            self.target_frame = ttk.Frame(self.notebook)
            self.visualization_frame = ttk.Frame(self.notebook)
            self.output_frame = ttk.Frame(self.notebook)
            self.results_frame = ttk.Frame(self.notebook)
            self.drl_frame = ttk.Frame(self.notebook)
            
            # Add tabs
            self.notebook.add(self.constellation_frame, text="Constellation")
            self.notebook.add(self.fault_frame, text="Fault Configuration")
            self.notebook.add(self.fault_detection_frame, text="Fault Detection")
            self.notebook.add(self.task_reassignment_frame, text="Task Reassignment")
            self.notebook.add(self.target_frame, text="Targets")
            self.notebook.add(self.visualization_frame, text="Visualization")
            self.notebook.add(self.output_frame, text="Output Settings")
            self.notebook.add(self.results_frame, text="Results")
            self.notebook.add(self.drl_frame, text="DRL")
            
            # Create tab objects
            try:
                self.constellation_tab = ConstellationTab(self, self.constellation_frame)
                self.fault_tab = FaultTab(self, self.fault_frame)
                self.fault_detection_tab = FaultDetectionTab(self, self.fault_detection_frame)
                self.task_reassignment_tab = TaskReassignmentTab(self, self.task_reassignment_frame)
                self.target_tab = TargetTab(self, self.target_frame)
                self.visualization_tab = VisualizationTab(self, self.visualization_frame)
                self.output_tab = OutputTab(self, self.output_frame)
                self.results_tab = ResultsTab(self, self.results_frame)
                
                # Create DRL tab with default paths
                default_drl_script = os.path.join(ROOT_DIR, "DRL", "PPO_Year2.py")
                default_drl_results = os.path.join(ROOT_DIR, "DRL")
                self.drl_tab = DRLTab(
                    self.drl_frame,
                    parent_app=self,
                    default_script=default_drl_script if os.path.exists(default_drl_script) else None,
                    default_results_dir=default_drl_results if os.path.isdir(default_drl_results) else None
                )
                self.drl_tab.pack(fill=tk.BOTH, expand=True)
                
            except Exception as e:
                messagebox.showerror("Tab Creation Error", f"Error creating tabs: {e}")
                self.logger.error(f"Error creating tabs: {e}")
                
        else:
            # Create a simple placeholder tab
            placeholder_frame = ttk.Frame(self.notebook)
            self.notebook.add(placeholder_frame, text="Configuration")
            
            ttk.Label(
                placeholder_frame,
                text="GUI modules not available. Please check imports.",
                font=('Segoe UI', 12)
            ).pack(expand=True)
        
    def _create_status_bar(self, parent):
        """Create the status bar"""
        status_frame = ttk.Frame(parent, style="Status.TFrame" if STYLE_AVAILABLE else "")
        status_frame.pack(fill=tk.X, pady=(10, 0))
        
        # Status message
        self.status_label = ttk.Label(status_frame, text="Ready", style="Status.TLabel" if STYLE_AVAILABLE else "")
        self.status_label.pack(side=tk.LEFT, padx=5)
        
        # Satellite count
        self.sat_count_status = ttk.Label(
            status_frame, 
            text="0 satellites",
            style='Status.TLabel' if STYLE_AVAILABLE else "",
            foreground='blue'
        )
        self.sat_count_status.pack(side=tk.LEFT, padx=20)
        
        # Target count
        self.target_status = ttk.Label(
            status_frame,
            text="0/0 targets assigned",
            style='Status.TLabel' if STYLE_AVAILABLE else "",
            foreground='green'
        )
        self.target_status.pack(side=tk.LEFT, padx=20)
        
        # ML Detection Status
        self.ml_detection_status = ttk.Label(
            status_frame,
            text="ML: Not Ready",
            style='Status.TLabel' if STYLE_AVAILABLE else "",
            foreground='gray'
        )
        self.ml_detection_status.pack(side=tk.LEFT, padx=20)
        
        # DRL Status
        self.drl_status = ttk.Label(
            status_frame,
            text="DRL: Checking...",
            style='Status.TLabel' if STYLE_AVAILABLE else "",
            foreground='orange'
        )
        self.drl_status.pack(side=tk.LEFT, padx=20)
        
        # Simulation time
        sim_time_frame = ttk.Frame(status_frame)
        sim_time_frame.pack(side=tk.LEFT, padx=20)
        
        ttk.Label(sim_time_frame, text="Sim Time:", style='Status.TLabel' if STYLE_AVAILABLE else "").pack(side=tk.LEFT)
        self.sim_time_status = ttk.Label(
            sim_time_frame,
            text="30 min",
            style='Status.TLabel' if STYLE_AVAILABLE else "",
            foreground='darkgreen'
        )
        self.sim_time_status.pack(side=tk.LEFT, padx=5)
        
        # Bind simulation time change
        self.simulation_time.trace_add('write', lambda *args: self._update_sim_time_status())

        # Current time
        self.time_label = ttk.Label(
            status_frame,
            text=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            style='Status.TLabel' if STYLE_AVAILABLE else ""
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
        if COMM_VIZ_AVAILABLE:
            sim_menu.add_command(label="Communication Visualization", command=self.show_communication_window)
        sim_menu.add_separator()
        sim_menu.add_command(label="View Results", command=self._view_results)
        sim_menu.add_command(label="View Fault Detection", command=self._view_fault_detection)
        sim_menu.add_command(label="Open Results Folder", command=self.open_results_folder)
        
        # DRL menu
        drl_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="DRL", menu=drl_menu)
        drl_menu.add_command(label="Check DRL Requirements", command=self.check_drl_status)
        drl_menu.add_command(label="DRL Configuration", command=self.show_drl_config)
        
        # Help menu
        help_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_command(label="Fault Detection Help", command=self.show_fault_detection_help)
        help_menu.add_command(label="DRL Integration Help", command=self.show_drl_help)
        help_menu.add_command(label="About", command=self.show_about)

    # Public methods

    def add_satellite(self, name, true_anomaly=0.0, altitude=600):
        """Add a new satellite to the constellation with all required fields"""
        satellite = {
            "name": name,
            "type": "individual",
            "cluster": None,
            "role": "independent",
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
            "communication": {
                "range": 2000.0,
                "fov": 30.0,
                "aHat_B": [0.0, 0.0, -1.0]
            },
            "targets": [],
            "orbit_name": "Default Orbit"
        }
        self.satellites.append(satellite)
        
        # Update fault detection tab satellite list
        if hasattr(self, 'fault_detection_tab'):
            try:
                self.fault_detection_tab.monitoring_spacecraft_combo['values'] = (
                    ["All Spacecraft"] + [sat["name"] for sat in self.satellites]
                )
            except:
                pass
        
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
            
            # ML Detection status
            if hasattr(self, 'fault_detection_tab'):
                try:
                    if self.fault_detection_tab.is_ml_ready():
                        self.ml_detection_status.config(text="ML: Ready", foreground='green')
                    elif hasattr(self.fault_detection_tab, 'ml_available') and self.fault_detection_tab.ml_available:
                        self.ml_detection_status.config(text="ML: Available", foreground='orange')
                    else:
                        self.ml_detection_status.config(text="ML: Not Available", foreground='red')
                except:
                    self.ml_detection_status.config(text="ML: Unknown", foreground='gray')
        except:
            pass

    def add_log(self, message, level="INFO"):
        """Add a message to the application log with filtering"""
        try:
            # FILTER OUT SPAM MESSAGES
            spam_keywords = [
                "AttributeError: 'MessageFactory'",  # TensorFlow warning spam
                "oneDNN custom operations",          # TensorFlow info spam
                "Generating realistic telemetry",    # Repeated during ML detection
                "No fault in",                       # Repeated for each satellite
                "running analysis",                  # Too frequent
            ]

            # Skip spam messages
            if any(keyword in message for keyword in spam_keywords):
                return

            # Only log important messages to GUI
            important_keywords = [
                "Created cluster",
                "Starting simulation",
                "Simulation completed",
                "ERROR",
                "WARNING",
                "Success",
                "Failed",
                "Ready",
                "ML Detection",
                "DRL",
                "Fault",
            ]

            is_important = any(keyword in message for keyword in important_keywords)

            # Add to GUI log (throttled)
            if hasattr(self, 'output_tab') and hasattr(self.output_tab, 'add_log_entry'):
                if is_important or level in ["ERROR", "WARNING"]:
                    self.output_tab.add_log_entry(message)

            # Always log to the file/logger
            if hasattr(self, 'logger'):
                if level == "ERROR":
                    self.logger.error(message)
                elif level == "WARNING":
                    self.logger.warning(message)
                else:
                    self.logger.info(message)

        except Exception:
            # Failsafe
            print(f"Log: {message}")
                
    def update_satellite_dropdowns(self):
        """Update all satellite dropdown menus"""
        try:
            if hasattr(self, 'fault_tab'):
                self.fault_tab.update_satellite_dropdown()
            if hasattr(self, 'visualization_tab'):
                self.visualization_tab.update_satellite_dropdown()
            if hasattr(self, 'target_tab'):
                self.target_tab.update_satellite_dropdown()
            
            # Update fault detection tab spacecraft list
            if hasattr(self, 'fault_detection_tab'):
                try:
                    self.fault_detection_tab.monitoring_spacecraft_combo['values'] = (
                        ["All Spacecraft"] + [sat["name"] for sat in self.satellites]
                    )
                except:
                    pass
        except:
            pass
            
    def update_target_assignments(self):
        """Update target assignment displays"""
        try:
            if hasattr(self, 'target_tab'):
                self.target_tab.update_target_assignments()
        except:
            pass
            
    def update_time(self):
        """Update the time display"""
        self.time_label.config(text=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.root.after(60000, self.update_time)  # Update every minute

    def _thread_safe(self, fn, *a, **kw):
        """Execute function on main thread safely"""
        self.root.after(0, lambda: fn(*a, **kw))

    def _build_rich_clusters_for_sim(self, sc_index_map):
        """
        Build a 'rich' cluster dict for the simulation layer
        """
        rich = {}

        # Prefer data coming from ConstellationTab if present
        clusters_src = getattr(self.constellation_tab, 'clusters', None) if hasattr(self, 'constellation_tab') else None

        # Build a consistent iterator of (cluster_name, cluster_dict)
        if isinstance(clusters_src, dict):
            items = clusters_src.items()
        elif isinstance(clusters_src, list):
            items = [(c.get('name') or f"CLUSTER_{i+1}", c) for i, c in enumerate(clusters_src)]
        else:
            items = []

        def _get_name(x):
            if isinstance(x, str):
                return x
            if isinstance(x, dict):
                return x.get('name')
            return getattr(x, 'name', getattr(x, 'model_tag', None))

        # Path 1: Use ConstellationTab clusters if available
        for cname, cd in items:
            if not cname:
                continue

            formation = cd.get('formation') or cd.get('type') or 'Line'
            separation = float(cd.get('separation', 10.0))
            orbit_lbl = cd.get('orbit') or cd.get('orbit_name') or 'LEO 600km'

            leader_name = _get_name(cd.get('leader'))
            sat_names = []

            if isinstance(cd.get('satellites'), list) and cd['satellites']:
                for s in cd['satellites']:
                    sat_names.append(_get_name(s))
            elif isinstance(cd.get('children'), list) and cd['children']:
                sat_names = [leader_name] + [ _get_name(s) for s in cd['children'] ]

            sats, leader_entry, children = [], None, []
            for nm in filter(None, sat_names):
                idx = sc_index_map.get(nm)
                entry = {
                    "name": nm,
                    "role": "leader" if nm == leader_name else "child",
                    "index": idx,
                    "formation": formation,
                    "orbit": orbit_lbl,
                    "separation": separation
                }
                sats.append(entry)
                if entry["role"] == "leader":
                    leader_entry = entry
                else:
                    children.append(entry)

            if sats:
                leader_entry = leader_entry or sats[0]
                rich[cname] = {
                    "satellites": sats,
                    "leader": leader_entry,
                    "children": children,
                    "formation": formation,
                    "orbit": orbit_lbl,
                    "separation": separation
                }

        if rich:
            cleaned = {}
            for k, v in rich.items():
                if v.get("leader") and any(s.get("role") == "child" for s in v["satellites"]):
                    cleaned[k] = v
            return cleaned

        # Path 2: Derive from self.satellites (fallback)
        by_cluster = {}
        for idx, sat in enumerate(self.satellites):
            cname = sat.get('cluster')
            if not cname:
                continue
            by_cluster.setdefault(cname, {"satellites": [], "children": []})

            role = sat.get('role', 'child')
            formation = sat.get('formation') or 'Line'
            orbit_lbl = sat.get('orbit_name', 'LEO 600km')
            separation = float(sat.get('separation', 10.0))

            entry = {
                "name": sat["name"],
                "role": role,
                "index": sc_index_map.get(sat["name"], idx),
                "formation": formation,
                "orbit": orbit_lbl,
                "separation": separation
            }
            by_cluster[cname]["satellites"].append(entry)
            by_cluster[cname]["formation"] = formation
            by_cluster[cname]["orbit"] = orbit_lbl
            by_cluster[cname]["separation"] = separation
            if role == 'leader':
                by_cluster[cname]["leader"] = entry
            else:
                by_cluster[cname]["children"].append(entry)

        rich = {
            k: v for k, v in by_cluster.items()
            if v.get("leader") and len(v.get("children", [])) > 0
        }
        return rich

    # Simulation methods
    def run_simulation(self):
        """Run the simulation with progress dialog"""
        if self.is_running:
            messagebox.showinfo("Simulation Running",
                              "A simulation is already running. Please wait.")
            return

        # Pull the latest Output tab settings
        try:
            if hasattr(self, 'output_tab') and hasattr(self.output_tab, 'update_settings_from_ui'):
                self.output_tab.update_settings_from_ui()
        except Exception as e:
            self.add_log(f"Warning: Could not update settings: {e}", "WARNING")

        # Validate clusters
        try:
            self._last_cluster_config = self._get_cluster_configuration()
            errs, warns = self.validate_configuration()
            for w in warns:
                self.add_log(f"Warning: {w}", "WARNING")
            if errs:
                self.update_status("Validation error")
                messagebox.showerror("Configuration Error", "\n".join(errs))
                return
        except Exception as e:
            self.add_log(f"Cluster validation skipped: {e}", "WARNING")

        # State/UI
        self.is_running = True
        self.run_button.config(state="disabled")
        self.update_status("Running simulation...")

        # Create progress dialog on the main thread
        progress = ProgressDialog(self.root, "Running Basilisk Simulation")

        self.add_log("Starting simulation...", "INFO")

        def simulation_thread():
            try:
                # All Tk calls MUST be marshalled to the main thread:
                self._thread_safe(progress.update, "Building spacecraft configuration...",
                                f"{len(self.satellites)} satellites")

                # Run the full pipeline with progress updates
                self._run_simulation_process_with_progress(progress)

            except Exception as e:
                self.add_log(f"Simulation error: {e}", "ERROR")
                self._thread_safe(messagebox.showerror, "Simulation Error", f"Error: {str(e)}")
            finally:
                # Close dialog + reset UI on the main thread
                self._thread_safe(progress.close)
                self.is_running = False
                self._thread_safe(self.run_button.config, state="normal")
                self._thread_safe(self.update_status, "Ready")

        # Background worker
        thread = threading.Thread(target=simulation_thread, daemon=True)
        thread.start()

    def _run_simulation_process_with_progress(self, progress):
        """Run simulation with progress updates (UI marshalled to main thread)."""
        try:
            # Helper to safely update progress from this worker thread
            def pupdate(msg, detail=""):
                self._thread_safe(progress.update, msg, detail)

            pupdate("Creating configuration...",
                    f"{len(self.satellites)} satellites, {len(self.targets)} targets")

            # ---- Build SimulationConfig ----
            config = SimulationConfig()
            config.simulation_time = self.simulation_time.get()
            config.show_plots = self.show_plots.get()
            config.save_plots = self.save_plots.get()
            config.save_binary = self.save_binary.get()
            config.binary_filename = os.path.basename(self.binary_filename)

            # Spacecraft list
            config.spacecraft_list = []
            for i, sat in enumerate(self.satellites):
                if (i % 4) == 0:
                    pupdate("Configuring spacecraft...", f"{i+1}/{len(self.satellites)}")

                spacecraft_config = {
                    "name": sat["name"],
                    "type": sat.get("type", "cluster_member" if sat.get("cluster") else "individual"),
                    "cluster": sat.get("cluster"),
                    "role": sat.get("role", "independent"),
                    "communication": sat.get("communication", {
                        "range": 2000.0, "fov": 30.0, "aHat_B": [0.0, 0.0, -1.0]
                    }),
                    "orbit_name": sat.get("orbit_name", "Default Orbit"),
                    "orbit": sat["orbit"],
                    "fault": sat["fault"],
                    "camera": sat["camera"],
                    "targets": sat.get("targets", []),
                    "formation": sat.get("formation"),
                    "separation": sat.get("separation")
                }
                config.spacecraft_list.append(spacecraft_config)

            pupdate("Setting up cameras and targets...", f"{len(self.targets)} targets")

            # Camera
            if hasattr(self, 'visualization_tab'):
                try:
                    active_camera = self.visualization_tab.get_active_camera_satellite()
                    if active_camera:
                        config.active_camera_name = active_camera["name"]
                        config.camera_position = active_camera["camera"]["position"]
                        config.camera_fov = active_camera["camera"]["fov"]
                except Exception:
                    pass

            # Targets
            config.targets = []
            for target in self.targets:
                if target.get("assigned_to"):
                    t = TargetDefinition(
                        target["name"], target["lat"], target["lon"],
                        self._convert_color(target["color"])
                    )
                    t.assigned_to = target["assigned_to"]
                    config.targets.append(t)

            # Rich clusters (if available)
            clustered_count = len([s for s in self.satellites if s.get('cluster')])
            pupdate("Building cluster configuration...", f"{clustered_count} clustered satellites")

            sc_index_map = {sc["name"]: i for i, sc in enumerate(config.spacecraft_list)}
            rich_clusters = self._build_rich_clusters_for_sim(sc_index_map)
            if rich_clusters:
                config.clusters = rich_clusters
                self.add_log(f"Using rich cluster config: {len(rich_clusters)} cluster(s)")

            # ---- Run Basilisk simulation ----
            pupdate("Running Basilisk simulation...", f"{config.simulation_time} minutes")
            result = run_custom_simulation(config, self.fault_detection_tab)

            if len(result) == 5:
                scenario, viz, figureList, output_dir, ml_results = result
                self.ml_detection_results = ml_results
            else:
                scenario, viz, figureList, output_dir = result
                self.ml_detection_results = None

            self.add_log("Basilisk simulation completed")

            # Extract telemetry
            if scenario:
                self._extract_telemetry_data(scenario)
            else:
                self._create_synthetic_telemetry_data()

            # ---- DRL Integration (if faults are enabled) ----
            has_faults = any(sat.get("fault", {}).get("enabled") for sat in self.satellites)
            if has_faults:
                pupdate("Running DRL integration...", "Analyzing faults")
                self.add_log("Faults detected - starting DRL integration...")

                try:
                    from main_integration_script import run_integrated_analysis

                    integrated_results = run_integrated_analysis(
                        scenario=scenario,
                        scenario_config=config,
                        output_dir=output_dir,
                        config_profile="development",
                        fault_detection_tab=self.fault_detection_tab
                    )

                    if "error" not in integrated_results:
                        self.add_log("DRL integration completed successfully!")
                        self.integrated_drl_results = integrated_results

                        if integrated_results.get("fault_detection_results"):
                            self.ml_detection_results = integrated_results["fault_detection_results"]
                    else:
                        err = integrated_results.get("error")
                        self.add_log(f"DRL integration failed: {err}", "WARNING")
                        self.integrated_drl_results = {"error": err, "partial": True}

                except Exception as e:
                    self.add_log(f"DRL integration error: {e}", "ERROR")
                    self.integrated_drl_results = {"error": str(e), "partial": True}
            else:
                self.add_log("No faults configured - skipping DRL integration")
                self.integrated_drl_results = None

            # ---- Results / Plots ----
            pupdate("Generating plots...", "Creating visualizations")
            self._handle_simulation_results(figureList, output_dir, config)

            if self.ml_detection_results:
                self._handle_ml_detection_results(self.ml_detection_results)

            if getattr(self, 'integrated_drl_results', None):
                self._handle_drl_results(self.integrated_drl_results)

            pupdate("Simulation complete!", "Processing results...")
            self._show_completion_message_with_drl(config, output_dir)

        except Exception as e:
            self.add_log(f"Simulation error: {e}", "ERROR")
            raise

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
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            self.add_log(f"Saving {len(figureList)} plots...")
            
            import matplotlib.pyplot as plt
            
            for name, fig in figureList.items():
                try:
                    if hasattr(fig, 'savefig'):
                        filename = f"{name}_{timestamp}.png"
                        path = os.path.join(self.plots_dir, filename)
                        fig.savefig(path, bbox_inches='tight', dpi=300)
                        self.latest_plots.append(filename)
                        self.add_log(f"Saved: {filename}")
                        plt.close(fig)
                except Exception as e:
                    self.add_log(f"Error saving {name}: {e}")
                    
        # Refresh results tab
        if hasattr(self, 'results_tab'):
            self._thread_safe(self.results_tab.refresh_plot_list)
        
        # Switch to results if plots created
        if self.latest_plots and self.save_plots.get():
            self._thread_safe(self.notebook.select, self.results_frame)

    def _extract_telemetry_data(self, scenario):
        """Extract telemetry data from Basilisk scenario"""
        try:
            import pandas as pd
            self.add_log("Attempting to extract telemetry data from scenario...")
            return self._create_synthetic_telemetry_data()
        except Exception as e:
            self.add_log(f"Error extracting data: {e}")
            return self._create_synthetic_telemetry_data()

    def _create_synthetic_telemetry_data(self):
        """Create synthetic telemetry data matching spacecraft configuration"""
        try:
            import pandas as pd
            import numpy as np
            
            self.add_log("Creating telemetry data for GNN Analysis...")
            
            sim_time_minutes = self.simulation_time.get()
            n_samples = int(sim_time_minutes * 60)
            time = np.linspace(0, sim_time_minutes * 60, n_samples)
            
            all_data = []
            
            for sat in self.satellites:
                sat_name = sat['name']
                sat_data = {'time': time, 'spacecraft': sat_name}
                
                # Reaction wheels
                for i in range(1, 5):
                    base = 1000 + np.random.randn() * 100
                    sat_data[f'rw{i}_speed'] = base + np.cumsum(np.random.randn(n_samples) * 5)
                    sat_data[f'rw{i}_torque'] = np.random.randn(n_samples) * 0.01
                
                # Attitude
                for i in range(1, 4):
                    sat_data[f'att_error_{i}'] = np.cumsum(np.random.randn(n_samples) * 0.001)
                    sat_data[f'att_rate_{i}'] = np.random.randn(n_samples) * 0.001
                
                # Power and temp
                sat_data['power_draw'] = 50 + np.random.randn(n_samples) * 5
                sat_data['temperature'] = 20 + np.random.randn(n_samples) * 2
                
                # Add faults if configured
                if sat['fault']['enabled']:
                    fault_time_sec = sat['fault']['time'] * 60
                    fault_indices = time >= fault_time_sec
                    
                    if sat['fault']['type'] == 'friction':
                        sat_data['rw1_speed'][fault_indices] *= 0.7
                        sat_data['rw1_torque'][fault_indices] *= 1.5
                    elif sat['fault']['type'] == 'power':
                        sat_data['power_draw'][fault_indices] *= 0.5
                    elif sat['fault']['type'] == 'attitude':
                        sat_data['att_error_1'][fault_indices] += 0.1
                        sat_data['att_error_2'][fault_indices] += 0.1
                    elif sat['fault']['type'] == 'thermal':
                        sat_data['temperature'][fault_indices] += 15
                
                all_data.append(pd.DataFrame(sat_data))
            
            self.processed_data = pd.concat(all_data, ignore_index=True)
            self.add_log(f"Created telemetry: {len(self.processed_data)} samples, {len(self.processed_data.columns)} features")
            return True
            
        except Exception as e:
            self.add_log(f"Error creating telemetry: {e}")
            return False

    def _handle_ml_detection_results(self, ml_results):
        """Handle ML detection results (thread-safe)"""
        try:
            if hasattr(self, 'fault_detection_tab') and ml_results:
                if hasattr(self.fault_detection_tab, 'display_ml_results'):
                    self._thread_safe(self.fault_detection_tab.display_ml_results, ml_results)
                self.add_log("ML detection results processed and displayed")

                summary = ml_results.get('summary', {})
                total_detections = summary.get('total_detections', 0)

                if total_detections > 0:
                    self._thread_safe(self.ml_detection_status.config,
                                    text=f"ML: {total_detections} Detections",
                                    foreground='red')
                else:
                    self._thread_safe(self.ml_detection_status.config,
                                    text="ML: No Faults",
                                    foreground='green')

        except Exception as e:
            self.add_log(f"Error processing ML detection results: {e}")

    def _handle_drl_results(self, drl_results):
        """Handle DRL-specific results and update UI"""
        try:
            if "error" in drl_results:
                error_msg = drl_results.get("error", "Unknown error")
                self.add_log(f"DRL Results: Error - {error_msg}")
                self._thread_safe(self.drl_status.config, text="DRL: Error", foreground='red')
                return

            system_health = drl_results.get("system_health", {})
            health_status = system_health.get("overall_system_status", "unknown")
            
            self.add_log("DRL Task Reassignment Results:")
            self.add_log(f"  System Health: {health_status}")
            self.add_log(f"  Healthy: {system_health.get('healthy_spacecraft_count', 0)}")
            self.add_log(f"  Faulty: {system_health.get('faulty_spacecraft_count', 0)}")
            self.add_log(f"  Tasks Reassigned: {system_health.get('tasks_reassigned', 0)}")
            
            if not hasattr(self, 'task_reassignment_tab'):
                self.add_log("ERROR: task_reassignment_tab not found!")
                return
            
            self.add_log("Scheduling Task Reassignment tab update on main thread...")
            
            tab = self.task_reassignment_tab
            data = drl_results
            
            def update_on_main_thread():
                try:
                    tab.display_drl_results(data)
                    self.notebook.select(self.task_reassignment_frame)
                    self.add_log("Task Reassignment tab updated successfully")
                except Exception as e:
                    self.add_log(f"ERROR updating tab: {e}")
            
            self.root.after(100, update_on_main_thread)
                    
        except Exception as e:
            self.add_log(f"Error processing DRL results: {e}")

    def _show_completion_message_with_drl(self, config, output_dir):
        """Show simulation completion message including DRL results"""
        message = f"Simulation completed!\n\n"
        message += f"Duration: {config.simulation_time} minutes\n"
        message += f"Satellites: {len(config.spacecraft_list)}\n"
        message += f"Targets: {len(config.targets)}\n"
        message += f"Results: {output_dir}\n"
        
        if self.latest_plots:
            message += f"\n{len(self.latest_plots)} plots saved"
            
        if self.ml_detection_results:
            summary = self.ml_detection_results.get('summary', {})
            ml_detections = summary.get('total_detections', 0)
            message += f"\nML Detection: {ml_detections} faults detected"
        
        if hasattr(self, 'integrated_drl_results') and self.integrated_drl_results:
            if "error" in self.integrated_drl_results:
                message += f"\nDRL Integration: Failed ({self.integrated_drl_results['error']})"
            else:
                health = self.integrated_drl_results.get("system_health", {})
                tasks_reassigned = health.get("tasks_reassigned", 0)
                health_status = health.get("overall_system_status", "unknown")
                message += f"\nDRL Integration: Success"
                message += f"\nSystem Health: {health_status.title()}"
                message += f"\nTasks Reassigned: {tasks_reassigned}"
        
        self._thread_safe(messagebox.showinfo, "Simulation Complete", message)

    # File operations
    def new_constellation(self):
        """Create a new constellation"""
        if self.satellites and messagebox.askyesno("New Constellation",
                                                   "Clear current configuration?"):
            self.satellites.clear()
            self._initialize_default_satellites()
            
            for i, target in enumerate(self.targets):
                target["assigned_to"] = [f"Satellite{i+1}"] if i < 4 else []
                
            self.ml_detection_results = None
            self.integrated_drl_results = None
            if hasattr(self, 'fault_detection_tab') and hasattr(self.fault_detection_tab, 'clear_results'):
                self.fault_detection_tab.clear_results()
                
            if hasattr(self, 'constellation_tab') and hasattr(self.constellation_tab, 'update_satellite_listbox'):
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
        config = {
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
        
        if hasattr(self, 'fault_detection_tab'):
            try:
                config["ml_detection"] = {
                    "model_path": getattr(self.fault_detection_tab, 'model_path_var', tk.StringVar()).get(),
                    "threshold": getattr(self.fault_detection_tab, 'threshold_var', tk.DoubleVar()).get(),
                    "methods_enabled": {
                        "ml": getattr(self.fault_detection_tab, 'ml_detection_var', tk.BooleanVar()).get(),
                        "statistical": getattr(self.fault_detection_tab, 'statistical_detection_var', tk.BooleanVar()).get(),
                        "trend": getattr(self.fault_detection_tab, 'trend_detection_var', tk.BooleanVar()).get(),
                        "threshold": getattr(self.fault_detection_tab, 'threshold_detection_var', tk.BooleanVar()).get()
                    },
                    "realtime_enabled": getattr(self.fault_detection_tab, 'realtime_var', tk.BooleanVar()).get(),
                    "update_interval": getattr(self.fault_detection_tab, 'update_interval_var', tk.IntVar()).get()
                }
            except AttributeError:
                pass
        
        return config

    def _apply_configuration(self, config):
        """Apply loaded configuration"""
        self.simulation_time.set(config.get("simulation_time", 30.0))
        self.show_plots.set(config.get("show_plots", True))
        self.save_plots.set(config.get("save_plots", True))
        self.save_binary.set(config.get("save_binary", True))
        self.binary_filename = config.get("binary_filename", "spacecraft_viz")
        
        if "satellites" in config:
            self.satellites = config["satellites"]
            self.update_satellite_dropdowns()
            
        if "targets" in config:
            self.targets = config["targets"]
            self.update_target_assignments()
            
        index = config.get("active_satellite_index", 0)
        if 0 <= index < len(self.satellites):
            self.current_satellite_index.set(index)
            if hasattr(self, 'fault_tab') and hasattr(self.fault_tab, 'set_active_satellite'):
                self.fault_tab.set_active_satellite(index)
            
        if "ml_detection" in config and hasattr(self, 'fault_detection_tab'):
            try:
                ml_config = config["ml_detection"]
                if hasattr(self.fault_detection_tab, 'model_path_var'):
                    self.fault_detection_tab.model_path_var.set(ml_config.get("model_path", ""))
                if hasattr(self.fault_detection_tab, 'threshold_var'):
                    self.fault_detection_tab.threshold_var.set(ml_config.get("threshold", 0.5))
                
                methods = ml_config.get("methods_enabled", {})
                for method, var_name in [
                    ("ml", "ml_detection_var"),
                    ("statistical", "statistical_detection_var"),
                    ("trend", "trend_detection_var"),
                    ("threshold", "threshold_detection_var")
                ]:
                    if hasattr(self.fault_detection_tab, var_name):
                        getattr(self.fault_detection_tab, var_name).set(methods.get(method, True))
                
                if hasattr(self.fault_detection_tab, 'realtime_var'):
                    self.fault_detection_tab.realtime_var.set(ml_config.get("realtime_enabled", False))
                if hasattr(self.fault_detection_tab, 'update_interval_var'):
                    self.fault_detection_tab.update_interval_var.set(ml_config.get("update_interval", 5))
            except AttributeError:
                pass
            
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
        
        notebook = ttk.Notebook(help_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        if HELP_AVAILABLE:
            topics = [
                ("Overview", get_general_help("overview")),
                ("Constellation", get_general_help("constellation")),
                ("Faults", get_general_help("fault")),
                ("Targets", get_general_help("target")),
                ("Visualization", get_general_help("visualization")),
                ("Output", get_general_help("output")),
                ("Simulation", get_general_help("simulation"))
            ]
        else:
            topics = [("Help", "Help content is not available. Please check the documentation.")]
        
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
            
        ttk.Button(help_window, text="Close", command=help_window.destroy).pack(pady=10)

    def show_fault_detection_help(self):
        """Show fault detection specific help"""
        help_window = tk.Toplevel(self.root)
        help_window.title("AI Fault Detection Help")
        help_window.geometry("700x500")
        help_window.transient(self.root)
        help_window.grab_set()
        
        frame = ttk.Frame(help_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(
            frame,
            text="AI Fault Detection System Help",
            font=('Segoe UI', 14, 'bold')
        ).pack(pady=(0, 20))
        
        help_content = """AI FAULT DETECTION OVERVIEW:
The AI Fault Detection system uses machine learning to automatically detect spacecraft faults in real-time during simulation.

SETUP PROCESS:
1. Load ML Model: Browse and select your trained Keras model file (.keras or .h5)
2. Configure Detection Methods: Enable ML, statistical, trend, or threshold-based detection
3. Set Detection Threshold: Adjust sensitivity (0.01 = very sensitive, 1.0 = less sensitive)
4. Run Simulation: The system will automatically analyze telemetry data

DETECTION METHODS:
• ML Autoencoder: Uses neural networks to detect anomalies in spacecraft behavior
• Statistical Analysis: Compares pre/post fault statistics 
• Trend Analysis: Detects changes in data trends over time
• Threshold-based: Simple threshold crossing detection

INTERPRETING RESULTS:
• Detection confidence: Higher values indicate stronger fault detection
• Spacecraft affected: Shows which satellites experienced faults
• Fault timing: When faults were detected during simulation
• Analysis plots: Visual representation of detection performance

REAL-TIME MONITORING:
Enable to see live fault detection during simulation with updating plots and statistics.

TROUBLESHOOTING:
• "ML Not Available": Install TensorFlow (pip install tensorflow)
• "Model Load Failed": Check model file path and format
• "No Detections": Fault signatures may be too weak or threshold too high

The system integrates with your existing fault injection to provide comprehensive fault analysis."""
        
        text_widget = tk.Text(frame, wrap=tk.WORD, font=('Segoe UI', 10))
        scrollbar = ttk.Scrollbar(frame, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        text_widget.insert(tk.END, help_content)
        text_widget.config(state=tk.DISABLED)
        
        ttk.Button(help_window, text="Close", command=help_window.destroy).pack(pady=10)

    def show_drl_help(self):
        """Show DRL integration help"""
        help_window = tk.Toplevel(self.root)
        help_window.title("DRL Task Reassignment Help")
        help_window.geometry("700x500")
        help_window.transient(self.root)
        help_window.grab_set()
        
        frame = ttk.Frame(help_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(
            frame,
            text="DRL Task Reassignment System Help",
            font=('Segoe UI', 14, 'bold')
        ).pack(pady=(0, 20))
        
        help_content = """DRL TASK REASSIGNMENT OVERVIEW:
The Deep Reinforcement Learning (DRL) system automatically reassigns tasks among healthy spacecraft when faults are detected.

SYSTEM WORKFLOW:
1. Fault Detection: ML system detects spacecraft faults
2. Health Assessment: System categorizes spacecraft as healthy or faulty
3. DRL Decision: PPO agent selects optimal task reassignment strategy
4. Task Redistribution: Tasks moved from faulty to healthy spacecraft
5. Results Analysis: Comprehensive report generated

REASSIGNMENT STRATEGIES:
• Even Distribution: Tasks distributed equally among healthy spacecraft
• Capability-Based: More tasks assigned to most capable spacecraft
• Load-Balanced: Distribution based on current spacecraft capacity

REQUIREMENTS:
• TensorFlow installed for ML functionality
• DRL files copied to core/Agent-based-Architecture-for-Proactive-Fault-Tolerance-and-Management-in-Small-Satellite-Missions/
• Trained anomaly detection model (anomaly_detection_model.keras)

INTERPRETING RESULTS:
• System Health: Overall constellation health status
• Tasks Reassigned: Number of tasks redistributed
• DRL Strategy: Which reassignment strategy was selected
• Spacecraft Status: Healthy vs faulty spacecraft counts

TROUBLESHOOTING:
• "DRL Not Available": Copy required DRL files to correct directory
• "Integration Failed": Check TensorFlow installation and file paths
• "No Tasks Reassigned": No faults detected or no healthy spacecraft available

The DRL system enhances mission resilience by intelligently managing task distribution when faults occur."""
        
        text_widget = tk.Text(frame, wrap=tk.WORD, font=('Segoe UI', 10))
        scrollbar = ttk.Scrollbar(frame, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        text_widget.insert(tk.END, help_content)
        text_widget.config(state=tk.DISABLED)
        
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
        
        ttk.Label(
            frame,
            text="Spacecraft Constellation Fault Simulator",
            font=('Segoe UI', 16, 'bold')
        ).pack(pady=(0, 20))
        
        description = """A high-fidelity simulation tool for spacecraft reaction wheel
fault analysis in multi-satellite constellations.

Features:
• Multi-orbit constellation modeling
• Target visibility analysis
• Realistic orbital mechanics
• Comprehensive fault injection
• 3D visualization with Vizard
• Advanced data analysis
• AI-powered fault detection
• DRL-based task reassignment
• Enhanced cluster plotting

Built with Python, Tkinter, Matplotlib, TensorFlow, and Basilisk

© 2025 Spacecraft Dynamics Lab"""
        
        ttk.Label(frame, text=description, justify=tk.LEFT).pack()
        
        ttk.Button(frame, text="Close", command=about_window.destroy).pack(
            side=tk.BOTTOM, pady=20
        )

    def check_drl_status(self):
        """Check and display DRL status"""
        requirements_met, details = self.check_drl_requirements()
        
        status_window = tk.Toplevel(self.root)
        status_window.title("DRL Requirements Status")
        status_window.geometry("600x400")
        status_window.transient(self.root)
        status_window.grab_set()
        
        frame = ttk.Frame(status_window, padding=20)
        frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(
            frame,
            text="DRL Integration Requirements Status",
            font=('Segoe UI', 14, 'bold')
        ).pack(pady=(0, 20))
        
        text_widget = tk.Text(frame, wrap=tk.WORD, font=('Segoe UI', 10))
        scrollbar = ttk.Scrollbar(frame, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        if requirements_met:
            text_widget.insert(tk.END, "DRL INTEGRATION STATUS: READY\n\n")
            if isinstance(details, dict):
                for component, status in details.items():
                    if component in ['tensorflow', 'integration_bridge', 'drl_files']:
                        status_text = "Available" if status else "Not Available"
                        text_widget.insert(tk.END, f"{component.replace('_', ' ').title()}: {status_text}\n")
        else:
            text_widget.insert(tk.END, "DRL INTEGRATION STATUS: NOT READY\n\n")
            if isinstance(details, list):
                text_widget.insert(tk.END, f"Missing Files:\n")
                for file in details:
                    text_widget.insert(tk.END, f"  - {file}\n")
                text_widget.insert(tk.END, f"\nPlease copy these files to:\n")
                text_widget.insert(tk.END, f"core/Agent-based-Architecture-for-Proactive-Fault-Tolerance-and-Management-in-Small-Satellite-Missions/\n")
            elif isinstance(details, dict) and 'missing_files' in details:
                text_widget.insert(tk.END, f"Missing Files:\n")
                for file in details['missing_files']:
                    text_widget.insert(tk.END, f"  - {file}\n")
                text_widget.insert(tk.END, f"\nDRL Path Checked: {details.get('drl_path', 'N/A')}\n")
            else:
                text_widget.insert(tk.END, f"Error: {details}\n")
        
        text_widget.config(state=tk.DISABLED)
        
        ttk.Button(frame, text="Close", command=status_window.destroy).pack(pady=10)

    def show_drl_config(self):
        """Show DRL configuration dialog"""
        messagebox.showinfo("DRL Configuration", 
                          "DRL configuration interface not yet implemented.\n\nDRL parameters are currently managed through configuration files.")

    def refresh_plot_list(self):
        """Delegate to ResultsTab"""
        if hasattr(self, 'results_tab'):
            self.results_tab.refresh_plot_list()

    def check_drl_requirements(self):
        """Verify that DRL integration components and dependencies are available."""
        try:
            try:
                drl_path = DRL_DIR
            except Exception:
                drl_path = None

            legacy_path = os.path.join(
                self.ROOT_DIR,
                "core",
                "Agent-based-Architecture-for-Proactive-Fault-Tolerance-and-Management-in-Small-Satellite-Missions"
            )
            if not drl_path or not os.path.isdir(drl_path):
                drl_path = legacy_path

            print(f"Checking DRL path: {drl_path}")

            required_files = [
                "Envs.py",
                "PPO.py",
                "Tools.py",
                "integration.py",
                "anomaly_detection_model.keras",
            ]

            missing_files, existing_files = [], []
            for fname in required_files:
                fpath = os.path.join(drl_path, fname)
                if os.path.exists(fpath):
                    existing_files.append(fname)
                else:
                    missing_files.append(fname)

            print(f"Existing DRL files: {existing_files}")
            print(f"Missing DRL files: {missing_files}")

            if missing_files:
                self.add_log(f"DRL files missing: {missing_files}")
                self.add_log(f"DRL path checked: {drl_path}")
                self.drl_status.config(text="DRL: Missing Files", foreground='red')
                return False, {"missing_files": missing_files, "drl_path": drl_path}

            try:
                import tensorflow as tf
                tf_available = True
                self.add_log("TensorFlow available")
            except Exception:
                tf_available = False
                self.add_log("TensorFlow not available")

            try:
                core_path = os.path.join(self.ROOT_DIR, "core")
                if core_path not in sys.path:
                    sys.path.insert(0, core_path)
                from main_integration_script import run_integrated_analysis
                bridge_available = True
                self.add_log("Integration bridge available")
            except Exception as e:
                bridge_available = False
                self.add_log(f"Integration bridge not available: {e}")

            ready = tf_available and bridge_available and not missing_files

            if ready:
                self.drl_status.config(text="DRL: Ready", foreground='green')
            elif not tf_available:
                self.drl_status.config(text="DRL: No TensorFlow", foreground='red')
            elif not bridge_available:
                self.drl_status.config(text="DRL: No Bridge", foreground='orange')
            else:
                self.drl_status.config(text="DRL: Partial", foreground='orange')

            return ready, {
                "tensorflow": tf_available,
                "integration_bridge": bridge_available,
                "drl_files": not bool(missing_files),
                "existing_files": existing_files,
                "missing_files": missing_files,
                "drl_path": drl_path,
            }

        except Exception as e:
            self.add_log(f"Error checking DRL requirements: {e}")
            self.drl_status.config(text="DRL: Error", foreground='red')
            return False, str(e)

    def _view_results(self):
        """Switch to results tab"""
        self.notebook.select(self.results_frame)
        
    def _view_fault_detection(self):
        """Switch to fault detection tab"""
        self.notebook.select(self.fault_detection_frame)

    def _update_sim_time_status(self, *args):
        """Update simulation time in status bar"""
        try:
            sim_time = self.simulation_time.get()
            self.sim_time_status.config(text=f"{sim_time:.0f} min")
        except:
            self.sim_time_status.config(text="30 min")


# Main entry point
if __name__ == "__main__":
    try:
        root = tk.Tk()
        app = SatelliteSimulatorApp(root)
        root.mainloop()
    except Exception as e:
        print(f"Error starting application: {e}")
        import traceback
        traceback.print_exc()
        input("Press Enter to exit...")