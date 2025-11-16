#!/usr/bin/env python3
"""
drl_tab_enhanced.py — FIXED DRL Tab with Proper GUI Integration

KEY FIXES:
1. Pulls actual spacecraft names from constellation_tab
2. Passes cluster configuration to training
3. Passes fault configuration to training
4. Updates Integration Status tab with real data
5. Shows actual names in training visualizations
"""

import os
import sys
import threading
import queue
import subprocess
import time
import datetime
import platform
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np

# Add paths
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
root_dir = os.path.dirname(parent_dir)
drl_dir = os.path.join(root_dir, "DRL")

for path in [parent_dir, root_dir, drl_dir]:
    if path not in sys.path:
        sys.path.insert(0, path)

# Optional imports
try:
    import pandas as pd
    PD_AVAILABLE = True
except Exception:
    pd = None
    PD_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    MPL_AVAILABLE = True
except Exception:
    Figure = None
    FigureCanvasTkAgg = None
    MPL_AVAILABLE = False


class DRLTab(ttk.Frame):
    """
    ENHANCED DRL Tab with proper GUI integration
    """

    def __init__(self, parent, parent_app=None, default_script=None, default_results_dir=None):
        super().__init__(parent)
        self.parent_app = parent_app
        self.proc = None
        self.worker_thread = None
        self.queue = queue.Queue()
        self._stop_reader = threading.Event()

        # State
        self.script_path = tk.StringVar(value=default_script or "")
        self.results_dir = tk.StringVar(value=default_results_dir or os.path.join(drl_dir, "results"))
        rd = os.path.normpath(self.results_dir.get())
        if rd == os.path.normpath(drl_dir):
            rd = os.path.join(drl_dir, "results")
            self.results_dir.set(rd)

        os.makedirs(self.results_dir.get(), exist_ok=True)
        self.status_text = tk.StringVar(value="Idle")
        self._log_lines = 0
        
        # Integration state
        self.last_scenario = None
        self.last_config = None
        self.fault_detection_results = None
        self.drl_results = None

        # Layout
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=5, pady=5)
        self._nb = nb

        self.overview = ttk.Frame(nb)
        self.training = ttk.Frame(nb)
        self.integration = ttk.Frame(nb)
        self.results = ttk.Frame(nb)
        
        nb.add(self.overview, text="Overview")
        nb.add(self.training, text="PPO Training")
        nb.add(self.integration, text="Integration Status")
        nb.add(self.results, text="Results")

        self._build_overview(self.overview)
        self._build_training(self.training)
        self._build_integration(self.integration)
        self._build_results(self.results)

        # Poll queue
        self.after(150, self._drain_queue)

    # ============================= OVERVIEW TAB =============================
    def _build_overview(self, root):
        """Build overview tab"""
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root, text="DRL Task Reassignment System",
                  font=("TkDefaultFont", 14, "bold")).pack(anchor="w", **pad)
        
        # Integration status
        status_frame = ttk.LabelFrame(root, text="Integration Status", padding=10)
        status_frame.pack(fill="x", **pad)
        
        components = [
            ("Constellation Tab", hasattr(self.parent_app, 'constellation_tab') if self.parent_app else False),
            ("Fault Tab", hasattr(self.parent_app, 'fault_tab') if self.parent_app else False),
            ("Fault Detection Tab", hasattr(self.parent_app, 'fault_detection_tab') if self.parent_app else False),
            ("Task Reassignment Tab", hasattr(self.parent_app, 'task_reassignment_tab') if self.parent_app else False),
        ]
        
        for name, available in components:
            row = ttk.Frame(status_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"{name}:", width=25).pack(side="left")
            status = "✓ Connected" if available else "✗ Not Available"
            color = "green" if available else "red"
            lbl = ttk.Label(row, text=status, foreground=color)
            lbl.pack(side="left")
        
        # Configuration summary
        config_frame = ttk.LabelFrame(root, text="Current Configuration", padding=10)
        config_frame.pack(fill="x", **pad)
        
        row = ttk.Frame(config_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Spacecraft:", width=25).pack(side="left")
        self.spacecraft_count_label = ttk.Label(row, text="Not configured")
        self.spacecraft_count_label.pack(side="left")
        
        row = ttk.Frame(config_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Spacecraft Names:", width=25).pack(side="left")
        self.spacecraft_names_label = ttk.Label(row, text="Not configured")
        self.spacecraft_names_label.pack(side="left")
        
        row = ttk.Frame(config_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Faults Configured:", width=25).pack(side="left")
        self.faults_count_label = ttk.Label(row, text="Not configured")
        self.faults_count_label.pack(side="left")
        
        row = ttk.Frame(config_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Clusters:", width=25).pack(side="left")
        self.clusters_count_label = ttk.Label(row, text="Not configured")
        self.clusters_count_label.pack(side="left")
        
        ttk.Button(config_frame, text="Refresh Configuration", 
                  command=self._refresh_config_summary).pack(pady=(10, 0))
        
        # Info text
        ttk.Separator(root, orient="horizontal").pack(fill="x", pady=10)
        
        info_text = """ENHANCED DRL TASK REASSIGNMENT:

✓ Now pulls actual spacecraft names from Constellation tab
✓ Uses real cluster configuration
✓ Uses real fault scenarios
✓ Training visualizations show YOUR spacecraft names

WORKFLOW:
1. Configure spacecraft in Constellation tab
2. Configure faults in Fault tab
3. Click "Refresh Configuration" to pull latest config
4. Start Training to train DRL with your configuration
5. Run Simulation to see DRL in action
6. View results in Task Reassignment tab

WHAT'S FIXED:
• Training now uses actual spacecraft names (not Sat1-Sat8)
• Load distribution charts show real names
• Cluster configuration is passed to training
• Fault scenarios are used in training environment
        """
        
        info_box = tk.Text(root, height=15, wrap="word", font=("TkDefaultFont", 9))
        info_box.pack(fill="both", expand=True, **pad)
        info_box.insert("1.0", info_text)
        info_box.configure(state="disabled", bg="#f0f0f0")
        
        # Auto-refresh
        self._refresh_config_summary()

    # ============================= TRAINING TAB =============================
    def _build_training(self, root):
        """Build training tab"""
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root, text="PPO Training Configuration",
                  font=("TkDefaultFont", 12, "bold")).pack(anchor="w", **pad)

        # Configuration source
        config_frame = ttk.LabelFrame(root, text="Configuration Source", padding=10)
        config_frame.pack(fill="x", **pad)
        
        ttk.Label(config_frame, text="✓ Training will use current GUI configuration",
                 font=("TkDefaultFont", 9, "bold"), foreground="green").pack(anchor="w", pady=(0, 5))
        
        ttk.Label(config_frame, text="This includes:",
                 font=("TkDefaultFont", 9)).pack(anchor="w", pady=(5, 2))
        ttk.Label(config_frame, text="  • Actual spacecraft names from Constellation tab",
                 font=("TkDefaultFont", 9)).pack(anchor="w")
        ttk.Label(config_frame, text="  • Cluster configurations",
                 font=("TkDefaultFont", 9)).pack(anchor="w")
        ttk.Label(config_frame, text="  • Fault scenarios from Fault tab",
                 font=("TkDefaultFont", 9)).pack(anchor="w")
        
        # Training parameters
        param_frame = ttk.LabelFrame(root, text="Training Parameters", padding=10)
        param_frame.pack(fill="x", **pad)
        
        row = ttk.Frame(param_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Training Episodes:", width=20).pack(side="left")
        self.episodes_var = tk.IntVar(value=100)
        ttk.Spinbox(row, from_=10, to=1000, textvariable=self.episodes_var, width=10).pack(side="left")
        
        row = ttk.Frame(param_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Max Steps/Episode:", width=20).pack(side="left")
        self.max_steps_var = tk.IntVar(value=20)
        ttk.Spinbox(row, from_=5, to=100, textvariable=self.max_steps_var, width=10).pack(side="left")

        # Control buttons
        ctrl_frame = ttk.Frame(root)
        ctrl_frame.pack(fill="x", **pad)
        
        self.train_btn = ttk.Button(ctrl_frame, text="Start Training", command=self._start_training)
        self.train_btn.pack(side="left", padx=(0, 5))
        
        ttk.Button(ctrl_frame, text="Stop", command=self._stop_training).pack(side="left", padx=(0, 10))
        ttk.Label(ctrl_frame, textvariable=self.status_text).pack(side="left")

        # Log
        log_frame = ttk.Frame(root)
        log_frame.pack(fill="both", expand=True, **pad)
        
        ttk.Label(log_frame, text="Training Log:").pack(anchor="w")
        
        text_frame = ttk.Frame(log_frame)
        text_frame.pack(fill="both", expand=True, pady=(2, 0))
        
        self.train_log = tk.Text(text_frame, height=12, state="disabled", wrap="word")
        scrollbar = ttk.Scrollbar(text_frame, command=self.train_log.yview)
        self.train_log.configure(yscrollcommand=scrollbar.set)
        
        self.train_log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    # ============================= INTEGRATION TAB =============================
    def _build_integration(self, root):
        """Build integration status tab"""
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root, text="Integration with Simulation Pipeline",
                  font=("TkDefaultFont", 12, "bold")).pack(anchor="w", **pad)

        # Last simulation info
        sim_frame = ttk.LabelFrame(root, text="Last Simulation Run", padding=10)
        sim_frame.pack(fill="x", **pad)
        
        row = ttk.Frame(sim_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Status:", width=20).pack(side="left")
        self.sim_status_label = ttk.Label(row, text="No simulation run yet")
        self.sim_status_label.pack(side="left")
        
        row = ttk.Frame(sim_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Spacecraft:", width=20).pack(side="left")
        self.sim_spacecraft_label = ttk.Label(row, text="-")
        self.sim_spacecraft_label.pack(side="left")
        
        row = ttk.Frame(sim_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Spacecraft Names:", width=20).pack(side="left")
        self.sim_names_label = ttk.Label(row, text="-")
        self.sim_names_label.pack(side="left")
        
        row = ttk.Frame(sim_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Faults Detected:", width=20).pack(side="left")
        self.sim_faults_label = ttk.Label(row, text="-")
        self.sim_faults_label.pack(side="left")
        
        # DRL results summary
        drl_frame = ttk.LabelFrame(root, text="DRL Task Reassignment Results", padding=10)
        drl_frame.pack(fill="x", **pad)
        
        row = ttk.Frame(drl_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="System Health:", width=20).pack(side="left")
        self.health_label = ttk.Label(row, text="-")
        self.health_label.pack(side="left")
        
        row = ttk.Frame(drl_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Tasks Reassigned:", width=20).pack(side="left")
        self.tasks_label = ttk.Label(row, text="-")
        self.tasks_label.pack(side="left")
        
        row = ttk.Frame(drl_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Healthy Spacecraft:", width=20).pack(side="left")
        self.healthy_label = ttk.Label(row, text="-")
        self.healthy_label.pack(side="left")
        
        row = ttk.Frame(drl_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Faulty Spacecraft:", width=20).pack(side="left")
        self.faulty_label = ttk.Label(row, text="-")
        self.faulty_label.pack(side="left")
        
        ttk.Button(drl_frame, text="View Detailed Results in Task Reassignment Tab",
                  command=self._view_task_reassignment).pack(pady=(10, 0))
        
        # Integration log
        log_frame = ttk.LabelFrame(root, text="Integration Log", padding=10)
        log_frame.pack(fill="both", expand=True, **pad)
        
        text_frame = ttk.Frame(log_frame)
        text_frame.pack(fill="both", expand=True)
        
        self.integration_log = tk.Text(text_frame, height=10, state="disabled", wrap="word")
        scrollbar = ttk.Scrollbar(text_frame, command=self.integration_log.yview)
        self.integration_log.configure(yscrollcommand=scrollbar.set)
        
        self.integration_log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    # ============================= RESULTS TAB =============================
    def _build_results(self, root):
        """Build results tab"""
        pad = {"padx": 10, "pady": 6}

        hdr = ttk.Frame(root)
        hdr.pack(fill="x", **pad)
        
        ttk.Label(hdr, text="Results Folder:").pack(side="left")
        self.res_entry = ttk.Entry(hdr, textvariable=self.results_dir)
        self.res_entry.pack(side="left", fill="x", expand=True, padx=(6,6))
        ttk.Button(hdr, text="Refresh", command=self._refresh_results).pack(side="left", padx=(0, 5))
        ttk.Button(hdr, text="Open Folder", command=self._open_results_folder).pack(side="left")

        # Treeview
        cols = ("name","size","modified")
        tvf = ttk.Frame(root)
        tvf.pack(fill="x", expand=False, **pad)
        
        self.tree = ttk.Treeview(tvf, columns=cols, show="headings", selectmode="browse", height=10)
        for c, w in zip(cols, (360, 100, 180)):
            self.tree.heading(c, text=c.title())
            self.tree.column(c, width=w, anchor="w")
        
        tree_scroll = ttk.Scrollbar(tvf, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        
        self.tree.pack(side="left", fill="both", expand=True)
        tree_scroll.pack(side="right", fill="y")

        # Preview
        self.preview_frame = ttk.LabelFrame(root, text="Preview")
        self.preview_frame.pack(fill="both", expand=True, **pad)

        if MPL_AVAILABLE and Figure and FigureCanvasTkAgg:
            self.fig = Figure(figsize=(20, 10.2), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
            self.canvas.get_tk_widget().pack(fill="x", expand=True)
        else:
            self.fig = None
            ttk.Label(self.preview_frame, text="Matplotlib not available").pack(pady=20)

        self.tree.bind("<<TreeviewSelect>>", self._on_select_artifact)
        
        if self.results_dir.get() and os.path.isdir(self.results_dir.get()):
            self._refresh_results()

    # ============================= CORE METHODS =============================
    def _refresh_config_summary(self):
        """Refresh configuration summary - PULL FROM GUI"""
        if not self.parent_app:
            return
        
        try:
            # Spacecraft count and names
            if hasattr(self.parent_app, 'satellites'):
                count = len(self.parent_app.satellites)
                names = [sat.get('name', f'Sat{i+1}') for i, sat in enumerate(self.parent_app.satellites)]
                self.spacecraft_count_label.config(text=f"{count} spacecraft configured")
                self.spacecraft_names_label.config(text=f"{', '.join(names[:3])}{'...' if len(names) > 3 else ''}")
            
            # Faults
            if hasattr(self.parent_app, 'satellites'):
                faults = sum(1 for sat in self.parent_app.satellites 
                           if sat.get('fault', {}).get('enabled', False))
                self.faults_count_label.config(text=f"{faults} faults configured")
            
            # Clusters
            if hasattr(self.parent_app, 'constellation_tab') and hasattr(self.parent_app.constellation_tab, 'clusters'):
                clusters = len(self.parent_app.constellation_tab.clusters)
                self.clusters_count_label.config(text=f"{clusters} clusters configured")
        except Exception as e:
            print(f"Error refreshing config: {e}")

    def _start_training(self):
        """Start PPO training WITH GUI CONFIGURATION"""
        if not self.parent_app or not hasattr(self.parent_app, 'satellites'):
            messagebox.showerror("Configuration Error", 
                            "No spacecraft configuration available.\nPlease configure constellation first.")
            return
        
        if len(self.parent_app.satellites) == 0:
            messagebox.showerror("Configuration Error", 
                            "No spacecraft configured.\nPlease add satellites in Constellation tab.")
            return
        
        # Find enhanced PPO script
        ppo_script = os.path.join(drl_dir, "PPO_simple_enhanced.py")
        if not os.path.exists(ppo_script):
            messagebox.showerror("Script Not Found", 
                            f"Enhanced PPO script not found at:\n{ppo_script}\n\n"
                            f"Please ensure PPO_simple_enhanced.py is in the DRL/ folder.")
            return

        # Extract configuration from GUI
        spacecraft_names = [sat.get('name') for sat in self.parent_app.satellites]
        
        # Cluster configuration
        cluster_config = {}
        if hasattr(self.parent_app, 'constellation_tab') and hasattr(self.parent_app.constellation_tab, 'clusters'):
            for cluster in self.parent_app.constellation_tab.clusters:
                cluster_config[cluster['name']] = {
                    'formation': cluster.get('formation', 'Leader-Follower'),
                    'satellites': cluster.get('satellites', [])
                }
        
        # Fault configuration
        fault_config = {}
        for sat in self.parent_app.satellites:
            if sat.get('fault', {}).get('enabled'):
                fault_config[sat['name']] = {
                    'enabled': True,
                    'type': sat['fault'].get('type'),
                    'time': sat['fault'].get('time'),
                    'magnitude': sat['fault'].get('magnitude')
                }

        self._append_train_log(f"\n=== Starting Enhanced PPO Training ===\n")
        self._append_train_log(f"Configuration:\n")
        self._append_train_log(f"  • Spacecraft: {len(spacecraft_names)}\n")
        self._append_train_log(f"  • Names: {', '.join(spacecraft_names)}\n")
        self._append_train_log(f"  • Clusters: {len(cluster_config)}\n")
        self._append_train_log(f"  • Faults: {len(fault_config)}\n")
        self._append_train_log(f"  • Episodes: {self.episodes_var.get()}\n")
        self._append_train_log(f"\nScript: {ppo_script}\n\n")
        
        # Run with configuration
        self._run_script_with_config(ppo_script, spacecraft_names, cluster_config, fault_config)

    def _run_script_with_config(self, script_path, spacecraft_names, cluster_config, fault_config):
        """Run script with GUI configuration passed via environment variables"""
        python = sys.executable
        cmd = [python, script_path]
        
        # Add environment variables with GUI configuration
        env = os.environ.copy()
        env['DRL_EPISODES'] = str(self.episodes_var.get())
        env['DRL_MAX_STEPS'] = str(self.max_steps_var.get())
        env['DRL_SPACECRAFT_NAMES'] = json.dumps(spacecraft_names)
        env['DRL_CLUSTER_CONFIG'] = json.dumps(cluster_config)
        env['DRL_FAULT_CONFIG'] = json.dumps(fault_config)

        try:
            self.status_text.set("Running...")
            self.train_btn.configure(state="disabled")
            self._stop_reader.clear()

            self.proc = subprocess.Popen(
                cmd,
                cwd=os.path.dirname(script_path) or None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True,
                env=env
            )
        except Exception as e:
            messagebox.showerror("Error", f"Failed to start:\n{e}")
            self.status_text.set("Failed")
            self.train_btn.configure(state="normal")
            return

        self.worker_thread = threading.Thread(target=self._reader_worker, daemon=True)
        self.worker_thread.start()

    def _reader_worker(self):
        """Read subprocess output"""
        try:
            for line in self.proc.stdout:
                if self._stop_reader.is_set():
                    break
                self.queue.put(("log", line))
        except Exception as e:
            self.queue.put(("log", f"[reader] error: {e}\n"))
        finally:
            rc = None
            try:
                rc = self.proc.wait(timeout=1.0)
            except Exception:
                pass
            self.queue.put(("log", f"\n[Process finished with code: {rc}]\n"))
            self.queue.put(("done", None))

    def _drain_queue(self):
        """Drain output queue"""
        drained = False
        while True:
            try:
                msg_type, item = self.queue.get_nowait()
            except queue.Empty:
                break
            drained = True
            
            if msg_type == "done":
                self.status_text.set("Idle")
                self.train_btn.configure(state="normal")
                self._refresh_results()
                continue
            
            if msg_type == "log":
                self._append_train_log(item)

        self.after(150 if drained else 250, self._drain_queue)

    def _stop_training(self):
        """Stop training"""
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
                self._append_train_log("\n[Training stopped by user]\n")
            except Exception as e:
                self._append_train_log(f"\n[Stop failed: {e}]\n")
        self.status_text.set("Stopped")
        self._stop_reader.set()

    def update_simulation_results(self, scenario, config, ml_results=None, drl_results=None):
        """Called by main GUI after simulation - UPDATE INTEGRATION STATUS"""
        self.last_scenario = scenario
        self.last_config = config
        self.fault_detection_results = ml_results
        self.drl_results = drl_results
        
        # Update integration display
        self._update_integration_display()
        
        # Log
        self._log_integration(f"\n[{datetime.datetime.now().strftime('%H:%M:%S')}] Simulation completed")
        
        if ml_results:
            summary = ml_results.get('summary', {})
            detections = summary.get('total_detections', 0)
            self._log_integration(f"  • Fault Detection: {detections} faults detected")
        
        if drl_results:
            if "error" not in drl_results:
                health = drl_results.get("system_health", {})
                self._log_integration(f"  • DRL: {health.get('tasks_reassigned', 0)} tasks reassigned")
                self._log_integration(f"  • System Health: {health.get('overall_system_status', 'unknown')}")
            else:
                self._log_integration(f"  • DRL: Error - {drl_results.get('error')}")

    def _update_integration_display(self):
        """Update integration status display WITH SPACECRAFT NAMES"""
        if self.last_scenario:
            spacecraft_count = len(getattr(self.last_scenario, 'sc_objects', []))
            self.sim_status_label.config(text="✓ Completed", foreground="green")
            self.sim_spacecraft_label.config(text=f"{spacecraft_count} spacecraft")
            
            # Show spacecraft names
            if self.last_config and hasattr(self.last_config, 'spacecraft_list'):
                names = [sc.get('name') for sc in self.last_config.spacecraft_list]
                self.sim_names_label.config(text=f"{', '.join(names[:3])}{'...' if len(names) > 3 else ''}")
        
        if self.fault_detection_results:
            summary = self.fault_detection_results.get('summary', {})
            detections = summary.get('total_detections', 0)
            self.sim_faults_label.config(text=f"{detections} faults detected")
        
        if self.drl_results and "error" not in self.drl_results:
            health = self.drl_results.get("system_health", {})
            
            status = health.get("overall_system_status", "unknown")
            self.health_label.config(text=status.title())
            
            tasks = health.get("tasks_reassigned", 0)
            self.tasks_label.config(text=f"{tasks} tasks")
            
            healthy = health.get("healthy_spacecraft_count", 0)
            self.healthy_label.config(text=f"{healthy} spacecraft")
            
            faulty = health.get("faulty_spacecraft_count", 0)
            self.faulty_label.config(text=f"{faulty} spacecraft")

    def _view_task_reassignment(self):
        """Switch to Task Reassignment tab"""
        if self.parent_app and hasattr(self.parent_app, 'notebook'):
            try:
                for i in range(self.parent_app.notebook.index("end")):
                    if self.parent_app.notebook.tab(i, "text") == "Task Reassignment":
                        self.parent_app.notebook.select(i)
                        break
            except Exception as e:
                print(f"Error switching tabs: {e}")

    def _log_integration(self, message):
        """Log to integration log"""
        self.integration_log.configure(state="normal")
        self.integration_log.insert("end", message + "\n")
        self.integration_log.see("end")
        self.integration_log.configure(state="disabled")

    def _append_train_log(self, s: str):
        """Append to training log"""
        self.train_log.configure(state="normal")
        self.train_log.insert("end", s)
        self._log_lines += s.count("\n")
        if self._log_lines > 3000:
            self.train_log.delete("1.0", "500.0")
            self._log_lines -= 500
        self.train_log.see("end")
        self.train_log.configure(state="disabled")

    def _refresh_results(self):
        """Refresh results list"""
        folder = self.results_dir.get().strip()
        self.tree.delete(*self.tree.get_children())
        if not folder or not os.path.isdir(folder):
            return
        
        allow = {".xlsx", ".csv", ".png", ".jpg", ".json", ".txt", ".pth", ".pt", ".zip"}
        rows = []
        
        try:
            for name in os.listdir(folder):
                ext = os.path.splitext(name)[1].lower()
                if ext in allow:
                    p = os.path.join(folder, name)
                    try:
                        stat = os.stat(p)
                        size_kb = int(stat.st_size / 1024)
                        mtime = datetime.datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
                        rows.append((name, f"{size_kb} KB", mtime))
                    except Exception:
                        rows.append((name, "-", "-"))
        except Exception as e:
            print(f"Error refreshing results: {e}")
            return
        
        rows.sort(key=lambda r: r[2], reverse=True)
        for r in rows:
            self.tree.insert("", "end", values=r)

    def _on_select_artifact(self, event=None):
        """Handle artifact selection"""
        sel = self.tree.selection()
        if not sel:
            return
        name = self.tree.item(sel[0], "values")[0]
        path = os.path.join(self.results_dir.get().strip(), name)
        
        if name.endswith(('.png', '.jpg')):
            self._preview_image(path)
        elif name.endswith('.json'):
            self._preview_json(path)

    def _preview_image(self, path):
        """Preview image"""
        if self.fig is None:
            return
        try:
            from PIL import Image as PILImage
            img = PILImage.open(path)
            self.ax.clear()
            self.ax.imshow(img)
            self.ax.axis('off')
            self.ax.set_title(os.path.basename(path))
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error previewing image: {e}")

    def _preview_json(self, path):
        """Preview JSON"""
        if self.fig is None:
            return
        try:
            with open(path, 'r') as f:
                data = json.load(f)
            
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"JSON File\n{os.path.basename(path)}\n{len(str(data))} characters", 
                        ha="center", va="center", fontsize=10)
            self.ax.set_title("JSON Preview")
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error previewing JSON: {e}")

    def _open_results_folder(self):
        """Open results folder"""
        folder = self.results_dir.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showinfo("DRL", "No folder selected or folder not found.")
            return
        
        try:
            if platform.system() == "Windows":
                os.startfile(folder)
            elif platform.system() == "Darwin":
                subprocess.call(["open", folder])
            else:
                subprocess.call(["xdg-open", folder])
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open folder:\n{e}")