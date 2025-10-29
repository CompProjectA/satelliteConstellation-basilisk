#!/usr/bin/env python3
"""
drl_tab.py — Fully Integrated DRL Tab

Integrates with:
- GUI system (spacecraft_simulator_gui.py)
- DRL algorithms (PPO, TDHD, DQN from DRL/ folder)
- Fault detection (real_ml_fault_detection.py)
- Task reassignment (drl_integration_bridge.py)
- Spacecraft simulation (spacecraft_simulation.py)
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

# ---------- Path setup ----------
current_dir = os.path.dirname(os.path.abspath(__file__))   # .../BskSim/core
root_dir = os.path.dirname(current_dir)                    # .../BskSim
drl_dir = os.path.join(root_dir, "DRL")


# Ensure import paths
for path in [current_dir, root_dir, drl_dir]:
    if path not in sys.path:
        sys.path.insert(0, path)

# Pull canonical paths from drl_config (single source of truth)
try:
    from drl_config import RESULT_DIR, DRL_DIR, LOGS_DIR
except Exception:
    # Fallback if importing constants fails; keep previous logic
    RESULT_DIR = os.path.join(drl_dir, "result")
    DRL_DIR = drl_dir
    LOGS_DIR = os.path.join(root_dir, "logs")

# Make sure these key dirs exist
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

# ---------- Optional imports ----------
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
    from matplotlib.backends_backend_tkagg import FigureCanvasTkAgg  # type: ignore
    MPL_AVAILABLE = True
except Exception:
    # Some environments use the older module path:
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # type: ignore
        MPL_AVAILABLE = True
    except Exception:
        Figure = None
        FigureCanvasTkAgg = None
        MPL_AVAILABLE = False

# ---------- Integration components ----------
try:
    from drl_integration_bridge import (
        DRLTaskReassignmentSystem,
        integrate_fault_detection_with_drl,
    )
    DRL_BRIDGE_AVAILABLE = True
except Exception as e:
    DRL_BRIDGE_AVAILABLE = False
    print(f"DRL bridge not available: {e}")

try:
    from real_ml_fault_detection import run_real_ml_detection_on_scenario
    FAULT_DETECTION_AVAILABLE = True
except Exception:
    FAULT_DETECTION_AVAILABLE = False

try:
    from main_integration_script import IntegratedFaultDetectionDRLSystem
    INTEGRATED_SYSTEM_AVAILABLE = True
except Exception:
    INTEGRATED_SYSTEM_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("TkAgg")  # interactive in Tk
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    MPL_AVAILABLE = True
except Exception:
    Figure = None
    FigureCanvasTkAgg = None
    MPL_AVAILABLE = False

# Model manager import
try:
    from drl_model_manager import get_model_manager
    MODEL_MANAGER_AVAILABLE = True
except Exception as e:
    MODEL_MANAGER_AVAILABLE = False
    print(f"Model manager not available: {e}")


class DRLTab(ttk.Frame):
    """
    Fully integrated DRL tab with training, fault detection, and task reassignment
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
        # Default results dir → canonical DRL/result
        default_res = default_results_dir or RESULT_DIR
        self.results_dir = tk.StringVar(value=default_res)
        self.status_text = tk.StringVar(value="Idle")
        self.algorithm_var = tk.StringVar(value="PPO")
        self.training_mode = tk.BooleanVar(value=False)
        self.integration_mode = tk.BooleanVar(value=True)
        self._log_lines = 0

        # Integration state
        self.fault_detection_results = None
        self.drl_reassignment_results = None
        self.integrated_system = None

        # Ensure results dir exists early
        self._ensure_results_dir(self.results_dir.get())

        # Layout: nested notebook
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=5, pady=5)
        self._nb = nb

        self.overview = ttk.Frame(nb)
        self.training = ttk.Frame(nb)
        self.integration = ttk.Frame(nb)
        self.results = ttk.Frame(nb)

        nb.add(self.overview, text="Overview")
        nb.add(self.training, text="Training")
        nb.add(self.integration, text="Integration")
        nb.add(self.results, text="Results")

        self._build_overview(self.overview)
        self._build_training(self.training)
        self._build_integration(self.integration)
        self._build_results(self.results)

        # Poll the queue for subprocess output
        self.after(150, self._drain_queue)

    # ------------------------------ Overview ------------------------------
    def _build_overview(self, root):
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root, text="Deep Reinforcement Learning System",
                  font=("TkDefaultFont", 14, "bold")).pack(anchor="w", **pad)

        status_frame = ttk.LabelFrame(root, text="System Status", padding=10)
        status_frame.pack(fill="x", **pad)

        components = [
            ("DRL Bridge", DRL_BRIDGE_AVAILABLE),
            ("Fault Detection", FAULT_DETECTION_AVAILABLE),
            ("Integrated System", INTEGRATED_SYSTEM_AVAILABLE),
            ("Model Manager", MODEL_MANAGER_AVAILABLE),
            ("Matplotlib", MPL_AVAILABLE),
            ("Pandas", PD_AVAILABLE),
        ]

        for name, available in components:
            row = ttk.Frame(status_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"{name}:", width=20).pack(side="left")
            status = "✓ Available" if available else "✗ Not Available"
            color = "green" if available else "red"
            ttk.Label(row, text=status, foreground=color).pack(side="left")

        ttk.Separator(root, orient="horizontal").pack(fill="x", pady=10)

        info_text = """DRL System Components:

1. Training: Train PPO, TDHD, or DQN algorithms for task reassignment
2. Integration: Combine fault detection with DRL-based task reassignment
3. Results: View training metrics, fault detection outputs, and reassignment plans

Workflow:
- Use Training tab to train DRL agents (or load pre-trained models)
- Use Integration tab to run complete fault detection + task reassignment
- View results in the Results tab with plots and data files

Getting Started:
1. Click on the Training tab to configure and train DRL models
2. Or click Integration to run fault detection with existing models
3. Check Results tab for output files and visualizations
        """

        info_box = tk.Text(root, height=15, wrap="word", font=("TkDefaultFont", 9))
        info_box.pack(fill="both", expand=True, **pad)
        info_box.insert("1.0", info_text)
        info_box.configure(state="disabled", bg="#f0f0f0")

    # ------------------------------ Training ------------------------------
    def _build_training(self, root):
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root, text="DRL Training Configuration",
                  font=("TkDefaultFont", 12, "bold")).pack(anchor="w", **pad)

        # Algorithm selection
        algo_frame = ttk.LabelFrame(root, text="Algorithm", padding=10)
        algo_frame.pack(fill="x", **pad)

        ttk.Radiobutton(algo_frame, text="PPO (Proximal Policy Optimization)",
                        variable=self.algorithm_var, value="PPO").pack(anchor="w", pady=2)
        ttk.Radiobutton(algo_frame, text="TD3HD (Twin Delayed DDPG with Hindsight)",
                        variable=self.algorithm_var, value="TDHD").pack(anchor="w", pady=2)
        ttk.Radiobutton(algo_frame, text="DQN (Deep Q-Network)",
                        variable=self.algorithm_var, value="DQN").pack(anchor="w", pady=2)

        # Training parameters
        param_frame = ttk.LabelFrame(root, text="Training Parameters", padding=10)
        param_frame.pack(fill="x", **pad)

        self.iterations_var = tk.IntVar(value=100)
        self.satellites_var = tk.IntVar(value=4)
        self.targets_var = tk.IntVar(value=40)

        row1 = ttk.Frame(param_frame); row1.pack(fill="x", pady=2)
        ttk.Label(row1, text="Training Iterations:", width=20).pack(side="left")
        ttk.Entry(row1, textvariable=self.iterations_var, width=10).pack(side="left", padx=5)

        row2 = ttk.Frame(param_frame); row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Number of Satellites:", width=20).pack(side="left")
        ttk.Entry(row2, textvariable=self.satellites_var, width=10).pack(side="left", padx=5)

        row3 = ttk.Frame(param_frame); row3.pack(fill="x", pady=2)
        ttk.Label(row3, text="Number of Targets:", width=20).pack(side="left")
        ttk.Entry(row3, textvariable=self.targets_var, width=10).pack(side="left", padx=5)

        # Control buttons
        ctrl_frame = ttk.Frame(root); ctrl_frame.pack(fill="x", **pad)
        self.train_btn = ttk.Button(ctrl_frame, text="Start Training", command=self._start_training)
        self.train_btn.pack(side="left", padx=(0, 5))
        ttk.Button(ctrl_frame, text="Stop", command=self._stop_training).pack(side="left", padx=(0, 10))
        ttk.Label(ctrl_frame, textvariable=self.status_text).pack(side="left")

        # Log
        log_frame = ttk.Frame(root); log_frame.pack(fill="both", expand=True, **pad)
        ttk.Label(log_frame, text="Training Log:").pack(anchor="w")

        text_frame = ttk.Frame(log_frame); text_frame.pack(fill="both", expand=True, pady=(2, 0))
        self.train_log = tk.Text(text_frame, height=12, state="disabled", wrap="word")
        scrollbar = ttk.Scrollbar(text_frame, command=self.train_log.yview)
        self.train_log.configure(yscrollcommand=scrollbar.set)
        self.train_log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    # ------------------------------ Integration ------------------------------
    def _build_integration(self, root):
        pad = {"padx": 10, "pady": 6}

        ttk.Label(root,
                  text="Integrated Fault Detection + DRL Task Reassignment",
                  font=("TkDefaultFont", 12, "bold")).pack(anchor="w", **pad)

        config_frame = ttk.LabelFrame(root, text="Configuration", padding=10)
        config_frame.pack(fill="x", **pad)

        ttk.Checkbutton(
            config_frame,
            text="Enable Fault Detection",
            variable=tk.BooleanVar(value=True),
            state="disabled"
        ).pack(anchor="w")

        ttk.Checkbutton(
            config_frame,
            text="Enable DRL Task Reassignment",
            variable=self.integration_mode
        ).pack(anchor="w")

        fault_frame = ttk.LabelFrame(root, text="Simulation Configuration", padding=10)
        fault_frame.pack(fill="x", **pad)

        ttk.Label(
            fault_frame,
            text="Faults and configuration will be loaded from the main GUI tabs.",
            wraplength=600
        ).pack(anchor="w")
        ttk.Button(
            fault_frame,
            text="Go to Fault Configuration Tab",
            command=self._switch_to_fault_tab
        ).pack(anchor="w", pady=5)

        # -------- Model Management Section --------
        if MODEL_MANAGER_AVAILABLE:
            model_frame = ttk.LabelFrame(root, text="Available DRL Models", padding=10)
            model_frame.pack(fill="x", **pad)

            try:
                manager = get_model_manager()
                models = manager.list_models()

                if models:
                    canvas = tk.Canvas(model_frame, height=160)
                    scrollbar = ttk.Scrollbar(model_frame, orient="vertical", command=canvas.yview)
                    scrollable_frame = ttk.Frame(canvas)

                    def _on_configure(_):
                        canvas.configure(scrollregion=canvas.bbox("all"))
                    scrollable_frame.bind("<Configure>", _on_configure)

                    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
                    canvas.configure(yscrollcommand=scrollbar.set)

                    hdr = ttk.Frame(scrollable_frame); hdr.pack(fill="x", pady=(0, 4))
                    ttk.Label(hdr, text="", width=3).pack(side="left")
                    ttk.Label(hdr, text="Name", width=28, font=("TkDefaultFont", 9, "bold")).pack(side="left")
                    ttk.Label(hdr, text="Purpose", width=36, font=("TkDefaultFont", 9, "bold")).pack(side="left")
                    ttk.Label(hdr, text="Type", width=12, font=("TkDefaultFont", 9, "bold")).pack(side="left")

                    for model in models:
                        exists = os.path.exists(model.get("path", ""))
                        status_icon = "✓" if exists else "✗"
                        status_color = "green" if exists else "red"

                        row = ttk.Frame(scrollable_frame); row.pack(fill="x", pady=1)
                        ttk.Label(row, text=status_icon, foreground=status_color, width=3).pack(side="left")
                        ttk.Label(row, text=model.get("name", "—"), width=28).pack(side="left")
                        ttk.Label(row, text=model.get("purpose", "—"), width=36, foreground="gray").pack(side="left")
                        ttk.Label(row, text=f"[{model.get('type','—')}]", width=12, foreground="blue").pack(side="left")

                    canvas.pack(side="left", fill="both", expand=True)
                    scrollbar.pack(side="right", fill="y")

                    ttk.Button(model_frame, text="Refresh Models", command=self._refresh_model_list).pack(pady=6)
                else:
                    ttk.Label(model_frame, text="No models found. Train a model first.", foreground="orange").pack(pady=10)
            except Exception as e:
                ttk.Label(model_frame, text=f"Error loading models: {e}", foreground="red").pack(pady=10)
        else:
            warning_frame = ttk.LabelFrame(root, text="Model Manager Not Available", padding=10)
            warning_frame.pack(fill="x", **pad)
            ttk.Label(
                warning_frame,
                text=("⚠️ Model manager is not available.\n"
                      "Make sure drl_model_manager.py is in the core/ directory."),
                foreground="orange"
            ).pack(pady=10)

            # Buttons must attach to warning_frame (fixed: was model_frame)
            btns = ttk.Frame(warning_frame); btns.pack(pady=(6, 0))
            def _load_latest(prefix):
                try:
                    from drl_model_manager import get_model_manager
                    m = get_model_manager()
                    model = m.load_latest_by_prefix(prefix)
                    messagebox.showinfo("DRL", f"Loaded latest {prefix} model OK.")
                except Exception as e:
                    messagebox.showerror("DRL", f"Load failed: {e}")

            ttk.Button(btns, text="Load Latest PPO",  command=lambda: _load_latest("PPO")).pack(side="left", padx=4)
            ttk.Button(btns, text="Load Latest TDHD", command=lambda: _load_latest("TDHD")).pack(side="left", padx=4)
            ttk.Button(btns, text="Load Latest DQN",  command=lambda: _load_latest("DQN")).pack(side="left", padx=4)

        # Run integration
        run_frame = ttk.Frame(root); run_frame.pack(fill="x", **pad)
        self.integrate_btn = ttk.Button(run_frame, text="Run Integrated System", command=self._run_integration)
        self.integrate_btn.pack(side="left", padx=(0, 5))
        ttk.Button(run_frame, text="Load Previous Results", command=self._load_integration_results).pack(side="left")

        # Results summary
        summary_frame = ttk.LabelFrame(root, text="Integration Results", padding=5)
        summary_frame.pack(fill="both", expand=True, **pad)

        text_frame = ttk.Frame(summary_frame); text_frame.pack(fill="both", expand=True)
        self.integration_text = tk.Text(text_frame, height=10, state="disabled", wrap="word")
        int_scrollbar = ttk.Scrollbar(text_frame, command=self.integration_text.yview)
        self.integration_text.configure(yscrollcommand=int_scrollbar.set)
        self.integration_text.pack(side="left", fill="both", expand=True)
        int_scrollbar.pack(side="right", fill="y")

    # ------------------------------ Results ------------------------------
    def _build_results(self, root):
        pad = {"padx": 10, "pady": 6}

        hdr = ttk.Frame(root); hdr.pack(fill="x", **pad)
        ttk.Label(hdr, text="Results Folder:").pack(side="left")

        self.res_entry = ttk.Entry(hdr, textvariable=self.results_dir)
        self.res_entry.pack(side="left", fill="x", expand=True, padx=(6, 6))

        ttk.Button(hdr, text="Use Default", command=lambda: self._set_results_dir(RESULT_DIR)).pack(side="left", padx=(0, 5))
        ttk.Button(hdr, text="Refresh", command=self._refresh_results).pack(side="left", padx=(0, 5))
        ttk.Button(hdr, text="Open Folder", command=self._open_results_folder).pack(side="left")

        cols = ("name", "size", "modified")
        tvf = ttk.Frame(root); tvf.pack(fill="both", expand=True, **pad)

        self.tree = ttk.Treeview(tvf, columns=cols, show="headings", selectmode="browse", height=10)
        for c, w in zip(cols, (360, 100, 180)):
            self.tree.heading(c, text=c.title())
            self.tree.column(c, width=w, anchor="w")

        tree_scroll = ttk.Scrollbar(tvf, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        self.tree.pack(side="left", fill="both", expand=True)
        tree_scroll.pack(side="right", fill="y")

        self.preview_frame = ttk.LabelFrame(root, text="Preview")
        self.preview_frame.pack(fill="both", expand=True, **pad)

        if MPL_AVAILABLE and Figure and FigureCanvasTkAgg:
            self.fig = Figure(figsize=(5, 2.6), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
        else:
            self.fig = None
            ttk.Label(self.preview_frame, text="Matplotlib not available - preview disabled").pack(pady=20)

        self.tree.bind("<<TreeviewSelect>>", self._on_select_artifact)

        # Initial refresh if folder exists
        if self.results_dir.get() and os.path.isdir(self.results_dir.get()):
            self._refresh_results()

    # ------------------------------ Training Actions ------------------------------
    def _start_training(self):
        """Start DRL training"""
        algo = self.algorithm_var.get()

        # Use the scripts that exist in your DRL folder
        script_map = {
            "PPO": "PPO.py",
            "TDHD": "TDHDYear2.py",
            "DQN": "DQNYear2.py",
        }

        script_name = script_map.get(algo)
        if not script_name:
            messagebox.showerror("Error", f"Unknown algorithm: {algo}")
            return

        script_path = os.path.join(DRL_DIR, script_name)
        if not os.path.exists(script_path):
            messagebox.showerror(
                "Error",
                f"Script not found:\n{script_path}\n\nPlease ensure DRL scripts are in the DRL/ folder.",
            )
            return

        self._append_train_log(f"Starting {algo} training...\n")
        self._append_train_log(f"Script: {script_path}\n")
        self._append_train_log(f"Iterations: {self.iterations_var.get()}\n")
        self._append_train_log(f"Satellites: {self.satellites_var.get()}\n")
        self._append_train_log(f"Targets: {self.targets_var.get()}\n\n")

        self._run_script(script_path)

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

    # ------------------------------ Integration Actions ------------------------------
    def _run_integration(self):
        """Run integrated fault detection + DRL task reassignment"""
        self._append_integration_log("=" * 60 + "\n")
        self._append_integration_log("INTEGRATED FAULT DETECTION + DRL TASK REASSIGNMENT\n")
        self._append_integration_log("=" * 60 + "\n\n")
        self._append_integration_log("Note: This feature requires a simulation to be run first.\n")
        self._append_integration_log("Please run a simulation from the main GUI with faults enabled.\n\n")

        messagebox.showinfo(
            "DRL Integration",
            "DRL Integration is triggered automatically when you:\n\n"
            "1. Configure faults in the Fault Configuration tab\n"
            "2. Run a simulation from the main GUI\n"
            "3. Results will appear in the Task Reassignment tab\n\n"
            "This tab is for standalone DRL training and result browsing."
        )

    def _integration_worker(self):
        """Worker thread for integration - placeholder"""
        pass

    def _get_simulation_config(self):
        """Get simulation configuration from parent GUI"""
        if not self.parent_app:
            return None

        config = {
            "satellites": getattr(self.parent_app, 'satellites', []),
            "targets": getattr(self.parent_app, 'targets', []),
            "duration": getattr(self.parent_app, 'simulation_time', tk.DoubleVar(value=30)).get()
        }
        return config if config["satellites"] else None

    def _switch_to_fault_tab(self):
        """Switch to fault configuration tab in main GUI"""
        if self.parent_app and hasattr(self.parent_app, 'notebook'):
            for i in range(self.parent_app.notebook.index("end")):
                if self.parent_app.notebook.tab(i, "text") == "Fault Configuration":
                    self.parent_app.notebook.select(i)
                    break

    def _load_integration_results(self):
        """Load previous integration results"""
        file_path = filedialog.askopenfilename(
            title="Load Integration Results",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not file_path:
            return

        try:
            with open(file_path, 'r') as f:
                results = json.load(f)

            self._integration_text_clear()
            self._append_integration_log("Loaded results from: " + os.path.basename(file_path) + "\n\n")
            self._append_integration_log(json.dumps(results, indent=2))

            messagebox.showinfo("Success", "Results loaded successfully")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load results:\n{e}")

    def _refresh_model_list(self):
        """Refresh the model list display in the Integration tab."""
        if not MODEL_MANAGER_AVAILABLE:
            messagebox.showwarning("Model Manager", "Model manager is not available.")
            return
        try:
            manager = get_model_manager()
            if hasattr(manager, "_scan_for_models"):
                manager._scan_for_models()
            messagebox.showinfo("Success", "Model list refreshed!")
            # Rebuild just the Integration tab to reflect changes
            for child in self.integration.winfo_children():
                child.destroy()
            self._build_integration(self.integration)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to refresh models:\n{e}")

    # ------------------------------ Script Execution ------------------------------
    def _run_script(self, script_path):
        """Run a Python script in a subprocess"""
        python = sys.executable
        cmd = [python, script_path]

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
                universal_newlines=True
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
        """Drain the output queue"""
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

    # ------------------------------ Logging Helpers ------------------------------
    def _append_train_log(self, s: str):
        self.train_log.configure(state="normal")
        self.train_log.insert("end", s)
        self._log_lines += s.count("\n")
        if self._log_lines > 3000:
            self.train_log.delete("1.0", "500.0")
            self._log_lines -= 500
        self.train_log.see("end")
        self.train_log.configure(state="disabled")

    def _append_integration_log(self, s: str):
        self.integration_text.configure(state="normal")
        self.integration_text.insert("end", s)
        self.integration_text.see("end")
        self.integration_text.configure(state="disabled")

    def _integration_text_clear(self):
        self.integration_text.configure(state="normal")
        self.integration_text.delete("1.0", "end")
        self.integration_text.configure(state="disabled")

    # ------------------------------ Results Handling ------------------------------
    def _ensure_results_dir(self, folder: str):
        try:
            if folder and not os.path.isdir(folder):
                os.makedirs(folder, exist_ok=True)
        except Exception as e:
            print(f"Could not create results dir '{folder}': {e}")

    def _set_results_dir(self, folder: str):
        if not folder:
            return
        self._ensure_results_dir(folder)
        self.results_dir.set(folder)
        self._refresh_results()

    def _refresh_results(self):
        folder = self.results_dir.get().strip()
        self.tree.delete(*self.tree.get_children())
        if not folder or not os.path.isdir(folder):
            return

        allow = {".xlsx", ".csv", ".png", ".jpg", ".json", ".txt", ".pth", ".pt", ".zip", ".h5", ".keras"}
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
        sel = self.tree.selection()
        if not sel:
            return
        name = self.tree.item(sel[0], "values")[0]
        path = os.path.join(self.results_dir.get().strip(), name)

        if name.lower().endswith(('.png', '.jpg', '.jpeg')):
            self._preview_image(path)
        elif name.lower().endswith('.json'):
            self._preview_json(path)
        else:
            # Clear preview text/axes for unsupported types
            if self.fig is not None:
                self.ax.clear()
                self.ax.text(0.5, 0.5, f"{name}\n(no preview)", ha="center", va="center", fontsize=10)
                self.ax.set_xticks([]); self.ax.set_yticks([])
                self.canvas.draw_idle()

    def _preview_image(self, path):
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
        if self.fig is None:
            return
        try:
            with open(path, 'r') as f:
                data = json.load(f)
            self.ax.clear()
            self.ax.text(0.5, 0.5, f"JSON File\n{os.path.basename(path)}\n{len(str(data))} characters",
                         ha="center", va="center", fontsize=10)
            self.ax.set_title("JSON Preview")
            self.ax.set_xticks([]); self.ax.set_yticks([])
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error previewing JSON: {e}")

    def _open_results_folder(self):
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
