#!/usr/bin/env python3
"""
drl_tab.py — Fully Integrated DRL Tab

Integrates with:
- GUI system (spacecraft_simulator_gui.py)
- DRL algorithms (PPO, DQN from DRL/ folder)
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

# Add paths for integration
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
        self.results_dir = tk.StringVar(value=default_results_dir or os.path.join(drl_dir, "results"))
        rd = os.path.normpath(self.results_dir.get())
        if rd == os.path.normpath(drl_dir):
            rd = os.path.join(drl_dir, "results")
            self.results_dir.set(rd)

        os.makedirs(self.results_dir.get(), exist_ok=True)  # ensure it exists
        self.status_text = tk.StringVar(value="Idle")
        self.algorithm_var = tk.StringVar(value="PPO")
        self.training_mode = tk.BooleanVar(value=False)
        self._log_lines = 0


        # Layout: nested notebook
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=5, pady=5)
        self._nb = nb

        self.overview = ttk.Frame(nb)
        self.training = ttk.Frame(nb)
        self.results = ttk.Frame(nb)
        
        nb.add(self.overview, text="Overview")
        nb.add(self.training, text="Training")
        nb.add(self.results, text="Results")

        self._build_overview(self.overview)
        self._build_training(self.training)
        self._build_results(self.results)

        # Poll the queue for subprocess output
        self.after(150, self._drain_queue)

    # ------------------------------ Overview ------------------------------
    def _build_overview(self, root):
        pad = {"padx": 10, "pady": 6}

        # Header
        ttk.Label(root, text="Deep Reinforcement Learning System",
                  font=("TkDefaultFont", 14, "bold")).pack(anchor="w", **pad)
        
        # Status frame
        status_frame = ttk.LabelFrame(root, text="System Status", padding=10)
        status_frame.pack(fill="x", **pad)
        
        # Component availability
        components = [
            ("Matplotlib", MPL_AVAILABLE),
            ("Pandas", PD_AVAILABLE)
        ]
        
        for name, available in components:
            row = ttk.Frame(status_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"{name}:", width=20).pack(side="left")
            status = "✓ Available" if available else "✗ Not Available"
            color = "green" if available else "red"
            lbl = ttk.Label(row, text=status, foreground=color)
            lbl.pack(side="left")
        
        # Separator
        ttk.Separator(root, orient="horizontal").pack(fill="x", pady=10)
        
        # Info text
        info_text = """DRL System Components:

1. Training: Train PPO, or DQN algorithms for task reassignment
2. Results: View training metrics, fault detection outputs, and reassignment plans

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
        ttk.Radiobutton(algo_frame, text="DQN (Deep Q-Network)", 
                       variable=self.algorithm_var, value="DQN").pack(anchor="w", pady=2)

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
        
        # Text widget with scrollbar
        text_frame = ttk.Frame(log_frame)
        text_frame.pack(fill="both", expand=True, pady=(2, 0))
        
        self.train_log = tk.Text(text_frame, height=12, state="disabled", wrap="word")
        scrollbar = ttk.Scrollbar(text_frame, command=self.train_log.yview)
        self.train_log.configure(yscrollcommand=scrollbar.set)
        
        self.train_log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

   # ------------------------------ Results ------------------------------
    def _build_results(self, root):
        pad = {"padx": 10, "pady": 6}

        # Header
        hdr = ttk.Frame(root)
        hdr.pack(fill="x", **pad)
        
        ttk.Label(hdr, text="Results Folder:").pack(side="left")
        self.res_entry = ttk.Entry(hdr, textvariable=self.results_dir)
        self.res_entry.pack(side="left", fill="x", expand=True, padx=(6,6))
        ttk.Button(hdr, text="Refresh", command=self._refresh_results).pack(side="left", padx=(0, 5))
        ttk.Button(hdr, text="Open Folder", command=self._open_results_folder).pack(side="left")

        # Treeview for files
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

        # Preview area
        self.preview_frame = ttk.LabelFrame(root, text="Preview")
        self.preview_frame.pack(fill="both", expand=True, **pad)

        if MPL_AVAILABLE and Figure and FigureCanvasTkAgg:
            self.fig = Figure(figsize=(20, 10.2), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
            self.canvas.get_tk_widget().pack(fill="x", expand=True)
            
        else:
            self.fig = None
            ttk.Label(self.preview_frame, text="Matplotlib not available - preview disabled").pack(pady=20)

        self.tree.bind("<<TreeviewSelect>>", self._on_select_artifact)
        
        # Auto-refresh if results dir exists
        if self.results_dir.get() and os.path.isdir(self.results_dir.get()):
            self._refresh_results()

    # ------------------------------ Training Actions ------------------------------
    def _start_training(self):
        """Start DRL training"""
        algo = self.algorithm_var.get()
        
        # Map algorithm to script
        script_map = {
            "PPO": "PPO.py",
            "DQN": "DQN.py"
        }
        
        script_name = script_map.get(algo)
        if not script_name:
            messagebox.showerror("Error", f"Unknown algorithm: {algo}")
            return
        
        script_path = os.path.join(drl_dir, script_name)
        if not os.path.exists(script_path):
            messagebox.showerror("Error", f"Script not found:\n{script_path}\n\nPlease ensure DRL scripts are in the DRL/ folder.")
            return

        self._append_train_log(f"Starting {algo} training...\n")
        self._append_train_log(f"Script: {script_path}\n")
        
        # Run the script
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
        """Append to training log"""
        self.train_log.configure(state="normal")
        self.train_log.insert("end", s)
        self._log_lines += s.count("\n")
        if self._log_lines > 3000:
            self.train_log.delete("1.0", "500.0")
            self._log_lines -= 500
        self.train_log.see("end")
        self.train_log.configure(state="disabled")

    # ------------------------------ Results Handling ------------------------------
    def _refresh_results(self):
        """Refresh results file list"""
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
        """Preview an image file"""
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
        """Preview JSON file"""
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
        """Open results folder in file explorer"""
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