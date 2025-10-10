#!/usr/bin/env python3
"""
drl_tab.py — a plug-in tab for Tkinter GUIs that adds:
- "DRL" tab with a nested Notebook: "Overview" and "Results"
- Lets the user pick a DRL script (e.g., "PPO2 3.py"), run it in a background
  thread/subprocess, stream logs, and browse result artifacts (xlsx/csv/png).
No external GUI framework required; pure Tkinter + ttk.
Matplotlib is optional for previewing plots if available.
"""

import os
import sys
import threading
import queue
import subprocess
import time
import datetime
import platform
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Optional imports (safe fallback if not present)
try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib
    matplotlib.use("Agg")  # headless-safe
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
except Exception:
    Figure = None
    FigureCanvasTkAgg = None


class DRLTab(ttk.Frame):
    """
    Drop-in widget:
        drl_tab = DRLTab(parent=self.notebook, parent_app=self)
        self.notebook.add(drl_tab, text="DRL")
    Expectation: parent_app may optionally provide add_log(msg: str).
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
        self.results_dir = tk.StringVar(value=default_results_dir or "")
        self.status_text = tk.StringVar(value="Idle")
        self._log_lines = 0

        # Layout: nested notebook
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)
        self._nb = nb

        self.overview = ttk.Frame(nb)
        self.results = ttk.Frame(nb)
        nb.add(self.overview, text="Overview")
        nb.add(self.results, text="Results")

        self._build_overview(self.overview)
        self._build_results(self.results)

        # Poll the queue for subprocess output
        self.after(150, self._drain_queue)

    # ------------------------------ Overview ------------------------------
    def _build_overview(self, root):
        pad = {"padx": 10, "pady": 6}

        # Header
        ttk.Label(root, text="Deep Reinforcement Learning — Runner",
                  font=("TkDefaultFont", 12, "bold")).pack(anchor="w", **pad)

        # Script picker
        f1 = ttk.Frame(root); f1.pack(fill="x", **pad)
        ttk.Label(f1, text="DRL script:").pack(side="left")
        e1 = ttk.Entry(f1, textvariable=self.script_path)
        e1.pack(side="left", fill="x", expand=True, padx=(6,6))
        ttk.Button(f1, text="Browse", command=self._browse_script).pack(side="left")

        # Results folder picker
        f2 = ttk.Frame(root); f2.pack(fill="x", **pad)
        ttk.Label(f2, text="Results folder:").pack(side="left")
        e2 = ttk.Entry(f2, textvariable=self.results_dir)
        e2.pack(side="left", fill="x", expand=True, padx=(6,6))
        ttk.Button(f2, text="Choose", command=self._choose_results_dir).pack(side="left")

        # Controls
        f3 = ttk.Frame(root); f3.pack(fill="x", **pad)
        self.run_btn = ttk.Button(f3, text="Run DRL", command=self._run_drl)
        self.run_btn.pack(side="left")
        ttk.Button(f3, text="Stop", command=self._stop_drl).pack(side="left", padx=(6,0))
        ttk.Label(f3, textvariable=self.status_text).pack(side="right")

        # Log box
        ttk.Label(root, text="Run log:").pack(anchor="w", padx=10, pady=(12, 0))
        self.log = tk.Text(root, height=18, state="disabled")
        self.log.pack(fill="both", expand=True, padx=10, pady=(2, 10))

        # Helper text
        helper = (
            "Tip: This tab simply runs your Python script in a background process and streams stdout/stderr.\n"
            "If your script writes .xlsx/.csv/.png to the results folder, use the Results tab to browse/open them."
        )
        ttk.Label(root, text=helper, foreground="#555").pack(anchor="w", padx=10, pady=(0,10))

    # ------------------------------ Results ------------------------------
    def _build_results(self, root):
        pad = {"padx": 10, "pady": 6}

        hdr = ttk.Frame(root); hdr.pack(fill="x", **pad)
        ttk.Label(hdr, text="Artifacts (xlsx/csv/png) in folder:").pack(side="left")
        self.res_entry = ttk.Entry(hdr, textvariable=self.results_dir)
        self.res_entry.pack(side="left", fill="x", expand=True, padx=(6,6))
        ttk.Button(hdr, text="Refresh", command=self._refresh_results).pack(side="left")
        ttk.Button(hdr, text="Open Folder", command=self._open_results_folder).pack(side="left", padx=(6,0))

        # Treeview for files
        cols = ("name","size","modified")
        tvf = ttk.Frame(root); tvf.pack(fill="both", expand=True, **pad)
        self.tree = ttk.Treeview(tvf, columns=cols, show="headings", selectmode="browse", height=12)
        for c, w in zip(cols, (360, 100, 180)):
            self.tree.heading(c, text=c.title())
            self.tree.column(c, width=w, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(tvf, orient="vertical", command=self.tree.yview)
        sb.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=sb.set)

        # Preview area (optional matplotlib)
        self.preview_frame = ttk.LabelFrame(root, text="Preview")
        self.preview_frame.pack(fill="both", expand=True, **pad)

        if Figure and FigureCanvasTkAgg:
            self.fig = Figure(figsize=(5, 2.6), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
        else:
            self.fig = None
            ttk.Label(self.preview_frame, text="Matplotlib not available; preview disabled.").pack(anchor="w", padx=10, pady=10)

        # Bind selection
        self.tree.bind("<<TreeviewSelect>>", self._on_select_artifact)

    # ------------------------------ Actions ------------------------------
    def _browse_script(self):
        path = filedialog.askopenfilename(
            title="Select DRL script",
            filetypes=[("Python files","*.py"), ("All files","*.*")]
        )
        if path:
            self.script_path.set(path)

    def _choose_results_dir(self):
        path = filedialog.askdirectory(title="Select results folder")
        if path:
            self.results_dir.set(path)
            self._refresh_results()

    def _run_drl(self):
        script = self.script_path.get().strip()
        if not script:
            messagebox.showerror("DRL", "Please choose a DRL script to run.")
            return
        if not os.path.exists(script):
            messagebox.showerror("DRL", f"Script not found:\n{script}")
            return

        # Use script's folder as working dir, but write artifacts to chosen results dir if env var supported
        env = os.environ.copy()
        if self.results_dir.get().strip():
            env["DRL_RESULTS_DIR"] = self.results_dir.get().strip()

        python = sys.executable
        cmd = [python, script]

        try:
            self._append_log(f"$ {' '.join(cmd)}\n")
            self.status_text.set("Running…")
            self.run_btn.configure(state="disabled")
            self._stop_reader.clear()

            self.proc = subprocess.Popen(
                cmd,
                cwd=os.path.dirname(script) or None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True,
                env=env
            )
        except Exception as e:
            messagebox.showerror("DRL", f"Failed to start:\n{e}")
            self.status_text.set("Failed")
            self.run_btn.configure(state="normal")
            return

        # Start reader thread
        self.worker_thread = threading.Thread(target=self._reader_worker, daemon=True)
        self.worker_thread.start()

    def _stop_drl(self):
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
                self._append_log("[DRL] Terminated process.\n")
            except Exception as e:
                self._append_log(f"[DRL] Terminate failed: {e}\n")
        self.status_text.set("Stopping…")
        self._stop_reader.set()

    def _reader_worker(self):
        """Read subprocess stdout lines and enqueue for the Tk thread"""
        try:
            for line in self.proc.stdout:
                if self._stop_reader.is_set():
                    break
                self.queue.put(line)
        except Exception as e:
            self.queue.put(f"[reader] error: {e}\n")
        finally:
            # Drain any remainder and signal close
            rc = None
            try:
                rc = self.proc.wait(timeout=1.0)
            except Exception:
                pass
            self.queue.put(f"\n[DRL] Finished with return code: {rc}\n")
            self.queue.put(None)  # sentinel

    def _drain_queue(self):
        """Periodic Tk callback to drain the queue and write to the log Text widget"""
        drained = False
        while True:
            try:
                item = self.queue.get_nowait()
            except queue.Empty:
                break
            drained = True
            if item is None:
                # Process finished
                self.status_text.set("Idle")
                self.run_btn.configure(state="normal")
                # Refresh after run
                try:
                    self._refresh_results()
                except Exception:
                    pass
                continue
            # Append line
            self._append_log(item)

        # Reschedule
        self.after(150 if drained else 250, self._drain_queue)

    def _append_log(self, s: str):
        self.log.configure(state="normal")
        self.log.insert("end", s)
        self._log_lines += s.count("\n")
        # cap lines to avoid unbounded growth
        if self._log_lines > 3000:
            self.log.delete("1.0", "500.0")  # drop first ~500 lines
            self._log_lines -= 500
        self.log.see("end")
        self.log.configure(state="disabled")
        # bubble up
        if getattr(self.parent_app, "add_log", None):
            try:
                self.parent_app.add_log(s.rstrip("\n"))
            except Exception:
                pass

    # ------------------------------ Results handling ------------------------------
    def _refresh_results(self):
        folder = self.results_dir.get().strip()
        self.tree.delete(*self.tree.get_children())
        if not folder or not os.path.isdir(folder):
            return
        allow = {".xlsx",".csv",".png",".jpg",".jpeg",".pdf",".txt"}
        rows = []
        for name in os.listdir(folder):
            ext = os.path.splitext(name)[1].lower()
            if ext in allow:
                p = os.path.join(folder, name)
                try:
                    stat = os.stat(p)
                    size_kb = int(stat.st_size/1024)
                    mtime = datetime.datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
                    rows.append((name, size_kb, mtime))
                except Exception:
                    rows.append((name, "-", "-"))
        rows.sort(key=lambda r: r[0].lower())
        for r in rows:
            self.tree.insert("", "end", values=r)

    def _on_select_artifact(self, event=None):
        sel = self.tree.selection()
        if not sel:
            return
        name, _, _ = self.tree.item(sel[0], "values")
        path = os.path.join(self.results_dir.get().strip(), name)
        self._preview_artifact(path)

    def _open_results_folder(self):
        folder = self.results_dir.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showinfo("DRL", "No folder selected or folder not found.")
            return
        self._open_path(folder)

    def _open_path(self, path):
        try:
            if platform.system() == "Windows":
                os.startfile(path)  # type: ignore[attr-defined]
            elif platform.system() == "Darwin":
                subprocess.call(["open", path])
            else:
                subprocess.call(["xdg-open", path])
        except Exception as e:
            messagebox.showerror("DRL", f"Open failed:\n{e}")

    def _preview_artifact(self, path):
        if not os.path.exists(path):
            return
        if self.fig is None:
            return  # matplotlib not available
        # Basic CSV preview: look for columns with 'reward' or a single numeric series
        import csv
        data = []
        try:
            if path.lower().endswith(".csv"):
                with open(path, "r", newline="") as f:
                    reader = csv.DictReader(f)
                    cols = reader.fieldnames or []
                    reward_col = None
                    for c in cols:
                        if "reward" in c.lower():
                            reward_col = c
                            break
                    for row in reader:
                        if reward_col and row.get(reward_col):
                            try:
                                data.append(float(row[reward_col]))
                            except Exception:
                                pass
            # Fallback: try to parse a simple one-number-per-line text file
            if not data and path.lower().endswith((".txt",".csv")):
                with open(path, "r") as f:
                    for line in f:
                        try:
                            data.append(float(line.strip()))
                        except Exception:
                            pass
        except Exception:
            data = []

        self.ax.clear()
        if data:
            self.ax.plot(range(len(data)), data)
            self.ax.set_title(f"Preview: {os.path.basename(path)}")
            self.ax.set_xlabel("Step / Episode")
            self.ax.set_ylabel("Reward")
        else:
            self.ax.text(0.5, 0.5, "No preview available", ha="center", va="center")
            self.ax.set_title(f"Preview: {os.path.basename(path)}")
            self.ax.set_xticks([]); self.ax.set_yticks([])
        self.canvas.draw_idle()
