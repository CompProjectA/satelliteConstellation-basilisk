#!/usr/bin/env python
"""
fault_tab.py - FIXED VERSION (hardened)
Enhanced Fault Configuration tab with proper scrolling and cluster support
"""
import tkinter as tk
from tkinter import ttk, messagebox
from .base_tab import BaseTab

# Visualization for cluster faults
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np


class FaultTab(BaseTab):
    """Fixed Fault Configuration tab with scrolling and cluster support"""

    def __init__(self, parent_app, parent_frame):
        super().__init__(parent_app, parent_frame)

        # Data handles
        self.satellites = parent_app.satellites

        # UI
        self.create_tab_ui()

        # React to fault type changes (Tk 8.6+ uses trace_add)
        try:
            self.fault_type_var.trace_add('write', self.on_fault_type_change)
        except Exception:
            self.fault_type_var.trace('w', self.on_fault_type_change)

        # Safe initial load
        self.update_satellite_dropdown()
        if self.satellites:
            self.load_fault_config(0)
        self.update_fault_summary()

    # -----------------------------
    # Selection helpers
    # -----------------------------
    def get_active_satellite_index(self):
        name = self.fault_satellite_var.get()
        for i, sat in enumerate(self.satellites):
            if sat["name"] == name:
                return i
        return None

    def set_active_satellite(self, index):
        if 0 <= index < len(self.satellites):
            self.fault_satellite_combo.current(index)
            self.load_fault_config(index)

    # -----------------------------
    # UI build
    # -----------------------------
    def create_tab_ui(self):
        self.fault_notebook = ttk.Notebook(self.parent_frame)
        self.fault_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        config_frame = ttk.Frame(self.fault_notebook)
        summary_frame = ttk.Frame(self.fault_notebook)
        cluster_frame = ttk.Frame(self.fault_notebook)

        self.fault_notebook.add(config_frame, text="Fault Configuration")
        self.fault_notebook.add(summary_frame, text="Fault Summary")
        self.fault_notebook.add(cluster_frame, text="Cluster Faults")

        self._create_fault_config_tab(config_frame)
        self._create_fault_summary_tab(summary_frame)
        self._create_cluster_fault_tab(cluster_frame)

    def _create_fault_config_tab(self, parent):
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True)

        # Scroll setup
        # Scroll setup
        canvas = tk.Canvas(main_container, highlightthickness=0)
        scrollbar = ttk.Scrollbar(main_container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        scrollable_frame = ttk.Frame(canvas)

        # 1) Keep a handle to the window item
        self._scroll_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

        # 2) Update scrollregion when inner frame changes
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        # 3) CRUCIAL: stretch inner frame to canvas width
        def _resize_inner(event):
            canvas.itemconfigure(self._scroll_window, width=event.width)

        canvas.bind("<Configure>", _resize_inner)

        # Mousewheel
        def _on_wheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", _on_wheel))
        canvas.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))


        # ----- Satellite Selection -----
        select_frame = ttk.LabelFrame(scrollable_frame, text="Satellite Selection", padding=10)
        select_frame.pack(fill=tk.X, padx=10, pady=10)

        row1 = ttk.Frame(select_frame)
        row1.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(row1, text="Select Cluster:").pack(side=tk.LEFT)
        self.cluster_var = tk.StringVar()
        self.cluster_combo = ttk.Combobox(row1, textvariable=self.cluster_var, width=20, state="readonly")
        self.cluster_combo.pack(side=tk.LEFT, padx=5)
        self.cluster_combo.bind('<<ComboboxSelected>>', self.on_cluster_selected)
        ttk.Button(row1, text="Show All", command=self.show_all_satellites).pack(side=tk.LEFT, padx=10)

        row2 = ttk.Frame(select_frame)
        row2.pack(fill=tk.X)
        ttk.Label(row2, text="Select Satellite:").pack(side=tk.LEFT)
        self.fault_satellite_var = tk.StringVar()
        self.fault_satellite_combo = ttk.Combobox(row2, textvariable=self.fault_satellite_var, width=20, state="readonly")
        self.fault_satellite_combo.pack(side=tk.LEFT, padx=5)
        self.fault_satellite_combo.bind('<<ComboboxSelected>>', self.on_fault_satellite_changed)

        self.fault_status_label = ttk.Label(row2, text="No fault enabled")
        self.fault_status_label.pack(side=tk.RIGHT, padx=5)

        self.sat_info_label = ttk.Label(select_frame, text="", font=('TkDefaultFont', 9, 'italic'))
        self.sat_info_label.pack(anchor=tk.W, pady=2)

        # ----- Fault Config -----
        fault_config_frame = ttk.LabelFrame(scrollable_frame, text="Fault Configuration", padding=10)
        fault_config_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fault_enabled_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            fault_config_frame, text="Enable Fault", variable=self.fault_enabled_var,
            command=self.update_fault_config
        ).pack(anchor=tk.W, pady=5)

        type_frame = ttk.Frame(fault_config_frame)
        type_frame.pack(fill=tk.X, pady=5)
        ttk.Label(type_frame, text="Fault Type:").pack(side=tk.LEFT)

        self.fault_type_var = tk.StringVar(value="friction")
        fault_types = [("Friction", "friction"), ("Power Limit", "power_limit"),
                       ("Encoder", "encoder"), ("Battery", "battery")]
        tgrp = ttk.Frame(type_frame); tgrp.pack(side=tk.LEFT, padx=5)
        for txt, val in fault_types:
            ttk.Radiobutton(tgrp, text=txt, value=val, variable=self.fault_type_var,
                            command=self.update_fault_config).pack(side=tk.LEFT, padx=5)

        self.fault_description_frame = ttk.LabelFrame(fault_config_frame, text="Fault Description", padding=10)
        self.fault_description_frame.pack(fill=tk.X, pady=5)
        self.fault_description_label = ttk.Label(
            self.fault_description_frame,
            text="Select a fault type to see its description",
            wraplength=500, justify=tk.LEFT
        )
        self.fault_description_label.pack(fill=tk.X, pady=5)

        self.params_frame = ttk.LabelFrame(fault_config_frame, text="Fault Parameters", padding=10)
        self.params_frame.pack(fill=tk.X, pady=5)

        self.fault_mag = tk.DoubleVar(value=0.0005)
        self.fault_time = tk.DoubleVar(value=10.0)
        self.fault_wheel = tk.IntVar(value=4)
        self.create_parameter_widgets()

        # ----- Periodic -----
        periodic_frame = ttk.LabelFrame(fault_config_frame, text="Periodic Fault", padding=10)
        periodic_frame.pack(fill=tk.X, pady=5)

        self.periodic_enabled_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(periodic_frame, text="Enable Periodic Fault",
                        variable=self.periodic_enabled_var,
                        command=self.toggle_periodic).pack(anchor=tk.W, pady=5)

        pf = ttk.Frame(periodic_frame); pf.pack(fill=tk.X)
        ttk.Label(pf, text="Interval (sec):").grid(row=0, column=0, sticky=tk.W)
        self.periodic_interval_var = tk.DoubleVar(value=360.0)
        self.periodic_interval_entry = ttk.Entry(pf, textvariable=self.periodic_interval_var, width=10, state="disabled")
        self.periodic_interval_entry.grid(row=0, column=1, padx=5)

        ttk.Label(pf, text="Magnitude:").grid(row=1, column=0, sticky=tk.W)
        self.periodic_magnitude_var = tk.DoubleVar(value=0.1)
        self.periodic_magnitude_entry = ttk.Entry(pf, textvariable=self.periodic_magnitude_var, width=10, state="disabled")
        self.periodic_magnitude_entry.grid(row=1, column=1, padx=5)

        ttk.Label(pf, text="Wheel Number:").grid(row=2, column=0, sticky=tk.W)
        self.periodic_wheel_var = tk.IntVar(value=1)
        self.periodic_wheel_spinbox = ttk.Spinbox(pf, from_=1, to=4, textvariable=self.periodic_wheel_var,
                                                  width=5, state="disabled")
        self.periodic_wheel_spinbox.grid(row=2, column=1, padx=5)

        # ----- Buttons -----
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.pack(fill=tk.X, padx=10, pady=20)
        ttk.Button(button_frame, text="Apply Fault Configuration",
                   command=self.apply_fault_config).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Apply to All Satellites",
                   command=self.apply_to_all_satellites).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Apply to Cluster",
                   command=self.apply_to_cluster).pack(side=tk.LEFT, padx=5)

        ttk.Frame(scrollable_frame, height=20).pack()  # bottom padding

    def _create_fault_summary_tab(self, parent):
        main = ttk.Frame(parent); main.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        hdr = ttk.Frame(main); hdr.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(hdr, text="Fault Summary Overview", font=('Segoe UI', 12, 'bold')).pack(side=tk.LEFT)

        ctrls = ttk.Frame(main); ctrls.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(ctrls, text="Refresh Summary", command=self.update_fault_summary).pack(side=tk.LEFT, padx=5)
        ttk.Button(ctrls, text="Clear All Faults", command=self.clear_all_faults).pack(side=tk.LEFT, padx=5)

        stats = ttk.LabelFrame(ctrls, text="Statistics", padding=5)
        stats.pack(side=tk.RIGHT, padx=10)
        self.fault_stats_label = ttk.Label(stats, text="Total: 0 satellites | Faults: 0 enabled")
        self.fault_stats_label.pack()

        tree_frame = ttk.LabelFrame(main, text="Satellite Fault Status", padding=10)
        tree_frame.pack(fill=tk.BOTH, expand=True)

        cols = ('Satellite', 'Cluster', 'Role', 'Fault Status', 'Type', 'Wheel', 'Time (min)', 'Magnitude')
        self.fault_summary_tree = ttk.Treeview(tree_frame, columns=cols, show='headings', height=15)
        for c in cols:
            self.fault_summary_tree.heading(c, text=c)
            self.fault_summary_tree.column(c, width=110 if c != 'Satellite' else 140)
        sb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.fault_summary_tree.yview)
        self.fault_summary_tree.configure(yscrollcommand=sb.set)
        self.fault_summary_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

    def _create_cluster_fault_tab(self, parent):
        main = ttk.Frame(parent); main.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        head = ttk.Frame(main); head.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(head, text="Cluster Fault Analysis", font=('Segoe UI', 12, 'bold')).pack(side=tk.LEFT)
        ctrls = ttk.Frame(main); ctrls.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(ctrls, text="Update Visualization", command=self._update_cluster_visualization).pack(side=tk.LEFT, padx=5)

        self.cluster_fig = Figure(figsize=(12, 6))
        self.cluster_canvas = FigureCanvasTkAgg(self.cluster_fig, main)
        self.cluster_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._update_cluster_visualization()

    # -----------------------------
    # Cluster / satellite selection
    # -----------------------------
    def on_cluster_selected(self, event=None):
        cname = self.cluster_var.get()
        if not cname or not hasattr(self.parent_app, 'constellation_tab'):
            return
        cluster = next((c for c in self.parent_app.constellation_tab.clusters if c['name'] == cname), None)
        if not cluster:
            return
        sats = cluster.get('satellites', [])
        self.fault_satellite_combo['values'] = sats
        if sats:
            self.fault_satellite_combo.current(0)
            self.on_fault_satellite_changed()

    def show_all_satellites(self):
        self.update_satellite_dropdown()
        self.cluster_var.set("All Satellites")
        if self.satellites:
            self.fault_satellite_combo.current(0)
            self.on_fault_satellite_changed()

    # -----------------------------
    # Apply config
    # -----------------------------
    def apply_to_cluster(self):
        cname = self.cluster_var.get()
        if not cname or cname == "All Satellites":
            messagebox.showwarning("No Cluster", "Please select a specific cluster")
            return
        if not hasattr(self.parent_app, 'constellation_tab'):
            return
        cluster = next((c for c in self.parent_app.constellation_tab.clusters if c['name'] == cname), None)
        if not cluster:
            return

        cfg = self._current_fault_config_payload()
        count = 0
        for sname in cluster.get('satellites', []):
            sat = next((s for s in self.satellites if s['name'] == sname), None)
            if sat:
                sat["fault"] = cfg.copy()
                count += 1

        self.parent_app.add_log(f"Applied fault to {count} satellites in cluster {cname}")
        self.update_fault_summary()
        self._update_cluster_visualization()
        messagebox.showinfo("Cluster Fault", f"Fault applied to {count} satellites in {cname}")

    def _update_cluster_visualization(self):
        if not hasattr(self, 'cluster_fig'):
            return
        self.cluster_fig.clear()

        clusters = getattr(getattr(self.parent_app, 'constellation_tab', None), 'clusters', [])
        ax = self.cluster_fig.add_subplot(121)
        ax2 = self.cluster_fig.add_subplot(122)

        if not clusters:
            ax.text(0.5, 0.5, 'No clusters configured', ha='center', va='center', fontsize=14)
            self.cluster_canvas.draw()
            return

        # Plot 1: totals vs faulted
        names, totals, faulted = [], [], []
        for c in clusters:
            names.append(c['name'])
            csats = c.get('satellites', [])
            totals.append(len(csats))
            faulted.append(sum(1 for n in csats
                               for s in self.satellites if s['name'] == n and s.get('fault', {}).get('enabled')))
        x = np.arange(len(names)); w = 0.35
        b1 = ax.bar(x - w/2, totals, w, label='Total', alpha=0.7)
        b2 = ax.bar(x + w/2, faulted, w, label='Faulted', alpha=0.7)
        ax.set_xlabel('Cluster'); ax.set_ylabel('Number of Satellites')
        ax.set_title('Fault Distribution Across Clusters')
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha='right')
        ax.legend(); ax.grid(True, alpha=0.3)
        for b in (*b1, *b2):
            h = b.get_height()
            if h > 0:
                ax.text(b.get_x() + b.get_width()/2., h, f'{int(h)}', ha='center', va='bottom')

        # Plot 2: fault type share
        types = {}
        for s in self.satellites:
            if s.get('fault', {}).get('enabled'):
                types[s['fault']['type']] = types.get(s['fault']['type'], 0) + 1
        if types:
            labels, sizes = list(types.keys()), list(types.values())
            ax2.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
            ax2.set_title('Fault Type Distribution')
        else:
            ax2.text(0.5, 0.5, 'No faults configured', ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Fault Type Distribution')

        self.cluster_fig.tight_layout()
        self.cluster_canvas.draw()

    def clear_all_faults(self):
        if not self.satellites:
            return
        if not messagebox.askyesno("Clear All Faults", "Disable all faults on all satellites?"):
            return
        for sat in self.satellites:
            sat["fault"]["enabled"] = False
            if "periodic" in sat["fault"]:
                sat["fault"]["periodic"]["enabled"] = False
        idx = self.fault_satellite_combo.current()
        if idx is not None and idx >= 0:
            self.load_fault_config(idx)
        self.update_fault_summary()
        self._update_cluster_visualization()
        self.parent_app.add_log("Cleared all faults from all satellites")
        messagebox.showinfo("Success", "All faults have been disabled")

    def update_satellite_dropdown(self):
        self.fault_satellite_combo['values'] = [s['name'] for s in self.satellites]
        if self.satellites:
            self.fault_satellite_combo.current(0)

        # clusters list
        clusters = []
        if hasattr(self.parent_app, 'constellation_tab'):
            clusters = [c['name'] for c in getattr(self.parent_app.constellation_tab, 'clusters', [])]
        self.cluster_combo['values'] = ['All Satellites'] + clusters
        self.cluster_combo.set('All Satellites')

        if hasattr(self, 'cluster_fig'):
            self._update_cluster_visualization()

    def on_fault_satellite_changed(self, event=None):
        name = self.fault_satellite_var.get()
        for i, s in enumerate(self.satellites):
            if s["name"] == name:
                self.load_fault_config(i)
                cluster = s.get('cluster', 'Individual')
                role = s.get('role', 'independent')
                orbit = s.get('orbit_name', 'Unknown')
                self.sat_info_label.config(text=f"Cluster: {cluster} | Role: {role} | Orbit: {orbit}")
                break

    # -----------------------------
    # Fault config / params
    # -----------------------------
    def load_fault_config(self, index):
        if not (0 <= index < len(self.satellites)):
            return
        self._loading_config = True
        f = self.satellites[index]["fault"]

        self.fault_enabled_var.set(f.get("enabled", False))
        self.fault_type_var.set(f.get("type", "friction"))
        self.fault_mag.set(f.get("magnitude", 0.0005))
        self.fault_wheel.set(f.get("wheel", 0) + 1)  # show 1..4
        self.fault_time.set(f.get("time", 10.0))

        p = f.get("periodic", {})
        self.periodic_enabled_var.set(p.get("enabled", False))
        self.periodic_interval_var.set(p.get("interval", 360.0))
        self.periodic_magnitude_var.set(p.get("magnitude", 0.1))
        self.periodic_wheel_var.set(p.get("wheel", 0) + 1)

        self.update_fault_config()
        self.update_fault_status()
        self._loading_config = False

    def on_fault_type_change(self, *args):
        ft = self.fault_type_var.get()
        defaults = {"friction": 0.0005, "power_limit": 0.5, "encoder": 20.0, "battery": 50.0}
        if ft in defaults:
            self.fault_mag.set(defaults[ft])
            if not getattr(self, "_loading_config", False):
                self.parent_app.add_log(f"Updated {ft} default magnitude → {defaults[ft]}")

    def update_fault_config(self):
        for w in self.params_frame.winfo_children():
            w.destroy()
        ft = self.fault_type_var.get()
        self.update_fault_description()

        if ft == "friction":
            self.create_friction_params()
        elif ft == "power_limit":
            self.create_power_limit_params()
        elif ft == "encoder":
            self.create_encoder_params()
        elif ft == "battery":
            self.create_battery_params()

        self.toggle_periodic()
        self.update_fault_status()
        self.update_fault_summary()

    def update_fault_description(self):
        ft = self.fault_type_var.get()
        desc = {
            "friction": "Friction fault increases the Coulomb friction in the reaction wheel. Default: 0.0005 N⋅m.",
            "power_limit": "Power limit fault restricts maximum power. Default: 0.5 W.",
            "encoder": "Encoder fault causes measurement errors. Default: 20% error.",
            "battery": "Battery fault simulates increased power drain. Default: 50 W additional drain."
        }.get(ft, "No description available")
        self.fault_description_label.config(text=desc)

    def update_fault_status(self):
        if self.fault_enabled_var.get():
            self.fault_status_label.config(
                text=f"Fault: {self.fault_type_var.get().replace('_',' ').title()} on RW{self.fault_wheel.get()}"
            )
        else:
            self.fault_status_label.config(text="No fault enabled")

    def create_parameter_widgets(self):
        self.create_friction_params()

    def create_friction_params(self):
        f = ttk.Frame(self.params_frame); f.pack(fill=tk.X, pady=2)
        ttk.Label(f, text="Magnitude (N⋅m):").pack(side=tk.LEFT)
        ttk.Entry(f, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)

        f2 = ttk.Frame(self.params_frame); f2.pack(fill=tk.X, pady=2)
        ttk.Label(f2, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(f2, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)

        f3 = ttk.Frame(self.params_frame); f3.pack(fill=tk.X, pady=2)
        ttk.Label(f3, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(f3, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)

    def create_power_limit_params(self):
        f = ttk.Frame(self.params_frame); f.pack(fill=tk.X, pady=2)
        ttk.Label(f, text="Power Limit (W):").pack(side=tk.LEFT)
        ttk.Entry(f, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)

        f2 = ttk.Frame(self.params_frame); f2.pack(fill=tk.X, pady=2)
        ttk.Label(f2, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(f2, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)

        f3 = ttk.Frame(self.params_frame); f3.pack(fill=tk.X, pady=2)
        ttk.Label(f3, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(f3, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)

    def create_encoder_params(self):
        f = ttk.Frame(self.params_frame); f.pack(fill=tk.X, pady=2)
        ttk.Label(f, text="Error (％):").pack(side=tk.LEFT)
        ttk.Entry(f, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)

        f2 = ttk.Frame(self.params_frame); f2.pack(fill=tk.X, pady=2)
        ttk.Label(f2, text="Wheel Number:").pack(side=tk.LEFT)
        ttk.Spinbox(f2, from_=1, to=4, textvariable=self.fault_wheel, width=5).pack(side=tk.LEFT, padx=5)

        f3 = ttk.Frame(self.params_frame); f3.pack(fill=tk.X, pady=2)
        ttk.Label(f3, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(f3, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)

    def create_battery_params(self):
        f = ttk.Frame(self.params_frame); f.pack(fill=tk.X, pady=2)
        ttk.Label(f, text="Power Drain (W):").pack(side=tk.LEFT)
        ttk.Entry(f, textvariable=self.fault_mag, width=10).pack(side=tk.LEFT, padx=5)

        f3 = ttk.Frame(self.params_frame); f3.pack(fill=tk.X, pady=2)
        ttk.Label(f3, text="Fault Time (min):").pack(side=tk.LEFT)
        ttk.Entry(f3, textvariable=self.fault_time, width=10).pack(side=tk.LEFT, padx=5)

    def toggle_periodic(self):
        state = "normal" if self.periodic_enabled_var.get() else "disabled"
        self.periodic_interval_entry.config(state=state)
        self.periodic_magnitude_entry.config(state=state)
        self.periodic_wheel_spinbox.config(state=state)

    def _current_fault_config_payload(self):
        return {
            "enabled": self.fault_enabled_var.get(),
            "type": self.fault_type_var.get(),
            "magnitude": self.fault_mag.get(),
            "wheel": self.fault_wheel.get() - 1,
            "time": self.fault_time.get(),
            "periodic": {
                "enabled": self.periodic_enabled_var.get(),
                "interval": self.periodic_interval_var.get(),
                "magnitude": self.periodic_magnitude_var.get(),
                "wheel": self.periodic_wheel_var.get() - 1
            }
        }

    def apply_to_all_satellites(self):
        if not self.satellites:
            return
        if not messagebox.askyesno("Confirm", "Apply the current fault configuration to ALL satellites?"):
            return
        cfg = self._current_fault_config_payload()
        for s in self.satellites:
            s["fault"] = cfg.copy()
        self.parent_app.add_log("Applied fault configuration to all satellites")
        messagebox.showinfo("Success", "Fault configuration applied to all satellites")
        self.update_fault_summary()

    def apply_fault_config(self):
        sname = self.fault_satellite_var.get()
        idx = next((i for i, s in enumerate(self.satellites) if s["name"] == sname), -1)
        if idx < 0:
            return
        cfg = self._current_fault_config_payload()
        self.satellites[idx]["fault"] = cfg
        self.update_fault_status()
        self.parent_app.add_log(f"{'Applied' if cfg['enabled'] else 'Disabled'} {cfg['type']} fault for {sname}")
        self.update_fault_summary()

    def update_fault_summary(self):
        if not hasattr(self, 'fault_summary_tree'):
            return
        tree = self.fault_summary_tree
        for i in tree.get_children():
            tree.delete(i)

        total = len(self.satellites)
        enabled = 0
        for sat in self.satellites:
            f = sat["fault"]
            if f.get("enabled"):
                enabled += 1
            status = "ENABLED" if f.get("enabled") else "Disabled"
            wheel_display = f"RW{(f.get('wheel', 0))+1}"
            magnitude_str = str(f.get("magnitude"))
            t = f.get("type", "")
            if t == "friction":
                magnitude_str = f"{f['magnitude']} N⋅m"
            elif t in ("power_limit", "battery"):
                magnitude_str = f"{f['magnitude']} W"
            elif t == "encoder":
                magnitude_str = f"{f['magnitude']}%"

            item = tree.insert('', 'end', values=(
                sat["name"],
                sat.get('cluster', '—'),
                sat.get('role', '—'),
                status,
                t.replace('_', ' ').title(),
                wheel_display,
                f"{f.get('time', 0)} min",
                magnitude_str
            ))
            tree.item(item, tags=('enabled',) if f.get("enabled") else ('disabled',))

        tree.tag_configure('enabled', foreground='darkred', font=('Segoe UI', 10, 'bold'))
        tree.tag_configure('disabled', foreground='gray')

        if hasattr(self, 'fault_stats_label'):
            self.fault_stats_label.config(text=f"Total: {total} satellites | Faults: {enabled} enabled")
        self.parent_app.add_log(f"Fault summary updated: {enabled}/{total} faults enabled")
