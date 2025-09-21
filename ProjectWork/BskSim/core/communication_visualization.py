#!/usr/bin/env python
# core/communication_visualization.py
"""
Integrated Communication Visualization

 real-time Timeline + Network Graph + live simulation + text log behavior

Works with either:
- a `cluster_manager` that exposes .clusters and .satellites (list OR dict shape)
- or a pre-baked `clusters_data` dict (via update_clusters)

Notes
- Range/FOV controls are on top (affect Windows + Matrix).
- Start/Stop/Clear controls drive the real-time Timeline + Network Graph and also
  add messages into the structured Message Log.
"""

import os
import time
import tkinter as tk
from tkinter import ttk, messagebox
from datetime import datetime
import numpy as np

import matplotlib
# Using TkAgg by default for embedded canvases
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D


class CommunicationVisualizer:
    """Unified communication visualization (Windows, Matrix, Timeline, Network, Log)."""

    # Palette used to tag clusters consistently across views
    PALETTE = ["#2ecc71", "#3498db", "#e67e22", "#9b59b6", "#e74c3c", "#1abc9c", "#f1c40f"]

    def __init__(self, parent_window, cluster_manager=None, clusters_data=None):
        self.parent_window = parent_window
        self.cluster_manager = cluster_manager

        # Window
        if self.parent_window is None:
            self.window = tk.Tk()
            self.window.title("Communication Visualization")
        else:
            self.window = tk.Toplevel(self.parent_window)
            self.window.title("Communication Visualization")
        self.window.geometry("1200x850")
        self.window.protocol("WM_DELETE_WINDOW", self._on_close)

        # Settings / state
        self.range_km = 2000.0
        self.fov_deg = 30.0
        self.simulation_time_min = 30.0

        # Data holders
        self._clusters_data = clusters_data  # canonical dict form (can be None; we derive when needed)
        self.comm_history = []  # for Timeline arrows
        self.messages = []      # structured log table
        self._matrix_cbar = None

        # Animation control
        self.animation_running = False
        self.animation = None

        # UI
        self._build_ui()

        # Initial draw
        self.refresh_all_views()
        self.start_animation()

    # --------------------------------------------------------------------- UI
    def _build_ui(self):
        main = ttk.Frame(self.window, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        # === Control panel ===
        ctrl = ttk.LabelFrame(main, text="Communication Controls", padding=10)
        ctrl.pack(fill=tk.X, pady=(0, 10))

        # Left controls: range/fov + refresh/help
        ttk.Label(ctrl, text="Range (km):").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.range_var = tk.DoubleVar(value=self.range_km)
        ttk.Entry(ctrl, textvariable=self.range_var, width=10).grid(row=0, column=1, padx=2)

        ttk.Label(ctrl, text="FOV (deg):").grid(row=0, column=2, sticky=tk.W, padx=8)
        self.fov_var = tk.DoubleVar(value=self.fov_deg)
        ttk.Entry(ctrl, textvariable=self.fov_var, width=10).grid(row=0, column=3, padx=2)

        ttk.Button(ctrl, text="Refresh Views", command=self.refresh_all_views).grid(row=0, column=4, padx=12)
        ttk.Button(ctrl, text="Help", command=self.show_help).grid(row=0, column=5, padx=4)

        ttk.Label(ctrl, text="(Bars = comm windows, dots = messages)", foreground="gray").grid(
            row=0, column=6, padx=10
        )

        # Right side: sim controls
        ttk.Button(ctrl, text="Start Simulation", command=self.start_simulation).grid(row=0, column=7, padx=10)
        ttk.Button(ctrl, text="Stop Simulation", command=self.stop_simulation).grid(row=0, column=8, padx=4)
        ttk.Button(ctrl, text="Clear History", command=self.clear_history).grid(row=0, column=9, padx=4)

        self.time_label = ttk.Label(ctrl, text="Time: 0.0 min", font=("Arial", 11, "bold"))
        self.time_label.grid(row=0, column=10, padx=10, sticky=tk.E)

        # === Notebook (tabs) ===
        self.nb = ttk.Notebook(main)
        self.nb.pack(fill=tk.BOTH, expand=True)

        self.windows_frame = ttk.Frame(self.nb)   # Gantt
        self.matrix_frame = ttk.Frame(self.nb)    # Heatmap matrix
        self.timeline_frame = ttk.Frame(self.nb)  # Arrows over time
        self.network_frame = ttk.Frame(self.nb)   # Topology
        self.log_frame = ttk.Frame(self.nb)       # Structured log

        self.nb.add(self.windows_frame, text="Windows (Prototype)")
        self.nb.add(self.matrix_frame, text="Links (Matrix)")
        self.nb.add(self.timeline_frame, text="Timeline View")
        self.nb.add(self.network_frame, text="Network Graph")
        self.nb.add(self.log_frame, text="Message Log")

        # Build tab contents
        self._build_windows_tab(self.windows_frame)
        self._build_matrix_tab(self.matrix_frame)
        self._build_timeline_tab(self.timeline_frame)
        self._build_network_tab(self.network_frame)
        self._build_log_tab(self.log_frame)

    # ------------------------------- Windows (Gantt)
    def _build_windows_tab(self, parent):
        self.windows_fig = Figure(figsize=(12, 8), facecolor="white")
        self.windows_ax = self.windows_fig.add_subplot(111)
        self.windows_canvas = FigureCanvasTkAgg(self.windows_fig, parent)
        self.windows_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        bar = ttk.Frame(parent); bar.pack(fill=tk.X, pady=6)
        ttk.Button(bar, text="Save Windows Plot", command=lambda: self._save_plot("windows")).pack(side=tk.LEFT, padx=5)

    # ------------------------------- Matrix
    def _build_matrix_tab(self, parent):
        top = ttk.Frame(parent); top.pack(fill=tk.X, pady=4)
        self.matrix_scope_var = tk.StringVar(value="intra")
        ttk.Radiobutton(top, text="Intra-cluster only", variable=self.matrix_scope_var, value="intra").pack(
            side=tk.LEFT, padx=5
        )
        ttk.Radiobutton(top, text="All satellites", variable=self.matrix_scope_var, value="all").pack(
            side=tk.LEFT, padx=5
        )
        ttk.Button(top, text="Update Matrix", command=self.update_links_matrix).pack(side=tk.LEFT, padx=12)
        ttk.Button(top, text="Save Matrix Plot", command=lambda: self._save_plot("matrix")).pack(side=tk.LEFT, padx=5)

        self.links_fig = Figure(figsize=(10, 8), facecolor="white")
        self.links_ax = self.links_fig.add_subplot(111)
        self.links_canvas = FigureCanvasTkAgg(self.links_fig, parent)
        self.links_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------------------------- Timeline
    def _build_timeline_tab(self, parent):
        self.timeline_fig = Figure(figsize=(10, 6), dpi=100)
        self.timeline_ax = self.timeline_fig.add_subplot(111)
        self.timeline_canvas = FigureCanvasTkAgg(self.timeline_fig, parent)
        self.timeline_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------------------------- Network
    def _build_network_tab(self, parent):
        self.network_fig = Figure(figsize=(8, 8), dpi=100)
        self.network_ax = self.network_fig.add_subplot(111)
        self.network_canvas = FigureCanvasTkAgg(self.network_fig, parent)
        self.network_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------------------------- Log (structured)
    def _build_log_tab(self, parent):
        # Treeview with columns (Claude style)
        cols = ("Time (min)", "From", "To", "Type", "Size (KB)", "Status")
        self.log_tree = ttk.Treeview(parent, columns=cols, show="headings", height=20)
        for c in cols:
            self.log_tree.heading(c, text=c)
            self.log_tree.column(c, width=120 if c != "Type" else 100)
        self.log_tree.pack(fill=tk.BOTH, expand=True, padx=6, pady=(6, 0))

        sb = ttk.Scrollbar(parent, orient="vertical", command=self.log_tree.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_tree.configure(yscrollcommand=sb.set)

        btns = ttk.Frame(parent); btns.pack(fill=tk.X, pady=6)
        ttk.Button(btns, text="Clear Log", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(btns, text="Export Log", command=self.export_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(btns, text="Generate Test Messages", command=self.generate_test_messages).pack(side=tk.LEFT, padx=5)

        # Tag colors
        self.log_tree.tag_configure("failed", foreground="red")
        self.log_tree.tag_configure("pending", foreground="orange")

    # =============================================================== Views API
    def refresh_all_views(self):
        self.range_km = float(self.range_var.get())
        self.fov_deg = float(self.fov_var.get())

        self.update_windows_gantt()
        self.update_links_matrix()
        if not self.messages:
            self.generate_test_messages()

        # Also refresh realtime views once
        self.update_timeline_plot()
        self.update_network_plot()

    # --------------------------------------------------------- Windows (Gantt)
    def update_windows_gantt(self):
        self.windows_ax.clear()
        clusters = self._get_clusters_data()
        if not clusters:
            self.windows_ax.text(0.5, 0.5, "No clusters available", ha="center", va="center")
            self.windows_canvas.draw(); return

        self.windows_ax.set_title(
            f"Communication Windows for {len(clusters)} Clusters", fontsize=14, fontweight="bold"
        )

        y = 0
        y_positions = []
        y_labels = []

        for idx, (cname, cdata) in enumerate(clusters.items()):
            color = self.PALETTE[idx % len(self.PALETTE)]
            leader = cdata.get("leader", {})
            children = cdata.get("children", [])
            for child in children:
                # Bars = windows
                for (t0, t1) in self._calculate_comm_windows(leader, child):
                    r = Rectangle((t0, y - 0.3), t1 - t0, 0.6, facecolor=color, alpha=0.30, edgecolor=color, linewidth=1)
                    self.windows_ax.add_patch(r)
                # Dots = actual messages for this link
                for mt in self._get_message_times_for_link(leader.get("name",""), child.get("name","")):
                    self.windows_ax.plot(mt, y, "o", color=color, markersize=6)

                # Left pill with cluster tag
                tag = Rectangle((-2.0, y - 0.35), 1.5, 0.7, facecolor=color, alpha=0.7)
                self.windows_ax.add_patch(tag)
                self.windows_ax.text(-1.25, y, cname[:3].upper(), color="white", fontsize=8,
                                     ha="center", va="center", fontweight="bold")

                y_positions.append(y)
                y_labels.append(f"{leader.get('name','')[:3]} → {child.get('name','')[-4:]}")
                y += 1

        self.windows_ax.set_xlim(-3, self.simulation_time_min)
        self.windows_ax.set_ylim(-0.5, max(0.5, y - 0.5))
        self.windows_ax.set_xlabel("Time (minutes)", fontsize=11)
        self.windows_ax.set_ylabel("Communication Links", fontsize=11)
        self.windows_ax.set_yticks(y_positions)
        self.windows_ax.set_yticklabels(y_labels)
        self.windows_ax.grid(True, axis="x", alpha=0.3)
        self.windows_ax.axvline(0.0, color="k", linewidth=0.5)

        # Legend
        legends = [Patch(facecolor=self.PALETTE[i % len(self.PALETTE)], label=f"{name} cluster", alpha=0.30)
                   for i, name in enumerate(clusters.keys())]
        legends.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                              markersize=8, label="Message sent"))
        self.windows_ax.legend(handles=legends, loc="upper right", fontsize=9)

        self.windows_ax.text(
            0.5, -0.06,
            "Colored bars = communication possible | Dots = messages sent",
            transform=self.windows_ax.transAxes, ha="center", fontsize=9, style="italic", color="gray"
        )
        self.windows_fig.tight_layout()
        self.windows_canvas.draw()

    # ----------------------------------------------------------- Links Matrix
    def update_links_matrix(self):
        self.links_ax.clear()
        clusters = self._get_clusters_data()
        if not clusters:
            self.links_ax.text(0.5, 0.5, "No clusters available", ha="center", va="center")
            self.links_canvas.draw(); return

        # Flatten satellites, keep cluster names and positions
        all_sats = []
        for cname, c in clusters.items():
            for s in c.get("satellites", []):
                s2 = dict(s)  # copy
                s2["cluster"] = cname
                all_sats.append(s2)

        n = len(all_sats)
        if n == 0:
            self.links_ax.text(0.5, 0.5, "No satellites present", ha="center", va="center")
            self.links_canvas.draw(); return

        adj = np.zeros((n, n), dtype=float)
        for i, s1 in enumerate(all_sats):
            for j, s2 in enumerate(all_sats):
                if i == j:
                    continue
                if self.matrix_scope_var.get() == "intra" and s1["cluster"] != s2["cluster"]:
                    continue
                d = self._distance_km(np.array(s1.get("position", [0, 0, 0])),
                                      np.array(s2.get("position", [0, 0, 0])))
                if d < self.range_km:
                    adj[i, j] = max(0.0, 1.0 - d / self.range_km)

        im = self.links_ax.imshow(adj, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")

        # Replace old colorbar
        if self._matrix_cbar:
            try:
                self._matrix_cbar.remove()
            except Exception:
                pass
            self._matrix_cbar = None
        self._matrix_cbar = self.links_fig.colorbar(im, ax=self.links_ax)
        self._matrix_cbar.set_label("Link Quality", rotation=270, labelpad=12)

        labels = [s["name"] for s in all_sats]
        self.links_ax.set_xticks(range(n)); self.links_ax.set_yticks(range(n))
        self.links_ax.set_xticklabels(labels, rotation=90, fontsize=8)
        self.links_ax.set_yticklabels(labels, fontsize=8)

        # Cluster separators
        sep_positions = []
        cur = None
        for i, s in enumerate(all_sats):
            if cur != s["cluster"]:
                if i > 0:
                    sep_positions.append(i - 0.5)
                cur = s["cluster"]
        for p in sep_positions:
            self.links_ax.axhline(p, color="black", linewidth=1.5)
            self.links_ax.axvline(p, color="black", linewidth=1.5)

        scope_text = "All Satellites" if self.matrix_scope_var.get() == "all" else "Intra-cluster Only"
        self.links_ax.set_title(f"Communication Links Matrix — {scope_text}", fontsize=12, fontweight="bold")
        self.links_ax.set_xlabel("Satellite (To)"); self.links_ax.set_ylabel("Satellite (From)")

        self.links_fig.tight_layout()
        self.links_canvas.draw()

    # ----------------------------------------------------------- Timeline view
    def update_timeline_plot(self):
        ax = self.timeline_ax
        ax.clear()

        if not self.comm_history:
            ax.text(0.5, 0.5, "No communication events yet", ha="center", va="center", fontsize=13)
            ax.set_xlim(0, self.simulation_time_min); ax.set_ylim(0, 1)
            self.timeline_canvas.draw(); return

        times = [e["time"] for e in self.comm_history]
        ids = sorted(set([e["sender"] for e in self.comm_history] + [e["receiver"] for e in self.comm_history]))
        ymap = {name: i for i, name in enumerate(ids)}

        for e in self.comm_history:
            y1, y2 = ymap[e["sender"]], ymap[e["receiver"]]
            t = e["time"]
            color = e.get("color", "tab:blue")
            ax.annotate("", xy=(t, y2), xytext=(t, y1),
                        arrowprops=dict(arrowstyle="->", color=color, lw=2, alpha=0.8))
            ax.scatter(t, y2, s=45, color=color, zorder=5, alpha=0.9)

        ax.set_yticks(range(len(ids))); ax.set_yticklabels(ids, fontsize=8)
        ax.set_xlabel("Time (minutes)"); ax.set_ylabel("Satellites")
        ax.set_title("Communication Timeline"); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max(max(times) + 3, self.simulation_time_min)); ax.set_ylim(-0.5, len(ids) - 0.5)

        self.timeline_fig.tight_layout()
        self.timeline_canvas.draw()

    # ----------------------------------------------------------- Network graph
    def update_network_plot(self):
        ax = self.network_ax
        ax.clear()

        clusters = self._get_clusters_data()
        if not clusters:
            ax.text(0.5, 0.5, "No network data available", ha="center", va="center", fontsize=13)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            self.network_canvas.draw(); return

        nC = len(clusters)
        if nC == 0:
            ax.text(0.5, 0.5, "No clusters configured", ha="center", va="center", fontsize=13)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            self.network_canvas.draw(); return

        # Layout clusters on a circle; place leader at center of its mini-cluster, kids around
        node_pos = {}
        node_col = {}
        edges = []

        R_cluster = 5.0
        for i, (cname, c) in enumerate(clusters.items()):
            theta = 2 * np.pi * i / max(1, nC)
            cx, cy = R_cluster * np.cos(theta), R_cluster * np.sin(theta)
            color = self.PALETTE[i % len(self.PALETTE)]

            leader = c.get("leader", {}).get("name")
            if leader:
                node_pos[leader] = (cx, cy)
                node_col[leader] = color

                kids = c.get("children", [])
                m = len(kids)
                if m > 0:
                    r = 1.2
                    for j, ch in enumerate(kids):
                        nm = ch.get("name", f"{cname}_C{j+1}")
                        phi = 2 * np.pi * j / m
                        node_pos[nm] = (cx + r * np.cos(phi), cy + r * np.sin(phi))
                        node_col[nm] = color
                        edges.append((leader, nm))

        # Inter-cluster leader links (dashed)
        leaders = [c.get("leader", {}).get("name") for c in clusters.values() if c.get("leader")]
        for i in range(len(leaders)):
            for j in range(i + 1, len(leaders)):
                a, b = leaders[i], leaders[j]
                if a in node_pos and b in node_pos:
                    x1, y1 = node_pos[a]; x2, y2 = node_pos[b]
                    ax.plot([x1, x2], [y1, y2], linestyle="--", color="tab:red", alpha=0.25, linewidth=2)

        # Draw edges
        for a, b in edges:
            if a in node_pos and b in node_pos:
                x1, y1 = node_pos[a]; x2, y2 = node_pos[b]
                ax.plot([x1, x2], [y1, y2], "k-", alpha=0.35, linewidth=1)

        # Draw nodes
        leader_set = set([nm for nm in leaders if nm])
        for nm, (x, y) in node_pos.items():
            is_leader = nm in leader_set
            size = 280 if is_leader else 160
            marker = "^" if is_leader else "o"
            ax.scatter(x, y, s=size, c=node_col.get(nm, "#95a5a6"), marker=marker,
                       edgecolor="black", linewidth=1.4, zorder=5)
            ax.text(x, y - 0.35, nm[:14], ha="center", fontsize=8)

        ax.set_title("Satellite Network Topology")
        ax.set_aspect("equal"); ax.axis("off")
        self.network_fig.tight_layout()
        self.network_canvas.draw()

    # -------------------------------------------------------------- Log table
    def update_message_log(self):
        # Clear
        for it in self.log_tree.get_children():
            self.log_tree.delete(it)
        # Fill
        for msg in self.messages:
            vals = (f"{msg['time']:.2f}", msg["from"], msg["to"], msg["type"],
                    f"{msg['size']:.1f}", msg["status"])
            tag = "failed" if msg["status"] == "Failed" else ("pending" if msg["status"] == "Pending" else "")
            self.log_tree.insert("", "end", values=vals, tags=(tag,))

    def clear_log(self):
        self.messages = []
        self.update_message_log()

    def export_log(self):
        if not self.messages:
            messagebox.showwarning("No Data", "No messages to export")
            return
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        os.makedirs("plots", exist_ok=True)
        fname = f"message_log_{ts}.txt"
        with open(fname, "w", encoding="utf-8") as f:
            f.write("Communication Message Log\n")
            f.write(f"Generated: {datetime.now()}\n")
            f.write(f"Total Messages: {len(self.messages)}\n\n")
            f.write("Time(min)\tFrom\tTo\tType\tSize(KB)\tStatus\n")
            f.write("-" * 60 + "\n")
            for m in self.messages:
                f.write(f"{m['time']:.2f}\t{m['from']}\t{m['to']}\t{m['type']}\t{m['size']:.1f}\t{m['status']}\n")
        messagebox.showinfo("Export Complete", f"Log exported to {fname}")

    def generate_test_messages(self):
        clusters = self._get_clusters_data()
        if not clusters:
            return
        kinds = ["Telemetry", "Command", "Data", "Status", "Health"]
        stats = ["Delivered", "Delivered", "Delivered", "Pending", "Failed"]
        out = []
        for c in clusters.values():
            leader = c.get("leader", {}).get("name")
            for child in c.get("children", []):
                cname = child.get("name")
                if not leader or not cname:
                    continue
                k = np.random.randint(3, 8)
                for _ in range(k):
                    t = np.random.uniform(0, self.simulation_time_min)
                    size = np.random.uniform(1, 100)
                    typ = np.random.choice(kinds)
                    st = np.random.choice(stats)
                    if np.random.random() > 0.5:
                        frm, to = leader, cname
                    else:
                        frm, to = cname, leader
                    out.append({"time": float(t), "from": frm, "to": to, "type": typ, "size": float(size), "status": st})
        out.sort(key=lambda x: x["time"])
        self.messages = out
        self.update_message_log()

    # ====================================================== Simulation / anim
    def simulate_communication(self, current_time):
        """Append simulated events (and mirror them into the log table)."""
        clusters = self._get_clusters_data()
        if not clusters:
            return

        # Intra-cluster leader<->child
        if np.random.random() > 0.72:
            cnames = list(clusters.keys())
            c = clusters[np.random.choice(cnames)]
            leader = c.get("leader", {}).get("name")
            kids = c.get("children", [])
            if leader and kids:
                child = kids[np.random.randint(0, len(kids))]["name"]
                # leader -> child
                self.comm_history.append({
                    "time": current_time, "sender": leader, "receiver": child,
                    "message": f"Status request from {leader}", "color": "tab:blue"
                })
                self._append_log_entry(current_time, leader, child, "Command", 2.5, "Delivered")

                if np.random.random() > 0.5:
                    # child -> leader
                    self.comm_history.append({
                        "time": current_time + 0.1, "sender": child, "receiver": leader,
                        "message": f"Status OK from {child}", "color": "tab:green"
                    })
                    self._append_log_entry(current_time + 0.1, child, leader, "Status", 1.2, "Delivered")

        # Inter-cluster leader↔leader
        if np.random.random() > 0.90:
            leaders = [c.get("leader", {}).get("name") for c in clusters.values() if c.get("leader")]
            if len(leaders) >= 2:
                a, b = np.random.choice(leaders, 2, replace=False)
                self.comm_history.append({
                    "time": current_time, "sender": a, "receiver": b,
                    "message": f"Inter-cluster sync {a}->{b}", "color": "tab:red"
                })
                self._append_log_entry(current_time, a, b, "Sync", 5.0, "Delivered")

        # Bound history + keep log compact
        self.comm_history = self.comm_history[-600:]
        self.messages = self.messages[-1200:]

    def _append_log_entry(self, t, frm, to, typ, size_kb, status):
        self.messages.append({"time": float(t), "from": frm, "to": to, "type": typ, "size": float(size_kb), "status": status})
        # Lightweight live update of the last few rows for smoothness
        if len(self.messages) % 5 == 0:
            self.update_message_log()

    def update_animation(self, frame):
        if self.animation_running:
            current_time = frame * 0.1  # 0.1 min/step
            self.time_label.config(text=f"Time: {current_time:.1f} min")

            self.simulate_communication(current_time)
            self.update_timeline_plot()

            if frame % 10 == 0:
                self.update_links_matrix()
                self.update_network_plot()

        return []

    def start_animation(self):
        self.animation_running = True
        self.animation = FuncAnimation(self.timeline_fig, self.update_animation, interval=100, blit=False)
        self.timeline_canvas.draw()

    def start_simulation(self):
        self.animation_running = True

    def stop_simulation(self):
        self.animation_running = False

    def clear_history(self):
        self.comm_history.clear()
        self.messages.clear()
        self.update_timeline_plot()
        self.update_message_log()
        self.update_links_matrix()
        self.update_network_plot()

    # ================================================================= Utils
    def _get_clusters_data(self):
        """
        Returns canonical dict:
        {
          "ClusterA": {
             "leader": {"name": "...", "position": [x,y,z]},
             "children": [{"name": "...", "position":[...]}, ...],
             "satellites": [{"name":"...", "position":[...]} ...]
          }, ...
        }
        Uses (in priority): self._clusters_data, else derives from cluster_manager.
        """
        if self._clusters_data:
            return self._clusters_data

        if not self.cluster_manager:
            return {}

        clusters_raw = getattr(self.cluster_manager, "clusters", None)
        satellites_raw = getattr(self.cluster_manager, "satellites", None)

        if clusters_raw is None or satellites_raw is None:
            return {}

        # Build satellite lookup
        sat_by_name = {}
        try:
            for s in satellites_raw:
                nm = s.get("name") if isinstance(s, dict) else getattr(s, "name", str(s))
                if not nm:
                    continue
                # Position heuristic:
                #   use stored 'position' or 'position_offset_km', else compute from orbit a/f
                pos = None
                if isinstance(s, dict):
                    pos = s.get("position")
                    if pos is None:
                        off = s.get("position_offset_km")
                        if off is not None:
                            pos = [float(off[0]), float(off[1]), float(off[2] if len(off) > 2 else 0.0)]
                if pos is None and isinstance(s, dict):
                    try:
                        a = float(s.get("orbit", {}).get("a", 6771))  # km from Earth's center
                        f = np.deg2rad(float(s.get("orbit", {}).get("f", 0.0)))
                        r = max(100.0, a - 6371.0)  # altitude-based radius
                        dx, dy = 0.0, 0.0
                        off = s.get("position_offset_km")
                        if off is not None:
                            dx += float(off[0]); dy += float(off[1])
                        pos = [r * np.cos(f) + dx, r * np.sin(f) + dy, 0.0]
                    except Exception:
                        pos = [0.0, 0.0, 0.0]
                sat_by_name[nm] = {"name": nm, "position": pos or [0.0, 0.0, 0.0]}
        except Exception:
            pass

        clusters = {}

        # clusters might be a list (from ConstellationTab) or dict (other manager)
        if isinstance(clusters_raw, list):
            for idx, c in enumerate(clusters_raw):
                cname = c.get("name", f"Cluster{idx+1}")
                leader_nm = c.get("leader")
                kids_nms = c.get("children", [])
                sat_nms = c.get("satellites", [])

                leader = {"name": leader_nm, "position": sat_by_name.get(leader_nm, {}).get("position", [0, 0, 0])} if leader_nm else {}
                children = [{"name": nm, "position": sat_by_name.get(nm, {}).get("position", [0, 0, 0])} for nm in kids_nms]
                sats = [{"name": nm, "position": sat_by_name.get(nm, {}).get("position", [0, 0, 0])} for nm in sat_nms]
                clusters[cname] = {"leader": leader, "children": children, "satellites": sats}

        elif isinstance(clusters_raw, dict):
            for cname, c in clusters_raw.items():
                # Try common shapes
                leader = c.get("leader")
                if isinstance(leader, dict):
                    leader_nm = leader.get("name")
                else:
                    leader_nm = getattr(leader, "name", leader) if leader is not None else None
                kids = c.get("children", [])
                satlist = c.get("satellites", [])

                def _nm(x):
                    return x.get("name") if isinstance(x, dict) else getattr(x, "name", x)

                leader_out = {"name": leader_nm, "position": sat_by_name.get(leader_nm, {}).get("position", [0, 0, 0])} if leader_nm else {}
                children_out = [{"name": _nm(k), "position": sat_by_name.get(_nm(k), {}).get("position", [0, 0, 0])} for k in kids]
                sats_out = [{"name": _nm(s), "position": sat_by_name.get(_nm(s), {}).get("position", [0, 0, 0])} for s in satlist]
                clusters[cname] = {"leader": leader_out, "children": children_out, "satellites": sats_out}

        return clusters

    def _calculate_comm_windows(self, sat1, sat2):
        """
        Simplified synthetic windows over the simulation horizon.
        Returns list of (start, end) minutes.
        """
        windows = []
        num = np.random.randint(3, 5)
        span = self.simulation_time_min / max(1, num)
        for i in range(num):
            t0 = i * span
            dur = 3.0 + np.random.random() * 2.0
            t1 = min(t0 + dur, self.simulation_time_min)
            windows.append((t0, t1))
        return windows

    def _get_message_times_for_link(self, frm, to):
        """Extract times from current structured log for frm→to; fallback to a few random dots."""
        ts = [m["time"] for m in self.messages if m["from"] == frm and m["to"] == to]
        if ts:
            return sorted(ts)
        # fallback sparse dots
        return sorted(list(np.random.uniform(0, self.simulation_time_min, size=np.random.randint(2, 4))))

    @staticmethod
    def _distance_km(p1, p2):
        try:
            return float(np.linalg.norm(np.array(p2, dtype=float) - np.array(p1, dtype=float)))
        except Exception:
            return 0.0

    # ------------------------------------------------------------- Save/Help/UI
    def _save_plot(self, which):
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        os.makedirs("plots", exist_ok=True)
        if which == "windows":
            fn = f"plots/{ts}_CommunicationWindows_AllClusters.png"
            self.windows_fig.savefig(fn, dpi=150, bbox_inches="tight")
        elif which == "matrix":
            scope = "All" if self.matrix_scope_var.get() == "all" else "IntraCluster"
            fn = f"plots/{ts}_CommunicationMatrix_{scope}.png"
            self.links_fig.savefig(fn, dpi=150, bbox_inches="tight")
        else:
            return
        messagebox.showinfo("Saved", f"Plot saved to {fn}")

    def show_help(self):
        messagebox.showinfo(
            "Communication Visualization Help",
            """Range (km): Max communication distance
FOV (deg): Antenna field of view

Windows tab: bars show when links are available; dots show messages.
Links tab: heatmap of link quality (green=strong); switch intra/all.
Timeline: arrows at the time messages are sent.
Network: leaders and children, with dashed inter-leader links.
Message Log: sortable table of messages (time, from, to, type, size, status)."""
        )

    def update_clusters(self, clusters_data):
        """External API to push pre-baked clusters dict (overrides cluster_manager)."""
        self._clusters_data = clusters_data
        self.refresh_all_views()

    # --------------------------------------------------------------- Lifecycle
    def _on_close(self):
        self.stop_simulation()
        self.animation = None
        self.window.destroy()


def show_communication_visualization(parent_window, cluster_manager=None, clusters_data=None):
    """Helper to open the visualization window."""
    return CommunicationVisualizer(parent_window, cluster_manager, clusters_data)
