#!/usr/bin/env python
"""
results_tab.py

Results viewing and analysis tab for the Spacecraft Constellation Fault Simulator.
"""

import os
import sys
import subprocess
import shutil
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from datetime import datetime

# Matplotlib (embed in Tk; avoid pyplot to prevent backend clashes)
import matplotlib
if "agg" in matplotlib.get_backend().lower():
    try:
        matplotlib.use("TkAgg")  # best-effort switch for embedded viewing
    except Exception:
        pass

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from PIL import Image
import numpy as np


class ResultsTab:
    def __init__(self, app, parent):
        self.app = app
        self.parent = parent

        self.plot_figure = None
        self.plot_canvas = None
        self.plot_toolbar = None
        self.plot_listbox = None
        self.plot_info_text = None
        self.current_plot_path = None

        self.create_ui()

    # ---------- UI ----------

    def create_ui(self):
        main_frame = ttk.Frame(self.parent, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Top controls
        controls_frame = ttk.Frame(main_frame)
        controls_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(controls_frame, text="Simulation Results", style="Title.TLabel").pack(side=tk.LEFT)

        ttk.Button(controls_frame, text="Refresh", command=self.refresh_plot_list).pack(side=tk.RIGHT, padx=5)
        ttk.Button(controls_frame, text="Open Results Folder", command=self.app.open_results_folder).pack(side=tk.RIGHT, padx=5)
        ttk.Button(controls_frame, text="Export Plot", command=self.export_current_plot).pack(side=tk.RIGHT, padx=5)
        ttk.Button(controls_frame, text="Open Plot", command=self._open_in_os).pack(side=tk.RIGHT, padx=5)

        # Main content area
        split_frame = ttk.Frame(main_frame)
        split_frame.pack(fill=tk.BOTH, expand=True)

        # Left: plot list
        list_frame = ttk.LabelFrame(split_frame, text="Available Plots", padding=10, width=320)
        list_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        list_frame.pack_propagate(False)

        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.plot_listbox = tk.Listbox(
            list_container,
            selectmode=tk.SINGLE,
            yscrollcommand=scrollbar.set,
            font=("Segoe UI", 10),
            exportselection=False
        )
        self.plot_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.plot_listbox.yview)

        self.plot_listbox.bind("<<ListboxSelect>>", self.on_plot_selected)
        self.plot_listbox.bind("<Double-1>", lambda e: self._open_in_os())

        # Right: plot display
        display_frame = ttk.LabelFrame(split_frame, text="Plot View", padding=10)
        display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))

        plot_container = ttk.Frame(display_frame)
        plot_container.pack(fill=tk.BOTH, expand=True)

        self.plot_figure = Figure(figsize=(10, 6), dpi=100)
        self.plot_canvas = FigureCanvasTkAgg(self.plot_figure, plot_container)

        toolbar_frame = ttk.Frame(plot_container)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        self.plot_toolbar = NavigationToolbar2Tk(self.plot_canvas, toolbar_frame)
        self.plot_toolbar.update()

        self.plot_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Bottom info
        info_frame = ttk.LabelFrame(main_frame, text="Plot Information", padding=10)
        info_frame.pack(fill=tk.X, pady=(10, 0))

        self.plot_info_text = scrolledtext.ScrolledText(info_frame, height=4, wrap=tk.WORD)
        self.plot_info_text.pack(fill=tk.BOTH, expand=True)
        self.plot_info_text.config(state=tk.DISABLED)

        # Initial state
        self.refresh_plot_list()

    # ---------- Data / Listing ----------

    def refresh_plot_list(self):
        """Refresh the list of available plots from the app's plots directory."""
        self.plot_listbox.delete(0, tk.END)

        # Reset info
        self._set_info_text("Select a plot to view details and analysis")

        # Reset canvas content
        self._show_placeholder("No plot selected\nSelect a plot from the list to view")

        if not os.path.exists(self.app.plots_dir):
            return

        valid_exts = {".png", ".jpg", ".jpeg", ".pdf", ".svg"}
        plot_files = [
            f for f in os.listdir(self.app.plots_dir)
            if os.path.splitext(f)[1].lower() in valid_exts
        ]
        plot_files.sort(
            key=lambda x: os.path.getmtime(os.path.join(self.app.plots_dir, x)),
            reverse=True
        )

        for filename in plot_files:
            self.plot_listbox.insert(tk.END, filename)

        if plot_files:
            self.plot_listbox.selection_set(0)
            self.on_plot_selected(None)

    # ---------- Selection / Display ----------

    def on_plot_selected(self, event):
        """Update info panel and render the selected plot."""
        selection = self.plot_listbox.curselection()
        if not selection:
            return

        index = selection[0]
        filename = self.plot_listbox.get(index)
        full_path = os.path.join(self.app.plots_dir, filename)
        self.current_plot_path = full_path

        # Update info panel
        self._update_info_for_file(filename, full_path)

        # Render plot/image
        ext = os.path.splitext(filename)[1].lower()
        self.plot_figure.clear()

        try:
            if ext in (".png", ".jpg", ".jpeg"):
                # Detach array to avoid locking files (esp. on Windows)
                with Image.open(full_path) as img:
                    arr = np.array(img)  # auto-handles modes incl. RGBA
                ax = self.plot_figure.add_subplot(111)
                ax.imshow(arr)
                ax.axis("off")

                title = filename.replace("_", " ")
                self.plot_figure.suptitle(title, fontsize=12, y=0.98)

            elif ext in (".pdf", ".svg"):
                # Placeholder for vector formats (no rasterization here)
                ax = self.plot_figure.add_subplot(111)
                ax.text(
                    0.5, 0.5,
                    f"Vector file selected: {filename}\nUse 'Open Plot' to view externally\nor 'Export Plot' to copy.",
                    ha="center", va="center", transform=ax.transAxes, fontsize=12
                )
                ax.axis("off")

            # Layout & draw
            self.plot_figure.tight_layout(rect=[0, 0, 1, 0.95])
            self.plot_canvas.draw()

        except Exception as e:
            self._show_error(f"Error displaying plot:\n{e}")

    # ---------- Export / Open ----------

    def export_current_plot(self):
        """Export the currently selected plot to a chosen location/format."""
        if not self.current_plot_path:
            messagebox.showwarning("No Plot Selected", "Please select a plot first.")
            return

        default_name = os.path.basename(self.current_plot_path)
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            initialfile=default_name,
            filetypes=[
                ("PNG files", "*.png"),
                ("JPEG files", "*.jpg"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("All files", "*.*")
            ],
            title="Export Plot"
        )

        if not file_path:
            return

        try:
            src_ext = os.path.splitext(self.current_plot_path)[1].lower()
            dst_ext = os.path.splitext(file_path)[1].lower()

            # Raster -> raster (copy or re-save)
            if src_ext in (".png", ".jpg", ".jpeg"):
                if src_ext == dst_ext:
                    shutil.copy2(self.current_plot_path, file_path)
                else:
                    # Re-save the currently displayed figure content in requested format
                    self.plot_figure.savefig(file_path, dpi=300, bbox_inches="tight")

            # Vector -> vector (copy only). No conversion in this viewer.
            elif src_ext in (".pdf", ".svg"):
                if src_ext == dst_ext:
                    shutil.copy2(self.current_plot_path, file_path)
                else:
                    messagebox.showwarning(
                        "Export Limitation",
                        "Converting PDF/SVG to another format isn’t supported here.\n"
                        "Choose the same extension to copy the original file."
                    )
                    return

            messagebox.showinfo("Export Successful", f"Plot exported to:\n{file_path}")

        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not export plot:\n{e}")

    def _open_in_os(self):
        """Open the currently selected plot in the system viewer."""
        if not self.current_plot_path:
            messagebox.showwarning("No Plot Selected", "Please select a plot first.")
            return
        try:
            if os.name == "nt":
                os.startfile(self.current_plot_path)  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.run(["open", self.current_plot_path], check=False)
            else:
                subprocess.run(["xdg-open", self.current_plot_path], check=False)
        except Exception as e:
            messagebox.showerror("Open Failed", f"Could not open file:\n{e}")

    # ---------- Helpers ----------

    def _set_info_text(self, text: str):
        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)
        self.plot_info_text.insert(tk.END, text)
        self.plot_info_text.config(state=tk.DISABLED)

    def _update_info_for_file(self, filename: str, full_path: str):
        try:
            file_size_kb = os.path.getsize(full_path) / 1024.0
            mod_time = datetime.fromtimestamp(os.path.getmtime(full_path))
            info_text = [
                f"File: {filename}",
                f"Size: {file_size_kb:.1f} KB",
                f"Created: {mod_time.strftime('%Y-%m-%d %H:%M:%S')}",
                self._infer_plot_type(filename)
            ]
            self._set_info_text("\n".join([line for line in info_text if line]))
        except Exception as e:
            self._set_info_text(f"Error reading file info: {e}")

    def _infer_plot_type(self, filename: str) -> str:
        name = filename.lower()
        if "clusters" in name or "clustercommunication" in name:
            return "Type: Cluster Communication Analysis"
        if "constellationoverview" in name or "constellation" in name:
            return "Type: Constellation Overview"
        if "intersatellitedistances" in name or "distance" in name:
            return "Type: Inter-Satellite Distance Analysis"
        if "real_" in name:
            return "Type: Real Fault Simulation Plot"
        if "friction" in name:
            return "Type: Friction Fault Analysis"
        if "encoder" in name:
            return "Type: Encoder Fault Analysis"
        if "battery" in name:
            return "Type: Battery Fault Analysis"
        if "power" in name:
            return "Type: Power-Limit Fault Analysis"
        return ""  # Unknown/unspecified

    def _show_placeholder(self, text: str):
        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, fontsize=14)
        ax.axis("off")
        self.plot_canvas.draw()

    def _show_error(self, text: str):
        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, fontsize=12, color="red")
        ax.axis("off")
        self.plot_canvas.draw()
