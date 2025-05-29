#!/usr/bin/env python
"""
results_tab.py

Results viewing and analysis tab for the Spacecraft Constellation Fault Simulator.
"""
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from datetime import datetime
import matplotlib.pyplot as plt
from PIL import Image

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

        # Main content area
        split_frame = ttk.Frame(main_frame)
        split_frame.pack(fill=tk.BOTH, expand=True)

        # Left panel - plot list
        list_frame = ttk.LabelFrame(split_frame, text="Available Plots", padding=10, width=320)
        list_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        list_frame.pack_propagate(False)  # Maintain width

        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.plot_listbox = tk.Listbox(list_container, selectmode=tk.SINGLE, 
                                      yscrollcommand=scrollbar.set, font=('Segoe UI', 10), 
                                      exportselection=False)
        self.plot_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.plot_listbox.yview)

        self.plot_listbox.bind('<<ListboxSelect>>', self.on_plot_selected)

        # Right panel - plot display
        display_frame = ttk.LabelFrame(split_frame, text="Plot View", padding=10)
        display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))

        # Create a frame for the canvas and toolbar
        plot_container = ttk.Frame(display_frame)
        plot_container.pack(fill=tk.BOTH, expand=True)

        # Create the matplotlib figure and canvas
        self.plot_figure = Figure(figsize=(10, 6), dpi=100)
        self.plot_canvas = FigureCanvasTkAgg(self.plot_figure, plot_container)
        
        # IMPORTANT: Create toolbar BEFORE packing canvas
        toolbar_frame = ttk.Frame(plot_container)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        self.plot_toolbar = NavigationToolbar2Tk(self.plot_canvas, toolbar_frame)
        self.plot_toolbar.update()
        
        # Now pack the canvas
        self.plot_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Bottom info panel
        info_frame = ttk.LabelFrame(main_frame, text="Plot Information", padding=10)
        info_frame.pack(fill=tk.X, pady=(10, 0))

        self.plot_info_text = scrolledtext.ScrolledText(info_frame, height=4, wrap=tk.WORD)
        self.plot_info_text.pack(fill=tk.BOTH, expand=True)
        self.plot_info_text.config(state=tk.DISABLED)

        # Initial refresh
        self.refresh_plot_list()

    def refresh_plot_list(self):
        """Refresh the list of available plots"""
        self.plot_listbox.delete(0, tk.END)

        # Clear info text
        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)
        self.plot_info_text.insert(tk.END, "Select a plot to view details and analysis")
        self.plot_info_text.config(state=tk.DISABLED)

        # Clear plot display
        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        ax.text(0.5, 0.5, "No plot selected\nSelect a plot from the list to view", 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.axis('off')
        self.plot_canvas.draw()

        if not os.path.exists(self.app.plots_dir):
            return

        # Get plot files
        plot_files = [f for f in os.listdir(self.app.plots_dir) 
                     if f.endswith(('.png', '.jpg', '.jpeg', '.pdf', '.svg'))]
        plot_files.sort(key=lambda x: os.path.getmtime(os.path.join(self.app.plots_dir, x)), 
                       reverse=True)

        # Add to listbox
        for filename in plot_files:
            self.plot_listbox.insert(tk.END, filename)

        # Select first plot if available
        if plot_files:
            self.plot_listbox.selection_set(0)
            self.on_plot_selected(None)

    def on_plot_selected(self, event):
        """Handle plot selection"""
        selection = self.plot_listbox.curselection()
        if not selection:
            return

        index = selection[0]
        filename = self.plot_listbox.get(index)
        full_path = os.path.join(self.app.plots_dir, filename)
        self.current_plot_path = full_path

        # Update info text
        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)

        try:
            file_size = os.path.getsize(full_path) / 1024
            mod_time = datetime.fromtimestamp(os.path.getmtime(full_path))
            info_text = f"File: {filename}\n"
            info_text += f"Size: {file_size:.1f} KB\n"
            info_text += f"Created: {mod_time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            
            # Add plot type info
            if "FrictionFault" in filename:
                info_text += "Type: Friction Fault Analysis"
            elif "InterSatelliteDistances" in filename:
                info_text += "Type: Inter-Satellite Distance Analysis"
            elif "ConstellationOrbits" in filename:
                info_text += "Type: Constellation Orbit Visualization"
            
            self.plot_info_text.insert(tk.END, info_text)
        except Exception as e:
            self.plot_info_text.insert(tk.END, f"Error reading file info: {e}")
        
        self.plot_info_text.config(state=tk.DISABLED)

        # Display the plot
        try:
            self.plot_figure.clear()
            
            if filename.endswith(('.png', '.jpg', '.jpeg')):
                # Load and display image
                img = Image.open(full_path)
                ax = self.plot_figure.add_subplot(111)
                ax.imshow(img)
                ax.axis('off')
                
                # Add title
                title = filename.replace('_', ' ').replace('.png', '').replace('.jpg', '')
                self.plot_figure.suptitle(title, fontsize=12, y=0.98)
                
            else:
                # For non-image files
                ax = self.plot_figure.add_subplot(111)
                ax.text(0.5, 0.5, f"Plot: {filename}\nUse 'Export Plot' to save", 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.axis('off')
            
            # Refresh canvas
            self.plot_figure.tight_layout()
            self.plot_canvas.draw()
            
        except Exception as e:
            self.plot_figure.clear()
            ax = self.plot_figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Error displaying plot:\n{str(e)}", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12, color='red')
            ax.axis('off')
            self.plot_canvas.draw()

    def export_current_plot(self):
        """Export the currently selected plot"""
        if not self.current_plot_path:
            messagebox.showwarning("No Plot Selected", "Please select a plot first.")
            return

        # Get save location
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
        
        if file_path:
            try:
                # If it's the same format, just copy
                if self.current_plot_path.endswith(os.path.splitext(file_path)[1]):
                    import shutil
                    shutil.copy2(self.current_plot_path, file_path)
                else:
                    # Otherwise, save the figure in the requested format
                    self.plot_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                
                messagebox.showinfo("Export Successful", f"Plot exported to:\n{file_path}")
                
            except Exception as e:
                messagebox.showerror("Export Failed", f"Could not export plot:\n{str(e)}")