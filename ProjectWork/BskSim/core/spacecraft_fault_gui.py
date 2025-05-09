#!/usr/bin/env python
"""
spacecraft_fault_gui.py with modern, styled Tkinter UI.
"""
import tkinter as tk
from tkinter import ttk, messagebox, colorchooser, filedialog
from PIL import ImageTk, Image
import sys
import os
import json

# Fix path resolution to work with new project structure
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CORE_DIR = os.path.join(ROOT_DIR, 'core')
PLOTS_DIR = os.path.join(ROOT_DIR, 'plots')

sys.path.extend([ROOT_DIR, CORE_DIR])

# Import CLI core functionality
from core.spacecraft_fault_cli import SimulationConfig, TargetDefinition, run_custom_simulation

class FaultSimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Spacecraft RW Fault Simulator")
        self.root.geometry("1000x750")
        self.config = SimulationConfig()
        self.setup_style()
        self.create_widgets()

    def setup_style(self):
        style = ttk.Style(self.root)
        style.theme_use('clam')
        style.configure('TNotebook.Tab', padding=(12, 8), font=('Segoe UI', 10, 'bold'))
        style.configure('TFrame', background='#f4f4f4')
        style.configure('TLabel', background='#f4f4f4', font=('Segoe UI', 10))
        style.configure('TButton', font=('Segoe UI', 10))
        style.configure('Header.TLabel', font=('Segoe UI', 12, 'bold'))
        style.configure('Treeview.Heading', font=('Segoe UI', 10, 'bold'))

    def create_widgets(self):
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True)

        # Create each tab
        self.sim_tab = ttk.Frame(notebook)
        self.fault_tab = ttk.Frame(notebook)
        self.viz_tab = ttk.Frame(notebook)

        notebook.add(self.sim_tab, text="Simulation Settings")
        notebook.add(self.fault_tab, text="Fault Configuration")
        notebook.add(self.viz_tab, text="Visualization")

        self.create_simulation_tab()
        self.create_fault_tab()
        self.create_viz_tab()

        # Bottom controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=10)

        self.status_label = ttk.Label(control_frame, text="Ready", style='Header.TLabel')
        self.status_label.pack(side=tk.LEFT, padx=10, fill=tk.X, expand=True)

        ttk.Button(control_frame, text="Export Config", command=self.export_config).pack(side=tk.RIGHT, padx=5)
        ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation).pack(side=tk.RIGHT, padx=5)

    def create_simulation_tab(self):
        f = self.sim_tab

        ttk.Label(f, text="Simulation Time (min):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.sim_time = tk.DoubleVar(value=self.config.simulation_time)
        ttk.Entry(f, textvariable=self.sim_time, width=10).grid(row=0, column=1, padx=5, pady=5)

        out_frame = ttk.LabelFrame(f, text="Output")
        out_frame.grid(row=1, column=0, columnspan=2, sticky=tk.EW, padx=5, pady=10)
        self.show_plots = tk.BooleanVar(value=self.config.show_plots)
        self.save_binary = tk.BooleanVar(value=self.config.save_binary)
        self.binary_filename = tk.StringVar(value=self.config.binary_filename)

        ttk.Checkbutton(out_frame, text="Show Plots", variable=self.show_plots).grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Checkbutton(out_frame, text="Save Binary", variable=self.save_binary).grid(row=0, column=1, sticky=tk.W, padx=5)

        ttk.Label(out_frame, text="Filename:").grid(row=1, column=0, sticky=tk.W, padx=5)
        ttk.Entry(out_frame, textvariable=self.binary_filename, width=20).grid(row=1, column=1, sticky=tk.W, padx=5)

    def create_fault_tab(self):
        f = self.fault_tab
        self.fault_mag = tk.DoubleVar(value=self.config.fault_magnitude)
        self.fault_time = tk.DoubleVar(value=self.config.fault_time)
        self.fault_wheel = tk.IntVar(value=self.config.fault_wheel_number)

        fault_frame = ttk.LabelFrame(f, text="One-Time Fault")
        fault_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(fault_frame, text="Magnitude:").grid(row=0, column=0, padx=5, pady=5)
        ttk.Entry(fault_frame, textvariable=self.fault_mag, width=10).grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(fault_frame, text="Wheel (0-3):").grid(row=0, column=2, padx=5, pady=5)
        ttk.Spinbox(fault_frame, from_=0, to=3, textvariable=self.fault_wheel, width=4).grid(row=0, column=3, padx=5, pady=5)

        ttk.Label(fault_frame, text="Time (min):").grid(row=0, column=4, padx=5, pady=5)
        ttk.Entry(fault_frame, textvariable=self.fault_time, width=10).grid(row=0, column=5, padx=5, pady=5)

    def create_viz_tab(self):
        f = self.viz_tab

        ttk.Label(f, text="Target Management").pack(anchor=tk.W, padx=10, pady=5)

        self.target_list = ttk.Treeview(f, columns=('name', 'lat', 'lon', 'color'), show='headings', height=10)
        for col in ('name', 'lat', 'lon', 'color'):
            self.target_list.heading(col, text=col.capitalize())
            self.target_list.column(col, width=100, anchor=tk.CENTER)
        self.target_list.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        for target in self.config.targets:
            self.target_list.insert('', tk.END, values=(target.name, target.latitude, target.longitude, target.color))

    def run_simulation(self):
        try:
            self.config.simulation_time = self.sim_time.get()
            self.config.show_plots = self.show_plots.get()
            self.config.save_binary = self.save_binary.get()
            self.config.binary_filename = self.binary_filename.get()
            self.config.fault_magnitude = self.fault_mag.get()
            self.config.fault_time = self.fault_time.get()
            self.config.fault_wheel_number = self.fault_wheel.get()
            self.config.targets = []
            for item in self.target_list.get_children():
                vals = self.target_list.item(item)['values']
                self.config.targets.append(TargetDefinition(vals[0], float(vals[1]), float(vals[2]), vals[3]))

            self.status_label.config(text="Running...")
            self.root.update()

            scenario, viz, figureList, output_dir = run_custom_simulation(self.config)

            for name, fig in figureList.items():
                fig_path = os.path.join(PLOTS_DIR, f"{name}.png")
                fig.savefig(fig_path)

            self.status_label.config(text="Completed")
            messagebox.showinfo("Simulation Complete", f"Results saved to {output_dir}\nPlots saved to {PLOTS_DIR}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text="Error")

    def export_config(self):
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".json",
                                                     filetypes=[("JSON files", "*.json")],
                                                     title="Save Configuration")
            if file_path:
                config_dict = {
                    "simulation_time": self.sim_time.get(),
                    "show_plots": self.show_plots.get(),
                    "save_binary": self.save_binary.get(),
                    "binary_filename": self.binary_filename.get(),
                    "fault_magnitude": self.fault_mag.get(),
                    "fault_time": self.fault_time.get(),
                    "fault_wheel_number": self.fault_wheel.get(),
                    "targets": [self.target_list.item(i)['values'] for i in self.target_list.get_children()]
                }
                with open(file_path, 'w') as f:
                    json.dump(config_dict, f, indent=4)
                messagebox.showinfo("Export Success", f"Configuration saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = FaultSimulationGUI(root)
    root.mainloop()
