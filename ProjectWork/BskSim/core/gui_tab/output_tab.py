#!/usr/bin/env python
"""
output_tab.py

Implements the Output Settings tab for the Spacecraft Constellation Fault Simulator.
"""
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
import os
import datetime
from .base_tab import BaseTab

class OutputTab(BaseTab):
    """Output Settings tab implementation"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the output tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.output_dir = parent_app.output_dir
        self.plots_dir = parent_app.plots_dir
        self.bin_dir = parent_app.bin_dir
        self.simulation_time = parent_app.simulation_time
        self.show_plots = parent_app.show_plots
        self.save_plots = parent_app.save_plots
        self.save_binary = parent_app.save_binary
        self.binary_filename = parent_app.binary_filename
        
        # Create the tab UI
        self.create_tab_ui()
        
    def create_tab_ui(self):
        """Create the Output Settings tab UI"""
        # Create main container with proper scrolling
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Top section for settings
        settings_container = ttk.Frame(main_container)
        settings_container.pack(fill=tk.X, pady=(0, 10))
        
        # Simulation settings
        sim_frame = ttk.LabelFrame(settings_container, text="Simulation Settings", padding=10)
        sim_frame.pack(fill=tk.X, pady=5)
        
        # Simulation time with proper 30-minute default and presets
        time_frame = ttk.Frame(sim_frame)
        time_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(time_frame, text="Simulation Time (minutes):").pack(side=tk.LEFT)
        time_entry = ttk.Entry(time_frame, textvariable=self.simulation_time, width=10)
        time_entry.pack(side=tk.LEFT, padx=5)
        
        # Add time presets for convenience
        preset_frame = ttk.Frame(time_frame)
        preset_frame.pack(side=tk.LEFT, padx=10)
        
        ttk.Label(preset_frame, text="Quick Presets:").pack(side=tk.LEFT)
        
        preset_10min_btn = ttk.Button(preset_frame, text="10 min", 
                                     command=lambda: self.simulation_time.set(10.0))
        preset_10min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_30min_btn = ttk.Button(preset_frame, text="30 min", 
                                     command=lambda: self.simulation_time.set(30.0))
        preset_30min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_60min_btn = ttk.Button(preset_frame, text="60 min", 
                                     command=lambda: self.simulation_time.set(60.0))
        preset_60min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_90min_btn = ttk.Button(preset_frame, text="90 min", 
                                     command=lambda: self.simulation_time.set(90.0))
        preset_90min_btn.pack(side=tk.LEFT, padx=2)
        
        # Output options
        output_options_frame = ttk.LabelFrame(settings_container, text="Output Options", padding=10)
        output_options_frame.pack(fill=tk.X, pady=5)
        
        # Show plots checkbox
        plots_check = ttk.Checkbutton(output_options_frame, text="Show Plots (display during simulation)", 
                                    variable=self.show_plots)
        plots_check.pack(anchor=tk.W, pady=2)
        
        # Save plots checkbox
        save_plots_check = ttk.Checkbutton(output_options_frame, text="Save Plots (to plots directory)", 
                                        variable=self.save_plots)
        save_plots_check.pack(anchor=tk.W, pady=2)
        
        # Save binary checkbox
        binary_check = ttk.Checkbutton(output_options_frame, text="Save Binary for Vizard (3D visualization)", 
                                    variable=self.save_binary)
        binary_check.pack(anchor=tk.W, pady=2)
        
        # Binary filename
        filename_frame = ttk.Frame(output_options_frame)
        filename_frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(filename_frame, text="Binary Filename:").pack(side=tk.LEFT)
        self.binary_filename_var = tk.StringVar(value=self.binary_filename)
        filename_entry = ttk.Entry(filename_frame, textvariable=self.binary_filename_var)
        filename_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        # Output directories section
        directories_frame = ttk.LabelFrame(settings_container, text="Output Directories", padding=10)
        directories_frame.pack(fill=tk.X, pady=5)
        
        # Output directory
        dir_frame = ttk.LabelFrame(directories_frame, text="Logs Directory", padding=5)
        dir_frame.pack(fill=tk.X, pady=2)
        
        self.output_dir_var = tk.StringVar(value=self.output_dir)
        dir_entry = ttk.Entry(dir_frame, textvariable=self.output_dir_var)
        dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        browse_btn = ttk.Button(dir_frame, text="Browse...", 
                            command=self.browse_output_dir)
        browse_btn.pack(side=tk.RIGHT)
        
        # Plot directory
        plots_dir_frame = ttk.LabelFrame(directories_frame, text="Plots Directory", padding=5)
        plots_dir_frame.pack(fill=tk.X, pady=2)
        
        self.plots_dir_var = tk.StringVar(value=self.plots_dir)
        plots_dir_entry = ttk.Entry(plots_dir_frame, textvariable=self.plots_dir_var)
        plots_dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        plots_browse_btn = ttk.Button(plots_dir_frame, text="Browse...", 
                                    command=self.browse_plots_dir)
        plots_browse_btn.pack(side=tk.RIGHT)
        
        # Binary file directory
        bin_dir_frame = ttk.LabelFrame(directories_frame, text="Vizard Files Directory", padding=5)
        bin_dir_frame.pack(fill=tk.X, pady=2)
        
        self.bin_dir_var = tk.StringVar(value=self.bin_dir)
        bin_dir_entry = ttk.Entry(bin_dir_frame, textvariable=self.bin_dir_var)
        bin_dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        bin_browse_btn = ttk.Button(bin_dir_frame, text="Browse...", 
                                command=self.browse_bin_dir)
        bin_browse_btn.pack(side=tk.RIGHT)
        
        # Directory reset button
        reset_dirs_btn = ttk.Button(directories_frame, text="Reset to Defaults", 
                                   command=self.reset_directories)
        reset_dirs_btn.pack(pady=5)
        
        # Simulation log - FIXED: Now properly sized and positioned
        log_frame = ttk.LabelFrame(main_container, text="Simulation Log", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        
        # Create text widget with scrollbar
        self.log_text = scrolledtext.ScrolledText(
            log_frame,
            height=12,  # Fixed height
            wrap=tk.WORD,
            font=('Consolas', 9),
            state=tk.DISABLED
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
        # Log control buttons
        log_controls_frame = ttk.Frame(log_frame)
        log_controls_frame.pack(fill=tk.X, pady=(5, 0))
        
        clear_btn = ttk.Button(log_controls_frame, text="Clear Log", 
                            command=self.clear_log)
        clear_btn.pack(side=tk.LEFT, padx=5)
        
        save_log_btn = ttk.Button(log_controls_frame, text="Save Log", 
                                 command=self.save_log)
        save_log_btn.pack(side=tk.LEFT, padx=5)
        
        # Auto-scroll checkbox
        self.auto_scroll_var = tk.BooleanVar(value=True)
        auto_scroll_check = ttk.Checkbutton(log_controls_frame, text="Auto-scroll", 
                                          variable=self.auto_scroll_var)
        auto_scroll_check.pack(side=tk.RIGHT, padx=5)
        
        # Add initial log entry
        self.add_initial_log_entry()
        
    def add_initial_log_entry(self):
        """Add an initial log entry to show the log is working"""
        self.log_text.config(state=tk.NORMAL)
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        initial_text = f"[{timestamp}] Spacecraft Constellation Fault Simulator ready\n"
        initial_text += f"[{timestamp}] Configure settings and run simulation\n"
        self.log_text.insert(tk.END, initial_text)
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)
        
    def browse_output_dir(self):
        """Open a directory browser to select output directory"""
        directory = filedialog.askdirectory(initialdir=self.output_dir_var.get())
        if directory:
            self.output_dir_var.set(directory)
            self.parent_app.output_dir = directory
            
    def browse_plots_dir(self):
        """Open a directory browser to select plots directory"""
        directory = filedialog.askdirectory(initialdir=self.plots_dir_var.get())
        if directory:
            self.plots_dir_var.set(directory)
            self.parent_app.plots_dir = directory
            
    def browse_bin_dir(self):
        """Open a directory browser to select binary files directory"""
        directory = filedialog.askdirectory(initialdir=self.bin_dir_var.get())
        if directory:
            self.bin_dir_var.set(directory)
            self.parent_app.bin_dir = directory
    
    def reset_directories(self):
        """Reset directories to default values"""
        # Reset to project structure defaults
        root_dir = self.parent_app.ROOT_DIR
        
        self.output_dir_var.set(os.path.join(root_dir, "logs"))
        self.plots_dir_var.set(os.path.join(root_dir, "plots"))
        self.bin_dir_var.set(os.path.join(root_dir, "Vizfile"))
        
        # Update parent app
        self.parent_app.output_dir = self.output_dir_var.get()
        self.parent_app.plots_dir = self.plots_dir_var.get()
        self.parent_app.bin_dir = self.bin_dir_var.get()
        
        self.add_log_entry("Reset output directories to defaults")
    
    def add_log_entry(self, message):
        """Add an entry to the simulation log"""
        if hasattr(self, 'log_text'):
            self.log_text.config(state=tk.NORMAL)
            timestamp = datetime.datetime.now().strftime("%H:%M:%S")
            log_message = f"[{timestamp}] {message}\n"
            self.log_text.insert(tk.END, log_message)
            
            if self.auto_scroll_var.get():
                self.log_text.see(tk.END)
            
            self.log_text.config(state=tk.DISABLED)
            
    def clear_log(self):
        """Clear the log text"""
        if self.log_text:
            self.log_text.config(state=tk.NORMAL)
            self.log_text.delete(1.0, tk.END)
            self.log_text.config(state=tk.DISABLED)
            # Add a cleared message
            self.add_log_entry("Log cleared")
    
    def save_log(self):
        """Save the current log to a file"""
        if not self.log_text:
            return
            
        # Ask user for save location
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("Log files", "*.log"),
                ("All files", "*.*")
            ],
            title="Save Simulation Log"
        )
        
        if file_path:
            try:
                # Get log content
                self.log_text.config(state=tk.NORMAL)
                log_content = self.log_text.get(1.0, tk.END)
                self.log_text.config(state=tk.DISABLED)
                
                # Write to file
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(f"Spacecraft Fault Simulator Log\n")
                    f.write(f"Saved: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(log_content)
                
                self.show_message("Save Successful", f"Log saved to:\n{file_path}")
                self.add_log_entry(f"Log saved to {os.path.basename(file_path)}")
                
            except Exception as e:
                self.show_message("Save Failed", f"Could not save log:\n{str(e)}", "error")
            
    def update_settings_from_ui(self):
        """Update parent app settings from UI values"""
        # Update settings from UI values
        self.parent_app.binary_filename = self.binary_filename_var.get()
        self.parent_app.output_dir = self.output_dir_var.get()
        self.parent_app.plots_dir = self.plots_dir_var.get()
        self.parent_app.bin_dir = self.bin_dir_var.get()
        
        # Create directories if they don't exist
        for directory in [self.parent_app.output_dir, self.parent_app.plots_dir, self.parent_app.bin_dir]:
            try:
                os.makedirs(directory, exist_ok=True)
            except Exception as e:
                self.add_log_entry(f"Warning: Could not create directory {directory}: {e}")