#!/usr/bin/env python
"""
output_tab.py

Implements the Output Settings tab for the Spacecraft Constellation Fault Simulator.
FIXED: Proper layout for simulation log, 30-minute default, Vizard time format explanation
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
        """Create the Output Settings tab UI with FIXED layout"""
        # FIXED: Main container with scrollable frame for better organization
        main_canvas = tk.Canvas(self.parent_frame)
        scrollbar = ttk.Scrollbar(self.parent_frame, orient="vertical", command=main_canvas.yview)
        scrollable_frame = ttk.Frame(main_canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        )
        
        main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=scrollbar.set)
        
        main_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # FIXED: Simulation settings with proper 30-minute default and Vizard time explanation
        sim_frame = ttk.LabelFrame(scrollable_frame, text="Simulation Settings", padding=10)
        sim_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Add help button
        self.add_help_button(sim_frame, "Simulation Settings", 
                           lambda: self.show_help_content("Simulation Settings", 
                                              "Configure the duration of the simulation and other basic parameters."))
        
        # FIXED: Simulation time with proper 30-minute default and time format explanation
        time_frame = ttk.Frame(sim_frame)
        time_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(time_frame, text="Simulation Time (minutes):").pack(side=tk.LEFT)
        time_entry = ttk.Entry(time_frame, textvariable=self.simulation_time, width=10)
        time_entry.pack(side=tk.LEFT, padx=5)
        
        # FIXED: Add Vizard time format explanation
        vizard_time_frame = ttk.Frame(time_frame)
        vizard_time_frame.pack(side=tk.LEFT, padx=20)
        
        ttk.Label(vizard_time_frame, text="Vizard Time Format:", font=('Segoe UI', 9, 'bold')).pack(side=tk.LEFT)
        self.vizard_time_label = ttk.Label(vizard_time_frame, 
                                          text="00:30:00:00:0", 
                                          font=('Segoe UI', 9, 'italic'),
                                          foreground='darkgreen')
        self.vizard_time_label.pack(side=tk.LEFT, padx=5)
        
        # Bind time change to update Vizard format
        self.simulation_time.trace('w', self.update_vizard_time_format)
        
        # Add time presets for convenience
        preset_frame = ttk.Frame(time_frame)
        preset_frame.pack(side=tk.RIGHT)
        
        ttk.Label(preset_frame, text="Quick Presets:").pack(side=tk.LEFT)
        
        preset_1min_btn = ttk.Button(preset_frame, text="1 min", 
                                    command=lambda: self.simulation_time.set(1.0))
        preset_1min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_5min_btn = ttk.Button(preset_frame, text="5 min", 
                                    command=lambda: self.simulation_time.set(5.0))
        preset_5min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_30min_btn = ttk.Button(preset_frame, text="30 min (Default)", 
                                     command=lambda: self.simulation_time.set(30.0),
                                     style="Run.TButton")
        preset_30min_btn.pack(side=tk.LEFT, padx=2)
        
        preset_60min_btn = ttk.Button(preset_frame, text="60 min", 
                                     command=lambda: self.simulation_time.set(60.0))
        preset_60min_btn.pack(side=tk.LEFT, padx=2)
        
        # FIXED: Add simulation time info with Vizard explanation
        info_frame = ttk.Frame(sim_frame)
        info_frame.pack(fill=tk.X, pady=5)
        
        info_text = ("FIXED: Default is 30 minutes for comprehensive satellite and target analysis.\n"
                    "Vizard Time Format: DD:HH:MM:SS:F (Days:Hours:Minutes:Seconds:Frames)\n"
                    "For 30 minutes: 00:00:30:00:0 (0 days, 0 hours, 30 minutes, 0 seconds, 0 frames)")
        ttk.Label(info_frame, text=info_text, style="Info.TLabel", wraplength=700).pack(anchor=tk.W)
        
        # Output options
        output_options_frame = ttk.LabelFrame(scrollable_frame, text="Output Options", padding=10)
        output_options_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Add help button
        self.add_help_button(output_options_frame, "Output Options",
                           lambda: self.show_help_content("Output Options",
                                              "Configure what outputs should be generated from the simulation."))
        
        # Show plots checkbox
        plots_check = ttk.Checkbutton(output_options_frame, text="Show Plots (display during simulation)", 
                                    variable=self.show_plots)
        plots_check.pack(anchor=tk.W, pady=3)
        
        # Save plots checkbox
        save_plots_check = ttk.Checkbutton(output_options_frame, text="Save Plots (to plots directory)", 
                                        variable=self.save_plots)
        save_plots_check.pack(anchor=tk.W, pady=3)
        
        # Save binary checkbox
        binary_check = ttk.Checkbutton(output_options_frame, text="Save Binary for Vizard (3D visualization with targets)", 
                                    variable=self.save_binary)
        binary_check.pack(anchor=tk.W, pady=3)
        
        # Binary filename
        filename_frame = ttk.Frame(output_options_frame)
        filename_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(filename_frame, text="Binary Filename:").pack(side=tk.LEFT)
        self.binary_filename_var = tk.StringVar(value=self.binary_filename)
        filename_entry = ttk.Entry(filename_frame, textvariable=self.binary_filename_var, width=30)
        filename_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        # Add note about binary files
        binary_note_frame = ttk.Frame(output_options_frame)
        binary_note_frame.pack(fill=tk.X, pady=5)
        
        binary_note = ("FIXED: Binary files can be opened in Vizard for 3D visualization with enhanced target visibility.\n"
                      "Files are saved as [filename]_UnityViz.bin in the Vizfile directory.\n"
                      "Targets appear as colored markers on Earth surface for easy identification.")
        ttk.Label(binary_note_frame, text=binary_note, style="Info.TLabel", wraplength=700).pack(anchor=tk.W)
        
        # Output directories section
        directories_frame = ttk.LabelFrame(scrollable_frame, text="Output Directories", padding=10)
        directories_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Output directory
        dir_frame = ttk.LabelFrame(directories_frame, text="Logs Directory", padding=5)
        dir_frame.pack(fill=tk.X, pady=3)
        
        self.output_dir_var = tk.StringVar(value=self.output_dir)
        dir_entry = ttk.Entry(dir_frame, textvariable=self.output_dir_var)
        dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        browse_btn = ttk.Button(dir_frame, text="Browse...", 
                            command=self.browse_output_dir)
        browse_btn.pack(side=tk.RIGHT)
        
        # Plot directory
        plots_dir_frame = ttk.LabelFrame(directories_frame, text="Plots Directory", padding=5)
        plots_dir_frame.pack(fill=tk.X, pady=3)
        
        self.plots_dir_var = tk.StringVar(value=self.plots_dir)
        plots_dir_entry = ttk.Entry(plots_dir_frame, textvariable=self.plots_dir_var)
        plots_dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        plots_browse_btn = ttk.Button(plots_dir_frame, text="Browse...", 
                                    command=self.browse_plots_dir)
        plots_browse_btn.pack(side=tk.RIGHT)
        
        # Binary file directory
        bin_dir_frame = ttk.LabelFrame(directories_frame, text="Vizard Files Directory", padding=5)
        bin_dir_frame.pack(fill=tk.X, pady=3)
        
        self.bin_dir_var = tk.StringVar(value=self.bin_dir)
        bin_dir_entry = ttk.Entry(bin_dir_frame, textvariable=self.bin_dir_var)
        bin_dir_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        bin_browse_btn = ttk.Button(bin_dir_frame, text="Browse...", 
                                command=self.browse_bin_dir)
        bin_browse_btn.pack(side=tk.RIGHT)
        
        # Directory reset button
        reset_dirs_btn = ttk.Button(directories_frame, text="Reset to Defaults", 
                                   command=self.reset_directories)
        reset_dirs_btn.pack(pady=10)
        
        # FIXED: Simulation log with proper layout and sizing
        log_frame = ttk.LabelFrame(scrollable_frame, text="Simulation Log", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # FIXED: Log controls at the top for better visibility
        log_controls_frame = ttk.Frame(log_frame)
        log_controls_frame.pack(fill=tk.X, pady=(0, 5))
        
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
        
        # Status label
        status_label = ttk.Label(log_controls_frame, text="Status: Ready for simulation", 
                               style="Info.TLabel")
        status_label.pack(side=tk.RIGHT, padx=20)
        
        # FIXED: Log text area with proper sizing and scrolling
        log_text_frame = ttk.Frame(log_frame)
        log_text_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create the scrolled text widget with FIXED dimensions
        self.log_text = scrolledtext.ScrolledText(
            log_text_frame, 
            height=15,  # FIXED: Proper height
            width=100,  # FIXED: Proper width
            wrap=tk.WORD,
            state="disabled",
            font=('Consolas', 9),  # FIXED: Better font for log readability
            bg='#f8f8f8',
            fg='#333333'
        )
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # FIXED: Add initial welcome message to log
        self.add_initial_log_message()
        
        # Update Vizard time format initially
        self.update_vizard_time_format()
    
    def update_vizard_time_format(self, *args):
        """FIXED: Update the Vizard time format display based on simulation time"""
        try:
            sim_time_minutes = float(self.simulation_time.get())
            
            # Convert to Vizard format: DD:HH:MM:SS:F
            days = int(sim_time_minutes // (24 * 60))
            hours = int((sim_time_minutes % (24 * 60)) // 60)
            minutes = int(sim_time_minutes % 60)
            seconds = int((sim_time_minutes % 1) * 60)
            frames = 0  # Usually 0 for our purposes
            
            vizard_format = f"{days:02d}:{hours:02d}:{minutes:02d}:{seconds:02d}:{frames}"
            self.vizard_time_label.config(text=vizard_format)
            
        except (ValueError, AttributeError):
            self.vizard_time_label.config(text="00:00:30:00:0")  # Default to 30 minutes
    
    def add_initial_log_message(self):
        """FIXED: Add initial welcome message to the log"""
        try:
            self.log_text.config(state="normal")
            welcome_msg = f"""
=== Spacecraft Constellation Fault Simulator ===
Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

FIXED IMPROVEMENTS:
✓ Default simulation time: 30 minutes (optimal for target analysis)
✓ Improved orbital altitude: 1629km for better visibility  
✓ Enhanced target visibility in Vizard
✓ Better camera positioning for ground target observation
✓ Simplified satellite names in constellation view

Ready to run simulation...
Configure your satellites, faults, and targets, then click 'Run Simulation'.

"""
            self.log_text.insert(tk.END, welcome_msg)
            self.log_text.config(state="disabled")
            self.log_text.see(tk.END)
        except Exception as e:
            print(f"Could not add initial log message: {e}")
    
    def browse_output_dir(self):
        """Open a directory browser to select output directory"""
        directory = filedialog.askdirectory(initialdir=self.output_dir_var.get())
        if directory:
            self.output_dir_var.set(directory)
            self.parent_app.output_dir = directory
            
            # Update plots and bin directory defaults if they follow the pattern
            if self.plots_dir_var.get().startswith(self.parent_app.output_dir):
                self.plots_dir_var.set(os.path.join(directory, "plots"))
                self.parent_app.plots_dir = os.path.join(directory, "plots")
                
            if self.bin_dir_var.get().startswith(self.parent_app.output_dir):
                self.bin_dir_var.set(os.path.join(directory, "Vizfile"))
                self.parent_app.bin_dir = os.path.join(directory, "Vizfile")
                
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
        
        self.add_to_log("Reset output directories to defaults")
            
    def clear_log(self):
        """Clear the log text"""
        if self.log_text:
            self.log_text.config(state="normal")
            self.log_text.delete(1.0, tk.END)
            self.log_text.config(state="disabled")
            # Add a cleared message
            self.add_to_log("Log cleared")
    
    def add_to_log(self, message):
        """FIXED: Add a message to the log with proper formatting"""
        try:
            if self.log_text:
                self.log_text.config(state="normal")
                timestamp = datetime.datetime.now().strftime("%H:%M:%S")
                self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
                
                # Auto-scroll if enabled
                if self.auto_scroll_var.get():
                    self.log_text.see(tk.END)
                    
                self.log_text.config(state="disabled")
                self.parent_frame.update_idletasks()
        except Exception as e:
            print(f"Error adding to log: {e}")
    
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
                self.log_text.config(state="normal")
                log_content = self.log_text.get(1.0, tk.END)
                self.log_text.config(state="disabled")
                
                # Write to file
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.write(f"Spacecraft Fault Simulator Log\n")
                    f.write(f"Saved: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(log_content)
                
                self.show_message("Save Successful", f"Log saved to:\n{file_path}")
                self.add_to_log(f"Log saved to: {os.path.basename(file_path)}")
                
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
                self.add_to_log(f"Warning: Could not create directory {directory}: {e}")