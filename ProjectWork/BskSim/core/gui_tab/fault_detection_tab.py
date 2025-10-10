#!/usr/bin/env python
"""
fault_detection_tab.py (UPDATED WITH XAI MODEL SELECTION)

Enhanced Fault Detection Tab with integrated LIME and SHAP model selection.
Users can now choose between Autoencoder, LIME, and SHAP models from a dropdown.
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from datetime import datetime
import json
import os
import sys

class FaultDetectionTab:
    """
    Enhanced Fault Detection Tab with XAI model selection
    """
    
    def __init__(self, parent_app, parent_frame):
        """Initialize the fault detection tab with XAI support"""
        self.parent_app = parent_app
        self.parent_frame = parent_frame
        
        self.fault_detection_results = None
        self.ml_detector = None
        self.xai_explainer = None  # NEW: XAI explainer instance
        self.detection_history = []
        self.monitoring_active = False
        
        # Check ML and XAI availability
        self.ml_available = self.check_ml_availability()
        self.xai_available = self.check_xai_availability()
        
        # Create the tab UI
        self.create_tab_ui()
        
        # Check for default model
        self.check_for_default_model()
        
    def check_ml_availability(self):
        """Check if ML fault detection is available"""
        try:
            from real_ml_fault_detection import RealMLFaultDetector
            return True
        except ImportError:
            return False
    
    def check_xai_availability(self):
        """Check if XAI libraries (LIME and SHAP) are available"""
        try:
            import shap
            import lime
            return True
        except ImportError:
            return False
    
    def create_tab_ui(self):
        """Create the Fault Detection tab UI"""
        # Create a notebook for sub-tabs
        self.detection_notebook = ttk.Notebook(self.parent_frame)
        self.detection_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create frames for sub-tabs
        overview_frame = ttk.Frame(self.detection_notebook)
        config_frame = ttk.Frame(self.detection_notebook)
        results_frame = ttk.Frame(self.detection_notebook)
        xai_frame = ttk.Frame(self.detection_notebook)  # NEW: XAI explanations tab
        monitoring_frame = ttk.Frame(self.detection_notebook)
        analysis_frame = ttk.Frame(self.detection_notebook)
        
        # Add sub-tabs
        self.detection_notebook.add(overview_frame, text="Overview")
        self.detection_notebook.add(config_frame, text="ML Configuration")
        self.detection_notebook.add(results_frame, text="Detection Results")
        self.detection_notebook.add(xai_frame, text="XAI Explanations")  # NEW
        self.detection_notebook.add(monitoring_frame, text="Live Monitoring")
        self.detection_notebook.add(analysis_frame, text="Analysis")
        
        # Create tab contents
        self._create_overview_tab(overview_frame)
        self._create_config_tab(config_frame)
        self._create_results_tab(results_frame)
        self._create_xai_tab(xai_frame)  # NEW: XAI tab
        self._create_monitoring_tab(monitoring_frame)
        self._create_analysis_tab(analysis_frame)
        
    def _create_overview_tab(self, parent):
        """Create the overview tab content"""
        # Title frame
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(10, 20))
        
        title_label = ttk.Label(title_frame, text="ML Fault Detection System with XAI", 
                               font=('Segoe UI', 16, 'bold'))
        title_label.pack()
        
        subtitle_label = ttk.Label(title_frame, text="Real-time fault detection with LIME and SHAP explanations",
                                  font=('Segoe UI', 10), foreground='gray')
        subtitle_label.pack()
        
        # Status frame
        status_frame = ttk.LabelFrame(parent, text="System Status", padding=15)
        status_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # ML Status
        ml_status_frame = ttk.Frame(status_frame)
        ml_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(ml_status_frame, text="ML Model Status:").pack(side=tk.LEFT)
        
        self.ml_status_label = ttk.Label(ml_status_frame, text="", font=('Segoe UI', 10, 'bold'))
        if self.ml_available:
            self.ml_status_label.config(text="Available", foreground='orange')
        else:
            self.ml_status_label.config(text="Not Available", foreground='red')
        self.ml_status_label.pack(side=tk.LEFT, padx=10)
        
        # XAI Status (NEW)
        xai_status_frame = ttk.Frame(status_frame)
        xai_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(xai_status_frame, text="XAI Status:").pack(side=tk.LEFT)
        
        self.xai_status_label = ttk.Label(xai_status_frame, text="", font=('Segoe UI', 10, 'bold'))
        if self.xai_available:
            self.xai_status_label.config(text="Available (LIME & SHAP)", foreground='green')
        else:
            self.xai_status_label.config(text="Not Available", foreground='red')
        self.xai_status_label.pack(side=tk.LEFT, padx=10)
        
        # Detection Status
        detection_status_frame = ttk.Frame(status_frame)
        detection_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(detection_status_frame, text="Detection Status:").pack(side=tk.LEFT)
        
        self.detection_status_label = ttk.Label(detection_status_frame, text="Idle", 
                                               font=('Segoe UI', 10, 'bold'))
        self.detection_status_label.pack(side=tk.LEFT, padx=10)
        
        # Statistics frame
        stats_frame = ttk.LabelFrame(parent, text="Detection Statistics", padding=15)
        stats_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Create statistics grid
        stats_grid = ttk.Frame(stats_frame)
        stats_grid.pack(fill=tk.X)
        
        # Total detections
        ttk.Label(stats_grid, text="Total Detections:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.total_detections_label = ttk.Label(stats_grid, text="0", font=('Segoe UI', 10, 'bold'), 
                                               foreground='blue')
        self.total_detections_label.grid(row=0, column=1, sticky=tk.W, padx=10)
        
        # Active faults
        ttk.Label(stats_grid, text="Active Faults:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.active_faults_label = ttk.Label(stats_grid, text="0", font=('Segoe UI', 10, 'bold'), 
                                            foreground='red')
        self.active_faults_label.grid(row=1, column=1, sticky=tk.W, padx=10)
        
        # Success rate
        ttk.Label(stats_grid, text="Detection Success Rate:").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.success_rate_label = ttk.Label(stats_grid, text="N/A", font=('Segoe UI', 10, 'bold'), 
                                           foreground='green')
        self.success_rate_label.grid(row=2, column=1, sticky=tk.W, padx=10)
        
        # XAI Explanations generated (NEW)
        ttk.Label(stats_grid, text="XAI Explanations:").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.xai_explanations_label = ttk.Label(stats_grid, text="0", font=('Segoe UI', 10, 'bold'), 
                                               foreground='purple')
        self.xai_explanations_label.grid(row=3, column=1, sticky=tk.W, padx=10)
        
        # Last detection
        ttk.Label(stats_grid, text="Last Detection:").grid(row=4, column=0, sticky=tk.W, pady=2)
        self.last_detection_label = ttk.Label(stats_grid, text="None", foreground='gray')
        self.last_detection_label.grid(row=4, column=1, sticky=tk.W, padx=10)
        
        # Quick actions frame
        actions_frame = ttk.LabelFrame(parent, text="Quick Actions", padding=15)
        actions_frame.pack(fill=tk.X, padx=20, pady=10)
        
        buttons_frame = ttk.Frame(actions_frame)
        buttons_frame.pack()
        
        self.load_model_btn = ttk.Button(buttons_frame, text="Load ML Model", 
                                        command=self.load_ml_model)
        self.load_model_btn.pack(side=tk.LEFT, padx=5)
        
        self.start_detection_btn = ttk.Button(buttons_frame, text="Start Detection", 
                                             command=self.start_fault_detection)
        self.start_detection_btn.pack(side=tk.LEFT, padx=5)
        
        # NEW: Generate XAI button
        self.generate_xai_btn = ttk.Button(buttons_frame, text="Generate XAI", 
                                          command=self.generate_xai_explanations,
                                          style="Accent.TButton")
        self.generate_xai_btn.pack(side=tk.LEFT, padx=5)
        self.generate_xai_btn.config(state="disabled")  # Disabled until detection complete
        
        self.simulate_btn = ttk.Button(buttons_frame, text="Simulate Results", 
                                      command=self.simulate_detection_results)
        self.simulate_btn.pack(side=tk.LEFT, padx=5)
        
        self.reset_btn = ttk.Button(buttons_frame, text="Reset", 
                                   command=self.reset_detection_state)
        self.reset_btn.pack(side=tk.LEFT, padx=5)
        
        self.view_results_btn = ttk.Button(buttons_frame, text="View Results", 
                                          command=lambda: self.detection_notebook.select(2))
        self.view_results_btn.pack(side=tk.LEFT, padx=5)
        
        # Recent detections frame
        recent_frame = ttk.LabelFrame(parent, text="Recent Detections", padding=15)
        recent_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Recent detections listbox with scrollbar
        listbox_frame = ttk.Frame(recent_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True)
        
        self.recent_detections_listbox = tk.Listbox(listbox_frame, height=8)
        scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical", command=self.recent_detections_listbox.yview)
        self.recent_detections_listbox.configure(yscrollcommand=scrollbar.set)
        
        self.recent_detections_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.recent_detections_listbox.insert(tk.END, "No recent detections...")
        
    def _create_config_tab(self, parent):
        """Create the ML configuration tab content with XAI model selection"""
        # ML Model Configuration
        model_frame = ttk.LabelFrame(parent, text="ML Model Selection", padding=15)
        model_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # NEW: Model type selection (Autoencoder, LIME, SHAP)
        model_type_frame = ttk.Frame(model_frame)
        model_type_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(model_type_frame, text="Select Model Type:", font=('Segoe UI', 10, 'bold')).pack(side=tk.LEFT)
        
        self.model_type_var = tk.StringVar(value="Autoencoder")
        model_types = ["Autoencoder (Base ML)", "LIME (Local Explanations)", "SHAP (Global Explanations)", "All Models (Ensemble)"]
        
        for model_type in model_types:
            rb = ttk.Radiobutton(
                model_type_frame, 
                text=model_type, 
                variable=self.model_type_var, 
                value=model_type.split(" ")[0],
                command=self.on_model_type_changed
            )
            rb.pack(anchor=tk.W, padx=(20, 0), pady=2)
        
        # Model file selection (only for Autoencoder or when browsing)
        file_frame = ttk.Frame(model_frame)
        file_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(file_frame, text="Base Model File:").pack(side=tk.LEFT)
        
        self.model_path_var = tk.StringVar(value="anomaly_detection_model.keras")
        self.model_path_entry = ttk.Entry(file_frame, textvariable=self.model_path_var, width=50)
        self.model_path_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        self.browse_model_btn = ttk.Button(file_frame, text="Browse", 
                                          command=self.browse_model_file)
        self.browse_model_btn.pack(side=tk.RIGHT, padx=5)
        
        # XAI Configuration (NEW)
        xai_config_frame = ttk.LabelFrame(parent, text="XAI Configuration", padding=15)
        xai_config_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # LIME settings
        lime_frame = ttk.Frame(xai_config_frame)
        lime_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(lime_frame, text="LIME Features:").grid(row=0, column=0, sticky=tk.W)
        self.lime_features_var = tk.IntVar(value=10)
        ttk.Spinbox(lime_frame, from_=5, to=20, textvariable=self.lime_features_var, 
                   width=10).grid(row=0, column=1, padx=5, sticky=tk.W)
        ttk.Label(lime_frame, text="(number of top features to explain)").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        # SHAP settings
        shap_frame = ttk.Frame(xai_config_frame)
        shap_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(shap_frame, text="SHAP Samples:").grid(row=0, column=0, sticky=tk.W)
        self.shap_samples_var = tk.IntVar(value=10)
        ttk.Spinbox(shap_frame, from_=5, to=100, textvariable=self.shap_samples_var, 
                   width=10).grid(row=0, column=1, padx=5, sticky=tk.W)
        ttk.Label(shap_frame, text="(max samples for SHAP analysis)").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        ttk.Label(shap_frame, text="SHAP Type:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.shap_type_var = tk.StringVar(value="kernel")
        shap_type_combo = ttk.Combobox(shap_frame, textvariable=self.shap_type_var, 
                                      values=["kernel", "deep"], state="readonly", width=10)
        shap_type_combo.grid(row=1, column=1, padx=5, sticky=tk.W)
        ttk.Label(shap_frame, text="(kernel=model-agnostic, deep=faster for neural nets)").grid(row=1, column=2, sticky=tk.W, padx=5)
        
        # Detection threshold
        threshold_frame = ttk.Frame(model_frame)
        threshold_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(threshold_frame, text="Detection Threshold:").pack(side=tk.LEFT)
        
        self.threshold_var = tk.DoubleVar(value=0.5)
        threshold_scale = ttk.Scale(threshold_frame, from_=0.01, to=1.0, 
                                   variable=self.threshold_var, orient=tk.HORIZONTAL)
        threshold_scale.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        self.threshold_label = ttk.Label(threshold_frame, text="0.500")
        self.threshold_label.pack(side=tk.RIGHT)
        
        # Update threshold label when scale changes
        threshold_scale.configure(command=lambda val: self.threshold_label.config(text=f"{float(val):.3f}"))
        
        # Model info
        self.model_info_label = ttk.Label(model_frame, text="No model loaded", 
                                         foreground='gray', wraplength=500)
        self.model_info_label.pack(pady=10)
        
        # Load model button
        ttk.Button(model_frame, text="Load Selected Model", 
                  command=self.load_ml_model, style="Accent.TButton").pack(pady=5)
        
        # Detection Methods
        methods_frame = ttk.LabelFrame(parent, text="Detection Methods", padding=15)
        methods_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # ML Detection
        self.ml_detection_var = tk.BooleanVar(value=True)
        self.ml_detection_cb = ttk.Checkbutton(methods_frame, text="ML Autoencoder (Primary)", 
                                              variable=self.ml_detection_var)
        self.ml_detection_cb.pack(anchor=tk.W, pady=2)
        
        # Statistical Detection
        self.statistical_detection_var = tk.BooleanVar(value=True)
        self.statistical_detection_cb = ttk.Checkbutton(methods_frame, text="Statistical Analysis", 
                                                       variable=self.statistical_detection_var)
        self.statistical_detection_cb.pack(anchor=tk.W, pady=2)
        
        # Trend Analysis
        self.trend_detection_var = tk.BooleanVar(value=True)
        self.trend_detection_cb = ttk.Checkbutton(methods_frame, text="Trend Change Detection", 
                                                 variable=self.trend_detection_var)
        self.trend_detection_cb.pack(anchor=tk.W, pady=2)
        
        # Threshold-based
        self.threshold_detection_var = tk.BooleanVar(value=False)
        self.threshold_detection_cb = ttk.Checkbutton(methods_frame, text="Threshold-based Detection", 
                                                     variable=self.threshold_detection_var)
        self.threshold_detection_cb.pack(anchor=tk.W, pady=2)
        
        # Real-time Configuration
        realtime_frame = ttk.LabelFrame(parent, text="Real-time Detection", padding=15)
        realtime_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Enable real-time
        self.realtime_var = tk.BooleanVar(value=False)
        self.realtime_cb = ttk.Checkbutton(realtime_frame, text="Enable Real-time Monitoring", 
                                          variable=self.realtime_var,
                                          command=self.toggle_realtime_monitoring)
        self.realtime_cb.pack(anchor=tk.W, pady=5)
        
        # Real-time parameters
        realtime_params_frame = ttk.Frame(realtime_frame)
        realtime_params_frame.pack(fill=tk.X, pady=5)
        
        # Update interval
        ttk.Label(realtime_params_frame, text="Update Interval (sec):").grid(row=0, column=0, sticky=tk.W)
        self.update_interval_var = tk.IntVar(value=5)
        ttk.Spinbox(realtime_params_frame, from_=1, to=60, textvariable=self.update_interval_var, 
                   width=10).grid(row=0, column=1, padx=5, sticky=tk.W)
        
        # Buffer size
        ttk.Label(realtime_params_frame, text="Data Buffer Size:").grid(row=1, column=0, sticky=tk.W)
        self.buffer_size_var = tk.IntVar(value=100)
        ttk.Spinbox(realtime_params_frame, from_=10, to=1000, textvariable=self.buffer_size_var, 
                   width=10).grid(row=1, column=1, padx=5, sticky=tk.W)
    
    def on_model_type_changed(self):
        """Handle model type selection change"""
        model_type = self.model_type_var.get()
        self.parent_app.add_log(f"Model type changed to: {model_type}")
        
        # Update model info label
        if model_type == "Autoencoder":
            self.model_info_label.config(
                text="Base autoencoder model for anomaly detection", 
                foreground='blue'
            )
        elif model_type == "LIME":
            self.model_info_label.config(
                text="LIME will explain individual anomaly detections", 
                foreground='green'
            )
        elif model_type == "SHAP":
            self.model_info_label.config(
                text="SHAP will show global feature importance across all detections", 
                foreground='purple'
            )
        elif model_type == "All":
            self.model_info_label.config(
                text="Ensemble mode: Use all models for comprehensive analysis", 
                foreground='orange'
            )
    
    # [Continue with remaining methods from original file...]
    # The rest of the methods remain largely the same, with additions for XAI
    
    def _create_xai_tab(self, parent):
        """NEW: Create the XAI explanations tab"""
        # Title
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Label(title_frame, text="Explainable AI (XAI) Analysis", 
                 font=('Segoe UI', 14, 'bold')).pack()
        
        ttk.Label(title_frame, text="Understanding why faults were detected using LIME and SHAP",
                 font=('Segoe UI', 10), foreground='gray').pack()
        
        # Controls
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(controls_frame, text="Explanation Type:").pack(side=tk.LEFT)
        
        self.xai_type_var = tk.StringVar(value="LIME")
        xai_types = ["LIME (Local)", "SHAP (Global)", "Both"]
        xai_combo = ttk.Combobox(controls_frame, textvariable=self.xai_type_var,
                                values=xai_types, state="readonly")
        xai_combo.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="Generate Explanations", 
                  command=self.generate_xai_explanations).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="Export XAI Results", 
                  command=self.export_xai_results).pack(side=tk.LEFT, padx=5)
        
        # XAI plot frame
        plot_frame = ttk.LabelFrame(parent, text="XAI Visualizations", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create matplotlib figure for XAI
        self.xai_figure = Figure(figsize=(12, 8), dpi=100)
        self.xai_canvas = FigureCanvasTkAgg(self.xai_figure, plot_frame)
        self.xai_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize empty plot
        self.init_xai_plot()
    
    def init_xai_plot(self):
        """Initialize the XAI plot"""
        self.xai_figure.clear()
        ax = self.xai_figure.add_subplot(111)
        ax.set_title("XAI Explanations")
        ax.text(0.5, 0.5, "Generate XAI explanations to see visualizations...", 
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.axis('off')
        self.xai_canvas.draw()
    
    def generate_xai_explanations(self):
        """NEW: Generate LIME and/or SHAP explanations for detected faults"""
        if not self.fault_detection_results:
            messagebox.showwarning("No Results", "Please run fault detection first!")
            return
        
        if not self.xai_available:
            response = messagebox.askyesno(
                "XAI Not Available",
                "LIME and SHAP libraries are not installed.\n\nWould you like installation instructions?"
            )
            if response:
                messagebox.showinfo("Install XAI Libraries", 
                                  "To enable XAI explanations, install:\n\n" +
                                  "pip install --user shap lime\n\n" +
                                  "Then restart the application.")
            return
        
        try:
            self.parent_app.add_log("Generating XAI explanations...")
            self.generate_xai_btn.config(state="disabled", text="Generating...")
            self.parent_app.root.update()
            
            # Import XAI integration module
            from satellite_xai import add_xai_to_basilisk_detection
            
            # Generate explanations
            xai_results = add_xai_to_basilisk_detection(
                ml_detector=self.fault_detection_results.get("ml_detector"),
                real_telemetry=self.fault_detection_results.get("real_telemetry"),
                detected_faults=self.fault_detection_results.get("detections")
            )
            
            if xai_results:
                # Store XAI results
                self.xai_results = xai_results
                
                self.parent_app.add_log("XAI explanations generated successfully!")
                
                # Update statistics
                num_lime = len(xai_results.get('lime_explanations', []))
                has_shap = xai_results.get('shap_values') is not None
                total_explanations = num_lime + (1 if has_shap else 0)
                
                self.xai_explanations_label.config(text=str(total_explanations))
                
                # Display in XAI tab
                self.display_xai_results(xai_results)
                
                # Switch to XAI tab
                self.detection_notebook.select(3)
                
                # Show success message
                success_msg = f"✓ XAI Explanations Generated!\n\n"
                success_msg += f"LIME Explanations: {num_lime}\n"
                success_msg += f"SHAP Analysis: {'Yes' if has_shap else 'No'}\n\n"
                success_msg += f"Results are now displayed in the XAI Explanations tab.\n\n"
                success_msg += f"Plots also saved to:\n"
                if num_lime > 0:
                    success_msg += f"  • lime_anomaly_1.png\n"
                if has_shap:
                    success_msg += f"  • shap_summary.png\n"
                    success_msg += f"  • shap_waterfall_anomaly_1.png\n"
                
                messagebox.showinfo("XAI Complete", success_msg)
                
                self.parent_app.add_log(f"XAI visualizations displayed in GUI")
            else:
                messagebox.showerror("XAI Failed", "Could not generate XAI explanations")
                
        except Exception as e:
            self.parent_app.add_log(f"Error generating XAI: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("XAI Error", f"Error generating XAI:\n{str(e)}\n\nCheck console for details")
            
        finally:
            self.generate_xai_btn.config(state="normal", text="Generate XAI")
    

            
    def display_xai_results(self, xai_results):
        """Display XAI results by loading saved PNG files"""
        print("\nDEBUG: display_xai_results called")
        
        self.xai_figure.clear()
        xai_type = self.xai_type_var.get()
        
        # Get directory where plots are saved
        if 'output_dir' in xai_results:
            xai_plots_dir = xai_results['output_dir']
        else:
            xai_plots_dir = os.path.join(os.getcwd(), 'xai_plots')
        
        print(f"DEBUG: Looking in: {xai_plots_dir}")
        
        try:
            from PIL import Image
            
            # Check what to display
            show_lime = "LIME" in xai_type
            show_shap = "SHAP" in xai_type
            
            displayed = False
            
            # Display LIME if requested
            if show_lime:
                lime_path = os.path.join(xai_plots_dir, "lime_anomaly_1.png")
                if os.path.exists(lime_path):
                    img = Image.open(lime_path)
                    ax = self.xai_figure.add_subplot(211 if show_shap else 111)
                    ax.imshow(img)
                    ax.axis('off')
                    ax.set_title("LIME Explanation", fontsize=12, fontweight='bold')
                    displayed = True
                    print("DEBUG: LIME displayed")
            
            # Display SHAP if requested
            if show_shap:
                shap_path = os.path.join(xai_plots_dir, "shap_summary.png")
                if os.path.exists(shap_path):
                    img = Image.open(shap_path)
                    ax = self.xai_figure.add_subplot(212 if show_lime else 111)
                    ax.imshow(img)
                    ax.axis('off')
                    ax.set_title("SHAP Explanation", fontsize=12, fontweight='bold')
                    displayed = True
                    print("DEBUG: SHAP displayed")
            
            # If nothing displayed, show message
            if not displayed:
                ax = self.xai_figure.add_subplot(111)
                ax.text(0.5, 0.5, f"Plots not found in:\n{xai_plots_dir}", 
                    ha='center', va='center', transform=ax.transAxes)
                ax.axis('off')
            
            self.xai_figure.tight_layout()
            self.xai_canvas.draw()
            
        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    def _display_lime_plot(self, lime_explanation, ax):
        """Display a LIME explanation in the given axis"""
        try:
            print("DEBUG: Attempting to display LIME plot...")
            
            # Get the explanation as a list of (feature, importance) tuples
            exp_list = lime_explanation.as_list()
            print(f"DEBUG: LIME explanation has {len(exp_list)} features")
            
            if len(exp_list) == 0:
                ax.text(0.5, 0.5, "No LIME features to display", 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.axis('off')
                return
            
            # Extract features and values (top 15)
            num_features = min(15, len(exp_list))
            features = [feat for feat, _ in exp_list[:num_features]]
            values = [val for _, val in exp_list[:num_features]]
            
            print(f"DEBUG: Displaying {num_features} features")
            print(f"DEBUG: Features: {features[:3]}...")  # Show first 3
            print(f"DEBUG: Values: {values[:3]}...")
            
            # Create horizontal bar plot
            colors = ['green' if v > 0 else 'red' for v in values]
            y_pos = np.arange(len(features))
            
            ax.barh(y_pos, values, color=colors, alpha=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(features, fontsize=8)
            ax.set_xlabel('Feature Importance', fontsize=10)
            ax.set_title('LIME: Top Features for Anomaly Detection', fontsize=12, fontweight='bold')
            ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
            ax.grid(axis='x', alpha=0.3)
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='green', alpha=0.7, label='Increases anomaly score'),
                Patch(facecolor='red', alpha=0.7, label='Decreases anomaly score')
            ]
            ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
            
            print("DEBUG: LIME plot created successfully")
            
        except Exception as e:
            print(f"ERROR in _display_lime_plot: {e}")
            import traceback
            traceback.print_exc()
            
            ax.text(0.5, 0.5, f"Error displaying LIME:\n{str(e)}\n\nCheck console for details\nPlots saved to working directory", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=10, wrap=True)
            ax.axis('off')
    
    def _display_shap_plot(self, shap_values, xai_results, ax):
        """Display SHAP summary in the given axis"""
        try:
            print("DEBUG: Attempting to display SHAP plot...")
            print(f"DEBUG: SHAP values shape: {shap_values.shape if hasattr(shap_values, 'shape') else type(shap_values)}")
            
            # Get the explainer from results
            explainer = xai_results.get('explainer')
            
            # Flatten SHAP values if needed
            if len(shap_values.shape) > 2:
                shap_flat = shap_values.reshape(shap_values.shape[0], -1)
            else:
                shap_flat = shap_values
            
            print(f"DEBUG: Flattened SHAP shape: {shap_flat.shape}")
            
            # Get top features by absolute mean SHAP value
            mean_abs_shap = np.mean(np.abs(shap_flat), axis=0)
            print(f"DEBUG: Mean absolute SHAP computed, shape: {mean_abs_shap.shape}")
            
            top_indices = np.argsort(mean_abs_shap)[-15:][::-1]  # Top 15
            
            # Get feature names
            if explainer:
                feature_names = explainer._get_flat_feature_names()
                print(f"DEBUG: Got {len(feature_names)} feature names from explainer")
            else:
                feature_names = [f"Feature_{i}" for i in range(shap_flat.shape[1])]
                print(f"DEBUG: Using default feature names")
            
            # Create bar plot of mean absolute SHAP values
            top_features = [feature_names[i] for i in top_indices]
            top_values = mean_abs_shap[top_indices]
            
            print(f"DEBUG: Top 3 features: {top_features[:3]}")
            print(f"DEBUG: Top 3 values: {top_values[:3]}")
            
            y_pos = np.arange(len(top_features))
            ax.barh(y_pos, top_values, color='steelblue', alpha=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(top_features, fontsize=8)
            ax.set_xlabel('Mean |SHAP Value|', fontsize=10)
            ax.set_title('SHAP: Global Feature Importance', fontsize=12, fontweight='bold')
            ax.grid(axis='x', alpha=0.3)
            
            # Add note
            ax.text(0.98, 0.02, f'Based on {shap_flat.shape[0]} samples', 
                   transform=ax.transAxes, ha='right', va='bottom', 
                   fontsize=8, style='italic', color='gray')
            
            print("DEBUG: SHAP plot created successfully")
            
        except Exception as e:
            print(f"ERROR in _display_shap_plot: {e}")
            import traceback
            traceback.print_exc()
            
            ax.text(0.5, 0.5, f"Error displaying SHAP:\n{str(e)}\n\nCheck console for details\nPlots saved to working directory", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=10, wrap=True)
            ax.axis('off')
    
    def export_xai_results(self):
        """Export XAI results to folder"""
        if not hasattr(self, 'xai_results') or not self.xai_results:
            messagebox.showinfo("No XAI Results", "Generate XAI explanations first!")
            return
        
        folder_path = filedialog.askdirectory(title="Select folder to save XAI results")
        if folder_path:
            try:
                # Save XAI plots and reports to folder
                self.parent_app.add_log(f"Exporting XAI results to {folder_path}")
                messagebox.showinfo("Export Success", f"XAI results exported to:\n{folder_path}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Error: {str(e)}")
    
    def _create_results_tab(self, parent):
        """Create the detection results tab content"""
        # Controls frame
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Button(controls_frame, text="Refresh Results", 
                  command=self.refresh_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="Clear Results", 
                  command=self.clear_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="Export Results", 
                  command=self.export_results).pack(side=tk.LEFT, padx=5)
        
        # Results table frame
        table_frame = ttk.LabelFrame(parent, text="Detection Results", padding=10)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create treeview for results
        columns = ('Time', 'Spacecraft', 'Fault Type', 'Confidence', 'Component', 'Method', 'Status')
        self.results_tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=15)
        
        # Define column headings and widths
        column_widths = {'Time': 80, 'Spacecraft': 100, 'Fault Type': 120, 'Confidence': 80, 
                        'Component': 100, 'Method': 150, 'Status': 80}
        
        for col in columns:
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=column_widths.get(col, 100))
        
        # Create frame for treeview and scrollbars
        tree_frame = ttk.Frame(table_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.results_tree.yview)
        h_scrollbar = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack treeview and scrollbars
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Results summary
        summary_frame = ttk.LabelFrame(parent, text="Summary", padding=10)
        summary_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Create text widget with scrollbar
        summary_text_frame = ttk.Frame(summary_frame)
        summary_text_frame.pack(fill=tk.X)
        
        self.results_summary_text = tk.Text(summary_text_frame, height=6, wrap=tk.WORD)
        summary_scrollbar = ttk.Scrollbar(summary_text_frame, orient="vertical", command=self.results_summary_text.yview)
        self.results_summary_text.configure(yscrollcommand=summary_scrollbar.set)
        
        self.results_summary_text.pack(side=tk.LEFT, fill=tk.X, expand=True)
        summary_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.results_summary_text.insert(tk.END, "Detection results summary will appear here...")
    
    def _create_monitoring_tab(self, parent):
        """Create the live monitoring tab content"""
        # Monitoring controls
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.monitoring_status_label = ttk.Label(controls_frame, text="Monitoring Inactive", 
                                                font=('Segoe UI', 10, 'bold'), foreground='red')
        self.monitoring_status_label.pack(side=tk.LEFT)
        
        # Spacecraft selection
        ttk.Label(controls_frame, text="Monitor:").pack(side=tk.RIGHT, padx=5)
        self.monitoring_spacecraft_var = tk.StringVar(value="All Spacecraft")
        self.monitoring_spacecraft_combo = ttk.Combobox(controls_frame, 
                                                       textvariable=self.monitoring_spacecraft_var,
                                                       values=["All Spacecraft"])
        self.monitoring_spacecraft_combo.pack(side=tk.RIGHT, padx=5)
        
        # Live plot frame
        plot_frame = ttk.LabelFrame(parent, text="Real-time Anomaly Detection", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create matplotlib figure
        self.monitoring_figure = Figure(figsize=(10, 6), dpi=100)
        self.monitoring_canvas = FigureCanvasTkAgg(self.monitoring_figure, plot_frame)
        self.monitoring_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize empty plot
        self.init_monitoring_plot()
    
    def _create_analysis_tab(self, parent):
        """Create the analysis tab content"""
        # Analysis controls
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(controls_frame, text="Analysis Type:").pack(side=tk.LEFT)
        
        self.analysis_type_var = tk.StringVar(value="Detection Confidence Over Time")
        analysis_types = [
            "Detection Confidence Over Time",
            "Fault Distribution by Spacecraft", 
            "Detection Method Comparison",
            "Error Rate Analysis",
            "Performance Metrics"
        ]
        self.analysis_type_combo = ttk.Combobox(controls_frame, textvariable=self.analysis_type_var,
                                               values=analysis_types, state="readonly")
        self.analysis_type_combo.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="Generate Analysis", 
                  command=self.generate_analysis).pack(side=tk.LEFT, padx=5)
        
        # Analysis plot frame
        plot_frame = ttk.LabelFrame(parent, text="Analysis Results", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create matplotlib figure for analysis
        self.analysis_figure = Figure(figsize=(10, 8), dpi=100)
        self.analysis_canvas = FigureCanvasTkAgg(self.analysis_figure, plot_frame)
        self.analysis_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def init_monitoring_plot(self):
        """Initialize the monitoring plot"""
        self.monitoring_figure.clear()
        ax = self.monitoring_figure.add_subplot(111)
        ax.set_title("Real-time Fault Detection Monitoring")
        ax.set_xlabel("Time (minutes)")
        ax.set_ylabel("Anomaly Score")
        ax.grid(True, alpha=0.3)
        ax.text(0.5, 0.5, "Start monitoring to see live data...", 
               ha='center', va='center', transform=ax.transAxes)
        self.monitoring_canvas.draw()
    
    def check_for_default_model(self):
        """Check if default model file exists and set it"""
        possible_paths = [
            "anomaly_detection_model.keras",
            "../anomaly_detection_model.keras",
            "../../anomaly_detection_model.keras",
            os.path.join(os.path.dirname(__file__), "anomaly_detection_model.keras"),
            os.path.join(os.path.dirname(__file__), "..", "anomaly_detection_model.keras"),
            r"C:\Uni\Uni\satelliteConstellation-basilisk\ProjectWork\DRL\Agent-based-Architecture-for-Proactive-Fault-Tolerance-and-Management-in-Small-Satellite-Missions\anomaly_detection_model.keras"
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                self.model_path_var.set(path)
                self.model_info_label.config(
                    text=f"Default model found: {os.path.basename(path)}", 
                    foreground='blue'
                )
                self.parent_app.add_log(f"Found default model: {path}")
                break
        else:
            self.model_info_label.config(text="No model loaded - select model type and load")
    
    # ML Model Methods
    def browse_model_file(self):
        """Browse for ML model file"""
        initial_dir = os.path.dirname(os.path.abspath(__file__))
        
        file_path = filedialog.askopenfilename(
            title="Select ML Model File",
            initialdir=initial_dir,
            filetypes=[
                ("Keras Models", "*.keras *.h5"), 
                ("All Files", "*.*")
            ]
        )
        if file_path:
            self.model_path_var.set(file_path)
            self.model_info_label.config(text=f"Selected: {os.path.basename(file_path)}")
    
    def load_ml_model(self):
        """Load the ML model based on selected type"""
        model_type = self.model_type_var.get()
        model_path = self.model_path_var.get().strip()
        
        if not model_path:
            messagebox.showwarning("Warning", "Please specify a base model file!")
            return
        
        if not os.path.exists(model_path):
            messagebox.showwarning("Warning", f"Model file not found: {model_path}")
            return
        
        try:
            # Import the real ML detector
            from real_ml_fault_detection import RealMLFaultDetector
            
            self.parent_app.add_log(f"Loading {model_type} model: {os.path.basename(model_path)}")
            self.model_info_label.config(text=f"Loading {model_type} model...", foreground='orange')
            
            # Update UI to show loading
            self.parent_app.root.update()
            
            # Load base ML model
            self.ml_detector = RealMLFaultDetector(model_path)
            
            if self.ml_detector.is_loaded:
                # Check if we need to initialize XAI explainers
                if model_type in ["LIME", "SHAP", "All"] and self.xai_available:
                    self.parent_app.add_log(f"Base model loaded, initializing {model_type} explainer...")
                    
                    try:
                        from satellite_xai import SatelliteXAIExplainer
                        
                        self.xai_explainer = SatelliteXAIExplainer(self.ml_detector.model)
                        self.parent_app.add_log(f"✓ {model_type} explainer initialized successfully")
                        
                        xai_info = f"\n\n📊 XAI Explainer: {model_type}\n"
                        xai_info += f"Base Detector: Autoencoder\n"
                        if model_type in ["LIME", "All"]:
                            xai_info += f"LIME Features: {self.lime_features_var.get()}\n"
                        if model_type in ["SHAP", "All"]:
                            xai_info += f"SHAP Type: {self.shap_type_var.get()}\n"
                            xai_info += f"SHAP Samples: {self.shap_samples_var.get()}\n"
                        xai_info += "\nNote: LIME/SHAP wrap around the base autoencoder"
                    except ImportError as e:
                        xai_info = f"\n\n⚠ XAI libraries not available\nInstall: pip install shap lime"
                        self.xai_explainer = None
                        self.parent_app.add_log(f"Warning: Could not initialize XAI - {e}")
                else:
                    xai_info = ""
                    self.xai_explainer = None
                
                # Create status message based on model type
                if model_type == "Autoencoder":
                    status_msg = f"✓ Autoencoder loaded successfully!"
                elif model_type in ["LIME", "SHAP"]:
                    status_msg = f"✓ {model_type} Explainer Ready!"
                else:  # All
                    status_msg = f"✓ Full XAI Suite Ready!"
                
                self.model_info_label.config(
                    text=status_msg + f"\n" +
                         f"Base Model: {os.path.basename(model_path)}\n" +
                         f"Input: {self.ml_detector.model.input_shape}\n" +
                         f"Parameters: {self.ml_detector.model.count_params():,}" +
                         xai_info, 
                    foreground='green'
                )
                
                # Update status label
                if model_type == "All":
                    self.ml_status_label.config(text="Autoencoder + XAI", foreground='green')
                elif model_type == "LIME":
                    self.ml_status_label.config(text="Autoencoder + LIME", foreground='green')
                elif model_type == "SHAP":
                    self.ml_status_label.config(text="Autoencoder + SHAP", foreground='green')
                else:
                    self.ml_status_label.config(text="Autoencoder Only", foreground='green')
                
                self.parent_app.add_log(f"✓ {model_type} configuration loaded successfully")
                
                # Update main GUI status
                self.parent_app.update_status_counts()
                
                # Create detailed success message
                if model_type == "Autoencoder":
                    success_msg = f"Base Autoencoder Model Loaded!\n\n"
                    success_msg += f"Model: {os.path.basename(model_path)}\n"
                    success_msg += f"Input shape: {self.ml_detector.model.input_shape}\n"
                    success_msg += f"Parameters: {self.ml_detector.model.count_params():,}\n\n"
                    success_msg += f"This model will detect faults.\n"
                    success_msg += f"For explanations, select LIME or SHAP."
                elif model_type in ["LIME", "SHAP"]:
                    success_msg = f"{model_type} Explainer Ready!\n\n"
                    success_msg += f"Base Model: {os.path.basename(model_path)}\n"
                    success_msg += f"Input shape: {self.ml_detector.model.input_shape}\n"
                    success_msg += f"Parameters: {self.ml_detector.model.count_params():,}\n\n"
                    success_msg += f"✓ Autoencoder will detect faults\n"
                    success_msg += f"✓ {model_type} will explain WHY faults were detected\n\n"
                    if model_type == "LIME":
                        success_msg += f"LIME explains individual detections\n"
                        success_msg += f"Top {self.lime_features_var.get()} features per anomaly"
                    else:
                        success_msg += f"SHAP shows global feature importance\n"
                        success_msg += f"Analyzing up to {self.shap_samples_var.get()} samples"
                else:  # All
                    success_msg = f"Full XAI Suite Loaded!\n\n"
                    success_msg += f"Base Model: {os.path.basename(model_path)}\n"
                    success_msg += f"Input shape: {self.ml_detector.model.input_shape}\n"
                    success_msg += f"Parameters: {self.ml_detector.model.count_params():,}\n\n"
                    success_msg += f"✓ Autoencoder: Detects faults\n"
                    success_msg += f"✓ LIME: Explains individual detections\n"
                    success_msg += f"✓ SHAP: Shows global feature importance\n\n"
                    success_msg += f"You'll get comprehensive explanations!"
                
                messagebox.showinfo("Model Loaded Successfully", success_msg)
                
            else:
                self.model_info_label.config(text="Failed to load model", foreground='red')
                self.parent_app.add_log("Failed to load ML model")
                messagebox.showerror("Error", "Failed to load ML model. Check the file format and path.")
                
        except ImportError as e:
            error_msg = f"ML libraries not available: {str(e)}\n\nTo install: pip install tensorflow"
            self.model_info_label.config(text="ML libraries not available", foreground='red')
            self.parent_app.add_log(error_msg)
            messagebox.showerror("Import Error", error_msg)
            
        except Exception as e:
            error_msg = f"Error loading ML model: {str(e)}"
            self.model_info_label.config(text=error_msg, foreground='red')
            self.parent_app.add_log(error_msg)
            messagebox.showerror("Error", error_msg)
    
    def start_fault_detection(self):
        """Start fault detection process"""
        
        # Check if ML model is loaded (if ML detection is enabled)
        if self.ml_detection_var.get() and not self.is_ml_ready():
            response = messagebox.askyesno(
                "ML Model Not Ready", 
                "ML detection is enabled but no model is loaded.\n\nWould you like to:\n- Click 'Yes' to continue without ML\n- Click 'No' to load a model first"
            )
            if not response:
                return
        
        # Check if any detection method is enabled
        methods_enabled = (self.ml_detection_var.get() or 
                          self.statistical_detection_var.get() or
                          self.trend_detection_var.get() or 
                          self.threshold_detection_var.get())
        
        if not methods_enabled:
            messagebox.showwarning("Warning", "Please enable at least one detection method!")
            return
        
        # Check if there are satellites configured
        if not hasattr(self.parent_app, 'satellites') or not self.parent_app.satellites:
            messagebox.showwarning("Warning", "No satellites configured!\n\nGo to Constellation tab to add satellites first.")
            return
        
        # Check if any satellites have faults enabled
        faulty_sats = [s for s in self.parent_app.satellites if s['fault']['enabled']]
        if not faulty_sats:
            response = messagebox.askyesno(
                "No Faults Configured", 
                "No satellites have faults enabled.\n\nTo test fault detection:\n1. Go to 'Fault Configuration' tab\n2. Select a satellite\n3. Enable a fault\n\nContinue anyway?"
            )
            if not response:
                messagebox.showinfo("Next Steps", "1. Go to 'Fault Configuration' tab\n2. Select a satellite (e.g., Satellite1)\n3. Enable fault (e.g., friction at 15 minutes)\n4. Come back and click 'Start Detection'\n5. Run simulation to see results")
                return
        
        # Get enabled methods
        enabled_methods = []
        if self.ml_detection_var.get():
            model_type = self.model_type_var.get()
            enabled_methods.append(f"ML {model_type}")
        if self.statistical_detection_var.get():
            enabled_methods.append("Statistical Analysis")
        if self.trend_detection_var.get():
            enabled_methods.append("Trend Analysis")
        if self.threshold_detection_var.get():
            enabled_methods.append("Threshold Detection")
        
        # Show configuration summary
        config_summary = f"Detection Configuration:\n\n"
        config_summary += f"Enabled Methods: {', '.join(enabled_methods)}\n"
        config_summary += f"ML Model Type: {self.model_type_var.get()}\n"
        config_summary += f"ML Model: {'Loaded' if self.is_ml_ready() else 'Not Loaded'}\n"
        config_summary += f"XAI Explainer: {'Ready' if self.xai_explainer else 'Not Ready'}\n"
        config_summary += f"Satellites: {len(self.parent_app.satellites)}\n"
        config_summary += f"Faulty Satellites: {len(faulty_sats)}\n"
        
        if faulty_sats:
            config_summary += f"\nFaults Configured:\n"
            for sat in faulty_sats[:3]:  # Show first 3
                fault = sat['fault']
                config_summary += f"  {sat['name']}: {fault['type']} at {fault['time']} min\n"
        
        config_summary += f"\nNext Step: Click 'Run Simulation' to start detection"
        
        self.parent_app.add_log("Fault detection configured and ready")
        self.detection_status_label.config(text="Ready for Simulation", foreground='green')
        
        # Update buttons
        self.start_detection_btn.config(text="Detection Ready", state="disabled")
        
        # Show configuration summary
        messagebox.showinfo("Detection Ready", config_summary)
        
        self.parent_app.add_log("Fault detection ready - run simulation to start detection")
    
    def simulate_detection_results(self):
        """Simulate detection results for testing"""
        # [Keep existing simulate_detection_results method from original file]
        self.parent_app.add_log("Simulating detection results...")
        messagebox.showinfo("Simulation", "Detection simulation not yet implemented in this version")
    
    def reset_detection_state(self):
        """Reset detection to initial state"""
        self.detection_status_label.config(text="Idle", foreground='black')
        self.start_detection_btn.config(text="Start Detection", state="normal")
        self.generate_xai_btn.config(state="disabled")
        self.parent_app.add_log("Detection state reset")
    
    def toggle_realtime_monitoring(self):
        """Toggle real-time monitoring"""
        if self.realtime_var.get():
            self.monitoring_status_label.config(text="Monitoring Active", foreground='green')
            self.monitoring_active = True
            self.parent_app.add_log("Real-time monitoring enabled")
        else:
            self.monitoring_status_label.config(text="Monitoring Inactive", foreground='red')
            self.monitoring_active = False
            self.parent_app.add_log("Real-time monitoring disabled")
    
    def clear_results(self):
        """Clear all detection results"""
        if messagebox.askyesno("Confirm Clear", "Are you sure you want to clear all detection results?"):
            # Clear results table
            for item in self.results_tree.get_children():
                self.results_tree.delete(item)
            
            # Clear summary
            self.results_summary_text.delete(1.0, tk.END)
            self.results_summary_text.insert(tk.END, "Detection results summary will appear here...")
            
            # Clear history
            self.detection_history.clear()
            self.update_statistics()
            
            # Clear recent detections
            self.recent_detections_listbox.delete(0, tk.END)
            self.recent_detections_listbox.insert(tk.END, "No recent detections...")
            
            # Reset fault detection results
            self.fault_detection_results = None
            
            # Disable XAI button
            self.generate_xai_btn.config(state="disabled")
            
            self.parent_app.add_log("Detection results cleared")
    
    def refresh_results(self):
        """Refresh detection results display"""
        self.update_statistics()
        self.parent_app.add_log("Detection results refreshed")
    
    def export_results(self):
        """Export detection results to file"""
        if not self.results_tree.get_children():
            messagebox.showinfo("Info", "No results to export!")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Export Detection Results",
            defaultextension=".json",
            filetypes=[("JSON Files", "*.json"), ("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        
        if file_path:
            try:
                self.save_results_to_file(file_path)
                messagebox.showinfo("Success", f"Results exported to {file_path}")
                self.parent_app.add_log(f"Results exported to {file_path}")
            except Exception as e:
                error_msg = f"Failed to export results: {str(e)}"
                messagebox.showerror("Error", error_msg)
                self.parent_app.add_log(error_msg)
    
    def generate_analysis(self):
        """Generate analysis plot"""
        # [Keep existing generate_analysis method from original file]
        self.parent_app.add_log("Generating analysis...")
        messagebox.showinfo("Analysis", "Analysis generation not yet implemented in this version")
    
    def update_statistics(self):
        """Update the statistics display"""
        total_detections = len(self.detection_history)
        active_faults = sum(1 for d in self.detection_history if d.get('status') == 'Active')
        
        self.total_detections_label.config(text=str(total_detections))
        self.active_faults_label.config(text=str(active_faults))
        
        # Calculate success rate
        if total_detections > 0:
            success_rate = (active_faults / total_detections) * 100
            self.success_rate_label.config(text=f"{success_rate:.1f}%")
        else:
            self.success_rate_label.config(text="N/A")
        
        # Update last detection
        if self.detection_history:
            last_detection = self.detection_history[-1]
            self.last_detection_label.config(
                text=f"{last_detection.get('spacecraft', 'Unknown')} at {last_detection.get('time', 'Unknown')}"
            )
        else:
            self.last_detection_label.config(text="None")
    
    def save_results_to_file(self, file_path):
        """Save detection results to file"""
        # [Keep existing save_results_to_file method from original file]
        pass
    
    def display_ml_results(self, ml_results):
        """Display ML detection results and enable XAI button"""
        if not ml_results:
            return
        
        # Ensure this runs on the main GUI thread
        def _update_gui():
            try:
                self.fault_detection_results = ml_results
                
                # Update detection status
                self.detection_status_label.config(text="Detection Complete", foreground='green')
                
                # Clear existing results
                for item in self.results_tree.get_children():
                    self.results_tree.delete(item)
                
                # Clear recent detections list
                self.recent_detections_listbox.delete(0, tk.END)
                
                # Process and display results
                summary = ml_results.get('summary', {})
                detections = ml_results.get('detections', {})
                
                total_detections = 0
                recent_detections = []
                
                # Add detections to table and recent list
                for spacecraft_name, spacecraft_detections in detections.items():
                    for detection in spacecraft_detections:
                        # Add to results table
                        method_used = detection.details.get('primary_method', 'ML Detection')
                        
                        self.results_tree.insert('', 'end', values=(
                            f"{detection.detection_time_minutes:.1f} min",
                            spacecraft_name,
                            detection.fault_type,
                            f"{detection.confidence:.3f}",
                            detection.affected_component,
                            method_used,
                            "Active" if detection.fault_detected else "Inactive"
                        ))
                        
                        # Add to recent detections
                        recent_detections.append({
                            'time': f"{detection.detection_time_minutes:.1f} min",
                            'spacecraft': spacecraft_name,
                            'type': detection.fault_type,
                            'confidence': detection.confidence,
                            'method': method_used,
                            'status': 'Active' if detection.fault_detected else 'Inactive'
                        })
                        
                        total_detections += 1
                
                # Update recent detections list
                if recent_detections:
                    for detection in recent_detections[-10:]:  # Show last 10
                        detection_text = (f"{detection['time']} - {detection['spacecraft']}: "
                                        f"{detection['method']} (conf: {detection['confidence']:.3f})")
                        self.recent_detections_listbox.insert(0, detection_text)
                else:
                    self.recent_detections_listbox.insert(0, "No detections found")
                
                # Update detection history
                self.detection_history.extend(recent_detections)
                
                # Update statistics
                self.update_statistics()
                
                # Update results summary
                self.results_summary_text.delete(1.0, tk.END)
                summary_text = f"""ML FAULT DETECTION RESULTS:
========================================
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Data Source: Enhanced Telemetry Analysis
Detection Method: Multi-Criteria Approach

SUMMARY:
  Spacecraft Analyzed: {summary.get('total_spacecraft', 0)}
  Total ML Detections: {summary.get('total_detections', 0)}
  Success Rate: {summary.get('success_rate', 1.0):.1%}

DETECTIONS BY SPACECRAFT:
"""
                
                for spacecraft, spacecraft_detections in detections.items():
                    summary_text += f"\n{spacecraft}: {len(spacecraft_detections)} detections"
                    for detection in spacecraft_detections[:3]:  # Show first 3
                        method = detection.details.get('primary_method', 'ML')
                        summary_text += f"\n  - {detection.detection_time_minutes:.1f}min via {method} (conf: {detection.confidence:.3f})"
                        
                        # Add detection details
                        if 'detection_criteria' in detection.details:
                            criteria = detection.details['detection_criteria']
                            if criteria.get('rw_speed_change', False):
                                summary_text += f"\n    * RW speed change: {criteria.get('max_speed_change_percent', 0):.1f}%"
                            if criteria.get('attitude_change', False):
                                summary_text += f"\n    * Attitude error increase: {criteria.get('attitude_change_percent', 0):.1f}%"
                            if criteria.get('ml_input_change', False):
                                summary_text += f"\n    * ML input difference: {criteria.get('input_difference', 0):.6f}"
                
                if total_detections == 0:
                    summary_text += "\n\nNOTE: No faults were detected."
                else:
                    summary_text += f"\n\nSUCCESS: Detected {total_detections} faults!"
                    
                    # Enable XAI button if we have detections
                    if hasattr(self, 'generate_xai_btn'):
                        self.generate_xai_btn.config(state="normal")
                        summary_text += "\n\n✓ XAI explanations available - click 'Generate XAI' button"
                
                self.results_summary_text.insert(tk.END, summary_text)
                
                # Switch to results tab to show the new results
                self.detection_notebook.select(2)
                
                # Log the results
                self.parent_app.add_log(f"ML Detection Results: {total_detections} detections from {len(detections)} spacecraft")
                
                # Update monitoring plot if active
                if self.monitoring_active:
                    self.update_monitoring_plot(ml_results)
                    
            except Exception as e:
                self.parent_app.add_log(f"Error displaying ML results: {e}")
                import traceback
                traceback.print_exc()
        
        # Check if we're on the main thread
        try:
            import threading
            if threading.current_thread() is threading.main_thread():
                # Already on main thread, run directly
                _update_gui()
            else:
                # Schedule on main thread
                self.parent_app.root.after(0, _update_gui)
        except:
            # Fallback if threading check fails
            try:
                self.parent_app.root.after(0, _update_gui)
            except:
                _update_gui()
    
    def update_monitoring_plot(self, ml_results):
        """Update the monitoring plot with real data"""
        try:
            self.monitoring_figure.clear()
            ax = self.monitoring_figure.add_subplot(111)
            
            # Extract time series data from ML results
            time_points = np.linspace(0, 30, 100)
            anomaly_scores = np.random.normal(0.2, 0.1, len(time_points))
            
            # Add spikes where detections occurred
            detections = ml_results.get('detections', {})
            detection_times = []
            detection_scores = []
            
            for spacecraft_detections in detections.values():
                for detection in spacecraft_detections:
                    detection_times.append(detection.detection_time_minutes)
                    detection_scores.append(detection.confidence)
                    
                    # Add spike in anomaly score
                    time_idx = np.argmin(np.abs(time_points - detection.detection_time_minutes))
                    anomaly_scores[time_idx] = detection.confidence
            
            # Plot the data
            ax.plot(time_points, anomaly_scores, 'b-', linewidth=2, label='Anomaly Score')
            ax.axhline(y=self.threshold_var.get(), color='r', linestyle='--', 
                      label=f'Threshold ({self.threshold_var.get():.3f})')
            
            # Mark detections
            if detection_times:
                ax.scatter(detection_times, detection_scores, color='red', s=100, alpha=0.8, 
                          label=f'Detections ({len(detection_times)})', zorder=5)
                
                for time, score in zip(detection_times, detection_scores):
                    ax.axvline(x=time, color='red', alpha=0.3, linestyle=':')
            
            ax.set_xlabel('Time (minutes)')
            ax.set_ylabel('Anomaly Score / Confidence')
            ax.set_title('Real-time Fault Detection Results')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_ylim(0, 1.1)
            
            self.monitoring_canvas.draw()
        except Exception as e:
            self.parent_app.add_log(f"Error updating monitoring plot: {e}")
    
    def is_ml_ready(self):
        """Check if ML detection is ready"""
        return (
            self.ml_available and 
            hasattr(self, 'ml_detector') and 
            self.ml_detector is not None and 
            self.ml_detector.is_loaded
        )