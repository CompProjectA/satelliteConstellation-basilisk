#!/usr/bin/env python
"""
fault_detection_tab.py (UPDATED WITH ISOLATION FOREST)

Enhanced Fault Detection Tab with integrated LIME, SHAP, and Isolation Forest model selection.
Users can now choose between Autoencoder, Isolation Forest, LIME, and SHAP models from a dropdown.
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
    Enhanced Fault Detection Tab with XAI and Isolation Forest support
    """
    
    def __init__(self, parent_app, parent_frame):
        """Initialize the fault detection tab with XAI and Isolation Forest support"""
        self.parent_app = parent_app
        self.parent_frame = parent_frame
        
        self.fault_detection_results = None
        self.ml_detector = None
        self.isolation_forest_detector = None  # NEW: Isolation Forest detector
        self.xai_explainer = None
        self.gnn_model = None  # NEW: GNN Autoencoder model
        self.detection_history = []
        self.monitoring_active = False
        
        # Check ML and XAI availability
        self.ml_available = self.check_ml_availability()
        self.xai_available = self.check_xai_availability()
        self.sklearn_available = self.check_sklearn_availability()  # NEW: Check for sklearn
        
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
    
    def check_sklearn_availability(self):
        """Check if scikit-learn (for Isolation Forest) is available"""
        try:
            from sklearn.ensemble import IsolationForest
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
        xai_frame = ttk.Frame(self.detection_notebook)
        analysis_frame = ttk.Frame(self.detection_notebook)
        
        # Add sub-tabs
        self.detection_notebook.add(overview_frame, text="Overview")
        self.detection_notebook.add(config_frame, text="ML Configuration")
        self.detection_notebook.add(results_frame, text="Detection Results")
        self.detection_notebook.add(xai_frame, text="XAI Explanations")
        self.detection_notebook.add(analysis_frame, text="GNN Analysis")
        
        # Create tab contents
        self._create_overview_tab(overview_frame)
        self._create_config_tab(config_frame)
        self._create_results_tab(results_frame)
        self._create_xai_tab(xai_frame)
        self._create_analysis_tab(analysis_frame)
        
    def _create_overview_tab(self, parent):
        """Create the overview tab content"""
        # Title frame
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(10, 20))
        
        title_label = ttk.Label(title_frame, text="ML Fault Detection System with XAI & Isolation Forest", 
                               font=('Segoe UI', 16, 'bold'))
        title_label.pack()
        
        subtitle_label = ttk.Label(title_frame, text="Real-time fault detection with LSTM, Isolation Forest, LIME and SHAP explanations",
                                  font=('Segoe UI', 10), foreground='gray')
        subtitle_label.pack()
        
        # Status frame
        status_frame = ttk.LabelFrame(parent, text="System Status", padding=15)
        status_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # ML Status
        ml_status_frame = ttk.Frame(status_frame)
        ml_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(ml_status_frame, text="LSTM Model Status:").pack(side=tk.LEFT)
        
        self.ml_status_label = ttk.Label(ml_status_frame, text="", font=('Segoe UI', 10, 'bold'))
        if self.ml_available:
            self.ml_status_label.config(text="Available", foreground='orange')
        else:
            self.ml_status_label.config(text="Not Available", foreground='red')
        self.ml_status_label.pack(side=tk.LEFT, padx=10)
        
        # Isolation Forest Status (NEW)
        if_status_frame = ttk.Frame(status_frame)
        if_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(if_status_frame, text="Isolation Forest Status:").pack(side=tk.LEFT)
        
        self.if_status_label = ttk.Label(if_status_frame, text="", font=('Segoe UI', 10, 'bold'))
        if self.sklearn_available:
            self.if_status_label.config(text="Available", foreground='green')
        else:
            self.if_status_label.config(text="Not Available (install scikit-learn)", foreground='red')
        self.if_status_label.pack(side=tk.LEFT, padx=10)
        
        # XAI Status
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
        
        # XAI Explanations generated
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
        
        # Generate XAI button
        self.generate_xai_btn = ttk.Button(buttons_frame, text="Generate XAI", 
                                          command=self.generate_xai_explanations,
                                          style="Accent.TButton")
        self.generate_xai_btn.pack(side=tk.LEFT, padx=5)
        self.generate_xai_btn.config(state="disabled")
        
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
        """Create the ML configuration tab content with Isolation Forest support"""
        # ML Model Configuration
        model_frame = ttk.LabelFrame(parent, text="ML Model Selection", padding=15)
        model_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Model type selection (Autoencoder and Isolation Forest only)
        model_type_frame = ttk.Frame(model_frame)
        model_type_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(model_type_frame, text="Select Model Type:", font=('Segoe UI', 10, 'bold')).pack(anchor=tk.W)
        
        self.model_type_var = tk.StringVar(value="Autoencoder")
        model_types = [
            "Autoencoder (LSTM/Deep Learning)",
            "Isolation Forest (Unsupervised)"
        ]
        
        for model_type in model_types:
            rb = ttk.Radiobutton(
                model_type_frame, 
                text=model_type, 
                variable=self.model_type_var, 
                value=model_type.split(" ")[0],
                command=self.on_model_type_changed
            )
            rb.pack(anchor=tk.W, padx=(20, 0), pady=2)
        
        # Model file selection (only for Autoencoder)
        file_frame = ttk.Frame(model_frame)
        file_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(file_frame, text="Base Model File (for LSTM):").pack(side=tk.LEFT)
        
        self.model_path_var = tk.StringVar(value="anomaly_detection_model.keras")
        self.model_path_entry = ttk.Entry(file_frame, textvariable=self.model_path_var, width=50)
        self.model_path_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        self.browse_model_btn = ttk.Button(file_frame, text="Browse", 
                                          command=self.browse_model_file)
        self.browse_model_btn.pack(side=tk.RIGHT, padx=5)
        
        # Isolation Forest Configuration (NEW)
        if_config_frame = ttk.LabelFrame(parent, text="Isolation Forest Configuration", padding=15)
        if_config_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Contamination parameter
        contam_frame = ttk.Frame(if_config_frame)
        contam_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(contam_frame, text="Contamination:").grid(row=0, column=0, sticky=tk.W)
        self.if_contamination_var = tk.DoubleVar(value=0.1)
        ttk.Scale(contam_frame, from_=0.01, to=0.5, variable=self.if_contamination_var, 
                 orient=tk.HORIZONTAL, length=200).grid(row=0, column=1, padx=5, sticky=tk.W)
        self.if_contamination_label = ttk.Label(contam_frame, text="10.0%")
        self.if_contamination_label.grid(row=0, column=2, sticky=tk.W, padx=5)
        ttk.Label(contam_frame, text="(expected % of anomalies)").grid(row=0, column=3, sticky=tk.W, padx=5)
        
        # Update contamination label
        def update_contam_label(val):
            self.if_contamination_label.config(text=f"{float(val)*100:.1f}%")
        
        contam_frame.grid_columnconfigure(1, weight=1)
        self.if_contamination_var.trace('w', lambda *args: update_contam_label(self.if_contamination_var.get()))
        
        # Number of estimators
        estimators_frame = ttk.Frame(if_config_frame)
        estimators_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(estimators_frame, text="N Estimators:").grid(row=0, column=0, sticky=tk.W)
        self.if_estimators_var = tk.IntVar(value=150)
        ttk.Spinbox(estimators_frame, from_=50, to=500, textvariable=self.if_estimators_var, 
                   width=10, increment=50).grid(row=0, column=1, padx=5, sticky=tk.W)
        ttk.Label(estimators_frame, text="(number of isolation trees)").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        # XAI Configuration
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
    
    def on_model_type_changed(self):
        """Handle model type selection change"""
        model_type = self.model_type_var.get()
        self.parent_app.add_log(f"Model type changed to: {model_type}")
        
        # Update model info label
        if model_type == "Autoencoder":
            self.model_info_label.config(
                text="LSTM-based autoencoder model for supervised anomaly detection\nRequires .keras model file", 
                foreground='blue'
            )
            self.model_path_entry.config(state="normal")
            self.browse_model_btn.config(state="normal")
        elif model_type == "Isolation":
            self.model_info_label.config(
                text="Unsupervised Isolation Forest - no model file needed\nTrained on normal telemetry data", 
                foreground='green'
            )
            self.model_path_entry.config(state="disabled")
            self.browse_model_btn.config(state="disabled")
    
    def _create_xai_tab(self, parent):
        """Create the XAI explanations tab"""
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
        """Generate LIME and/or SHAP explanations for detected faults"""
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
                                  "pip install shap lime\n\n" +
                                  "Then restart the application.")
            return
        
        try:
            self.parent_app.add_log("Generating XAI explanations...")
            self.generate_xai_btn.config(state="disabled", text="Generating...")
            self.parent_app.root.update()
            
            # Import XAI integration module
            from satellite_xai import add_xai_to_basilisk_detection
            real_telemetry = self.fault_detection_results.get("real_telemetry")
            if real_telemetry is None:
                self.parent_app.add_log("XAI aborted: real_telemetry is None (mock results).")
                messagebox.showwarning(
                    "No Telemetry for XAI",
                    "XAI needs real telemetry, but this run used mock results.\n"
                    "Run a simulation with the Fault Detection tab attached to the pipeline, then try again."
                )
                self.generate_xai_btn.config(state="normal", text="Generate XAI")
                return
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
                success_msg = f"XAI Explanations Generated!\n\n"
                success_msg += f"LIME Explanations: {num_lime}\n"
                success_msg += f"SHAP Analysis: {'Yes' if has_shap else 'No'}\n\n"
                success_msg += f"Results are now displayed in the XAI Explanations tab.\n\n"
                success_msg += f"Plots also saved to:\n"
                if num_lime > 0:
                    success_msg += f"  - lime_anomaly_1.png\n"
                if has_shap:
                    success_msg += f"  - shap_summary.png\n"
                    success_msg += f"  - shap_waterfall_anomaly_1.png\n"
                
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
            # Check what to display
            xai_type_upper = (xai_type or "").upper()
            show_lime = ("LIME" in xai_type_upper) or (xai_type_upper == "BOTH")
            show_shap = ("SHAP" in xai_type_upper) or (xai_type_upper == "BOTH")

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
    
    def _create_analysis_tab(self, parent):
        """Create the GNN Autoencoder Analysis tab with OSHAD-CSL implementation"""
        # Title
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(title_frame, text="GNN Autoencoder Analysis (OSHAD-CSL)", 
                 font=('Segoe UI', 14, 'bold')).pack()
        ttk.Label(title_frame, text="Online Spacecraft Health Anomaly Detection with Causal Structure Learning",
                 font=('Segoe UI', 10), foreground='gray').pack()
        
        # Configuration frame
        config_frame = ttk.LabelFrame(parent, text="GNN Configuration", padding=10)
        config_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Parameters
        param_grid = ttk.Frame(config_frame)
        param_grid.pack(fill=tk.X)
        
        ttk.Label(param_grid, text="Sequence Length:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.gnn_seq_len_var = tk.IntVar(value=10)
        ttk.Spinbox(param_grid, from_=5, to=50, textvariable=self.gnn_seq_len_var, 
                   width=10).grid(row=0, column=1, sticky='w', padx=5, pady=5)
        
        ttk.Label(param_grid, text="Hidden Dimensions:").grid(row=0, column=2, sticky='w', padx=5, pady=5)
        self.gnn_hidden_dim_var = tk.IntVar(value=16)
        ttk.Spinbox(param_grid, from_=8, to=64, textvariable=self.gnn_hidden_dim_var, 
                   width=10).grid(row=0, column=3, sticky='w', padx=5, pady=5)
        
        ttk.Label(param_grid, text="Epochs:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        self.gnn_epochs_var = tk.IntVar(value=300)
        ttk.Spinbox(param_grid, from_=100, to=1000, textvariable=self.gnn_epochs_var, 
                   width=10, increment=50).grid(row=1, column=1, sticky='w', padx=5, pady=5)
        
        ttk.Label(param_grid, text="Learning Rate:").grid(row=1, column=2, sticky='w', padx=5, pady=5)
        self.gnn_lr_var = tk.DoubleVar(value=0.001)
        ttk.Entry(param_grid, textvariable=self.gnn_lr_var, width=10).grid(row=1, column=3, sticky='w', padx=5, pady=5)
        
        ttk.Label(param_grid, text="POT Quantile:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        self.gnn_pot_q_var = tk.DoubleVar(value=0.01)
        ttk.Entry(param_grid, textvariable=self.gnn_pot_q_var, width=10).grid(row=2, column=1, sticky='w', padx=5, pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(config_frame)
        button_frame.pack(pady=10)
        
        ttk.Button(button_frame, text="Train GNN Model", 
                  command=self.train_gnn_model).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Run Analysis", 
                  command=self.run_gnn_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Results", 
                  command=self.export_gnn_results).pack(side=tk.LEFT, padx=5)
        
        # Status label
        self.gnn_status_label = ttk.Label(config_frame, text="Ready to train", foreground='gray')
        self.gnn_status_label.pack()
        
        # Create notebook for visualizations
        viz_notebook = ttk.Notebook(parent)
        viz_notebook.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Anomaly Detection Plot
        anomaly_frame = ttk.Frame(viz_notebook)
        viz_notebook.add(anomaly_frame, text="Anomaly Detection")
        
        self.gnn_anomaly_figure = Figure(figsize=(10, 5))
        self.gnn_anomaly_canvas = FigureCanvasTkAgg(self.gnn_anomaly_figure, anomaly_frame)
        self.gnn_anomaly_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Causal Graph Plot
        causal_frame = ttk.Frame(viz_notebook)
        viz_notebook.add(causal_frame, text="Causal Graph")
        
        self.gnn_causal_figure = Figure(figsize=(10, 6))
        self.gnn_causal_canvas = FigureCanvasTkAgg(self.gnn_causal_figure, causal_frame)
        self.gnn_causal_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Training History Plot
        training_frame = ttk.Frame(viz_notebook)
        viz_notebook.add(training_frame, text="Training History")
        
        self.gnn_training_figure = Figure(figsize=(10, 5))
        self.gnn_training_canvas = FigureCanvasTkAgg(self.gnn_training_figure, training_frame)
        self.gnn_training_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Results text
        results_frame = ttk.Frame(viz_notebook)
        viz_notebook.add(results_frame, text="Detailed Results")
        
        text_scroll = ttk.Scrollbar(results_frame)
        text_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.gnn_results_text = tk.Text(results_frame, wrap=tk.WORD, yscrollcommand=text_scroll.set)
        self.gnn_results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        text_scroll.config(command=self.gnn_results_text.yview)
    
    def check_for_default_model(self):
        """Check if default model file exists and set it"""
        possible_paths = [
            "anomaly_detection_model.keras",
            "../anomaly_detection_model.keras",
            "../../anomaly_detection_model.keras",
            os.path.join(os.path.dirname(__file__), "anomaly_detection_model.keras"),
            os.path.join(os.path.dirname(__file__), "..", "anomaly_detection_model.keras"),
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
        """Load the ML model based on selected type - UPDATED WITH ISOLATION FOREST"""
        model_type = self.model_type_var.get()
        
        # Handle Isolation Forest separately
        if model_type == "Isolation":
            if not self.sklearn_available:
                messagebox.showerror(
                    "scikit-learn Not Available",
                    "Isolation Forest requires scikit-learn.\n\n" +
                    "Install with: pip install scikit-learn\n\n" +
                    "Then restart the application."
                )
                return
            
            try:
                from isolation_forest_fault_detection import SatelliteIsolationForestDetector
                
                self.parent_app.add_log(f"Initializing Isolation Forest detector...")
                self.model_info_label.config(text="Initializing Isolation Forest...", foreground='orange')
                self.parent_app.root.update()
                
                # Create Isolation Forest detector with configured parameters
                self.isolation_forest_detector = SatelliteIsolationForestDetector(
                    contamination=self.if_contamination_var.get(),
                    n_estimators=self.if_estimators_var.get(),
                    random_state=42
                )
                
                self.parent_app.add_log("Isolation Forest detector initialized successfully")
                
                # Update UI
                status_msg = "Isolation Forest Ready!"
                self.model_info_label.config(
                    text=status_msg + f"\n" +
                         f"Algorithm: Unsupervised Anomaly Detection\n" +
                         f"Contamination: {self.if_contamination_var.get()*100:.1f}%\n" +
                         f"N Estimators: {self.if_estimators_var.get()}\n" +
                         f"Note: Will train on normal spacecraft telemetry", 
                    foreground='green'
                )
                
                self.if_status_label.config(text="Loaded", foreground='green')
                self.parent_app.add_log("Isolation Forest configuration complete")
                self.parent_app.update_status_counts()
                
                success_msg = f"Isolation Forest Initialized!\n\n"
                success_msg += f"Algorithm: Unsupervised Anomaly Detection\n"
                success_msg += f"Contamination: {self.if_contamination_var.get()*100:.1f}% (expected anomalies)\n"
                success_msg += f"Number of Trees: {self.if_estimators_var.get()}\n\n"
                success_msg += f"How it works:\n"
                success_msg += f"1. Trains on normal (fault-free) spacecraft telemetry\n"
                success_msg += f"2. Detects anomalies by measuring isolation difficulty\n"
                success_msg += f"3. No labeled data required\n"
                success_msg += f"4. Extracts 22 telemetry features automatically\n\n"
                success_msg += f"Features include:\n"
                success_msg += f"- RW speeds, torques, derivatives\n"
                success_msg += f"- Attitude error and rates\n"
                success_msg += f"- Cross-wheel correlations\n"
                success_msg += f"- System momentum and power\n\n"
                success_msg += f"Ready for detection!"
                
                messagebox.showinfo("Isolation Forest Ready", success_msg)
                return
                
            except ImportError as e:
                error_msg = f"Isolation Forest module not found: {str(e)}\n\nMake sure isolation_forest_fault_detection.py is in the same directory."
                self.model_info_label.config(text="Isolation Forest module not found", foreground='red')
                self.parent_app.add_log(error_msg)
                messagebox.showerror("Import Error", error_msg)
                return
                
            except Exception as e:
                error_msg = f"Error initializing Isolation Forest: {str(e)}"
                self.model_info_label.config(text=error_msg, foreground='red')
                self.parent_app.add_log(error_msg)
                messagebox.showerror("Error", error_msg)
                return
        
        # For LSTM/Autoencoder models, check for model file
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
                if model_type in ["LIME", "SHAP", "Ensemble"] and self.xai_available:
                    self.parent_app.add_log(f"Base model loaded, initializing {model_type} explainer...")
                    
                    try:
                        from satellite_xai import SatelliteXAIExplainer
                        
                        self.xai_explainer = SatelliteXAIExplainer(self.ml_detector.model)
                        self.parent_app.add_log(f"{model_type} explainer initialized successfully")
                        
                        xai_info = f"\n\nXAI Explainer: {model_type}\n"
                        xai_info += f"Base Detector: Autoencoder\n"
                        if model_type in ["LIME", "Ensemble"]:
                            xai_info += f"LIME Features: {self.lime_features_var.get()}\n"
                        if model_type in ["SHAP", "Ensemble"]:
                            xai_info += f"SHAP Type: {self.shap_type_var.get()}\n"
                            xai_info += f"SHAP Samples: {self.shap_samples_var.get()}\n"
                        xai_info += "\nNote: LIME/SHAP wrap around the base autoencoder"
                    except ImportError as e:
                        xai_info = f"\n\nXAI libraries not available\nInstall: pip install shap lime"
                        self.xai_explainer = None
                        self.parent_app.add_log(f"Warning: Could not initialize XAI - {e}")
                else:
                    xai_info = ""
                    self.xai_explainer = None
                
                # Create status message based on model type
                if model_type == "Autoencoder":
                    status_msg = f"Autoencoder loaded successfully!"
                elif model_type in ["LIME", "SHAP"]:
                    status_msg = f"{model_type} Explainer Ready!"
                else:  # Ensemble
                    status_msg = f"Full XAI Suite Ready!"
                
                self.model_info_label.config(
                    text=status_msg + f"\n" +
                         f"Base Model: {os.path.basename(model_path)}\n" +
                         f"Input: {self.ml_detector.model.input_shape}\n" +
                         f"Parameters: {self.ml_detector.model.count_params():,}" +
                         xai_info, 
                    foreground='green'
                )
                
                # Update status label
                if model_type == "Ensemble":
                    self.ml_status_label.config(text="Autoencoder + XAI", foreground='green')
                elif model_type == "LIME":
                    self.ml_status_label.config(text="Autoencoder + LIME", foreground='green')
                elif model_type == "SHAP":
                    self.ml_status_label.config(text="Autoencoder + SHAP", foreground='green')
                else:
                    self.ml_status_label.config(text="Autoencoder Only", foreground='green')
                
                self.parent_app.add_log(f"{model_type} configuration loaded successfully")
                
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
                    success_msg += f"Autoencoder will detect faults\n"
                    success_msg += f"{model_type} will explain WHY faults were detected\n\n"
                    if model_type == "LIME":
                        success_msg += f"LIME explains individual detections\n"
                        success_msg += f"Top {self.lime_features_var.get()} features per anomaly"
                    else:
                        success_msg += f"SHAP shows global feature importance\n"
                        success_msg += f"Analyzing up to {self.shap_samples_var.get()} samples"
                else:  # Ensemble
                    success_msg = f"Full XAI Suite Loaded!\n\n"
                    success_msg += f"Base Model: {os.path.basename(model_path)}\n"
                    success_msg += f"Input shape: {self.ml_detector.model.input_shape}\n"
                    success_msg += f"Parameters: {self.ml_detector.model.count_params():,}\n\n"
                    success_msg += f"Autoencoder: Detects faults\n"
                    success_msg += f"LIME: Explains individual detections\n"
                    success_msg += f"SHAP: Shows global feature importance\n\n"
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
        """Start fault detection process - UPDATED TO SUPPORT ISOLATION FOREST"""
        
        model_type = self.model_type_var.get()
        
        # Check if appropriate model is loaded
        if model_type == "Isolation":
            if not hasattr(self, 'isolation_forest_detector') or self.isolation_forest_detector is None:
                response = messagebox.askyesno(
                    "Isolation Forest Not Ready", 
                    "Isolation Forest detector is not initialized.\n\nLoad it first?"
                )
                if response:
                    self.load_ml_model()
                return
        else:
            # Check LSTM/Autoencoder
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
                          self.threshold_detection_var.get() or
                          model_type == "Isolation")
        
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
        if model_type == "Isolation":
            enabled_methods.append("Isolation Forest (Unsupervised)")
        elif self.ml_detection_var.get():
            enabled_methods.append(f"ML {model_type}")
        if self.statistical_detection_var.get():
            enabled_methods.append("Statistical Analysis")
        if self.trend_detection_var.get():
            enabled_methods.append("Trend Analysis")
        if self.threshold_detection_var.get():
            enabled_methods.append("Threshold Detection")
        
        # Show configuration summary
        config_summary = f"Detection Configuration:\n\n"
        config_summary += f"Primary Method: {model_type}\n"
        config_summary += f"Enabled Methods: {', '.join(enabled_methods)}\n"
        
        if model_type == "Isolation":
            config_summary += f"IF Contamination: {self.if_contamination_var.get()*100:.1f}%\n"
            config_summary += f"IF Estimators: {self.if_estimators_var.get()}\n"
            config_summary += f"Status: {'Ready' if self.isolation_forest_detector else 'Not Loaded'}\n"
        else:
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
        
        self.parent_app.add_log(f"Fault detection configured ({model_type}) and ready")
        self.detection_status_label.config(text="Ready for Simulation", foreground='green')
        
        # Update buttons
        self.start_detection_btn.config(text="Detection Ready", state="disabled")
        
        # Show configuration summary
        messagebox.showinfo("Detection Ready", config_summary)
        
        self.parent_app.add_log("Fault detection ready - run simulation to start detection")
    
    def simulate_detection_results(self):
        """Simulate detection results for testing"""
        self.parent_app.add_log("Simulating detection results...")
        messagebox.showinfo("Simulation", "Detection simulation not yet implemented in this version")
    
    def reset_detection_state(self):
        """Reset detection to initial state"""
        self.detection_status_label.config(text="Idle", foreground='black')
        self.start_detection_btn.config(text="Start Detection", state="normal")
        self.generate_xai_btn.config(state="disabled")
        self.parent_app.add_log("Detection state reset")
    
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
        pass
    
    def display_ml_results(self, ml_results):
        """Display ML detection results (LSTM or Isolation Forest) and enable XAI button"""
        if not ml_results:
            return

        def _update_gui():
            try:
                self.fault_detection_results = ml_results

                # Which detector?
                detector_type = "Unknown"
                if 'detector' in ml_results:
                    detector_type = "Isolation Forest"
                elif 'ml_detector' in ml_results:
                    detector_type = "LSTM/Autoencoder"

                # Status
                self.detection_status_label.config(
                    text=f"Detection Complete ({detector_type})", foreground='green'
                )

                # Clear tables/lists
                for item in self.results_tree.get_children():
                    self.results_tree.delete(item)
                self.recent_detections_listbox.delete(0, tk.END)

                # Pull data
                summary = ml_results.get('summary', {}) or {}
                detections = ml_results.get('detections', {}) or {}

                total_detections = 0
                recent_detections = []

                # Rows for table + recent list
                for spacecraft_name, spacecraft_detections in detections.items():
                    for detection in spacecraft_detections:
                        # Build method label per-detector safely
                        if detector_type == "Isolation Forest":
                            method_used = f"Isolation Forest ({getattr(detection, 'fault_type', 'unknown')})"
                        else:
                            details = getattr(detection, 'details', {}) or {}
                            safe_method = details.get('primary_method', 'ML') if isinstance(details, dict) else 'ML'
                            method_used = f"{safe_method}"

                        # Insert into table
                        self.results_tree.insert(
                            '',
                            'end',
                            values=(
                                f"{getattr(detection, 'detection_time_minutes', 0.0):.1f} min",
                                spacecraft_name,
                                getattr(detection, 'fault_type', 'unknown'),
                                f"{getattr(detection, 'confidence', 0.0):.3f}",
                                getattr(detection, 'affected_component', '—'),
                                method_used,
                                "Active" if getattr(detection, 'fault_detected', False) else "Inactive",
                            ),
                        )

                        # Track for "Recent Detections"
                        recent_detections.append(
                            {
                                'time': f"{getattr(detection, 'detection_time_minutes', 0.0):.1f} min",
                                'spacecraft': spacecraft_name,
                                'type': getattr(detection, 'fault_type', 'unknown'),
                                'confidence': getattr(detection, 'confidence', 0.0),
                                'method': method_used,
                                'status': 'Active' if getattr(detection, 'fault_detected', False) else 'Inactive',
                            }
                        )
                        total_detections += 1

                # Recent list UI
                if recent_detections:
                    for d in recent_detections[-10:]:
                        txt = f"{d['time']} - {d['spacecraft']}: {d['method']} (conf: {d['confidence']:.3f})"
                        self.recent_detections_listbox.insert(0, txt)
                else:
                    self.recent_detections_listbox.insert(0, "No detections found")

                # History + stats
                self.detection_history.extend(recent_detections)
                self.update_statistics()

                # Build summary AFTER we've done the table, so summary_text exists before use
                self.results_summary_text.delete(1.0, tk.END)
                summary_text = (
                    f"FAULT DETECTION RESULTS ({detector_type}):\n"
                    f"========================================\n"
                    f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                    f"Detection Method: {detector_type}\n\n"
                    f"SUMMARY:\n"
                    f"  Spacecraft Analyzed: {summary.get('total_spacecraft', 0)}\n"
                    f"  Total Detections: {summary.get('total_detections', 0)}\n"
                    f"  Success Rate: {summary.get('success_rate', 1.0):.1%}\n\n"
                )

                if detector_type == "Isolation Forest":
                    summary_text += "ALGORITHM DETAILS:\n"
                    summary_text += "  Method: Unsupervised Anomaly Detection\n"
                    summary_text += "  Features: 22 telemetry parameters\n"
                    if 'confidence_scores' in summary:
                        try:
                            import numpy as _np
                            summary_text += f"  Avg Confidence: {_np.mean(summary['confidence_scores']):.3f}\n"
                        except Exception:
                            pass
                    if 'isolation_scores' in summary:
                        try:
                            import numpy as _np
                            summary_text += f"  Avg Isolation Score: {_np.mean(summary['isolation_scores']):.3f}\n"
                        except Exception:
                            pass
                    summary_text += "\n"

                summary_text += "DETECTIONS BY SPACECRAFT:\n"

                for spacecraft, spacecraft_detections in detections.items():
                    summary_text += f"\n{spacecraft}: {len(spacecraft_detections)} detections"
                    for detection in spacecraft_detections[:3]:
                        if detector_type == "Isolation Forest":
                            iso = getattr(detection, 'isolation_score', 0.0)
                            summary_text += (
                                f"\n  - {getattr(detection, 'detection_time_minutes', 0.0):.1f}min: "
                                f"{getattr(detection, 'fault_type', 'unknown')} "
                                f"(conf: {getattr(detection, 'confidence', 0.0):.3f}, iso: {iso:.3f})"
                            )
                            # Optional per-feature importance if provided
                            feat_imp = getattr(detection, 'feature_importance', None)
                            if isinstance(feat_imp, dict) and feat_imp:
                                top_features = list(feat_imp.items())[:3]
                                summary_text += f"\n    Top features: {', '.join([k for k, _ in top_features])}"
                        else:
                            details = getattr(detection, 'details', {}) or {}
                            safe_method = details.get('primary_method', 'ML') if isinstance(details, dict) else 'ML'
                            summary_text += (
                                f"\n  - {getattr(detection, 'detection_time_minutes', 0.0):.1f}min via {safe_method} "
                                f"(conf: {getattr(detection, 'confidence', 0.0):.3f})"
                            )
                            # Optional criteria lines
                            if isinstance(details, dict):
                                criteria = details.get('detection_criteria', {}) or {}
                                if criteria.get('rw_speed_change'):
                                    summary_text += f"\n    * RW speed change: {criteria.get('max_speed_change_percent', 0):.1f}%"
                                if criteria.get('attitude_change'):
                                    summary_text += f"\n    * Attitude error increase: {criteria.get('attitude_change_percent', 0):.1f}%"

                if total_detections == 0:
                    summary_text += "\n\nNOTE: No faults were detected."
                else:
                    summary_text += f"\n\nSUCCESS: Detected {total_detections} faults using {detector_type}!"
                    # Enable XAI only if we have detections
                    if hasattr(self, 'generate_xai_btn'):
                        self.generate_xai_btn.config(state="normal")
                        summary_text += "\n\nXAI explanations available — click 'Generate XAI'"

                self.results_summary_text.insert(tk.END, summary_text)

                # Show results tab
                self.detection_notebook.select(2)

                # Log
                self.parent_app.add_log(
                    f"{detector_type} Detection Results: {total_detections} detections from {len(detections)} spacecraft"
                )

            except Exception as e:
                self.parent_app.add_log(f"Error displaying results: {e}")
                import traceback
                traceback.print_exc()

        # Ensure main thread
        try:
            import threading
            if threading.current_thread() is threading.main_thread():
                _update_gui()
            else:
                self.parent_app.root.after(0, _update_gui)
        except Exception:
            try:
                self.parent_app.root.after(0, _update_gui)
            except Exception:
                _update_gui()

    
    def train_gnn_model(self):
        """Train the GNN Autoencoder model"""
        try:
            import torch
            import torch.nn as nn
            import torch.optim as optim
            from sklearn.preprocessing import StandardScaler
            from scipy.stats import genpareto
            import networkx as nx
        except ImportError as e:
            messagebox.showerror("Error", f"Required library not available: {e}\n\nPlease install: torch, scikit-learn, scipy, networkx")
            return
        
        # Check if we have data
        if not hasattr(self.parent_app, 'processed_data') or self.parent_app.processed_data is None:
            messagebox.showwarning("No Data", "Please load and process spacecraft data first.")
            return
        
        try:
            self.gnn_status_label.config(text="Training GNN model...", foreground='orange')
            self.parent_app.add_log("Starting GNN Autoencoder training...")
            self.parent_app.root.update()
            
            # Get parameters
            seq_len = self.gnn_seq_len_var.get()
            hidden_dim = self.gnn_hidden_dim_var.get()
            epochs = self.gnn_epochs_var.get()
            lr = self.gnn_lr_var.get()
            
            # Prepare data
            df = self.parent_app.processed_data
            
            # Select numerical columns
            numerical_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            if 'spacecraft_id' in numerical_cols:
                numerical_cols.remove('spacecraft_id')
            
            # Normalize data
            scaler = StandardScaler()
            data = scaler.fit_transform(df[numerical_cols].values)
            
            n_features = data.shape[1]
            device = "cuda" if torch.cuda.is_available() else "cpu"
            
            self.parent_app.add_log(f"Data prepared: {len(data)} samples, {n_features} features")
            self.parent_app.add_log(f"Using device: {device}")
            
            # Convert to sliding windows
            X = []
            y = []
            for i in range(seq_len, len(data)):
                X.append(data[i - 1])
                y.append(data[i])
            
            X = np.array(X)
            y = np.array(y)
            X_tensor = torch.tensor(X, dtype=torch.float32).to(device)
            y_tensor = torch.tensor(y, dtype=torch.float32).to(device)
            
            # Define GNN Autoencoder classes
            class GCNLayer(nn.Module):
                def __init__(self, in_features, out_features):
                    super().__init__()
                    self.linear = nn.Linear(in_features, out_features)

                def forward(self, X, A):
                    XW = self.linear(X)
                    return torch.relu(torch.matmul(XW, A))

            class GNNEncoder(nn.Module):
                def __init__(self, n_features, hidden_dim):
                    super().__init__()
                    self.gcn1 = GCNLayer(n_features, hidden_dim)
                    self.gcn2 = GCNLayer(hidden_dim, hidden_dim)

                def forward(self, X, A):
                    X = self.gcn1(X, A)
                    return self.gcn2(X, A)

            class AutoEncoder(nn.Module):
                def __init__(self, n_features, hidden_dim):
                    super().__init__()
                    self.encoder = GNNEncoder(n_features, hidden_dim)
                    self.decoder = nn.Sequential(
                        nn.Linear(hidden_dim, hidden_dim),
                        nn.ReLU(),
                        nn.Linear(hidden_dim, n_features)
                    )

                def forward(self, X, A):
                    H = self.encoder(X, A)
                    return self.decoder(H)
            
            # Initialize model
            A = torch.eye(hidden_dim, requires_grad=True, device=device)
            model = AutoEncoder(n_features, hidden_dim).to(device)
            optimizer = optim.Adam(list(model.parameters()) + [A], lr=lr)
            loss_fn = nn.MSELoss()
            
            self.parent_app.add_log("Model initialized. Starting training...")
            
            # Training loop
            training_losses = []
            for epoch in range(epochs):
                model.train()
                optimizer.zero_grad()
                output = model(X_tensor, A)
                loss = loss_fn(output, y_tensor)
                loss.backward()
                optimizer.step()
                
                training_losses.append(loss.item())
                
                if epoch % 50 == 0:
                    self.parent_app.add_log(f"Epoch {epoch}, Loss: {loss.item():.4f}")
                    self.gnn_status_label.config(text=f"Training... Epoch {epoch}/{epochs}, Loss: {loss.item():.4f}")
                    self.parent_app.root.update()
            
            # Plot training history
            self.gnn_training_figure.clear()
            ax = self.gnn_training_figure.add_subplot(111)
            ax.plot(training_losses, linewidth=2)
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Loss (MSE)')
            ax.set_title('GNN Autoencoder Training History')
            ax.grid(True, alpha=0.3)
            self.gnn_training_canvas.draw()
            
            # Save trained model
            self.gnn_model = {
                'model': model,
                'A': A,
                'scaler': scaler,
                'device': device,
                'n_features': n_features,
                'hidden_dim': hidden_dim,
                'seq_len': seq_len,
                'numerical_cols': numerical_cols
            }
            
            self.gnn_status_label.config(text="Training completed successfully!", foreground='green')
            self.parent_app.add_log("GNN Autoencoder training completed successfully!")
            
            messagebox.showinfo("Success", f"GNN model trained successfully!\nFinal Loss: {training_losses[-1]:.4f}")
            
        except Exception as e:
            self.gnn_status_label.config(text="Training failed", foreground='red')
            self.parent_app.add_log(f"Error training GNN model: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Training Error", f"Failed to train GNN model:\n{str(e)}")
    
    def run_gnn_analysis(self):
        """Run GNN Autoencoder analysis on the data"""
        try:
            import torch
            import networkx as nx
            from scipy.stats import genpareto
        except ImportError as e:
            messagebox.showerror("Error", f"Required library not available: {e}")
            return
        
        if not hasattr(self, 'gnn_model') or self.gnn_model is None:
            messagebox.showwarning("No Model", "Please train the GNN model first.")
            return
        
        if not hasattr(self.parent_app, 'processed_data') or self.parent_app.processed_data is None:
            messagebox.showwarning("No Data", "Please load and process spacecraft data first.")
            return
        
        try:
            self.gnn_status_label.config(text="Running analysis...", foreground='orange')
            self.parent_app.add_log("Running GNN Autoencoder analysis...")
            self.parent_app.root.update()
            
            # Get model components
            model = self.gnn_model['model']
            A = self.gnn_model['A']
            scaler = self.gnn_model['scaler']
            device = self.gnn_model['device']
            seq_len = self.gnn_model['seq_len']
            numerical_cols = self.gnn_model['numerical_cols']
            
            # Prepare data
            df = self.parent_app.processed_data
            data = scaler.transform(df[numerical_cols].values)
            
            # Convert to sliding windows
            X = []
            y = []
            for i in range(seq_len, len(data)):
                X.append(data[i - 1])
                y.append(data[i])
            
            X = np.array(X)
            y = np.array(y)
            X_tensor = torch.tensor(X, dtype=torch.float32).to(device)
            y_tensor = torch.tensor(y, dtype=torch.float32).to(device)
            
            # Predict and calculate errors
            model.eval()
            with torch.no_grad():
                predictions = model(X_tensor, A).cpu().numpy()
                targets = y_tensor.cpu().numpy()
                errors = np.abs(predictions - targets).mean(axis=1)
            
            # Dynamic Threshold using POT
            def dynamic_threshold(errors, q=0.01, n=100):
                fit_data = errors[:n]
                if len(fit_data) == 0 or all(fit_data <= 0):
                    return np.mean(errors) + 3 * np.std(errors)
                threshold = np.percentile(fit_data, 100 * (1 - q))
                excesses = fit_data[fit_data > threshold] - threshold
                if len(excesses) == 0:
                    return threshold
                try:
                    shape, loc, scale = genpareto.fit(excesses)
                except:
                    pass
                return threshold
            
            pot_q = self.gnn_pot_q_var.get()
            thresh = dynamic_threshold(errors, q=pot_q)
            labels = (errors > thresh).astype(int)
            
            anomaly_count = np.sum(labels)
            anomaly_rate = anomaly_count / len(labels) * 100
            
            # Plot anomaly detection results
            self.gnn_anomaly_figure.clear()
            ax = self.gnn_anomaly_figure.add_subplot(111)
            ax.plot(errors, label="Reconstruction Error", linewidth=2)
            ax.axhline(y=thresh, color='r', linestyle='--', label=f"Dynamic Threshold ({thresh:.4f})", linewidth=2)
            
            # Highlight anomalies
            anomaly_indices = np.where(labels == 1)[0]
            if len(anomaly_indices) > 0:
                ax.scatter(anomaly_indices, errors[anomaly_indices], color='red', s=50, 
                          alpha=0.6, label=f'Anomalies ({anomaly_count})', zorder=5)
            
            ax.set_xlabel('Time Step')
            ax.set_ylabel('Reconstruction Error')
            ax.set_title('GNN Autoencoder Anomaly Detection')
            ax.legend()
            ax.grid(True, alpha=0.3)
            self.gnn_anomaly_canvas.draw()
            
            # Visualize Learned Causal Graph
            adj_matrix = A.detach().cpu().numpy()
            G = nx.DiGraph()
            
            for i in range(adj_matrix.shape[0]):
                for j in range(adj_matrix.shape[1]):
                    if abs(adj_matrix[i, j]) > 0.05:
                        G.add_edge(f"H{i}", f"H{j}", weight=round(adj_matrix[i, j], 2))
            
            self.gnn_causal_figure.clear()
            ax2 = self.gnn_causal_figure.add_subplot(111)
            
            pos = nx.spring_layout(G, seed=42)
            nx.draw(G, pos, ax=ax2, with_labels=True, node_color='lightblue', 
                   node_size=500, font_size=9, edge_color='gray', arrows=True, 
                   arrowsize=15, arrowstyle='->')
            
            edge_labels = nx.get_edge_attributes(G, 'weight')
            nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax2, font_size=7)
            
            ax2.set_title(f"Learned Causal Graph ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
            ax2.axis('off')
            self.gnn_causal_canvas.draw()
            
            # Display detailed results
            self.gnn_results_text.delete(1.0, tk.END)
            results_summary = f"""GNN AUTOENCODER ANALYSIS RESULTS
{'=' * 60}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MODEL CONFIGURATION:
  Sequence Length: {seq_len}
  Hidden Dimensions: {self.gnn_model['hidden_dim']}
  Number of Features: {self.gnn_model['n_features']}
  Device: {device}

ANOMALY DETECTION SUMMARY:
  Total Time Steps: {len(errors)}
  Anomalies Detected: {anomaly_count}
  Anomaly Rate: {anomaly_rate:.2f}%
  Dynamic Threshold: {thresh:.6f}
  POT Quantile: {pot_q}

ERROR STATISTICS:
  Mean Error: {np.mean(errors):.6f}
  Std Error: {np.std(errors):.6f}
  Min Error: {np.min(errors):.6f}
  Max Error: {np.max(errors):.6f}

CAUSAL GRAPH STATISTICS:
  Nodes: {G.number_of_nodes()}
  Edges: {G.number_of_edges()}
  Avg Degree: {2 * G.number_of_edges() / G.number_of_nodes() if G.number_of_nodes() > 0 else 0:.2f}

TOP ANOMALOUS TIME STEPS:
"""
            # Get top 10 anomalies
            top_anomaly_indices = np.argsort(errors)[-10:][::-1]
            for idx in top_anomaly_indices:
                results_summary += f"  Step {idx}: Error = {errors[idx]:.6f} {'[ANOMALY]' if labels[idx] == 1 else ''}\n"
            
            results_summary += f"\nFEATURES ANALYZED:\n"
            for i, col in enumerate(numerical_cols):
                results_summary += f"  {i+1}. {col}\n"
            
            self.gnn_results_text.insert(tk.END, results_summary)
            
            self.gnn_status_label.config(text="Analysis completed successfully!", foreground='green')
            self.parent_app.add_log(f"GNN Analysis completed: {anomaly_count} anomalies detected ({anomaly_rate:.2f}%)")
            
            messagebox.showinfo("Analysis Complete", 
                              f"GNN Autoencoder analysis completed!\n\n"
                              f"Anomalies Detected: {anomaly_count}\n"
                              f"Anomaly Rate: {anomaly_rate:.2f}%\n"
                              f"Causal Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
            
        except Exception as e:
            self.gnn_status_label.config(text="Analysis failed", foreground='red')
            self.parent_app.add_log(f"Error during GNN analysis: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Analysis Error", f"Failed to run GNN analysis:\n{str(e)}")
    
    def export_gnn_results(self):
        """Export GNN analysis results"""
        if not hasattr(self, 'gnn_model') or self.gnn_model is None:
            messagebox.showwarning("No Results", "Please run GNN analysis first.")
            return
        
        try:
            filename = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
                title="Export GNN Results"
            )
            
            if filename:
                with open(filename, 'w') as f:
                    f.write(self.gnn_results_text.get(1.0, tk.END))
                
                messagebox.showinfo("Success", f"Results exported to:\n{filename}")
                self.parent_app.add_log(f"GNN results exported to {filename}")
        
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export results:\n{str(e)}")
    
    def is_ml_ready(self):
        """Check if ML detection is ready"""
        return (
            self.ml_available and 
            hasattr(self, 'ml_detector') and 
            self.ml_detector is not None and 
            self.ml_detector.is_loaded
        )
    
    def is_isolation_forest_ready(self):
        """Check if Isolation Forest detection is ready"""
        return (
            self.sklearn_available and
            hasattr(self, 'isolation_forest_detector') and
            self.isolation_forest_detector is not None
        )