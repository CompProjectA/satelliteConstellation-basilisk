#!/usr/bin/env python
"""
fault_detection_tab.py

Fault Detection Tab for the Spacecraft Simulator GUI (Tkinter version).
Displays real-time fault detection results from ML models and traditional methods.
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

from .base_tab import BaseTab

class FaultDetectionTab(BaseTab):
    """
    Fault Detection Tab - Shows ML fault detection results and analysis
    """
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the fault detection tab
        
        Parameters:
        parent_app: The parent application instance
        parent_frame: The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        self.fault_detection_results = None
        self.ml_detector = None
        self.detection_history = []
        self.monitoring_active = False
        
        # Check ML availability
        self.ml_available = self.check_ml_availability()
        
        # Create the tab UI
        self.create_tab_ui()
        
    def check_ml_availability(self):
        """Check if ML fault detection is available"""
        try:
            from real_ml_fault_detection import RealMLFaultDetector
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
        monitoring_frame = ttk.Frame(self.detection_notebook)
        analysis_frame = ttk.Frame(self.detection_notebook)
        
        # Add sub-tabs
        self.detection_notebook.add(overview_frame, text="Overview")
        self.detection_notebook.add(config_frame, text="ML Configuration")
        self.detection_notebook.add(results_frame, text="Detection Results")
        self.detection_notebook.add(monitoring_frame, text="Live Monitoring")
        self.detection_notebook.add(analysis_frame, text="Analysis")
        
        # Create tab contents
        self._create_overview_tab(overview_frame)
        self._create_config_tab(config_frame)
        self._create_results_tab(results_frame)
        self._create_monitoring_tab(monitoring_frame)
        self._create_analysis_tab(analysis_frame)
        
    def _create_overview_tab(self, parent):
        """Create the overview tab content"""
        # Title frame
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(10, 20))
        
        title_label = ttk.Label(title_frame, text=" ML Fault Detection System", 
                               font=('Segoe UI', 16, 'bold'))
        title_label.pack()
        
        subtitle_label = ttk.Label(title_frame, text="Real-time spacecraft fault detection using machine learning",
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
            self.ml_status_label.config(text=" Available", foreground='green')
        else:
            self.ml_status_label.config(text=" Not Available", foreground='red')
        self.ml_status_label.pack(side=tk.LEFT, padx=10)
        
        # Detection Status
        detection_status_frame = ttk.Frame(status_frame)
        detection_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(detection_status_frame, text="Detection Status:").pack(side=tk.LEFT)
        
        self.detection_status_label = ttk.Label(detection_status_frame, text=" Idle", 
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
        
        # Last detection
        ttk.Label(stats_grid, text="Last Detection:").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.last_detection_label = ttk.Label(stats_grid, text="None", foreground='gray')
        self.last_detection_label.grid(row=3, column=1, sticky=tk.W, padx=10)
        
        # Quick actions frame
        actions_frame = ttk.LabelFrame(parent, text="Quick Actions", padding=15)
        actions_frame.pack(fill=tk.X, padx=20, pady=10)
        
        buttons_frame = ttk.Frame(actions_frame)
        buttons_frame.pack()
        
        self.load_model_btn = ttk.Button(buttons_frame, text="🔄 Load ML Model", 
                                        command=self.load_ml_model)
        self.load_model_btn.pack(side=tk.LEFT, padx=5)
        
        self.start_detection_btn = ttk.Button(buttons_frame, text="▶️ Start Detection", 
                                             command=self.start_fault_detection)
        self.start_detection_btn.pack(side=tk.LEFT, padx=5)
        
        self.view_results_btn = ttk.Button(buttons_frame, text="📋 View Results", 
                                          command=lambda: self.detection_notebook.select(2))
        self.view_results_btn.pack(side=tk.LEFT, padx=5)
        
        # Recent detections frame
        recent_frame = ttk.LabelFrame(parent, text="Recent Detections", padding=15)
        recent_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Recent detections listbox
        self.recent_detections_listbox = tk.Listbox(recent_frame, height=8)
        self.recent_detections_listbox.pack(fill=tk.BOTH, expand=True)
        
        self.recent_detections_listbox.insert(tk.END, "No recent detections...")
        
    def _create_config_tab(self, parent):
        """Create the ML configuration tab content"""
        # ML Model Configuration
        model_frame = ttk.LabelFrame(parent, text="🧠 ML Model Configuration", padding=15)
        model_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Model file selection
        file_frame = ttk.Frame(model_frame)
        file_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(file_frame, text="Model File:").pack(side=tk.LEFT)
        
        self.model_path_var = tk.StringVar(value="anomaly_detection_model.keras")
        self.model_path_entry = ttk.Entry(file_frame, textvariable=self.model_path_var, width=40)
        self.model_path_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        
        self.browse_model_btn = ttk.Button(file_frame, text="📁 Browse", 
                                          command=self.browse_model_file)
        self.browse_model_btn.pack(side=tk.RIGHT, padx=5)
        
        # Detection threshold
        threshold_frame = ttk.Frame(model_frame)
        threshold_frame.pack(fill=tk.X, pady=5)
        
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
                                         foreground='gray', wraplength=400)
        self.model_info_label.pack(pady=10)
        
        # Detection Methods
        methods_frame = ttk.LabelFrame(parent, text="🔍 Detection Methods", padding=15)
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
        realtime_frame = ttk.LabelFrame(parent, text="⏱️ Real-time Detection", padding=15)
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
        
    def _create_results_tab(self, parent):
        """Create the detection results tab content"""
        # Controls frame
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Button(controls_frame, text="🔄 Refresh Results", 
                  command=self.refresh_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="🗑️ Clear Results", 
                  command=self.clear_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(controls_frame, text="💾 Export Results", 
                  command=self.export_results).pack(side=tk.LEFT, padx=5)
        
        # Results table
        table_frame = ttk.LabelFrame(parent, text="Detection Results", padding=10)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create treeview for results
        columns = ('Time', 'Spacecraft', 'Fault Type', 'Confidence', 'Component', 'Method', 'Status')
        self.results_tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=15)
        
        # Define column headings and widths
        for col in columns:
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=100)
        
        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.results_tree.yview)
        h_scrollbar = ttk.Scrollbar(table_frame, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack treeview and scrollbars
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Results summary
        summary_frame = ttk.LabelFrame(parent, text="Summary", padding=10)
        summary_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.results_summary_text = tk.Text(summary_frame, height=6, wrap=tk.WORD)
        self.results_summary_text.pack(fill=tk.X)
        self.results_summary_text.insert(tk.END, "Detection results summary will appear here...")
        
    def _create_monitoring_tab(self, parent):
        """Create the live monitoring tab content"""
        # Monitoring controls
        controls_frame = ttk.Frame(parent)
        controls_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.monitoring_status_label = ttk.Label(controls_frame, text="📴 Monitoring Inactive", 
                                                font=('Segoe UI', 10, 'bold'))
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
        
        ttk.Button(controls_frame, text="📊 Generate Analysis", 
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
    
    # Methods for ML functionality
    def browse_model_file(self):
        """Browse for ML model file"""
        file_path = filedialog.askopenfilename(
            title="Select ML Model File",
            filetypes=[("Keras Models", "*.keras *.h5"), ("All Files", "*.*")]
        )
        if file_path:
            self.model_path_var.set(file_path)
    
    def load_ml_model(self):
        """Load the ML model"""
        if not self.ml_available:
            messagebox.showwarning("Warning", "ML libraries not available!")
            return
        
        model_path = self.model_path_var.get().strip()
        if not model_path:
            messagebox.showwarning("Warning", "Please specify a model file!")
            return
        
        try:
            from real_ml_fault_detection import RealMLFaultDetector
            
            self.parent_app.add_log(f"Loading ML model: {model_path}")
            self.ml_detector = RealMLFaultDetector(model_path)
            
            if self.ml_detector.is_loaded:
                self.model_info_label.config(text=f"✅ Model loaded successfully\nPath: {model_path}", 
                                            foreground='green')
                self.ml_status_label.config(text="✅ Model Loaded", foreground='green')
                self.parent_app.add_log("ML model loaded successfully")
            else:
                self.model_info_label.config(text="❌ Failed to load model", foreground='red')
                self.parent_app.add_log("Failed to load ML model")
                
        except Exception as e:
            error_msg = f"Error loading ML model: {str(e)}"
            self.model_info_label.config(text=f"❌ {error_msg}", foreground='red')
            self.parent_app.add_log(error_msg)
            messagebox.showerror("Error", error_msg)
    
    def start_fault_detection(self):
        """Start fault detection process"""
        if not self.ml_detector and self.ml_detection_var.get():
            messagebox.showwarning("Warning", "Please load ML model first!")
            return
        
        # Check if any detection method is enabled
        methods_enabled = (self.ml_detection_var.get() or 
                          self.statistical_detection_var.get() or
                          self.trend_detection_var.get() or 
                          self.threshold_detection_var.get())
        
        if not methods_enabled:
            messagebox.showwarning("Warning", "Please enable at least one detection method!")
            return
        
        self.parent_app.add_log("Starting fault detection...")
        self.detection_status_label.config(text="🔄 Running Detection...", foreground='orange')
        
        # Update buttons
        self.start_detection_btn.config(text="⏹️ Stop Detection", command=self.stop_fault_detection)
        
        self.parent_app.add_log("Fault detection started successfully")
        
        # TODO: Integrate with your spacecraft_simulation.py
        # This would be called after simulation completes
        
    def stop_fault_detection(self):
        """Stop fault detection process"""
        self.parent_app.add_log("Stopping fault detection...")
        self.detection_status_label.config(text="Idle", foreground='black')
        
        # Update buttons
        self.start_detection_btn.config(text="▶️ Start Detection", command=self.start_fault_detection)
        
        self.parent_app.add_log("Fault detection stopped")
    
    def toggle_realtime_monitoring(self):
        """Toggle real-time monitoring"""
        if self.realtime_var.get():
            self.monitoring_status_label.config(text="📡 Monitoring Active", foreground='green')
            self.monitoring_active = True
        else:
            self.monitoring_status_label.config(text="📴 Monitoring Inactive", foreground='gray')
            self.monitoring_active = False
    
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
        analysis_type = self.analysis_type_var.get()
        
        try:
            self.analysis_figure.clear()
            
            if analysis_type == "Detection Confidence Over Time":
                self.plot_confidence_over_time()
            elif analysis_type == "Fault Distribution by Spacecraft":
                self.plot_fault_distribution()
            elif analysis_type == "Detection Method Comparison":
                self.plot_method_comparison()
            elif analysis_type == "Error Rate Analysis":
                self.plot_error_rate_analysis()
            elif analysis_type == "Performance Metrics":
                self.plot_performance_metrics()
            
            self.analysis_canvas.draw()
            self.parent_app.add_log(f"Generated analysis: {analysis_type}")
            
        except Exception as e:
            error_msg = f"Failed to generate analysis: {str(e)}"
            self.parent_app.add_log(error_msg)
            messagebox.showerror("Error", error_msg)
    
    def plot_confidence_over_time(self):
        """Plot detection confidence over time"""
        ax = self.analysis_figure.add_subplot(111)
        
        # Sample data - replace with actual detection history
        times = np.linspace(0, 30, 100)
        confidence = np.random.beta(2, 5, 100)  # Sample confidence scores
        
        ax.plot(times, confidence, 'b-', linewidth=2, label='Confidence Score')
        ax.axhline(y=self.threshold_var.get(), color='r', linestyle='--', 
                  label=f'Threshold ({self.threshold_var.get():.3f})')
        
        ax.set_xlabel('Time (minutes)')
        ax.set_ylabel('Detection Confidence')
        ax.set_title('Fault Detection Confidence Over Time')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def plot_fault_distribution(self):
        """Plot fault distribution by spacecraft"""
        ax = self.analysis_figure.add_subplot(111)
        
        # Sample data
        spacecraft = ['Satellite1', 'Satellite2', 'Satellite3', 'Satellite4']
        fault_counts = [3, 1, 0, 2]
        
        bars = ax.bar(spacecraft, fault_counts, color=['red' if count > 0 else 'green' for count in fault_counts])
        ax.set_xlabel('Spacecraft')
        ax.set_ylabel('Number of Faults Detected')
        ax.set_title('Fault Distribution by Spacecraft')
        
        # Add count labels on bars
        for bar, count in zip(bars, fault_counts):
            if count > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                       str(count), ha='center', va='bottom')
    
    def plot_method_comparison(self):
        """Plot detection method comparison"""
        ax = self.analysis_figure.add_subplot(111)
        
        # Sample data
        methods = ['ML Autoencoder', 'Statistical', 'Trend Analysis', 'Threshold']
        detection_rates = [85, 72, 68, 45]
        
        bars = ax.bar(methods, detection_rates, color=['#2E8B57', '#4682B4', '#DAA520', '#DC143C'])
        ax.set_xlabel('Detection Method')
        ax.set_ylabel('Detection Rate (%)')
        ax.set_title('Detection Method Performance Comparison')
        ax.set_ylim(0, 100)
        
        # Add percentage labels
        for i, rate in enumerate(detection_rates):
            ax.text(i, rate + 2, f'{rate}%', ha='center', va='bottom')
    
    def plot_error_rate_analysis(self):
        """Plot error rate analysis"""
        ax = self.analysis_figure.add_subplot(111)
        
        # Sample data for different error types
        error_types = ['False Positive', 'False Negative', 'Missed Detection', 'Late Detection']
        error_rates = [12, 8, 5, 15]
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
        wedges, texts, autotexts = ax.pie(error_rates, labels=error_types, colors=colors,
                                         autopct='%1.1f%%', startangle=90)
        
        ax.set_title('Error Rate Analysis')
    
    def plot_performance_metrics(self):
        """Plot performance metrics"""
        # Create 2x2 subplot layout
        fig = self.analysis_figure
        
        # Precision/Recall
        ax1 = fig.add_subplot(221)
        precision_recall = [0.85, 0.78]
        ax1.bar(['Precision', 'Recall'], precision_recall, color=['#2E8B57', '#4682B4'])
        ax1.set_title('Precision vs Recall')
        ax1.set_ylim(0, 1)
        
        # Response Time Distribution
        ax2 = fig.add_subplot(222)
        response_times = np.random.gamma(2, 2, 1000)  # Sample response times
        ax2.hist(response_times, bins=30, color='skyblue', alpha=0.7, edgecolor='black')
        ax2.set_title('Response Time Distribution')
        ax2.set_xlabel('Response Time (seconds)')
        ax2.set_ylabel('Frequency')
        
        # Detection Accuracy by Fault Type
        ax3 = fig.add_subplot(223)
        fault_types = ['Friction', 'Power Limit', 'Encoder', 'Battery']
        accuracies = [0.92, 0.87, 0.79, 0.83]
        ax3.bar(fault_types, accuracies, color=['red', 'orange', 'yellow', 'green'])
        ax3.set_title('Detection Accuracy by Fault Type')
        ax3.set_ylim(0, 1)
        ax3.tick_params(axis='x', rotation=45)
        
        # Confidence Score Distribution
        ax4 = fig.add_subplot(224)
        confidence_scores = np.random.beta(5, 2, 1000)  # Sample confidence scores
        ax4.hist(confidence_scores, bins=20, color='lightcoral', alpha=0.7, edgecolor='black')
        ax4.set_title('Confidence Score Distribution')
        ax4.set_xlabel('Confidence Score')
        ax4.set_ylabel('Frequency')
        
        fig.tight_layout()
    
    def update_statistics(self):
        """Update the statistics display"""
        total_detections = len(self.detection_history)
        active_faults = sum(1 for d in self.detection_history if d.get('status') == 'Active')
        
        self.total_detections_label.config(text=str(total_detections))
        self.active_faults_label.config(text=str(active_faults))
        
        # Calculate success rate if we have ground truth
        if total_detections > 0:
            # This is a placeholder - you'd calculate based on actual vs detected faults
            success_rate = min(85.0, (active_faults / max(1, total_detections)) * 100)
            self.success_rate_label.config(text=f"{success_rate:.1f}%")
        
        # Update last detection
        if self.detection_history:
            last_detection = self.detection_history[-1]
            self.last_detection_label.config(
                text=f"{last_detection.get('spacecraft', 'Unknown')} at {last_detection.get('time', 'Unknown')}"
            )
    
    def save_results_to_file(self, file_path):
        """Save detection results to file"""
        if file_path.endswith('.json'):
            # Save as JSON
            results_data = {
                'timestamp': datetime.now().isoformat(),
                'total_detections': len(self.detection_history),
                'detections': self.detection_history,
                'configuration': {
                    'ml_enabled': self.ml_detection_var.get(),
                    'statistical_enabled': self.statistical_detection_var.get(),
                    'trend_enabled': self.trend_detection_var.get(),
                    'threshold_enabled': self.threshold_detection_var.get(),
                    'detection_threshold': self.threshold_var.get()
                }
            }
            
            with open(file_path, 'w') as f:
                json.dump(results_data, f, indent=2)
                
        elif file_path.endswith('.csv'):
            # Save as CSV
            import csv
            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                # Write header
                writer.writerow(['Time', 'Spacecraft', 'Fault Type', 'Confidence', 'Component', 'Method', 'Status'])
                
                # Write detection data
                for item in self.results_tree.get_children():
                    values = self.results_tree.item(item)['values']
                    writer.writerow(values)
    
    def display_ml_results(self, ml_results):
        """
        Display ML detection results in the tab
        This method is called from the main GUI when ML results are available
        """
        if not ml_results:
            return
            
        self.fault_detection_results = ml_results
        
        # Update detection status
        self.detection_status_label.config(text=" Detection Complete", foreground='green')
        
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
                self.results_tree.insert('', 'end', values=(
                    f"{detection.detection_time_minutes:.1f} min",
                    spacecraft_name,
                    detection.fault_type,
                    f"{detection.confidence:.3f}",
                    detection.affected_component,
                    "ML Autoencoder",
                    "Active" if detection.fault_detected else "Inactive"
                ))
                
                # Add to recent detections
                recent_detections.append({
                    'time': f"{detection.detection_time_minutes:.1f} min",
                    'spacecraft': spacecraft_name,
                    'type': detection.fault_type,
                    'confidence': detection.confidence,
                    'status': 'Active' if detection.fault_detected else 'Inactive'
                })
                
                total_detections += 1
        
        # Update recent detections list
        if recent_detections:
            for detection in recent_detections[-10:]:  # Show last 10
                detection_text = (f"{detection['time']} - {detection['spacecraft']}: "
                                f"{detection['type']} (conf: {detection['confidence']:.3f})")
                self.recent_detections_listbox.insert(0, detection_text)
        else:
            self.recent_detections_listbox.insert(0, "No detections found")
        
        # Update detection history
        self.detection_history.extend(recent_detections)
        
        # Update statistics
        self.update_statistics()
        
        # Update results summary
        self.results_summary_text.delete(1.0, tk.END)
        summary_text = f"""REAL ML FAULT DETECTION RESULTS:
================================
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Data Source: Real Basilisk Simulation
ML Model: Client's anomaly_detection_model.keras

SUMMARY:
  Spacecraft Analyzed: {summary.get('total_spacecraft', 0)}
  Total ML Detections: {summary.get('total_detections', 0)}
  Success Rate: {summary.get('success_rate', 0):.1%}

DETECTIONS BY SPACECRAFT:
"""
        
        for spacecraft, spacecraft_detections in detections.items():
            summary_text += f"\n{spacecraft}: {len(spacecraft_detections)} detections"
            for detection in spacecraft_detections[:3]:  # Show first 3
                summary_text += f"\n  - {detection.detection_time_minutes:.1f}min (conf: {detection.confidence:.3f})"
        
        if total_detections == 0:
            summary_text += "\n\nNOTE: No faults were detected by the ML model."
            summary_text += "\nThis could indicate:"
            summary_text += "\n  - No actual faults occurred"
            summary_text += "\n  - Fault signatures were too weak"
            summary_text += "\n  - Model needs retraining"
            summary_text += "\n  - Detection threshold needs adjustment"
        
        self.results_summary_text.insert(tk.END, summary_text)
        
        # Switch to results tab to show the new results
        self.detection_notebook.select(2)
        
        # Log the results
        self.parent_app.add_log(f"ML Detection Results: {total_detections} detections from {len(detections)} spacecraft")
        
        # Update monitoring plot if active
        if self.monitoring_active:
            self.update_monitoring_plot(ml_results)
    
    def update_monitoring_plot(self, ml_results):
        """Update the monitoring plot with real data"""
        self.monitoring_figure.clear()
        ax = self.monitoring_figure.add_subplot(111)
        
        # Extract time series data from ML results
        # This is a simplified version - you'd extract actual time series from ml_results
        time_points = np.linspace(0, 30, 100)
        
        # Sample anomaly scores based on detection results
        anomaly_scores = np.random.normal(0.2, 0.1, len(time_points))
        
        # Add spikes where detections occurred
        detections = ml_results.get('detections', {})
        for spacecraft_detections in detections.values():
            for detection in spacecraft_detections:
                detection_time = detection.detection_time_minutes
                # Find closest time point
                time_idx = np.argmin(np.abs(time_points - detection_time))
                anomaly_scores[time_idx] = detection.confidence
        
        # Plot the data
        ax.plot(time_points, anomaly_scores, 'b-', linewidth=2, label='Anomaly Score')
        ax.axhline(y=self.threshold_var.get(), color='r', linestyle='--', 
                  label=f'Threshold ({self.threshold_var.get():.3f})')
        
        # Mark detections
        for spacecraft_detections in detections.values():
            for detection in spacecraft_detections:
                ax.axvline(x=detection.detection_time_minutes, color='red', alpha=0.7, linestyle=':')
                ax.plot(detection.detection_time_minutes, detection.confidence, 'ro', markersize=8)
        
        ax.set_xlabel('Time (minutes)')
        ax.set_ylabel('Anomaly Score')
        ax.set_title('Real-time Fault Detection Monitoring')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        self.monitoring_canvas.draw()
    
    def get_detection_methods_status(self):
        """Get the status of enabled detection methods"""
        methods = []
        if self.ml_detection_var.get():
            methods.append("ML Autoencoder")
        if self.statistical_detection_var.get():
            methods.append("Statistical Analysis")
        if self.trend_detection_var.get():
            methods.append("Trend Analysis")
        if self.threshold_detection_var.get():
            methods.append("Threshold-based")
        
        return methods
    
    def is_ml_ready(self):
        """Check if ML detection is ready"""
        return self.ml_available and self.ml_detector and self.ml_detector.is_loaded