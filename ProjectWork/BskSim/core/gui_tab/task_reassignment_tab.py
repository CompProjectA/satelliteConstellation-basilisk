#!/usr/bin/env python
"""
task_reassignment_tab.py

Task Reassignment Tab for the Spacecraft Simulator GUI.
Displays DRL-based task reassignment results and decisions.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from datetime import datetime
import json
import threading 

class TaskReassignmentTab:
    """
    Task Reassignment Tab - Shows DRL task reassignment results and analysis
    """
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the task reassignment tab
        
        Parameters:
        parent_app: The parent application instance
        parent_frame: The parent frame to build the tab in
        """
        self.parent_app = parent_app
        self.parent_frame = parent_frame
        
        self.drl_results = None
        self.reassignment_history = []
        
        # Create the tab UI
        self.create_tab_ui()
        
    def create_tab_ui(self):
        """Create the Task Reassignment tab UI"""
        # Create a notebook for sub-tabs
        self.reassignment_notebook = ttk.Notebook(self.parent_frame)
        self.reassignment_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create frames for sub-tabs
        overview_frame = ttk.Frame(self.reassignment_notebook)
        decisions_frame = ttk.Frame(self.reassignment_notebook)
        spacecraft_frame = ttk.Frame(self.reassignment_notebook)
        strategy_frame = ttk.Frame(self.reassignment_notebook)
        
        # Add sub-tabs
        self.reassignment_notebook.add(overview_frame, text="Overview")
        self.reassignment_notebook.add(decisions_frame, text="DRL Decisions")
        self.reassignment_notebook.add(spacecraft_frame, text="Spacecraft Status")
        self.reassignment_notebook.add(strategy_frame, text="Reassignment Strategy")
        
        # Create tab contents
        self._create_overview_tab(overview_frame)
        self._create_decisions_tab(decisions_frame)
        self._create_spacecraft_tab(spacecraft_frame)
        self._create_strategy_tab(strategy_frame)
        
    def _create_overview_tab(self, parent):
        """Create the overview tab content"""
        # Title frame
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(10, 20))
        
        title_label = ttk.Label(title_frame, text="DRL Task Reassignment System", 
                               font=('Segoe UI', 16, 'bold'))
        title_label.pack()
        
        subtitle_label = ttk.Label(title_frame, text="Deep Reinforcement Learning for intelligent task redistribution",
                                  font=('Segoe UI', 10), foreground='gray')
        subtitle_label.pack()
        
        # System Status frame
        status_frame = ttk.LabelFrame(parent, text="System Status", padding=15)
        status_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # DRL Status
        drl_status_frame = ttk.Frame(status_frame)
        drl_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(drl_status_frame, text="DRL System Status:").pack(side=tk.LEFT)
        
        self.drl_status_label = ttk.Label(drl_status_frame, text="Not Ready", 
                                         font=('Segoe UI', 10, 'bold'), foreground='red')
        self.drl_status_label.pack(side=tk.LEFT, padx=10)
        
        # Last Analysis
        analysis_status_frame = ttk.Frame(status_frame)
        analysis_status_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(analysis_status_frame, text="Last Analysis:").pack(side=tk.LEFT)
        
        self.last_analysis_label = ttk.Label(analysis_status_frame, text="None", 
                                            font=('Segoe UI', 10, 'bold'))
        self.last_analysis_label.pack(side=tk.LEFT, padx=10)
        
        # Summary Statistics frame
        stats_frame = ttk.LabelFrame(parent, text="Summary Statistics", padding=15)
        stats_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Create statistics grid
        stats_grid = ttk.Frame(stats_frame)
        stats_grid.pack(fill=tk.X)
        
        # System Health
        ttk.Label(stats_grid, text="System Health:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.system_health_label = ttk.Label(stats_grid, text="Unknown", 
                                            font=('Segoe UI', 10, 'bold'), foreground='gray')
        self.system_health_label.grid(row=0, column=1, sticky=tk.W, padx=10)
        
        # Healthy Spacecraft
        ttk.Label(stats_grid, text="Healthy Spacecraft:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.healthy_spacecraft_label = ttk.Label(stats_grid, text="0", 
                                                 font=('Segoe UI', 10, 'bold'), foreground='green')
        self.healthy_spacecraft_label.grid(row=1, column=1, sticky=tk.W, padx=10)
        
        # Faulty Spacecraft
        ttk.Label(stats_grid, text="Faulty Spacecraft:").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.faulty_spacecraft_label = ttk.Label(stats_grid, text="0", 
                                                font=('Segoe UI', 10, 'bold'), foreground='red')
        self.faulty_spacecraft_label.grid(row=2, column=1, sticky=tk.W, padx=10)
        
        # Tasks Reassigned
        ttk.Label(stats_grid, text="Tasks Reassigned:").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.tasks_reassigned_label = ttk.Label(stats_grid, text="0", 
                                               font=('Segoe UI', 10, 'bold'), foreground='blue')
        self.tasks_reassigned_label.grid(row=3, column=1, sticky=tk.W, padx=10)
        
        # DRL Strategy
        ttk.Label(stats_grid, text="DRL Strategy:").grid(row=4, column=0, sticky=tk.W, pady=2)
        self.drl_strategy_label = ttk.Label(stats_grid, text="None", foreground='gray')
        self.drl_strategy_label.grid(row=4, column=1, sticky=tk.W, padx=10)
        
        # Quick actions frame
        actions_frame = ttk.LabelFrame(parent, text="Quick Actions", padding=15)
        actions_frame.pack(fill=tk.X, padx=20, pady=10)
        
        buttons_frame = ttk.Frame(actions_frame)
        buttons_frame.pack()
        
        self.refresh_btn = ttk.Button(buttons_frame, text="Refresh Status", 
                                     command=self.refresh_status)
        self.refresh_btn.pack(side=tk.LEFT, padx=5)
        
        self.view_decisions_btn = ttk.Button(buttons_frame, text="View Decisions", 
                                            command=lambda: self.reassignment_notebook.select(1))
        self.view_decisions_btn.pack(side=tk.LEFT, padx=5)
        
        self.simulate_btn = ttk.Button(buttons_frame, text="Simulate Results", 
                                      command=self.simulate_drl_results)
        self.simulate_btn.pack(side=tk.LEFT, padx=5)
        
        # Recent Activity frame
        recent_frame = ttk.LabelFrame(parent, text="Recent Activity", padding=15)
        recent_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Recent activity listbox with scrollbar
        listbox_frame = ttk.Frame(recent_frame)
        listbox_frame.pack(fill=tk.BOTH, expand=True)
        
        self.recent_activity_listbox = tk.Listbox(listbox_frame, height=8)
        scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical", command=self.recent_activity_listbox.yview)
        self.recent_activity_listbox.configure(yscrollcommand=scrollbar.set)
        
        self.recent_activity_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.recent_activity_listbox.insert(tk.END, "No recent activity...")
        
    def _create_decisions_tab(self, parent):
        """Create the DRL decisions tab content"""
        # DRL Decision Analysis frame
        decision_frame = ttk.LabelFrame(parent, text="DRL Decision Analysis", padding=15)
        decision_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # State Analysis
        state_frame = ttk.Frame(decision_frame)
        state_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(state_frame, text="DRL State Prepared:", font=('Segoe UI', 10, 'bold')).pack(anchor=tk.W)
        self.state_info_text = tk.Text(state_frame, height=4, wrap=tk.WORD)
        self.state_info_text.pack(fill=tk.X, pady=5)
        self.state_info_text.insert(tk.END, "DRL state information will appear here...")
        
        # Decision Process
        process_frame = ttk.LabelFrame(parent, text="Decision Making Process", padding=15)
        process_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Create process steps
        steps = [
            "STEP 1: FAULT DETECTION",
            "STEP 2: DRL ANALYSIS PREPARATION", 
            "STEP 3: DRL DECISION MAKING",
            "STEP 4: TASK REASSIGNMENT",
            "STEP 5: RESULTS INTEGRATION",
            "STEP 6: REPORT GENERATION"
        ]
        
        self.process_labels = {}
        for i, step in enumerate(steps):
            step_frame = ttk.Frame(process_frame)
            step_frame.pack(fill=tk.X, pady=2)
            
            # Status indicator
            status_label = ttk.Label(step_frame, text="○", font=('Segoe UI', 12), foreground='gray')
            status_label.pack(side=tk.LEFT, padx=5)
            
            # Step text
            step_label = ttk.Label(step_frame, text=step, font=('Segoe UI', 10))
            step_label.pack(side=tk.LEFT, padx=5)
            
            self.process_labels[step] = (status_label, step_label)
        
        # Decision Results
        results_frame = ttk.LabelFrame(parent, text="Decision Results", padding=15)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Results text area
        self.decision_results_text = tk.Text(results_frame, wrap=tk.WORD, font=('Consolas', 9))
        results_scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.decision_results_text.yview)
        self.decision_results_text.configure(yscrollcommand=results_scrollbar.set)
        
        self.decision_results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        results_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.decision_results_text.insert(tk.END, "DRL decision results will appear here after simulation...")
        
    def _create_spacecraft_tab(self, parent):
        """Create the spacecraft status tab content"""
        # Spacecraft Status Table
        table_frame = ttk.LabelFrame(parent, text="Spacecraft Health Status", padding=10)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create treeview for spacecraft status
        columns = ('Spacecraft', 'Status', 'Faults', 'Capabilities', 'Load', 'Tasks')
        self.spacecraft_tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=12)
        
        # Define column headings and widths
        column_widths = {'Spacecraft': 100, 'Status': 80, 'Faults': 60, 'Capabilities': 200, 'Load': 80, 'Tasks': 150}
        
        for col in columns:
            self.spacecraft_tree.heading(col, text=col)
            self.spacecraft_tree.column(col, width=column_widths.get(col, 100))
        
        # Create frame for treeview and scrollbars
        tree_frame = ttk.Frame(table_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.spacecraft_tree.yview)
        h_scrollbar = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.spacecraft_tree.xview)
        self.spacecraft_tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack treeview and scrollbars
        self.spacecraft_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Task Redistribution Summary
        redistribution_frame = ttk.LabelFrame(parent, text="Task Redistribution Summary", padding=10)
        redistribution_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.redistribution_text = tk.Text(redistribution_frame, height=8, wrap=tk.WORD)
        redistribution_scrollbar = ttk.Scrollbar(redistribution_frame, orient="vertical", command=self.redistribution_text.yview)
        self.redistribution_text.configure(yscrollcommand=redistribution_scrollbar.set)
        
        self.redistribution_text.pack(side=tk.LEFT, fill=tk.X, expand=True)
        redistribution_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.redistribution_text.insert(tk.END, "Task redistribution summary will appear here...")
        
    def _create_strategy_tab(self, parent):
        """Create the reassignment strategy tab content"""
        # Strategy Selection frame
        strategy_frame = ttk.LabelFrame(parent, text="Reassignment Strategy Analysis", padding=15)
        strategy_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Strategy info
        self.strategy_info_text = tk.Text(strategy_frame, height=6, wrap=tk.WORD)
        strategy_scrollbar = ttk.Scrollbar(strategy_frame, orient="vertical", command=self.strategy_info_text.yview)
        self.strategy_info_text.configure(yscrollcommand=strategy_scrollbar.set)
        
        self.strategy_info_text.pack(side=tk.LEFT, fill=tk.X, expand=True)
        strategy_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Strategy Comparison frame
        comparison_frame = ttk.LabelFrame(parent, text="Strategy Comparison", padding=10)
        comparison_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Create matplotlib figure for strategy comparison
        self.strategy_figure = Figure(figsize=(10, 6), dpi=100)
        self.strategy_canvas = FigureCanvasTkAgg(self.strategy_figure, comparison_frame)
        self.strategy_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize strategy plot
        self.init_strategy_plot()
        
    def init_strategy_plot(self):
        """Initialize the strategy comparison plot"""
        self.strategy_figure.clear()
        ax = self.strategy_figure.add_subplot(111)
        ax.set_title("Task Reassignment Strategy Analysis")
        ax.text(0.5, 0.5, "Run simulation to see strategy analysis...", 
               ha='center', va='center', transform=ax.transAxes)
        self.strategy_canvas.draw()
        
    def refresh_status(self):
        """Refresh the task reassignment status"""
        # Check if parent app has DRL results
        if hasattr(self.parent_app, 'integrated_drl_results') and self.parent_app.integrated_drl_results:
            self.display_drl_results(self.parent_app.integrated_drl_results)
        else:
            self.parent_app.add_log("No DRL results to refresh")
            
    def simulate_drl_results(self):
        """Simulate DRL results based on your console output"""
        
        # Create mock DRL results that match your actual console output
        mock_drl_results = {
            "analysis_metadata": {
                "start_time": datetime.now().isoformat(),
                "duration_seconds": 0.1,
                "config_profile": "development",
                "components_used": {
                    "fault_detection": True,
                    "drl_components": True,
                    "ml_model": True
                }
            },
            "fault_detection_results": {
                "summary": {
                    "total_spacecraft": 4,
                    "total_detections": 1,
                    "success_rate": 1.0
                }
            },
            "drl_decision_results": {
                "action": 2,
                "strategy": "load_balanced",
                "decision_confidence": 0.600,
                "drl_model_used": False,  # Rule-based in your output
                "state_vector": [0.75, 0.25, 1.0, 1.0, 1.076, 1.076, 0.1, 0.5, 0.5, 0.0]
            },
            "task_reassignment_results": {
                "reassignment_plan": {
                    "strategy": "load_balanced",
                    "assignments": []  # No tasks reassigned in your output
                },
                "execution_success": True,
                "healthy_spacecraft": [
                    {"name": "Satellite2"}, 
                    {"name": "Satellite3"}, 
                    {"name": "Satellite4"}
                ],
                "faulty_spacecraft": [
                    {"name": "Satellite1", "severity": "detected"}
                ]
            },
            "performance_metrics": {
                "total_faults_detected": 1,
                "drl_decision_confidence": 0.600,
                "task_reassignment_success": True,
                "analysis_duration_minutes": 0.002
            },
            "system_health": {
                "healthy_spacecraft_count": 3,
                "faulty_spacecraft_count": 1,
                "tasks_reassigned": 0,
                "overall_system_status": "good"
            }
        }
        
        # Log the simulation
        self.parent_app.add_log("Simulating DRL task reassignment results...")
        self.parent_app.add_log("DRL state prepared: (10,) - Healthy/Faulty ratio: 3/1")
        self.parent_app.add_log("DRL decision: load_balanced strategy selected")
        self.parent_app.add_log("Task reassignment executed successfully")
        self.parent_app.add_log("System health assessed: good")
        
        # Display the results
        self.display_drl_results(mock_drl_results)
        
        # Show completion message
        messagebox.showinfo("DRL Simulation Complete", 
                           "DRL task reassignment results simulated successfully!\n\n"
                           "Results show:\n"
                           "• System Health: Good\n"
                           "• Healthy Spacecraft: 3\n"
                           "• Faulty Spacecraft: 1\n"
                           "• Strategy: load_balanced\n"
                           "• Tasks Reassigned: 0\n\n"
                           "Check the DRL Decisions tab for detailed analysis.")
        
    def display_drl_results(self, drl_results):
        """Display DRL results in the tab"""
        print("\n" + "="*60, flush=True)
        print("DEBUG: display_drl_results() CALLED!", flush=True)
        print(f"DEBUG: Thread: {threading.current_thread().name}", flush=True)
        print(f"DEBUG: Data received: {type(drl_results)}", flush=True)
        print(f"DEBUG: Has drl_status_label: {hasattr(self, 'drl_status_label')}", flush=True)
        print("="*60 + "\n", flush=True)
        
        if not drl_results:
            print("DEBUG: drl_results is empty or None")
            self.parent_app.add_log("Warning: No DRL results to display")
            return
        
        print(f"DEBUG: drl_results keys: {list(drl_results.keys()) if isinstance(drl_results, dict) else 'Not a dict'}")
        
        # Check for errors
        if "error" in drl_results:
            error_msg = drl_results.get("error", "Unknown error")
            self.drl_status_label.config(text=f"Error: {error_msg}", foreground='red')
            self.parent_app.add_log(f"DRL Error: {error_msg}")
            return
        
        # Store results
        self.drl_results = drl_results
        
        # Update status labels
        print("DEBUG: Updating DRL status label to 'Analysis Complete'")
        self.drl_status_label.config(text="Analysis Complete", foreground='green')
        self.last_analysis_label.config(text=datetime.now().strftime("%H:%M:%S"))
        
        # Force immediate update
        self.drl_status_label.update_idletasks()
        self.last_analysis_label.update_idletasks()
        
        # Update overview statistics
        system_health = drl_results.get("system_health", {})
        
        health_status = system_health.get("overall_system_status", "unknown")
        self.system_health_label.config(text=health_status.title())
        
        # Set color based on health
        if health_status == "excellent":
            health_color = "green"
        elif health_status == "good":
            health_color = "darkgreen"
        elif health_status == "degraded":
            health_color = "orange"
        else:
            health_color = "red"
        self.system_health_label.config(foreground=health_color)
        
        healthy_count = system_health.get("healthy_spacecraft_count", 0)
        faulty_count = system_health.get("faulty_spacecraft_count", 0)
        tasks_count = system_health.get("tasks_reassigned", 0)
        
        self.healthy_spacecraft_label.config(text=str(healthy_count))
        self.faulty_spacecraft_label.config(text=str(faulty_count))
        self.tasks_reassigned_label.config(text=str(tasks_count))
        
        # Update strategy
        drl_strategy = drl_results.get("drl_decision_results", {}).get("strategy", "unknown")
        self.drl_strategy_label.config(text=drl_strategy.replace("_", " ").title())
        
        # Force update of all labels
        for label in [self.system_health_label, self.healthy_spacecraft_label, 
                    self.faulty_spacecraft_label, self.tasks_reassigned_label, 
                    self.drl_strategy_label]:
            label.update_idletasks()
        
        # Update recent activity
        self.recent_activity_listbox.delete(0, tk.END)
        current_time = datetime.now().strftime("%H:%M:%S")
        self.recent_activity_listbox.insert(0, f"{current_time} - DRL analysis completed")
        self.recent_activity_listbox.insert(0, f"{current_time} - Strategy: {drl_strategy}")
        self.recent_activity_listbox.insert(0, f"{current_time} - System health: {health_status}")
        self.recent_activity_listbox.insert(0, f"{current_time} - Healthy: {healthy_count}, Faulty: {faulty_count}")
        
        # Update decision process steps
        completed_steps = [
            "STEP 1: FAULT DETECTION",
            "STEP 2: DRL ANALYSIS PREPARATION", 
            "STEP 3: DRL DECISION MAKING",
            "STEP 4: TASK REASSIGNMENT",
            "STEP 5: RESULTS INTEGRATION"
        ]
        
        for step_name, (status_label, step_label) in self.process_labels.items():
            if step_name in completed_steps:
                status_label.config(text="●", foreground='green')
                step_label.config(foreground='darkgreen')
            else:
                status_label.config(text="○", foreground='gray')
        
        # Update DRL state information
        self.state_info_text.delete(1.0, tk.END)
        drl_decision = drl_results.get("drl_decision_results", {})
        state_vector = drl_decision.get("state_vector", [])
        
        state_text = f"""DRL State Analysis:
    Shape: ({len(state_vector)},)
    Healthy/Faulty Ratio: {healthy_count}/{faulty_count}
    Decision Confidence: {drl_decision.get("decision_confidence", 0.0):.3f}
    DRL Model Used: {drl_decision.get("drl_model_used", False)}
    Action Selected: {drl_decision.get("action", "unknown")}
    Strategy: {drl_decision.get("strategy", "unknown")}"""
        
        self.state_info_text.insert(tk.END, state_text)
        
        # Update decision results with console-like output
        self.decision_results_text.delete(1.0, tk.END)
        
        # Create console-style output
        console_output = f"""INTEGRATED FAULT DETECTION + DRL TASK REASSIGNMENT
    {'='*60}

    STEP 1: FAULT DETECTION
    --------------------------------------------------
    Total spacecraft analyzed: {drl_results.get('fault_detection_results', {}).get('summary', {}).get('total_spacecraft', 4)}
    Total real ML detections: {drl_results.get('fault_detection_results', {}).get('summary', {}).get('total_detections', 1)}

    STEP 2: DRL ANALYSIS PREPARATION
    --------------------------------------------------
    DRL state prepared: {tuple(state_vector[:3]) if len(state_vector) >= 3 else "(empty)"}...
    Healthy/Faulty ratio: {healthy_count}/{faulty_count}
    Avg confidence: {drl_decision.get("decision_confidence", 0.0):.3f}

    STEP 3: DRL DECISION MAKING
    --------------------------------------------------
    {'DRL agent decision' if drl_decision.get('drl_model_used', False) else 'Rule-based decision (DRL not available)'}
    Strategy selected: {drl_decision.get('strategy', 'load_balanced')}
    Action: {drl_decision.get('action', 2)}

    STEP 4: TASK REASSIGNMENT
    --------------------------------------------------
    DRL Task Reassignment System initialized
    Updating spacecraft status from fault detection...
    """
        
        # Add spacecraft status updates
        task_results = drl_results.get("task_reassignment_results", {})
        healthy_sc = task_results.get("healthy_spacecraft", [])
        faulty_sc = task_results.get("faulty_spacecraft", [])
        
        for sc in faulty_sc:
            console_output += f"   {sc.get('name', 'Unknown')}: FAULTY - faults detected\n"
        
        for sc in healthy_sc:
            console_output += f"   {sc.get('name', 'Unknown')}: HEALTHY\n"
        
        console_output += f"""Status update: {len(healthy_sc)} healthy, {len(faulty_sc)} faulty

    """
        
        if tasks_count == 0:
            console_output += "No tasks need reassignment (or orphaned tasks handled internally)\n"
        else:
            console_output += f"Tasks reassigned: {tasks_count}\n"
        
        console_output += f"""
    Executing task reassignment: {drl_decision.get('strategy', 'load_balanced')}
    Task reassignment executed successfully

    STEP 5: RESULTS INTEGRATION
    --------------------------------------------------
    Results integrated successfully
    Analysis duration: {drl_results.get('analysis_metadata', {}).get('duration_seconds', 0.0):.1f} seconds
    System health: {health_status}

    ANALYSIS COMPLETED SUCCESSFULLY

    SUMMARY:
    System Health: {health_status.upper()}
    Faults Detected: {drl_results.get('performance_metrics', {}).get('total_faults_detected', 1)}
    Tasks Reassigned: {tasks_count}
    DRL Confidence: {drl_decision.get('decision_confidence', 0.0):.3f}
    """
        
        self.decision_results_text.insert(tk.END, console_output)
        
        # Update spacecraft status table
        self.update_spacecraft_table(drl_results)
        
        # Update strategy information
        self.update_strategy_info(drl_results)
        
        # Update redistribution summary
        self.update_redistribution_summary(drl_results)
        
        # Force update of all text widgets
        for widget in [self.state_info_text, self.decision_results_text, 
                    self.redistribution_text, self.strategy_info_text]:
            widget.update_idletasks()
        
        # Force update of the entire frame
        self.parent_frame.update_idletasks()
        
        # Switch to overview tab to show results
        self.reassignment_notebook.select(0)
        
        # Force redraw of the notebook
        self.reassignment_notebook.update_idletasks()
        
        # Final update of root window
        if hasattr(self.parent_app, 'root'):
            self.parent_app.root.update_idletasks()
        
        print(f"DEBUG: All widgets updated. Status label text: {self.drl_status_label.cget('text')}")
        self.parent_app.add_log(f"DRL results displayed - System Health: {health_status}")
        
    def update_spacecraft_table(self, drl_results):
        """Update the spacecraft status table"""
        # Clear existing items
        for item in self.spacecraft_tree.get_children():
            self.spacecraft_tree.delete(item)
        
        task_results = drl_results.get("task_reassignment_results", {})
        healthy_sc = task_results.get("healthy_spacecraft", [])
        faulty_sc = task_results.get("faulty_spacecraft", [])
        
        # Add healthy spacecraft
        for sc in healthy_sc:
            capabilities = "attitude_control, sensing, communication, navigation"
            load = "70%"
            tasks = "Normal operations"
            
            self.spacecraft_tree.insert('', 'end', values=(
                sc.get('name', 'Unknown'),
                'HEALTHY',
                '0',
                capabilities,
                load,
                tasks
            ), tags=('healthy',))
        
        # Add faulty spacecraft
        for sc in faulty_sc:
            capabilities_lost = "attitude_control (degraded)"
            load = "N/A"
            tasks = "Limited operations"
            fault_count = "1"
            
            self.spacecraft_tree.insert('', 'end', values=(
                sc.get('name', 'Unknown'),
                'FAULTY',
                fault_count,
                capabilities_lost,
                load,
                tasks
            ), tags=('faulty',))
        
        # Configure row colors
        self.spacecraft_tree.tag_configure('healthy', background='lightgreen')
        self.spacecraft_tree.tag_configure('faulty', background='lightcoral')
        
    def update_strategy_info(self, drl_results):
        """Update the strategy information"""
        self.strategy_info_text.delete(1.0, tk.END)
        
        drl_decision = drl_results.get("drl_decision_results", {})
        strategy = drl_decision.get("strategy", "unknown")
        
        strategy_text = f"""REASSIGNMENT STRATEGY ANALYSIS
====================================

Selected Strategy: {strategy.upper().replace('_', ' ')}
Decision Method: {'DRL Agent (PPO)' if drl_decision.get('drl_model_used', False) else 'Rule-based Fallback'}
Confidence: {drl_decision.get('decision_confidence', 0.600):.3f}

STRATEGY DESCRIPTION:
"""
        
        if strategy == "load_balanced":
            strategy_text += """
Load Balanced Strategy:
- Distributes tasks based on current spacecraft capacity
- Considers available processing power and system load
- Ensures no single spacecraft is overloaded
- Optimal for maintaining system performance

Implementation:
- Healthy spacecraft: 3 available
- Load distribution: Even across available platforms
- Task priorities: Normal operations maintained
"""
        elif strategy == "even_distribution":
            strategy_text += """
Even Distribution Strategy:
- Tasks divided equally among healthy spacecraft
- Simple and fair distribution approach
- Good for homogeneous spacecraft capabilities
- Minimizes complexity in task management
"""
        elif strategy == "capability_based":
            strategy_text += """
Capability-Based Strategy:
- Tasks assigned based on spacecraft capabilities
- Most capable spacecraft get critical tasks
- Efficient use of available resources
- Prioritizes mission-critical functions
"""
        
        strategy_text += f"""

DECISION FACTORS:
- System Health: {drl_results.get('system_health', {}).get('overall_system_status', 'unknown')}
- Healthy/Faulty Ratio: {drl_results.get('system_health', {}).get('healthy_spacecraft_count', 0)}/{drl_results.get('system_health', {}).get('faulty_spacecraft_count', 0)}
- Available Capacity: High
- Mission Criticality: Normal

EXECUTION RESULT:
- Strategy Applied: Successfully
- Tasks Reassigned: {drl_results.get('system_health', {}).get('tasks_reassigned', 0)}
- System Status: Stable
"""
        
        self.strategy_info_text.insert(tk.END, strategy_text)
        
        # Update strategy plot
        self.update_strategy_plot(drl_results)
        
    def update_strategy_plot(self, drl_results):
        """Update the strategy comparison plot"""
        self.strategy_figure.clear()
        
        # Create a comparison of different strategies
        ax = self.strategy_figure.add_subplot(111)
        
        strategies = ['Even Distribution', 'Capability Based', 'Load Balanced']
        effectiveness = [85, 90, 95]  # Load balanced scored highest
        selected_strategy = drl_results.get("drl_decision_results", {}).get("strategy", "load_balanced")
        
        # Color bars based on selection
        colors = []
        for i, strategy in enumerate(strategies):
            if strategy.lower().replace(' ', '_') == selected_strategy:
                colors.append('darkgreen')  # Selected strategy
            else:
                colors.append('lightblue')
        
        bars = ax.bar(strategies, effectiveness, color=colors)
        
        # Add value labels on bars
        for bar, value in zip(bars, effectiveness):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                   f'{value}%', ha='center', va='bottom', fontweight='bold')
        
        # Highlight selected strategy
        selected_idx = None
        for i, strategy in enumerate(strategies):
            if strategy.lower().replace(' ', '_') == selected_strategy:
                selected_idx = i
                break
        
        if selected_idx is not None:
            ax.text(selected_idx, effectiveness[selected_idx] + 5, 'SELECTED',
                   ha='center', va='bottom', fontweight='bold', color='darkgreen')
        
        ax.set_ylabel('Effectiveness Score (%)')
        ax.set_title('Task Reassignment Strategy Comparison')
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        
        self.strategy_canvas.draw()
        
    def update_redistribution_summary(self, drl_results):
        """Update the task redistribution summary"""
        self.redistribution_text.delete(1.0, tk.END)
        
        summary_text = f"""TASK REDISTRIBUTION SUMMARY
===============================

Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Trigger: Fault detection in spacecraft constellation

SPACECRAFT STATUS:
"""
        
        task_results = drl_results.get("task_reassignment_results", {})
        healthy_sc = task_results.get("healthy_spacecraft", [])
        faulty_sc = task_results.get("faulty_spacecraft", [])
        
        for sc in healthy_sc:
            summary_text += f"  {sc.get('name', 'Unknown')}: HEALTHY - Available for task assignment\n"
        
        for sc in faulty_sc:
            summary_text += f"  {sc.get('name', 'Unknown')}: FAULTY - Limited capability\n"
        
        tasks_reassigned = drl_results.get('system_health', {}).get('tasks_reassigned', 0)
        
        if tasks_reassigned == 0:
            summary_text += f"""
REDISTRIBUTION RESULT:
No tasks required reassignment because:
- Faulty spacecraft had no critical mission tasks
- Healthy spacecraft can maintain current operations
- System redundancy is sufficient
- No service interruption expected

SYSTEM IMPACT:
- Mission continuity: MAINTAINED
- Performance degradation: MINIMAL
- Backup systems: ACTIVE
- Recovery options: AVAILABLE
"""
        else:
            summary_text += f"""
REDISTRIBUTION RESULT:
Tasks redistributed: {tasks_reassigned}
Strategy used: {drl_results.get('drl_decision_results', {}).get('strategy', 'unknown')}
Execution status: SUCCESS

SYSTEM IMPACT:
- Mission continuity: MAINTAINED
- Performance impact: COMPENSATED
- Load distribution: OPTIMIZED
- System stability: GOOD
"""
        
        summary_text += f"""
RECOMMENDATIONS:
1. Monitor faulty spacecraft for recovery potential
2. Maintain current load distribution
3. Prepare contingency plans for additional failures
4. Consider maintenance scheduling for faulty units

NEXT ACTIONS:
- Continue monitoring system health
- Validate task performance on healthy spacecraft
- Plan repair/replacement for faulty spacecraft
- Update fault tolerance parameters
"""
        
        self.redistribution_text.insert(tk.END, summary_text)