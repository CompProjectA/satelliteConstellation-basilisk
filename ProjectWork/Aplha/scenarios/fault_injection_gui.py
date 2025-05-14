import tkinter as tk
from tkinter import ttk, messagebox
import threading
import time
import sys
import os
import numpy as np

# Add paths to ensure Basilisk modules can be found
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../models')
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../plotting')

from Basilisk.utilities import vizSupport, macros
import BSK_Plotting as BSK_plt
from Basilisk.architecture import messaging

class SimulationThread(threading.Thread):
    def __init__(self, gui):
        threading.Thread.__init__(self)
        self.gui = gui
        self.scenario = None
        self.running = False
        self.daemon = True  # Make thread close when main program exits
        self.vizInterface = None
        
# Replace the run() method in SimulationThread class with this fixed version
def run(self):
    self.running = True
    
    # Import scenario here to avoid circular imports
    try:
        from scenario_AddRWFault import scenario_AddRWFault
    except ImportError:
        self.gui.update_status("ERROR: Could not import scenario_AddRWFault!")
        return
        
    self.gui.update_status("Initializing simulation...")
    
    try:
        # Create scenario
        self.scenario = scenario_AddRWFault()
        
        # Set the mode request before initialization
        self.gui.update_status("Initializing simulation - this may take a few seconds...")
        self.scenario.modeRequest = "hillPoint"
        
        # IMPORTANT: Initialize the simulation first
        self.scenario.InitializeSimulation()

        # FIXED: Setup visualization AFTER initialization
        self.gui.update_status("Setting up visualization...")
        
        # Get the process name and task name directly from scenario
        processName = self.scenario.DynModels.processName
        taskName = self.scenario.DynModels.taskName
        
        # Make sure the spacecraft object is properly set up
        if hasattr(self.scenario.DynModels, 'scObject'):
            # Enable visualization with correct parameters
            self.vizInterface = vizSupport.enableUnityVisualization(
                self.scenario.TotalSim,        # Pass the simulation object
                self.scenario.DynModels.taskName,  # Pass the task name
                self.scenario.DynModels.scObject,  # Pass the spacecraft object
                saveFile=False,                # Don't save visualization data
                liveStream=True               # Enable live streaming
            )
            
            # Configure visualization settings explicitly
            if self.vizInterface and hasattr(self.vizInterface, 'settings'):
                self.vizInterface.settings.orbitLinesOn = 2
                self.vizInterface.settings.showSpacecraftLabels = True
            
            # Print connection information more clearly
            self.gui.update_status("Simulation running - Connect to Basilisk Sim in Vizard now")
            print("\n[SIMULATION] =====================================================")
            print("[SIMULATION] Connect to Basilisk Sim in Vizard now")
            print("[SIMULATION] Socket Address: tcp://localhost:5556")
            print("[SIMULATION] =====================================================\n")
            
            # Update the connection status periodically
            self.gui.root.after(1000, self.gui.update_connection_status)
        else:
            self.gui.update_status("ERROR: Spacecraft object not found!")
            return
        
        # Main simulation loop - keep running as long as the thread is active
        while self.running:
            self.scenario.TotalSim.SingleStepProcesses()
            current_time = self.scenario.TotalSim.CurrentNanos * macros.NANO2MIN
            self.gui.update_time(f"Simulation Time: {current_time:.2f} min")
            
            # Check if it's time for the automatic RW fault (if configured)
            if hasattr(self.scenario, 'oneTimeRWFaultFlag') and self.scenario.oneTimeRWFaultFlag:
                if self.scenario.TotalSim.CurrentNanos >= self.scenario.oneTimeFaultTime:
                    self.scenario.DynModels.AddRWFault('friction', 0.05, 1, self.scenario.TotalSim.CurrentNanos)
                    self.scenario.oneTimeRWFaultFlag = 0
                    self.gui.update_status(f"Automatic RW Fault triggered at simulation time " + 
                                           f"{self.scenario.TotalSim.CurrentNanos * macros.NANO2MIN:.2f} min")
            
            time.sleep(0.01)  # About 100Hz real-time
            
    except Exception as e:
        print(f"[ERROR] Simulation error: {str(e)}")
        import traceback
        traceback.print_exc()
        self.gui.update_status(f"Error: {str(e)[:50]}...")
        
    # Replace the update_connection_status method in SpacecraftFaultInjectionGUI class
def update_connection_status(self):
    if self.sim_thread and self.sim_thread.is_alive():
        if self.sim_thread.vizInterface:
            # Try to get the connection status
            try:
                # Check if the server is running
                if hasattr(self.sim_thread.vizInterface, '_server') and self.sim_thread.vizInterface._server:
                    # Check if there are any active clients
                    if hasattr(self.sim_thread.vizInterface._server, '_clients') and self.sim_thread.vizInterface._server._clients:
                        self.connection_status.config(text="Connection Status: Connected ✅", fg="green")
                        self.log_message("Vizard connected successfully!")
                    else:
                        self.connection_status.config(text="Connection Status: Server Running, Waiting for Vizard ⏳", fg="orange")
                        self.log_message("Server running, waiting for Vizard to connect...")
                        # Check again after 2 seconds
                        self.root.after(2000, self.update_connection_status)
                else:
                    self.connection_status.config(text="Connection Status: Server Ready, Open Vizard ⏳", fg="orange")
                    self.log_message("Server ready, waiting for Vizard...")
                    # Check again after 2 seconds
                    self.root.after(2000, self.update_connection_status)
            except Exception as e:
                self.log_message(f"Error checking connection: {str(e)}")
                self.connection_status.config(text="Connection Status: Error Checking Connection ⚠️", fg="red")
                # Check again after 3 seconds
                self.root.after(3000, self.update_connection_status)
        else:
            self.connection_status.config(text="Connection Status: Initializing Visualization ⏳", fg="orange")
            # Check again after 2 seconds
            self.root.after(2000, self.update_connection_status)
            
    # Replace the inject_fault method in SpacecraftFaultInjectionGUI class
def inject_fault(self):
    if not self.sim_thread or not self.sim_thread.is_alive() or not self.sim_thread.scenario:
        messagebox.showwarning("Warning", "Start simulation first!")
        return
        
    selected_fault = self.fault_type_var.get()
    if selected_fault == "" or selected_fault == "Select Fault Type":
        messagebox.showwarning("Warning", "Select a fault type!")
        return
    
    scenario = self.sim_thread.scenario
    severity = self.severity_var.get()
    
    try:
        # Get current simulation time for logging
        current_time = scenario.TotalSim.CurrentNanos * macros.NANO2MIN
        
        if selected_fault == "Reaction Wheel Friction Fault":
            # Target wheel ID (0-3 in BSK_Dynamics but documentation shows 1-4)
            # Use wheel_id-1 if BSK_Dynamics uses 0-based indexing
            wheel_id = 1  # You can make this configurable if needed
            scenario.DynModels.AddRWFault('friction', severity, wheel_id, scenario.TotalSim.CurrentNanos)
            self.log_message(f"Injected friction fault on RW {wheel_id} with severity {severity:.3f}")
            
        elif selected_fault == "Encoder Fault":
            scenario.DynModels.AddSensorNoiseFault('encoder', severity)
            self.log_message(f"Injected encoder fault with severity {severity:.3f}")
            
        elif selected_fault == "Battery Fault":
            scenario.DynModels.SimulateBatteryDrop()
            self.log_message("Injected battery voltage drop")
            
        elif selected_fault == "Thruster Fault":
            scenario.DynModels.AddThrusterFault(severity)
            self.log_message(f"Injected thruster fault with severity {severity:.3f}")
            
        elif selected_fault == "Power Loss Fault":
            scenario.DynModels.CutAllPower()
            self.log_message("Injected complete power loss")
        
        self.update_status(f"{selected_fault} Injected at {current_time:.2f} min ⚡")
        
    except Exception as e:
        error_msg = f"Error injecting fault: {str(e)}"
        self.log_message(error_msg)
        import traceback
        traceback.print_exc()
        messagebox.showerror("Fault Injection Error", error_msg)
        
    def stop(self):
        self.running = False
        if self.scenario:
            # Close any open plot windows
            import matplotlib.pyplot as plt
            plt.close('all')

class SpacecraftFaultInjectionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🛰️ Spacecraft Fault Injection GUI")
        self.root.geometry("550x650")  # Slightly taller to accommodate more information
        self.root.resizable(False, False)
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        self.sim_thread = None
        
        self.setup_ui()
        
    def setup_ui(self):
        title_label = tk.Label(self.root, text="Spacecraft Fault Injection GUI", font=("Helvetica", 18, "bold"))
        title_label.pack(pady=10)
        
        # Vizard connection frame
        vizard_frame = tk.LabelFrame(self.root, text="⚠️ IMPORTANT: Vizard Connection", font=("Helvetica", 12, "bold"), bg="lightyellow")
        vizard_frame.pack(fill="x", padx=20, pady=10)
        
        viz_label = tk.Label(vizard_frame, 
                            text="STEPS TO CONNECT:\n"
                            "1. Start simulation by clicking 'Start' below\n"
                            "2. Open Vizard application\n"
                            "3. In Vizard: Click 'Connect to Basilisk Sim'\n"
                            "4. Enter: tcp://localhost:5556\n"
                            "5. Select 'Live Display' mode\n"
                            "6. Click 'Start Visualization'",
                            font=("Helvetica", 11), bg="lightyellow", justify="left")
        viz_label.pack(pady=10)
        
        # Connection status indicator
        self.connection_status = tk.Label(vizard_frame, text="Connection Status: Not Connected", 
                                         font=("Helvetica", 10, "bold"), fg="red", bg="lightyellow")
        self.connection_status.pack(pady=5)
        
        # Simulation control frame
        sim_frame = tk.LabelFrame(self.root, text="Simulation Control", font=("Helvetica", 12))
        sim_frame.pack(fill="x", padx=20, pady=10)
        
        button_frame = tk.Frame(sim_frame)
        button_frame.pack(pady=10)
        
        self.start_button = tk.Button(button_frame, text="Start ▶️", font=("Helvetica", 12), 
                                     bg="lightgreen", width=12, command=self.start_simulation)
        self.start_button.grid(row=0, column=0, padx=5)
        
        self.stop_button = tk.Button(button_frame, text="Stop ⏹️", font=("Helvetica", 12),
                                    bg="lightcoral", width=12, command=self.stop_simulation)
        self.stop_button.grid(row=0, column=1, padx=5)
        self.stop_button.config(state=tk.DISABLED)
        
        self.show_plots_button = tk.Button(sim_frame, text="Show Telemetry Plots 📊", 
                                          font=("Helvetica", 12), width=25, command=self.show_plots)
        self.show_plots_button.pack(pady=10)
        self.show_plots_button.config(state=tk.DISABLED)
        
        # Fault injection frame
        fault_frame = tk.LabelFrame(self.root, text="Fault Injection", font=("Helvetica", 12))
        fault_frame.pack(fill="x", padx=20, pady=10)
        
        self.fault_type_var = tk.StringVar()
        fault_type_dropdown = ttk.Combobox(fault_frame, textvariable=self.fault_type_var, 
                                          state="readonly", font=("Helvetica", 12), width=30)
        fault_type_dropdown['values'] = ("Reaction Wheel Friction Fault", 
                                        "Encoder Fault", 
                                        "Battery Fault", 
                                        "Thruster Fault", 
                                        "Power Loss Fault")
        fault_type_dropdown.pack(pady=10)
        fault_type_dropdown.set("Select Fault Type")
        
        # Add severity slider for faults
        severity_frame = tk.Frame(fault_frame)
        severity_frame.pack(fill="x", pady=5)
        
        tk.Label(severity_frame, text="Fault Severity:", font=("Helvetica", 11)).pack(side=tk.LEFT, padx=5)
        
        self.severity_var = tk.DoubleVar(value=0.05)
        severity_slider = ttk.Scale(severity_frame, from_=0.01, to=0.2, 
                                   variable=self.severity_var, length=250)
        severity_slider.pack(side=tk.LEFT, padx=5)
        
        self.severity_label = tk.Label(severity_frame, text="0.05", font=("Helvetica", 11))
        self.severity_label.pack(side=tk.LEFT, padx=5)
        
        # Update severity label when slider changes
        def update_severity_label(*args):
            self.severity_label.config(text=f"{self.severity_var.get():.2f}")
        
        self.severity_var.trace("w", update_severity_label)
        
        self.inject_button = tk.Button(fault_frame, text="Inject Fault ⚡", font=("Helvetica", 13), 
                                      bg="orange", width=25, command=self.inject_fault)
        self.inject_button.pack(pady=10)
        self.inject_button.config(state=tk.DISABLED)
        
        # Status frame
        status_frame = tk.LabelFrame(self.root, text="Status", font=("Helvetica", 12))
        status_frame.pack(fill="x", padx=20, pady=10)
        
        self.status_label = tk.Label(status_frame, text="Status: Ready 🔄", font=("Helvetica", 12))
        self.status_label.pack(pady=5)
        
        self.time_label = tk.Label(status_frame, text="Simulation Time: 0.00 min", font=("Helvetica", 12))
        self.time_label.pack(pady=5)
        
        # Add debug log
        log_frame = tk.LabelFrame(self.root, text="Debug Log", font=("Helvetica", 10))
        log_frame.pack(fill="both", expand=True, padx=20, pady=10)
        
        self.log_text = tk.Text(log_frame, height=5, width=50, font=("Courier", 9))
        self.log_text.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Add a scrollbar
        scrollbar = ttk.Scrollbar(self.log_text, command=self.log_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_text.config(yscrollcommand=scrollbar.set)
        
        self.log_message("GUI initialized. Ready to start simulation.")
        
    def log_message(self, message):
        timestamp = time.strftime("%H:%M:%S", time.localtime())
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)  # Auto-scroll to the end
        
    def start_simulation(self):
        if self.sim_thread and self.sim_thread.is_alive():
            messagebox.showinfo("Info", "Simulation already running!")
            return
        
        # Start the simulation in a separate thread
        self.sim_thread = SimulationThread(self)
        self.sim_thread.start()
        
        # Update button states
        self.start_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)
        self.inject_button.config(state=tk.NORMAL)
        self.show_plots_button.config(state=tk.NORMAL)
        
        self.log_message("Starting simulation thread...")
        
        # Start a timer to check connection status periodically
        self.root.after(3000, self.update_connection_status)
        
    def update_connection_status(self):
        if self.sim_thread and self.sim_thread.is_alive():
            if self.sim_thread.vizInterface:
                # Check if vizInterface has connected property and it's True
                if hasattr(self.sim_thread.vizInterface, 'connected') and self.sim_thread.vizInterface.connected:
                    self.connection_status.config(text="Connection Status: Connected ✅", fg="green")
                    self.log_message("Vizard connected successfully!")
                else:
                    self.connection_status.config(text="Connection Status: Waiting for Vizard ⏳", fg="orange")
                    self.log_message("Still waiting for Vizard connection...")
                    # Check again after 2 seconds
                    self.root.after(2000, self.update_connection_status)
            else:
                self.connection_status.config(text="Connection Status: Ready, Open Vizard ⏳", fg="orange")
                # Check again after 2 seconds
                self.root.after(2000, self.update_connection_status)
        
    def stop_simulation(self):
        if self.sim_thread and self.sim_thread.is_alive():
            self.log_message("Stopping simulation...")
            self.sim_thread.stop()
            self.sim_thread.join(timeout=1.0)
            self.update_status("Simulation stopped")
            self.connection_status.config(text="Connection Status: Not Connected", fg="red")
        
        # Update button states
        self.start_button.config(state=tk.NORMAL)
        self.stop_button.config(state=tk.DISABLED)
        self.inject_button.config(state=tk.DISABLED)
        
    def inject_fault(self):
        if not self.sim_thread or not self.sim_thread.is_alive() or not self.sim_thread.scenario:
            messagebox.showwarning("Warning", "Start simulation first!")
            return
            
        selected_fault = self.fault_type_var.get()
        if selected_fault == "" or selected_fault == "Select Fault Type":
            messagebox.showwarning("Warning", "Select a fault type!")
            return
        
        scenario = self.sim_thread.scenario
        severity = self.severity_var.get()
        
        try:
            # Get current simulation time for logging
            current_time = scenario.TotalSim.CurrentNanos * macros.NANO2MIN
            
            if selected_fault == "Reaction Wheel Friction Fault":
                # Target wheel ID 0-3
                wheel_id = 1  # You can make this configurable if needed
                scenario.DynModels.AddRWFault('friction', severity, wheel_id, scenario.TotalSim.CurrentNanos)
                self.log_message(f"Injected friction fault on RW {wheel_id} with severity {severity:.3f}")
                
            elif selected_fault == "Encoder Fault":
                scenario.DynModels.AddSensorNoiseFault('encoder', severity)
                self.log_message(f"Injected encoder fault with severity {severity:.3f}")
                
            elif selected_fault == "Battery Fault":
                scenario.DynModels.SimulateBatteryDrop()
                self.log_message("Injected battery voltage drop")
                
            elif selected_fault == "Thruster Fault":
                scenario.DynModels.AddThrusterFault(severity)
                self.log_message(f"Injected thruster fault with severity {severity:.3f}")
                
            elif selected_fault == "Power Loss Fault":
                scenario.DynModels.CutAllPower()
                self.log_message("Injected complete power loss")
            
            self.update_status(f"{selected_fault} Injected at {current_time:.2f} min ⚡")
            
        except Exception as e:
            error_msg = f"Error injecting fault: {str(e)}"
            self.log_message(error_msg)
            messagebox.showerror("Fault Injection Error", error_msg)
        
    def show_plots(self):
        if not self.sim_thread or not self.sim_thread.is_alive() or not self.sim_thread.scenario:
            messagebox.showwarning("Warning", "Start simulation first!")
            return
            
        # Create a new thread for displaying plots to avoid freezing GUI
        def display_plots():
            try:
                self.log_message("Generating telemetry plots...")
                self.sim_thread.scenario.pull_outputs(True)
                self.log_message("Plots generated successfully")
            except Exception as e:
                error_msg = f"Error displaying plots: {str(e)}"
                self.log_message(error_msg)
                messagebox.showerror("Plot Error", error_msg)
                
        plot_thread = threading.Thread(target=display_plots)
        plot_thread.daemon = True
        plot_thread.start()
        
    def update_status(self, message):
        self.status_label.config(text=f"Status: {message}")
        self.log_message(message)
        
    def update_time(self, time_text):
        self.time_label.config(text=time_text)
        
    def on_closing(self):
        if self.sim_thread and self.sim_thread.is_alive():
            self.log_message("Shutting down simulation...")
            self.sim_thread.stop()
            self.sim_thread.join(timeout=1.0)
        self.root.destroy()

if __name__ == "__main__":
    # Make sure we have the required modules
    try:
        import Basilisk
        from Basilisk.utilities import vizSupport
        if not vizSupport.vizFound:
            print("WARNING: Vizard visualization support not found!")
    except ImportError:
        print("ERROR: Basilisk modules not found. Check your Python environment.")
        sys.exit(1)
        
    root = tk.Tk()
    app = SpacecraftFaultInjectionGUI(root)
    root.mainloop()