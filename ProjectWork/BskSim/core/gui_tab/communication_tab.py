#!/usr/bin/env python
"""
communication_tab.py - Enhanced Communication Tab with full cluster integration
Handles intra-cluster (leader-child) and inter-cluster (leader-leader) communication
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from datetime import datetime
import threading

class CommunicationTab:
    def __init__(self, parent, parent_frame):
        self.parent = parent
        self.parent_frame = parent_frame
        self.message_history = []
        self.access_matrix = {}
        self.communication_active = False
        self.sim_time = 0.0
        
        self._create_widgets()
        
    def _create_widgets(self):
        """Create communication tab widgets"""
        # Main container with grid
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        main_container.grid_columnconfigure(0, weight=1)
        main_container.grid_columnconfigure(1, weight=1)
        main_container.grid_rowconfigure(2, weight=1)
        
        # === CLUSTER COMMUNICATION CONTROL ===
        comm_frame = ttk.LabelFrame(main_container, text="Cluster Communication Control", padding=10)
        comm_frame.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 10))
        
        # Communication type selection
        type_frame = ttk.Frame(comm_frame)
        type_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(type_frame, text="Communication Type:").pack(side=tk.LEFT, padx=5)
        self.comm_type = tk.StringVar(value="intra_cluster")
        ttk.Radiobutton(type_frame, text="Intra-Cluster (Leaderâ†”Child)", 
                       variable=self.comm_type, value="intra_cluster",
                       command=self._update_comm_options).pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(type_frame, text="Inter-Cluster (Leaderâ†”Leader)", 
                       variable=self.comm_type, value="inter_cluster",
                       command=self._update_comm_options).pack(side=tk.LEFT, padx=10)
        
        # Source/destination selection
        select_frame = ttk.Frame(comm_frame)
        select_frame.pack(fill=tk.X, pady=5)
        
        # Source cluster/satellite
        ttk.Label(select_frame, text="From:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.source_cluster = tk.StringVar()
        self.source_combo = ttk.Combobox(select_frame, textvariable=self.source_cluster, 
                                         width=20, state="readonly")
        self.source_combo.grid(row=0, column=1, padx=5)
        self.source_combo.bind('<<ComboboxSelected>>', self._update_target_options)
        
        # Target selection
        ttk.Label(select_frame, text="To:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.target_sat = tk.StringVar()
        self.target_combo = ttk.Combobox(select_frame, textvariable=self.target_sat, 
                                        width=20, state="readonly")
        self.target_combo.grid(row=0, column=3, padx=5)
        
        # Message content
        ttk.Label(select_frame, text="Message:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.message_text = tk.StringVar(value="Hello from cluster")
        msg_entry = ttk.Entry(select_frame, textvariable=self.message_text, width=50)
        msg_entry.grid(row=1, column=1, columnspan=3, padx=5, pady=5)
        
        # Control buttons
        btn_frame = ttk.Frame(comm_frame)
        btn_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(btn_frame, text="Send Message", 
                  command=self._send_message).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Check Access", 
                  command=self._check_access).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Start Auto Communication", 
                  command=self._toggle_auto_comm).pack(side=tk.LEFT, padx=5)
        
        self.comm_status = ttk.Label(btn_frame, text="Status: Ready", foreground="green")
        self.comm_status.pack(side=tk.LEFT, padx=20)

        # Add a button to launch this in your Communication tab:
        ttk.Button(comm_frame, text="Launch Visualization Window", 
                 command=self.launch_communication_visualization).pack(side=tk.LEFT, padx=5)  
        
        # === ACCESS STATUS DISPLAY ===
        access_frame = ttk.LabelFrame(main_container, text="Communication Access Status", padding=10)
        access_frame.grid(row=1, column=0, sticky='nsew', padx=(0, 5))
        
        # Access matrix display
        self.access_tree = ttk.Treeview(access_frame, columns=('Link', 'Status', 'Quality', 'Last Message'), 
                                        show='tree headings', height=8)
        self.access_tree.heading('#0', text='')
        self.access_tree.heading('Link', text='Communication Link')
        self.access_tree.heading('Status', text='Status')
        self.access_tree.heading('Quality', text='Link Quality')
        self.access_tree.heading('Last Message', text='Last Message')
        
        self.access_tree.column('#0', width=0, stretch=False)
        self.access_tree.column('Link', width=200)
        self.access_tree.column('Status', width=80)
        self.access_tree.column('Quality', width=80)
        self.access_tree.column('Last Message', width=150)
        
        self.access_tree.pack(fill=tk.BOTH, expand=True)
        
        # === MESSAGE HISTORY ===
        history_frame = ttk.LabelFrame(main_container, text="Message History", padding=10)
        history_frame.grid(row=1, column=1, sticky='nsew', padx=(5, 0))
        
        # Message history display
        self.history_tree = ttk.Treeview(history_frame, 
                                         columns=('Time', 'From', 'To', 'Message', 'Status'), 
                                         show='tree headings', height=8)
        self.history_tree.heading('#0', text='')
        self.history_tree.heading('Time', text='Time (min)')
        self.history_tree.heading('From', text='Sender')
        self.history_tree.heading('To', text='Receiver')
        self.history_tree.heading('Message', text='Message')
        self.history_tree.heading('Status', text='Status')
        
        self.history_tree.column('#0', width=0, stretch=False)
        self.history_tree.column('Time', width=80)
        self.history_tree.column('From', width=100)
        self.history_tree.column('To', width=100)
        self.history_tree.column('Message', width=150)
        self.history_tree.column('Status', width=80)
        
        self.history_tree.pack(fill=tk.BOTH, expand=True)
        
        # === COMMUNICATION PLOTS ===
        plot_frame = ttk.LabelFrame(main_container, text="Communication Timeline & Analysis", padding=10)
        plot_frame.grid(row=2, column=0, columnspan=2, sticky='nsew', pady=(10, 0))
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(12, 4))
        self.fig.suptitle('Cluster Communication Analysis', fontsize=12)
        
        # Create subplots
        self.ax1 = self.fig.add_subplot(131)  # Access timeline
        self.ax2 = self.fig.add_subplot(132)  # Communication heatmap
        self.ax3 = self.fig.add_subplot(133)  # Message throughput
        
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Control buttons for plots
        plot_btn_frame = ttk.Frame(plot_frame)
        plot_btn_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(plot_btn_frame, text="Update Plots", 
                  command=self._update_plots).pack(side=tk.LEFT, padx=5)
        ttk.Button(plot_btn_frame, text="Export Data", 
                  command=self._export_comm_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(plot_btn_frame, text="Clear History", 
                  command=self._clear_history).pack(side=tk.LEFT, padx=5)
        
        # Initialize with empty plots
        self._initialize_plots()
        
    def _initialize_plots(self):
        """Initialize empty plots"""
        # Access timeline
        self.ax1.set_xlabel('Time (min)')
        self.ax1.set_ylabel('Communication Links')
        self.ax1.set_title('Access Timeline')
        self.ax1.grid(True, alpha=0.3)
        
        # Communication heatmap
        self.ax2.set_title('Link Quality Heatmap')
        self.ax2.set_xlabel('Receiver')
        self.ax2.set_ylabel('Sender')
        
        # Message throughput
        self.ax3.set_xlabel('Time (min)')
        self.ax3.set_ylabel('Messages/min')
        self.ax3.set_title('Message Throughput')
        self.ax3.grid(True, alpha=0.3)
        
        self.canvas.draw()
        
    def _update_comm_options(self):
        """Update communication options based on type"""
        if not hasattr(self.parent, 'constellation_tab'):
            return
            
        clusters = self.parent.constellation_tab.clusters
        
        if self.comm_type.get() == "intra_cluster":
            # Show clusters for source
            cluster_names = [c['name'] for c in clusters]
            self.source_combo['values'] = cluster_names
            if cluster_names:
                self.source_combo.current(0)
            self._update_target_options()
        else:
            # Show leaders for inter-cluster
            leaders = [c['leader'] for c in clusters if c.get('leader')]
            self.source_combo['values'] = leaders
            self.target_combo['values'] = [l for l in leaders if l != self.source_cluster.get()]
            
    def _update_target_options(self, event=None):
        """Update target options based on source selection"""
        if self.comm_type.get() == "intra_cluster":
            # Get children of selected cluster
            cluster_name = self.source_cluster.get()
            if hasattr(self.parent, 'constellation_tab'):
                cluster = next((c for c in self.parent.constellation_tab.clusters 
                              if c['name'] == cluster_name), None)
                if cluster:
                    # Show all satellites in cluster
                    self.target_combo['values'] = cluster.get('satellites', [])
                    if cluster.get('satellites'):
                        self.target_combo.current(0)
                        
    def _send_message(self):
        """Send a message based on current selection"""
        source = self.source_cluster.get()
        target = self.target_sat.get()
        message = self.message_text.get()
        
        if not source or not target:
            messagebox.showwarning("Selection Required", "Please select source and target")
            return
            
        # Check access first
        has_access = self._check_10s_access(source, target)
        
        if not has_access:
            self.comm_status.config(text="Status: No Access - 10s window required", foreground="red")
            return
            
        # Add to message history
        self.message_history.append({
            'time': self.sim_time,
            'from': source,
            'to': target,
            'message': message,
            'status': 'Sent' if has_access else 'Failed'
        })
        
        # Update history display
        self._update_history_display()
        
        # Update status
        self.comm_status.config(text=f"Status: Message sent to {target}", foreground="green")
        
        # Log the message
        self.parent.add_log(f"Message sent from {source} to {target}: {message}")

    # Add this method to your main GUI class or communication_tab.py:

    def launch_communication_visualization(self):
        """Launch the communication visualization window"""
        try:
            from communication_visualization import CommunicationVisualizer
            
            # Build cluster manager data from current configuration
            class SimpleClusterManager:
                def __init__(self):
                    self.clusters = {}
            
            cluster_mgr = SimpleClusterManager()
            
            # Get clusters from constellation tab
            if hasattr(self, 'constellation_tab'):
                for cluster in self.constellation_tab.clusters:
                    cluster_mgr.clusters[cluster['name']] = {
                        'leader': None,
                        'children': []
                    }
                    
                    # Find leader and children satellites
                    for sat in self.constellation_tab.satellites:
                        if sat.get('cluster') == cluster['name']:
                            if sat.get('role') == 'leader':
                                # Create a simple object with model_tag
                                leader_obj = type('obj', (object,), {'model_tag': sat['name']})()
                                cluster_mgr.clusters[cluster['name']]['leader'] = leader_obj
                            elif sat.get('role') == 'child':
                                child_obj = type('obj', (object,), {'model_tag': sat['name']})()
                                cluster_mgr.clusters[cluster['name']]['children'].append(child_obj)
            
            # Launch visualization window
            viz_window = CommunicationVisualizer(
                parent_window=self.root if hasattr(self, 'root') else None,
                cluster_manager=cluster_mgr
            )
            
            self.add_log("Launched communication visualization window")
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not launch visualization: {e}")

  
    def _check_10s_access(self, source, target):
        """Check if 10-second continuous access is available"""
        # Simulate access check - in real implementation, this would check actual access data
        # For demo, use random with higher probability for same cluster
        if self.comm_type.get() == "intra_cluster":
            return np.random.random() > 0.2  # 80% success for intra-cluster
        else:
            return np.random.random() > 0.5  # 50% success for inter-cluster
            
    def _check_access(self):
        """Check current access status for all links"""
        if not hasattr(self.parent, 'constellation_tab'):
            return
            
        # Clear current display
        for item in self.access_tree.get_children():
            self.access_tree.delete(item)
            
        clusters = self.parent.constellation_tab.clusters
        
        for cluster in clusters:
            leader = cluster.get('leader')
            if not leader:
                continue
                
            # Check intra-cluster links
            for child in cluster.get('children', []):
                has_access = np.random.random() > 0.3  # Simulate access
                quality = np.random.randint(60, 100) if has_access else 0
                
                self.access_tree.insert('', 'end', values=(
                    f"{leader} â†’ {child}",
                    "Connected" if has_access else "No Access",
                    f"{quality}%",
                    "10 min ago" if has_access else "Never"
                ))
                
        # Check inter-cluster links
        if len(clusters) > 1:
            for i, cluster1 in enumerate(clusters[:-1]):
                for cluster2 in clusters[i+1:]:
                    leader1 = cluster1.get('leader')
                    leader2 = cluster2.get('leader')
                    if leader1 and leader2:
                        has_access = np.random.random() > 0.6  # Lower probability
                        quality = np.random.randint(40, 80) if has_access else 0
                        
                        self.access_tree.insert('', 'end', values=(
                            f"{leader1} â†” {leader2}",
                            "Connected" if has_access else "No Access",
                            f"{quality}%",
                            "15 min ago" if has_access else "Never"
                        ))
                        
    def _toggle_auto_comm(self):
        """Toggle automatic communication simulation"""
        self.communication_active = not self.communication_active
        
        if self.communication_active:
            self.comm_status.config(text="Status: Auto-communication active", foreground="blue")
            # Start communication thread
            threading.Thread(target=self._auto_communicate, daemon=True).start()
        else:
            self.comm_status.config(text="Status: Auto-communication stopped", foreground="orange")
            
    def _auto_communicate(self):
        """Automatic communication simulation"""
        while self.communication_active:
            # Simulate time progression
            self.sim_time += 0.5
            
            # Random communication events
            if np.random.random() > 0.7:  # 30% chance per cycle
                self._simulate_random_message()
                
            # Update displays
            self.parent.root.after(0, self._update_plots)
            self.parent.root.after(0, self._check_access)
            
            # Sleep for simulation
            import time
            time.sleep(2)
            
    def _simulate_random_message(self):
        """Simulate a random message between satellites"""
        if not hasattr(self.parent, 'constellation_tab'):
            return
            
        clusters = self.parent.constellation_tab.clusters
        if not clusters:
            return
            
        # Pick random cluster
        cluster = np.random.choice(clusters)
        leader = cluster.get('leader')
        children = cluster.get('children', [])
        
        if leader and children:
            # Random intra-cluster message
            if np.random.random() > 0.5:
                # Leader to child
                child = np.random.choice(children)
                self.message_history.append({
                    'time': self.sim_time,
                    'from': leader,
                    'to': child,
                    'message': f"Command #{len(self.message_history)}",
                    'status': 'Sent'
                })
            else:
                # Child to leader
                child = np.random.choice(children)
                self.message_history.append({
                    'time': self.sim_time,
                    'from': child,
                    'to': leader,
                    'message': f"Telemetry #{len(self.message_history)}",
                    'status': 'Sent'
                })
                
        self._update_history_display()
        
    def _update_history_display(self):
        """Update message history display"""
        # Clear current display
        for item in self.history_tree.get_children():
            self.history_tree.delete(item)
            
        # Add recent messages (last 20)
        for msg in self.message_history[-20:]:
            self.history_tree.insert('', 0, values=(
                f"{msg['time']:.1f}",
                msg['from'],
                msg['to'],
                msg['message'],
                msg['status']
            ))
            
    def _update_plots(self):
        """Update communication plots"""
        # Clear plots
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        
        if not self.message_history:
            self._initialize_plots()
            return
            
        # === Plot 1: Access Timeline ===
        self.ax1.set_title('Communication Access Timeline')
        
        # Create timeline for each link
        links = {}
        for msg in self.message_history:
            link = f"{msg['from']}â†’{msg['to']}"
            if link not in links:
                links[link] = []
            links[link].append(msg['time'])
            
        y_pos = 0
        for link, times in links.items():
            # Plot access windows
            for t in times:
                self.ax1.plot([t-0.5, t+0.5], [y_pos, y_pos], 'g-', linewidth=3, alpha=0.6)
                self.ax1.plot(t, y_pos, 'ro', markersize=6)  # Message marker
            self.ax1.text(-1, y_pos, link[:15], fontsize=8, va='center')
            y_pos += 1
            
        self.ax1.set_xlabel('Time (min)')
        self.ax1.set_ylabel('Communication Links')
        self.ax1.set_xlim(max(0, self.sim_time - 30), self.sim_time + 1)
        self.ax1.set_ylim(-0.5, max(0.5, y_pos - 0.5))
        self.ax1.grid(True, alpha=0.3)
        
        # === Plot 2: Communication Heatmap ===
        self.ax2.set_title('Communication Frequency Heatmap')
        
        # Build communication matrix
        all_sats = set()
        for msg in self.message_history:
            all_sats.add(msg['from'])
            all_sats.add(msg['to'])
        all_sats = sorted(list(all_sats))
        
        if len(all_sats) > 1:
            comm_matrix = np.zeros((len(all_sats), len(all_sats)))
            for msg in self.message_history:
                i = all_sats.index(msg['from'])
                j = all_sats.index(msg['to'])
                comm_matrix[i, j] += 1
                
            im = self.ax2.imshow(comm_matrix, cmap='YlOrRd', aspect='auto')
            self.ax2.set_xticks(range(len(all_sats)))
            self.ax2.set_yticks(range(len(all_sats)))
            self.ax2.set_xticklabels([s[:8] for s in all_sats], rotation=45, ha='right')
            self.ax2.set_yticklabels([s[:8] for s in all_sats])
            self.ax2.set_xlabel('Receiver')
            self.ax2.set_ylabel('Sender')
            
            # Add colorbar
            self.fig.colorbar(im, ax=self.ax2, fraction=0.046, pad=0.04)
            
        # === Plot 3: Message Throughput ===
        self.ax3.set_title('Message Throughput Over Time')
        
        # Calculate throughput in bins
        if self.message_history:
            times = [msg['time'] for msg in self.message_history]
            bins = np.arange(0, max(times) + 5, 5)  # 5-minute bins
            hist, bin_edges = np.histogram(times, bins=bins)
            throughput = hist / 5.0  # Messages per minute
            
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            self.ax3.bar(bin_centers, throughput, width=4, alpha=0.7, color='blue')
            self.ax3.plot(bin_centers, throughput, 'r-', linewidth=2)
            
        self.ax3.set_xlabel('Time (min)')
        self.ax3.set_ylabel('Messages/min')
        self.ax3.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.canvas.draw()
        
    def _export_comm_data(self):
        """Export communication data to CSV"""
        from tkinter import filedialog
        import csv
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            with open(filename, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['time', 'from', 'to', 'message', 'status'])
                writer.writeheader()
                writer.writerows(self.message_history)
            messagebox.showinfo("Export Complete", f"Data exported to {filename}")
            
    def _clear_history(self):
        """Clear message history"""
        if messagebox.askyesno("Clear History", "Clear all message history?"):
            self.message_history.clear()
            self._update_history_display()
            self._update_plots()
            self.parent.add_log("Communication history cleared")
            
    def update_for_clusters(self):
        """Update communication tab when clusters change"""
        self._update_comm_options()
        self._check_access()
        self._update_plots()