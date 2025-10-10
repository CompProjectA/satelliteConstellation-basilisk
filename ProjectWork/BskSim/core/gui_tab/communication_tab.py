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
        
        
        self.cm = getattr(parent, "cluster_manager", None)

        self._create_widgets()
        self.refresh_from_gui_clusters()
        

        
    def _create_widgets(self):
        main_container = ttk.Frame(self.parent_frame)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        main_container.grid_columnconfigure(0, weight=1)
        main_container.grid_columnconfigure(1, weight=1)
        main_container.grid_rowconfigure(2, weight=1)

        # === CLUSTER COMMUNICATION CONTROL (top row) ===
        comm_frame = ttk.LabelFrame(main_container, text="Cluster Communication Control", padding=10)
        comm_frame.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 10))

        type_frame = ttk.Frame(comm_frame)
        type_frame.pack(fill=tk.X, pady=5)
        ttk.Label(type_frame, text="Communication Type:").pack(side=tk.LEFT, padx=5)
        self.comm_type = tk.StringVar(value="intra_cluster")
        ttk.Radiobutton(
            type_frame, text="Intra-Cluster (Leader→Child)",
            variable=self.comm_type, value="intra_cluster",
            command=self._update_comm_options
        ).pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(
            type_frame, text="Inter-Cluster (Leader↔Leader)",
            variable=self.comm_type, value="inter_cluster",
            command=self._update_comm_options
        ).pack(side=tk.LEFT, padx=10)

        select_frame = ttk.Frame(comm_frame)
        select_frame.pack(fill=tk.X, pady=5)
        ttk.Label(select_frame, text="From:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.source_cluster = tk.StringVar()
        self.source_combo = ttk.Combobox(select_frame, textvariable=self.source_cluster, width=20, state="readonly")
        self.source_combo.grid(row=0, column=1, padx=5)
        self.source_combo.bind('<<ComboboxSelected>>', self._update_target_options)

        ttk.Label(select_frame, text="To:").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.target_sat = tk.StringVar()
        self.target_combo = ttk.Combobox(select_frame, textvariable=self.target_sat, width=20, state="readonly")
        self.target_combo.grid(row=0, column=3, padx=5)

        ttk.Label(select_frame, text="Message:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.message_text = tk.StringVar(value="Hello from cluster")
        msg_entry = ttk.Entry(select_frame, textvariable=self.message_text, width=50)
        msg_entry.grid(row=1, column=1, columnspan=3, padx=5, pady=5)

        btn_frame = ttk.Frame(comm_frame)
        btn_frame.pack(fill=tk.X, pady=10)
        ttk.Button(btn_frame, text="Send Message", command=self._send_message).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Check Access", command=self._check_access).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Start Auto Communication", command=self._toggle_auto_comm).pack(side=tk.LEFT, padx=5)
        self.comm_status = ttk.Label(btn_frame, text="Status: Ready", foreground="green")
        self.comm_status.pack(side=tk.LEFT, padx=20)
        ttk.Button(comm_frame, text="Launch Visualization Window", command=self.launch_communication_visualization)\
        .pack(side=tk.LEFT, padx=5)

        # === NOTEBOOK (row 1) ===
        self.subtabs = ttk.Notebook(main_container)
        self.subtabs.grid(row=1, column=0, columnspan=2, sticky='nsew', pady=(5, 0))
        console = ttk.Frame(self.subtabs)
        viz = ttk.Frame(self.subtabs)
        self.subtabs.add(console, text="Console")
        self.subtabs.add(viz, text="Visualization")

        # --- Console: Access + History ---
        access_frame = ttk.LabelFrame(console, text="Communication Access Status", padding=10)
        access_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
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

        history_frame = ttk.LabelFrame(console, text="Message History", padding=10)
        history_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(5, 0))
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

        # --- Visualization: Plots ---
        plot_frame = ttk.LabelFrame(viz, text="Communication Timeline & Analysis", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(12, 4))
        self.fig.suptitle('Cluster Communication Analysis', fontsize=12)
        self.ax1 = self.fig.add_subplot(131)
        self.ax2 = self.fig.add_subplot(132)
        self.ax3 = self.fig.add_subplot(133)
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        plot_btn_frame = ttk.Frame(plot_frame)
        plot_btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(plot_btn_frame, text="Update Plots", command=self._update_plots).pack(side=tk.LEFT, padx=5)
        ttk.Button(plot_btn_frame, text="Export Data", command=self._export_comm_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(plot_btn_frame, text="Clear History", command=self._clear_history).pack(side=tk.LEFT, padx=5)

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
        ct = getattr(self.parent, "constellation_tab", None)
        if not ct:
            return
        clusters = ct.clusters

        if self.comm_type.get() == "intra_cluster":
            names = [c['name'] for c in clusters]
            self.source_combo['values'] = names
            if names:
                self.source_combo.current(0)
            self._update_target_options()
        else:
            # inter-cluster: show leaders
            leaders = [c['leader'] for c in clusters if c.get('leader')]
            self.source_combo['values'] = leaders
            if leaders:
                self.source_combo.current(0)
                self._update_target_options()


    def _update_target_options(self, *_):
        ct = getattr(self.parent, "constellation_tab", None)
        if not ct: 
            self.target_combo['values'] = []
            return
        clusters = ct.clusters

        if self.comm_type.get() == "intra_cluster":
            src_cluster_name = self.source_cluster.get()
            c = next((x for x in clusters if x['name'] == src_cluster_name), None)
            vals = [c['leader']] + c['children'] if c else []
        else:
            src_leader = self.source_cluster.get()
            vals = [c['leader'] for c in clusters if c.get('leader') and c['leader'] != src_leader]

        self.target_combo['values'] = vals
        if vals:
            self.target_sat.set(vals[0])

  

    def _send_message(self):
        source = self.source_cluster.get()
        target = self.target_sat.get()
        message = self.message_text.get().strip()
        if not source or not target or not message:
            messagebox.showwarning("Selection Required", "Please select source/target and enter a message")
            return

        # Gate on access if intra
        if not self._check_10s_access(source, target):
            self.comm_status.config(text="Status: No Access - 10s window required", foreground="red")
            return

        sent = False
        if self.cm:
            ct = getattr(self.parent, "constellation_tab", None)
            if ct:
                if self.comm_type.get() == "intra_cluster":
                    # source = cluster name; target is child name (or leader)
                    cluster = next((c for c in ct.clusters if c['name'] == source), None)
                    if cluster:
                        if target in cluster['children']:
                            idx = cluster['children'].index(target)
                            sent = self.cm.send_message_in_cluster(
                                cluster_name=cluster['name'],
                                message_content=message,
                                time_min=self.sim_time,
                                from_leader=True,
                                to_child_index=idx,
                                require_access=True
                            )
                        else:
                            # child -> leader (quick pass; extend UI later if needed)
                            sent = self.cm.send_message_in_cluster(
                                cluster_name=cluster['name'],
                                message_content=message,
                                time_min=self.sim_time,
                                from_leader=False,
                                to_child_index=0,   # pick first child by default
                                require_access=False
                            )
                else:
                    # inter-cluster: source/target are leader names; map to cluster names
                    def cluster_by_leader(name):
                        for c in ct.clusters:
                            if c.get('leader') == name:
                                return c['name']
                        return None
                    a = cluster_by_leader(source)
                    b = cluster_by_leader(target)
                    if a and b:
                        sent = self.cm.send_inter_cluster_message(a, b, message, self.sim_time)
        else:
            # fallback simulated send
            sent = True

        status_text = "Sent" if sent else "Failed"
        self.message_history.append({'time': self.sim_time, 'from': source, 'to': target, 'message': message, 'status': status_text})
        self._update_history_display()
        color = "green" if sent else "red"
        self.comm_status.config(text=f"Status: {status_text} to {target}", foreground=color)
        self.parent.add_log(f"Message {status_text.lower()} from {source} to {target}: {message}")




    def update_cluster_communication_display(self):
        """Update communication display based on current clusters"""
        # Clear existing display
        for item in self.access_tree.get_children():
            self.access_tree.delete(item)
        
        if not hasattr(self.parent, 'constellation_tab'):
            return
        
        clusters = self.parent.constellation_tab.clusters
        satellites = self.parent.constellation_tab.satellites
        
        # Build cluster index mapping
        sat_to_index = {sat['name']: i for i, sat in enumerate(satellites)}
        
        for cluster in clusters:
            cluster_name = cluster['name']
            leader_name = cluster.get('leader')
            children_names = cluster.get('children', [])
            
            if not leader_name:
                continue
            
            # Add cluster header
            cluster_item = self.access_tree.insert('', 'end', values=(
                f"[{cluster_name}]",
                "CLUSTER",
                "",
                ""
            ))
            
            # Add intra-cluster links
            for child_name in children_names:
                # Simulate communication status
                has_access = np.random.random() > 0.2  # 80% chance
                quality = np.random.randint(70, 100) if has_access else 0
                
                self.access_tree.insert(cluster_item, 'end', values=(
                    f"  {leader_name} ↔ {child_name}",
                    "Active" if has_access else "No Link",
                    f"{quality}%" if has_access else "0%",
                    "Real-time" if has_access else "N/A"
                ))

    def launch_communication_visualization(self):
        try:
            from communication_visualization import CommunicationVisualizer
            ct = getattr(self.parent, "constellation_tab", None)
            if not ct:
                messagebox.showerror("Error", "No constellation tab available.")
                return

            class SimpleClusterManager:
                def __init__(self): self.clusters = {}

            cluster_mgr = SimpleClusterManager()
            for cluster in ct.clusters:
                cluster_mgr.clusters[cluster['name']] = {'leader': None, 'children': []}
                for sat in ct.satellites:
                    if sat.get('cluster') == cluster['name']:
                        stub = type('obj', (object,), {'model_tag': sat['name']})()
                        if sat.get('role') == 'leader':
                            cluster_mgr.clusters[cluster['name']]['leader'] = stub
                        elif sat.get('role') == 'child':
                            cluster_mgr.clusters[cluster['name']]['children'].append(stub)

            CommunicationVisualizer(parent_window=self.parent.root, cluster_manager=cluster_mgr)
            if hasattr(self.parent, 'add_log'):
                self.parent.add_log("Launched communication visualization window")

        except Exception as e:
            messagebox.showerror("Error", f"Could not launch visualization: {e}")



    def refresh_from_gui_clusters(self):
        ct = getattr(self.parent, "constellation_tab", None)
        clusters = ct.clusters if ct else []
        names = [c['name'] for c in clusters]
        self.source_combo['values'] = names
        if names and not self.source_cluster.get():
            self.source_cluster.set(names[0])
            self._update_comm_options()


    def _check_10s_access(self, source, target):
        # Intra: source = cluster name, target = child/leader name
        if self.cm and self.comm_type.get() == "intra_cluster":
            ct = getattr(self.parent, "constellation_tab", None)
            if not ct: return False
            cluster = next((c for c in ct.clusters if c['name'] == source), None)
            if not cluster: return False
            # only gate leader->child by access (you can extend for child->leader later)
            if target in cluster['children']:
                idx = cluster['children'].index(target)
                return self.cm.has_access_child(cluster['name'], idx)
            return True  # leader selected as target -> treat as local
        # Inter: keep as probabilistic unless you add a manager-level access API for leaders
        return np.random.random() > 0.5

            
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
                    f"{leader} → {child}",
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
                            f"{leader1} ↔ {leader2}",
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
            link = f"{msg['from']}→{msg['to']}"
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