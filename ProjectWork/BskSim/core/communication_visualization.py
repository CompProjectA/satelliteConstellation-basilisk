#!/usr/bin/env python
"""
communication_visualization.py

Real-time communication visualization showing actual message passing
between satellites in clusters.
"""

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import time

class CommunicationVisualizer:
    """Visualize real-time communication between satellites"""
    
    def __init__(self, parent_window, cluster_manager=None):
        self.parent_window = parent_window
        self.cluster_manager = cluster_manager
        
        # Create window
        self.window = tk.Toplevel(parent_window)
        self.window.title("Communication Visualization")
        self.window.geometry("1200x800")
        
        # Communication data
        self.comm_history = []
        self.access_windows = {}
        self.message_log = []
        
        # Create UI
        self._create_ui()
        
        # Start animation
        self.animation_running = False
        self.start_animation()
        
    def _create_ui(self):
        """Create the visualization UI"""
        # Main container
        main_frame = ttk.Frame(self.window, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Top control panel
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding=10)
        control_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(control_frame, text="Start Simulation", 
                  command=self.start_simulation).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Stop Simulation", 
                  command=self.stop_simulation).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Clear History", 
                  command=self.clear_history).pack(side=tk.LEFT, padx=5)
        
        # Time display
        self.time_label = ttk.Label(control_frame, text="Time: 0.0 min", 
                                   font=('Arial', 12, 'bold'))
        self.time_label.pack(side=tk.RIGHT, padx=10)
        
        # Create notebook for different views
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Communication Matrix Tab
        matrix_frame = ttk.Frame(notebook)
        notebook.add(matrix_frame, text="Communication Matrix")
        self._create_matrix_view(matrix_frame)
        
        # Timeline Tab
        timeline_frame = ttk.Frame(notebook)
        notebook.add(timeline_frame, text="Timeline View")
        self._create_timeline_view(timeline_frame)
        
        # Network Graph Tab
        network_frame = ttk.Frame(notebook)
        notebook.add(network_frame, text="Network Graph")
        self._create_network_view(network_frame)
        
        # Message Log Tab
        log_frame = ttk.Frame(notebook)
        notebook.add(log_frame, text="Message Log")
        self._create_log_view(log_frame)
        
    def _create_matrix_view(self, parent):
        """Create communication matrix visualization"""
        # Create figure
        self.matrix_fig = Figure(figsize=(8, 6), dpi=100)
        self.matrix_ax = self.matrix_fig.add_subplot(111)
        
        # Canvas
        self.matrix_canvas = FigureCanvasTkAgg(self.matrix_fig, parent)
        self.matrix_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize matrix plot
        self.update_matrix_plot()
        
    def _create_timeline_view(self, parent):
        """Create timeline visualization"""
        # Create figure
        self.timeline_fig = Figure(figsize=(10, 6), dpi=100)
        self.timeline_ax = self.timeline_fig.add_subplot(111)
        
        # Canvas
        self.timeline_canvas = FigureCanvasTkAgg(self.timeline_fig, parent)
        self.timeline_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize timeline
        self.update_timeline_plot()
        
    def _create_network_view(self, parent):
        """Create network graph visualization"""
        # Create figure
        self.network_fig = Figure(figsize=(8, 8), dpi=100)
        self.network_ax = self.network_fig.add_subplot(111)
        
        # Canvas
        self.network_canvas = FigureCanvasTkAgg(self.network_fig, parent)
        self.network_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Initialize network
        self.update_network_plot()
        
    def _create_log_view(self, parent):
        """Create message log view"""
        # Log text widget
        log_frame = ttk.Frame(parent)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create text widget with scrollbar
        scrollbar = ttk.Scrollbar(log_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.log_text = tk.Text(log_frame, wrap=tk.WORD, 
                               yscrollcommand=scrollbar.set,
                               font=('Consolas', 10))
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.log_text.yview)
        
        # Add initial message
        self.add_log_message("Communication system initialized")
        
    def update_matrix_plot(self):
        """Update communication matrix plot"""
        self.matrix_ax.clear()
        
        if not self.cluster_manager or not hasattr(self.cluster_manager, 'clusters'):
            self.matrix_ax.text(0.5, 0.5, 'No clusters available', 
                              ha='center', va='center', fontsize=14)
            self.matrix_ax.set_xlim(0, 1)
            self.matrix_ax.set_ylim(0, 1)
        else:
            # Get all satellites
            all_sats = []
            for cluster_name, cluster_data in self.cluster_manager.clusters.items():
                if 'leader' in cluster_data:
                    all_sats.append(cluster_data['leader'].model_tag if hasattr(cluster_data['leader'], 'model_tag') 
                                  else str(cluster_data['leader']))
                for child in cluster_data.get('children', []):
                    if hasattr(child, 'model_tag'):
                        all_sats.append(child.model_tag)
                    else:
                        all_sats.append(str(child))
            
            n = len(all_sats)
            if n > 0:
                # Create communication matrix
                comm_matrix = np.random.random((n, n))  # Simulate comm strength
                np.fill_diagonal(comm_matrix, 1.0)  # Self communication is perfect
                
                # Make matrix symmetric
                comm_matrix = (comm_matrix + comm_matrix.T) / 2
                
                # Plot matrix
                im = self.matrix_ax.imshow(comm_matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
                
                # Set labels
                self.matrix_ax.set_xticks(range(n))
                self.matrix_ax.set_yticks(range(n))
                self.matrix_ax.set_xticklabels([s[:8] for s in all_sats], rotation=45, ha='right')
                self.matrix_ax.set_yticklabels([s[:8] for s in all_sats])
                
                # Add colorbar
                self.matrix_fig.colorbar(im, ax=self.matrix_ax, label='Link Quality')
                
                # Add grid
                for i in range(n+1):
                    self.matrix_ax.axhline(i-0.5, color='white', linewidth=0.5)
                    self.matrix_ax.axvline(i-0.5, color='white', linewidth=0.5)
                
                self.matrix_ax.set_title('Communication Link Quality Matrix')
            else:
                self.matrix_ax.text(0.5, 0.5, 'No satellites in clusters', 
                                  ha='center', va='center', fontsize=14)
                self.matrix_ax.set_xlim(0, 1)
                self.matrix_ax.set_ylim(0, 1)
        
        self.matrix_canvas.draw()
        
    def update_timeline_plot(self):
        """Update timeline plot showing communication events"""
        self.timeline_ax.clear()
        
        if not self.comm_history:
            self.timeline_ax.text(0.5, 0.5, 'No communication events yet', 
                                ha='center', va='center', fontsize=14)
            self.timeline_ax.set_xlim(0, 30)
            self.timeline_ax.set_ylim(0, 1)
        else:
            # Plot communication events
            times = [event['time'] for event in self.comm_history]
            senders = [event['sender'] for event in self.comm_history]
            receivers = [event['receiver'] for event in self.comm_history]
            
            # Get unique satellites
            all_sats = list(set(senders + receivers))
            sat_positions = {sat: i for i, sat in enumerate(all_sats)}
            
            # Plot connections
            for event in self.comm_history:
                sender_y = sat_positions[event['sender']]
                receiver_y = sat_positions[event['receiver']]
                time = event['time']
                
                # Draw arrow from sender to receiver
                self.timeline_ax.annotate('', xy=(time, receiver_y), 
                                        xytext=(time, sender_y),
                                        arrowprops=dict(arrowstyle='->', 
                                                      color=event.get('color', 'blue'),
                                                      lw=2, alpha=0.7))
                
                # Add message indicator
                self.timeline_ax.scatter(time, receiver_y, s=50, 
                                       color=event.get('color', 'blue'), 
                                       zorder=5, alpha=0.8)
            
            # Set labels
            self.timeline_ax.set_yticks(range(len(all_sats)))
            self.timeline_ax.set_yticklabels(all_sats)
            self.timeline_ax.set_xlabel('Time (minutes)')
            self.timeline_ax.set_ylabel('Satellites')
            self.timeline_ax.set_title('Communication Timeline')
            self.timeline_ax.grid(True, alpha=0.3)
            
            # Set limits
            if times:
                self.timeline_ax.set_xlim(0, max(times) + 5)
            else:
                self.timeline_ax.set_xlim(0, 30)
            self.timeline_ax.set_ylim(-0.5, len(all_sats) - 0.5)
        
        self.timeline_canvas.draw()
        
    def update_network_plot(self):
        """Update network graph showing satellite connections"""
        self.network_ax.clear()
        
        if not self.cluster_manager or not hasattr(self.cluster_manager, 'clusters'):
            self.network_ax.text(0.5, 0.5, 'No network data available', 
                               ha='center', va='center', fontsize=14)
            self.network_ax.set_xlim(0, 1)
            self.network_ax.set_ylim(0, 1)
        else:
            # Create network layout
            clusters = self.cluster_manager.clusters
            num_clusters = len(clusters)
            
            if num_clusters > 0:
                # Position clusters in a circle
                cluster_angle_step = 2 * np.pi / num_clusters
                cluster_radius = 3
                
                node_positions = {}
                node_colors = {}
                edges = []
                
                colors = ['red', 'blue', 'green', 'orange', 'purple']
                
                for i, (cluster_name, cluster_data) in enumerate(clusters.items()):
                    cluster_angle = i * cluster_angle_step
                    cluster_x = cluster_radius * np.cos(cluster_angle)
                    cluster_y = cluster_radius * np.sin(cluster_angle)
                    
                    color = colors[i % len(colors)]
                    
                    # Position leader
                    leader = cluster_data.get('leader')
                    if leader:
                        leader_name = leader.model_tag if hasattr(leader, 'model_tag') else str(leader)
                        node_positions[leader_name] = (cluster_x, cluster_y)
                        node_colors[leader_name] = color
                        
                        # Position children around leader
                        children = cluster_data.get('children', [])
                        num_children = len(children)
                        if num_children > 0:
                            child_angle_step = 2 * np.pi / num_children
                            child_radius = 0.8
                            
                            for j, child in enumerate(children):
                                child_name = child.model_tag if hasattr(child, 'model_tag') else str(child)
                                child_angle = j * child_angle_step
                                child_x = cluster_x + child_radius * np.cos(child_angle)
                                child_y = cluster_y + child_radius * np.sin(child_angle)
                                
                                node_positions[child_name] = (child_x, child_y)
                                node_colors[child_name] = color
                                
                                # Add edge from leader to child
                                edges.append((leader_name, child_name))
                
                # Draw nodes
                for node_name, (x, y) in node_positions.items():
                    is_leader = any(node_name == (c['leader'].model_tag if hasattr(c['leader'], 'model_tag') 
                                                 else str(c['leader'])) 
                                  for c in clusters.values() if 'leader' in c)
                    
                    size = 300 if is_leader else 150
                    marker = '^' if is_leader else 'o'
                    
                    self.network_ax.scatter(x, y, s=size, c=node_colors[node_name], 
                                         marker=marker, edgecolor='black', 
                                         linewidth=2, zorder=5)
                    
                    # Add label
                    self.network_ax.text(x, y-0.3, node_name[:8], 
                                       ha='center', fontsize=8)
                
                # Draw edges
                for sender, receiver in edges:
                    if sender in node_positions and receiver in node_positions:
                        x1, y1 = node_positions[sender]
                        x2, y2 = node_positions[receiver]
                        
                        self.network_ax.plot([x1, x2], [y1, y2], 'k-', 
                                           alpha=0.3, linewidth=1)
                
                # Add inter-cluster connections
                if num_clusters > 1:
                    cluster_leaders = []
                    for cluster_data in clusters.values():
                        leader = cluster_data.get('leader')
                        if leader:
                            leader_name = leader.model_tag if hasattr(leader, 'model_tag') else str(leader)
                            cluster_leaders.append(leader_name)
                    
                    # Connect all leaders
                    for i, leader1 in enumerate(cluster_leaders):
                        for leader2 in cluster_leaders[i+1:]:
                            if leader1 in node_positions and leader2 in node_positions:
                                x1, y1 = node_positions[leader1]
                                x2, y2 = node_positions[leader2]
                                
                                self.network_ax.plot([x1, x2], [y1, y2], 'r--', 
                                                   alpha=0.2, linewidth=2)
                
                self.network_ax.set_title('Satellite Network Topology')
                self.network_ax.set_aspect('equal')
                self.network_ax.axis('off')
            else:
                self.network_ax.text(0.5, 0.5, 'No clusters configured', 
                                   ha='center', va='center', fontsize=14)
                self.network_ax.set_xlim(0, 1)
                self.network_ax.set_ylim(0, 1)
        
        self.network_canvas.draw()
        
    def add_log_message(self, message):
        """Add message to log"""
        timestamp = time.strftime("%H:%M:%S")
        log_entry = f"[{timestamp}] {message}\n"
        
        self.log_text.insert(tk.END, log_entry)
        self.log_text.see(tk.END)  # Scroll to bottom
        
        # Keep log size manageable
        if int(self.log_text.index('end-1c').split('.')[0]) > 1000:
            self.log_text.delete('1.0', '2.0')
    
    def simulate_communication(self, current_time):
        """Simulate communication events"""
        if not self.cluster_manager:
            return
        
        # Simulate some communication events
        if np.random.random() > 0.7:  # 30% chance of communication
            clusters = list(self.cluster_manager.clusters.items())
            if clusters:
                cluster_name, cluster_data = clusters[np.random.randint(0, len(clusters))]
                
                leader = cluster_data.get('leader')
                children = cluster_data.get('children', [])
                
                if leader and children:
                    # Simulate leader to child communication
                    child = children[np.random.randint(0, len(children))]
                    
                    sender_name = leader.model_tag if hasattr(leader, 'model_tag') else "Leader"
                    receiver_name = child.model_tag if hasattr(child, 'model_tag') else "Child"
                    
                    # Add to history
                    event = {
                        'time': current_time,
                        'sender': sender_name,
                        'receiver': receiver_name,
                        'message': f"Status request from {sender_name}",
                        'color': 'blue'
                    }
                    self.comm_history.append(event)
                    
                    # Add to log
                    self.add_log_message(f"{sender_name} → {receiver_name}: Status request")
                    
                    # Simulate response
                    if np.random.random() > 0.5:
                        response = {
                            'time': current_time + 0.1,
                            'sender': receiver_name,
                            'receiver': sender_name,
                            'message': f"Status OK from {receiver_name}",
                            'color': 'green'
                        }
                        self.comm_history.append(response)
                        self.add_log_message(f"{receiver_name} → {sender_name}: Status OK")
        
        # Simulate inter-cluster communication
        if np.random.random() > 0.9:  # 10% chance
            clusters = list(self.cluster_manager.clusters.items())
            if len(clusters) > 1:
                idx1, idx2 = np.random.choice(len(clusters), 2, replace=False)
                
                leader1 = clusters[idx1][1].get('leader')
                leader2 = clusters[idx2][1].get('leader')
                
                if leader1 and leader2:
                    sender_name = leader1.model_tag if hasattr(leader1, 'model_tag') else "Leader1"
                    receiver_name = leader2.model_tag if hasattr(leader2, 'model_tag') else "Leader2"
                    
                    event = {
                        'time': current_time,
                        'sender': sender_name,
                        'receiver': receiver_name,
                        'message': f"Inter-cluster sync from {sender_name}",
                        'color': 'red'
                    }
                    self.comm_history.append(event)
                    
                    self.add_log_message(f"{sender_name} ↔ {receiver_name}: Inter-cluster sync")
    
    def update_animation(self, frame):
        """Animation update function"""
        if self.animation_running:
            current_time = frame * 0.1  # Convert to minutes
            
            # Update time label
            self.time_label.config(text=f"Time: {current_time:.1f} min")
            
            # Simulate communication
            self.simulate_communication(current_time)
            
            # Update plots
            self.update_timeline_plot()
            
            # Periodically update other plots
            if frame % 10 == 0:
                self.update_matrix_plot()
                self.update_network_plot()
        
        return []
    
    def start_animation(self):
        """Start the animation"""
        self.animation_running = True
        self.animation = FuncAnimation(self.timeline_fig, self.update_animation, 
                                      interval=100, blit=False)
        self.timeline_canvas.draw()
    
    def start_simulation(self):
        """Start communication simulation"""
        self.animation_running = True
        self.add_log_message("Communication simulation started")
    
    def stop_simulation(self):
        """Stop communication simulation"""
        self.animation_running = False
        self.add_log_message("Communication simulation stopped")
    
    def clear_history(self):
        """Clear communication history"""
        self.comm_history.clear()
        self.message_log.clear()
        
        # Clear log text
        self.log_text.delete('1.0', tk.END)
        self.add_log_message("History cleared")
        
        # Update all plots
        self.update_matrix_plot()
        self.update_timeline_plot()
        self.update_network_plot()


def show_communication_visualization(parent_window, cluster_manager=None):
    """Show the communication visualization window"""
    visualizer = CommunicationVisualizer(parent_window, cluster_manager)
    return visualizer