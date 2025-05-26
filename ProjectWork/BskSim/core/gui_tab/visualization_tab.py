#!/usr/bin/env python
"""
visualization_tab.py

Implements the Visualization tab for the Spacecraft Constellation Fault Simulator.
FIXED: Better layout, enhanced camera positioning, improved target visibility
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from .base_tab import BaseTab

class VisualizationTab(BaseTab):
    """Visualization tab implementation"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the visualization tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        super().__init__(parent_app, parent_frame)
        
        # Store references to parent app data
        self.satellites = parent_app.satellites
        
        # Create the tab UI
        self.create_tab_ui()

    def get_camera_position(self):
        """
        Get the camera position from the UI
        
        Returns:
        list: [x, y, z] position of the camera, or None if not available
        """
        try:
            return [self.camera_x_var.get(), 
                    self.camera_y_var.get(), 
                    self.camera_z_var.get()]
        except:
            # Return default position if values are not available
            return [0.0, 0.0, 5.0]  # FIXED: Higher above spacecraft for better view
    
    def get_active_camera_satellite(self):
        """
        Get the satellite that has the active camera
        
        Returns:
        dict: Satellite object with active camera, or None if no active camera
        """
        for sat in self.satellites:
            if sat["camera"]["enabled"]:
                return sat
        return None
        
    def create_tab_ui(self):
        """Create the Visualization tab UI with FIXED layout"""
        # FIXED: Main container with proper padding
        main_container = ttk.Frame(self.parent_frame, padding=10)
        main_container.pack(fill=tk.BOTH, expand=True)
        
        # Top frame for satellite selection and camera settings
        top_frame = ttk.Frame(main_container)
        top_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Satellite selection at the top
        select_frame = ttk.Frame(top_frame)
        select_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(select_frame, text="Select Satellite:").pack(side=tk.LEFT)
        
        self.viz_satellite_var = tk.StringVar()
        self.viz_satellite_combo = ttk.Combobox(select_frame, textvariable=self.viz_satellite_var, width=20)
        self.viz_satellite_combo.pack(side=tk.LEFT, padx=10)
        self.update_satellite_dropdown()
        
        # Bind selection change event
        self.viz_satellite_combo.bind('<<ComboboxSelected>>', self.on_viz_satellite_changed)
        
        # FIXED: Camera configuration with better layout
        camera_frame = ttk.LabelFrame(top_frame, text="Camera Settings (for Vizard)", padding=10)
        camera_frame.pack(fill=tk.X, pady=5)
        
        # First row: Enable camera and active camera indicator
        row1_frame = ttk.Frame(camera_frame)
        row1_frame.pack(fill=tk.X, pady=5)
        
        self.camera_enabled_var = tk.BooleanVar(value=False)
        camera_enabled_check = ttk.Checkbutton(row1_frame, text="Enable Camera for this Satellite", 
                                             variable=self.camera_enabled_var)
        camera_enabled_check.pack(side=tk.LEFT)
        
        # Active camera indicator
        active_frame = ttk.Frame(row1_frame)
        active_frame.pack(side=tk.RIGHT, padx=10)
        
        ttk.Label(active_frame, text="Active Camera:").pack(side=tk.LEFT)
        self.active_camera_label = ttk.Label(active_frame, text="None", foreground="red", font=('Segoe UI', 9, 'bold'))
        self.active_camera_label.pack(side=tk.LEFT, padx=5)
        
        # Update indicator
        self.update_active_camera_indicator()
        
        # Second row: Camera position controls
        pos_frame = ttk.Frame(camera_frame)
        pos_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(pos_frame, text="Camera Position (body frame):").pack(side=tk.LEFT)
        
        # X, Y, Z position controls in a compact layout
        xyz_frame = ttk.Frame(pos_frame)
        xyz_frame.pack(side=tk.LEFT, padx=20)
        
        ttk.Label(xyz_frame, text="X:").pack(side=tk.LEFT)
        self.camera_x_var = tk.DoubleVar(value=0.0)
        x_entry = ttk.Entry(xyz_frame, textvariable=self.camera_x_var, width=8)
        x_entry.pack(side=tk.LEFT, padx=2)
        
        ttk.Label(xyz_frame, text="Y:").pack(side=tk.LEFT, padx=(10, 0))
        self.camera_y_var = tk.DoubleVar(value=0.0)
        y_entry = ttk.Entry(xyz_frame, textvariable=self.camera_y_var, width=8)
        y_entry.pack(side=tk.LEFT, padx=2)
        
        ttk.Label(xyz_frame, text="Z:").pack(side=tk.LEFT, padx=(10, 0))
        self.camera_z_var = tk.DoubleVar(value=5.0)  # FIXED: Higher default for better view
        z_entry = ttk.Entry(xyz_frame, textvariable=self.camera_z_var, width=8)
        z_entry.pack(side=tk.LEFT, padx=2)
        
        # Field of view control
        fov_frame = ttk.Frame(pos_frame)
        fov_frame.pack(side=tk.RIGHT, padx=10)
        
        ttk.Label(fov_frame, text="FOV (deg):").pack(side=tk.LEFT)
        self.camera_fov_var = tk.DoubleVar(value=70.0)
        fov_entry = ttk.Entry(fov_frame, textvariable=self.camera_fov_var, width=8)
        fov_entry.pack(side=tk.LEFT, padx=5)
        
        # Third row: Preset buttons and controls
        controls_frame = ttk.Frame(camera_frame)
        controls_frame.pack(fill=tk.X, pady=5)
        
        # FIXED: Enhanced preset buttons with better positions for target visibility
        preset_frame = ttk.Frame(controls_frame)
        preset_frame.pack(side=tk.LEFT)
        
        ttk.Label(preset_frame, text="Position Presets:").pack(side=tk.LEFT)
        
        above_btn = ttk.Button(preset_frame, text="Above (Best for Targets)", 
                              command=lambda: self.set_camera_preset([0.0, 0.0, 5.0]))
        above_btn.pack(side=tk.LEFT, padx=2)
        
        high_btn = ttk.Button(preset_frame, text="High Above",
                             command=lambda: self.set_camera_preset([0.0, 0.0, 10.0]))
        high_btn.pack(side=tk.LEFT, padx=2)
        
        side_btn = ttk.Button(preset_frame, text="Side View",
                             command=lambda: self.set_camera_preset([5.0, 0.0, 0.0]))
        side_btn.pack(side=tk.LEFT, padx=2)
        
        earth_btn = ttk.Button(preset_frame, text="Earth View",
                             command=lambda: self.set_camera_preset([0.0, 0.0, 15.0]))
        earth_btn.pack(side=tk.LEFT, padx=2)
        
        # Control buttons
        button_frame = ttk.Frame(controls_frame)
        button_frame.pack(side=tk.RIGHT)
        
        # Update button
        update_btn = ttk.Button(button_frame, text="Apply Settings", 
                               command=self.apply_camera_settings,
                               style="Run.TButton")
        update_btn.pack(side=tk.LEFT, padx=5)
        
        # Auto-enable camera button for convenience
        auto_enable_btn = ttk.Button(button_frame, text="Auto-Enable Best Camera", 
                                    command=self.auto_enable_camera)
        auto_enable_btn.pack(side=tk.LEFT, padx=5)
        
        # FIXED: Add enhanced note about camera parameters
        note_frame = ttk.Frame(camera_frame)
        note_frame.pack(fill=tk.X, pady=(10, 0))
        note_text = ("FIXED: 'Above' position (Z=5.0) provides excellent view of spacecraft and targets.\n"
                    "Higher Z values give wider Earth view. Only one camera can be active at a time.\n"
                    "Targets will be visible as colored markers on Earth surface in Vizard.")
        ttk.Label(note_frame, text=note_text, style="Info.TLabel", wraplength=800).pack(anchor=tk.W)
        
        # FIXED: Visualization preview with proper layout
        preview_frame = ttk.LabelFrame(main_container, text="Visualization Preview", padding=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        
        # FIXED: Create a properly centered container for the matplotlib figure
        canvas_container = ttk.Frame(preview_frame)
        canvas_container.pack(fill=tk.BOTH, expand=True)
        
        # Simple satellite visualization using matplotlib with FIXED size and positioning
        self.viz_figure = Figure(figsize=(12, 8), dpi=100, facecolor='white')
        self.viz_subplot = self.viz_figure.add_subplot(111, projection='3d')
        
        # FIXED: Properly pack the canvas
        self.viz_canvas = FigureCanvasTkAgg(self.viz_figure, canvas_container)
        self.viz_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Add navigation toolbar for zoom/pan functionality
        toolbar_frame = ttk.Frame(preview_frame)
        toolbar_frame.pack(fill=tk.X, pady=(5, 0))
        
        self.viz_toolbar = NavigationToolbar2Tk(self.viz_canvas, toolbar_frame)
        self.viz_toolbar.update()
        
        # Add view control buttons
        view_controls_frame = ttk.Frame(toolbar_frame)
        view_controls_frame.pack(side=tk.RIGHT, padx=10)
        
        reset_view_btn = ttk.Button(view_controls_frame, text="Reset View", 
                                   command=self.reset_3d_view)
        reset_view_btn.pack(side=tk.LEFT, padx=2)
        
        top_view_btn = ttk.Button(view_controls_frame, text="Top View", 
                                 command=self.set_top_view)
        top_view_btn.pack(side=tk.LEFT, padx=2)
        
        side_view_btn = ttk.Button(view_controls_frame, text="Side View", 
                                  command=self.set_side_view)
        side_view_btn.pack(side=tk.LEFT, padx=2)
        
        # Draw a basic satellite preview
        self.update_visualization_preview()

    def auto_enable_camera(self):
        """Auto-enable camera for the first satellite with optimal settings"""
        try:
            # Check if any cameras are enabled
            camera_enabled = any(sat["camera"]["enabled"] for sat in self.satellites)
            
            if not camera_enabled and self.satellites:
                # Enable camera for first satellite with FIXED optimal settings
                self.satellites[0]["camera"]["enabled"] = True
                self.satellites[0]["camera"]["position"] = [0.0, 0.0, 5.0]  # FIXED: Better height
                self.satellites[0]["camera"]["fov"] = 70.0
                self.parent_app.add_log(f"Auto-enabled camera for {self.satellites[0]['name']} with optimal target visibility")
                
                # Update the UI to reflect the change
                if self.viz_satellite_var.get() == self.satellites[0]["name"]:
                    self.camera_enabled_var.set(True)
                    self.camera_x_var.set(0.0)
                    self.camera_y_var.set(0.0)
                    self.camera_z_var.set(5.0)  # FIXED: Better height
                    self.camera_fov_var.set(70.0)
                    
                # Update the active camera indicator
                self.update_active_camera_indicator()
                
                # Update constellation tab
                try:
                    self.parent_app.constellation_tab.update_satellite_listbox()
                except:
                    pass
                    
                messagebox.showinfo("Camera Enabled", 
                                  f"Auto-enabled camera for {self.satellites[0]['name']}\n"
                                  f"Position: Above spacecraft (0, 0, 5) - OPTIMAL FOR TARGETS\n"
                                  f"This provides excellent view for Vizard visualization.\n"
                                  f"Targets will be visible as colored markers on Earth.")
            else:
                messagebox.showinfo("Camera Status", 
                                  "A camera is already enabled, or no satellites are available.")
                    
        except Exception as e:
            print(f"Could not auto-enable camera: {e}")
        
    def update_satellite_dropdown(self):
        """Update the satellite dropdown with current satellites"""
        self.viz_satellite_combo['values'] = [sat['name'] for sat in self.satellites]
        if self.satellites:
            self.viz_satellite_combo.current(0)
            
    def on_viz_satellite_changed(self, event):
        """Handle visualization satellite selection change"""
        sat_name = self.viz_satellite_var.get()
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                self.load_camera_settings(i)
                break
                
    def load_camera_settings(self, index):
        """Load camera settings for the specified satellite"""
        if 0 <= index < len(self.satellites):
            camera = self.satellites[index]["camera"]
            
            # Load camera enabled state
            self.camera_enabled_var.set(camera.get("enabled", False))
            
            # Load camera position
            self.camera_x_var.set(camera["position"][0])
            self.camera_y_var.set(camera["position"][1])
            self.camera_z_var.set(camera["position"][2])
            
            # Load field of view
            self.camera_fov_var.set(camera.get("fov", 70.0))
            
            self.update_visualization_preview()
            
    def set_camera_preset(self, position):
        """Set a camera position preset"""
        self.camera_x_var.set(position[0])
        self.camera_y_var.set(position[1])
        self.camera_z_var.set(position[2])
        
        # Update the preview with new camera position
        self.update_visualization_preview()
        
        # Log the preset selection
        preset_name = "Custom"
        if position == [0.0, 0.0, 5.0]:
            preset_name = "Above (Best for Targets)"
        elif position == [0.0, 0.0, 10.0]:
            preset_name = "High Above"
        elif position == [5.0, 0.0, 0.0]:
            preset_name = "Side View"
        elif position == [0.0, 0.0, 15.0]:
            preset_name = "Earth View"
            
        self.parent_app.add_log(f"Set camera preset: {preset_name} at {position}")
        
    def update_active_camera_indicator(self):
        """Update the active camera indicator"""
        active_sat = self.get_active_camera_satellite()
        if active_sat:
            self.active_camera_label.config(text=f"{active_sat['name']}", foreground="green")
        else:
            self.active_camera_label.config(text="None", foreground="red")
            
    def apply_camera_settings(self):
        """Apply camera settings to the selected satellite"""
        sat_name = self.viz_satellite_var.get()
        sat_index = -1
        
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                sat_index = i
                break
                
        if sat_index >= 0:
            # Check if we're enabling a camera
            enabling_camera = (not self.satellites[sat_index]["camera"].get("enabled", False)) and self.camera_enabled_var.get()
            
            if enabling_camera and self.get_active_camera_satellite():
                # If another satellite already has an enabled camera, confirm with user
                confirm = messagebox.askyesno("Camera Already Active", 
                                           f"Another satellite already has an active camera. Disable that camera and enable this one?")
                if confirm:
                    # Disable all other cameras
                    for sat in self.satellites:
                        if sat["name"] != sat_name:
                            sat["camera"]["enabled"] = False
                else:
                    # User cancelled, don't enable this camera
                    self.camera_enabled_var.set(False)
            
            # Update camera settings
            camera = self.satellites[sat_index]["camera"]
            
            # Set enabled state
            camera["enabled"] = self.camera_enabled_var.get()
            
            # Update position
            camera["position"] = [
                self.camera_x_var.get(),
                self.camera_y_var.get(),
                self.camera_z_var.get()
            ]
            
            # Update field of view
            camera["fov"] = self.camera_fov_var.get()
            
            # Update active camera indicator
            self.update_active_camera_indicator()
            
            # Update visualization preview
            self.update_visualization_preview()
            
            # Update constellation tab to show camera indicator
            try:
                self.parent_app.constellation_tab.update_satellite_listbox()
            except:
                pass
            
            self.parent_app.add_log(f"Applied camera settings to {sat_name} - Position: {camera['position']}")
            
    def reset_3d_view(self):
        """Reset the 3D view to default angle"""
        self.viz_subplot.view_init(elev=20, azim=45)
        self.viz_canvas.draw()
    
    def set_top_view(self):
        """Set 3D view to top-down"""
        self.viz_subplot.view_init(elev=90, azim=0)
        self.viz_canvas.draw()
    
    def set_side_view(self):
        """Set 3D view to side view"""
        self.viz_subplot.view_init(elev=0, azim=0)
        self.viz_canvas.draw()
            
    def update_visualization_preview(self):
        """FIXED: Update the visualization preview with enhanced graphics and better target visibility"""
        sat_name = self.viz_satellite_var.get()
        if not sat_name:
            return
            
        # Clear the previous plot
        self.viz_subplot.clear()
        
        # FIXED: Set up the plot with enhanced styling and proper background
        self.viz_subplot.set_facecolor('lightblue')  # Sky-like background
        self.viz_subplot.set_xlabel('X (m)', fontsize=12, color='darkblue')
        self.viz_subplot.set_ylabel('Y (m)', fontsize=12, color='darkblue')
        self.viz_subplot.set_zlabel('Z (m)', fontsize=12, color='darkblue')
        self.viz_subplot.set_title(f'Enhanced Preview: {sat_name} - Camera & Target View', 
                                  fontsize=14, fontweight='bold', color='darkgreen')
        
        # FIXED: Draw a more detailed spacecraft body with enhanced visibility
        self.viz_subplot.scatter([0], [0], [0], color='black', s=200, marker='s', 
                                label='Spacecraft Body', zorder=10, edgecolor='white', linewidth=2)
        
        # FIXED: Draw the reaction wheels with enhanced visibility
        wheel_colors = ['red', 'green', 'blue', 'orange']
        wheel_labels = ['Wheel 0', 'Wheel 1', 'Wheel 2', 'Wheel 3']
        wheel_positions = [
            [0.5, 0.5, 0], [0.5, -0.5, 0], 
            [-0.5, 0.5, 0], [-0.5, -0.5, 0]
        ]
        
        for i in range(4):
            pos = wheel_positions[i]
            self.viz_subplot.scatter([pos[0]], [pos[1]], [pos[2]], 
                                   color=wheel_colors[i], s=100, marker='o',
                                   label=wheel_labels[i], zorder=8, 
                                   edgecolor='black', linewidth=1)
            
            # Draw connection lines to show wheel mounting
            self.viz_subplot.plot([0, pos[0]], [0, pos[1]], [0, pos[2]], 
                                 color=wheel_colors[i], alpha=0.6, linewidth=2, zorder=6)
        
        # FIXED: Draw camera position if enabled with enhanced visualization
        if self.camera_enabled_var.get():
            camera_pos = [self.camera_x_var.get(), self.camera_y_var.get(), self.camera_z_var.get()]
            self.viz_subplot.scatter([camera_pos[0]], [camera_pos[1]], [camera_pos[2]], 
                                   color='cyan', s=150, marker='^', 
                                   label='Camera (ENHANCED)', zorder=10,
                                   edgecolor='darkblue', linewidth=2)
            
            # Draw a line from origin to camera
            self.viz_subplot.plot([0, camera_pos[0]], [0, camera_pos[1]], [0, camera_pos[2]], 
                                 color='cyan', linestyle='--', linewidth=3, alpha=0.8, zorder=7)
            
            # FIXED: Draw enhanced field of view cone for target visibility
            fov_rad = np.radians(self.camera_fov_var.get() / 2)
            dist = np.sqrt(camera_pos[0]**2 + camera_pos[1]**2 + camera_pos[2]**2)
            
            # Normalize camera position vector
            norm = np.sqrt(sum(p**2 for p in camera_pos))
            if norm > 0:
                camera_dir = [p/norm for p in camera_pos]
                
                # Create an enhanced cone to represent field of view
                cone_radius = max(dist * np.tan(fov_rad), 1.0)  # Larger minimum visible size
                cone_tip = camera_pos
                
                # FIXED: Draw more lines to represent the cone for enhanced target visibility
                num_lines = 12  # More lines for better visibility
                for i in range(num_lines):
                    angle = i * (2 * np.pi / num_lines)
                    # Calculate a point on the cone base
                    dx = cone_radius * np.cos(angle)
                    dy = cone_radius * np.sin(angle)
                    
                    # Project perpendicular to camera direction (enhanced)
                    if abs(camera_dir[2]) < 0.9:  # Not pointing straight up/down
                        base_x = cone_tip[0] - camera_dir[0] * dist * 0.6 + dx * 0.6
                        base_y = cone_tip[1] - camera_dir[1] * dist * 0.6 + dy * 0.6
                        base_z = cone_tip[2] - camera_dir[2] * dist * 0.6
                    else:  # Pointing straight up/down
                        base_x = cone_tip[0] + dx * 0.6
                        base_y = cone_tip[1] + dy * 0.6
                        base_z = cone_tip[2] - camera_dir[2] * dist * 0.6
                    
                    self.viz_subplot.plot([cone_tip[0], base_x], 
                                         [cone_tip[1], base_y],
                                         [cone_tip[2], base_z], 
                                         color='cyan', linestyle=':', alpha=0.7, linewidth=1.5, zorder=5)
        
        # FIXED: Add targets if any are assigned to this satellite with ENHANCED visibility
        target_count = 0
        earth_targets = []
        for target in self.parent_app.targets:
            if sat_name in target["assigned_to"]:
                # FIXED: Position targets to simulate Earth surface visibility
                angle = target_count * (2 * np.pi / max(len(target["assigned_to"]), 1))
                # Simulate Earth surface at larger distance for visibility
                earth_radius = 8.0  # Simulated Earth surface radius for preview
                target_x = earth_radius * np.cos(angle) * np.cos(target_count * 0.3)
                target_y = earth_radius * np.sin(angle) * np.cos(target_count * 0.3)
                target_z = -earth_radius * 0.8 + target_count * 0.5  # Distributed on Earth surface
                
                self.viz_subplot.scatter([target_x], [target_y], [target_z], 
                                       color=target["color"], s=200, marker='*',
                                       label=f'{target["name"]} (TARGET)', zorder=9,
                                       edgecolor='black', linewidth=1)
                
                # Draw line to show visibility from camera
                if self.camera_enabled_var.get():
                    camera_pos = [self.camera_x_var.get(), self.camera_y_var.get(), self.camera_z_var.get()]
                    self.viz_subplot.plot([camera_pos[0], target_x], 
                                         [camera_pos[1], target_y], 
                                         [camera_pos[2], target_z], 
                                         color=target["color"], alpha=0.4, linewidth=2, 
                                         linestyle=':', zorder=4)
                
                earth_targets.append((target_x, target_y, target_z))
                target_count += 1
        
        # FIXED: Draw a simplified Earth representation
        if earth_targets:
            # Draw Earth as a sphere mesh for context
            u = np.linspace(0, 2 * np.pi, 10)
            v = np.linspace(0, np.pi, 10)
            earth_radius = 7.0
            x_earth = earth_radius * np.outer(np.cos(u), np.sin(v))
            y_earth = earth_radius * np.outer(np.sin(u), np.sin(v))
            z_earth = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v)) - 8.0
            
            self.viz_subplot.plot_surface(x_earth, y_earth, z_earth, alpha=0.3, color='lightgreen', zorder=1)
        
        # FIXED: Add enhanced coordinate axes for reference
        axis_length = 2.0
        self.viz_subplot.plot([0, axis_length], [0, 0], [0, 0], 'r-', alpha=0.8, linewidth=4, zorder=3)
        self.viz_subplot.plot([0, 0], [0, axis_length], [0, 0], 'g-', alpha=0.8, linewidth=4, zorder=3)
        self.viz_subplot.plot([0, 0], [0, 0], [0, axis_length], 'b-', alpha=0.8, linewidth=4, zorder=3)
        
        # Add enhanced axis labels
        self.viz_subplot.text(axis_length*1.3, 0, 0, 'X', color='red', fontsize=14, fontweight='bold')
        self.viz_subplot.text(0, axis_length*1.3, 0, 'Y', color='green', fontsize=14, fontweight='bold')
        self.viz_subplot.text(0, 0, axis_length*1.3, 'Z', color='blue', fontsize=14, fontweight='bold')
        
        # FIXED: Add enhanced legend with better positioning
        self.viz_subplot.legend(loc='upper left', bbox_to_anchor=(0, 1), fontsize=9, framealpha=0.9)
        
        # FIXED: Set enhanced aspect ratio and appropriate limits for target visibility
        camera_pos = [self.camera_x_var.get(), self.camera_y_var.get(), self.camera_z_var.get()]
        max_coord = max(10.0, max(abs(val) for val in camera_pos) * 1.5)
        
        self.viz_subplot.set_xlim(-max_coord, max_coord)
        self.viz_subplot.set_ylim(-max_coord, max_coord)
        self.viz_subplot.set_zlim(-max_coord, max_coord)
        
        # Set a good viewing angle for enhanced visibility
        self.viz_subplot.view_init(elev=25, azim=45)
        
        # FIXED: Update the canvas with tight layout
        try:
            self.viz_figure.tight_layout(pad=2.0)
        except:
            pass  # Ignore layout warnings
        self.viz_canvas.draw()