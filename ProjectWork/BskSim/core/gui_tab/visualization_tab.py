#!/usr/bin/env python
"""
visualization_tab.py

Visualization settings and preview tab for the Spacecraft Constellation Fault Simulator.
Provides camera configuration, orbit display controls, and 3D constellation preview.
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set backend before importing pyplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from .base_tab import BaseTab


class VisualizationTab(BaseTab):
    """Tab for visualization settings and constellation preview"""
    
    def __init__(self, parent_app, parent_frame):
        """Initialize the visualization tab"""
        super().__init__(parent_app, parent_frame)
        
        # References
        self.satellites = parent_app.satellites
        self.targets = parent_app.targets
        
        # Visualization settings
        self.orbit_visibility = {}
        self.satellite_colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 'orange', 'purple']
        
        # Camera settings
        self.camera_enabled_var = tk.BooleanVar(value=False)
        self.camera_x_var = tk.DoubleVar(value=0.0)
        self.camera_y_var = tk.DoubleVar(value=0.0)
        self.camera_z_var = tk.DoubleVar(value=15.0)
        self.camera_fov_var = tk.DoubleVar(value=80.0)
        
        # Display settings
        self.show_satellite_labels_var = tk.BooleanVar(value=True)
        self.show_assignment_lines_var = tk.BooleanVar(value=True)
        self.satellite_size_var = tk.DoubleVar(value=8.0)
        self.orbit_line_width_var = tk.DoubleVar(value=2.0)
        self.orbit_transparency_var = tk.DoubleVar(value=0.8)
        self.target_altitude_var = tk.DoubleVar(value=50.0)
        self.target_marker_size_var = tk.DoubleVar(value=2000.0)
        self.color_scheme_var = tk.StringVar(value="distinct")
        
        # Satellite selection
        self.viz_satellite_var = tk.StringVar()
        
        # Create UI
        self.create_ui()
        
    def create_ui(self):
        """Create the main UI layout"""
        # Main container
        main_container = ttk.Frame(self.parent_frame, padding=5)
        main_container.pack(fill=tk.BOTH, expand=True)

        # Create notebook for organization (settings tabs go here)
        self.notebook = ttk.Notebook(main_container)
        self.notebook.pack(fill=tk.BOTH, expand=True, pady=(0, 10))  # Add padding below if needed

        # Create setting tabs
        self._create_camera_tab()
        self._create_orbit_tab()
        self._create_display_tab()
        self._create_target_tab()
        self._create_preview_panel_tab()

        # Initial update
        self.update_satellite_dropdown()
        self.update_preview()

    def _create_preview_panel_tab(self):
        """Add preview as a separate tab in the notebook"""
        preview_frame = ttk.Frame(self.notebook)
        self.notebook.add(preview_frame, text="Visualization Preview")
        self._create_preview_panel(preview_frame)

        
    def _create_camera_tab(self):
        """Create camera settings tab"""
        camera_frame = ttk.Frame(self.notebook)
        self.notebook.add(camera_frame, text="Camera & View")
        
        # Satellite selection
        select_frame = ttk.Frame(camera_frame, padding=10)
        select_frame.pack(fill=tk.X)
        
        ttk.Label(select_frame, text="Satellite for Camera:").pack(side=tk.LEFT, padx=(0, 10))
        self.viz_satellite_combo = ttk.Combobox(select_frame, textvariable=self.viz_satellite_var, 
                                              state="readonly", width=20)
        self.viz_satellite_combo.pack(side=tk.LEFT)
        self.viz_satellite_combo.bind('<<ComboboxSelected>>', self._on_satellite_changed)
        
        # Active camera indicator
        self.active_camera_label = ttk.Label(select_frame, text="No active camera", 
                                           foreground="red", font=('Segoe UI', 10, 'bold'))
        self.active_camera_label.pack(side=tk.RIGHT, padx=10)
        
        # Camera settings
        settings_frame = ttk.LabelFrame(camera_frame, text="Camera Settings", padding=10)
        settings_frame.pack(fill=tk.X, padx=10, pady=10)
        
        # Enable checkbox
        ttk.Checkbutton(settings_frame, text="Enable Camera for Vizard", 
                       variable=self.camera_enabled_var).pack(anchor=tk.W)
        
        # Position controls
        pos_frame = ttk.Frame(settings_frame)
        pos_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(pos_frame, text="Position (Body Frame):").grid(row=0, column=0, sticky=tk.W, pady=5)
        
        # X, Y, Z entries
        coords = [("X:", self.camera_x_var), ("Y:", self.camera_y_var), ("Z:", self.camera_z_var)]
        for i, (label, var) in enumerate(coords):
            ttk.Label(pos_frame, text=label).grid(row=1, column=i*2, padx=(10, 5))
            ttk.Entry(pos_frame, textvariable=var, width=10).grid(row=1, column=i*2+1, padx=(0, 10))
        
        # FOV control
        ttk.Label(pos_frame, text="FOV (deg):").grid(row=1, column=6, padx=(20, 5))
        ttk.Entry(pos_frame, textvariable=self.camera_fov_var, width=10).grid(row=1, column=7)
        
        # Preset buttons
        preset_frame = ttk.Frame(settings_frame)
        preset_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(preset_frame, text="View Presets:").pack(side=tk.LEFT, padx=(0, 10))
        
        presets = [
            ("Earth View", [0.0, 0.0, 15.0]),
            ("Target View", [0.0, 0.0, 12.0]),
            ("Side View", [8.0, 0.0, 0.0]),
            ("Front View", [0.0, 8.0, 0.0]),
            ("Close", [0.0, 0.0, 8.0]),
            ("Far", [0.0, 0.0, 20.0])
        ]
        
        for name, pos in presets:
            ttk.Button(preset_frame, text=name, 
                      command=lambda p=pos: self._set_camera_preset(p)).pack(side=tk.LEFT, padx=2)
        
        # Action buttons
        action_frame = ttk.Frame(settings_frame)
        action_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(action_frame, text="Apply Settings", 
                  command=self._apply_camera_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Auto Setup", 
                  command=self._auto_setup_camera).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Disable All", 
                  command=self._disable_all_cameras).pack(side=tk.RIGHT, padx=5)
        
        
        
    def _create_orbit_tab(self):
        """Create orbit control tab"""
        orbit_frame = ttk.Frame(self.notebook)
        self.notebook.add(orbit_frame, text="Orbit Control")
        
        # Master controls
        master_frame = ttk.LabelFrame(orbit_frame, text="Master Controls", padding=10)
        master_frame.pack(fill=tk.X, padx=10, pady=10)
        
        control_row = ttk.Frame(master_frame)
        control_row.pack(fill=tk.X)
        
        ttk.Button(control_row, text="Show All", 
                  command=self._show_all_orbits).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_row, text="Hide All", 
                  command=self._hide_all_orbits).pack(side=tk.LEFT, padx=5)
        
        # Orbit line settings
        settings_row = ttk.Frame(master_frame)
        settings_row.pack(fill=tk.X, pady=10)
        
        ttk.Label(settings_row, text="Line Width:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Spinbox(settings_row, from_=1.0, to=5.0, increment=0.5,
                   textvariable=self.orbit_line_width_var, width=8).pack(side=tk.LEFT, padx=(0, 20))
        
        ttk.Label(settings_row, text="Opacity:").pack(side=tk.LEFT, padx=(0, 5))
        ttk.Spinbox(settings_row, from_=0.1, to=1.0, increment=0.1,
                   textvariable=self.orbit_transparency_var, width=8).pack(side=tk.LEFT)
        
        # Individual orbit controls
        orbit_list_frame = ttk.LabelFrame(orbit_frame, text="Individual Orbits", padding=10)
        orbit_list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Scrollable frame
        canvas = tk.Canvas(orbit_list_frame, height=200)
        scrollbar = ttk.Scrollbar(orbit_list_frame, orient="vertical", command=canvas.yview)
        self.orbit_controls_frame = ttk.Frame(canvas)
        
        self.orbit_controls_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.orbit_controls_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
    def _create_display_tab(self):
        """Create satellite display settings tab"""
        display_frame = ttk.Frame(self.notebook)
        self.notebook.add(display_frame, text="Satellite Display")
        
        # Display settings
        settings_frame = ttk.LabelFrame(display_frame, text="Display Settings", padding=10)
        settings_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Checkbutton(settings_frame, text="Show satellite names", 
                       variable=self.show_satellite_labels_var).pack(anchor=tk.W, pady=5)
        
        size_frame = ttk.Frame(settings_frame)
        size_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(size_frame, text="Satellite Size:").pack(side=tk.LEFT, padx=(0, 10))
        ttk.Spinbox(size_frame, from_=1.0, to=15.0, increment=0.5,
                   textvariable=self.satellite_size_var, width=10).pack(side=tk.LEFT)
        ttk.Label(size_frame, text="(multiplier for visibility)").pack(side=tk.LEFT, padx=10)
        
        # Color scheme
        color_frame = ttk.LabelFrame(display_frame, text="Color Scheme", padding=10)
        color_frame.pack(fill=tk.X, padx=10, pady=10)
        
        schemes = [("Distinct", "distinct"), ("Rainbow", "rainbow"), 
                  ("Cool", "cool"), ("Warm", "warm")]
        
        for text, value in schemes:
            ttk.Radiobutton(color_frame, text=text, value=value,
                          variable=self.color_scheme_var,
                          command=self._update_color_scheme).pack(anchor=tk.W, pady=2)
        
        # Satellite list
        list_frame = ttk.LabelFrame(display_frame, text="Satellite Information", padding=10)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        columns = ('Name', 'Color', 'Camera', 'Fault', 'Targets', 'Altitude')
        self.satellite_tree = ttk.Treeview(list_frame, columns=columns, show='headings', height=8)
        
        for col in columns:
            self.satellite_tree.heading(col, text=col)
            self.satellite_tree.column(col, width=100)
        
        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.satellite_tree.yview)
        self.satellite_tree.configure(yscrollcommand=scrollbar.set)
        
        self.satellite_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
    def _create_target_tab(self):
        """Create target display settings tab"""
        target_frame = ttk.Frame(self.notebook)
        self.notebook.add(target_frame, text="Target Settings")
        
        # Target visualization settings
        settings_frame = ttk.LabelFrame(target_frame, text="Target Display", padding=10)
        settings_frame.pack(fill=tk.X, padx=10, pady=10)
        
        alt_frame = ttk.Frame(settings_frame)
        alt_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(alt_frame, text="Altitude (km):").pack(side=tk.LEFT, padx=(0, 10))
        ttk.Spinbox(alt_frame, from_=10.0, to=100.0, increment=10.0,
                   textvariable=self.target_altitude_var, width=10).pack(side=tk.LEFT)
        ttk.Label(alt_frame, text="(above Earth surface)").pack(side=tk.LEFT, padx=10)
        
        range_frame = ttk.Frame(settings_frame)
        range_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(range_frame, text="Marker Range (km):").pack(side=tk.LEFT, padx=(0, 10))
        ttk.Spinbox(range_frame, from_=500.0, to=5000.0, increment=500.0,
                   textvariable=self.target_marker_size_var, width=10).pack(side=tk.LEFT)
        
        ttk.Checkbutton(settings_frame, text="Show assignment connections",
                       variable=self.show_assignment_lines_var).pack(anchor=tk.W, pady=10)
        
        # Assignment summary
        summary_frame = ttk.LabelFrame(target_frame, text="Target Assignments", padding=10)
        summary_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.assignment_text = tk.Text(summary_frame, height=10, wrap=tk.WORD, state="disabled")
        scrollbar = ttk.Scrollbar(summary_frame, command=self.assignment_text.yview)
        self.assignment_text.config(yscrollcommand=scrollbar.set)
        
        self.assignment_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
    def _create_preview_panel(self, main_container):
        """Create the 3D preview panel with tabs"""
        preview_frame = ttk.LabelFrame(main_container, text="Visualization Preview", padding=10)

        preview_frame.pack(fill=tk.BOTH, expand=True, pady=10)

        # Create preview tab notebook
        preview_notebook = ttk.Notebook(preview_frame)
        preview_notebook.pack(fill=tk.BOTH, expand=True)

        # ----- Constellation Preview Tab -----
        constellation_frame = ttk.Frame(preview_notebook)
        preview_notebook.add(constellation_frame, text="Constellation Preview")

        # Create matplotlib figure
        self.figure = Figure(figsize=(12, 4), dpi=100)
        self.ax = self.figure.add_subplot(111, projection='3d')

        self.canvas = FigureCanvasTkAgg(self.figure, constellation_frame)
        

        # Pack canvas first
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Pack toolbar directly (full row, not inside another nested frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, constellation_frame)
        self.toolbar.update()

        # Add custom view buttons below the toolbar
        button_row = ttk.Frame(constellation_frame)
        button_row.pack(fill=tk.X, pady=5)

        ttk.Button(button_row, text="Reset View", command=self._reset_view).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_row, text="Constellation View", command=self._constellation_view).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_row, text="Earth View", command=self._earth_view).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_row, text="Update", command=self.update_preview).pack(side=tk.RIGHT, padx=10)


        # --- Satellite Preview Tab ---
        satellite_frame = ttk.Frame(preview_notebook)
        preview_notebook.add(satellite_frame, text="Satellite Preview")

        # Dropdown to choose satellite
        select_frame = ttk.Frame(satellite_frame)
        select_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(select_frame, text="Select Satellite:").pack(side=tk.LEFT)
        self.sat_preview_var = tk.StringVar()
        self.sat_preview_combo = ttk.Combobox(
            select_frame,
            textvariable=self.sat_preview_var,
            values=[s["name"] for s in self.satellites],
            state="readonly",
            width=20
        )
        self.sat_preview_combo.pack(side=tk.LEFT, padx=10)
        self.sat_preview_combo.bind("<<ComboboxSelected>>", self._update_satellite_preview)

        #  3D Figure
        self.sat_fig = Figure(figsize=(6, 5), dpi=100)
        self.sat_ax = self.sat_fig.add_subplot(111, projection='3d')
        self.sat_canvas = FigureCanvasTkAgg(self.sat_fig, master=satellite_frame)
        self.sat_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Optional: Add toolbar
        toolbar = NavigationToolbar2Tk(self.sat_canvas, satellite_frame)
        toolbar.update()



    def _update_satellite_preview(self, event=None):
        """Update satellite-specific preview in Satellite Preview tab."""
        sat_name = self.sat_preview_var.get()
        sat = next((s for s in self.satellites if s["name"] == sat_name), None)
        
        if not sat:
            return

        self.sat_ax.clear()
        self.sat_ax.set_title(f"Preview: {sat_name}", fontsize=12, fontweight='bold')
        self.sat_ax.set_xlabel("X (m)")
        self.sat_ax.set_ylabel("Y (m)")
        self.sat_ax.set_zlabel("Z (m)")

        # Spacecraft Body
        self.sat_ax.scatter(0, 0, 0, color='black', s=150, label="Spacecraft Body", zorder=3)

        # Reaction wheels
        wheel_colors = ['red', 'green', 'blue', 'orange']
        wheel_offsets = [[0.5, 0.5, 0], [0.5, -0.5, 0], [-0.5, 0.5, 0], [-0.5, -0.5, 0]]
        for i in range(4):
            x, y, z = wheel_offsets[i]
            self.sat_ax.scatter(x, y, z, color=wheel_colors[i], s=100, label=f"Wheel {i}")
            self.sat_ax.plot([0, x], [0, y], [0, z], color=wheel_colors[i], alpha=0.5)

        # Camera
        if sat["camera"]["enabled"]:
            cam_pos = sat["camera"]["position"]
            self.sat_ax.scatter(*cam_pos, color='cyan', s=120, marker='^', label='Camera')
            self.sat_ax.plot([0, cam_pos[0]], [0, cam_pos[1]], [0, cam_pos[2]], 'c--')

        # Targets
        for i, target in enumerate(self.parent_app.targets):
            if sat["name"] in target.get("assigned_to", []):
                angle = i * (2 * np.pi / max(len(self.parent_app.targets), 1))
                x = 3.0 * np.cos(angle)
                y = 3.0 * np.sin(angle)
                z = -2.0 - 0.2 * i
                self.sat_ax.scatter(x, y, z, s=180, marker='*', color=target["color"], label=target["name"])
                self.sat_ax.plot([0, x], [0, y], [0, z], linestyle=':', color=target["color"], alpha=0.6)

        # Axes
        self.sat_ax.plot([0, 2], [0, 0], [0, 0], 'r-', label='X')
        self.sat_ax.plot([0, 0], [0, 2], [0, 0], 'g-', label='Y')
        self.sat_ax.plot([0, 0], [0, 0], [0, 2], 'b-', label='Z')

        self.sat_ax.legend(loc='upper left', fontsize=8)
        self.sat_ax.view_init(elev=25, azim=45)
        self.sat_ax.set_xlim(-4, 4)
        self.sat_ax.set_ylim(-4, 4)
        self.sat_ax.set_zlim(-4, 4)

        self.sat_canvas.draw()


        
    def update_preview(self):
        """Update the 3D constellation preview"""
        try:
            self.ax.clear()
            
            # Set labels and title
            self.ax.set_xlabel('X (km)', fontsize=10)
            self.ax.set_ylabel('Y (km)', fontsize=10)
            self.ax.set_zlabel('Z (km)', fontsize=10)
            
            active_sat = self.viz_satellite_var.get() or "None"
            self.ax.set_title(f'Constellation Preview - Camera: {active_sat}', fontsize=12)
            
            if not self.satellites:
                self.ax.text(0, 0, 0, 'No satellites configured\nAdd satellites in Constellation tab',
                           ha='center', va='center', fontsize=12)
                self.ax.set_xlim(-1, 1)
                self.ax.set_ylim(-1, 1)
                self.ax.set_zlim(-1, 1)
                self.canvas.draw()
                return
            
            max_coord = 10.0
            
            # Draw satellites
            for i, sat in enumerate(self.satellites):
                color = self.satellite_colors[i % len(self.satellite_colors)]
                
                # Position satellites in a pattern
                if i == 0:
                    pos = [0, 0, 0]
                else:
                    angle = 2 * np.pi * i / len(self.satellites)
                    radius = 5.0
                    pos = [radius * np.cos(angle), radius * np.sin(angle), 2 * np.sin(angle * 2)]
                
                # Draw satellite
                size = 200 if sat["name"] == active_sat else 150
                self.ax.scatter(*pos, color=color, s=size, marker='o', 
                              edgecolor='black', linewidth=2, alpha=0.9, label=sat["name"])
                
                # Label
                if self.show_satellite_labels_var.get():
                    self.ax.text(pos[0], pos[1], pos[2] + 1.5, sat["name"],
                               fontsize=9, fontweight='bold', color=color, ha='center')
                
                # Camera indicator
                if sat["camera"]["enabled"]:
                    cam_pos = [pos[j] + sat["camera"]["position"][j] for j in range(3)]
                    self.ax.scatter(*cam_pos, color='cyan', s=100, marker='^', 
                                  edgecolor='black', linewidth=2)
                    self.ax.plot([pos[0], cam_pos[0]], [pos[1], cam_pos[1]], 
                               [pos[2], cam_pos[2]], 'c--', linewidth=2, alpha=0.8)
                
                # Orbit line
                if self.orbit_visibility.get(sat["name"], True):
                    theta = np.linspace(0, 2*np.pi, 50)
                    orbit_r = 8.0 + i * 0.5
                    orbit_x = orbit_r * np.cos(theta) + pos[0]
                    orbit_y = orbit_r * np.sin(theta) + pos[1]
                    orbit_z = np.zeros_like(orbit_x) + pos[2]
                    
                    self.ax.plot(orbit_x, orbit_y, orbit_z, color=color,
                               linewidth=self.orbit_line_width_var.get(),
                               alpha=self.orbit_transparency_var.get())
            
            # Draw coordinate axes
            axis_len = max_coord * 0.3
            self.ax.plot([0, axis_len], [0, 0], [0, 0], 'r-', linewidth=3, alpha=0.8)
            self.ax.plot([0, 0], [0, axis_len], [0, 0], 'g-', linewidth=3, alpha=0.8)
            self.ax.plot([0, 0], [0, 0], [0, axis_len], 'b-', linewidth=3, alpha=0.8)
            
            self.ax.text(axis_len*1.1, 0, 0, 'X', color='red', fontsize=10, fontweight='bold')
            self.ax.text(0, axis_len*1.1, 0, 'Y', color='green', fontsize=10, fontweight='bold')
            self.ax.text(0, 0, axis_len*1.1, 'Z', color='blue', fontsize=10, fontweight='bold')
            
            # Draw targets
            assigned_targets = [t for t in self.parent_app.targets if t.get("assigned_to", [])]
            for target in assigned_targets:
                # Simple target visualization
                self.ax.scatter(0, 0, -3, color=target["color"], s=300, marker='*',
                              edgecolor='black', linewidth=2, alpha=0.9)
                self.ax.text(0, 0, -4, target["name"], fontsize=8, ha='center',
                           color=target["color"], fontweight='bold')
            
            # Set limits and legend
            self.ax.set_xlim(-max_coord, max_coord)
            self.ax.set_ylim(-max_coord, max_coord)
            self.ax.set_zlim(-max_coord, max_coord)
            
            self.ax.legend(loc='upper left', fontsize=8)
            self.ax.view_init(elev=20, azim=45)
            
            self.figure.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            self.add_log(f"Error updating visualization preview: {e}")
            self.ax.clear()
            self.ax.text(0, 0, 0, f'Error: {str(e)[:50]}...\nCheck console',
                       ha='center', va='center', fontsize=10)
            self.ax.set_xlim(-1, 1)
            self.ax.set_ylim(-1, 1)
            self.ax.set_zlim(-1, 1)
            self.canvas.draw()
    
    # Camera methods
    def get_camera_position(self):
        """Get current camera position"""
        return [self.camera_x_var.get(), self.camera_y_var.get(), self.camera_z_var.get()]
    
    def get_active_camera_satellite(self):
        """Get the satellite with active camera"""
        for sat in self.satellites:
            if sat["camera"]["enabled"]:
                return sat
        return None
    
    def update_satellite_dropdown(self):
        """Update satellite dropdown list"""
        sat_names = [sat['name'] for sat in self.satellites]
        self.viz_satellite_combo['values'] = sat_names
        
        if sat_names and not self.viz_satellite_var.get():
            self.viz_satellite_combo.current(0)
            self._load_camera_settings(0)
            
        self._update_orbit_controls()
        self._update_satellite_list()
        self._update_active_camera_indicator()
        self._update_assignment_summary()
    
    # Private methods
    def _on_satellite_changed(self, event):
        """Handle satellite selection change"""
        sat_name = self.viz_satellite_var.get()
        for i, sat in enumerate(self.satellites):
            if sat["name"] == sat_name:
                self._load_camera_settings(i)
                break
    
    def _load_camera_settings(self, index):
        """Load camera settings for satellite"""
        if 0 <= index < len(self.satellites):
            camera = self.satellites[index]["camera"]
            self.camera_enabled_var.set(camera.get("enabled", False))
            self.camera_x_var.set(camera["position"][0])
            self.camera_y_var.set(camera["position"][1])
            self.camera_z_var.set(camera["position"][2])
            self.camera_fov_var.set(camera.get("fov", 80.0))
    
    def _set_camera_preset(self, position):
        """Set camera to preset position"""
        self.camera_x_var.set(position[0])
        self.camera_y_var.set(position[1])
        self.camera_z_var.set(position[2])
        self.add_log(f"Camera preset applied: {position}")
    
    def _apply_camera_settings(self):
        """Apply current camera settings"""
        sat_name = self.viz_satellite_var.get()
        for sat in self.satellites:
            if sat["name"] == sat_name:
                # Disable other cameras if enabling this one
                if self.camera_enabled_var.get() and not sat["camera"]["enabled"]:
                    for other_sat in self.satellites:
                        if other_sat["name"] != sat_name:
                            other_sat["camera"]["enabled"] = False
                
                # Apply settings
                sat["camera"]["enabled"] = self.camera_enabled_var.get()
                sat["camera"]["position"] = self.get_camera_position()
                sat["camera"]["fov"] = self.camera_fov_var.get()
                
                self._update_active_camera_indicator()
                self._update_satellite_list()
                self.update_preview()
                
                status = "enabled" if sat["camera"]["enabled"] else "disabled"
                self.add_log(f"Camera {status} for {sat_name}")
                break
    
    def _auto_setup_camera(self):
        """Automatically setup camera for optimal viewing"""
        if not self.satellites:
            messagebox.showinfo("No Satellites", "Add satellites first")
            return
            
        # Disable all cameras
        for sat in self.satellites:
            sat["camera"]["enabled"] = False
            
        # Enable first satellite with optimal settings
        self.satellites[0]["camera"]["enabled"] = True
        self.satellites[0]["camera"]["position"] = [0.0, 0.0, 15.0]
        self.satellites[0]["camera"]["fov"] = 80.0
        
        # Update UI
        self.viz_satellite_combo.current(0)
        self._load_camera_settings(0)
        self._update_active_camera_indicator()
        self._update_satellite_list()
        self.update_preview()
        
        self.add_log("Auto camera setup complete")
        messagebox.showinfo("Camera Setup", 
                          f"Camera enabled for {self.satellites[0]['name']}\n"
                          "Position: Earth View (0, 0, 15)")
    
    def _disable_all_cameras(self):
        """Disable all satellite cameras"""
        for sat in self.satellites:
            sat["camera"]["enabled"] = False
            
        self.camera_enabled_var.set(False)
        self._update_active_camera_indicator()
        self._update_satellite_list()
        self.update_preview()
        
        self.add_log("All cameras disabled")
    
    def _update_active_camera_indicator(self):
        """Update the active camera label"""
        active_sat = self.get_active_camera_satellite()
        if active_sat:
            self.active_camera_label.config(text=active_sat['name'], foreground="green")
        else:
            self.active_camera_label.config(text="No active camera", foreground="red")
    
    # Orbit control methods
    def _show_all_orbits(self):
        """Show all satellite orbits"""
        for var in self.orbit_visibility.values():
            if isinstance(var, tk.BooleanVar):
                var.set(True)
        self.update_preview()
        self.add_log("All orbits visible")
    
    def _hide_all_orbits(self):
        """Hide all satellite orbits"""
        for var in self.orbit_visibility.values():
            if isinstance(var, tk.BooleanVar):
                var.set(False)
        self.update_preview()
        self.add_log("All orbits hidden")
    
    def _update_orbit_controls(self):
        """Update individual orbit visibility controls"""
        # Clear existing
        for widget in self.orbit_controls_frame.winfo_children():
            widget.destroy()
        
        # Create controls for each satellite
        for i, sat in enumerate(self.satellites):
            frame = ttk.Frame(self.orbit_controls_frame)
            frame.pack(fill=tk.X, pady=2, padx=5)
            
            # Color indicator
            color = self.satellite_colors[i % len(self.satellite_colors)]
            canvas = tk.Canvas(frame, width=20, height=20, bg=color, highlightthickness=1)
            canvas.pack(side=tk.LEFT, padx=(0, 10))
            
            # Name
            ttk.Label(frame, text=sat['name'], font=('Segoe UI', 10, 'bold')).pack(side=tk.LEFT, padx=(0, 20))
            
            # Visibility checkbox
            if sat['name'] not in self.orbit_visibility:
                self.orbit_visibility[sat['name']] = tk.BooleanVar(value=True)
                
            ttk.Checkbutton(frame, text="Show Orbit",
                          variable=self.orbit_visibility[sat['name']],
                          command=self.update_preview).pack(side=tk.LEFT)
            
            # Orbit info
            try:
                alt = sat['orbit']['a'] - 6371
                period = 2 * np.pi * np.sqrt((sat['orbit']['a'] * 1000)**3 / 3.986004418e14) / 60
                info = f"Alt: {alt:.0f}km, Period: {period:.1f}min"
                ttk.Label(frame, text=info, foreground="gray").pack(side=tk.RIGHT, padx=10)
            except:
                pass
    
    # Display methods
    def _update_color_scheme(self):
        """Update satellite color scheme"""
        scheme = self.color_scheme_var.get()
        
        schemes = {
            "distinct": ['red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 'orange', 'purple'],
            "rainbow": ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet', 'pink'],
            "cool": ['blue', 'cyan', 'teal', 'green', 'purple', 'indigo', 'navy', 'darkblue'],
            "warm": ['red', 'orange', 'yellow', 'gold', 'coral', 'salmon', 'pink', 'magenta']
        }
        
        self.satellite_colors = schemes.get(scheme, schemes["distinct"])
        self._update_orbit_controls()
        self._update_satellite_list()
        self.update_preview()
        self.add_log(f"Color scheme changed to: {scheme}")
    
    def _update_satellite_list(self):
        """Update satellite information list"""
        # Clear existing
        for item in self.satellite_tree.get_children():
            self.satellite_tree.delete(item)
        
        # Add satellites
        for i, sat in enumerate(self.satellites):
            color = self.satellite_colors[i % len(self.satellite_colors)]
            camera = "Yes" if sat.get("camera", {}).get("enabled", False) else "No"
            fault = "Yes" if sat.get("fault", {}).get("enabled", False) else "No"
            
            # Count targets
            target_count = sum(1 for t in self.parent_app.targets 
                             if sat['name'] in t.get('assigned_to', []))
            
            altitude = f"{sat['orbit']['a'] - 6371:.0f}"
            
            self.satellite_tree.insert('', 'end', values=(
                sat['name'], color, camera, fault, target_count, altitude
            ))
    
    def _update_assignment_summary(self):
        """Update target assignment summary"""
        self.assignment_text.config(state="normal")
        self.assignment_text.delete(1.0, tk.END)
        
        total = len(self.parent_app.targets)
        assigned = len([t for t in self.parent_app.targets if t.get('assigned_to', [])])
        
        summary = f"TARGETS: {assigned}/{total} assigned\n"
        summary += "=" * 40 + "\n\n"
        
        # By satellite
        for sat in self.satellites:
            targets = [t['name'] for t in self.parent_app.targets 
                      if sat['name'] in t.get('assigned_to', [])]
            if targets:
                summary += f"{sat['name']}:\n"
                for target in targets:
                    summary += f"  • {target}\n"
                summary += "\n"
        
        # Unassigned
        unassigned = [t['name'] for t in self.parent_app.targets if not t.get('assigned_to', [])]
        if unassigned:
            summary += "UNASSIGNED:\n"
            for target in unassigned:
                summary += f"  ✗ {target}\n"
        
        self.assignment_text.insert(tk.END, summary)
        self.assignment_text.config(state="disabled")
    
    # View control methods
    def _reset_view(self):
        """Reset 3D view to default"""
        self.ax.view_init(elev=20, azim=45)
        self.canvas.draw()
    
    def _constellation_view(self):
        """Set view for constellation overview"""
        self.ax.view_init(elev=30, azim=60)
        self.canvas.draw()
    
    def _earth_view(self):
        """Set view for Earth/target viewing"""
        self.ax.view_init(elev=10, azim=0)
        self.canvas.draw()