#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rw_fault_help.py

Contains help content for the Spacecraft Constellation Fault Simulator.
"""

def get_general_help(section="overview"):
    """Get general help content about the simulator"""
    help_sections = {
        "overview": """Spacecraft Constellation Fault Simulator
=======================================

The Spacecraft Constellation Fault Simulator allows you to model multiple satellites 
in formation, inject various faults, and analyze their effects on spacecraft attitude 
control and mission performance.

This simulator is built on the Basilisk astrodynamics framework, providing high-fidelity 
simulation of spacecraft dynamics in Earth orbit.

Use the different tabs to configure your constellation, set up fault scenarios, 
define target locations, and configure visualization options.

Key features:
- Configure multiple satellites in a constellation
- Inject different fault types into reaction wheels
- Assign observation targets to specific satellites
- Configure camera positions for visualization
- View simulation results in plots and 3D visualization""",

        "simulation": """Simulation Settings
==================

Simulation Parameters:
- Simulation Time: Duration in minutes (30-60 minutes recommended for a full orbit)
- Show Plots: Display results after simulation
- Save Binary: Generate visualization file for Vizard

The simulation uses numerical integration to solve the equations of motion 
for each spacecraft, accounting for orbital dynamics, attitude control, 
and injected faults.

Workflow:
1. Configure your satellite constellation in the Constellation tab
2. Set up fault parameters in the Fault Configuration tab
3. Define and assign targets in the Targets tab
4. Configure visualization in the Visualization tab
5. Set output options in the Output Settings tab
6. Click the Run Simulation button
7. View results in the Results tab or open Vizard for 3D visualization""",

        "visualization": """Visualization Settings
====================

The Visualization tab controls how the simulation appears in Vizard and the preview.

Camera Configuration:
- Select which satellite has the active camera
- Set camera position relative to the satellite body frame
- Adjust field of view (FOV) for wider or narrower viewing angle
- Use presets for common viewing angles (Earth View, Target View, etc.)

Orbit Control:
- Show/hide individual satellite orbits
- Adjust orbit line width and transparency
- Master controls to show/hide all orbits at once

Satellite Display:
- Toggle satellite name labels
- Adjust satellite size multiplier for visibility
- Choose color schemes (distinct, rainbow, cool, warm)
- View satellite information including altitude and assignments

Target Settings:
- Set target altitude above Earth surface
- Adjust target marker size
- Show/hide assignment connections between satellites and targets
- View target assignment summary

The visualization preview shows a simplified 3D view of your constellation.
The actual Vizard visualization will show full orbital dynamics and motion.""",

        "constellation": """Constellation Configuration
=========================

The Constellation tab manages satellite creation and orbital parameters.

Orbit Management:
- Create multiple orbit configurations (Default, MEO Navigation, High Coverage)
- Each orbit has altitude and inclination parameters
- Delete unused orbits (must be empty)

Quick Setup Buttons:
- 4-Sat Default: Creates 4 satellites at 600km altitude
- 6-Sat MEO: Creates 6 satellites at 1200km altitude  
- 2-Sat High: Creates 2 satellites at 2000km altitude

Satellite List:
- Shows all satellites with their orbit and altitude
- Add/remove individual satellites
- Displays orbit assignment and altitude for each satellite

Satellite Details:
- Edit satellite name
- Change orbit assignment (moves satellite to new orbit)
- Adjust position in orbit (true anomaly in degrees)
- View orbital parameters and coverage information

Orbital Parameters:
- Semi-major axis (a): Earth radius + altitude in km
- Eccentricity (e): Shape of orbit (kept small for circular orbits)
- Inclination (i): Angle between orbit and equator
- RAAN (Ω): Right ascension of ascending node
- Argument of periapsis (ω): Orientation of orbit ellipse
- True anomaly (f): Position along the orbit

All satellites in the same orbit share the same orbital plane.""",

        "fault": """Fault Configuration
=================

The Fault tab configures reaction wheel failures for individual satellites.

Fault Types:
1. Friction: Adds mechanical resistance to wheel rotation
   - Simulates bearing damage or contamination
   - Magnitude in N⋅m (higher = more severe)

2. Power Limit: Restricts maximum electrical power to wheel
   - Simulates power system degradation
   - Magnitude in Watts (lower = more severe)

3. Encoder: Causes measurement errors in wheel speed
   - Simulates sensor failures
   - Results in control oscillations

4. Battery: Increases power consumption
   - Simulates battery degradation
   - Magnitude in kW additional drain

Fault Parameters:
- Enable Fault: Turns fault on/off for selected satellite
- Wheel Number: Which reaction wheel (0-3) is affected
- Fault Time: When to inject the fault (minutes into simulation)
- Magnitude: Severity of the fault (units depend on fault type)

Periodic Faults:
- Enable recurring faults at set intervals
- Interval: Time between fault occurrences (seconds)
- Can affect different wheel than main fault

Apply Options:
- Apply to selected satellite only
- Apply same configuration to all satellites""",

        "target": """Target Configuration
==================

The Target tab manages ground locations for satellite observations.

Target List:
- Shows all defined targets with assignment status
- Add/remove targets individually
- Generate random targets at major cities

Target Details:
- Name: Identifier for the target location
- Latitude: -90 to 90 degrees
- Longitude: -180 to 180 degrees  
- Priority: 1-5 (higher = more important)
- Color: Visual marker color in Vizard

Target Assignment:
- Assign targets to specific satellites
- Multiple satellites can observe same target
- Only assigned targets appear in Vizard visualization
- Auto-assign distributes targets evenly

Coverage Map:
- Visual representation of target locations
- Shows assignment status with different markers
- VISIBLE: Target can be seen by assigned satellite
- ASSIGNED: Target assigned but may not be visible
- NOT VISIBLE: Unassigned targets

Coverage Analysis:
- Check Coverage button analyzes all targets
- Reports which targets are visible based on satellite altitude
- Satellites need >200km altitude for good coverage
- Provides recommendations for improving coverage""",

        "output": """Output Settings
=============

The Output tab controls simulation parameters and logging.

Simulation Settings:
- Simulation Time: Total duration in minutes
- Quick presets: 10, 30, 60, 90 minute buttons
- Typical simulations run 30-60 minutes

Output Options:
- Show Plots: Display graphs during/after simulation
- Save Plots: Store plot images to disk
- Save Binary: Create Vizard visualization file

File Directories:
- Logs Directory: Simulation summaries and data
- Plots Directory: Generated graphs and charts
- Vizard Files Directory: Binary files for 3D visualization
- Reset to Defaults: Restore standard directory paths

Binary Filename:
- Name for the Vizard visualization file
- Don't include extension (.bin added automatically)
- File saved in Vizfile/_VizFiles/ directory

Simulation Log:
- Real-time updates during simulation
- Shows configuration changes and progress
- Save log to file for record keeping
- Clear log to remove old messages
- Auto-scroll keeps latest messages visible

The simulation generates:
- Summary text files with configuration details
- Plot images showing spacecraft behavior
- Binary file for Vizard 3D visualization""",

        "results": """Results Tab
===========

The Results tab displays simulation outputs and analysis.

Plot List:
- Shows all generated plots from simulations
- Sorted by creation time (newest first)
- Click to view plot in main display

Plot Viewer:
- Displays selected plot with zoom/pan controls
- Navigation toolbar for plot interaction
- Export plot in various formats (PNG, PDF, SVG)

Available Plots:
- Constellation orbit plots
- Inter-satellite distance graphs
- Fault effect visualizations
- Reaction wheel behavior
- Attitude control performance

Plot Information:
- Filename and creation timestamp
- File size
- Analysis notes

Using Results:
1. Run simulation to generate plots
2. Select plot from list to view
3. Use toolbar to zoom/pan for details
4. Export plots for reports/presentations
5. Open results folder to access all files

The plots help analyze:
- Fault impact on spacecraft performance
- Constellation geometry over time
- Target visibility windows
- System behavior during failures"""
    }

    return help_sections.get(section, help_sections["overview"])