#!/usr/bin/env python
"""
rw_fault_help.py

Contains help content for the Spacecraft Constellation Fault Simulator.
"""

def get_general_help(section="overview"):
    """Get general help content about the simulator"""
    
    help_sections = {
        "overview": """
Spacecraft Constellation Fault Simulator
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
- View simulation results in plots and 3D visualization
""",
        
        "simulation": """
Simulation Settings
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
7. View results in the Results tab or open Vizard for 3D visualization
""",
        
        "visualization": """
Visualization Settings
====================

The simulator outputs visualization files that can be viewed in Vizard:

- Binary files contain time history of spacecraft state
- Camera position can be adjusted to change viewpoint
- Ground targets are shown on Earth's surface
- Reaction wheels are visualized with their spin axes

To view the visualization:
1. Run a simulation with 'Save Binary' enabled
2. Open the Vizard application
3. Load the generated binary file from the Vizfile directory

Camera Configuration:
- Each satellite can have one camera defined
- Only one satellite's camera can be active at a time
- Camera position is defined in the body frame of the satellite
- Field of view controls how much the camera can see
- Use the preset buttons for common viewing angles

For best constellation visibility:
- Use the Earth View preset with a Z value of 50000 or higher
- This positions the camera far above the spacecraft
- Increase field of view to see more of the environment

Vizard Controls:
- In Vizard, use the camera dropdown to switch between views
- Adjust the Field of View to see more or less of the scene
- Use the time slider at the bottom to control simulation playback
- Camera settings can be adjusted in the left panel
"""
    }
    
    return help_sections.get(section, help_sections["overview"])

def get_help_content(topic):
    """Get help content for a specific topic"""
    
    help_topics = {
        "constellation": """
Constellation Configuration
=========================

The Constellation tab allows you to configure multiple satellites in formation.

You can add, remove, and modify satellite orbital parameters:
- Semi-major axis (a): Distance from the center of orbit to the farthest point
- Eccentricity (e): How much the orbit deviates from a perfect circle
- Inclination (i): Tilt of the orbit with respect to the equator
- RAAN (Ω): Right Ascension of the Ascending Node
- Argument of Periapsis (ω): Defines the orientation of the orbit
- True Anomaly (f): Position of the satellite in its orbit

Adding Satellites:
1. Click "Add Satellite" to add individual satellites
2. Use "Add Multiple" to create a constellation with evenly spaced satellites
3. The "True Anomaly" parameter determines the satellite's position in orbit

Tips:
- For better visualization, use semi-major axis values above 7000 km
- Set different inclinations for more interesting orbital dynamics
- Space satellites evenly in true anomaly for symmetric constellations
""",
        
        "fault": """
Fault Configuration
=================

The Fault Configuration tab allows you to inject various faults into spacecraft components.

Supported fault types include:

Friction Fault:
- Adds additional friction to the reaction wheel, simulating mechanical issues like bearing damage
- Magnitude: 0.0001-0.001 typical (dimensionless friction coefficient)
- Effect: Reduces wheel acceleration capability and increases power consumption

Power Limit Fault:
- Restricts maximum power available to the reaction wheel, simulating power system limitations
- Magnitude: 0.2-1.0 typical (Watts)
- Effect: Limits maximum torque the wheel can provide

Encoder Fault:
- Causes measurement errors in the wheel speed feedback
- Effect: Results in attitude control errors and potential instability

Battery Fault:
- Simulates increased power consumption or battery degradation
- Magnitude: 0.01-0.1 typical (kW drain)
- Effect: Reduces available power for all spacecraft systems

Fault Configuration:
- Enable the fault using the checkbox
- Select the fault type from the radio buttons
- Set the magnitude appropriate for the fault type
- Select which reaction wheel to affect (0-3)
- Set the time when the fault occurs (in minutes)

Periodic Faults:
- Enable periodic faults to simulate intermittent issues
- Set interval (seconds), magnitude, and affected wheel
- Useful for simulating cyclical environmental effects or degradation
""",
        
        "target": """
Target Configuration
==================

The Targets tab allows you to define ground locations that the spacecraft will observe.

For each target, you can specify:
- Name: Identifier for the location
- Latitude: North-south position (-90° to 90°)
- Longitude: East-west position (-180° to 180°)
- Color: Visual representation in the simulation
- Priority: Importance level (1-5)

Target Assignment:
- Select a target from the list
- Choose a satellite from the dropdown
- Click "Assign Target" to link them
- Assignments are shown in the Current Assignments box

Auto-Assignment:
- The "Auto-Assign Targets" button will distribute targets among satellites
- Higher priority targets are assigned first
- Targets are distributed evenly among available satellites

Target Map:
- Shows a visual representation of all targets
- Assigned targets are shown with stars, unassigned with circles
- Larger markers indicate higher priority targets
- Colors match those used in the visualization

Random Target Generation:
- "Generate Random Targets" creates locations around the world
- Useful for quickly populating the simulation with observation points
- You can clear existing targets when generating new ones
""",
        
        "visualization": """
Visualization Configuration
=========================

The Visualization tab configures how the spacecraft simulation is rendered in 3D.

Camera Settings:
- Position: Location of the camera in the body frame
- Field of View: Width of the visible area (in degrees)
- Only one camera can be active at a time

Presets:
- Side View: Camera positioned to view from the side
- Top View: Camera positioned above the spacecraft
- Front View: Camera positioned in front of the spacecraft
- Earth View: Camera positioned far away to view the constellation

For Best Results:
- Use high Z values (15000+) to see the entire constellation
- Use wider field of view (90-120 degrees) for better visibility
- Enable only one camera to avoid conflicts
- The "Disable Other Cameras" button ensures only one camera is active

Camera Preview:
- Shows a 3D representation of the spacecraft and camera
- Helps visualize the field of view and camera positioning
- Blue cone represents the camera's field of view
- Colored dots represent assigned targets

Visualization Binary Files:
- Generated binary files can be opened in Vizard
- Binaries are saved in the Vizfile directory
- Each simulation creates a new binary file
- Target locations and spacecraft trajectories are included
""",
        
        "output": """
Output Settings
=============

The Output Settings tab controls how simulation results are saved and displayed.

Simulation Settings:
- Simulation Time: Duration of the simulation in minutes
- Show Plots: Display plots immediately after simulation
- Save Plots: Save plots as image files
- Save Binary: Generate visualization file for Vizard

Directory Settings:
- Logs Directory: Where simulation logs and summaries are saved
- Plots Directory: Where plot images are saved
- Vizard Files Directory: Where binary visualization files are saved

Binary Filename:
- Base name for the visualization binary file
- "_UnityViz.bin" will be appended automatically
- Use only alphanumeric characters and underscores

Simulation Log:
- Shows real-time messages during simulation
- Records errors and warnings
- Displays information about the simulation progress
- Can be cleared with the "Clear Log" button

Opening Results:
- Results can be viewed directly in the Results tab
- "Open Results Folder" button opens the plots directory
- Binary files can be opened in Vizard after simulation
""",
        
        "friction": """
Friction Fault Details
====================

The friction fault simulates increased friction in a reaction wheel, which can occur due to:
- Bearing damage or wear
- Lubricant degradation
- Contamination or debris
- Mechanical misalignment

Parameters:
- Magnitude: Friction coefficient (typical range: 0.0001-0.001)
- Wheel: Which reaction wheel to affect (0-3)
- Time: When the fault occurs during simulation (minutes)

Effects:
- Reduced wheel acceleration
- Increased power consumption
- Potential attitude control errors
- Possible wheel speed saturation

Characteristic Plots:
- Wheel speed shows reduced acceleration and deceleration
- Friction torque increases for the affected wheel
- Higher power consumption for the affected wheel
- Attitude errors may increase during maneuvers

Mitigation Strategies:
- Fault detection and isolation
- Controller retuning
- Wheel speed avoidance regions
- Switching to redundant wheels
""",
        
        "power_limit": """
Power Limit Fault Details
=======================

The power limit fault simulates restrictions in the maximum power available to a reaction wheel, which can occur due to:
- Battery degradation
- Solar panel damage
- Power distribution issues
- Thermal constraints

Parameters:
- Magnitude: Maximum power limit in Watts (typical range: 0.2-1.0)
- Wheel: Which reaction wheel to affect (0-3)
- Time: When the fault occurs during simulation (minutes)

Effects:
- Limited maximum torque
- Reduced maneuver capability
- Increased settling time
- Potential stability issues during high-demand maneuvers

Characteristic Plots:
- Wheel torque saturates at the power limit
- Wheel speed shows rate limiting
- Attitude errors increase during rapid maneuvers
- Control effort may oscillate near the limit

Mitigation Strategies:
- Reduced maneuver rates
- Modified controller gains
- Load balancing across wheels
- Alternate attitude control methods (e.g., magnetorquers)
""",
        
        "encoder": """
Encoder Fault Details
===================

The encoder fault simulates errors in the wheel speed measurements, which can occur due to:
- Sensor damage
- Electrical interference
- Communication errors
- Software glitches

Parameters:
- Wheel: Which reaction wheel encoder to affect (0-3)
- Time: When the fault occurs during simulation (minutes)

Effects:
- Inaccurate control feedback
- Oscillations in wheel speed
- Attitude pointing errors
- Potential instability

Characteristic Plots:
- Reported wheel speed differs from actual speed
- Wheel command shows oscillatory behavior
- Attitude errors exhibit high-frequency components
- Control effort may become erratic

Mitigation Strategies:
- Sensor filtering and fusion
- Robust control techniques
- Observer-based state estimation
- Fault detection and reconfiguration
""",
        
        "battery": """
Battery Fault Details
===================

The battery fault simulates increased power consumption or battery degradation, which can affect the entire spacecraft:
- Cell degradation
- Short circuits
- Charge controller failure
- Thermal runaway

Parameters:
- Magnitude: Additional power drain in kW (typical range: 0.01-0.1)
- Time: When the fault occurs during simulation (minutes)

Effects:
- Reduced available power for all systems
- Faster battery depletion
- Potential power cycling of non-critical components
- Limited operational capability

Characteristic Plots:
- Battery state of charge decreases more rapidly
- Power generation vs. consumption balance shifts
- System voltage may show drops during high-demand periods
- Power management system may activate load shedding

Mitigation Strategies:
- Power conservation modes
- Reduced duty cycles
- Prioritization of critical systems
- Adjusting orbit for better solar illumination
"""
    }
    
    return help_topics.get(topic, f"No help content available for: {topic}")