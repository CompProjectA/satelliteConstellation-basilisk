#!/usr/bin/env python
"""
rw_fault_help.py

Contains help content and documentation for the spacecraft reaction wheel fault simulator.
This module separates the help content from the GUI code for better organization.
"""

# Dictionary of help content for different topics
HELP_CONTENT = {
    "friction": """Friction Fault

Friction faults simulate increased mechanical friction in a reaction wheel, which requires more torque to overcome static friction.

Parameters:
• Magnitude: The strength of the friction effect (higher values = more friction)
• Wheel: Which reaction wheel to apply the fault to (0-3)
• Time: When to inject the fault during simulation (in minutes)

Effects:
Increased friction causes higher power consumption, reduced wheel speed response, and can lead to attitude control errors if severe.

Periodic Fault:
When enabled, adds cycling friction effects that can simulate intermittent mechanical issues or varying friction levels.
""",
    "power_limit": """Power Limit Fault

Power limit faults restrict the maximum power available to the reaction wheel, limiting its ability to generate torque.

Parameters:
• Power Limit (W): Maximum power the wheel can draw (in Watts)
• Wheel: Which reaction wheel to apply the limit to (0-3)
• Time: When to apply the power limit during simulation (in minutes)

Effects:
Limited power means the wheel cannot respond correctly to large torque commands, causing slower response and potentially instability in spacecraft attitude control.

Periodic Fault:
When enabled, creates cycling power limits that can simulate power system fluctuations or partial failures.
""",
    "encoder": """Encoder Fault

Encoder faults cause errors in the reaction wheel speed feedback to the control system, creating a mismatch between commanded and actual wheel speeds.

Parameters:
• Wheel: Which reaction wheel's encoder to fault (0-3)
• Time: When to inject the fault during simulation (in minutes)

Effects:
With inaccurate speed feedback, the control system cannot properly command the wheel, leading to oscillations, incorrect wheel speeds, and potentially severe attitude control errors.

Periodic Fault:
When enabled, creates intermittent encoder errors that can simulate partial encoder failure or electrical noise in the feedback loop.
""",
    "battery": """Battery Fault

Battery faults simulate increased power consumption or battery degradation in the spacecraft's power system.

Parameters:
• Power Drain (kW): The additional power consumption caused by the fault
• Time: When to inject the fault during simulation (in minutes)

Effects:
Increased power consumption causes faster battery depletion, which can lead to reduced spacecraft capabilities
or complete power failure if the battery level falls below critical thresholds.

Periodic Fault:
When enabled, adds cycling power drain effects that can simulate intermittent electrical issues or varying
power consumption patterns.
""",
    "targets": """Target Management

Targets represent ground locations that the spacecraft can observe or communicate with during its orbit.

Adding Targets:
1. Enter a name for the target location
2. Specify the latitude (-90 to 90 degrees)
3. Specify the longitude (-180 to 180 degrees)
4. Select a color for visualization
5. Click "Add Target"

Visibility:
Targets are visible to the spacecraft when they are above the horizon from the spacecraft's position. The visibility timeline is shown in the simulation results.

Examples:
Default targets include major cities around the world. You can add custom targets like ground stations, research locations, or areas of interest.
""",
    "camera": """Camera Settings

The camera settings control the viewing position for the 3D visualization in Vizard.

Parameters:
• X, Y, Z: Camera position coordinates in the spacecraft body frame

Recommended Settings:
• Side view: (0.0, 2.0, 0.0)
• Top view: (0.0, 0.0, 2.0)
• Front view: (2.0, 0.0, 0.0)

Tips:
Adjust these values to get the best view of the spacecraft and its reaction wheels. The Vizard visualization allows you to rotate and zoom during playback.
""",
    "periodic": """Periodic Fault Injection

Periodic faults cycle on and off during the simulation, creating time-varying effects on the spacecraft.

Parameters:
• Interval (sec): How often the fault repeats (in seconds)
• Magnitude: Strength of the periodic fault effect
• Wheel: Which wheel to apply the periodic fault to (0-3)

Use Cases:
Periodic faults are useful for simulating:
• Intermittent hardware issues
• Orbital thermal cycling effects
• Power system fluctuations
• Noise in sensor or actuator systems

Notes:
The periodic fault runs independently from the one-time fault, allowing you to simulate multiple fault conditions simultaneously.
"""
}

# General help content for different sections
GENERAL_HELP = {
    "overview": """The Spacecraft Reaction Wheel Fault Simulator is a specialized tool for modeling and analyzing the effects of various fault types on spacecraft attitude control systems.

Key Features:
• High-fidelity spacecraft dynamics simulation using the Basilisk astrodynamics framework
• Multiple fault types (friction, power limit, encoder, battery)
• Detailed visualization of spacecraft behavior
• Comprehensive data analysis through plots and visualizations
• Target visibility tracking for ground station analysis

This simulator helps engineers and researchers understand how spacecraft will respond to hardware failures in orbit, allowing for better fault detection, isolation, and recovery strategies.

Getting Started:
1. Configure basic simulation parameters in the "Simulation Settings" tab
2. Select a fault type and configure its parameters in the "Fault Configuration" tab
3. Add ground targets and set camera options in the "Visualization" tab
4. Click "Run Simulation" to execute the simulation
5. View the results in plots and/or the Vizard 3D visualization tool

The simulator creates detailed plots of spacecraft attitude, reaction wheel speeds, and fault-specific metrics to help analyze the effects of the fault on spacecraft performance.
""",
    "simulation": """Simulation Parameters:

• Simulation Time: Duration of the simulation in minutes. Longer simulations show more complete orbit effects but take longer to run. Typical values are 30-60 minutes for a full orbit.

• Output Options:
  - Show Plots: Display plots of simulation results after completion
  - Save Binary: Save a binary file for Vizard 3D visualization
  - Binary Filename: Name for the saved binary file (without extension)

Simulation Process:
1. The simulator initializes the spacecraft in a predefined orbit with reaction wheels for attitude control
2. The spacecraft is commanded to maintain a specific attitude profile
3. At the specified fault time, the selected fault is injected into the reaction wheel
4. The simulation continues, showing the spacecraft's response to the fault
5. Results are displayed in plots and saved for further analysis

Performance Tips:
• Shorter simulations (5-10 minutes) are good for quick tests and fault response analysis
• Longer simulations (30-60 minutes) show full orbital effects and long-term stabilization
• For periodic faults, ensure the simulation runs long enough to capture multiple fault cycles
• Use meaningful binary filenames to identify different test cases
""",
    "faults": """The simulator supports four main types of faults:

1. Friction Fault:
   • Simulates increased mechanical friction in the reaction wheel
   • Parameters: Magnitude (torque in N·m), Wheel (0-3), Time (minutes)
   • Effects: Higher power consumption, reduced response speed, potentially destabilized attitude
   • Real-world causes: Bearing damage, lubricant degradation, mechanical wear
   • Typical magnitude values: 0.0001 to 0.001 N·m

2. Power Limit Fault:
   • Restricts the maximum electrical power available to the wheel
   • Parameters: Power Limit (Watts), Wheel (0-3), Time (minutes)
   • Effects: Torque saturation, inability to respond to large commands
   • Real-world causes: Power system failures, overheating, current limiters
   • Typical power limit values: 0.2 to 1.0 W (small spacecraft)

3. Encoder Fault:
   • Creates measurement errors in the wheel speed feedback
   • Parameters: Wheel (0-3), Time (minutes)
   • Effects: Control oscillations, incorrect wheel speeds, attitude instability
   • Real-world causes: Sensor failure, wiring issues, electronic noise
   • Notes: Unlike other faults, encoder faults affect the control system rather than the hardware

4. Battery Fault:
   • Simulates increased power consumption or battery degradation
   • Parameters: Power Drain (kW), Time (minutes)
   • Effects: Faster battery depletion, potential power failure
   • Real-world causes: Short circuits, damaged cells, thermal issues
   • Typical power drain values: 0.01 to 0.1 kW (10-100W)

Periodic Faults:
   • Optional cycling of fault effects at regular intervals
   • Parameters: Interval (seconds), Magnitude, Wheel (0-3)
   • Useful for simulating intermittent issues or thermal cycling effects
   • Can be combined with any fault type for more complex scenarios

Selecting the appropriate fault parameters:
   • Start with small magnitude values and gradually increase to observe effects
   • Test different wheels to see how fault location affects spacecraft response
   • Try different fault times to see how the vehicle's state at fault injection affects response
   • For realistic scenarios, combine steady and periodic faults to simulate complex failures
""",
    "visualization": """The simulator provides two types of visualization:

1. Data Plots:
   • Attitude error plots show spacecraft pointing accuracy
   • Reaction wheel speed plots show wheel behavior
   • Fault-specific plots highlight the effects of each fault type:
     - Friction plots show increased friction torque
     - Power limit plots show power consumption vs. limits
     - Encoder plots show commanded vs. actual wheel speeds
     - Battery plots show charge level and power consumption
   • Target visibility plots show when ground locations are visible

2. 3D Visualization (Vizard):
   • Requires the Vizard application (included with Basilisk)
   • Shows detailed 3D model of the spacecraft
   • Animates spacecraft motion throughout the simulation
   • Displays Earth and configured target locations
   • Visualizes reaction wheel rotation speeds
   • For battery faults, shows battery charge level indicators

Target Management:
   • Add ground locations by name, latitude, and longitude
   • Each target will appear in the 3D visualization and visibility plots
   • Use targets to simulate ground stations, observation targets, or points of interest
   • The simulator calculates when each target is visible from the spacecraft

Camera Settings:
   • Control the viewpoint in the 3D visualization
   • X, Y, Z coordinates define camera position in spacecraft body frame
   • Presets provide common viewing angles:
     - Side View: See spacecraft from the side
     - Top View: Look down on the spacecraft
     - Front View: See spacecraft from the front
   • Custom positions allow precise control of visualization perspective

Using Vizard:
   1. Run a simulation with "Save Binary" enabled
   2. Open the Vizard application
   3. Load the binary file created by the simulation
   4. Use Vizard controls to play, pause, and navigate the visualization
   5. Observe spacecraft attitude changes in response to the fault
"""
}

# Fault type descriptions for the main interface
FAULT_DESCRIPTIONS = {
    "Friction": "Friction fault adds additional friction to the reaction wheel, simulating mechanical issues like bearing damage.",
    "Power Limit": "Power limit fault restricts the maximum power available to the reaction wheel, simulating power system failures.",
    "Encoder": "Encoder fault causes measurement errors in the reaction wheel speed feedback, leading to control issues.",
    "Battery": "Battery fault simulates increased power consumption or battery degradation, which can lead to power system failure."
}

# Detailed parameter descriptions
PARAMETER_DESCRIPTIONS = {
    "friction": """Friction fault adds a constant friction torque to the reaction wheel, requiring more motor torque to overcome.
Higher magnitude values create more severe effects. A typical spacecraft reaction wheel may experience
friction values in the range of 0.0001 to 0.001 N·m.""",
    "power_limit": """Power limit fault restricts how much electrical power is available to the reaction wheel motor.
When the wheel requires more power than the limit, it cannot produce the required torque.
A typical small spacecraft reaction wheel uses 0.5-5W, so limits below 1W create noticeable effects.""",
    "encoder": """Encoder fault creates errors in the speed measurements fed back to the control system.
The fault simulates a broken encoder, providing incorrect wheel speed information.
This can lead to severe control issues as the spacecraft thinks the wheel is at a different speed than its actual value.""",
    "battery": """Battery fault increases power consumption, simulating damaged or degraded power system components.
Higher values (in kW) create more severe battery drain. For example, a value of 0.05 represents
a 50W power drain, which can significantly impact small spacecraft with limited power generation capacity."""
}

def get_help_content(topic):
    """Get help content for a specific topic"""
    if topic in HELP_CONTENT:
        return HELP_CONTENT[topic]
    return f"No help content available for topic: {topic}"

def get_general_help(section):
    """Get general help content for a section"""
    if section in GENERAL_HELP:
        return GENERAL_HELP[section]
    return f"No general help available for section: {section}"

def get_fault_description(fault_type):
    """Get description for a fault type"""
    if fault_type in FAULT_DESCRIPTIONS:
        return FAULT_DESCRIPTIONS[fault_type]
    return "Unknown fault type"

def get_parameter_description(fault_type):
    """Get parameter description for a fault type"""
    key = fault_type.lower().replace(" ", "_")
    if key in PARAMETER_DESCRIPTIONS:
        return PARAMETER_DESCRIPTIONS[key]
    return ""