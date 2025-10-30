#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rw_fault_help.py

Contains help content for the Spacecraft Constellation Fault Simulator.
Updated to include Cluster Management, Fault Detection, Task Reassignment, and DRL tabs.
"""

def get_general_help(section="overview"):
    """Get general help content about the simulator"""
    help_sections = {
        "overview": """Spacecraft Constellation Fault Simulator
=======================================

The Spacecraft Constellation Fault Simulator allows you to model multiple satellites 
in formation, inject various faults, and analyze their effects on spacecraft attitude 
control and mission performance using advanced AI/ML techniques.

This simulator is built on the Basilisk astrodynamics framework, providing high-fidelity 
simulation of spacecraft dynamics in Earth orbit, enhanced with machine learning-based 
fault detection and deep reinforcement learning for autonomous task reassignment.

Available Tabs:
1. Constellation: Configure clusters, formations, and individual satellites
2. Fault Configuration: Inject various fault types into reaction wheels
3. Fault Detection: ML-based anomaly detection and fault identification
4. Task Reassignment: DRL-driven autonomous task redistribution
5. Targets: Define ground observation targets and coverage
6. Visualization: Configure 3D view settings for Vizard
7. Output Settings: Control simulation parameters and logging
8. Results: View plots, metrics, and analysis outputs
9. DRL: Train and manage DRL agents for task reassignment

Key Features:
- Cluster-based constellation management with multiple formation types
- Leader-follower satellite coordination with communication links
- AI-powered fault detection using autoencoders and statistical analysis
- Deep reinforcement learning for autonomous task reassignment
- Realistic orbital mechanics and attitude control simulation
- Comprehensive visualization with Vizard 3D playback
- Real-time telemetry monitoring and anomaly detection
- Export plots, logs, and binary files for detailed analysis""",

        "simulation": """Simulation Settings
==================

Simulation Parameters:
- Simulation Time: Duration in minutes (30-60 minutes recommended for a full orbit)
- Show Plots: Display results after simulation
- Save Binary: Generate visualization file for Vizard

The simulation uses numerical integration to solve the equations of motion 
for each spacecraft, accounting for orbital dynamics, attitude control, 
and injected faults.

Enhanced Simulation Pipeline:
1. Basilisk dynamics simulation with fault injection
2. Automated ML fault detection analysis
3. DRL-based task reassignment (if faults detected)
4. Enhanced cluster plotting and communication analysis
5. Results export with comprehensive metrics

Workflow:
1. Configure your satellite constellation in the Constellation tab
2. Set up fault parameters in the Fault Configuration tab
3. Configure ML detection in Fault Detection tab (optional)
4. Define and assign targets in the Targets tab
5. Configure visualization in the Visualization tab
6. Set output options in the Output Settings tab
7. Click the Run Simulation button
8. View results in the Results tab or open Vizard for 3D visualization
9. Check Task Reassignment tab for DRL analysis (if faults occurred)

Simulation generates:
- Telemetry data for all spacecraft
- ML detection results with confidence scores
- DRL reassignment recommendations
- Enhanced cluster communication plots
- Standard attitude and orbit plots
- Vizard binary file for 3D playback""",

        "constellation": """Constellation Configuration
=========================

The Constellation tab has four sub-tabs for complete satellite management:

1. CLUSTER MANAGEMENT
- Create New Cluster: Define cluster name, formation type, and number of satellites
- Formation Types: Leader-Follower, Line, Triangle, Diamond
- Separation Distance: Inter-satellite spacing (min 15km for Vizard visibility)
- Cluster Orbit: Select shared orbit for all cluster members
- Existing Clusters: View, modify, delete, or test communication

2. ORBIT MANAGEMENT
- Create multiple orbit configurations (LEO 600km, MEO 1200km, GEO, etc.)
- Each orbit defines altitude and inclination
- Delete unused orbits (must be empty)
- All satellites in same orbit share orbital plane

3. INDIVIDUAL SATELLITES
- Quick Setup: 4-Sat Default, 6-Sat MEO, 2-Sat High presets
- Satellite List: Shows all satellites with orbit and cluster assignment
- Add/Remove: Manage individual satellites
- Satellite Details: Edit name, orbit, position (true anomaly)
- Cluster Assignment: Change satellite cluster membership

4. COMMUNICATION
- View cluster communication topology
- Leader-follower link status
- Communication range and field-of-view settings
- Test communication links between satellites

Cluster Features:
- Up to 4 clusters supported per simulation
- Each cluster needs: 1 leader + 1+ children (followers)
- Formation geometry automatically calculated
- Communication links visualized in plots
- Enhanced cluster plots show: formation geometry, inter-satellite distances, 
  communication status, and relative positions

Orbital Parameters:
- Semi-major axis (a): Earth radius + altitude in km
- Eccentricity (e): Shape of orbit (kept small for circular orbits)
- Inclination (i): Angle between orbit and equator
- RAAN (Ω): Right ascension of ascending node
- Argument of periapsis (ω): Orientation of orbit ellipse
- True anomaly (f): Position along the orbit

Important Notes:
- All satellites in a cluster share the same orbit
- Separation must be ≥15km for Vizard visibility
- Formation geometry is calculated automatically based on type
- Leader coordinates cluster activities and communication""",

        "fault": """Fault Configuration
=================

The Fault tab configures reaction wheel failures for individual satellites.

Fault Types:
1. Friction: Adds mechanical resistance to wheel rotation
   - Simulates bearing damage or contamination
   - Magnitude in N·m (higher = more severe)
   - Example: 0.0005 N·m = mild friction, 0.005 N·m = severe

2. Power Limit: Restricts maximum electrical power to wheel
   - Simulates power system degradation
   - Magnitude in Watts (lower = more severe)
   - Example: 50W = significant limitation

3. Encoder: Causes measurement errors in wheel speed
   - Simulates sensor failures
   - Results in control oscillations
   - Magnitude: measurement error percentage

4. Battery: Increases power consumption
   - Simulates battery degradation
   - Magnitude in kW additional drain
   - Example: 0.5 kW = significant drain

5. RW (Reaction Wheel): Complete wheel failure
   - Simulates mechanical breakdown
   - Wheel becomes completely non-functional
   - Requires spacecraft to use redundant wheels

Fault Parameters:
- Enable Fault: Turns fault on/off for selected satellite
- Wheel Number: Which reaction wheel (0-3) is affected
- Fault Time: When to inject the fault (minutes into simulation)
- Magnitude: Severity of the fault (units depend on fault type)

Periodic Faults:
- Enable recurring faults at set intervals
- Interval: Time between fault occurrences (seconds)
- Can affect different wheel than main fault
- Useful for testing recovery and adaptation

Apply Options:
- Apply to Selected: Configure fault for current satellite only
- Apply to All: Apply same fault configuration to entire constellation
- Apply to Cluster: Apply to all satellites in a specific cluster

Integration with ML Detection:
- All configured faults are automatically detected during simulation
- ML system analyzes telemetry for anomalies
- Detection results appear in Fault Detection tab
- DRL system responds to detected faults automatically""",

        "fault_detection": """Fault Detection System
====================

The Fault Detection tab provides AI-powered anomaly detection and fault identification.

OVERVIEW:
Uses machine learning (ML) autoencoder models combined with statistical analysis 
to automatically detect spacecraft faults during and after simulation.

DETECTION METHODS:

1. ML Autoencoder Detection
   - Neural network trained on normal spacecraft behavior
   - Detects anomalies by reconstruction error
   - Confidence score: 0.0 (no fault) to 1.0 (definite fault)
   - Real-time analysis during simulation

2. Statistical Analysis
   - Compares pre-fault and post-fault statistics
   - Detects changes in mean, variance, trends
   - Works without pre-trained models
   - Good for unknown fault types

3. Trend Analysis
   - Identifies sudden changes in data trends
   - Uses moving averages and derivatives
   - Detects gradual degradation
   - Useful for predictive maintenance

4. Threshold-Based Detection
   - Simple rule-based detection
   - Compares values against fixed limits
   - Fast and reliable for known patterns
   - Requires manual threshold configuration

SETUP PROCESS:

1. Load ML Model:
   - Browse to your trained Keras model (.keras or .h5)
   - Default: core/anomaly_detection_model.keras
   - Model should be trained on normal telemetry

2. Configure Detection Methods:
   - Enable/disable: ML, Statistical, Trend, Threshold
   - Multiple methods can run simultaneously
   - Results are combined for robust detection

3. Set Detection Parameters:
   - ML Confidence Threshold: 0.01 (sensitive) to 1.0 (conservative)
   - Anomaly Score Threshold: Reconstruction error limit
   - Speed Change Threshold: % change to trigger detection
   - Attitude Error Threshold: Max acceptable error

4. Real-Time Monitoring (Optional):
   - Enable for live detection during simulation
   - Set update interval (seconds)
   - View plots updating in real-time
   - Higher intervals = better performance

RUNNING DETECTION:

Automatic Mode:
- Detection runs automatically when simulation completes
- Results appear in this tab after simulation
- ML and statistical analysis performed on full dataset
- No user interaction required

Manual Mode:
- Use "Run Detection" button to analyze current data
- Specify target spacecraft or analyze all
- Adjust thresholds and re-run as needed
- Export results to file

INTERPRETING RESULTS:

Detection Summary:
- Total Detections: Number of faults found
- Spacecraft Affected: Which satellites have faults
- Detection Confidence: Average confidence score (0-1)
- Fault Timing: When each fault was detected

Detailed View:
- Per-spacecraft detection status
- Confidence scores for each method
- Telemetry plots with fault annotations
- Anomaly score time series
- Threshold crossings highlighted

Plots Generated:
- Detection confidence over time
- Anomaly scores per spacecraft
- Telemetry with fault markers
- Statistical analysis results
- Method comparison charts

TROUBLESHOOTING:

"ML Not Available":
- Install TensorFlow: pip install tensorflow
- Check model file path and permissions
- Verify model format (.keras or .h5)

"Model Load Failed":
- Ensure model file exists at specified path
- Check model was trained with compatible TensorFlow version
- Verify model input shape matches telemetry features

"No Detections" (but faults injected):
- Decrease ML confidence threshold (try 0.3-0.5)
- Enable multiple detection methods
- Check that fault magnitude is significant
- Verify telemetry data was generated correctly

"Too Many False Positives":
- Increase confidence threshold (try 0.7-0.9)
- Use ML method only (disable statistical/threshold)
- Retrain model with more diverse normal data
- Adjust anomaly score threshold upward

INTEGRATION WITH DRL:
- Detection results automatically feed into DRL system
- Faulty spacecraft identified for task reassignment
- Healthy spacecraft identified as reassignment targets
- System health status computed from detection results
- Full integration pipeline runs automatically if faults detected""",

        "task_reassignment": """Task Reassignment System
======================

The Task Reassignment tab displays results from the Deep Reinforcement Learning (DRL) 
system that automatically redistributes tasks when faults are detected.

OVERVIEW:
When faults are detected during simulation, the DRL system automatically determines 
the best strategy to reassign tasks from faulty spacecraft to healthy ones, ensuring 
mission continuity and optimal performance.

SYSTEM WORKFLOW:

1. Fault Detection
   - ML system identifies faulty spacecraft
   - Confidence scores computed for each detection
   - Spacecraft categorized as healthy or faulty

2. Health Assessment
   - System health status: operational, degraded, or critical
   - Spacecraft capabilities evaluated
   - Current task load assessed

3. DRL Decision Making
   - Trained PPO agent analyzes system state
   - Evaluates multiple reassignment strategies
   - Selects optimal strategy based on:
     * Number of healthy vs. faulty spacecraft
     * Current task distribution
     * Spacecraft capabilities
     * Mission priorities

4. Task Redistribution
   - Tasks moved from faulty to healthy spacecraft
   - Load balancing applied
   - Communication links verified
   - Assignments updated in system

5. Results Display
   - Reassignment plan shown in this tab
   - System health metrics displayed
   - Before/after comparison provided

REASSIGNMENT STRATEGIES:

1. Even Distribution
   - Tasks distributed equally among healthy spacecraft
   - Simplest strategy, good for homogeneous constellation
   - Minimizes maximum load on any spacecraft

2. Capability-Based
   - More tasks to spacecraft with better capabilities
   - Considers power, pointing accuracy, sensors
   - Optimal for heterogeneous constellation

3. Load-Balanced
   - Distribution based on current spacecraft capacity
   - Accounts for existing task load
   - Prevents overloading any spacecraft

4. Priority-Based
   - High-priority tasks assigned to best spacecraft
   - Lower priorities distributed evenly
   - Ensures critical tasks maintained

DISPLAYED INFORMATION:

System Health Summary:
- Overall Status: Operational / Degraded / Critical
- Healthy Spacecraft Count
- Faulty Spacecraft Count
- Total Tasks in System
- Tasks Successfully Reassigned

DRL State Vector:
- Input state provided to DRL agent
- Includes: health ratios, confidence scores, capabilities
- Shows what information DRL used for decision

DRL Decision Details:
- Selected Strategy: Which reassignment approach
- Action Index: Numerical action from DRL agent
- Decision Confidence: How confident DRL is (0-1)
- Model Used: Whether DRL model was available

Spacecraft Status Table:
- Name: Spacecraft identifier
- Status: Healthy / Faulty
- Faults: Specific fault types detected
- Previous Tasks: Tasks before reassignment
- New Tasks: Tasks after reassignment
- Load Change: Increase/decrease in task count

Task Reassignment Plan:
- From/To: Which tasks moved where
- Priority: Task importance level
- Reason: Why task was reassigned
- Success: Whether reassignment succeeded

Console Output:
- Step-by-step log of reassignment process
- DRL decision making details
- Task transfer confirmations
- Fallback actions if needed

WHEN REASSIGNMENT OCCURS:

Automatic Trigger:
- Any fault detected during simulation
- ML confidence above threshold
- At least one healthy spacecraft available
- Tasks exist to reassign

Manual Override:
- Not currently supported in GUI
- Would require integration scripting
- Future feature: manual reassignment button

INTERPRETING RESULTS:

Good Reassignment:
- All tasks redistributed successfully
- Load balanced across healthy spacecraft
- High decision confidence (>0.7)
- No spacecraft overloaded

Degraded Reassignment:
- Some tasks couldn't be reassigned
- Load imbalanced
- Medium confidence (0.4-0.7)
- Some spacecraft near capacity limits

Failed Reassignment:
- No reassignment possible
- All spacecraft faulty or overloaded
- Low confidence (<0.4)
- System enters safe mode

FALLBACK BEHAVIOR:

If DRL Not Available:
- System falls back to rule-based reassignment
- Uses simple even distribution strategy
- Still functional but less optimal
- Message shown in results

If Reassignment Fails:
- Faulty spacecraft enters safe mode
- Tasks marked as unassigned
- Operator alerts generated
- Manual intervention required

REQUIREMENTS:

Software:
- TensorFlow installed (pip install tensorflow)
- DRL model files in DRL/ directory
- Integration bridge (main_integration_script.py)

Files:
- anomaly_detection_model.keras (for fault detection)
- PPO trained model (or uses rule-based fallback)
- Configuration files in core/

Simulation:
- At least one fault configured and injected
- Fault detection enabled
- Healthy spacecraft available for reassignment""",

        "drl": """Deep Reinforcement Learning (DRL) Training
======================================

The DRL tab provides tools for training and managing reinforcement learning agents 
for autonomous task reassignment.

OVERVIEW:
Train PPO (Proximal Policy Optimization), TDHD, or DQN agents to make intelligent 
decisions about task reassignment when spacecraft faults occur. Agents learn optimal 
strategies through trial and error in simulated fault scenarios.

TAB SECTIONS:

1. OVERVIEW
   - System status: Shows availability of DRL components
   - Component check: TensorFlow, DRL Bridge, Integration System
   - Getting started guide
   - Workflow explanation

2. TRAINING
   - Algorithm selection: PPO / TDHD / DQN
   - Training parameters configuration
   - Start/stop training controls
   - Real-time training log
   - Progress monitoring

3. INTEGRATION
   - Link to main simulation workflow
   - Explains automatic integration
   - Testing capabilities
   - Configuration summary

4. RESULTS
   - Browse training output files
   - View training metrics and plots
   - Load and compare model performance
   - Export results for analysis

TRAINING ALGORITHMS:

1. PPO (Proximal Policy Optimization)
   - Recommended for most users
   - Stable and reliable training
   - Good sample efficiency
   - Produces robust policies
   - Training script: DRL/PPO.py

2. TDHD (Twin Delayed DDPG with Hindsight)
   - Advanced continuous control
   - Uses hindsight experience replay
   - Better for complex state spaces
   - Requires more tuning
   - Training script: DRL/TDHD.py

3. DQN (Deep Q-Network)
   - Discrete action space
   - Simpler than policy gradient methods
   - Good for finite action sets
   - Fast training
   - Training script: DRL/DQN.py

TRAINING PARAMETERS:

Basic Settings:
- Iterations: Number of training episodes (100-1000)
- Satellites: Number of spacecraft in training scenarios (4-8)
- Targets: Number of observation targets (4-12)
- Learning Rate: Step size for weight updates (0.0001-0.01)
- Batch Size: Samples per training update (32-256)

Advanced Settings:
- Discount Factor (Gamma): Future reward weight (0.9-0.99)
- Exploration Rate (Epsilon): Random action probability (0.1-0.5)
- Hidden Layers: Neural network architecture (2-4 layers)
- Layer Size: Neurons per layer (64-512)
- Update Frequency: Steps between model updates (5-50)

Environment Settings:
- Fault Types: Which faults to include in training
- Fault Frequency: How often faults occur
- Task Complexity: Number and type of tasks
- Success Criteria: What constitutes successful reassignment

TRAINING PROCESS:

1. Configure Algorithm:
   - Select PPO/TDHD/DQN
   - Set training parameters
   - Choose output directory

2. Start Training:
   - Click "Start Training" button
   - Script runs in background
   - Log shows real-time progress
   - Can take minutes to hours

3. Monitor Progress:
   - Training log updates continuously
   - Shows episode rewards
   - Reports success rate
   - Indicates convergence

4. Training Completes:
   - Final model saved to DRL/result/
   - Training plots generated
   - Metrics logged to file
   - Model ready for use

USING TRAINED MODELS:

Automatic Integration:
- Place trained model in DRL/result/
- Update drl_config.py with model path
- Model loads automatically during simulation
- Used for task reassignment decisions

Manual Testing:
- Load model in Integration tab
- Run test scenarios
- Compare with rule-based baseline
- Evaluate decision quality

Model Selection:
- Multiple models can be saved
- Choose best based on training metrics
- Test different models in simulations
- Keep best performer

RESULTS AND ANALYSIS:

Training Plots:
- Episode rewards over time
- Success rate curve
- Loss function values
- Action distribution
- State visitation frequency

Model Files:
- Trained weights (.h5 or .keras)
- Training configuration (.json)
- Replay buffer (if applicable)
- Tensorboard logs

Performance Metrics:
- Average reward per episode
- Success rate (tasks successfully reassigned)
- Training time and convergence
- Final policy quality

Comparison Tools:
- Compare different algorithms
- Baseline vs. trained agent
- Before/after training performance
- Ablation studies

TROUBLESHOOTING:

"DRL Not Available":
- Install TensorFlow: pip install tensorflow
- Check DRL files in DRL/ directory
- Verify Python path includes DRL/
- Restart application

"Training Fails to Start":
- Check training script path
- Verify parameters are valid
- Ensure output directory writable
- Review log for specific errors

"Training Doesn't Converge":
- Reduce learning rate
- Increase training iterations
- Simplify environment
- Adjust reward function
- Check for bugs in environment code

"Model Performance Poor":
- Train for more iterations
- Tune hyperparameters
- Add more training scenarios
- Increase network size
- Use different algorithm

WORKFLOW INTEGRATION:

Standalone Training:
- Use this tab to train agents offline
- Test different configurations
- Optimize hyperparameters
- Generate multiple models

Integrated Usage:
- Trained models used automatically in main simulation
- No need to explicitly load model
- System selects best available model
- Falls back to rules if model unavailable

Best Practices:
- Train multiple models with different parameters
- Test models before deployment
- Keep training logs for comparison
- Backup successful models
- Document training conditions

REQUIREMENTS:

Software:
- Python 3.7+ with TensorFlow 2.x
- NumPy, Matplotlib, Pandas
- DRL training scripts in DRL/
- Integration bridge in core/

Hardware:
- CPU training: Works but slow
- GPU recommended for faster training
- At least 8GB RAM
- Sufficient disk space for models

Files:
- Training scripts: DRL/PPO.py, TDHD.py, DQN.py
- Environment: DRL/Envs.py
- Tools: DRL/Tools.py
- Config: core/drl_config.py""",

        "target": """Target Configuration
==================

The Target tab manages ground locations for satellite observations.

Target List:
- Shows all defined targets with assignment status
- Add/remove targets individually
- Generate random targets at major cities
- Import targets from file (CSV format)

Target Details:
- Name: Identifier for the target location
- Latitude: -90 to 90 degrees (negative = South)
- Longitude: -180 to 180 degrees (negative = West)  
- Priority: 1-5 (higher = more important)
- Color: Visual marker color in Vizard
- Altitude: Elevation above sea level (meters)

Target Assignment:
- Assign targets to specific satellites
- Multiple satellites can observe same target
- Only assigned targets appear in Vizard visualization
- Auto-assign distributes targets evenly across constellation
- Cluster-aware assignment: assign to cluster leaders

Coverage Map:
- Visual representation of target locations on world map
- Shows assignment status with different markers:
  * VISIBLE: Target can be seen by assigned satellite
  * ASSIGNED: Target assigned but may not be visible currently
  * NOT VISIBLE: Unassigned targets
- Real-time coverage updates during simulation

Coverage Analysis:
- "Check Coverage" button analyzes all targets
- Reports which targets are visible based on:
  * Satellite altitude and orbit inclination
  * Target latitude and longitude
  * Communication constraints
  * Minimum elevation angle
- Satellites need >200km altitude for good coverage
- Provides recommendations for improving coverage:
  * Adjust satellite orbits
  * Reassign targets
  * Add more satellites

Target Priorities in Task Reassignment:
- High-priority targets (4-5) assigned to healthy spacecraft first
- Medium-priority (2-3) distributed based on availability
- Low-priority (1) may be dropped if resources limited
- DRL system considers priorities when reassigning tasks

Integration with Fault Detection:
- If assigned satellite fails, target becomes unassigned
- DRL task reassignment automatically reassigns targets
- Coverage maintained through intelligent redistribution
- Priority determines reassignment urgency""",

        "visualization": """Visualization Settings
====================

The Visualization tab controls how the simulation appears in Vizard and plots.

Camera Configuration:
- Select which satellite has the active camera
- Set camera position relative to the satellite body frame
- Adjust field of view (FOV) for wider or narrower viewing angle
- Use presets for common viewing angles:
  * Earth View: Look toward Earth
  * Target View: Focus on assigned target
  * Formation View: See entire constellation
  * Close-up: Detailed satellite view
  * Wide Angle: Overview of orbit

Orbit Control:
- Show/hide individual satellite orbits
- Adjust orbit line width (1-5 pixels)
- Set orbit transparency (0-100%)
- Master controls to show/hide all orbits at once
- Color orbits by: satellite, cluster, orbit type

Satellite Display:
- Toggle satellite name labels on/off
- Adjust satellite size multiplier for visibility (1-10x)
- Choose color schemes:
  * Distinct: Each satellite different color
  * Rainbow: Spectrum gradient
  * Cool: Blues and greens
  * Warm: Reds and oranges
  * Cluster: Color by cluster membership
- View satellite information: altitude, speed, assignments
- Highlight faulty satellites with special markers

Cluster Visualization:
- Show/hide cluster formation geometry
- Display inter-satellite communication links
- Highlight leader-follower relationships
- Show separation distances
- Formation shape overlays (line, triangle, diamond)

Target Settings:
- Set target altitude above Earth surface (0-1000m)
- Adjust target marker size (0.5-5.0x)
- Show/hide assignment connections (lines from satellites to targets)
- View target assignment summary
- Color targets by priority level
- Display coverage zones (circles around targets)

Communication Links:
- Visualize satellite-to-satellite links
- Show communication range (cones or spheres)
- Display link quality/strength
- Highlight active vs. inactive links
- Leader-follower connections emphasized

Advanced Options:
- Time acceleration: Speed up/slow down playback
- Grid overlay: Show latitude/longitude grid
- Axis display: Show X, Y, Z reference axes
- Lighting: Configure sun position and shadows
- Background: Stars, Earth texture, space color

The visualization preview in the GUI shows a simplified 3D view of your constellation.
The actual Vizard visualization (from .bin file) will show:
- Full orbital dynamics and motion over time
- High-quality Earth rendering with textures
- Smooth camera transitions
- Time-stamped events (fault injection, target passes)
- Animated reaction wheel states
- Playback controls (play, pause, fast-forward, rewind)

Export Settings:
- Binary filename: Name for .bin file (no extension needed)
- Output directory: Where to save Vizard file (Vizfile/_VizFiles/)
- Include telemetry: Embed data in binary for playback analysis
- Compression: Reduce file size (may increase load time)

Opening in Vizard:
1. Run simulation with "Save Binary" enabled
2. Find .bin file in Vizfile/_VizFiles/ directory
3. Open Vizard application
4. File → Open → Select your .bin file
5. Use Vizard controls for playback and analysis
6. Can export frames, create videos, or take screenshots""",

        "output": """Output Settings
=============

The Output tab controls simulation parameters and logging.

Simulation Settings:
- Simulation Time: Total duration in minutes
  * Minimum: 5 minutes (partial orbit)
  * Recommended: 30-60 minutes (1-2 orbits)
  * Maximum: 120 minutes (multiple orbits)
- Quick presets: 10, 30, 60, 90 minute buttons
- Real-time checkbox: Run at actual time scale (slow)

Output Options:
- Show Plots: Display graphs during/after simulation
- Save Plots: Store plot images to disk (PNG format)
- Save Binary: Create Vizard visualization file (.bin)
- Auto-open Results: Switch to Results tab when complete
- Generate Report: Create PDF summary of simulation

File Directories:
- Logs Directory: Simulation summaries and data
  * Default: ../logs/
  * Contains: run logs, error logs, telemetry CSVs
- Plots Directory: Generated graphs and charts
  * Default: ../plots/
  * Contains: PNG images of all plots
- Vizard Files Directory: Binary files for 3D visualization
  * Default: ../Vizfile/_VizFiles/
  * Contains: .bin files for Vizard playback
- Results Directory: DRL and ML outputs
  * Default: ../DRL/result/
  * Contains: training logs, model files, metrics
- Reset to Defaults: Restore standard directory paths

Binary Filename:
- Name for the Vizard visualization file
- Don't include extension (.bin added automatically)
- Default: spacecraft_viz
- Timestamp automatically added to prevent overwrites
- File saved in Vizfile/_VizFiles/ directory

Simulation Log:
- Real-time updates during simulation
- Shows configuration changes and progress
- Displays fault injection events
- Reports ML detection results
- Shows DRL task reassignment decisions
- Save log to file for record keeping
- Clear log to remove old messages
- Auto-scroll keeps latest messages visible
- Filter by message type: INFO, WARNING, ERROR

The simulation generates:
- Summary text files with configuration details
- Telemetry CSV files for each spacecraft
- Plot images showing spacecraft behavior:
  * Attitude errors (roll, pitch, yaw)
  * Reaction wheel speeds
  * Power consumption
  * Orbit trajectories
  * Target visibility
  * Cluster formations
  * Communication links
  * Fault detection confidence
  * Task reassignment results
- Binary file for Vizard 3D visualization
- ML detection results (JSON format)
- DRL reassignment plan (JSON format)
- Training logs (if DRL training occurred)

Log Message Types:
- INFO: Normal operations, status updates
- WARNING: Non-critical issues, suboptimal conditions
- ERROR: Critical problems, simulation failures
- DEBUG: Detailed technical information (if enabled)

Export Options:
- Export Configuration: Save current setup to JSON file
- Export Results: Save all results in ZIP archive
- Export Plots: Batch export all plots as PNG
- Export Report: Generate PDF with plots and summary
- Export Telemetry: Save all telemetry as CSV""",

        "results": """Results Tab
===========

The Results tab displays simulation outputs and analysis.

Plot List:
- Shows all generated plots from simulations
- Sorted by creation time (newest first)
- Filtered by type: Attitude, Orbit, Fault, Cluster, ML, DRL
- Click to view plot in main display
- Double-click to open in external viewer
- Right-click for options: Export, Delete, Rename

Plot Categories:
1. Attitude Plots:
   - Attitude error history (roll, pitch, yaw)
   - Rate error (angular velocity)
   - Control torques applied
   - Reaction wheel speeds

2. Orbit Plots:
   - 3D orbit trajectories
   - Ground track (lat/lon over time)
   - Altitude profile
   - Orbital elements over time

3. Fault Plots:
   - Fault injection timing
   - Before/after fault comparison
   - Wheel performance degradation
   - Recovery behavior

4. Cluster Plots (NEW):
   - Constellation overview with all satellites
   - Formation geometry for each cluster
   - Inter-satellite distance analysis
   - Communication link status
   - Relative positions over time

5. ML Detection Plots:
   - Anomaly scores per spacecraft
   - Detection confidence time series
   - Telemetry with fault annotations
   - Statistical analysis results
   - Method comparison charts

6. DRL Plots:
   - Training reward curves
   - Success rate over episodes
   - Task reassignment plans
   - System health metrics
   - State-action distributions

Plot Viewer:
- Displays selected plot with zoom/pan controls
- Navigation toolbar for plot interaction:
  * Home: Reset view
  * Back/Forward: Navigate zoom history
  * Pan: Click and drag to move
  * Zoom: Click and drag to select region
  * Configure: Adjust plot settings
  * Save: Export plot in various formats
- Mouse wheel zoom
- Right-click for additional options

Export Formats:
- PNG: Raster image (default, good for documents)
- PDF: Vector format (scalable, good for publications)
- SVG: Vector format (editable in Inkscape, Illustrator)
- JPG: Compressed raster (smaller file size)
- EPS: Vector format (for LaTeX)

File Information Panel:
- Filename and full path
- Creation timestamp
- File size
- Image dimensions
- Associated simulation run
- Analysis notes (editable)

Plot Comparison:
- Select multiple plots to compare side-by-side
- Overlay plots for direct comparison
- Sync zoom/pan across plots
- Export comparison as composite image

Search and Filter:
- Search by filename or keywords
- Filter by date range
- Filter by simulation run ID
- Filter by spacecraft name
- Filter by fault type
- Show only favorites

Batch Operations:
- Export all plots in category
- Delete old plots (before date)
- Rename with consistent pattern
- Add notes to multiple plots
- Create plot collections

Plot Refresh:
- Auto-refresh: Update list when new plots created
- Manual refresh button
- Shows "New" badge on recently created plots
- Notification when simulation generates plots

Using Results:
1. Run simulation to generate plots
2. Plots appear automatically in list
3. Select plot from list to view
4. Use toolbar to zoom/pan for details
5. Export plots for reports/presentations
6. Open results folder to access all files
7. View telemetry CSV files for raw data
8. Load Vizard .bin file for 3D playback

The plots help analyze:
- Fault impact on spacecraft performance
- Constellation geometry over time
- Formation maintenance quality
- Target visibility windows
- Communication link availability
- System behavior during failures
- ML detection accuracy and timing
- DRL reassignment effectiveness
- Overall mission success metrics

Advanced Analysis:
- Load telemetry CSV into Python/MATLAB for custom analysis
- Compare multiple simulation runs
- Statistical analysis of fault effects
- Correlation between faults and task reassignment
- Training effectiveness visualization
- System resilience metrics"""
    }

    return help_sections.get(section, help_sections["overview"])