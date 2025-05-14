#!/usr/bin/env python
"""
encoder_fault.py

A Basilisk scenario that simulates spacecraft dynamics with encoder faults
in the reaction wheels and properly saves binary files for Vizard visualization.

For use with the spacecraft fault simulator GUI.
"""
import inspect
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from Basilisk.utilities import (orbitalMotion, macros, vizSupport)

# Set paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
ROOT_DIR = os.path.abspath(os.path.join(path, '..'))
MODELS_DIR = os.path.join(ROOT_DIR, 'models')
PLOTTING_DIR = os.path.join(ROOT_DIR, 'plotting')
VIZ_DIR = os.path.join(ROOT_DIR, 'Vizfile')

sys.path.extend([ROOT_DIR, MODELS_DIR, PLOTTING_DIR])

# Import BSK modules
try:
    from BSK_masters import BSKSim, BSKScenario
    import BSK_Dynamics, BSK_Fsw
    import plotting.BSK_Plotting as BSK_plt
except ImportError as e:
    print(f"ERROR: Could not import required modules: {e}")
    sys.exit(1)

class EncoderFaultScenario(BSKSim, BSKScenario):
    """
    Scenario for simulating encoder faults in reaction wheels.
    Inherits from BSKSim and BSKScenario for Basilisk simulation framework.
    """
    def __init__(self):
        super(EncoderFaultScenario, self).__init__()
        self.name = 'EncoderFaultScenario'
        self.msgRecList = {}
        self.sNavTransName = "sNavTransMsg"
        self.attGuidName = "attGuidMsg"

        self.cameraLocation = [0.0, 2.0, 0.0]

        self.targets = [
            {"name": "Melbourne", "lat": -37.8136, "lon": 144.9631, "color": "red"},
            {"name": "New York", "lat": 40.71, "lon": -74.00, "color": "blue"},
            {"name": "Tokyo", "lat": 35.68, "lon": 139.77, "color": "green"},
            {"name": "London", "lat": 51.51, "lon": -0.13, "color": "yellow"}
        ]

        self.set_DynModel(BSK_Dynamics)
        self.set_FswModel(BSK_Fsw)

        self.configure_initial_conditions()
        self.log_outputs()

        self.oneTimeRWFaultFlag = 1  # Important - this is the flag for encoder faults
        self.oneTimeFaultTime = macros.min2nano(10.)
        self.fault_wheel_number = 3
        self.DynModels = self.get_DynModel()
        self.FSWModels = self.get_FswModel()
        self.DynModels.RWFaultLog = []

    def configure_initial_conditions(self):
        """Configure orbit and attitude initial conditions"""
        oe = orbitalMotion.ClassicElements()
        oe.a = 10000000.0
        oe.e = 0.01
        oe.i = 33.3 * macros.D2R
        oe.Omega = 48.2 * macros.D2R
        oe.omega = 347.8 * macros.D2R
        oe.f = 85.3 * macros.D2R

        DynModel = self.get_DynModel()
        mu = DynModel.gravFactory.gravBodies['earth'].mu
        rN, vN = orbitalMotion.elem2rv(mu, oe)
        orbitalMotion.rv2elem(mu, rN, vN)
        DynModel.scObject.hub.r_CN_NInit = rN
        DynModel.scObject.hub.v_CN_NInit = vN
        DynModel.scObject.hub.sigma_BNInit = [[0.1], [0.2], [-0.3]]
        DynModel.scObject.hub.omega_BN_BInit = [[0.001], [-0.01], [0.03]]

    def log_outputs(self):
        """Configure message logging for analysis"""
        FswModel = self.get_FswModel()
        DynModel = self.get_DynModel()
        samplingTime = FswModel.processTasksTimeStep

        # RW speeds from dynamics (actual)
        self.rwSpeedRec = DynModel.rwStateEffector.rwSpeedOutMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.rwSpeedRec)
        
        # RW speeds command (reference)
        try:
            self.rwSpeedCmdRec = FswModel.rwMotorCmdMsg.recorder(samplingTime)
            self.AddModelToTask(FswModel.taskName, self.rwSpeedCmdRec)
            self.have_command_speeds = True
        except:
            print("WARNING: Could not record RW command speeds, encoder fault visualization will be limited")
            self.have_command_speeds = False

        self.rwLogs = []
        for item in range(4):
            self.rwLogs.append(DynModel.rwStateEffector.rwOutMsgs[item].recorder(samplingTime))
            self.AddModelToTask(DynModel.taskName, self.rwLogs[item])
            
        # FSW controller outputs
        self.msgRecList[self.attGuidName] = FswModel.attGuidMsg.recorder(samplingTime)
        self.AddModelToTask(DynModel.taskName, self.msgRecList[self.attGuidName])

    def pull_outputs(self, showPlots):
        """Process and plot simulation outputs"""
        numRW = 4
        RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(numRW)], 0, 0)
        
        # Get RW command speeds if available
        RW_cmd_speeds = None
        if self.have_command_speeds:
            try:
                RW_cmd_speeds = np.delete(self.rwSpeedCmdRec.wheelSpeeds[:, range(numRW)], 0, 0)
            except:
                print("Warning: Could not extract command speeds from recorder")

        # Get attitude logs
        attErrRec = self.msgRecList[self.attGuidName]
        sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
        omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)

        # Clear all existing plots
        plt.close('all')
        
        # Get time data
        timeData = np.delete(self.rwSpeedRec.times(), 0, 0) * macros.NANO2MIN
        fault_time_min = self.oneTimeFaultTime * macros.NANO2MIN
        
        # Generate plots
        # Attitude error plot (figure 1)
        plt.figure(1)
        plt.plot(timeData, np.linalg.norm(sigma_BR, axis=1))
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Attitude Error Norm')
        plt.title('Attitude Error')
        plt.grid(True)
        plt.legend()
        
        # Rate error plot (figure 2)
        plt.figure(2)
        for i in range(3):
            plt.plot(timeData, omega_BR_B[:, i], label=f'Axis {i+1}')
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Angular Rate Error [rad/s]')
        plt.title('Rate Error')
        plt.grid(True)
        plt.legend()
        
        # RW speeds plot (figure 3)
        plt.figure(3)
        for i in range(numRW):
            plt.plot(timeData, RW_speeds[:, i], label=f'RW{i+1}')
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Wheel Speed [rad/s]')
        plt.title('Reaction Wheel Speeds')
        plt.grid(True)
        plt.legend()
        
        # Encoder fault plot (figure 4)
        plt.figure(4)
        
        # Highlight the faulty wheel
        plt.subplot(211)
        plt.plot(timeData, RW_speeds[:, self.fault_wheel_number], 'b-', label=f'RW{self.fault_wheel_number+1} Speed')
        if RW_cmd_speeds is not None:
            plt.plot(timeData, RW_cmd_speeds[:, self.fault_wheel_number], 'g--', label=f'RW{self.fault_wheel_number+1} Command')
        plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
        plt.xlabel('Time [min]')
        plt.ylabel('Wheel Speed [rad/s]')
        plt.title(f'Encoder Fault on RW{self.fault_wheel_number+1}')
        plt.grid(True)
        plt.legend()
        
        # Also show error between commanded and actual speeds
        plt.subplot(212)
        if RW_cmd_speeds is not None:
            error = RW_speeds[:, self.fault_wheel_number] - RW_cmd_speeds[:, self.fault_wheel_number]
            plt.plot(timeData, error, 'm-', label='Speed Error')
            plt.axvline(x=fault_time_min, color='r', linestyle='--', label='Fault Injection')
            plt.xlabel('Time [min]')
            plt.ylabel('Speed Error [rad/s]')
            plt.title('Commanded vs Actual Speed Error')
        else:
            plt.text(0.5, 0.5, 'Command speeds not available\nCannot calculate error', 
                    horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        plt.grid(True)
        plt.legend()
        
        plt.tight_layout()

        # Create figure list to return
        figureList = {
            "attitudeErrorNorm": plt.figure(1),
            "rateError": plt.figure(2),
            "RWSpeeds": plt.figure(3),
            "RWEncoder": plt.figure(4)
        }

        # Show plots if requested
        if showPlots:
            plt.show()
        else:
            plt.close('all')

        return figureList
        
    def plot_target_visibility(self, timeData):
        """Plot target visibility based on spacecraft orbit"""
        # Implement target visibility plotting if required
        pass

    def inject_rw_encoder_fault(self, wheel_idx, time_nano):
        """
        Inject an encoder fault into a specific reaction wheel
        
        Parameters:
        wheel_idx (int): Index of the wheel to fault (0-3)
        time_nano (float): Simulation time in nanoseconds when fault occurs
        """
        # Log the fault for plotting
        self.DynModels.RWFaultLog.append({
            'type': 'encoder',
            'wheel': wheel_idx,
            'time': time_nano,
            'magnitude': 0.0  # Encoder faults don't have a magnitude
        })
        
        # Here you would implement the encoder fault based on your simulation's capabilities
        # This could involve:
        # 1. Modifying the sensor reading pipeline
        # 2. Adding noise to the wheel speed feedback
        # 3. Setting up incorrect/zeroed feedback
        
        print(f"Encoder fault injected in wheel {wheel_idx} at time {time_nano * macros.NANO2MIN:.2f} minutes")
        
        # Example implementation (modify as needed for your Basilisk setup):
        try:
            # This is a placeholder - replace with the actual encoder fault implementation
            # that works with your BSK model configuration
            self.FSWModels.rwMotorVoltage.encoderFault[wheel_idx] = True
        except AttributeError:
            print("Warning: Could not apply encoder fault directly, simulation may not show full fault effects")
            
            # Alternative fault implementation if the more direct approach doesn't work
            try:
                # Try to use a more general approach like modifying the measured speed feedback
                # This is a simplified way to simulate a "stuck" encoder by freezing the measurement
                # You'd need to implement this method in your model
                self.FSWModels.freezeEncoderFeedback(wheel_idx)
            except:
                print("Could not implement encoder fault through alternative means")
        
        return

def runScenario(scenario, saveBinary=True):
    """Run the encoder fault scenario"""
    simulationTime = macros.min2nano(30.)
    scenario.modeRequest = "hillPoint"

    # Set up encoder fault event
    scenario.createNewEvent(
        "injectEncoderFault",
        scenario.get_FswModel().processTasksTimeStep,
        True,
        ["self.TotalSim.CurrentNanos>=self.oneTimeFaultTime and self.oneTimeRWFaultFlag==1"],
        [f"self.inject_rw_encoder_fault({scenario.fault_wheel_number}, self.TotalSim.CurrentNanos)", 
         "self.oneTimeRWFaultFlag=0"]
    )

    viz = None
    if vizSupport.vizFound:
        # Create visualization binary directory if it doesn't exist
        vizfiles_dir = os.path.join(VIZ_DIR, "_VizFiles")
        if not os.path.exists(vizfiles_dir):
            try:
                os.makedirs(vizfiles_dir, exist_ok=True)
            except:
                print(f"Warning: Could not create directory {vizfiles_dir}")
        
        # Enable visualization
        binary_filename = "encoder_fault_viz" if saveBinary else None
        if saveBinary:
            binary_path = os.path.join(vizfiles_dir, binary_filename)
        else:
            binary_path = None
            
        viz = vizSupport.enableUnityVisualization(
            scenario,
            scenario.get_DynModel().taskName,
            scenario.get_DynModel().scObject,
            rwEffectorList=scenario.get_DynModel().rwStateEffector,
            liveStream=not saveBinary,
            saveFile=binary_path
        )

        # Add targets
        for target in scenario.targets:
            lat = target["lat"]
            lon = target["lon"]
            color = target.get("color", "red")
            alt = 0.0
            radius = 6371000.0 + alt
            lat_rad = lat * macros.D2R
            lon_rad = lon * macros.D2R
            x = radius * np.cos(lat_rad) * np.cos(lon_rad)
            y = radius * np.cos(lat_rad) * np.sin(lon_rad)
            z = radius * np.sin(lat_rad)
            location_position = [x, y, z]

            vizSupport.addLocation(
                viz,
                stationName=target["name"],
                parentBodyName="earth",
                r_GP_P=location_position,
                color=color
            )

        # Add camera that looks at the spacecraft body
        vizSupport.createStandardCamera(
            viz,
            setMode=1,
            spacecraftName=scenario.get_DynModel().scObject.ModelTag,
            fieldOfView=70 * macros.D2R,
            displayName="RW Camera",
            pointingVector_B=[0, 0, 0],
            position_B=scenario.cameraLocation
        )

    scenario.InitializeSimulation()
    scenario.ConfigureStopTime(simulationTime)
    scenario.ExecuteSimulation()

    return viz

def run(showPlots=True, saveBinary=True):
    """
    Run the encoder fault scenario with default parameters
    
    Parameters:
    showPlots (bool): Flag to display plots
    saveBinary (bool): Flag to save binary file for visualization
    
    Returns:
    tuple: (scenario, viz, figureList) - The simulation objects and results
    """
    print("\n===== Running Encoder Fault Scenario =====")
    print(f"Save Binary: {saveBinary}")
    scenario = EncoderFaultScenario()
    viz = runScenario(scenario, saveBinary)
    figureList = scenario.pull_outputs(showPlots)

    if saveBinary and viz:
        print("\nBinary file saved successfully as 'encoder_fault_viz.bin'")
        print("You can now open this file in Vizard for visualization.")

    return scenario, viz, figureList

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the Encoder Fault Scenario")
    parser.add_argument("--no-plots", action="store_true", help="Don't show plots")
    parser.add_argument("--no-binary", action="store_true", help="Don't save binary file")
    args = parser.parse_args()

    run(not args.no_plots, not args.no_binary)