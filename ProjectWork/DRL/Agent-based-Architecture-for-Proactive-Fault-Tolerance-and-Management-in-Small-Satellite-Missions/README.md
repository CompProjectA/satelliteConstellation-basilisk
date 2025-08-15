# Fault-Tolerance and Management architecture for Small Satellite Missions based on Basilisk simulator

testing

Designing a Fault-Tolerance and Management architecture for Small Satellite Missions by considering several functions including anomaly detection, Explainable Artificial Intelligence (XAI), decision making & recovery, coordinated adaptation via inter-satellite links, to enhance the resiliency and autonomy of satellite constellation, is the main goal of this project.
This work has been supported by the SmartSat CRC (under Project No. P2-52), whose activities are funded by the Australian Government’s CRC Program.

## Getting Started
At first step, [Basilisk](https://avslab.github.io/basilisk/), [Vizard](https://avslab.github.io/basilisk/Vizard/Vizard.html) and [BSK-RL](https://avslab.github.io/bsk_rl/) (Basilisk + Reinforcement Learning) should be installed on your device.

1- Basilisk, or BSK for short, is a software framework capable of both faster-than realtime spacecraft simulations, including repeatable Monte-Carlo simulation options, as well as providing real-time options for hardware-in-the-loop simulations.

Installation instructions, examples, and documentation can be found on the [Basilisk website](https://avslab.github.io/basilisk/) and [GitHub Page](https://github.com/AVSLab/basilisk).

2- The Vizard Unity-based Basilisk visualization is able to display in a three-dimensional view the Basilisk simulation data.

For installation, download Vizard from [website](https://avslab.github.io/basilisk/Vizard/VizardDownload.html), then extract the .zip file to a specified folder. Finally, inside the extracted Vizard folder, go to Windows\Vizard and run the Vizard.exe file

3- BSK-RL is a Python package for constructing Gymnasium environments for spacecraft tasking problems. It is built on top of Basilisk, making the simulation environments high-fidelity and computationally efficient.

Installation instructions, examples, and documentation can be found on the [BSK-RL website](https://avslab.github.io/bsk_rl/) and [GitHub Page](https://github.com/AVSLab/bsk_rl).

Note 1: BSK-RL should be cloned and installed as a folder in Basilisk. Basilisk v2.2.0, Vizard v2.1.5 and BSK-RL v1.0.1 is used for current simulation.

Note 2: It is recommended to install Visual Studio and dring installation, select the “Desktop development with C++” and tick the box.

## Prerequsites
Integration requires the following Python packages:
TensorFlow 2.17.0, Pandas 2.2.1, NumPy 1.26.4, Matplotlib 3.8.3 and sklearn 1.5.2 (ensure all mentioned packages are installed (by using **pip list** command))

## Repository Structure and Explanation

Envs.py: Represents the use case scenario with one satellite that orbit the earth and unresponsive reaction wheel is injected in predefined time that uses PPO agent for learning and decision making

PD.py: Represents the classic PD control logic for spacecraft attitude control.

PPO.py: Contains the PPO DRL agent used for intelligent decision-making.

Tools.py: Utility functions data plotting.

anomaly_detection_model.keras: Trained Keras model for real-time fault detection.

Integration.py: integrates the trained anomaly detection model with decision making & recovery. The log data is loaded and checks to detect the anomaly via loaded anomaly detection model. Whenever the anomaly is detected a true flag will be activated and sent to decision making agent. Finally, agent will take the recovey action.

*The excel files represent the logged data during the simulation.


These files should be copied in the bsk_rl folder and Integration.py is the main file that should be executed.


