from Basilisk.utilities import RigidBodyKinematics, macros
import numpy as np
import matplotlib.pyplot as plt

import Tools
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras.models import load_model

# Define PD Controller
def PDControllor(error, omega):
    Kp = 20.0
    Kd = 2.011
    torque = Kp * error - Kd * omega
    return np.array(torque)

def random_euler():
    phi = np.random.uniform(-np.pi, np.pi)
    theta = np.random.uniform(-np.pi, np.pi)
    psi = np.random.uniform(-np.pi, np.pi)
    return np.array([phi, theta, psi])


from Envs import BasiliskEnv
import Tools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Envs import BasiliskModel
from Basilisk.utilities import RigidBodyKinematics

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.preprocessing import MinMaxScaler

def PD_test():
    # Load the saved LSTM model
    model = load_model('anomaly_detection_mode.keras')
    print("Model loaded successfully")

    # Load reference data and fit the scaler on it (e.g., training data)
    training_data = pd.read_excel('simulation_log.xlsx')  # Reference dataset for scaling
    scaler = MinMaxScaler()
    scaler.fit(training_data.drop(columns=['Time']))

    ref_MRP = RigidBodyKinematics.euler1232MRP(random_euler())
    env = BasiliskEnv(faulty=True, torque_mode="wheel")
    b_model = BasiliskModel(I=[0.025, 0, 0, 0, 0.05, 0, 0, 0, 0.065], ref_MRP=ref_MRP, torque_mode="wheel")

    log_data = []
    reconstruction_errors = []
    time_intervals = []
    b_MRP_history = []
    b_omega_history = []
    b_error_MRP_history = []
    torque_history = [np.array([0, 0, 0, 0])]
    dataRW = []
    dataRwPower = []

    # Run the simulation loop
    for i in range(200):
        torque = PDControllor(b_model.cur_error_MRP, b_model.cur_omega)
        if b_model.faulty and b_model.step_count > b_model.fault_time and b_model.wheel_num != -1:
            torque[b_model.wheel_num] = 0
        torque = np.append(torque, 0)
        torque_history.append(torque)

        cur_MRP, cur_error_MRP, cur_error_angle, cur_omega, cur_omega_dot = b_model.step(macros.sec2nano((i+1)*0.08), torque)
        b_MRP_history.append(cur_MRP)
        b_omega_history.append(cur_omega)
        b_error_MRP_history.append(cur_error_MRP)
        # Log the data for each step
        step_data = {
            "Time": (i+1),
            "Torque_RW0": torque[0],
            "Torque_RW1": torque[1],
            "Torque_RW2": torque[2],
            "Torque_RW3": torque[3],
            "MRP_x": cur_MRP[0],
            "MRP_y": cur_MRP[1],
            "MRP_z": cur_MRP[2],
            "Error_MRP_x": cur_error_MRP[0],
            "Error_MRP_y": cur_error_MRP[1],
            "Error_MRP_z": cur_error_MRP[2],
            "Omega_x": cur_omega[0],
            "Omega_y": cur_omega[1],
            "Omega_z": cur_omega[2],
            "RW_Power_0": 0.0,
            "RW_Power_1": 0.0,
            "RW_Power_2": 0.0,
            "RW_Power_3": 0.0
        }
        log_data.append(step_data)

        # Process data after 10 timesteps
        if len(log_data) >= 10:
            X = np.array([scaler.transform(pd.DataFrame(log_data[-10:]).drop(columns=['Time']))])
            X = X.reshape((1, 10, X.shape[2]))  # Reshape for LSTM input
            X_pred = model.predict(X)
            reconstruction_error = np.mean(np.abs(X_pred - X), axis=2)
            reconstruction_error = np.mean(reconstruction_error, axis=1)[0]
            reconstruction_errors.append(reconstruction_error)
            time_intervals.append(step_data["Time"])

            # Define the threshold for anomaly detection
            threshold = 0.185  # Adjust threshold if necessary
           
        #    threshold_upper = 0.17 
         #   threshold_lower = 0.247 

# Identify anomalies where the reconstruction error is either greater than the upper threshold or less than the lower threshold 
#anomalies = (reconstruction_error > threshold_upper) | (reconstruction_error < threshold_lower)

            if reconstruction_error > threshold:
                print(f"Anomaly detected at time {step_data['Time']} with reconstruction error: {reconstruction_error}")
    b_MRP_history = np.array(b_MRP_history)
    b_omega_history = np.array(b_omega_history)
    b_error_MRP_history = np.array(b_error_MRP_history)
    torque_history = np.array(torque_history)
    # Plot reconstruction errors
    plt.figure(figsize=(10, 6))
    plt.plot(time_intervals, reconstruction_errors, label='Reconstruction Error')
    anomalies = np.array(reconstruction_errors) > threshold
    #anomalies = (np.array(reconstruction_errors) > threshold_upper) | (np.array(reconstruction_errors) < threshold_lower)
    plt.scatter(np.array(time_intervals)[anomalies], np.array(reconstruction_errors)[anomalies], color='red', label='Anomalies')
    plt.title('Real-Time Anomaly Detection')
    plt.xlabel('Time (s)')
    plt.ylabel('Reconstruction Error')
    plt.legend()
    plt.show()
    fig, axes = plt.subplots(2, 2)
    fig.suptitle(f"Fault time: {b_model.fault_time}, Wheel num: {b_model.wheel_num}")

    ref_MRP_array = np.array(ref_MRP)
    time_steps = len(b_MRP_history)
    ref_line = np.tile(ref_MRP_array, (time_steps, 1))
    
    ax = axes[0, 0]
    ax.plot(b_MRP_history, label=["x-actual", "y-actual", "z-actual"])
    ax.plot(ref_line, '--', label=["x-Desired", "y-Desired", "z-Desired"])
    ax.set_title("Desired Vs Actual Attitude")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(b_error_MRP_history, label=["error_x", "error_y", "error_z"])
    ax.set_title("Error History")
    ax.legend()

    error_angle_history = []
    for mrp in b_error_MRP_history:
        error_angle_history.append(4 * np.arctan(np.linalg.norm(mrp)))

    ax = axes[1, 0]
    ax.plot(b_omega_history, label=["omega_x-actual", "omega_y-actual", "omega_z-actual"])
    
    error_angle_history_array = np.array(error_angle_history)
    adjusted_omega_history = b_omega_history - error_angle_history_array[:, np.newaxis]
    ax.plot(adjusted_omega_history, '--', label=["omega_x-Desired", "omega_y-Desired", "omega_z-Desired"])
    ax.set_title("Angular Velocity")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(error_angle_history, label="error")
    ax.set_title("Angle Error")
    ax.legend()

    plt.show()

    # Plot Torque History
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))

    # RW0
    ax = axes[0, 0]
    ax.plot(torque_history[:, 0], label=["RW0"], color='red')
    ax.set_title("Torque History RW0")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Torque (N·m)")
    ax.legend()

    # RW1
    ax = axes[0, 1]
    ax.plot(torque_history[:, 1], label=["RW1"], color='blue')
    ax.set_title("Torque History RW1")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Torque (N·m)")
    ax.legend()

    # RW2
    ax = axes[1, 0]
    ax.plot(torque_history[:, 2], label=["RW2"], color='green')
    ax.set_title("Torque History RW2")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Torque (N·m)")
    ax.legend()

    # RW3
    ax = axes[1, 1]
    ax.plot(torque_history[:, 3], label=["RW3"], color='yellow')
    ax.set_title("Torque History RW3")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Torque (N·m)")
    ax.legend()


    fig.suptitle(f"Fault time: {b_model.fault_time}, Wheel num: {b_model.wheel_num}")
    plt.show()

if __name__ == "__main__":
    PD_test()
