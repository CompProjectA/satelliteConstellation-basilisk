from Basilisk.utilities import RigidBodyKinematics, macros
import numpy as np
import matplotlib.pyplot as plt

import Tools
import pandas as pd


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
from Envs import BasiliskModel, DynamicModel
from Basilisk.utilities import RigidBodyKinematics

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def PD_test():
    ref_MRP = RigidBodyKinematics.euler1232MRP(random_euler())
    env = BasiliskEnv(faulty=True, torque_mode="wheel")
    b_model = BasiliskModel(I=[0.025, 0, 0, 0, 0.05, 0, 0, 0, 0.065], ref_MRP=ref_MRP, torque_mode="wheel")

    b_MRP_history = []
    b_omega_history = []
    b_error_MRP_history = []
    torque_history = [np.array([0, 0, 0, 0])]
    dataRW = []
    dataRwPower = []

    # Run the simulation loop
    for i in range(3000):
        torque = PDControllor(b_model.cur_error_MRP, b_model.cur_omega)
        if b_model.faulty and b_model.step_count > b_model.fault_time and b_model.wheel_num != -1:
            torque[b_model.wheel_num] = 0
        torque = np.append(torque, 0)
        torque_history.append(torque)
        cur_MRP, cur_error_MRP, cur_error_angle, cur_omega, cur_omega_dot = b_model.step(macros.sec2nano((i+1)*0.08), torque)
        b_MRP_history.append(cur_MRP)
        b_omega_history.append(cur_omega)
        b_error_MRP_history.append(cur_error_MRP)

    for c in range(0, 4):
        dataRW.append(b_model.rwOutLog[c].u_current)
        dataRwPower.append(b_model.rwPowLog[c].netPower)

    # Convert history lists to numpy arrays
    b_MRP_history = np.array(b_MRP_history)
    b_omega_history = np.array(b_omega_history)
    b_error_MRP_history = np.array(b_error_MRP_history)
    torque_history = np.array(torque_history)
    dataRwPower = np.array(dataRwPower)

    # Adjust the length of arrays with 2001 elements by removing the last element
    dataRwPower = dataRwPower[:, :-1]
    torque_history = torque_history[:-1]
    
    # Continue with the print statements to check the lengths
    print('b_model.timeData', len(b_model.timeData))
    print('dataRwPower', len(dataRwPower[0]))
    print('torque_history', len(torque_history[:, 0]))
    print('torque_history', len(torque_history[:, 1]))
    print('torque_history', len(torque_history[:, 2]))
    print('torque_history', len(torque_history[:, 3]))
    print('b_MRP_history', len(b_MRP_history[:, 0]))
    print('b_MRP_history', len(b_MRP_history[:, 1]))
    print('b_MRP_history', len(b_MRP_history[:, 2]))
    print('b_error_MRP_history', len(b_error_MRP_history[:, 0]))
    print('b_error_MRP_history', len(b_error_MRP_history[:, 1]))
    print('b_error_MRP_history', len(b_error_MRP_history[:, 2]))
    print('b_omega_history', len(b_omega_history[:, 0]))
    print('b_omega_history', len(b_omega_history[:, 1]))
    print('b_omega_history', len(b_omega_history[:, 2]))
    print('dataRwPower', len(dataRwPower[ 0]))
    print('dataRwPower', len(dataRwPower[ 1]))
    print('dataRwPower', len(dataRwPower[ 2]))
    print('dataRwPower', len(dataRwPower[ 3]))
   
    timeData = (b_model.timeData * 1e-7)
    timeData=timeData[:-1]
    # Create a dictionary for storing all the data to be written to Excel
    data_dict = {
        "Time ": timeData,
        "Torque_RW0": torque_history[:, 0],
        "Torque_RW1": torque_history[:, 1],
        "Torque_RW2": torque_history[:, 2],
        "Torque_RW3": torque_history[:, 3],
        "MRP_x": b_MRP_history[:, 0],
        "MRP_y": b_MRP_history[:, 1],
        "MRP_z": b_MRP_history[:, 2],
        "Error_MRP_x": b_error_MRP_history[:, 0],
        "Error_MRP_y": b_error_MRP_history[:, 1],
        "Error_MRP_z": b_error_MRP_history[:, 2],
        "Omega_x": b_omega_history[:, 0],
        "Omega_y": b_omega_history[:, 1],
        "Omega_z": b_omega_history[:, 2],
        "RW_Power_0": dataRwPower[0],
        "RW_Power_1": dataRwPower[1],
        "RW_Power_2": dataRwPower[2],
        "RW_Power_3": dataRwPower[3]
    }

    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(data_dict)
    #filename = f"simulation_results_{count}.xlsx"
    # Write DataFrame to an Excel file
    #df.to_excel(filename, index=False)

    # Plot RW Power
    for idx in range(4):
        plt.plot(timeData, dataRwPower[idx], label='$p_{rw,' + str(idx) + '}$')
    
    plt.legend(loc='lower right')
    plt.xlabel('Time')
    plt.ylabel('RW Power (W)')
    #plt.show()

    # Continue with your original plotting...
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
