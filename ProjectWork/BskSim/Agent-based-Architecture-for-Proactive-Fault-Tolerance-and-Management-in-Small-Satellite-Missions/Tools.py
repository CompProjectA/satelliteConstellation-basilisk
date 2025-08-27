import os

os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
from Basilisk.utilities import macros, unitTestSupport
import numpy as np
import matplotlib.pyplot as plt
from Envs import BasiliskEnv
from stable_baselines3.common.env_util import make_vec_env
from stable_baselines3.common.vec_env import DummyVecEnv, VecMonitor, SubprocVecEnv





def eulerXYZ_to_quat(euler: np.ndarray) -> np.ndarray:
    c1 = np.cos(euler[0] / 2)
    s1 = np.sin(euler[0] / 2)
    c2 = np.cos(euler[1] / 2)
    s2 = np.sin(euler[1] / 2)
    c3 = np.cos(euler[2] / 2)
    s3 = np.sin(euler[2] / 2)

    q = np.array([c1 * c2 * c3 + s1 * s2 * s3, s1 * c2 * c3 - c1 * s2 * s3, c1 * s2 * c3 + s1 * c2 * s3, c1 * c2 * s3 - s1 * s2 * c3])
    return q / np.linalg.norm(q)


def random_quar_ref() -> np.ndarray:
    psi = np.random.uniform(-np.pi, np.pi)  # yaw，偏航角，绕z轴
    theta = np.random.uniform(-np.pi / 2, np.pi / 2)  # pitch，俯仰角，绕y轴
    phi = np.random.uniform(-np.pi, np.pi)  # roll，横滚角，绕x轴
    # print(phi, theta, psi)
    ref = eulerXYZ_to_quat(np.array([phi, theta, psi]))
    ref = ref / np.linalg.norm(ref)
    return ref


def quat_mul(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return np.array(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ]
    )


def quat_error(q_cur: np.ndarray, q_des: np.ndarray) -> np.ndarray:
    q_cur_inv = np.array([q_cur[0], -q_cur[1], -q_cur[2], -q_cur[3]]) / np.linalg.norm(q_cur)
    q_cur_inv = q_cur_inv / np.linalg.norm(q_cur_inv)
    return quat_mul(q_des, q_cur_inv)


def mrp_to_quaternion(mrp: np.ndarray) -> np.ndarray:
    mrp_squared_norm = np.dot(mrp, mrp)
    scalar_part = (1 - mrp_squared_norm) / (1 + mrp_squared_norm)
    vector_part = (2 * mrp) / (1 + mrp_squared_norm)
    quaternion = np.hstack((scalar_part, vector_part))

    return quaternion


def random_euler():
    phi = np.random.uniform(-np.pi, np.pi)
    theta = np.random.uniform(-np.pi, np.pi)
    psi = np.random.uniform(-np.pi, np.pi)

    return np.array([phi, theta, psi])

import pandas as pd

def plot_history(dynamic_model, ref_MRP, torque_history, env):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    
    # Attitude Data (MRP)
    mrp_history = np.array(dynamic_model.MRP_history)
    ref_MRP_array = np.array(ref_MRP)
    time_steps = len(dynamic_model.MRP_history)
    ref_line = np.tile(ref_MRP_array, (time_steps, 1))
    
    ax = axes[0, 0]
    ax.plot(mrp_history, label=["x-actual", "y-actual", "z-actual"])
    ax.plot(ref_line, '--', label=["x-Desired", "y-Desired", "z-Desired"])  # Dashed lines for reference MRP
    ax.set_title("Desired vs Actual Attitude")
    ax.legend()

    # Error Data (MRP)
    error_mrp_history = np.array(dynamic_model.error_MRP_history)
    ax = axes[0, 1]
    ax.plot(error_mrp_history, label=["error_x", "error_y", "error_z"])
    ax.set_title("Attitude Error")
    ax.legend()

    # Angular Velocity Data
    omega_history = np.array(dynamic_model.omega_history)
    ax = axes[1, 0]
    ax.plot(omega_history, label=["omega_x-actual", "omega_y-actual", "omega_z-actual"])
    ax.set_title("Angular Velocity")
    
    error_angle_history_array = np.array(dynamic_model.error_angle_history)
    adjusted_omega_history = omega_history - error_angle_history_array[:, np.newaxis]
    ax.plot(adjusted_omega_history, '--', label=["omega_x-Desired", "omega_y-Desired", "omega_z-Desired"])
    ax.legend()

    # Error Angle Data
    ax = axes[1, 1]
    ax.plot(error_angle_history_array, label="error")
    ax.set_title("Error Angle")
    ax.legend()

    fig.suptitle(f"fault time: {env.fault_time}, wheel num: {env.wheel_num}")
    plt.show()

    # Torque History Data
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    
    torque_history = np.array(torque_history)

    ax = axes[0, 0]
    ax.plot(torque_history[:,0], label=["RW0"], color='red')
    ax.set_title("Torque History RW0")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(torque_history[:,1], label=["RW1"], color='blue')
    ax.set_title("Torque History RW1")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(torque_history[:,2], label=["RW2"], color='green')
    ax.set_title("Torque History RW2")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(torque_history[:,3], label=["RW3"], color='orange')
    ax.set_title("Torque History RW3")
    ax.legend()

    fig.suptitle(f"fault time: {env.fault_time}, wheel num: {env.wheel_num}")
    plt.show()
    file_path = 'plot_data.xlsx'
    mode = 'a' if os.path.exists(file_path) else 'w'
    # Log data to Excel
    with pd.ExcelWriter(file_path, mode=mode, engine='openpyxl') as writer:
        pd.DataFrame(mrp_history, columns=["x-actual", "y-actual", "z-actual"]).to_excel(writer, sheet_name="MRP History", index=False)
        pd.DataFrame(ref_line, columns=["x-Desired", "y-Desired", "z-Desired"]).to_excel(writer, sheet_name="Desired MRP", index=False)
        pd.DataFrame(error_mrp_history, columns=["error_x", "error_y", "error_z"]).to_excel(writer, sheet_name="Attitude Error", index=False)
        pd.DataFrame(omega_history, columns=["omega_x-actual", "omega_y-actual", "omega_z-actual"]).to_excel(writer, sheet_name="Angular Velocity", index=False)
        pd.DataFrame(adjusted_omega_history, columns=["omega_x-Desired", "omega_y-Desired", "omega_z-Desired"]).to_excel(writer, sheet_name="Desired Angular Velocity", index=False)
        pd.DataFrame(error_angle_history_array, columns=["Error Angle"]).to_excel(writer, sheet_name="Error Angle", index=False)
        pd.DataFrame(torque_history, columns=["RW0", "RW1", "RW2", "RW3"]).to_excel(writer, sheet_name="Torque History", index=False)

    return fig, axes





def plot_individual(dynamic_model, ref_MRP, torque_history, env):
    # Attitude Data (MRP)
    mrp_history = np.array(dynamic_model.MRP_history)
    ref_MRP_array = np.array(ref_MRP)
    time_steps = len(dynamic_model.MRP_history)
    ref_line = np.tile(ref_MRP_array, (time_steps, 1))
    
    # Plot Desired vs Actual Attitude (MRP)
    plt.figure(figsize=(5, 3))
    plt.plot(mrp_history, label=["x-actual", "y-actual", "z-actual"])
    plt.plot(ref_line, '--', label=["x-Desired", "y-Desired", "z-Desired"])  # Dashed lines for reference MRP
    plt.title("Desired vs Actual Attitude")
    plt.xlabel("Time (s)")
    plt.ylabel("Attitude (MRP)")
    plt.legend()
    plt.show()
    
    # Plot Attitude Error (MRP)
    error_mrp_history = np.array(dynamic_model.error_MRP_history)
    plt.figure(figsize=(5, 3))
    plt.plot(error_mrp_history, label=["error_x", "error_y", "error_z"])
    plt.title("Attitude Error")
    plt.xlabel("Time (s)")
    plt.ylabel("Attitude Error (deg)")
    plt.legend()
    plt.show()

    # Plot Angular Velocity
    omega_history = np.array(dynamic_model.omega_history)
    error_angle_history_array = np.array(dynamic_model.error_angle_history)
    adjusted_omega_history = omega_history - error_angle_history_array[:, np.newaxis]

    plt.figure(figsize=(5, 3))
    plt.plot(omega_history, label=["omega_x-actual", "omega_y-actual", "omega_z-actual"])
    plt.plot(adjusted_omega_history, '--', label=["omega_x-Desired", "omega_y-Desired", "omega_z-Desired"])
    plt.title("Angular Velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Angular Velocity (deg/s)")
    plt.legend()
    plt.show()

    # Plot Error Angle
    plt.figure(figsize=(5, 3))
    plt.plot(error_angle_history_array, label="error")
    plt.title("Error Hitory")
    plt.xlabel("Time (s)")
    plt.ylabel("Angular Error (deg)")
    plt.legend()
    plt.show()

    # Plot Torque History for each RW
    torque_history = np.array(torque_history)

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


    fig.suptitle(f"Fault time: {env.fault_time}, Wheel num: {env.wheel_num}")
    plt.show()
    
    # Optional: save to Excel (if needed)






def plot_rw_power0(timeData, dataRwPower, numRW):
    """Plot the RW actual motor torques and log data to Excel."""
    print(f"Original dataRwPower slice: {dataRwPower}")  # Print first few elements for inspection

    # Convert time to minutes and slice time data
    timeData_in_minutes = (timeData * 1e-7)[:2000]
    
    # Ensure consistent lengths by flattening, padding, or trimming each RW's data
    dataRwPower_trimmed = []
    for idx in range(numRW):
        rw_slice = dataRwPower[idx*2000:(idx+1)*2000]
        
        # Flatten the rw_slice in case it's a list of arrays
        rw_slice_flattened = np.hstack(rw_slice) if isinstance(rw_slice[0], np.ndarray) else rw_slice
        
        # If rw_slice is shorter than 2000, pad with the last value
        if len(rw_slice_flattened) < 2000:
            pad_size = 2000 - len(rw_slice_flattened)
            rw_slice_flattened = np.pad(rw_slice_flattened, (0, pad_size), mode='edge')
        elif len(rw_slice_flattened) > 2000:
            rw_slice_flattened = rw_slice_flattened[:2000]  # Trim if longer than 2000
        
        dataRwPower_trimmed.append(rw_slice_flattened)

    # Convert dataRwPower_trimmed to a NumPy array for consistency
    dataRwPower_trimmed = np.array(dataRwPower_trimmed)

    # Debugging: Check the shapes and lengths after trimming/padding
    print(f"Shape of timeData_in_minutes: {timeData_in_minutes.shape}")
    for idx in range(numRW):
        print(f"Shape of dataRwPower_trimmed[{idx}]: {dataRwPower_trimmed[idx].shape}")

    # Ensure lengths match between time and power data
    assert len(timeData_in_minutes) == len(dataRwPower_trimmed[0]), "Length mismatch between time data and RW power data"

    # Plotting
    plt.figure(3)
    for idx in range(numRW):
        plt.plot(timeData_in_minutes, dataRwPower_trimmed[idx], label=f'$p_{{rw,{idx}}}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Power (W)')
    plt.show()

# Example usage
# plot_rw_power0(b_model.timeData, dataRwPower, 4)


# Example usage
# plot_rw_power0(b_model.timeData, dataRwPower, 4)

# Example usage
# plot_rw_power0(b_model.timeData, dataRwPower, 4)

# Example usage
# plot_rw_power0(b_model.timeData, dataRwPower, 4)

def plot_rw_power(timeData, dataRwPower, numRW):
    """Plot the RW actual motor torques."""
    #timeData = timeData[:-1]  # Remove the last item from timeData if necessary
    print(timeData)
    timeData = (timeData * 1e-7)
    for idx in range(4):
        # Flatten the data if needed
        #flat_data = np.ravel(dataRwPower[idx])  # Convert to 1D array if it's not already
        
        plt.plot(timeData, dataRwPower[idx],
               
                 label='$p_{rw,' + str(idx) + '}$')
    
    plt.legend(loc='lower right')
    plt.xlabel('Time')
    plt.ylabel('RW Power (W)')
    plt.show()

def plot_rw_power2(timeData, dataRwPower, numRW):
    """Plot the RW actual motor torques over the time span."""
    
    # Check if timeData spans from 0 to 8000 seconds
    if timeData is not None and len(timeData) > 0:
        total_time = (timeData[-1] - timeData[0]) * 1e-9  # Convert nanoseconds to seconds
        if total_time != 8000:
            print(f"Warning: total simulation time is {total_time} seconds, expected 8000 seconds.")
    
    # Now plot the RW power data for each RW
    for idx in range(numRW):
        # Concatenate the sub-arrays for RW power data into a single flat array
        flat_data = np.concatenate(dataRwPower[idx])
        
        # Ensure timeData matches the length of the flattened data
        if len(timeData) < len(flat_data):
            print(f"Warning: timeData length {len(timeData)} is shorter than flattened data length {len(flat_data)}.")
            # Trimming flat_data to match timeData's length if necessary
            flat_data = flat_data[:len(timeData)]
        elif len(timeData) > len(flat_data):
            print(f"Warning: timeData length {len(timeData)} is longer than flattened data length {len(flat_data)}.")
            # Trimming timeData to match flat_data's length if necessary
            timeData = timeData[:len(flat_data)]
        
        # Convert timeData to minutes for better readability (assuming it's in nanoseconds)
        timeData_in_minutes = (timeData * 1e-7) #/ 60.0  # Convert to seconds and then to minutes
        print(len(timeData_in_minutes))
        # Plot the RW power over time
        plt.plot(timeData_in_minutes, flat_data,
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label='$p_{rw,' + str(idx) + '}$')
    print(f"Start time: {timeData[0] * 1e-7} seconds")
    print(f"End time: {timeData[-1] * 1e-7} seconds")

    plt.legend(loc='lower right')
   # plt.xlabel('Time [min]')
    plt.ylabel('RW Power (W)')
    plt.show()



def plot_rw_speeds(timeData, dataOmegaRW, numRW, figID=None):
    """Plot the RW spin rates."""
    plt.figure(figID, figsize=(5, 2.75))
    for idx in range(numRW):
        plt.plot(timeData, dataOmegaRW[:, idx] / macros.RPM,
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label=r'$\Omega_{' + str(idx + 1) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [hours]')
    plt.ylabel('RW Speed (RPM) ')
    plt.show()
from scipy.ndimage import uniform_filter1d

def plot_rw_power3(timeData, dataRwPower, numRW):
    """Plot the RW actual motor torques with improved visualization."""
    
    timeData = timeData[:-1]  # Removing the last element to match lengths
    timeData_scaled = timeData / 1e10  # Rescale the time axis for better readability
    
    fig, axs = plt.subplots(numRW, 1, figsize=(12, 10), sharex=True)

    for idx in range(numRW):
        # Concatenate the sub-arrays to get a single flat array
        flat_data = np.concatenate(dataRwPower[idx])
        
        # Apply smoothing to the flat_data
        smoothed_data = uniform_filter1d(flat_data, size=50)
        
        # Ensure timeData matches the length of the flattened data
        if len(timeData_scaled) < len(smoothed_data):
            print(f"Warning: timeData length {len(timeData_scaled)} is shorter than flattened data length {len(smoothed_data)}.")
            # Trimming smoothed_data to match timeData's length if necessary
            smoothed_data = smoothed_data[:len(timeData_scaled)]
        elif len(timeData_scaled) > len(smoothed_data):
            print(f"Warning: timeData length {len(timeData_scaled)} is longer than flattened data length {len(smoothed_data)}.")
            # Trimming timeData to match smoothed_data's length if necessary
            timeData_scaled = timeData_scaled[:len(smoothed_data)]
        
        # Now, plot the trimmed and smoothed data
        axs[idx].plot(timeData_scaled, smoothed_data,
                      color=unitTestSupport.getLineColor(idx, numRW),
                      label=f'$p_{{rw,{idx}}}$')
        
        axs[idx].legend(loc='lower right')
        axs[idx].set_ylabel('RW Power (W)')
    
    axs[-1].set_xlabel('Time [min]')
    plt.tight_layout()
    plt.show()
def make_env(env_name: str, env_num: int, faulty: bool, torque_mode: str, vec_env_cls=DummyVecEnv):
    env_classes = {
        "benv01": BEnv1,
        "benv02": BEnv2,
        "benv03": BEnv3,
        "benv04": BEnv4,
        "benv05": BEnv5,
        "benv06": BEnv6,
    
    }
    if env_name not in env_classes.keys():
        raise ValueError("env name not found")

    if env_num == 1:
        env = env_classes[env_name](faulty=faulty, torque_mode=torque_mode)
        env.reset()
    else:
        env = make_vec_env(
            env_classes[env_name], n_envs=env_num, env_kwargs={"faulty": faulty, "torque_mode": torque_mode}, vec_env_cls=vec_env_cls,seed=42
        )

    return env
