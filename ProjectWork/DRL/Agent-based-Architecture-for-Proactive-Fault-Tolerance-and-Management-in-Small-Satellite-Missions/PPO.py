import os

import pandas as pd

os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
import numpy as np
from stable_baselines3 import PPO
import matplotlib.pyplot as plt
import Tools





def evaluate(path: str, episode_num=10):
    model = PPO.load(path)
    reward_history = []
    ss_error_history = []
    for k in range(episode_num):
        episode_reward = 0
        env = Tools.make_env("benv01", 1, False, "wheel")
        obs, _ = env.reset()
        for i in range(6000):
            action, _states = model.predict(obs, deterministic=True)
            obs, reward, terminated, truncated, info = env.step(action)
            episode_reward += reward
        ss_error = env.model.error_angle_history[-1]
        reward_history.append(episode_reward)
        ss_error_history.append(ss_error)
        print(f"episode: {k}/{episode_num}, reward: {episode_reward}, error: {ss_error}")
    plt.plot(reward_history)
    plt.show()
    plt.plot(ss_error_history)
    plt.show()
import os

os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
from stable_baselines3 import PPO
from Envs import BasiliskEnv
import Tools
from Basilisk.utilities import unitTestSupport
def plot_rw_power(timeData, dataRwPower, numRW):
    """Plot the RW actual motor torques."""
    plt.figure(3)
    for idx in range(3):
        plt.plot(timeData, dataRwPower[idx],
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label='$p_{rw,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Power (W)')
def load(path: str):
    model = PPO.load(path)
    env = BasiliskEnv(faulty=True, torque_mode="wheel")
    obs, _ = env.reset()

    for i in range(8000):
        action, _states = model.predict(obs)
        
        obs, reward, terminated, truncated, info = env.step(action)
        
        #print(f"action: {action}")
        
    dataRW = []
    dataRwPower = []
    for c in range(0, 4):
        dataRW.append(env.model.rwOutLog[c].u_current)
        dataRwPower.append(env.model.rwPowLog[c].netPower)

    




    timeData=env.model.timeData[:-1]
    torque_history=np.array( env.torque_history)
    b_MRP_history=np.array(env.model.MRP_history)
    b_error_MRP_history=np.array(env.model.error_MRP_history)
    b_omega_history=np.array(env.model.omega_history)
    dataRwPower = np.array(dataRwPower)
    
    
    
    
    dataRwPower = dataRwPower[:, :-1]
    torque_history = torque_history[:-1]
    b_MRP_history=b_MRP_history[: -1]
    b_error_MRP_history=b_error_MRP_history[: -1]
    b_omega_history=b_omega_history[:-1]
    print('timeData', len(timeData))
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
    # Create a dictionary for storing all the data to be written to Excel
    data_dict = {
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

    }


        # Set the desired display precision for floats
    pd.options.display.float_format = '{:.10f}'.format  # Adjust the precision as needed

    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(data_dict)

    # Format specific columns (if needed)
    for col in df.columns:
        df[col] = df[col].map(lambda x: '{:.10f}'.format(x))

    # Save to Excel
    filename = "simulation_results_PPO3.xlsx"
    df.to_excel(filename, index=False)













    #Tools.plot_rw_power0(env.model.timeData,dataRwPower, 4)
    #Tools.plot_rw_speeds(env.model.timeData, env.model.dataOmegaRW, 4, figID=3)

    #Tools.plot_rw_power2(env.model.timeData, env.model.dataRwPower, 4)
    
    
    Tools.plot_individual(env.model,env.ref_MRP,env.torque_history,env)
    fig, axes =Tools.plot_history(env.model,env.ref_MRP,env.torque_history,env)
    fig.suptitle(f"fault time: {env.fault_time}, wheel num: {env.wheel_num}")
    plt.show()



load("C:\\training\\train_Basilisk\\benv05_True_wheel\\model\\65536000.zip")
