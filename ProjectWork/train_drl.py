# train_fault_drl.py
# PPO training loop for the SingleSatFaultEnv

import os
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import DummyVecEnv
from fault_drl_env import SingleSatFaultEnv, EnvConfig

def make_env():
    cfg = EnvConfig(
        step_dt_s=1.0,
        episode_len_s=1200.0,
        single_fault_per_episode=True,      # set False to train on 4 faults in one episode
        include_no_fault_episodes=True,
        no_fault_prob=0.15,
        r_correct=1.0,
        r_false_alarm=-0.2,
        r_wrong_label=-0.5,
        r_missed_after_window=-1.0,
        detection_grace_s=5.0,
        seed=42
    )
    return SingleSatFaultEnv(cfg)

if __name__ == "__main__":
    env = DummyVecEnv([make_env])
    model = PPO(
        "MlpPolicy",
        env,
        verbose=1,
        tensorboard_log="./tb_fault/",
        n_steps=1024,
        batch_size=256,
        gae_lambda=0.95,
        gamma=0.995,
        learning_rate=3e-4,
        ent_coef=0.01,
        clip_range=0.2,
    )
    model.learn(total_timesteps=500_000)
    os.makedirs("models", exist_ok=True)
    model.save("models/ppo_fault_detector")

    print("Training complete. Model saved to models/ppo_fault_detector.zip")
