import os
import argparse

from stable_baselines3 import PPO
from stable_baselines3.common.monitor import Monitor

from friction_fault_env import FrictionFaultEnv


# Default directories for saving models and logs
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_LOG_DIR = os.path.join(THIS_DIR, "ppo_logs")
DEFAULT_MODEL_PATH = os.path.join(THIS_DIR, "friction_fault_ppo_model")


def train(model_path: str, timesteps: int = 10000, log_dir: str = DEFAULT_LOG_DIR) -> None:
    """Train a PPO agent on the FrictionFaultEnv."""
    os.makedirs(log_dir, exist_ok=True)

    env = FrictionFaultEnv()
    env = Monitor(env, log_dir)

    model = PPO(
        "MlpPolicy",
        env,
        verbose=1,
        tensorboard_log=log_dir,
    )
    model.learn(total_timesteps=timesteps)
    model.save(model_path)
    env.close()


def evaluate(model_path: str, episodes: int = 5) -> None:
    """Evaluate a trained PPO agent."""
    env = FrictionFaultEnv()
    model = PPO.load(model_path)

    for ep in range(episodes):
        obs, _ = env.reset()
        done = False
        total_reward = 0.0
        while not done:
            action, _ = model.predict(obs, deterministic=True)
            obs, reward, terminated, truncated, _ = env.step(action)
            done = terminated or truncated
            total_reward = reward
        print(f"Episode {ep , 1}: reward={total_reward}")
    env.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train or evaluate PPO on FrictionFaultEnv")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--train", action="store_true", help="Run training mode")
    group.add_argument("--eval", action="store_true", help="Run evaluation mode")
    parser.add_argument("--model-path", default=DEFAULT_MODEL_PATH, help="Path to save or load the PPO model")
    parser.add_argument("--timesteps", type=int, default=10000, help="Timesteps for training")
    parser.add_argument("--episodes", type=int, default=5, help="Episodes for evaluation")
    parser.add_argument("--log-dir", default=DEFAULT_LOG_DIR, help="Directory for tensorboard logs")
    args = parser.parse_args()

    if args.train:
        train(args.model_path, args.timesteps, args.log_dir)
    elif args.eval:
        evaluate(args.model_path, args.episodes)
 