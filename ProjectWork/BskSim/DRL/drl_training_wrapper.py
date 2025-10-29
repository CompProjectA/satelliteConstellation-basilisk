#!/usr/bin/env python
"""
drl_training_wrapper.py

Wrapper for training DRL algorithms from GUI.
Extracts core training logic from notebooks.
"""

import os
import sys
import numpy as np
from datetime import datetime
from typing import Dict, Optional, Callable
import json

# Add paths
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)


class DRLTrainingWrapper:
    """
    Unified training wrapper for PPO, DQN, and TD3HD algorithms
    """
    
    def __init__(self, algorithm="PPO"):
        self.algorithm = algorithm
        self.config = self._get_default_config()
        self.training_results = None
        self.model_save_path = None
        self.progress_callback = None
    
    def _get_default_config(self) -> Dict:
        """Get default configuration for algorithm"""
        
        base_config = {
            "n_satellites": 4,
            "n_targets": 40,
            "n_episodes": 100,
            "max_steps_per_episode": 1000,
            "orbit_altitude": 600,  # km
            "orbit_inclination": 97.0,  # degrees
            "target_lat_range": (-60, 60),
            "target_lon_range": (0, 360),
        }
        
        if self.algorithm == "PPO":
            base_config.update({
                "learning_rate": 5e-5,
                "train_batch_size": 4000,
                "sgd_minibatch_size": 128,
                "num_sgd_iter": 30,
                "gamma": 0.995,
                "lambda_": 0.95,
                "num_workers": 4
            })
        
        elif self.algorithm == "DQN":
            base_config.update({
                "learning_rate": 1e-4,
                "buffer_size": 10000,
                "batch_size": 32,
                "gamma": 0.99,
                "epsilon_start": 1.0,
                "epsilon_end": 0.01,
                "epsilon_decay": 0.995,
                "target_update_frequency": 100
            })
        
        elif self.algorithm == "TD3HD":
            base_config.update({
                "learning_rate": 3e-4,
                "buffer_size": 1000000,
                "batch_size": 256,
                "gamma": 0.99,
                "tau": 0.005,
                "policy_noise": 0.2,
                "noise_clip": 0.5,
                "policy_delay": 2
            })
        
        return base_config
    
    def set_config(self, **kwargs):
        """Update configuration parameters"""
        self.config.update(kwargs)
    
    def set_progress_callback(self, callback: Callable):
        """Set callback function for progress updates"""
        self.progress_callback = callback
    
    def _report_progress(self, episode: int, reward: float, message: str = ""):
        """Report training progress"""
        if self.progress_callback:
            self.progress_callback(episode, reward, message)
        else:
            print(f"Episode {episode}: Reward = {reward:.2f} | {message}")
    
    def train(self, save_dir="DRL/result"):
        """
        Train the selected algorithm
        """
        
        print(f"\n{'='*60}")
        print(f"Training {self.algorithm} Algorithm")
        print(f"{'='*60}")
        print(f"Configuration:")
        for key, value in self.config.items():
            print(f"  {key}: {value}")
        print(f"{'='*60}\n")
        
        os.makedirs(save_dir, exist_ok=True)
        
        try:
            if self.algorithm == "PPO":
                results = self._train_ppo(save_dir)
            elif self.algorithm == "DQN":
                results = self._train_dqn(save_dir)
            elif self.algorithm == "TD3HD":
                results = self._train_td3hd(save_dir)
            else:
                raise ValueError(f"Unknown algorithm: {self.algorithm}")
            
            self.training_results = results
            return results
            
        except Exception as e:
            print(f"Training failed: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def _train_ppo(self, save_dir):
        """Train PPO algorithm"""
        
        from ray.rllib.algorithms.ppo import PPOConfig
        from ray.tune.registry import register_env
        
        # Import environment (adjust path as needed)
        try:
            from PPO import TaskEnv  # Your custom environment
            register_env("TaskEnv", lambda config: TaskEnv(config))
        except ImportError:
            print("Warning: Could not import TaskEnv from PPO.py")
            print("Using placeholder environment")
        
        # Configure PPO
        config = PPOConfig()
        config = config.environment(
            env="TaskEnv",
            env_config={
                "n_satellites": self.config["n_satellites"],
                "n_targets": self.config["n_targets"]
            }
        )
        config = config.framework("torch")
        config = config.rollouts(num_rollout_workers=self.config["num_workers"])
        config = config.training(
            train_batch_size=self.config["train_batch_size"],
            sgd_minibatch_size=self.config["sgd_minibatch_size"],
            num_sgd_iter=self.config["num_sgd_iter"],
            lr=self.config["learning_rate"],
            gamma=self.config["gamma"],
            lambda_=self.config["lambda_"]
        )
        
        # Build algorithm
        algo = config.build()
        
        # Training loop
        results = {
            "episodes": [],
            "rewards": [],
            "episode_lengths": []
        }
        
        for episode in range(self.config["n_episodes"]):
            result = algo.train()
            
            reward = result["episode_reward_mean"]
            length = result["episode_len_mean"]
            
            results["episodes"].append(episode)
            results["rewards"].append(reward)
            results["episode_lengths"].append(length)
            
            self._report_progress(episode, reward, f"Length: {length:.0f}")
            
            # Save checkpoint every 10 episodes
            if (episode + 1) % 10 == 0:
                checkpoint_path = algo.save(save_dir)
                print(f"  Checkpoint saved: {checkpoint_path}")
        
        # Save final model
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        final_path = os.path.join(save_dir, f"ppo_final_{timestamp}")
        algo.save(final_path)
        
        self.model_save_path = final_path
        print(f"\n✓ Training complete! Model saved to: {final_path}")
        
        # Save training results
        results_file = os.path.join(save_dir, f"training_results_ppo_{timestamp}.json")
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def _train_dqn(self, save_dir):
        """Train DQN algorithm"""
        
        print("DQN training implementation...")
        
        # Import DQN implementation
        from DQNYear2 import DQNAgent, BasiliskEnv
        
        # Create environment
        env = BasiliskEnv(faulty=False)
        
        # Create agent
        agent = DQNAgent(
            state_dim=env.observation_space.shape[0],
            action_dim=env.action_space.n,
            learning_rate=self.config["learning_rate"],
            gamma=self.config["gamma"],
            epsilon=self.config["epsilon_start"],
            epsilon_min=self.config["epsilon_end"],
            epsilon_decay=self.config["epsilon_decay"]
        )
        
        # Training loop
        results = {
            "episodes": [],
            "rewards": [],
            "epsilon": []
        }
        
        for episode in range(self.config["n_episodes"]):
            state = env.reset()[0]
            episode_reward = 0
            done = False
            steps = 0
            
            while not done and steps < self.config["max_steps_per_episode"]:
                action = agent.select_action(state)
                next_state, reward, terminated, truncated, _ = env.step(action)
                done = terminated or truncated
                
                agent.store_transition(state, action, reward, next_state, done)
                agent.train()
                
                state = next_state
                episode_reward += reward
                steps += 1
            
            results["episodes"].append(episode)
            results["rewards"].append(episode_reward)
            results["epsilon"].append(agent.epsilon)
            
            self._report_progress(episode, episode_reward, f"Epsilon: {agent.epsilon:.3f}")
        
        # Save model
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        model_path = os.path.join(save_dir, f"dqn_model_{timestamp}.pth")
        agent.save(model_path)
        
        self.model_save_path = model_path
        print(f"\n✓ Training complete! Model saved to: {model_path}")
        
        return results
    
    def _train_td3hd(self, save_dir):
        """Train TD3HD algorithm"""
        
        print("TD3HD training implementation...")
        
        # Import TD3HD implementation
        from TDHDYear2 import TD3HDAgent, BasiliskEnv
        
        # Create environment
        env = BasiliskEnv(faulty=False)
        
        # Create agent
        agent = TD3HDAgent(
            state_dim=env.observation_space.shape[0],
            action_dim=env.action_space.shape[0],
            learning_rate=self.config["learning_rate"],
            gamma=self.config["gamma"]
        )
        
        # Training loop (similar to DQN)
        results = {
            "episodes": [],
            "rewards": []
        }
        
        for episode in range(self.config["n_episodes"]):
            state = env.reset()[0]
            episode_reward = 0
            done = False
            steps = 0
            
            while not done and steps < self.config["max_steps_per_episode"]:
                action = agent.select_action(state)
                next_state, reward, terminated, truncated, _ = env.step(action)
                done = terminated or truncated
                
                agent.store_transition(state, action, reward, next_state, done)
                agent.train()
                
                state = next_state
                episode_reward += reward
                steps += 1
            
            results["episodes"].append(episode)
            results["rewards"].append(episode_reward)
            
            self._report_progress(episode, episode_reward)
        
        # Save model
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        model_path = os.path.join(save_dir, f"td3hd_model_{timestamp}.pth")
        agent.save(model_path)
        
        self.model_save_path = model_path
        print(f"\n✓ Training complete! Model saved to: {model_path}")
        
        return results


# Convenience function for GUI
def train_from_gui(algorithm, config, progress_callback=None):
    """
    Train DRL algorithm with GUI integration
    
    Args:
        algorithm: "PPO", "DQN", or "TD3HD"
        config: Dictionary with training parameters
        progress_callback: Function(episode, reward, message) for progress updates
    
    Returns:
        (results, model_path): Training results and saved model path
    """
    
    trainer = DRLTrainingWrapper(algorithm)
    trainer.set_config(**config)
    
    if progress_callback:
        trainer.set_progress_callback(progress_callback)
    
    results = trainer.train()
    
    return results, trainer.model_save_path


if __name__ == "__main__":
    print("DRL Training Wrapper - Test Mode")
    
    # Test with small configuration
    test_config = {
        "n_satellites": 4,
        "n_targets": 10,
        "n_episodes": 5
    }
    
    trainer = DRLTrainingWrapper("PPO")
    trainer.set_config(**test_config)
    
    def test_callback(episode, reward, message):
        print(f"GUI Update: Episode {episode}, Reward {reward:.2f} - {message}")
    
    trainer.set_progress_callback(test_callback)
    
    try:
        results = trainer.train()
        print("\nTraining Results:")
        print(f"  Episodes: {len(results['episodes'])}")
        print(f"  Average Reward: {np.mean(results['rewards']):.2f}")
        print(f"  Model Path: {trainer.model_save_path}")
    except Exception as e:
        print(f"Training test failed: {e}")