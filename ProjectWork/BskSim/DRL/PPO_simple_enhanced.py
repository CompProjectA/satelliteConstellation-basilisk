#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPO_simple_enhanced.py

ENHANCED PPO Training for Task Reassignment
NOW PULLS CONFIGURATION FROM GUI!

Features:
- Uses actual spacecraft names from constellation
- Pulls cluster configuration
- Uses real fault scenarios
- Environment variables from GUI

FIXED: All Unicode symbols replaced with ASCII for Windows compatibility
"""

import os
import sys
import numpy as np
import json
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Stable-Baselines3
try:
    from stable_baselines3 import PPO
    from stable_baselines3.common.env_util import make_vec_env
    from stable_baselines3.common.callbacks import BaseCallback
    import gymnasium as gym
    SB3_AVAILABLE = True
except ImportError:
    SB3_AVAILABLE = False
    print("Install stable-baselines3: pip install stable-baselines3 gymnasium")
    sys.exit(1)

# Results directory
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Strategy names
STRATEGY_NAMES = ["even_distribution", "capability_based", "load_balanced"]


class TaskReassignmentEnv(gym.Env):
    """
    ENHANCED Task reassignment environment that uses real GUI configuration
    """
    
    def __init__(self, spacecraft_names=None, cluster_config=None, fault_config=None):
        super().__init__()
        
        # USE ACTUAL SPACECRAFT NAMES FROM GUI
        self.spacecraft_names = spacecraft_names or [f"Sat{i+1}" for i in range(8)]
        self.n_satellites = len(self.spacecraft_names)
        
        # Cluster configuration (if provided)
        self.cluster_config = cluster_config or {}
        self.fault_config = fault_config or {}
        
        print(f"Environment initialized with {self.n_satellites} spacecraft:")
        for i, name in enumerate(self.spacecraft_names):
            print(f"  [{i}] {name}")
        
        # Observation: health (0/1) + load (0-1) for each satellite
        self.observation_space = gym.spaces.Box(
            low=0, high=1, 
            shape=(self.n_satellites * 2,), 
            dtype=np.float32
        )
        
        # Action: 3 strategies
        self.action_space = gym.spaces.Discrete(3)
        
        # Tracking
        self.episode_count = 0
        self.strategy_counts = np.zeros(3)
        self.episode_rewards = []
        
        self.reset()
    
    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        
        # Initialize spacecraft states
        self.sat_health = np.ones(self.n_satellites, dtype=np.float32)
        self.sat_load = np.random.uniform(0.2, 0.6, self.n_satellites).astype(np.float32)
        
        # Store initial state
        self.initial_health = self.sat_health.copy()
        self.initial_load = self.sat_load.copy()
        self.initial_names = self.spacecraft_names.copy()
        
        # Inject faults based on configuration or random
        if self.fault_config:
            # Use actual fault configuration
            for sat_name, fault_info in self.fault_config.items():
                if sat_name in self.spacecraft_names:
                    idx = self.spacecraft_names.index(sat_name)
                    if fault_info.get('enabled'):
                        self.sat_health[idx] = 0.0
                        print(f"  Injected fault: {sat_name} (index {idx})")
        else:
            # Random faults (fallback)
            n_faults = np.random.randint(1, 3)
            faulty_indices = np.random.choice(self.n_satellites, n_faults, replace=False)
            self.sat_health[faulty_indices] = 0.0
        
        # Orphaned tasks from faulty satellites
        self.orphaned_tasks = sum(self.sat_load[self.sat_health == 0])
        
        self.step_count = 0
        self.episode_count += 1
        
        return self._get_obs(), {}
    
    def _get_obs(self):
        """Observation: [health1, load1, health2, load2, ...]"""
        obs = np.zeros(self.n_satellites * 2, dtype=np.float32)
        obs[0::2] = self.sat_health
        obs[1::2] = self.sat_load
        return obs
    
    def step(self, action):
        """Execute task reassignment strategy."""
        self.step_count += 1
        self.strategy_counts[action] += 1
        
        healthy_indices = np.where(self.sat_health > 0)[0]
        n_healthy = len(healthy_indices)
        
        if n_healthy == 0:
            reward = -10.0
            self.episode_rewards.append(reward)
            return self._get_obs(), reward, True, False, {
                'action': action, 
                'strategy': STRATEGY_NAMES[action],
                'spacecraft_names': self.spacecraft_names
            }
        
        # Apply strategy
        new_load = self.sat_load.copy()
        
        if action == 0:  # Even distribution
            per_sat = self.orphaned_tasks / n_healthy
            new_load[healthy_indices] += per_sat
            
        elif action == 1:  # Capability-based
            sorted_idx = healthy_indices[np.argsort(self.sat_load[healthy_indices])]
            weights = 1.0 / (self.sat_load[sorted_idx] + 0.1)
            weights /= weights.sum()
            new_load[sorted_idx] += self.orphaned_tasks * weights
            
        elif action == 2:  # Load-balanced
            current_loads = self.sat_load[healthy_indices]
            target_load = (current_loads.sum() + self.orphaned_tasks) / n_healthy
            new_load[healthy_indices] = target_load
        
        # Calculate reward
        max_load = new_load[healthy_indices].max()
        load_variance = new_load[healthy_indices].var()
        
        if max_load > 0.9:
            reward = -5.0
        elif load_variance < 0.05:
            reward = 10.0
        else:
            reward = 5.0 - load_variance * 10
        
        # Update state
        self.final_load = new_load.copy()
        self.final_health = self.sat_health.copy()
        self.final_names = self.spacecraft_names.copy()
        self.sat_load = new_load
        self.orphaned_tasks = 0
        
        self.episode_rewards.append(reward)
        
        done = True
        truncated = False
        
        return self._get_obs(), reward, done, truncated, {
            'action': action, 
            'strategy': STRATEGY_NAMES[action],
            'load_variance': load_variance,
            'max_load': max_load,
            'spacecraft_names': self.spacecraft_names
        }


class TrainingCallback(BaseCallback):
    """Callback to track training progress."""
    
    def __init__(self, verbose=0):
        super().__init__(verbose)
        self.episode_rewards = []
        self.episode_lengths = []
        self.strategy_counts = np.zeros(3)
    
    def _on_step(self) -> bool:
        if 'episode' in self.locals.get('infos', [{}])[0]:
            info = self.locals['infos'][0]['episode']
            self.episode_rewards.append(info['r'])
            self.episode_lengths.append(info['l'])
        return True


def create_training_plots(callback, test_results, timestamp, results_dir, spacecraft_names):
    """Create comprehensive training visualization plots WITH REAL NAMES."""
    
    print("\nGenerating training visualization plots...")
    
    sns.set_style("whitegrid")
    
    # 1. Training Reward Curve
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if len(callback.episode_rewards) > 0:
        episodes = range(1, len(callback.episode_rewards) + 1)
        ax.plot(episodes, callback.episode_rewards, linewidth=2, color='#2E86AB')
        
        if len(callback.episode_rewards) > 10:
            window = 10
            moving_avg = np.convolve(callback.episode_rewards, 
                                    np.ones(window)/window, mode='valid')
            ax.plot(range(window, len(callback.episode_rewards) + 1), 
                   moving_avg, linewidth=2, color='#A23B72', 
                   label='Moving Average (10 episodes)')
        
        ax.set_xlabel('Episode', fontsize=12)
        ax.set_ylabel('Reward', fontsize=12)
        ax.set_title('DRL Training Progress - Task Reassignment', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
    plot1_path = os.path.join(results_dir, f"training_reward_curve_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: training_reward_curve_{timestamp}.png")
    
    # 2. Strategy Selection Distribution
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    strategy_counts = test_results['strategy_counts']
    colors = ['#06A77D', '#F77F00', '#D62828']
    
    ax1.bar(STRATEGY_NAMES, strategy_counts, color=colors, alpha=0.8)
    ax1.set_ylabel('Selection Count', fontsize=12)
    ax1.set_title('Strategy Selection Frequency', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=15)
    
    ax2.pie(strategy_counts, labels=STRATEGY_NAMES, autopct='%1.1f%%',
           colors=colors, startangle=90)
    ax2.set_title('Strategy Distribution', fontsize=14, fontweight='bold')
    
    plot2_path = os.path.join(results_dir, f"strategy_distribution_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: strategy_distribution_{timestamp}.png")
    
    # 3. Load Distribution Before/After - WITH REAL NAMES!
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    example = test_results['example_episode']
    # USE REAL SPACECRAFT NAMES
    satellites = example.get('spacecraft_names', spacecraft_names)
    x = np.arange(len(satellites))
    
    # Before
    colors_before = ['red' if h == 0 else 'green' for h in example['initial_health']]
    ax1.bar(x, example['initial_load'], color=colors_before, alpha=0.7)
    ax1.axhline(y=0.9, color='red', linestyle='--', label='Overload Threshold')
    ax1.set_xlabel('Spacecraft', fontsize=12)
    ax1.set_ylabel('Task Load', fontsize=12)
    ax1.set_title('Before Task Reassignment', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(satellites, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # After
    colors_after = ['red' if h == 0 else 'blue' for h in example['final_health']]
    ax2.bar(x, example['final_load'], color=colors_after, alpha=0.7)
    ax2.axhline(y=0.9, color='red', linestyle='--', label='Overload Threshold')
    ax2.set_xlabel('Spacecraft', fontsize=12)
    ax2.set_ylabel('Task Load', fontsize=12)
    ax2.set_title(f'After Task Reassignment ({example["strategy"]})', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(satellites, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plot3_path = os.path.join(results_dir, f"load_distribution_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: load_distribution_{timestamp}.png")
    
    # 4. Performance Metrics
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    test_rewards = test_results['episode_rewards']
    episodes_test = range(1, len(test_rewards) + 1)
    ax1.plot(episodes_test, test_rewards, marker='o', linewidth=2, markersize=8, color='#06A77D')
    ax1.axhline(y=np.mean(test_rewards), color='red', linestyle='--', label=f'Mean: {np.mean(test_rewards):.2f}')
    ax1.set_xlabel('Test Episode', fontsize=12)
    ax1.set_ylabel('Reward', fontsize=12)
    ax1.set_title('Test Episode Rewards', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    load_vars = test_results['load_variances']
    ax2.plot(episodes_test, load_vars, marker='s', linewidth=2, markersize=8, color='#F77F00')
    ax2.axhline(y=0.05, color='green', linestyle='--', label='Target < 0.05')
    ax2.set_xlabel('Test Episode', fontsize=12)
    ax2.set_ylabel('Load Variance', fontsize=12)
    ax2.set_title('Load Balance Quality', fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    max_loads = test_results['max_loads']
    ax3.plot(episodes_test, max_loads, marker='^', linewidth=2, markersize=8, color='#D62828')
    ax3.axhline(y=0.9, color='red', linestyle='--', label='Overload Threshold')
    ax3.set_xlabel('Test Episode', fontsize=12)
    ax3.set_ylabel('Maximum Load', fontsize=12)
    ax3.set_title('Peak Spacecraft Load', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    stats_text = f"""
    TRAINING SUMMARY
    
    Configuration: {len(satellites)} spacecraft
    Spacecraft: {', '.join(satellites[:4])}{'...' if len(satellites) > 4 else ''}
    Total Episodes: {len(callback.episode_rewards)}
    Avg Reward: {np.mean(callback.episode_rewards):.2f}
    Best Reward: {np.max(callback.episode_rewards):.2f}
    
    TEST RESULTS
    
    Avg Reward: {np.mean(test_rewards):.2f}
    Avg Load Variance: {np.mean(load_vars):.4f}
    Avg Max Load: {np.mean(max_loads):.2f}
    
    STRATEGY PREFERENCE
    
    Even Dist: {strategy_counts[0]:.0f}
    Capability: {strategy_counts[1]:.0f}
    Load Balanced: {strategy_counts[2]:.0f}
    """
    
    ax4.text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
            verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax4.axis('off')
    ax4.set_title('Performance Summary', fontsize=13, fontweight='bold')
    
    plot4_path = os.path.join(results_dir, f"performance_metrics_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(plot4_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: performance_metrics_{timestamp}.png")
    
    print("All plots generated successfully!\n")
    
    return [plot1_path, plot2_path, plot3_path, plot4_path]


def load_gui_configuration():
    """Load configuration from GUI environment variables or config file"""
    
    # Try to load from environment variables (set by GUI)
    spacecraft_names = os.environ.get('DRL_SPACECRAFT_NAMES')
    if spacecraft_names:
        spacecraft_names = json.loads(spacecraft_names)
        print(f"[OK] Loaded spacecraft names from GUI: {spacecraft_names}")
    else:
        # Fallback to default
        spacecraft_names = [f"Sat{i+1}" for i in range(8)]
        print(f"[X] No GUI config found, using default: {spacecraft_names}")
    
    # Load cluster configuration
    cluster_config = os.environ.get('DRL_CLUSTER_CONFIG')
    if cluster_config:
        cluster_config = json.loads(cluster_config)
        print(f"[OK] Loaded cluster configuration: {len(cluster_config)} clusters")
    else:
        cluster_config = {}
    
    # Load fault configuration
    fault_config = os.environ.get('DRL_FAULT_CONFIG')
    if fault_config:
        fault_config = json.loads(fault_config)
        print(f"[OK] Loaded fault configuration: {len(fault_config)} faults")
    else:
        fault_config = {}
    
    return spacecraft_names, cluster_config, fault_config


def train_model(n_episodes=100, save_path=None):
    """Train PPO model using GUI configuration"""
    
    print("\n" + "="*60)
    print("TRAINING PPO FOR TASK REASSIGNMENT")
    print("="*60)
    
    # LOAD CONFIGURATION FROM GUI
    spacecraft_names, cluster_config, fault_config = load_gui_configuration()
    
    episodes = int(os.environ.get('DRL_EPISODES', n_episodes))
    
    print(f"\nConfiguration:")
    print(f"  Episodes: {episodes}")
    print(f"  Spacecraft: {len(spacecraft_names)}")
    print(f"  Names: {spacecraft_names}")
    print(f"  Clusters: {len(cluster_config)}")
    print(f"  Faults: {len(fault_config)}")
    
    # Create vectorized environment WITH GUI CONFIG
    env = make_vec_env(
        lambda: TaskReassignmentEnv(
            spacecraft_names=spacecraft_names,
            cluster_config=cluster_config,
            fault_config=fault_config
        ), 
        n_envs=4
    )
    
    # Create PPO model
    model = PPO(
        "MlpPolicy",
        env,
        learning_rate=3e-4,
        n_steps=2048,
        batch_size=64,
        n_epochs=10,
        gamma=0.99,
        gae_lambda=0.95,
        clip_range=0.2,
        verbose=1
    )
    
    # Training callback
    callback = TrainingCallback()
    
    # Train
    print("\nStarting training...")
    total_timesteps = episodes * 10
    model.learn(total_timesteps=total_timesteps, callback=callback)
    
    # Save model
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if save_path is None:
        save_path = os.path.join(RESULTS_DIR, f"ppo_task_reassignment_{timestamp}.zip")
    
    model.save(save_path)
    print(f"\n[OK] Model saved: {save_path}")
    
    # Test WITH ACTUAL NAMES
    print("\nTesting trained model...")
    test_env = TaskReassignmentEnv(
        spacecraft_names=spacecraft_names,
        cluster_config=cluster_config,
        fault_config=fault_config
    )
    
    test_results = {
        'episode_rewards': [],
        'strategy_counts': np.zeros(3),
        'load_variances': [],
        'max_loads': [],
        'example_episode': None
    }
    
    for episode in range(10):
        obs, _ = test_env.reset()
        action, _ = model.predict(obs, deterministic=True)
        obs, reward, done, truncated, info = test_env.step(action)
        
        test_results['episode_rewards'].append(reward)
        test_results['strategy_counts'][action] += 1
        test_results['load_variances'].append(info['load_variance'])
        test_results['max_loads'].append(info['max_load'])
        
        # Save first episode for visualization WITH NAMES
        if episode == 0:
            test_results['example_episode'] = {
                'initial_health': test_env.initial_health.copy(),
                'initial_load': test_env.initial_load.copy(),
                'final_health': test_env.final_health.copy(),
                'final_load': test_env.final_load.copy(),
                'spacecraft_names': test_env.spacecraft_names.copy(),
                'strategy': STRATEGY_NAMES[action],
                'reward': reward
            }
        
        print(f"  Episode {episode+1}: Strategy={STRATEGY_NAMES[action]}, Reward={reward:.2f}, Variance={info['load_variance']:.4f}")
    
    avg_reward = np.mean(test_results['episode_rewards'])
    print(f"\nAverage test reward: {avg_reward:.2f}")
    
    # Generate plots WITH REAL NAMES
    plot_paths = create_training_plots(callback, test_results, timestamp, RESULTS_DIR, spacecraft_names)
    
    # Save results WITH SPACECRAFT NAMES
    results = {
        "timestamp": datetime.now().isoformat(),
        "episodes": episodes,
        "spacecraft": spacecraft_names,
        "n_satellites": len(spacecraft_names),
        "clusters": cluster_config,
        "faults": fault_config,
        "total_timesteps": total_timesteps,
        "avg_test_reward": float(avg_reward),
        "test_rewards": [float(r) for r in test_results['episode_rewards']],
        "strategy_distribution": test_results['strategy_counts'].tolist(),
        "avg_load_variance": float(np.mean(test_results['load_variances'])),
        "avg_max_load": float(np.mean(test_results['max_loads'])),
        "model_path": save_path,
        "plot_paths": plot_paths
    }
    
    results_json = os.path.join(RESULTS_DIR, f"training_results_{timestamp}.json")
    with open(results_json, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"[OK] Results saved: {results_json}")
    print("="*60)
    print("\nTRAINING COMPLETE!")
    print(f"  Model: {os.path.basename(save_path)}")
    print(f"  Results: {os.path.basename(results_json)}")
    print(f"  Plots: {len(plot_paths)} visualization files")
    print(f"  Spacecraft: {', '.join(spacecraft_names)}")
    print("="*60)
    
    return model, results


if __name__ == "__main__":
    if not SB3_AVAILABLE:
        print("ERROR: stable-baselines3 not installed")
        print("Install with: pip install stable-baselines3 gymnasium seaborn")
        sys.exit(1)
    
    train_model()