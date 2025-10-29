#!/usr/bin/env python3
"""
DQNYear2.py - DQN with Ray RLlib for bsk_rl multi-agent satellite constellation

Based on working DQN.ipynb code.
Uses Ray RLlib + bsk_rl for multi-agent constellation tasking.
Saves checkpoints to DRL/result/ directory.
"""

import time
import os
import sys
from datetime import datetime
from typing import Dict, Any

# Path setup
DRL_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(DRL_DIR)
RESULT_DIR = os.path.join(DRL_DIR, "result")
os.makedirs(RESULT_DIR, exist_ok=True)

for path in [DRL_DIR, ROOT_DIR]:
    if path not in sys.path:
        sys.path.insert(0, path)

# Third-party imports
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gymnasium as gym

try:
    import pandas as pd
except Exception:
    pd = None

from collections import Counter

# Ray RLlib imports
from ray.tune.registry import register_env
from ray.rllib.policy.policy import PolicySpec
from ray.rllib.env.multi_agent_env import MultiAgentEnv
from ray.rllib.algorithms.dqn import DQNConfig, DQN

# Basilisk imports
from Basilisk.architecture import messaging

# bsk_rl imports
from bsk_rl import ConstellationTasking
from bsk_rl.sats import ImagingSatellite
from bsk_rl.act import Action, Image
from bsk_rl.act.actions import ActionBuilder
from bsk_rl import obs
from bsk_rl.sim import dyn, fsw
from bsk_rl.scene.targets import UniformTargets, Target
from bsk_rl.data import UniqueImageReward
from bsk_rl.comm import LOSCommunication
from bsk_rl.utils.orbital import walker_delta_args


# =========================================================
# Target Area Definition
# =========================================================

class TargetAreas(UniformTargets):
    """Custom target area generator for specific geographic region"""
    def __init__(self, n_targets: int = 20, priority_distribution=None, radius=6378136.6):
        super().__init__(n_targets, priority_distribution, radius)
        self.lat_min, self.lat_max = -38.0, -25.0
        self.lon_min, self.lon_max = 129.0, 141.0

    def regenerate_targets(self) -> None:
        n_targets = self.n_targets if isinstance(self.n_targets, int) else np.random.randint(
            self.n_targets[0], self.n_targets[1] + 1
        )
        lats = np.random.uniform(self.lat_min, self.lat_max, n_targets)
        lons = np.random.uniform(self.lon_min, self.lon_max, n_targets)
        r = self.radius
        targets = []
        for i, (lat, lon) in enumerate(zip(lats, lons)):
            lat_rad = np.radians(lat)
            lon_rad = np.radians(lon)
            x = r * np.cos(lat_rad) * np.cos(lon_rad)
            y = r * np.cos(lat_rad) * np.sin(lon_rad)
            z = r * np.sin(lat_rad)
            priority = self.priority_distribution() if self.priority_distribution else np.random.uniform(0, 1)
            targets.append(Target(f"SA_Target_{i}", [x, y, z], priority))
        self.targets = targets


# =========================================================
# Custom Action: Reallocate
# =========================================================

class Reallocate(Action):
    """Task reallocation action for multi-agent coordination"""
    def __init__(self, n_sats=4):
        self.action_space = gym.spaces.Discrete(n_sats)
        self.n_actions = self.action_space.n
        self.option = None

    @property
    def builder_type(self):
        return Image.builder_type if isinstance(Image.builder_type, type) else ActionBuilder

    def set_action(self, option: int, **kwargs):
        self.option = option

    def action(self, satellite, state):
        return 0.0


# =========================================================
# Advanced Imaging Satellite
# =========================================================

class AdvancedImagingSatellite(ImagingSatellite):
    """Satellite with imaging and reallocation capabilities"""
    observation_spec = [
        obs.OpportunityProperties(
            dict(prop="priority"),
            dict(prop="opportunity_open", norm=5700.0),
            n_ahead_observe=5,
        )
    ]
    action_spec = [Image(n_ahead_image=5), Reallocate()]
    dyn_type = dyn.FullFeaturedDynModel
    fsw_type = fsw.SteeringImagerFSWModel


# =========================================================
# Opportunity Generator Tuning
# =========================================================

def _tune_access_generation(sat, initial=600.0, step=300.0, max_dur=1800.0):
    """Configure opportunity generation parameters for faster training"""
    candidates = [
        getattr(getattr(sat, "fsw", None), "opportunity_generator", None),
        getattr(getattr(sat, "dynamics", None), "accessGenerator", None),
        getattr(getattr(sat, "dynamics", None), "opportunityGenerator", None),
    ]
    og = next((c for c in candidates if c is not None), None)
    if og is None:
        return

    for name, val in [
        ("initial_generation_duration", float(initial)),
        ("generation_step", float(step)),
        ("max_generation_duration", float(max_dur)),
        ("max_lookahead", float(max_dur)),
    ]:
        if hasattr(og, name):
            setattr(og, name, val)

    # Optimized settings for training speed
    extras = {
        "retask_on_image_complete": False,
        "max_access_compute_time_s": 1.0,
    }
    for k, v in extras.items():
        if hasattr(og, k):
            setattr(og, k, v)


# =========================================================
# Custom Reward Function
# =========================================================

class CustomUniqueImageReward(UniqueImageReward):
    """Enhanced reward function tracking images per satellite"""
    
    DEBUG_PROBE = False
    DEBUG_PROBE_STEPS = 5

    def __init__(self):
        try:
            super().__init__(data_store_kwargs={"keys": ["imaged", "image", "image_complete"]})
        except TypeError:
            super().__init__()
        self.imaged_by_sat = {f"Sat-{i}": 0 for i in range(4)}
        self._probe_count = 0
        self.imaged_targets = set()

    def _event_key(self, evt):
        for k in ("key", "event", "type", "name"):
            if hasattr(evt, k):
                v = getattr(evt, k)
                if isinstance(v, str):
                    return v
        if isinstance(evt, dict):
            for k in ("key", "event", "type", "name"):
                if k in evt and isinstance(evt[k], str):
                    return evt[k]
        return None

    def _event_sat(self, evt):
        for k in ("agent", "satellite", "sat", "who"):
            if hasattr(evt, k):
                return getattr(evt, k)
        if isinstance(evt, dict):
            for k in ("agent", "satellite", "sat", "who"):
                if k in evt:
                    return evt[k]
        return None

    def _event_target(self, evt):
        for k in ("target", "tgt", "obj"):
            if hasattr(evt, k):
                return getattr(evt, k)
        if isinstance(evt, dict):
            for k in ("target", "tgt", "obj"):
                if k in evt:
                    return evt[k]
        return None

    def _target_id_and_priority(self, tgt):
        if tgt is None:
            return ("<none>", 0.0)
        tid = getattr(tgt, "name", None) or getattr(tgt, "id", None)
        pr = getattr(tgt, "priority", None)
        if tid is None and isinstance(tgt, dict):
            tid = tgt.get("name") or tgt.get("id") or tgt.get("uid")
            pr = tgt.get("priority", pr)
        if tid is None:
            tid = str(tgt)
        try:
            pr = float(pr) if pr is not None else 0.0
        except Exception:
            pr = 0.0
        return (tid, pr)

    def reward(self, new_data_dict):
        all_step_targets = []
        for data in new_data_dict.values():
            imgs = getattr(data, "imaged", []) or []
            all_step_targets.extend(imgs)
        occ = Counter(all_step_targets)

        rewards = {}
        for sat_id, data in new_data_dict.items():
            total = 0.0
            new_unique_count = 0
            for tgt in getattr(data, "imaged", []) or []:
                if tgt not in self.data.imaged:
                    prio = float(getattr(tgt, "priority", 0.0))
                    denom = occ.get(tgt, 1) or 1
                    total += (1.0 + prio) / denom
                    new_unique_count += 1
                    self.data.imaged.add(tgt)
                    self.imaged_targets.add(self._target_id_and_priority(tgt)[0])

            if new_unique_count > 0:
                self.imaged_by_sat[sat_id] = self.imaged_by_sat.get(sat_id, 0) + new_unique_count

            rewards[sat_id] = total
        return rewards


# =========================================================
# Custom Constellation Environment
# =========================================================

class CustomConstellationTasking(ConstellationTasking):
    """Extended constellation environment with task tracking"""
    def __init__(self, *args, max_episode_steps=32, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_episode_steps = int(max_episode_steps)
        self._step_count = 0
        self.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 3, "Sat-3": 3}

    def reset(self, *, seed=None, options=None):
        self._step_count = 0
        out = super().reset(seed=seed, options=options) if "seed" in super().reset.__code__.co_varnames else super().reset()
        if isinstance(out, tuple) and len(out) == 2:
            obs, info = out
        else:
            obs, info = out, {}
        self.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 3, "Sat-3": 3}
        return obs, info

    def step(self, action_dict: Dict[str, Any]):
        self._step_count += 1
        obs, rews, terms, truncs, infos = super().step(action_dict)

        for agent, r in rews.items():
            if r > 0:
                self.remaining_tasks[agent] = max(0, self.remaining_tasks.get(agent, 0) - 1)

        horizon_reached = self._step_count >= self.max_episode_steps
        no_agents_left = len(self.agents) == 0

        all_keys = set(rews.keys()) | set(terms.keys()) | set(truncs.keys()) | set(self.agents)
        all_done_per_agent = all(terms.get(a, False) or truncs.get(a, False) for a in all_keys) if all_keys else True

        terms["__all__"] = no_agents_left or all_done_per_agent
        truncs["__all__"] = (horizon_reached and not terms["__all__"])

        return obs, rews, terms, truncs, infos


# =========================================================
# Environment Creator
# =========================================================

sat_args = {
    "imageAttErrorRequirement": 0.01,
    "imageRateErrorRequirement": 0.01,
    "batteryStorageCapacity": 1e9,
    "storedCharge_Init": 1e9,
    "dataStorageCapacity": 1e12,
    "u_max": 0.4,
    "K1": 0.25,
    "K3": 3.0,
    "omega_max": 0.087,
    "servo_Ki": 5.0,
    "servo_P": 150 / 5,
}
sat_arg_randomizer = walker_delta_args(altitude=800.0, inc=60.0, n_planes=1)


def env_creator(env_config):
    """Create constellation environment"""
    max_episode_steps = env_config.get("max_episode_steps", 32)
    satellites = [AdvancedImagingSatellite(f"Sat-{i}", sat_args) for i in range(4)]

    for sat in satellites:
        _tune_access_generation(sat, initial=600.0, step=60.0, max_dur=900.0)

    return CustomConstellationTasking(
        satellites=satellites,
        scenario=TargetAreas(n_targets=20),
        rewarder=CustomUniqueImageReward(),
        communicator=LOSCommunication(),
        sat_arg_randomizer=sat_arg_randomizer,
        log_level="WARNING",
        max_episode_steps=max_episode_steps,
    )


# =========================================================
# Training Function
# =========================================================

def train_dqn(num_iterations=50, save_freq=10, verbose=True):
    """Train DQN agent using Ray RLlib"""
    
    print("=" * 60)
    print("DQN Training for Satellite Constellation")
    print("=" * 60)
    
    # Register environment
    register_env("constellation", env_creator)
    
    # Create sample environment to get specs
    sample_env = env_creator({"max_episode_steps": 32})
    obs_space = sample_env.observation_space(sample_env.agents[0])
    act_space = sample_env.action_space(sample_env.agents[0])
    sample_env.close()
    
    # Configure DQN
    config = (
        DQNConfig()
        .environment("constellation", env_config={"max_episode_steps": 32})
        .framework("torch")
        .training(
            gamma=0.99,
            lr=5e-4,
            train_batch_size=64,
            target_network_update_freq=500,
            replay_buffer_config={
                "type": "MultiAgentReplayBuffer",
                "capacity": 100000,
            },
        )
        .multi_agent(
            policies={
                "shared_policy": PolicySpec(
                    policy_class=None,
                    observation_space=obs_space,
                    action_space=act_space,
                )
            },
            policy_mapping_fn=lambda agent_id, *args, **kwargs: "shared_policy",
        )
        .resources(num_gpus=0)
        .rollouts(num_rollout_workers=0)
    )
    
    # Build algorithm
    algo = DQN(config=config)
    
    print(f"\nStarting training for {num_iterations} iterations...")
    print("=" * 60)
    
    rewards_history = []
    start_time = time.time()
    checkpoint_dir = None
    
    try:
        for iteration in range(num_iterations):
            result = algo.train()
            episode_reward_mean = result.get("episode_reward_mean", 0)
            rewards_history.append(episode_reward_mean)
            
            if verbose and (iteration + 1) % 5 == 0:
                elapsed = time.time() - start_time
                print(f"Iter {iteration+1}/{num_iterations}, "
                      f"Mean Reward: {episode_reward_mean:.2f}, "
                      f"Time: {elapsed:.1f}s")
            
            # Save checkpoint
            if (iteration + 1) % save_freq == 0:
                checkpoint_dir = algo.save(RESULT_DIR)
                if verbose:
                    print(f"  Saved checkpoint: {checkpoint_dir}")
        
        # Save final checkpoint
        checkpoint_dir = algo.save(RESULT_DIR)
        print(f"Final checkpoint saved at {checkpoint_dir}")
        
    except Exception as e:
        print(f"\nTraining failed with error: {e}")
        try:
            checkpoint_dir = algo.save(RESULT_DIR)
            print(f"Emergency checkpoint saved at {checkpoint_dir}")
        except:
            print("Could not save emergency checkpoint")
    
    finally:
        print("\nTraining session completed!")
        print(f"Total training time: {time.time() - start_time:.2f} seconds")
    
    # Save training history plot
    if rewards_history:
        plt.figure(figsize=(10, 6))
        plt.plot(rewards_history)
        plt.xlabel('Iteration')
        plt.ylabel('Mean Episode Reward')
        plt.title('DQN Training Progress')
        plt.grid(True)
        plt.savefig(os.path.join(RESULT_DIR, 'dqn_training_rewards.png'))
        plt.close()
        print(f"Saved training plot: {os.path.join(RESULT_DIR, 'dqn_training_rewards.png')}")
    
    return checkpoint_dir, rewards_history


# =========================================================
# Main
# =========================================================

def main():
    """Main entry point"""
    print("\nDQN Training Script for Satellite Constellation")
    print("Using Ray RLlib + bsk_rl")
    print(f"Results will be saved to: {RESULT_DIR}\n")
    
    checkpoint_dir, rewards_history = train_dqn(
        num_iterations=50,
        save_freq=10,
        verbose=True
    )
    
    if checkpoint_dir:
        print(f"\n✓ Training complete! Checkpoint saved to:")
        print(f"  {checkpoint_dir}")
    else:
        print("\n✗ Training finished but no checkpoint was saved")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())