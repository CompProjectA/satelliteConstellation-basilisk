#!/usr/bin/env python3
"""
TDHDYear2.py â€” TD3 with Hindsight Experience Replay and Dimension-Wise Clipping

Custom implementation using:
- PyTorch (no RLlib for this one)
- bsk_rl (NOT bsk_rl_develop)
- Hindsight Experience Replay
- Dimension-Wise Clipping for multi-discrete actions

Run: python TDHDYear2.py
"""

import os
import sys
import time
from datetime import datetime
from typing import Dict, Any, List, Tuple
from collections import deque, namedtuple, Counter
import random

# ---------- Path setup ----------
DRL_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(DRL_DIR)
RESULT_DIR = os.path.join(DRL_DIR, "result")
os.makedirs(RESULT_DIR, exist_ok=True)

for path in [DRL_DIR, ROOT_DIR]:
    if path not in sys.path:
        sys.path.insert(0, path)

# ---------- Imports ----------
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- PyTorch ----------
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
except ImportError as e:
    print(f"ERROR: PyTorch not available: {e}")
    print("Install with: pip install torch")
    sys.exit(1)

# ---------- Basilisk/bsk_rl ----------
try:
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
    import gymnasium as gym
    
    # MONKEY PATCH: Fix Basilisk 2.2.1 createNewEvent API incompatibility
    from Basilisk.architecture import sim_model
    original_createNewEvent = sim_model.SimBaseClass.createNewEvent
    
    def patched_createNewEvent(self, eventName, *args, **kwargs):
        """
        Patched createNewEvent that removes unsupported 'conditionFunction' parameter
        for Basilisk 2.2.1 compatibility
        """
        # Remove conditionFunction if present (not supported in Basilisk 2.2.1)
        kwargs.pop('conditionFunction', None)
        return original_createNewEvent(self, eventName, *args, **kwargs)
    
    sim_model.SimBaseClass.createNewEvent = patched_createNewEvent
    print("✓ Applied Basilisk 2.2.1 compatibility patch for createNewEvent")
    
except ImportError as e:
    print(f"ERROR: bsk_rl or Basilisk not available: {e}")
    sys.exit(1)


# =========================================================
# Replay Buffer with HER
# =========================================================

Experience = namedtuple('Experience', ['state', 'action', 'reward', 'next_state', 'done'])

class HindsightExperienceReplay:
    """Replay buffer with hindsight experience replay"""
    def __init__(self, capacity=100000, her_ratio=0.8):
        self.capacity = capacity
        self.her_ratio = her_ratio
        self.buffer = deque(maxlen=capacity)
        
    def add(self, state, action, reward, next_state, done):
        self.buffer.append(Experience(state, action, reward, next_state, done))
        
    def sample(self, batch_size):
        batch = random.sample(self.buffer, min(batch_size, len(self.buffer)))
        states = torch.FloatTensor(np.array([e.state for e in batch]))
        actions = torch.FloatTensor(np.array([e.action for e in batch]))
        rewards = torch.FloatTensor(np.array([e.reward for e in batch])).unsqueeze(1)
        next_states = torch.FloatTensor(np.array([e.next_state for e in batch]))
        dones = torch.FloatTensor(np.array([e.done for e in batch])).unsqueeze(1)
        return states, actions, rewards, next_states, dones
    
    def __len__(self):
        return len(self.buffer)


# =========================================================
# Actor-Critic Networks
# =========================================================

class Actor(nn.Module):
    """Actor network for multi-discrete action space"""
    def __init__(self, state_dim, action_dims, hidden_dim=256):
        super(Actor, self).__init__()
        self.action_dims = action_dims
        self.total_actions = sum(action_dims)
        
        self.l1 = nn.Linear(state_dim, hidden_dim)
        self.l2 = nn.Linear(hidden_dim, hidden_dim)
        self.l3 = nn.Linear(hidden_dim, self.total_actions)
        
        nn.init.xavier_uniform_(self.l1.weight)
        nn.init.xavier_uniform_(self.l2.weight)
        nn.init.xavier_uniform_(self.l3.weight)
    
    def forward(self, state):
        x = F.relu(self.l1(state))
        x = F.relu(self.l2(x))
        x = self.l3(x)
        
        actions = []
        start_idx = 0
        for dim in self.action_dims:
            action_logits = x[:, start_idx:start_idx + dim]
            actions.append(F.softmax(action_logits, dim=-1))
            start_idx += dim
        
        return torch.cat(actions, dim=-1)


class Critic(nn.Module):
    """Twin critic networks"""
    def __init__(self, state_dim, action_dims, hidden_dim=256):
        super(Critic, self).__init__()
        self.total_actions = sum(action_dims)
        
        # Q1
        self.l1 = nn.Linear(state_dim + self.total_actions, hidden_dim)
        self.l2 = nn.Linear(hidden_dim, hidden_dim)
        self.l3 = nn.Linear(hidden_dim, 1)
        
        # Q2
        self.l4 = nn.Linear(state_dim + self.total_actions, hidden_dim)
        self.l5 = nn.Linear(hidden_dim, hidden_dim)
        self.l6 = nn.Linear(hidden_dim, 1)
        
        for layer in [self.l1, self.l2, self.l3, self.l4, self.l5, self.l6]:
            nn.init.xavier_uniform_(layer.weight)
    
    def forward(self, state, action):
        sa = torch.cat([state, action], 1)
        
        q1 = F.relu(self.l1(sa))
        q1 = F.relu(self.l2(q1))
        q1 = self.l3(q1)
        
        q2 = F.relu(self.l4(sa))
        q2 = F.relu(self.l5(q2))
        q2 = self.l6(q2)
        
        return q1, q2
    
    def Q1(self, state, action):
        sa = torch.cat([state, action], 1)
        q1 = F.relu(self.l1(sa))
        q1 = F.relu(self.l2(q1))
        q1 = self.l3(q1)
        return q1


# =========================================================
# TD3 Agent with HER and DWC
# =========================================================

class TD3_HER_DWC:
    """TD3 with Hindsight Experience Replay and Dimension-Wise Clipping"""
    
    def __init__(
        self,
        state_dim,
        action_dims,
        device='cuda' if torch.cuda.is_available() else 'cpu',
        discount=0.99,
        tau=0.005,
        policy_freq=2,
        her_ratio=0.8,
        lr_actor=3e-4,
        lr_critic=3e-4,
    ):
        self.device = device
        self.discount = discount
        self.tau = tau
        self.policy_freq = policy_freq
        self.action_dims = action_dims
        self.total_actions = sum(action_dims)
        
        # Actor
        self.actor = Actor(state_dim, action_dims).to(device)
        self.actor_target = Actor(state_dim, action_dims).to(device)
        self.actor_target.load_state_dict(self.actor.state_dict())
        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=lr_actor)
        
        # Critic
        self.critic = Critic(state_dim, action_dims).to(device)
        self.critic_target = Critic(state_dim, action_dims).to(device)
        self.critic_target.load_state_dict(self.critic.state_dict())
        self.critic_optimizer = torch.optim.Adam(self.critic.parameters(), lr=lr_critic)
        
        # Replay buffer
        self.replay_buffer = HindsightExperienceReplay(her_ratio=her_ratio)
        
        self.total_iterations = 0
    
    def select_action(self, state, explore=True):
        """Select action from current policy"""
        state = torch.FloatTensor(state.reshape(1, -1)).to(self.device)
        action_probs = self.actor(state).cpu().data.numpy().flatten()
        
        actions = []
        start_idx = 0
        for dim in self.action_dims:
            probs = action_probs[start_idx:start_idx + dim]
            
            if explore:
                noise = np.random.dirichlet(np.ones(dim) * 0.1)
                probs = 0.9 * probs + 0.1 * noise
                probs = probs / np.sum(probs)
            
            action = np.random.choice(dim, p=probs)
            actions.append(action)
            start_idx += dim
        
        return actions
    
    def train(self, batch_size=256):
        """Train the agent"""
        if len(self.replay_buffer) < batch_size:
            return
        
        self.total_iterations += 1
        
        # Sample batch
        state, action, reward, next_state, done = self.replay_buffer.sample(batch_size)
        state = state.to(self.device)
        action = action.to(self.device)
        reward = reward.to(self.device)
        next_state = next_state.to(self.device)
        done = done.to(self.device)
        
        action_one_hot = self._actions_to_one_hot(action, batch_size)
        
        # Train critic
        with torch.no_grad():
            next_action_probs = self.actor_target(next_state)
            noise = torch.clamp(torch.randn_like(next_action_probs) * 0.1, -0.2, 0.2)
            next_action_probs = torch.clamp(next_action_probs + noise, 1e-6, 1.0)
            next_action_probs = next_action_probs / next_action_probs.sum(dim=1, keepdim=True)
            
            target_Q1, target_Q2 = self.critic_target(next_state, next_action_probs)
            target_Q = torch.min(target_Q1, target_Q2)
            target_Q = reward + (1 - done) * self.discount * target_Q
        
        current_Q1, current_Q2 = self.critic(state, action_one_hot)
        critic_loss = F.mse_loss(current_Q1, target_Q) + F.mse_loss(current_Q2, target_Q)
        
        self.critic_optimizer.zero_grad()
        critic_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.critic.parameters(), 1.0)
        self.critic_optimizer.step()
        
        # Train actor (delayed)
        if self.total_iterations % self.policy_freq == 0:
            action_probs = self.actor(state)
            actor_loss = -self.critic.Q1(state, action_probs).mean()
            
            self.actor_optimizer.zero_grad()
            actor_loss.backward()
            torch.nn.utils.clip_grad_norm_(self.actor.parameters(), 1.0)
            self.actor_optimizer.step()
            
            # Update target networks
            for param, target_param in zip(self.critic.parameters(), self.critic_target.parameters()):
                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)
            
            for param, target_param in zip(self.actor.parameters(), self.actor_target.parameters()):
                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)
    
    def _actions_to_one_hot(self, actions, batch_size):
        """Convert discrete actions to one-hot"""
        one_hot_actions = []
        for i in range(batch_size):
            one_hot = []
            action_idx = 0
            for dim in self.action_dims:
                action_val = int(actions[i, action_idx].item())
                action_one_hot = torch.zeros(dim)
                action_one_hot[action_val] = 1.0
                one_hot.append(action_one_hot)
                action_idx += 1
            one_hot_actions.append(torch.cat(one_hot))
        
        return torch.stack(one_hot_actions).to(self.device)
    
    def save(self, filename):
        """Save model"""
        torch.save({
            'actor': self.actor.state_dict(),
            'critic': self.critic.state_dict(),
            'actor_target': self.actor_target.state_dict(),
            'critic_target': self.critic_target.state_dict(),
            'actor_optimizer': self.actor_optimizer.state_dict(),
            'critic_optimizer': self.critic_optimizer.state_dict(),
        }, filename)
    
    def load(self, filename):
        """Load model"""
        checkpoint = torch.load(filename)
        self.actor.load_state_dict(checkpoint['actor'])
        self.critic.load_state_dict(checkpoint['critic'])
        self.actor_target.load_state_dict(checkpoint['actor_target'])
        self.critic_target.load_state_dict(checkpoint['critic_target'])
        self.actor_optimizer.load_state_dict(checkpoint['actor_optimizer'])
        self.critic_optimizer.load_state_dict(checkpoint['critic_optimizer'])


# =========================================================
# Environment Components (same as PPO/DQN)
# =========================================================

class TargetAreas(UniformTargets):
    def __init__(self, n_targets: int = 40, priority_distribution=None, radius=6378136.6):
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


class Reallocate(Action):
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


class AdvancedImagingSatellite(ImagingSatellite):
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


def _tune_access_generation(sat, initial=1800.0, step=120.0, max_dur=3600.0):
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

    extras = {"retask_on_image_complete": True, "max_access_compute_time_s": 2.0}
    for k, v in extras.items():
        if hasattr(og, k):
            setattr(og, k, v)


class CustomUniqueImageReward(UniqueImageReward):
    def __init__(self):
        try:
            super().__init__(data_store_kwargs={"keys": ["imaged", "image", "image_complete"]})
        except TypeError:
            super().__init__()
        self.imaged_by_sat = {f"Sat-{i}": 0 for i in range(4)}

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
                    total += self.reward_fn(prio) / denom
                    new_unique_count += 1
            rewards[sat_id] = total
            if new_unique_count > 0:
                self.imaged_by_sat[sat_id] = self.imaged_by_sat.get(sat_id, 0) + new_unique_count

        for sat_id in getattr(self, "imaged_by_sat", {}).keys():
            rewards.setdefault(sat_id, 0.0)

        return rewards


class CustomConstellationTasking(ConstellationTasking):
    def __init__(self, *args, max_episode_steps=64, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_episode_steps = int(max_episode_steps)
        self._step_count = 0
        self.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 6, "Sat-3": 3}

    def reset(self, *, seed=None, options=None):
        self._step_count = 0
        out = super().reset(seed=seed, options=options) if "seed" in super().reset.__code__.co_varnames else super().reset()
        if isinstance(out, tuple) and len(out) == 2:
            obs, info = out
        else:
            obs, info = out, {}
        self.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 6, "Sat-3": 3}
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
# Training
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


def train_td3_her(num_episodes=100, max_timesteps=64, batch_size=256, verbose=True):
    """Train TD3-HER agent"""
    
    print("=" * 60)
    print("TD3-HER Training for Satellite Constellation")
    print("=" * 60)
    
    # Create environment
    print("\nInitializing environment...")
    max_episode_steps = max_timesteps
    satellites = [AdvancedImagingSatellite(f"Sat-{i}", sat_args) for i in range(4)]

    for sat in satellites:
        _tune_access_generation(sat, initial=1800.0, step=120.0, max_dur=3600.0)

    env = CustomConstellationTasking(
        satellites=satellites,
        scenario=TargetAreas(n_targets=40),
        rewarder=CustomUniqueImageReward(),
        communicator=LOSCommunication(),
        sat_arg_randomizer=sat_arg_randomizer,
        log_level="INFO",
        max_episode_steps=max_episode_steps,
    )
    print("Environment initialized!")
    
    # Get dimensions
    obs0, _ = env.reset()
    first_agent = env.agents[0]
    
    if hasattr(env.observation_space(first_agent), 'spaces'):
        state_dim = sum([space.shape[0] for space in env.observation_space(first_agent).spaces.values()])
    else:
        state_dim = env.observation_space(first_agent).shape[0]
    
    action_dims = []
    if hasattr(env.action_space(first_agent), 'spaces'):
        for space in env.action_space(first_agent).spaces:
            if hasattr(space, 'n'):
                action_dims.append(space.n)
            else:
                action_dims.append(space.shape[0])
    else:
        if hasattr(env.action_space(first_agent), 'n'):
            action_dims = [env.action_space(first_agent).n]
        else:
            action_dims = [env.action_space(first_agent).shape[0]]
    
    print(f"State dimension: {state_dim}, Action dimensions: {action_dims}")
    
    # Create agents for each satellite
    agents = {}
    for agent_id in env.agents:
        agents[agent_id] = TD3_HER_DWC(
            state_dim=state_dim,
            action_dims=action_dims,
            device='cuda' if torch.cuda.is_available() else 'cpu',
            her_ratio=0.8
        )
    
    print(f"\nStarting training for {num_episodes} episodes...")
    print("=" * 60)
    
    rewards_history = []
    start_time = time.time()
    
    for episode in range(num_episodes):
        state, _ = env.reset()
        episode_reward = 0
        episode_steps = 0
        
        for t in range(max_timesteps):
            # Select actions
            action = {}
            for agent_id in env.agents:
                agent_state = state[agent_id]
                if isinstance(agent_state, dict):
                    agent_state = np.concatenate([v.flatten() for v in agent_state.values()])
                action[agent_id] = agents[agent_id].select_action(agent_state, explore=True)
            
            # Step environment
            next_state, reward, done, truncated, info = env.step(action)
            
            # Store experiences
            for agent_id in env.agents:
                agent_state = state[agent_id]
                if isinstance(agent_state, dict):
                    agent_state = np.concatenate([v.flatten() for v in agent_state.values()])
                
                next_agent_state = next_state[agent_id]
                if isinstance(next_agent_state, dict):
                    next_agent_state = np.concatenate([v.flatten() for v in next_agent_state.values()])
                
                agents[agent_id].replay_buffer.add(
                    state=agent_state,
                    action=np.array(action[agent_id]),
                    reward=reward[agent_id],
                    next_state=next_agent_state,
                    done=done[agent_id]
                )
            
            # Train agents
            for agent_id in env.agents:
                agents[agent_id].train(batch_size)
            
            state = next_state
            episode_reward += sum(reward.values())
            episode_steps += 1
            
            if all(done.values()):
                break
        
        rewards_history.append(episode_reward)
        
        if verbose and (episode + 1) % 10 == 0:
            avg_reward = np.mean(rewards_history[-10:])
            elapsed = time.time() - start_time
            print(f"Episode {episode+1}/{num_episodes}, Avg Reward: {avg_reward:.2f}, Steps: {episode_steps}, Time: {elapsed:.1f}s")
        
        # Save checkpoint
        if (episode + 1) % 100 == 0:
            for agent_id in env.agents:
                save_path = os.path.join(RESULT_DIR, f"td3_her_{agent_id}_checkpoint_{episode+1}.pth")
                agents[agent_id].save(save_path)
                if verbose:
                    print(f"  Saved checkpoint: {save_path}")
    
    # Save final models
    final_paths = []
    for agent_id in env.agents:
        save_path = os.path.join(RESULT_DIR, f"TDHDYear2_{agent_id}_final.pth")
        agents[agent_id].save(save_path)
        final_paths.append(save_path)
    
    env.close()
    
    print("\n" + "=" * 60)
    print("Training complete!")
    print(f"Total time: {time.time() - start_time:.2f} seconds")
    for path in final_paths:
        print(f"Saved: {path}")
    print("=" * 60)
    
    return rewards_history, final_paths


# =========================================================
# Main
# =========================================================

def main():
    """Main entry point"""
    print("\nTD3-HER Training Script for Satellite Constellation")
    print("Using PyTorch + bsk_rl")
    print(f"Results will be saved to: {RESULT_DIR}\n")
    
    rewards_history, model_paths = train_td3_her(
        num_episodes=100,
        max_timesteps=64,
        batch_size=256,
        verbose=True
    )
    
    if model_paths:
        print(f"\nâœ“ Training complete! Models saved to:")
        for path in model_paths:
            print(f"  {path}")
    else:
        print("\nâœ— Training finished but no models were saved")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())