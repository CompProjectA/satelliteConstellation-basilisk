

import os
import time
from datetime import datetime
from typing import Dict, Any, Tuple, List

import numpy as np
import gymnasium as gym
import matplotlib.pyplot as plt

try:
    import pandas as pd
except Exception:
    pd = None

# ------- Basilisk / bsk_rl_develop imports -------
from Basilisk.architecture import messaging
from collections import Counter

from bsk_rl_develop import ConstellationTasking
from bsk_rl_develop.sats import ImagingSatellite
from bsk_rl_develop.act import Action, Image
from bsk_rl_develop.act.actions import ActionBuilder
from bsk_rl_develop import obs
from bsk_rl_develop.sim import dyn, fsw
from bsk_rl_develop.scene.targets import UniformTargets, Target
from bsk_rl_develop.data import UniqueImageReward
from bsk_rl_develop.comm import LOSCommunication
from bsk_rl_develop.utils.orbital import walker_delta_args

# ------------------- Tianshou imports -------------------
import torch
import torch.nn as nn
from tianshou.env import DummyVectorEnv
from tianshou.data import Collector, VectorReplayBuffer, Batch
from tianshou.policy import TD3Policy
from tianshou.utils.net.common import Net
from tianshou.utils.net.continuous import Actor, Critic

# Exploration noise (handle old/new paths)
try:
    from tianshou.exploration import GaussianNoise, OUNoise
except Exception:
    try:
        from tianshou.exploration.gaussian import GaussianNoise  # type: ignore
        from tianshou.exploration.ou import OUNoise              # type: ignore
    except Exception:
        GaussianNoise = None
        OUNoise = None


# =========================================================
# Target region
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
            lat_rad = np.radians(lat); lon_rad = np.radians(lon)
            x = r * np.cos(lat_rad) * np.cos(lon_rad)
            y = r * np.cos(lat_rad) * np.sin(lon_rad)
            z = r * np.sin(lat_rad)
            priority = self.priority_distribution() if self.priority_distribution else np.random.uniform(0, 1)
            targets.append(Target(f"SA_Target_{i}", [x, y, z], priority))
        self.targets = targets


# =========================================================
# Custom action: Reallocate (fix abstract builder_type)
# =========================================================
class Reallocate(Action):
    # Some builds read builder_type as class attr:
    builder_type = ActionBuilder

    def __init__(self, n_sats=4):
        self.action_space = gym.spaces.Discrete(n_sats)
        self.n_actions = self.action_space.n
        self.option = None

    # Some builds require a property or method; provide property:
    @property
    def builder_type(self):
        bt = getattr(Image, "builder_type", None)
        return bt if isinstance(bt, type) else ActionBuilder

    def set_action(self, option: int, **kwargs):
        self.option = option

    def action(self, satellite, state):
        return 0.0


# =========================================================
# Satellite and base multi-agent env
# =========================================================
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
    for k, v in {"retask_on_image_complete": True, "max_access_compute_time_s": 2.0}.items():
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
# Joint single-agent wrapper for Tianshou + DWC
# =========================================================
class JointConstellationEnv(gym.Env):
    """
    Single-agent facade over the multi-agent constellation:
    - Observation: concatenation of per-agent obs (flattened).
    - Action: Box(-1,1)^D continuous joint vector, one dim per discrete head.
      Internally mapped to per-sat discrete choices (Image + Reallocate).
    - DWC: per-dim clipping of action delta between steps.
    - Reward: sum of per-agent rewards (change to avg if desired).
    - IMPORTANT: reset/step always return the same info keys for Tianshou.
    """

    metadata = {"render_modes": []}

    def __init__(self, max_episode_steps=64, dwc_delta=0.5):
        super().__init__()
        # Build base env
        satellites = [AdvancedImagingSatellite(f"Sat-{i}", sat_args) for i in range(4)]
        for sat in satellites:
            _tune_access_generation(sat, initial=1800.0, step=120.0, max_dur=3600.0)

        self.base = CustomConstellationTasking(
            satellites=satellites,
            scenario=TargetAreas(n_targets=40),
            rewarder=CustomUniqueImageReward(),
            communicator=LOSCommunication(),
            sat_arg_randomizer=sat_arg_randomizer,
            log_level="INFO",
            max_episode_steps=max_episode_steps,
        )
        self.max_episode_steps = max_episode_steps
        self._last_action = None
        self._dwc = float(dwc_delta)

        # Reset once to discover spaces and agent IDs
        obs0, _ = self.base.reset()
        self._agent_ids = list(self.base.agents)  # e.g., ["Sat-0", ..., "Sat-3"]
        # Make a stable order for per-agent vectors
        self._sat_order = [f"Sat-{i}" for i in range(4)]

        self._per_agent_obs_space = {aid: self.base.observation_space(aid) for aid in self._agent_ids}
        self._per_agent_act_space = {aid: self.base.action_space(aid) for aid in self._agent_ids}

        # Build flattened obs
        self._obs_slices: List[Tuple[int, int]] = []
        start = 0
        for aid in self._agent_ids:
            sp = self._per_agent_obs_space[aid]
            size = int(np.prod(sp.shape if hasattr(sp, "shape") else ()))
            _ = self._flatten(obs0[aid], size)  # probe flatten
            self._obs_slices.append((start, start + size))
            start += size
        self._obs_dim = start

        # Build joint action bins per sat
        self._per_agent_bins = []
        for aid in self._agent_ids:
            act_sp = self._per_agent_act_space[aid]
            bins = []
            if isinstance(act_sp, gym.spaces.Discrete):
                bins.append(act_sp.n)
            elif isinstance(act_sp, gym.spaces.Tuple):
                for s in act_sp.spaces:
                    assert isinstance(s, gym.spaces.Discrete), "Only Discrete in Tuple supported"
                    bins.append(s.n)
            else:
                raise ValueError(f"Unsupported per-agent action space: {act_sp}")
            self._per_agent_bins.append(bins)

        self._act_bins = np.array([b for bins in self._per_agent_bins for b in bins], dtype=np.int32)
        self._act_dim = len(self._act_bins)

        self.observation_space = gym.spaces.Box(low=-np.inf, high=np.inf, shape=(self._obs_dim,), dtype=np.float32)
        self.action_space = gym.spaces.Box(low=-1.0, high=1.0, shape=(self._act_dim,), dtype=np.float32)

        self._t = 0

    def _flatten(self, ob, size):
        arr = np.array(ob, dtype=np.float32).reshape(-1)
        if arr.size != size and isinstance(ob, dict):
            # try dict-like {"observation": ...} or {"obs": ...}
            key = "observation" if "observation" in ob else ("obs" if "obs" in ob else None)
            if key is not None:
                arr = np.array(ob[key], dtype=np.float32).reshape(-1)
        return arr

    def _concat_obs(self, obs_dict):
        vecs = []
        for (start, end), aid in zip(self._obs_slices, self._agent_ids):
            size = end - start
            vecs.append(self._flatten(obs_dict[aid], size))
        return np.concatenate(vecs, dtype=np.float32)

    def _apply_dwc(self, a: np.ndarray) -> np.ndarray:
        if self._last_action is None:
            self._last_action = np.zeros_like(a, dtype=np.float32)
        delta = a - self._last_action
        delta = np.clip(delta, -self._dwc, self._dwc)
        a_clipped = np.clip(self._last_action + delta, -1.0, 1.0)
        self._last_action = a_clipped.copy()
        return a_clipped

    def _cont_to_discrete(self, a_cont: np.ndarray) -> Dict[str, Any]:
        out = {}
        idx = 0
        for aid, bins in zip(self._agent_ids, self._per_agent_bins):
            choices = []
            for n in bins:
                if n <= 1:
                    choices.append(0)
                else:
                    v = a_cont[idx]
                    k = int(np.round((v + 1.0) * 0.5 * (n - 1)))
                    choices.append(np.clip(k, 0, n - 1))
                idx += 1
            out[aid] = choices[0] if len(choices) == 1 else tuple(choices)
        return out

    def _rews_array(self, rews: Dict[str, float]) -> np.ndarray:
        # Convert dict of per-agent rewards into a stable-order vector
        return np.array([float(rews.get(s, 0.0)) for s in self._sat_order], dtype=np.float32)

    def reset(self, *, seed=None, options=None):
        self._t = 0
        self._last_action = None
        obs_dict, _ = self.base.reset(seed=seed, options=options)
        obs = self._concat_obs(obs_dict)
        # IMPORTANT: return a stable info structure that will also appear in step()
        info = {"per_agent_reward": np.zeros(len(self._sat_order), dtype=np.float32)}
        return obs, info

    def step(self, action):
        self._t += 1
        a = np.asarray(action, dtype=np.float32)
        a = np.clip(a, -1.0, 1.0)
        a = self._apply_dwc(a)  # DWC
        act_dict = self._cont_to_discrete(a)

        obs_dict, rews, terms, truncs, _infos = self.base.step(act_dict)
        obs = self._concat_obs(obs_dict)

        # Aggregate reward (sum) + provide stable per-agent reward vector in info
        reward = float(sum(rews.values()))
        info = {"per_agent_reward": self._rews_array(rews)}

        terminated = bool(terms.get("__all__", False))
        truncated = bool(truncs.get("__all__", False))
        done = terminated or truncated
        return obs, reward, done, False, info

    @property
    def satellites(self):
        return self.base.unwrapped.satellites

    @property
    def agents(self):
        return self.base.agents

    @property
    def rewarder(self):
        return self.base.rewarder

    @property
    def remaining_tasks(self):
        return self.base.remaining_tasks

    @remaining_tasks.setter
    def remaining_tasks(self, v):
        self.base.remaining_tasks = v

# =========================================================
# Satellite args + orbit randomizer
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


# =========================================================
# Build Tianshou training objects
# =========================================================
def make_env():
    return JointConstellationEnv(max_episode_steps=64, dwc_delta=0.5)

train_envs = DummyVectorEnv([make_env for _ in range(1)])   # increase N for parallelism
test_envs  = DummyVectorEnv([make_env for _ in range(1)])

example_env = make_env()
state_shape = example_env.observation_space.shape
action_shape = example_env.action_space.shape

device = "cuda" if torch.cuda.is_available() else "cpu"

# Actor / Critic networks
net_a = Net(state_shape, hidden_sizes=[256, 256], device=device)
actor = Actor(net_a, action_shape, device=device, max_action=1.0).to(device)

net_c1 = Net(state_shape, action_shape, hidden_sizes=[256, 256], concat=True, device=device)
net_c2 = Net(state_shape, action_shape, hidden_sizes=[256, 256], concat=True, device=device)
critic1 = Critic(net_c1, device=device).to(device)
critic2 = Critic(net_c2, device=device).to(device)

actor_optim = torch.optim.Adam(actor.parameters(), lr=1e-4)
critic_optim = torch.optim.Adam(list(critic1.parameters()) + list(critic2.parameters()), lr=1e-3)

# Use an exploration noise OBJECT (not a float)
if GaussianNoise is None:
    raise ImportError("Could not import GaussianNoise from Tianshou. Please upgrade/install Tianshou properly.")
exp_noise = GaussianNoise(sigma=0.1)  # or OUNoise(sigma=0.2, theta=0.15)

policy = TD3Policy(
    actor=actor,
    actor_optim=actor_optim,
    critic1=critic1,
    critic1_optim=critic_optim,
    critic2=critic2,
    critic2_optim=critic_optim,
    tau=0.005,
    gamma=0.99,
    exploration_noise=exp_noise,   # <-- fixed: noise object
    policy_noise=0.2,
    update_actor_freq=2,
    noise_clip=0.5,
    estimation_step=1,
    action_space=example_env.action_space,
)

# Replay buffer (swap to HERReplayBuffer when you define goals)
buffer_size = int(1e6)
train_buffer = VectorReplayBuffer(buffer_size, len(train_envs))
# Example to enable HER later:
# from tianshou.data import HERReplayBuffer
# train_buffer = HERReplayBuffer(buffer_size, len(train_envs),
#                                reward_fn=your_her_reward_fn,
#                                sample_fn=your_her_sample_fn)

train_collector = Collector(policy, train_envs, train_buffer, exploration_noise=True)
test_collector  = Collector(policy, test_envs)

# Optional warm-up to stabilize replay stats
train_collector.collect(n_step=2048)


# =========================================================
# Training loop (off-policy style)
# =========================================================
def save_checkpoint(path_prefix="td3_hd"):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fname = f"{path_prefix}_{ts}.pth"
    torch.save(policy.state_dict(), fname)
    print(f"[TD3] Saved checkpoint: {fname}")
    return fname

epoch = 20
step_per_epoch = 10000
step_per_collect = 64
update_per_step = 1.0
batch_size = 512
test_num_episodes = 5

env_steps = 0
for ep in range(epoch):
    policy.train()
    steps_this_epoch = 0
    while steps_this_epoch < step_per_epoch:
        result = train_collector.collect(n_step=step_per_collect)
        steps_this_epoch += result["n/st"]
        env_steps += result["n/st"]
        # update TD3 (pulls minibatches from train_buffer internally)
        for _ in range(int(step_per_collect * update_per_step)):
            policy.update(batch_size, train_buffer)
    policy.eval()
    res = test_collector.collect(n_episode=test_num_episodes)
    print(f"[Epoch {ep}] env_steps={env_steps} "
          f"train_rew_mean={result['rews'].mean():.3f} "
          f"test_rew_mean={res['rews'].mean():.3f} "
          f"buffer={train_buffer.size}")
    if (ep + 1) % 5 == 0:
        save_checkpoint()

final_ckpt = save_checkpoint()


# =========================================================
# Test run with failure injections + plots
# =========================================================
test_env = make_env()

def inject_failure(sat_index: int):
    """Force a satellite 'dead' by zeroing power each step."""
    base = test_env.base
    sat = base.unwrapped.satellites[sat_index]
    def isnt_alive(log_failure=False, _sat=sat):
        death_message = messaging.PowerStorageStatusMsgPayload()
        death_message.storageLevel = 0.0
        _sat.dynamics.powerMonitor.batPowerOutMsg.write(death_message)
        return _sat.dynamics.is_alive(log_failure=log_failure) and _sat.fsw.is_alive(log_failure=log_failure)
    sat.is_alive = isnt_alive

max_steps = 40
max_wallclock_s = 120
fail_plan = {5: 0, 15: 2, 25: 3}

step_log = []
health_log = {f"Sat-{i}": [] for i in range(4)}
imaged_log = {f"Sat-{i}": [] for i in range(4)}
remaining_log = {f"Sat-{i}": [] for i in range(4)}
faulty_steps = {f"Sat-{i}": [] for i in range(4)}
remaining_matrix = np.zeros((4, max_steps))

obs, _ = test_env.reset()
test_env.base.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 6, "Sat-3": 3}
current_step = 0
episode_reward = 0.0
t0 = time.time()

policy.eval()
# turn off exploration noise for evaluation rollout
policy.set_exp_noise(None)

while current_step < max_steps and len(test_env.base.agents) > 0:
    if time.time() - t0 > max_wallclock_s:
        print(f"[Test] Wallclock timeout ({max_wallclock_s}s). Breaking to plots.")
        break

    # Tianshou policy forward -> Batch(act=...)
    with torch.no_grad():
        obs_tensor = torch.as_tensor(obs, dtype=torch.float32, device=device).unsqueeze(0)
        out = policy.forward(Batch(obs=obs_tensor), state=None, info=Batch(), exploration_noise=False)
        action = out.act[0].detach().cpu().numpy()

    obs, reward, terminated, truncated, info = test_env.step(action)
    episode_reward += float(reward)
    step_log.append(current_step)

    # Logging
    for i in range(4):
        sat_name = f"Sat-{i}"
        is_alive = 1 if sat_name in test_env.base.agents else 0
        health_log[sat_name].append(is_alive)
        imaged_count = getattr(test_env.base.rewarder, "imaged_by_sat", {}).get(sat_name, 0)
        imaged_log[sat_name].append(imaged_count)
        remaining = test_env.base.remaining_tasks.get(sat_name, 0)
        remaining_log[sat_name].append(remaining)
        remaining_matrix[i, current_step] = remaining
        if is_alive == 0:
            faulty_steps[sat_name].append(current_step)

    # Inject failures on schedule
    if current_step in fail_plan:
        idx = fail_plan[current_step]
        if 0 <= idx < 4 and (test_env.base.unwrapped.satellites[idx].name in test_env.base.agents):
            inject_failure(idx)

    current_step += 1
    if terminated or truncated:
        break

print(f"Test episode reward: {episode_reward:.3f}")

# ---------- Plots ----------
fig1, ax1 = plt.subplots()
for sat_name in health_log:
    ax1.plot(step_log, health_log[sat_name], label=sat_name)
ax1.set_xlabel('Step'); ax1.set_ylabel('Health (0=Faulty, 1=Healthy)')
ax1.set_title('Health Status Over Time'); ax1.legend()
plt.tight_layout(); plt.savefig('health_status.png'); plt.close(fig1)

fig2, axs = plt.subplots(4, 1, figsize=(9, 12), sharex=True)
for i, sat_name in enumerate(remaining_log):
    ax = axs[i]
    ax.plot(step_log, imaged_log[sat_name], label='Imaged')
    ax.plot(step_log, remaining_log[sat_name], label='Remaining')
    for fs in faulty_steps[sat_name]:
        ax.axvspan(fs - 0.5, fs + 0.5, alpha=0.3)
    ax.set_ylabel(f'{sat_name} (#)'); ax.legend()
axs[-1].set_xlabel('Step')
fig2.suptitle('Per-Satellite Task Progress (shaded = faulty steps)')
plt.tight_layout(); plt.savefig('task_progress.png'); plt.close(fig2)

fig3, ax3 = plt.subplots(figsize=(9, 3))
cax = ax3.imshow(remaining_matrix, aspect='auto', interpolation='nearest')
ax3.set_yticks(range(4)); ax3.set_yticklabels(['Sat-0', 'Sat-1', 'Sat-2', 'Sat-3'])
ax3.set_xlabel('Step'); ax3.set_title('Remaining Tasks Heatmap')
fig3.colorbar(cax, label='Remaining Tasks')
plt.tight_layout(); plt.savefig('remaining_tasks_heatmap.png'); plt.close(fig3)

print("Testing complete. Plots saved as:")
print(" - health_status.png")
print(" - task_progress.png")
print(" - remaining_tasks_heatmap.png")


# =========================================================
# Excel export
# =========================================================
def save_results_to_excel(step_log,
                          health_log,
                          imaged_log,
                          remaining_log,
                          remaining_matrix,
                          episode_reward,
                          algo_name="TD3_HD_Tianshou"):
    if pd is None:
        print("[Excel] pandas not installed; skipping Excel export."); return
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fname = f"results_{algo_name}_{ts}.xlsx"
    steps = len(step_log); sats = list(health_log.keys())
    rows = []
    for s in sats:
        for i, step in enumerate(step_log):
            rows.append({"step": step, "sat": s,
                         "health": health_log[s][i] if i < len(health_log[s]) else None,
                         "imaged": imaged_log[s][i] if i < len(imaged_log[s]) else None,
                         "remaining": remaining_log[s][i] if i < len(remaining_log[s]) else None})
    df_long = pd.DataFrame(rows)
    health_wide = pd.DataFrame({s: health_log[s][:steps] for s in sats}, index=step_log)
    imaged_wide = pd.DataFrame({s: imaged_log[s][:steps] for s in sats}, index=step_log)
    remaining_wide = pd.DataFrame({s: remaining_log[s][:steps] for s in sats}, index=step_log)
    hm = pd.DataFrame(remaining_matrix[:, :steps], index=sats, columns=step_log)
    summary = pd.DataFrame({"metric": ["episode_reward", "steps_logged", "sats"],
                            "value": [episode_reward, steps, len(sats)]})
    engine = "xlsxwriter"
    try:
        import xlsxwriter  # noqa: F401
    except Exception:
        engine = "openpyxl"
    with pd.ExcelWriter(fname, engine=engine) as writer:
        df_long.to_excel(writer, "log_long", index=False)
        health_wide.to_excel(writer, "health", index=True)
        imaged_wide.to_excel(writer, "imaged", index=True)
        remaining_wide.to_excel(writer, "remaining", index=True)
        hm.to_excel(writer, "heatmap", index=True)
        summary.to_excel(writer, "summary", index=False)
        if engine == "xlsxwriter":
            workbook = writer.book; ws = workbook.add_worksheet("plots")
            y = 2
            for img in ["health_status.png", "task_progress.png", "remaining_tasks_heatmap.png"]:
                if os.path.exists(img):
                    ws.insert_image(y, 2, img); y += 20
    print(f"Saved Excel: {fname}")

save_results_to_excel(
    step_log=step_log,
    health_log=health_log,
    imaged_log=imaged_log,
    remaining_log=remaining_log,
    remaining_matrix=remaining_matrix,
    episode_reward=episode_reward,
    algo_name="TD3_HD_Tianshou"
)
