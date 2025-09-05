

import time
import os
from datetime import datetime
from typing import Dict, Any

from bsk_rl.act.actions import ActionBuilder
import numpy as np
import matplotlib.pyplot as plt
import gymnasium as gym


# ---- Excel export deps (optional but recommended) ----
try:
    import pandas as pd
except Exception:
    pd = None

from ray.tune.registry import register_env
from ray.rllib.policy.policy import PolicySpec
from ray.rllib.env.multi_agent_env import MultiAgentEnv
from ray.rllib.algorithms.ppo import PPOConfig, PPO

from Basilisk.architecture import messaging
from collections import Counter  


from bsk_rl import ConstellationTasking
from bsk_rl.sats import ImagingSatellite
from bsk_rl.act import Action, Image
from bsk_rl import obs
from bsk_rl.sim import dyn, fsw
from bsk_rl.scene.targets import UniformTargets, Target
from bsk_rl.data import UniqueImageReward
from bsk_rl.comm import LOSCommunication
from bsk_rl.utils.orbital import walker_delta_args



class TargetAreas(UniformTargets):
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
# Custom action: Reallocate (kept simple/no-op reward shaping)
# =========================================================
class Reallocate(Action):
    def __init__(self, n_sats=4):
        self.action_space = gym.spaces.Discrete(n_sats)
        self.n_actions = self.action_space.n
        self.option = None

    @property
    def builder_type(self):
        # Match Image's builder; fallback to ActionBuilder if needed
        return Image.builder_type if isinstance(Image.builder_type, type) else ActionBuilder

    def set_action(self, option: int, **kwargs):
        self.option = option

    def action(self, satellite, state):
        # No direct side-effects here; environment bookkeeping handles tasks.
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


def _tune_access_generation(sat,
                            initial=600.0,   # ~30 min initial generation
                            step=60.0,       # 2-min stride
                            max_dur=1200.0):  # clamp lookahead to 1 hour
    """
    Tries common locations/attributes for the opportunity/access generator.
    Silently no-ops if attributes don't exist in your build.
    """
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

    # Optional extra limits (if present)
    extras = {"retask_on_image_complete": True, "max_access_compute_time_s": 2.0}
    for k, v in extras.items():
        if hasattr(og, k):
            setattr(og, k, v)


class CustomUniqueImageReward(UniqueImageReward):
    
    DEBUG_PROBE = False      # set True temporarily to inspect incoming events
    DEBUG_PROBE_STEPS = 5    # how many calls to print

    def __init__(self):
        # Many builds route imaging under one of these keys. Base class may
        # accept a subscription hint via data_store_kwargs (ignored if unknown).
        try:
            super().__init__(data_store_kwargs={"keys": ["imaged", "image", "image_complete"]})
        except TypeError:
            super().__init__()
        self.imaged_by_sat = {f"Sat-{i}": 0 for i in range(4)}
        self._probe_count = 0
        # Keep a set of target identifiers to prevent double-count
        self.imaged_targets = set()

    def _event_key(self, evt):
        # Try a handful of common field names
        for k in ("key", "event", "type", "name"):
            if hasattr(evt, k):
                v = getattr(evt, k)
                if isinstance(v, str):
                    return v
        # Sometimes events are dict-like
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
        # Derive a stable ID string and priority float from various shapes
        if tgt is None:
            return ("<none>", 0.0)
        # object with attributes
        tid = getattr(tgt, "name", None) or getattr(tgt, "id", None)
        pr = getattr(tgt, "priority", None)
        # dict-like
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
                    total += self.reward_fn(prio) / denom
                    new_unique_count += 1
            rewards[sat_id] = total
            if new_unique_count > 0:
                self.imaged_by_sat[sat_id] = self.imaged_by_sat.get(sat_id, 0) + new_unique_count

        # Ensure zeros for any known sat without entries this step
        for sat_id in getattr(self, "imaged_by_sat", {}).keys():
            rewards.setdefault(sat_id, 0.0)

        return rewards

class CustomConstellationTasking(ConstellationTasking, MultiAgentEnv):
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

        # Simple proxy: any positive reward => one task done
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
# Env factory + registration (includes speed tuning)
# =========================================================
def env_creator(env_config):
    max_episode_steps = env_config.get("max_episode_steps", 64)
    satellites = [AdvancedImagingSatellite(f"Sat-{i}", sat_args) for i in range(4)]

    # SPEED FIX: bound access lookahead and stride
    for sat in satellites:
        _tune_access_generation(sat, initial=1800.0, step=120.0, max_dur=3600.0)

    return CustomConstellationTasking(
        satellites=satellites,
        scenario=TargetAreas(n_targets=20),
        rewarder=CustomUniqueImageReward(),
        communicator=LOSCommunication(),
        sat_arg_randomizer=sat_arg_randomizer,
        log_level="INFO",
        max_episode_steps=max_episode_steps,
    )

register_env("custom_constellation", env_creator)


# =========================================================
# Probe spaces (for RLlib policy specs)
# =========================================================
_sample_env = env_creator({})
obs0, _ = _sample_env.reset()
first_agent = _sample_env.agents[0]
obs_space = _sample_env.observation_space(first_agent)
act_space = _sample_env.action_space(first_agent)
_sample_env.close()


# =========================================================
# PPO config + trainer (short horizons, small fragments)
# =========================================================
config = (
    PPOConfig()
    .api_stack(enable_rl_module_and_learner=False, enable_env_runner_and_connector_v2=False)
    .environment("custom_constellation", env_config={"max_episode_steps": 64})
    .multi_agent(
        policies={"shared_policy": PolicySpec(None, obs_space, act_space)},
        policy_mapping_fn=lambda agent_id, *args, **kwargs: "shared_policy",
    )
    .framework("torch")
    .resources(num_gpus=0)
    .env_runners(num_env_runners=1, rollout_fragment_length=16, sample_timeout_s=120.0)
    .training(
        gamma=0.99,
        lr=5e-4,
        train_batch_size=2048,
        vf_clip_param=10.0,
        clip_param=0.2,
        num_sgd_iter=10,
        lambda_=0.95,
    )
)

algo = PPO(config=config)



num_iterations = 10
for i in range(num_iterations):
    result = algo.train()
    print(f"Iteration {i}: mean reward = {result.get('episode_reward_mean', 'N/A')}")

checkpoint_dir = algo.save()
print(f"Checkpoint saved at {checkpoint_dir}")


# =========================================================
# Test run with failure injections + wall-clock failsafe + plotting
# =========================================================
test_env = env_creator({})

def inject_failure(sat_index: int):
    """Force a satellite 'dead' by zeroing power each step."""
    sat = test_env.unwrapped.satellites[sat_index]
    def isnt_alive(log_failure=False, _sat=sat):
        death_message = messaging.PowerStorageStatusMsgPayload()
        death_message.storageLevel = 0.0
        _sat.dynamics.powerMonitor.batPowerOutMsg.write(death_message)
        return _sat.dynamics.is_alive(log_failure=log_failure) and _sat.fsw.is_alive(log_failure=log_failure)
    sat.is_alive = isnt_alive

max_steps = 20
max_wallclock_s = 30  # hard cap to ensure we reach plots
fail_plan = {5: 0, 15: 2, 25: 3}  # step -> satellite index

# Logs
step_log = []
health_log = {sat.name: [] for sat in test_env.unwrapped.satellites}
imaged_log = {sat.name: [] for sat in test_env.unwrapped.satellites}
remaining_log = {sat.name: [] for sat in test_env.unwrapped.satellites}
faulty_steps = {sat.name: [] for sat in test_env.unwrapped.satellites}
remaining_matrix = np.zeros((4, max_steps))

# Reset
observations, _ = test_env.reset()
test_env.remaining_tasks = {"Sat-0": 3, "Sat-1": 3, "Sat-2": 6, "Sat-3": 3}

current_step = 0
episode_reward = 0.0
t0 = time.time()

while current_step < max_steps and test_env.agents:
    if time.time() - t0 > max_wallclock_s:
        print(f"[Test] Wallclock timeout ({max_wallclock_s}s). Breaking to plots.")
        break

    # Actions from policy (fine despite deprecation warning)
    actions = {}
    for agent in test_env.agents:
        pid = algo.config.policy_mapping_fn(agent)
        action = algo.compute_single_action(observations[agent], policy_id=pid)
        actions[agent] = action

    observations, rewards, terminations, truncations, infos = test_env.step(actions)
    episode_reward += sum(rewards.values())

    step_log.append(current_step)

    # Logging
    for i, sat in enumerate(test_env.unwrapped.satellites):
        sat_name = sat.name
        is_alive = 1 if sat_name in test_env.agents else 0
        health_log[sat_name].append(is_alive)

        imaged_count = getattr(test_env.rewarder, "imaged_by_sat", {}).get(sat_name, 0)
        imaged_log[sat_name].append(imaged_count)

        remaining = test_env.remaining_tasks.get(sat_name, 0)
        remaining_log[sat_name].append(remaining)
        remaining_matrix[i, current_step] = remaining

        if is_alive == 0:
            faulty_steps[sat_name].append(current_step)

    # Inject failures
    if current_step in fail_plan:
        idx = fail_plan[current_step]
        if 0 <= idx < len(test_env.unwrapped.satellites):
            if test_env.unwrapped.satellites[idx].name in test_env.agents:
                inject_failure(idx)

    current_step += 1

print(f"Test episode reward: {episode_reward}")

# =========================================================
# Plots
# =========================================================
# 1) Health over time
fig1, ax1 = plt.subplots()
for sat_name in health_log:
    ax1.plot(step_log, health_log[sat_name], label=sat_name)
ax1.set_xlabel('Step')
ax1.set_ylabel('Health (0=Faulty, 1=Healthy)')
ax1.set_title('Health Status Over Time')
ax1.legend()
plt.tight_layout()
plt.savefig('health_status.png')
plt.close(fig1)

# 2) Per-satellite task progress with shaded faulty steps
fig2, axs = plt.subplots(4, 1, figsize=(9, 12), sharex=True)
for i, sat_name in enumerate(remaining_log):
    ax = axs[i]
    ax.plot(step_log, imaged_log[sat_name], label='Imaged')
    ax.plot(step_log, remaining_log[sat_name], label='Remaining')
    for fs in faulty_steps[sat_name]:
        ax.axvspan(fs - 0.5, fs + 0.5, alpha=0.3)
    ax.set_ylabel(f'{sat_name} (#)')
    ax.legend()
axs[-1].set_xlabel('Step')
fig2.suptitle('Per-Satellite Task Progress (shaded = faulty steps)')
plt.tight_layout()
plt.savefig('task_progress.png')
plt.close(fig2)

# 3) Remaining tasks heatmap
fig3, ax3 = plt.subplots(figsize=(9, 3))
cax = ax3.imshow(remaining_matrix, aspect='auto', interpolation='nearest')
ax3.set_yticks(range(4))
ax3.set_yticklabels(['Sat-0', 'Sat-1', 'Sat-2', 'Sat-3'])
ax3.set_xlabel('Step')
ax3.set_title('Remaining Tasks Heatmap')
fig3.colorbar(cax, label='Remaining Tasks')
plt.tight_layout()
plt.savefig('remaining_tasks_heatmap.png')
plt.close(fig3)

print("Testing complete. Plots saved as:")
print(" - health_status.png")
print(" - task_progress.png")
print(" - remaining_tasks_heatmap.png")


# =========================================================
# Excel export (tidy + wide tables; embeds plots if xlsxwriter present)
# =========================================================
def save_results_to_excel(step_log,
                          health_log,
                          imaged_log,
                          remaining_log,
                          remaining_matrix,
                          episode_reward,
                          algo_name="PPO"):
    if pd is None:
        print("[Excel] pandas not installed; skipping Excel export.")
        return

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fname = f"results_{algo_name}_{ts}.xlsx"

    steps = len(step_log)
    sats = list(health_log.keys())

    # Long/tidy table
    rows = []
    for s in sats:
        for i, step in enumerate(step_log):
            rows.append({
                "step": step,
                "sat": s,
                "health": health_log[s][i] if i < len(health_log[s]) else None,
                "imaged": imaged_log[s][i] if i < len(imaged_log[s]) else None,
                "remaining": remaining_log[s][i] if i < len(remaining_log[s]) else None,
            })
    df_long = pd.DataFrame(rows)

    # Wide tables
    health_wide = pd.DataFrame({s: health_log[s][:steps] for s in sats}, index=step_log)
    imaged_wide = pd.DataFrame({s: imaged_log[s][:steps] for s in sats}, index=step_log)
    remaining_wide = pd.DataFrame({s: remaining_log[s][:steps] for s in sats}, index=step_log)

    # Heatmap matrix (rows=sats, cols=steps actually logged)
    hm = pd.DataFrame(remaining_matrix[:, :steps], index=sats, columns=step_log)

    summary = pd.DataFrame(
        {"metric": ["episode_reward", "steps_logged", "sats"],
         "value": [episode_reward, steps, len(sats)]}
    )

    # Choose engine; only xlsxwriter can embed images
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

        # Optionally embed the plots if engine supports it and files exist
        if engine == "xlsxwriter":
            workbook = writer.book
            ws = workbook.add_worksheet("plots")
            y = 2
            for img in ["health_status.png", "task_progress.png", "remaining_tasks_heatmap.png"]:
                if os.path.exists(img):
                    ws.insert_image(y, 2, img)
                    y += 20

    print(f"Saved Excel: {fname}")

# Export
save_results_to_excel(
    step_log=step_log,
    health_log=health_log,
    imaged_log=imaged_log,
    remaining_log=remaining_log,
    remaining_matrix=remaining_matrix,
    episode_reward=episode_reward,
    algo_name="PPO"
)

