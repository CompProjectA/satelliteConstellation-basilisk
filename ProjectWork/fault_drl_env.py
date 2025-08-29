# fault_drl_env.py
# Gym environment for DRL-based fault detection on a single Basilisk satellite.
# Action = {0: no_fault, 1: battery, 2: friction, 3: powerlimit, 4: encoder}
# Observation = compact telemetry vector (battery + nav/RW if you wire them)
# Reward = online classification: +1 correct after onset, penalties for false alarms & delay.

import math
import random
from dataclasses import dataclass
from typing import Dict, Any, Optional, Tuple, List

import gymnasium as gym
import numpy as np
from gymnasium import spaces

# === Basilisk imports (assumes Basilisk installed) ===
from Basilisk.utilities import (SimulationBaseClass, macros, orbitalMotion,
                                simIncludeGravBody)
from Basilisk.simulation import spacecraft, simpleNav, simpleBattery, simplePowerSink
from Basilisk.architecture import messaging
from Basilisk.architecture import bskLogging
bskLogging.setDefaultLogLevel(bskLogging.BSK_WARNING)
# === Your fault modules (you uploaded these) ===
# Wire these in the adapters below (see FaultRegistry).
from BskSim.faults import battery_fault as bf_mod
from BskSim.faults import friction_fault as ff_mod
from BskSim.faults import powerlimit_fault as pf_mod
from BskSim.faults import encoder_fault as ef_mod

# ----------------------------
# Config dataclasses
# ----------------------------
@dataclass
class FaultConfig:
    name: str
    code: int  # action label
    # Injection time window (seconds) for single-fault-per-episode mode:
    t_start_min_s: float = 60.0
    t_start_max_s: float = 600.0

@dataclass
class EnvConfig:
    step_dt_s: float = 1.0
    episode_len_s: float = 1200.0       # 20 minutes
    single_fault_per_episode: bool = True  # False => schedule all 4 faults in one ep
    include_no_fault_episodes: bool = True
    no_fault_prob: float = 0.15          # chance to run an episode with no faults
    # Reward shaping
    r_correct: float = 1.0               # per-step reward after onset if correct
    r_false_alarm: float = -0.20         # per-step penalty if predicts fault before onset
    r_wrong_label: float = -0.50         # per-step penalty after onset if wrong label
    r_missed_after_window: float = -1.0  # penalty every step after onset+grace if still wrong
    detection_grace_s: float = 5.0       # grace window after onset before harsh miss penalty
    seed: Optional[int] = None


# ----------------------------
# Fault Adapters (thin wrappers)
# ----------------------------
class FaultAdapter:
    """Base adapter. You will implement 'inject' using your module."""
    def __init__(self, name: str):
        self.name = name
        self.is_applied = False

    def schedule(self, t_start_s: float):
        self.t_start_s = t_start_s
        self.is_applied = False

    def maybe_apply(self, t_now_s: float, ctx: Dict[str, Any]):
        if (not self.is_applied) and (t_now_s >= self.t_start_s):
            self.inject(ctx)  # <-- call into your module here
            self.is_applied = True

    def inject(self, ctx: Dict[str, Any]):
        raise NotImplementedError


class BatteryFault(FaultAdapter):
    def inject(self, ctx: Dict[str, Any]):
        """
        Example strategies (pick ONE and comment out the others):
          1) Increase a fault sink draw (W -> kW):
             ctx['sink_fault'].nodePowerOut = -(50.0 / 1000.0)
          2) Degrade capacity multiplicatively:
             ctx['battery'].storageCapacity *= 0.9
          3) Call your own function:
             bf_mod.activate_fault(ctx['scSim'], ctx['battery'])
        """
        battery = ctx['battery']
        # If you already created a dedicated fault sink earlier:
        if 'sink_fault' in ctx and ctx['sink_fault'] is not None:
            ctx['sink_fault'].nodePowerOut = -(50.0 / 1000.0)  # +50 W
        else:
            # Capacity fade as a generic fallback
            battery.storageCapacity = max(1e3, float(battery.storageCapacity) * 0.9)
        # If your module provides a function, uncomment:
        # bf_mod.inject(ctx['scSim'], battery)


class FrictionFault(FaultAdapter):
    def inject(self, ctx: Dict[str, Any]):
        """
        Typical: introduce unmodeled external torque or increase RW friction.
        Replace the line below with a call into your friction_fault.py.
        Example:
            ff_mod.inject(ctx['scSim'], ctx['sc'], axis=0, tau=1e-4)  # Nm bias
        """
        if hasattr(ff_mod, "inject"):
            ff_mod.inject(ctx['scSim'], ctx['sc'])
        else:
            # Generic placeholder: do nothing (safe)
            pass


class PowerLimitFault(FaultAdapter):
    def inject(self, ctx: Dict[str, Any]):
        """
        Enforce an artificial bus power cap or reduce solar input.
        Example (if your module exposes it):
            pf_mod.inject(ctx['scSim'], ctx['battery'], maxPower_W=20.0)
        """
        if hasattr(pf_mod, "inject"):
            pf_mod.inject(ctx['scSim'], ctx['battery'])
        else:
            # Generic placeholder: reduce capacity hard to emulate limit effects
            cap = float(ctx['battery'].storageCapacity)
            ctx['battery'].storageCapacity = max(1e3, 0.5 * cap)


class EncoderFault(FaultAdapter):
    def inject(self, ctx: Dict[str, Any]):
        """
        Typical: corrupt RW encoder (speed/angle) or freeze one wheel’s encoder.
        Example:
            ef_mod.inject(ctx['scSim'], ctx['rw_model'], wheel=2, mode="stuck")
        """
        if hasattr(ef_mod, "inject"):
            ef_mod.inject(ctx['scSim'], ctx['sc'])
        else:
            pass


# Registry to map numeric label -> adapter
def build_fault_registry() -> Dict[int, FaultAdapter]:
    return {
        1: BatteryFault("battery"),
        2: FrictionFault("friction"),
        3: PowerLimitFault("powerlimit"),
        4: EncoderFault("encoder"),
    }


# ----------------------------
# The Gym environment
# ----------------------------
class SingleSatFaultEnv(gym.Env):
    metadata = {"render.modes": ["human"]}

    def __init__(self, cfg: EnvConfig = EnvConfig()):
        super().__init__()
        self.cfg = cfg
        if cfg.seed is not None:
            random.seed(cfg.seed)
            np.random.seed(cfg.seed)

        # --- Basilisk sim + single satellite + battery + optional fault sink ---
        self.scSim, self.sc, self.battery, self.sink_base, self.sink_fault, self.nav = self._build_sim()

        # --- Message readers for observations ---
        self.batt_reader = messaging.PowerStorageStatusMsgReader()
        self.batt_reader.subscribeTo(self.battery.batPowerOutMsg)

        self.nav_reader = messaging.NavAttMsgReader()
        self.nav_reader.subscribeTo(self.nav.attOutMsg)

        # --- Action/Observation spaces ---
        # action: predict current fault label [0..4]: 0=no-fault, 1=battery, 2=friction, 3=powerlimit, 4=encoder
        self.action_space = spaces.Discrete(5)

        # observation: [soc_pct, netPower_kW, storedCharge_kJ, omega_norm]
        low = np.array([0.0, -10.0, 0.0, 0.0], dtype=np.float32)
        high = np.array([100.0, 10.0, 1e6, 1e2], dtype=np.float32)
        self.observation_space = spaces.Box(low=low, high=high, dtype=np.float32)

        # --- Fault setup ---
        self.faults: Dict[int, FaultAdapter] = build_fault_registry()
        self.fault_list: List[FaultConfig] = [
            FaultConfig("battery", 1),
            FaultConfig("friction", 2),
            FaultConfig("powerlimit", 3),
            FaultConfig("encoder", 4),
        ]

        self._reset_episode_state()

    # ----------------------------
    # Basilisk construction
    # ----------------------------
    def _build_sim(self):
        scSim = SimulationBaseClass.SimBaseClass()
        simTaskName = "simTask"
        dynProcess = scSim.CreateNewProcess("simProcess")
        dt = macros.sec2nano(1.0)
        dynProcess.addTask(scSim.CreateNewTask(simTaskName, dt))

        # Gravity + a simple LEO orbit
        gravFactory = simIncludeGravBody.gravBodyFactory()
        earth = gravFactory.createEarth()
        earth.isCentralBody = True
        mu = earth.mu

        oe = orbitalMotion.ClassicElements()
        oe.a = 7000e3
        oe.e = 0.0
        oe.i = 56 * macros.D2R
        oe.Omega = 0.0
        oe.omega = 0.0
        oe.f = 45 * macros.D2R
        rN, vN = orbitalMotion.elem2rv(mu, oe)

        # Spacecraft
        sc = spacecraft.Spacecraft()
        sc.ModelTag = "SingleSat"
        sc.hub.mHubMass = 10.0
        sc.hub.mHubI = [0.1, 0.1, 0.1]
        sc.hub.r_CN_NInit = rN
        sc.hub.v_CN_NInit = vN
        scSim.AddModelToTask(simTaskName, sc)

        # Simple nav (for angular rates)
        nav = simpleNav.SimpleNav()
        nav.ModelTag = "nav"
        nav.scStateInMsg.subscribeTo(sc.scStateOutMsg)
        scSim.AddModelToTask(simTaskName, nav)

        # Battery + loads
        battery = simpleBattery.SimpleBattery()
        battery.ModelTag = "bat"
        battery.storageCapacity = 3.6e6      # 1 kWh (J)
        battery.storedCharge_Init = 1.8e6    # 50% SoC
        scSim.AddModelToTask(simTaskName, battery)

        sink_base = simplePowerSink.SimplePowerSink()
        sink_base.ModelTag = "sink_base"
        sink_base.nodePowerOut = -0.01      # 10W base
        scSim.AddModelToTask(simTaskName, sink_base)
        battery.addPowerNodeToModel(sink_base.nodePowerOutMsg)

        sink_fault = simplePowerSink.SimplePowerSink()
        sink_fault.ModelTag = "sink_fault"
        sink_fault.nodePowerOut = 0.0       # off initially
        scSim.AddModelToTask(simTaskName, sink_fault)
        battery.addPowerNodeToModel(sink_fault.nodePowerOutMsg)

        scSim.InitializeSimulation()

        return scSim, sc, battery, sink_base, sink_fault, nav

    # ----------------------------
    # Gym API
    # ----------------------------
    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        super().reset(seed=seed)
        if seed is not None:
            random.seed(seed); np.random.seed(seed)

        # Rebuild sim fresh each episode (simplest & robust)
        self.scSim, self.sc, self.battery, self.sink_base, self.sink_fault, self.nav = self._build_sim()
        self.batt_reader = messaging.PowerStorageStatusMsgReader()
        self.batt_reader.subscribeTo(self.battery.batPowerOutMsg)
        self.nav_reader = messaging.NavAttMsgReader()
        self.nav_reader.subscribeTo(self.nav.attOutMsg)

        self._reset_episode_state()
        return self._observe(), {}

    def step(self, action: int):
        # Advance simulation by cfg.step_dt_s
        self.t_s += self.cfg.step_dt_s
        nextStop = macros.sec2nano(self.t_s)
        self.scSim.ConfigureStopTime(nextStop)
        self.scSim.ExecuteSimulation()

        # Maybe apply scheduled faults
       # Maybe apply scheduled faults
        if self.cfg.single_fault_per_episode:
            target = self.episode_fault
            if target.code == 0:  # <-- no-fault episode
                true_label = 0
                onset_s = float("inf")  # or self.fault_onset_s
            else:
                adapter = self.faults[target.code]
                adapter.maybe_apply(self.t_s, self.ctx)
                true_label = target.code if adapter.is_applied else 0
                onset_s = self.fault_onset_s
        else:
            # Multiple faults in the same episode, label = latest active (or 0)
            true_label = 0
            onset_s = None
            for f, onset in zip(self.multi_faults, self.multi_onsets):
                ad = self.faults[f.code]
                ad.maybe_apply(self.t_s, self.ctx)
                if ad.is_applied:
                    true_label = f.code
                    onset_s = onset  # latest applied becomes the label


        # Compute reward
        reward = 0.0
        done = self.t_s >= self.cfg.episode_len_s

        if true_label == 0:
            # No fault active yet
            if action == 0:
                reward += 0.0   # neutral before onset
            else:
                reward += self.cfg.r_false_alarm
        else:
            # Fault active
            if action == true_label:
                reward += self.cfg.r_correct
                # small early-detection bonus for being quick
                if (self.t_s - onset_s) <= self.cfg.detection_grace_s:
                    reward += 0.25
            else:
                reward += self.cfg.r_wrong_label
                if (self.t_s - onset_s) > self.cfg.detection_grace_s:
                    reward += self.cfg.r_missed_after_window

        obs = self._observe()
        info = {
            "t_s": self.t_s,
            "true_label": true_label,
            "fault_onset_s": onset_s,
            "episode_fault": (self.episode_fault.name if self.cfg.single_fault_per_episode else "multi"),
        }
        return obs, reward, done, False, info

    # ----------------------------
    # Helpers
    # ----------------------------
    def _reset_episode_state(self):
        self.t_s = 0.0

        # Context shared with adapters
        self.ctx: Dict[str, Any] = {
            "scSim": self.scSim,
            "sc": self.sc,
            "battery": self.battery,
            "sink_fault": self.sink_fault,
            # add more here if your fault modules need other handles (RW model, controller, etc.)
        }
        # Reset all faults
        for ad in self.faults.values():
            ad.is_applied = False
            ad.t_start_s = float("inf")

        if self.cfg.single_fault_per_episode:
            # Maybe run a no-fault episode
            if self.cfg.include_no_fault_episodes and random.random() < self.cfg.no_fault_prob:
                self.episode_fault = FaultConfig("none", 0)
                self.fault_onset_s = float("inf")
            else:
                self.episode_fault = random.choice(self.fault_list)
                self.fault_onset_s = random.uniform(self.episode_fault.t_start_min_s,
                                                    self.episode_fault.t_start_max_s)
                self.faults[self.episode_fault.code].schedule(self.fault_onset_s)
            self.multi_faults = []
            self.multi_onsets = []
        else:
            # Schedule all 4 in one episode at spaced times
            start = 120.0
            gap = (self.cfg.episode_len_s - 240.0) / 4.0  # spread across ep
            self.multi_faults = self.fault_list[:]
            self.multi_onsets = []
            for i, f in enumerate(self.multi_faults):
                onset = start + i * gap
                self.faults[f.code].schedule(onset)
                self.multi_onsets.append(onset)

    def _observe(self) -> np.ndarray:
        # Battery features
        batt = self.batt_reader()
        # Common fields in PowerStorageStatus: storageLevel (%), netPower (kW), storedCharge (J)
        soc_pct = float(getattr(batt, "storageLevel", 100.0))
        net_kW = float(getattr(batt, "netPower", 0.0))
        stored_J = float(getattr(batt, "storedCharge", 0.0))
        stored_kJ = stored_J / 1000.0

        # Nav features (angular rate magnitude)
        nav = self.nav_reader()
        omega_BN_B = np.array([nav.omega_BN_B[0], nav.omega_BN_B[1], nav.omega_BN_B[2]], dtype=np.float64)
        omega_norm = float(np.linalg.norm(omega_BN_B))

        obs = np.array([soc_pct, net_kW, stored_kJ, omega_norm], dtype=np.float32)
        return obs

    # Optional: simple print render
    def render(self, mode="human"):
        pass

    def close(self):
        pass


# Quick smoke test
if __name__ == "__main__":
    env = SingleSatFaultEnv()
    o, _ = env.reset()
    total = 0
    for _ in range(50):
        a = env.action_space.sample()
        o, r, d, t, info = env.step(a)
        total += r
        if d:
            break
    print("smoke reward:", total)
