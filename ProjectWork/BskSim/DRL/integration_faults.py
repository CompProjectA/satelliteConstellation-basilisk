"""
Integration.py
--------------
Real-time fault detection + recovery for a Basilisk simulation using your
trained Keras anomaly detector (autoencoder).

Pipeline per step:
  1) Get current spacecraft state from Envs.BasiliskModel
  2) Compute baseline control torques (PD by default)
  3) Build feature window -> scale -> model.predict -> reconstruction error
  4) If anomaly: classify likely fault (battery | powerlimit | encoder | friction)
  5) Apply tailored recovery (safe torque clamp, zero affected wheels)
  6) Step the sim forward

Files this uses:
  - Envs.py .......... BasiliskModel wrapper (you already have)
  - PD.py ............ PDControllor(error_MRP, omega)
  - anomaly_detection_model.keras .... your trained Keras model
  - (optional) simulation_log.xlsx .... to fit MinMaxScaler; otherwise we auto-fit online

Config via environment variables:
  MODEL_PATH=anomaly_detection_model.keras
  SCALER_DATA=simulation_log.xlsx
  WINDOW=10
  THRESHOLD=0.185
  DT_SEC=0.08
  STEPS=5000
  SAFE_TORQUE_LIMIT=0.02
  RECOVERY_STEPS=200

Detection/classification thresholds (battery/power/encoder/friction) are also tunable via env vars below.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler
import os

# quiet noisy TF logs (optional)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# avoid protobuf C++ impl contention
os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"

# keep thread pools tiny to prevent abseil/grpc lock contention
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["TF_NUM_INTEROP_THREADS"] = "1"
os.environ["TF_NUM_INTRA_OP_THREADS"] = "1"
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"  # disable oneDNN thread pooling (often helps)

import tensorflow as tf
tf.config.threading.set_inter_op_parallelism_threads(1)
tf.config.threading.set_intra_op_parallelism_threads(1)
from tensorflow.keras.models import load_model

from Basilisk.utilities import macros, RigidBodyKinematics
from Envs import BasiliskModel
from PD import PDControllor

# ---------------------------- Config ----------------------------
MODEL_PATH =  "anomaly_detection_mode.keras"
SCALER_DATA = os.environ.get("SCALER_DATA", "simulation_log.xlsx")

WINDOW = int(os.environ.get("WINDOW", 10))
THRESHOLD = float(os.environ.get("THRESHOLD", 0.185))

DT_SEC = float(os.environ.get("DT_SEC", 0.08))
STEPS = int(os.environ.get("STEPS", 1000))

SAFE_TORQUE_LIMIT = float(os.environ.get("SAFE_TORQUE_LIMIT", 0.02))
RECOVERY_STEPS = int(os.environ.get("RECOVERY_STEPS", 200))

# Battery thresholds
BATTERY_MIN_LEVEL = float(os.environ.get("BATTERY_MIN_LEVEL", "300000"))         # W·s
BATTERY_NET_NEG_THRESH = float(os.environ.get("BATTERY_NET_NEG_THRESH", "-40")) # W
BATTERY_NEG_PERSIST_STEPS = int(os.environ.get("BATTERY_NEG_PERSIST_STEPS", "40"))

# Power-limit thresholds
PER_WHEEL_POWER_LIMIT = float(os.environ.get("PER_WHEEL_POWER_LIMIT", "60.0"))  # W
TOTAL_POWER_LIMIT = float(os.environ.get("TOTAL_POWER_LIMIT", "180.0"))         # W
POWER_PERSIST_STEPS = int(os.environ.get("POWER_PERSIST_STEPS", "10"))

# Encoder thresholds
TORQUE_MIN_FOR_ENCODER = float(os.environ.get("TORQUE_MIN_FOR_ENCODER", "0.02"))  # N·m
SPEED_DELTA_MIN = float(os.environ.get("SPEED_DELTA_MIN", "1.0"))                 # rad/s
ENCODER_STUCK_STEPS = int(os.environ.get("ENCODER_STUCK_STEPS", "20"))

# Friction thresholds (electrical vs mechanical power gap)
FRICTION_ALPHA = float(os.environ.get("FRICTION_ALPHA", "1.0"))   # expected P_elec ~ alpha*|τ|*|ω|
FRICTION_MARGIN = float(os.environ.get("FRICTION_MARGIN", "8.0")) # W extra over expected
FRICTION_PERSIST_STEPS = int(os.environ.get("FRICTION_PERSIST_STEPS", "20"))

# The feature vector your model expects (matches your earlier Integration.py)
FEATURE_COLUMNS = [
    "Torque_RW0", "Torque_RW1", "Torque_RW2", "Torque_RW3",
    "MRP_x", "MRP_y", "MRP_z",
    "Error_MRP_x", "Error_MRP_y", "Error_MRP_z",
    "Omega_x", "Omega_y", "Omega_z",
    "RW_Power_0", "RW_Power_1", "RW_Power_2", "RW_Power_3",
]

# ---------------------------- Helpers ----------------------------
def latest_row(x):
    """Return last row/element as 1D numpy array or None."""
    try:
        a = np.array(x)
        if a.ndim == 0:
            return np.array([float(a)])
        return np.array(a[-1]).reshape(-1)
    except Exception:
        return None
def vec4(x, fill=0.0):
    """Return a 1D float array of length 4 (pad or trim as needed)."""
    a = np.asarray(x, dtype=float).reshape(-1)
    if a.size < 4:
        a = np.pad(a, (0, 4 - a.size), mode="constant", constant_values=fill)
    elif a.size > 4:
        a = a[:4]
    return a

def get_rw_speeds(model):
    try:
        return latest_row(model.rwSpeedLog.wheelSpeeds)[:4].astype(float)
    except Exception:
        return np.zeros(4, float)

def get_rw_powers(model):
    vals = np.zeros(4)
    try:
        for i in range(4):
            vals[i] = float(latest_row(model.rwPowLog[i].netPower)[-1])
    except Exception:
        pass
    return vals

def get_battery_metrics(model):
    storage, net = None, None
    try:
        storage = float(latest_row(model.batPowLog.storageLevel)[-1])
        net = float(latest_row(model.batPowLog.currentNetPower)[-1])
    except Exception:
        pass
    return storage, net

def build_feature_row(torque4, cur_MRP3, err_MRP3, omega3, rw_power4):
    return np.concatenate([
        torque4[:4], cur_MRP3[:3], err_MRP3[:3], omega3[:3],
        (rw_power4[:4] if rw_power4 is not None else np.zeros(4))
    ]).astype(np.float32)

# ---------------------------- Rule-based classifiers ----------------------------
class Counters:
    def __init__(self):
        self.batt_neg = 0
        self.power = 0
        self.encoder = np.zeros(4, int)
        self.friction = np.zeros(4, int)

def classify_battery(storage_Ws, net_W, C: Counters):
    fault, reason = False, None
    if storage_Ws is not None and storage_Ws < BATTERY_MIN_LEVEL:
        return True, f"battery_low({storage_Ws:.0f}Ws<{BATTERY_MIN_LEVEL})"
    if net_W is not None and net_W < BATTERY_NET_NEG_THRESH:
        C.batt_neg += 1
    else:
        C.batt_neg = 0
    if C.batt_neg >= BATTERY_NEG_PERSIST_STEPS:
        fault, reason = True, f"battery_net_negative({net_W:.1f} W) for {C.batt_neg} steps"
    return fault, reason

def classify_powerlimit(powers4, C: Counters):
    total = float(np.sum(powers4))
    over = (powers4 > PER_WHEEL_POWER_LIMIT).any() or (total > TOTAL_POWER_LIMIT)
    C.power = C.power + 1 if over else 0
    if C.power >= POWER_PERSIST_STEPS:
        why = "per-wheel>limit" if (powers4 > PER_WHEEL_POWER_LIMIT).any() else "total>limit"
        return True, f"power_limit({why})"
    return False, None

def classify_encoder(torque4, speeds4, prev_speeds4, C: Counters):
    if prev_speeds4 is None:
        return [], []
    dω = speeds4 - prev_speeds4
    wheels, reasons = [], []
    for k in range(4):
        stuck = abs(torque4[k]) > TORQUE_MIN_FOR_ENCODER and abs(dω[k]) < SPEED_DELTA_MIN
        C.encoder[k] = C.encoder[k] + 1 if stuck else 0
        if C.encoder[k] >= ENCODER_STUCK_STEPS:
            wheels.append(k)
            reasons.append(f"encoder_stuck(w{k})")
    return wheels, reasons

def classify_friction(torque4, speeds4, powers4, C: Counters):
    mech = np.abs(torque4[:4] * speeds4[:4])
    expected = FRICTION_ALPHA * mech
    excess = powers4[:4] - expected
    wheels, reasons = [], []
    for k in range(4):
        high = excess[k] > FRICTION_MARGIN
        C.friction[k] = C.friction[k] + 1 if high else 0
        if C.friction[k] >= FRICTION_PERSIST_STEPS:
            wheels.append(k)
            reasons.append(f"friction_high_power(w{k})")
    return wheels, reasons

# ---------------------------- Recovery logic ----------------------------
class Recovery:
    def __init__(self):
        self.active = False
        self.steps = 0
        self.mode = None         # 'battery' | 'powerlimit' | 'encoder' | 'friction' | 'unknown'
        self.wheels = []         # affected wheels for encoder/friction

    def trigger(self, mode, wheels=None):
        self.active = True
        self.steps = RECOVERY_STEPS
        self.mode = mode
        self.wheels = list(wheels or [])

    def tick(self):
        if self.active:
            self.steps -= 1
            if self.steps <= 0:
                self.active = False
                self.mode = None
                self.wheels = []

def apply_recovery(mode, torque4, wheels):
    t = torque4.copy()
    if mode in ("battery", "powerlimit", "unknown"):
        # Global throttle
        t[:4] = np.clip(t[:4], -SAFE_TORQUE_LIMIT, SAFE_TORQUE_LIMIT)
    elif mode in ("encoder", "friction"):
        # Zero affected wheels; clamp the rest
        for w in wheels:
            if 0 <= w < 4:
                t[w] = 0.0
        t[:4] = np.clip(t[:4], -SAFE_TORQUE_LIMIT, SAFE_TORQUE_LIMIT)
    return t

# ---------------------------- Main integration ----------------------------
def main():

    # Load anomaly model
    model = load_model("anomaly_detection_model.keras")
    print(f"[Integration] Loaded anomaly model: {MODEL_PATH}")

    # Fit scaler
    scaler = MinMaxScaler()
    if os.path.exists(SCALER_DATA):
        try:
            df = pd.read_excel(SCALER_DATA)
            cols = [c for c in df.columns if c in FEATURE_COLUMNS]
            for c in FEATURE_COLUMNS:
                if c not in cols:
                    df[c] = 0.0
            scaler.fit(df[FEATURE_COLUMNS])
            print(f"[Integration] Scaler fitted on {SCALER_DATA} (columns={FEATURE_COLUMNS})")
        except Exception as e:
            print(f"[Integration] Warning: failed to fit scaler from {SCALER_DATA}: {e}\n"
                  f"            Falling back to online fit on first window.")
            scaler = None
    else:
        print("[Integration] No scaler data found; will online-fit on first window.")
        scaler = None

    # Init Basilisk
    ref_MRP = RigidBodyKinematics.euler1232MRP(np.array([0.4, 0.2, -0.3]))
    b_model = BasiliskModel(
    I=[0.025, 0,     0,
       0,     0.05,  0,
       0,     0,     0.065],
    ref_MRP=ref_MRP,
    torque_mode="wheel"
)


    # Buffers, logs
    window = []
    C = Counters()
    recovery = Recovery()
    prev_speeds = None

    times, scores, flags = [], [], []
    modes, trig_notes = [], []
    torque_hist, power_hist = [], []
    batt_hist = []

    for i in range(STEPS):
        # 1) baseline control (swap with your agent here if desired)
        raw_torque = PDControllor(b_model.cur_error_MRP, b_model.cur_omega)
        torque4 = vec4(raw_torque, fill=0.0)           # ensure length 4 for RW0..RW3


        # 2) read telemetry
        rw_speeds = vec4(get_rw_speeds(b_model), fill=0.0)    # rad/s
        rw_power  = vec4(get_rw_powers(b_model),  fill=0.0)   # W
        # W
        batt_storage, batt_net = get_battery_metrics(b_model)

        # 3) build feature window & anomaly score
        feat = build_feature_row(torque4, b_model.cur_MRP, b_model.cur_error_MRP, b_model.cur_omega, rw_power)
        window.append(feat)
        if len(window) > WINDOW:
            window.pop(0)

        anomaly = False
        score = np.nan
        if len(window) == WINDOW:
            X = np.array(window, dtype=np.float32)               # (T, F)
            if scaler is None:
                scaler = MinMaxScaler().fit(X)
            Xn = scaler.transform(X).reshape(1, WINDOW, -1)      # (1, T, F)
            Xp = model.predict(Xn, verbose=0)                    # (1, T, F)
            err_tf = np.mean(np.abs(Xp - Xn), axis=2)            # (1, T)
            score = float(np.mean(err_tf, axis=1)[0])
            anomaly = score > THRESHOLD

        # 4) fault classification (rule-based), ONLY when anomaly or to keep fast dominance faults responsive
        reason = None
        mode_to_trigger, wheels_to_trigger = None, []
        # Battery highest priority
        batt_fault, batt_reason = classify_battery(batt_storage, batt_net, C)
        power_fault, power_reason = classify_powerlimit(rw_power, C)
        enc_wheels, enc_reasons = classify_encoder(torque4, rw_speeds, prev_speeds, C)
        fric_wheels, fric_reasons = classify_friction(torque4, rw_speeds, rw_power, C)

        if anomaly:
            if batt_fault:
                mode_to_trigger, reason = "battery", batt_reason
            elif power_fault:
                mode_to_trigger, reason = "powerlimit", power_reason
            elif len(enc_wheels) > 0:
                mode_to_trigger, wheels_to_trigger = "encoder", enc_wheels
                reason = "; ".join(enc_reasons)
            elif len(fric_wheels) > 0:
                mode_to_trigger, wheels_to_trigger = "friction", fric_wheels
                reason = "; ".join(fric_reasons)
            else:
                mode_to_trigger, reason = "unknown", "anomaly>threshold"

        # 5) trigger recovery once per anomaly episode
        if mode_to_trigger and not recovery.active:
            recovery.trigger(mode_to_trigger, wheels=wheels_to_trigger)
            note = f"t={((i+1)*DT_SEC):.2f}s score={score:.3f} trigger={mode_to_trigger} {reason or ''}"
            print("[Integration] " + note)
            trig_notes.append(note)

        # 6) apply recovery action (if active)
        if recovery.active:
            torque4 = apply_recovery(recovery.mode, torque4, recovery.wheels)

        # 7) step sim
        b_model.step(macros.sec2nano((i+1)*DT_SEC), np.append(torque4, 0.0))


        # bookkeeping
        recovery.tick()
        prev_speeds = rw_speeds.copy()

        times.append((i+1)*DT_SEC/60.0)    # minutes
        scores.append(score)
        flags.append(bool(anomaly))
        modes.append(recovery.mode or "normal")
        torque_hist.append(torque4[:4].copy())
        power_hist.append(rw_power[:4].copy())
        batt_hist.append([batt_storage if batt_storage is not None else np.nan,
                          batt_net if batt_net is not None else np.nan])

    # ---------------------------- Plots ----------------------------
    times = np.array(times)
    torque_hist = np.array(torque_hist)
    power_hist = np.array(power_hist)
    batt_hist = np.array(batt_hist)

    fig, axes = plt.subplots(3, 2, figsize=(13, 12))
    # Score
    ax = axes[0][0]
    ax.plot(times, scores, label="Recon. error")
    ax.axhline(THRESHOLD, linestyle="--", label=f"Threshold {THRESHOLD:.3f}")
    ax.set_title("Anomaly score"); ax.set_xlabel("Time (min)"); ax.legend(); ax.grid(True, alpha=0.3)

    # Recovery timeline
    ax = axes[0][1]
    y = np.array([0 if m=="normal" else 1 for m in modes])
    ax.plot(times*60, y, drawstyle="steps-post", label="Recovery active")
    ax.set_title("Recovery state"); ax.set_xlabel("Time (s)"); ax.set_yticks([0,1],["No","Yes"]); ax.grid(True, alpha=0.3)

    # Battery
    ax = axes[1][0]
    ax.plot(times, batt_hist[:,0], label="Battery storage (Ws)")
    ax.axhline(BATTERY_MIN_LEVEL, linestyle="--", label="Min storage")
    ax.set_title("Battery storage"); ax.set_xlabel("Time (min)"); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1][1]
    ax.plot(times, batt_hist[:,1], label="Battery net power (W)")
    ax.axhline(BATTERY_NET_NEG_THRESH, linestyle="--", label="Neg power thresh")
    ax.set_title("Battery net power"); ax.set_xlabel("Time (min)"); ax.legend(); ax.grid(True, alpha=0.3)

    # Power and torques
    ax = axes[2][0]
    ax.plot(times*60, np.sum(power_hist, axis=1), label="Total RW power")
    ax.axhline(TOTAL_POWER_LIMIT, linestyle="--", label="Total limit")
    ax.set_title("RW power (total)"); ax.set_xlabel("Time (s)"); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2][1]
    ax.plot(times*60, torque_hist[:,0], label="τ0")
    ax.plot(times*60, torque_hist[:,1], label="τ1")
    ax.plot(times*60, torque_hist[:,2], label="τ2")
    ax.plot(times*60, torque_hist[:,3], label="τ3")
    ax.set_title("RW torques"); ax.set_xlabel("Time (s)"); ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
