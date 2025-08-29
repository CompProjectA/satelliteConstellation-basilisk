# eval_ppo.py
from stable_baselines3 import PPO
import numpy as np
from collections import defaultdict

from fault_drl_env import SingleSatFaultEnv, EnvConfig

# 1) Build the SAME env config you want to test
env = SingleSatFaultEnv(EnvConfig(
    single_fault_per_episode=True,          # set False to test 4 faults per episode
    include_no_fault_episodes=True,         # include no-fault episodes if desired
    seed=123
))

# 2) Load your trained model (the one you just saved)
model = PPO.load("models/ppo_fault_detector.zip", env=env)

# 3) Roll out a few episodes and print metrics
EPISODES = 10
conf = defaultdict(int)        # confusion counts (true, pred)
delays = []                    # detection latency (s), onset -> first correct prediction
correct_after_onset = 0
total_after_onset = 0

for ep in range(EPISODES):
    obs, _ = env.reset()
    done = False
    first_correct_time = None
    onset_s = None
    last_pred = None

    while not done:
        # use deterministic policy for evaluation
        action, _ = model.predict(obs, deterministic=True)
        obs, reward, done, truncated, info = env.step(int(action))
        done = done or truncated

        true = info["true_label"]
        t    = info["t_s"]
        if onset_s is None:
            onset_s = info["fault_onset_s"]

        # confusion accounting
        conf[(true, int(action))] += 1

        # latency: first time we get it right after onset
        if true != 0:
            total_after_onset += 1
            if int(action) == true and first_correct_time is None:
                first_correct_time = t
            if int(action) == true:
                correct_after_onset += 1

        last_pred = int(action)

    if onset_s not in (None, float("inf")) and first_correct_time is not None:
        delays.append(first_correct_time - onset_s)

    print(f"Episode {ep}: final true={true} pred={last_pred} "
          f"latency={None if not delays else round(delays[-1],2)}s")

# 4) Report summary
labels = ["no", "battery", "friction", "powerlimit", "encoder"]
print("\nConfusion (true,pred) counts:")
for i in range(5):
    row = [conf[(i, j)] for j in range(5)]
    print(f"{labels[i]:>10}: {row}")

if delays:
    print(f"\nMean detection latency: {np.mean(delays):.2f}s  (n={len(delays)})")
if total_after_onset > 0:
    print(f"Step-wise accuracy after onset: {correct_after_onset/total_after_onset:.3f}")
