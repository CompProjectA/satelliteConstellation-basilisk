import numpy as np
import matplotlib.pyplot as plt

# Stub: replace with actual simulator interaction
def run_simulation_with_fault(rw_fault):
    np.random.seed(hash(rw_fault) % 1000)
    x_requests = np.linspace(0, 10000, 50)
    x_orbits = np.linspace(0, 3.0, 50)

    return {
        "requests": x_requests,
        "imaged_success": np.clip(30 + np.random.randn(len(x_requests)) * 5, 0, 40),
        "reward": np.clip(20 + np.random.randn(len(x_requests)) * 3, 0, 30),
        "battery_ratio": np.clip(0.45 + np.random.randn(len(x_requests)) * 0.03, 0.3, 0.6),
        "orbits": x_orbits,
        "alive_percent": np.clip(90 - np.linspace(0, 40, len(x_orbits)) + np.random.randn(len(x_orbits)) * 3, 0, 100)
    }

# Run simulations for each RW fault
faults = ["RW1", "RW2", "RW3", "RW4"]
results = {rw: run_simulation_with_fault(rw) for rw in faults}

# Plotting helper
def plot_metric(title, ylabel, xkey, ykey, ylim, xlim, xlabel, filename):
    colors = {
        "RW1": "tab:blue",
        "RW2": "tab:orange",
        "RW3": "tab:green",
        "RW4": "tab:red"
    }
    plt.figure(figsize=(10, 6))
    for rw in faults:
        plt.plot(results[rw][xkey], results[rw][ykey], label=f"{rw} Fault", color=colors[rw])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(*ylim)
    plt.xlim(*xlim)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

# ✅ Graph 1: Imaged Success / Non-fault (%)
plot_metric(
    title="Graph 1: Imaged Success / Non-fault (%)",
    ylabel="Imaged Success (%)",
    xkey="requests",
    ykey="imaged_success",
    ylim=(0, 40),
    xlim=(0, 10000),
    xlabel="Requests |R|",
    filename="graph1_imaged_success.png"
)

# ✅ Graph 2: Reward / Non-fault (%)
plot_metric(
    title="Graph 2: Reward / Non-fault (%)",
    ylabel="Reward",
    xkey="requests",
    ykey="reward",
    ylim=(0, 30),
    xlim=(0, 10000),
    xlabel="Requests |R|",
    filename="graph2_reward.png"
)

# ✅ Graph 3: Average Battery Ratio
plot_metric(
    title="Graph 3: Average Battery Ratio",
    ylabel="Battery Ratio",
    xkey="requests",
    ykey="battery_ratio",
    ylim=(0.3, 0.6),
    xlim=(0, 10000),
    xlabel="Requests |R|",
    filename="graph3_battery_ratio.png"
)

# ✅ Graph 4: Alive / Non-fault (%)
plot_metric(
    title="Graph 4: Alive / Non-fault (%)",
    ylabel="Alive (%)",
    xkey="orbits",
    ykey="alive_percent",
    ylim=(0, 100),
    xlim=(0, 3.0),
    xlabel="Orbits Completed",
    filename="graph4_alive_percent.png"
)
