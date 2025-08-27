from Basilisk.architecture import messaging
import numpy as np
from Basilisk import __path__
from Basilisk.simulation import (
    reactionWheelStateEffector,
    simpleBattery,
    spacecraft,
    extForceTorque,
    ReactionWheelPower,
)
from Basilisk.utilities import (
    RigidBodyKinematics,
    unitTestSupport,
    macros,
    simIncludeRW,
    SimulationBaseClass,
    vizSupport,
)
import gymnasium as gym
from gymnasium import spaces
import os
import pandas as pd

bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])


class FrictionFaultModel:
    """Interface to the Basilisk friction fault simulation."""

    def __init__(self) -> None:
        self.rw_speed_in_msg = messaging.RWSpeedMsg()
        self.reset()

    def reset(self) -> None:
        self.cur_error_MRP = np.zeros(3, dtype=np.float32)
        self.cur_omega = np.zeros(3, dtype=np.float32)
        self.cur_rwSpeeds = np.zeros(4, dtype=np.float32)

    def step(self, torque_cmd: np.ndarray) -> np.ndarray:
        """Propagate the model by one step and return an observation vector."""

        # Read current wheel speeds from the message and keep only the first four
        rw_speed_msg = self.rw_speed_in_msg.read()
        self.cur_rwSpeeds = np.array(rw_speed_msg.wheelSpeeds[:4], dtype=np.float32)

        # Observation contains error MRP, angular rates and wheel speeds
        obs = np.hstack((self.cur_error_MRP, self.cur_omega, self.cur_rwSpeeds)).astype(
            np.float32
        )
        return obs


class FrictionFaultEnv(gym.Env):
    """Gymnasium environment exposing the friction fault model."""

    def __init__(self):
        super().__init__()
        self.model = FrictionFaultModel()

        # Observation is 10 elements: 3 error MRP, 3 body rates, 4 wheel speeds
        self.observation_space = spaces.Box(
            low=-np.inf, high=np.inf, shape=(10,), dtype=np.float32
        )
        # Action is 4 commanded wheel torques in the range [-1, 1]
        self.action_space = spaces.Box(low=-1.0, high=1.0, shape=(4,), dtype=np.float32)

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        self.model.reset()
        obs = np.hstack(
            (self.model.cur_error_MRP, self.model.cur_omega, self.model.cur_rwSpeeds)
        ).astype(np.float32)
        return obs, {}

    def step(self, action: np.ndarray):
        obs = self.model.step(action)
        # Placeholder reward: encourage small attitude error
        reward = -float(np.linalg.norm(self.model.cur_error_MRP))
        terminated = False
        truncated = False
        return obs, reward, terminated, truncated, {}

    def close(self) -> None:
        pass


if __name__ == "__main__":
    env = FrictionFaultEnv()
    obs, _ = env.reset()
    for _ in range(3):
        action = env.action_space.sample()
        obs, reward, terminated, truncated, _ = env.step(action)
        if terminated or truncated:
            break
    env.close()
