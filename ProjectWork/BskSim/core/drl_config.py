#!/usr/bin/env python
"""
drl_config.py

Configuration management for DRL-based task reassignment system.
Manages thresholds, model parameters, and integration settings.

This version normalises all paths to the current BskSim layout:
- ROOT_DIR = .../BskSim
- CORE_DIR = .../BskSim/core
- DRL_DIR  = .../BskSim/DRL
- Default model_path -> .../BskSim/core/anomaly_detection_model.keras
"""

import os
import json
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

# -------- Project paths (single source of truth) --------
# Absolute paths derived from this file location (core/)
CORE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CORE_DIR)                  # .../BskSim
DRL_DIR = os.path.join(ROOT_DIR, "DRL")
MODELS_DIR = os.path.join(ROOT_DIR, "models")
PLOTS_DIR = os.path.join(ROOT_DIR, "plots")
LOGS_DIR = os.path.join(ROOT_DIR, "logs")
RESULT_DIR = os.path.join(DRL_DIR, "result")
CONFIG_DIR = os.path.join(CORE_DIR, "config")
DEFAULT_CONFIG_FILE = os.path.join(CONFIG_DIR, "integrated_system_config.json")

# Known model file you actually have (see your dir dump)
DEFAULT_ANOMALY_MODEL = os.path.join(CORE_DIR, "anomaly_detection_model.keras")


@dataclass
class FaultDetectionConfig:
    """Configuration for fault detection thresholds"""
    ml_confidence_threshold: float = 0.6
    anomaly_score_threshold: float = 50.0
    speed_change_threshold_percent: float = 10.0
    attitude_error_threshold_percent: float = 50.0
    detection_window_minutes: float = 1.0
    minimum_data_points: int = 30


@dataclass
class DRLConfig:
    """Configuration for DRL agent parameters"""
    # Updated to real file in core/
   
    model_path: str = "core/anomaly_detection_model.keras"

    state_dim: int = 10
    action_dim: int = 3
    learning_rate: float = 0.0003
    training_episodes: int = 1000
    max_steps_per_episode: int = 500
    update_frequency: int = 20
    use_pretrained: bool = True


@dataclass
class TaskReassignmentConfig:
    """Configuration for task reassignment strategies"""
    reassignment_strategies: List[str] = None
    max_load_increase_per_spacecraft: float = 0.3
    priority_levels: List[str] = None
    task_redistribution_timeout_seconds: int = 30
    fallback_to_rules: bool = True

    def __post_init__(self):
        if self.reassignment_strategies is None:
            self.reassignment_strategies = [
                "even_distribution",
                "capability_based",
                "load_balanced",
            ]
        if self.priority_levels is None:
            self.priority_levels = ["low", "normal", "high", "critical"]


@dataclass
class IntegrationConfig:
    """Main configuration for fault detection + DRL integration"""
    fault_detection: FaultDetectionConfig = None
    drl: DRLConfig = None
    task_reassignment: TaskReassignmentConfig = None
    simulation_time_limit_minutes: float = 30.0
    enable_real_time_detection: bool = False
    enable_drl_logging: bool = True
    output_directory: str = LOGS_DIR  # normalised

    def __post_init__(self):
        if self.fault_detection is None:
            self.fault_detection = FaultDetectionConfig()
        if self.drl is None:
            self.drl = DRLConfig()
        if self.task_reassignment is None:
            self.task_reassignment = TaskReassignmentConfig()


class ConfigManager:
    """Manages configuration for the integrated system"""

    def __init__(self, config_file: Optional[str] = None):
        self.config_file = config_file or DEFAULT_CONFIG_FILE
        # Ensure config/log/result dirs exist
        for d in (CONFIG_DIR, LOGS_DIR, RESULT_DIR, PLOTS_DIR, MODELS_DIR):
            os.makedirs(d, exist_ok=True)

        self.config = IntegrationConfig()

        # Load existing config if available, else save defaults
        if os.path.exists(self.config_file):
            self.load_config()
        else:
            self.save_config()

    def load_config(self) -> bool:
        """Load configuration from JSON file"""
        try:
            with open(self.config_file, "r") as f:
                config_dict = json.load(f)

            # Parse nested configuration (falling back to defaults)
            fault_config = FaultDetectionConfig(**config_dict.get("fault_detection", {}))
            # Normalise model_path if missing/old
            drl_dict = config_dict.get("drl", {})
            model_path = drl_dict.get("model_path", DEFAULT_ANOMALY_MODEL)
            # If someone saved a stale path, hard-fix to known file if it exists
            if not os.path.isfile(model_path) and os.path.isfile(DEFAULT_ANOMALY_MODEL):
                model_path = DEFAULT_ANOMALY_MODEL
            drl_dict["model_path"] = model_path
            drl_config = DRLConfig(**drl_dict)

            task_config = TaskReassignmentConfig(**config_dict.get("task_reassignment", {}))

            self.config = IntegrationConfig(
                fault_detection=fault_config,
                drl=drl_config,
                task_reassignment=task_config,
                **{
                    k: v
                    for k, v in config_dict.items()
                    if k not in ["fault_detection", "drl", "task_reassignment"]
                },
            )

            print(f"Configuration loaded from {self.config_file}")
            return True

        except Exception as e:
            print(f"Error loading config: {e}")
            print("Using default configuration")
            # Ensure we still save a valid default for next time
            self.config.drl.model_path = DEFAULT_ANOMALY_MODEL
            self.save_config()
            return False

    def save_config(self) -> bool:
        """Save current configuration to JSON file"""
        try:
            config_dict = {
                "fault_detection": asdict(self.config.fault_detection),
                "drl": asdict(self.config.drl),
                "task_reassignment": asdict(self.config.task_reassignment),
                "simulation_time_limit_minutes": self.config.simulation_time_limit_minutes,
                "enable_real_time_detection": self.config.enable_real_time_detection,
                "enable_drl_logging": self.config.enable_drl_logging,
                "output_directory": self.config.output_directory,
            }

            with open(self.config_file, "w") as f:
                json.dump(config_dict, f, indent=2)

            print(f"Configuration saved to {self.config_file}")
            return True

        except Exception as e:
            print(f"Error saving config: {e}")
            return False

    def update_fault_thresholds(self, **kwargs):
        """Update fault detection thresholds"""
        for key, value in kwargs.items():
            if hasattr(self.config.fault_detection, key):
                setattr(self.config.fault_detection, key, value)
                print(f"Updated fault detection {key}: {value}")

    def update_drl_parameters(self, **kwargs):
        """Update DRL parameters"""
        for key, value in kwargs.items():
            if hasattr(self.config.drl, key):
                setattr(self.config.drl, key, value)
                print(f"Updated DRL {key}: {value}")

    def update_task_settings(self, **kwargs):
        """Update task reassignment settings"""
        for key, value in kwargs.items():
            if hasattr(self.config.task_reassignment, key):
                setattr(self.config.task_reassignment, key, value)
                print(f"Updated task reassignment {key}: {value}")

    def get_detection_thresholds(self) -> Dict:
        """Get current fault detection thresholds"""
        return {
            "ml_confidence": self.config.fault_detection.ml_confidence_threshold,
            "anomaly_score": self.config.fault_detection.anomaly_score_threshold,
            "speed_change": self.config.fault_detection.speed_change_threshold_percent,
            "attitude_error": self.config.fault_detection.attitude_error_threshold_percent,
        }

    def get_drl_parameters(self) -> Dict:
        """Get current DRL parameters"""
        return {
            "state_dim": self.config.drl.state_dim,
            "action_dim": self.config.drl.action_dim,
            "learning_rate": self.config.drl.learning_rate,
            "model_path": self.config.drl.model_path,
        }

    def validate_config(self) -> List[str]:
        """Validate configuration and return list of warnings/errors"""
        warnings = []

        # Fault thresholds
        if not (0 <= self.config.fault_detection.ml_confidence_threshold <= 1):
            warnings.append("ML confidence threshold should be between 0 and 1")

        if self.config.fault_detection.anomaly_score_threshold < 0:
            warnings.append("Anomaly score threshold should be positive")

        # DRL params
        if self.config.drl.state_dim <= 0:
            warnings.append("DRL state dimension should be positive")

        if self.config.drl.action_dim <= 0:
            warnings.append("DRL action dimension should be positive")

        # File paths
        model_dir = os.path.dirname(self.config.drl.model_path or "")
        if model_dir and not os.path.isdir(model_dir):
            warnings.append(f"DRL model directory does not exist: {model_dir}")

        if self.config.drl.model_path and not os.path.isfile(self.config.drl.model_path):
            warnings.append(f"DRL model file not found: {self.config.drl.model_path}")

        # Task settings
        if not (0 <= self.config.task_reassignment.max_load_increase_per_spacecraft <= 1):
            warnings.append("Max load increase should be between 0 and 1")

        return warnings

    def print_config_summary(self):
        """Print a summary of the current configuration"""
        print("\nCURRENT CONFIGURATION SUMMARY")
        print("=" * 50)

        print("\nPaths:")
        print(f"  ROOT_DIR : {ROOT_DIR}")
        print(f"  CORE_DIR : {CORE_DIR}")
        print(f"  DRL_DIR  : {DRL_DIR}")
        print(f"  RESULT   : {RESULT_DIR}")
        print(f"  LOGS     : {LOGS_DIR}")

        print("\nFault Detection:")
        print(f"  ML Confidence Threshold: {self.config.fault_detection.ml_confidence_threshold}")
        print(f"  Anomaly Score Threshold: {self.config.fault_detection.anomaly_score_threshold}")
        print(f"  Speed Change Threshold: {self.config.fault_detection.speed_change_threshold_percent}%")
        print(f"  Attitude Error Threshold: {self.config.fault_detection.attitude_error_threshold_percent}%")

        print("\nDRL Agent:")
        print(f"  Model Path: {self.config.drl.model_path}")
        print(f"  State Dimension: {self.config.drl.state_dim}")
        print(f"  Action Dimension: {self.config.drl.action_dim}")
        print(f"  Learning Rate: {self.config.drl.learning_rate}")
        print(f"  Use Pretrained: {self.config.drl.use_pretrained}")

        print("\nTask Reassignment:")
        print(f"  Strategies: {', '.join(self.config.task_reassignment.reassignment_strategies)}")
        print(f"  Max Load Increase: {self.config.task_reassignment.max_load_increase_per_spacecraft}")
        print(f"  Fallback to Rules: {self.config.task_reassignment.fallback_to_rules}")

        print("\nSystem Settings:")
        print(f"  Simulation Time Limit: {self.config.simulation_time_limit_minutes} min")
        print(f"  Real-time Detection: {self.config.enable_real_time_detection}")
        print(f"  DRL Logging: {self.config.enable_drl_logging}")
        print(f"  Output Directory: {self.config.output_directory}")

        warnings = self.validate_config()
        if warnings:
            print("\nCONFIGURATION WARNINGS:")
            for warning in warnings:
                print(f"  ⚠️ {warning}")


# Predefined configuration profiles for different scenarios
class ConfigProfiles:
    """Predefined configuration profiles for different use cases"""

    @staticmethod
    def get_development_config() -> IntegrationConfig:
        """Configuration for development and testing"""
        return IntegrationConfig(
            fault_detection=FaultDetectionConfig(
                ml_confidence_threshold=0.5,
                anomaly_score_threshold=30.0,
                speed_change_threshold_percent=5.0,
                attitude_error_threshold_percent=25.0,
            ),
            drl=DRLConfig(
                model_path=DEFAULT_ANOMALY_MODEL,
                learning_rate=0.001,
                training_episodes=500,
                use_pretrained=False,
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.5,
                fallback_to_rules=True,
            ),
            enable_drl_logging=True,
            enable_real_time_detection=False,
            output_directory=LOGS_DIR,
        )

    @staticmethod
    def get_production_config() -> IntegrationConfig:
        """Configuration for production deployment"""
        return IntegrationConfig(
            fault_detection=FaultDetectionConfig(
                ml_confidence_threshold=0.8,
                anomaly_score_threshold=70.0,
                speed_change_threshold_percent=15.0,
                attitude_error_threshold_percent=75.0,
            ),
            drl=DRLConfig(
                model_path=DEFAULT_ANOMALY_MODEL,
                learning_rate=0.0001,
                training_episodes=2000,
                use_pretrained=True,
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.2,
                fallback_to_rules=True,
            ),
            enable_drl_logging=False,
            enable_real_time_detection=True,
            output_directory=LOGS_DIR,
        )

    @staticmethod
    def get_research_config() -> IntegrationConfig:
        """Configuration for research and experimentation"""
        return IntegrationConfig(
            fault_detection=FaultDetectionConfig(
                ml_confidence_threshold=0.3,
                anomaly_score_threshold=20.0,
                speed_change_threshold_percent=2.0,
                attitude_error_threshold_percent=10.0,
            ),
            drl=DRLConfig(
                model_path=DEFAULT_ANOMALY_MODEL,
                learning_rate=0.0005,
                training_episodes=5000,
                use_pretrained=False,
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.8,
                fallback_to_rules=False,
            ),
            enable_drl_logging=True,
            enable_real_time_detection=True,
            output_directory=LOGS_DIR,
        )


def create_config_manager(profile: str = "development") -> ConfigManager:
    """Create a configuration manager with a specific profile"""
    cfg = ConfigManager()
    if profile == "development":
        cfg.config = ConfigProfiles.get_development_config()
    elif profile == "production":
        cfg.config = ConfigProfiles.get_production_config()
    elif profile == "research":
        cfg.config = ConfigProfiles.get_research_config()
    else:
        print(f"Unknown profile '{profile}', using default configuration")

    # Always normalise known paths on save
    if not cfg.config.drl.model_path or not os.path.isfile(cfg.config.drl.model_path):
        if os.path.isfile(DEFAULT_ANOMALY_MODEL):
            cfg.config.drl.model_path = DEFAULT_ANOMALY_MODEL

    cfg.save_config()
    return cfg


# Handy exports for other modules to import paths directly
__all__ = [
    "ConfigManager",
    "IntegrationConfig",
    "FaultDetectionConfig",
    "DRLConfig",
    "TaskReassignmentConfig",
    "ConfigProfiles",
    "create_config_manager",
    # paths
    "ROOT_DIR",
    "CORE_DIR",
    "DRL_DIR",
    "MODELS_DIR",
    "RESULT_DIR",
    "LOGS_DIR",
    "PLOTS_DIR",
    "DEFAULT_ANOMALY_MODEL",
]

if __name__ == "__main__":
    print("DRL Configuration Manager (normalised paths)")
    cfg = create_config_manager("development")
    cfg.print_config_summary()
