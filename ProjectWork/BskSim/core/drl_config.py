#!/usr/bin/env python
"""
drl_config.py

Configuration management for DRL-based task reassignment system.
Manages thresholds, model parameters, and integration settings.
"""

import os
import json
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict


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
    model_path: str = "core/Agent-based-Architecture-for-Proactive-Fault-Tolerance-and-Management-in-Small-Satellite-Missions/anomaly_detection_model.keras"
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
                "load_balanced"
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
    output_directory: str = "logs"
    
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
        self.config_file = config_file or "config/integrated_system_config.json"
        self.config = IntegrationConfig()
        
        # Create config directory if it doesn't exist
        os.makedirs(os.path.dirname(self.config_file), exist_ok=True)
        
        # Load existing config if available
        if os.path.exists(self.config_file):
            self.load_config()
        else:
            self.save_config()  # Save default config
    
    def load_config(self) -> bool:
        """Load configuration from JSON file"""
        try:
            with open(self.config_file, 'r') as f:
                config_dict = json.load(f)
            
            # Parse nested configuration
            fault_config = FaultDetectionConfig(**config_dict.get("fault_detection", {}))
            drl_config = DRLConfig(**config_dict.get("drl", {}))
            task_config = TaskReassignmentConfig(**config_dict.get("task_reassignment", {}))
            
            # Update main config
            self.config = IntegrationConfig(
                fault_detection=fault_config,
                drl=drl_config,
                task_reassignment=task_config,
                **{k: v for k, v in config_dict.items() 
                   if k not in ["fault_detection", "drl", "task_reassignment"]}
            )
            
            print(f"Configuration loaded from {self.config_file}")
            return True
            
        except Exception as e:
            print(f"Error loading config: {e}")
            print("Using default configuration")
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
                "output_directory": self.config.output_directory
            }
            
            with open(self.config_file, 'w') as f:
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
            "attitude_error": self.config.fault_detection.attitude_error_threshold_percent
        }
    
    def get_drl_parameters(self) -> Dict:
        """Get current DRL parameters"""
        return {
            "state_dim": self.config.drl.state_dim,
            "action_dim": self.config.drl.action_dim,
            "learning_rate": self.config.drl.learning_rate,
            "model_path": self.config.drl.model_path
        }
    
    def validate_config(self) -> List[str]:
        """Validate configuration and return list of warnings/errors"""
        warnings = []
        
        # Validate fault detection thresholds
        if self.config.fault_detection.ml_confidence_threshold < 0 or self.config.fault_detection.ml_confidence_threshold > 1:
            warnings.append("ML confidence threshold should be between 0 and 1")
        
        if self.config.fault_detection.anomaly_score_threshold < 0:
            warnings.append("Anomaly score threshold should be positive")
        
        # Validate DRL parameters
        if self.config.drl.state_dim <= 0:
            warnings.append("DRL state dimension should be positive")
        
        if self.config.drl.action_dim <= 0:
            warnings.append("DRL action dimension should be positive")
        
        # Validate file paths
        if not os.path.exists(os.path.dirname(self.config.drl.model_path)):
            warnings.append(f"DRL model directory does not exist: {os.path.dirname(self.config.drl.model_path)}")
        
        # Validate task reassignment settings
        if self.config.task_reassignment.max_load_increase_per_spacecraft < 0 or self.config.task_reassignment.max_load_increase_per_spacecraft > 1:
            warnings.append("Max load increase should be between 0 and 1")
        
        return warnings
    
    def print_config_summary(self):
        """Print a summary of the current configuration"""
        print("\nCURRENT CONFIGURATION SUMMARY")
        print("=" * 50)
        
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
        
        # Check for configuration warnings
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
                ml_confidence_threshold=0.5,  # Lower threshold for testing
                anomaly_score_threshold=30.0,
                speed_change_threshold_percent=5.0,
                attitude_error_threshold_percent=25.0
            ),
            drl=DRLConfig(
                learning_rate=0.001,  # Higher learning rate for faster training
                training_episodes=500,
                use_pretrained=False
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.5,  # More aggressive reassignment
                fallback_to_rules=True
            ),
            enable_drl_logging=True,
            enable_real_time_detection=False
        )
    
    @staticmethod
    def get_production_config() -> IntegrationConfig:
        """Configuration for production deployment"""
        return IntegrationConfig(
            fault_detection=FaultDetectionConfig(
                ml_confidence_threshold=0.8,  # Higher threshold for reliability
                anomaly_score_threshold=70.0,
                speed_change_threshold_percent=15.0,
                attitude_error_threshold_percent=75.0
            ),
            drl=DRLConfig(
                learning_rate=0.0001,  # Lower learning rate for stability
                training_episodes=2000,
                use_pretrained=True
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.2,  # Conservative reassignment
                fallback_to_rules=True
            ),
            enable_drl_logging=False,  # Disable for performance
            enable_real_time_detection=True
        )
    
    @staticmethod
    def get_research_config() -> IntegrationConfig:
        """Configuration for research and experimentation"""
        return IntegrationConfig(
            fault_detection=FaultDetectionConfig(
                ml_confidence_threshold=0.3,  # Very sensitive for research
                anomaly_score_threshold=20.0,
                speed_change_threshold_percent=2.0,
                attitude_error_threshold_percent=10.0
            ),
            drl=DRLConfig(
                learning_rate=0.0005,
                training_episodes=5000,  # Extended training
                use_pretrained=False
            ),
            task_reassignment=TaskReassignmentConfig(
                max_load_increase_per_spacecraft=0.8,  # Experimental limits
                fallback_to_rules=False  # Force DRL decisions
            ),
            enable_drl_logging=True,
            enable_real_time_detection=True
        )


def create_config_manager(profile: str = "development") -> ConfigManager:
    """Create a configuration manager with a specific profile"""
    
    config_manager = ConfigManager()
    
    if profile == "development":
        config_manager.config = ConfigProfiles.get_development_config()
    elif profile == "production":
        config_manager.config = ConfigProfiles.get_production_config()
    elif profile == "research":
        config_manager.config = ConfigProfiles.get_research_config()
    else:
        print(f"Unknown profile '{profile}', using default configuration")
    
    # Save the profile-based configuration
    config_manager.save_config()
    
    return config_manager


if __name__ == "__main__":
    print("DRL Configuration Manager")
    print("Creating development profile configuration...")
    
    config_manager = create_config_manager("development")
    config_manager.print_config_summary()
    
    # Example of updating configuration
    print("\nUpdating fault detection thresholds...")
    config_manager.update_fault_thresholds(
        ml_confidence_threshold=0.7,
        anomaly_score_threshold=60.0
    )
    
    # Validate configuration
    warnings = config_manager.validate_config()
    if warnings:
        print("\nConfiguration warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("\nConfiguration is valid!")
    
    # Save updated configuration
    config_manager.save_config()