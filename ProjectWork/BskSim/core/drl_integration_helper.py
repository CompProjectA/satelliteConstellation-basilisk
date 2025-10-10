#!/usr/bin/env python
"""
drl_integration_helper.py

Helper functions for integrating DRL with the GUI and backend systems.
Place this file in: core/drl_integration_helper.py
"""

import os
import sys
import json
import numpy as np
from typing import Dict, List, Optional, Any
from datetime import datetime

# Add paths
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
drl_dir = os.path.join(parent_dir, "DRL")

for path in [parent_dir, drl_dir]:
    if path not in sys.path:
        sys.path.insert(0, path)


class DRLIntegrationHelper:
    """
    Helper class to bridge GUI, DRL algorithms, and spacecraft simulation
    """
    
    def __init__(self, parent_app=None):
        self.parent_app = parent_app
        self.drl_algorithms = {
            "PPO": "PPO_Year2.py",
            "TDHD": "TDHDYear2.py",
            "DQN": "DQNYear2.py"
        }
    
    def get_available_drl_scripts(self) -> Dict[str, str]:
        """Get available DRL training scripts"""
        scripts = {}
        for name, filename in self.drl_algorithms.items():
            path = os.path.join(drl_dir, filename)
            if os.path.exists(path):
                scripts[name] = path
        return scripts
    
    def get_constellation_config(self) -> Dict[str, Any]:
        """Extract constellation configuration from GUI"""
        if not self.parent_app:
            return self._get_default_config()
        
        config = {
            "satellites": [],
            "clusters": [],
            "targets": [],
            "faults": [],
            "simulation_time": 30.0,
            "timestamp": datetime.now().isoformat()
        }
        
        # Extract satellites
        if hasattr(self.parent_app, 'satellites'):
            for sat in self.parent_app.satellites:
                config["satellites"].append({
                    "name": sat.get("name", "Unknown"),
                    "altitude": sat.get("altitude", 600),
                    "inclination": sat.get("inclination", 97.0),
                    "raan": sat.get("raan", 0.0),
                    "true_anomaly": sat.get("true_anomaly", 0.0),
                    "cluster_id": sat.get("cluster_id", None),
                    "is_leader": sat.get("is_leader", False)
                })
        
        # Extract clusters
        if hasattr(self.parent_app, 'constellation_tab'):
            clusters = getattr(self.parent_app.constellation_tab, 'clusters', {})
            for cluster_id, cluster_data in clusters.items():
                config["clusters"].append({
                    "id": cluster_id,
                    "name": cluster_data.get("name", f"Cluster {cluster_id}"),
                    "formation": cluster_data.get("formation", "Leader-Follower"),
                    "num_satellites": cluster_data.get("num_satellites", 4)
                })
        
        # Extract targets
        if hasattr(self.parent_app, 'targets'):
            for target in self.parent_app.targets:
                config["targets"].append({
                    "name": target.get("name", "Unknown"),
                    "latitude": target.get("latitude", 0.0),
                    "longitude": target.get("longitude", 0.0),
                    "priority": target.get("priority", 1.0),
                    "assigned_satellites": target.get("assigned_satellites", [])
                })
        
        # Extract faults
        if hasattr(self.parent_app, 'fault_tab'):
            fault_tab = self.parent_app.fault_tab
            if hasattr(fault_tab, 'get_configured_faults'):
                config["faults"] = fault_tab.get_configured_faults()
        
        # Get simulation time
        if hasattr(self.parent_app, 'simulation_time'):
            config["simulation_time"] = self.parent_app.simulation_time.get()
        
        return config
    
    def _get_default_config(self) -> Dict[str, Any]:
        """Get default configuration for testing"""
        return {
            "satellites": [
                {"name": f"Sat-{i}", "altitude": 600, "inclination": 97.0, 
                 "raan": i * 90, "true_anomaly": 0.0}
                for i in range(4)
            ],
            "clusters": [
                {"id": 0, "name": "Cluster 0", "formation": "Leader-Follower", "num_satellites": 4}
            ],
            "targets": [
                {"name": f"Target-{i}", "latitude": -30.0 + i * 2, "longitude": 130.0 + i * 2, 
                 "priority": 1.0}
                for i in range(10)
            ],
            "faults": [],
            "simulation_time": 30.0,
            "timestamp": datetime.now().isoformat()
        }
    
    def prepare_drl_environment_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare configuration for DRL environment"""
        return {
            "n_satellites": len(config["satellites"]),
            "n_targets": len(config["targets"]),
            "simulation_duration": config["simulation_time"],
            "satellite_params": config["satellites"],
            "target_params": config["targets"],
            "fault_config": config.get("faults", [])
        }
    
    def extract_fault_detection_results(self) -> Dict[str, Any]:
        """Extract fault detection results from GUI"""
        if not self.parent_app or not hasattr(self.parent_app, 'ml_detection_results'):
            return {"detections": {}, "summary": {"total_detections": 0}}
        
        results = self.parent_app.ml_detection_results
        if not results:
            # Try to get from fault detection tab
            if hasattr(self.parent_app, 'fault_detection_tab'):
                fd_tab = self.parent_app.fault_detection_tab
                if hasattr(fd_tab, 'detection_results'):
                    results = fd_tab.detection_results
        
        return results or {"detections": {}, "summary": {"total_detections": 0}}
    
    def convert_faults_to_capabilities_lost(self, faults: List[Dict]) -> Dict[str, List[str]]:
        """Convert fault information to lost capabilities per satellite"""
        capabilities_lost = {}
        
        capability_mapping = {
            "friction": ["attitude_control", "pointing_accuracy"],
            "power": ["power_generation", "system_operations"],
            "encoder": ["attitude_control", "sensing"],
            "battery": ["power_generation", "power_storage"],
            "communication": ["communication", "data_relay"],
            "sensor": ["sensing", "navigation"]
        }
        
        for fault in faults:
            sat_name = fault.get("satellite", "Unknown")
            fault_type = fault.get("type", "").lower()
            
            if sat_name not in capabilities_lost:
                capabilities_lost[sat_name] = []
            
            # Map fault type to capabilities
            for key, caps in capability_mapping.items():
                if key in fault_type:
                    capabilities_lost[sat_name].extend(caps)
                    break
            else:
                # Default: any fault affects attitude control
                capabilities_lost[sat_name].append("attitude_control")
        
        # Remove duplicates
        for sat_name in capabilities_lost:
            capabilities_lost[sat_name] = list(set(capabilities_lost[sat_name]))
        
        return capabilities_lost
    
    def create_drl_state_vector(self, config: Dict, fault_results: Dict) -> np.ndarray:
        """Create state vector for DRL agent"""
        n_sats = len(config["satellites"])
        n_faults = len(fault_results.get("detections", {}))
        
        healthy_ratio = (n_sats - n_faults) / n_sats if n_sats > 0 else 0.0
        fault_ratio = n_faults / n_sats if n_sats > 0 else 0.0
        
        # Calculate average anomaly score
        detections = fault_results.get("detections", {})
        anomaly_scores = [d.get("anomaly_score", 0.0) for d in detections.values()]
        avg_anomaly = np.mean(anomaly_scores) if anomaly_scores else 0.0
        max_anomaly = np.max(anomaly_scores) if anomaly_scores else 0.0
        
        # Task load metrics (placeholder)
        avg_task_load = 0.5
        max_task_load = 0.8
        
        state = np.array([
            healthy_ratio,      # 0: ratio of healthy spacecraft
            fault_ratio,        # 1: ratio of faulty spacecraft  
            n_sats / 10.0,      # 2: normalized number of satellites
            n_faults / 10.0,    # 3: normalized number of faults
            avg_anomaly,        # 4: average anomaly score
            max_anomaly,        # 5: maximum anomaly score
            avg_task_load,      # 6: average task load
            max_task_load       # 7: maximum task load
        ], dtype=np.float32)
        
        return state
    
    def interpret_drl_action(self, action: int, state: np.ndarray) -> Dict[str, Any]:
        """Interpret DRL action and convert to reassignment strategy"""
        strategies = {
            0: "even_distribution",
            1: "capability_based",
            2: "load_balanced",
            3: "priority_based"
        }
        
        strategy = strategies.get(action, "even_distribution")
        
        return {
            "action": action,
            "strategy": strategy,
            "state_vector": state.tolist(),
            "recommendations": self._generate_recommendations(strategy, state)
        }
    
    def _generate_recommendations(self, strategy: str, state: np.ndarray) -> List[str]:
        """Generate human-readable recommendations"""
        recommendations = []
        
        healthy_ratio = state[0]
        fault_ratio = state[1]
        avg_anomaly = state[4]
        
        if healthy_ratio < 0.3:
            recommendations.append("CRITICAL: Less than 30% of satellites are healthy")
            recommendations.append("Consider activating backup systems")
        
        if fault_ratio > 0.5:
            recommendations.append("WARNING: More than 50% of satellites have faults")
        
        if avg_anomaly > 0.7:
            recommendations.append("HIGH SEVERITY: Average anomaly score is high")
        
        if strategy == "even_distribution":
            recommendations.append("Distributing tasks evenly across all healthy spacecraft")
        elif strategy == "capability_based":
            recommendations.append("Prioritizing most capable spacecraft for critical tasks")
        elif strategy == "load_balanced":
            recommendations.append("Balancing workload based on current capacity")
        elif strategy == "priority_based":
            recommendations.append("Assigning tasks based on target priority")
        
        return recommendations
    
    def save_integration_results(self, results: Dict[str, Any], output_dir: str):
        """Save integration results to file"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save JSON
        json_path = os.path.join(output_dir, 
                                f"drl_integration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        try:
            with open(json_path, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            print(f"Integration results saved: {json_path}")
        except Exception as e:
            print(f"Error saving JSON: {e}")
        
        # Save text report
        report_path = os.path.join(output_dir,
                                  f"drl_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
        try:
            with open(report_path, 'w') as f:
                self._write_text_report(f, results)
            print(f"Text report saved: {report_path}")
        except Exception as e:
            print(f"Error saving report: {e}")
        
        return json_path, report_path
    
    def _write_text_report(self, file, results: Dict[str, Any]):
        """Write human-readable text report"""
        file.write("=" * 70 + "\n")
        file.write("DRL INTEGRATION REPORT\n")
        file.write("=" * 70 + "\n\n")
        file.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Configuration
        if "configuration" in results:
            config = results["configuration"]
            file.write("CONFIGURATION:\n")
            file.write(f"  Satellites: {len(config.get('satellites', []))}\n")
            file.write(f"  Targets: {len(config.get('targets', []))}\n")
            file.write(f"  Simulation Time: {config.get('simulation_time', 0)} min\n\n")
        
        # Fault Detection
        if "fault_detection" in results:
            fd = results["fault_detection"]
            file.write("FAULT DETECTION:\n")
            summary = fd.get("summary", {})
            file.write(f"  Total Detections: {summary.get('total_detections', 0)}\n")
            file.write(f"  Success Rate: {summary.get('success_rate', 0):.1%}\n\n")
        
        # DRL Decisions
        if "drl_decisions" in results:
            drl = results["drl_decisions"]
            file.write("DRL DECISIONS:\n")
            file.write(f"  Strategy: {drl.get('strategy', 'N/A')}\n")
            file.write(f"  Action: {drl.get('action', 'N/A')}\n")
            file.write(f"  Confidence: {drl.get('decision_confidence', 0):.2%}\n\n")
            
            recs = drl.get("recommendations", [])
            if recs:
                file.write("  Recommendations:\n")
                for rec in recs:
                    file.write(f"    - {rec}\n")
                file.write("\n")
        
        # Task Reassignment
        if "task_reassignment" in results:
            tr = results["task_reassignment"]
            file.write("TASK REASSIGNMENT:\n")
            file.write(f"  Execution Success: {tr.get('execution_success', False)}\n")
            file.write(f"  Healthy Spacecraft: {len(tr.get('healthy_spacecraft', []))}\n")
            file.write(f"  Faulty Spacecraft: {len(tr.get('faulty_spacecraft', []))}\n\n")


def create_drl_helper(parent_app=None) -> DRLIntegrationHelper:
    """Factory function to create DRL integration helper"""
    return DRLIntegrationHelper(parent_app)


# Quick test
if __name__ == "__main__":
    print("DRL Integration Helper")
    print("=" * 50)
    
    helper = create_drl_helper()
    
    print("\nAvailable DRL Scripts:")
    scripts = helper.get_available_drl_scripts()
    for name, path in scripts.items():
        exists = "✓" if os.path.exists(path) else "✗"
        print(f"  {exists} {name}: {path}")
    
    print("\nDefault Configuration:")
    config = helper._get_default_config()
    print(f"  Satellites: {len(config['satellites'])}")
    print(f"  Targets: {len(config['targets'])}")
    print(f"  Simulation Time: {config['simulation_time']} min")
    
    print("\nDRL Environment Config:")
    env_config = helper.prepare_drl_environment_config(config)
    print(f"  N Satellites: {env_config['n_satellites']}")
    print(f"  N Targets: {env_config['n_targets']}")
    
    print("\nTest State Vector:")
    fault_results = {"detections": {"Sat-0": {"anomaly_score": 0.8}}}
    state = helper.create_drl_state_vector(config, fault_results)
    print(f"  State Shape: {state.shape}")
    print(f"  State Values: {state}")
    
    print("\nTest Action Interpretation:")
    decision = helper.interpret_drl_action(1, state)
    print(f"  Strategy: {decision['strategy']}")
    print(f"  Recommendations:")
    for rec in decision['recommendations']:
        print(f"    - {rec}")