#!/usr/bin/env python
"""
drl_integration_bridge.py

Bridges your existing fault detection with the DRL task reassignment system.
This file connects real_ml_fault_detection.py with the DRL Integration.py system.
"""

import os
import sys
import numpy as np
from typing import Dict, List, Optional
import json
from datetime import datetime

# Add DRL folder to path
drl_path = os.path.join(os.path.dirname(__file__), 'bsk_rl')
if drl_path not in sys.path:
    sys.path.insert(0, drl_path)

# Import your existing fault detection
from real_ml_fault_detection import RealMLFaultDetector, RealDetectionResult

# Import DRL components (these will be from the provided files)
try:
    from Integration import main as drl_integration_main
    from PPO import PPOAgent
    from Envs import SatelliteEnv
    DRL_AVAILABLE = True
    print("DRL components loaded successfully")
except ImportError as e:
    DRL_AVAILABLE = False
    print(f"DRL components not available: {e}")


class DRLTaskReassignmentSystem:
    """
    Manages task reassignment using DRL after fault detection
    """
    
    def __init__(self, constellation_config):
        self.constellation_config = constellation_config
        self.healthy_spacecraft = []
        self.faulty_spacecraft = []
        self.task_assignments = {}
        self.drl_agent = None
        self.environment = None
        
        print("DRL Task Reassignment System initialized")
        
        if DRL_AVAILABLE:
            self.initialize_drl_components()
    
    def initialize_drl_components(self):
        """Initialize the DRL agent and environment"""
        try:
            # Initialize DRL environment for constellation management
            self.environment = SatelliteEnv()
            
            # Initialize PPO agent for decision making
            self.drl_agent = PPOAgent(
                state_dim=self.environment.observation_space.shape[0],
                action_dim=self.environment.action_space.n
            )
            
            # Load pre-trained model if available
            model_path = os.path.join(drl_path, 'trained_models', 'ppo_constellation.pth')
            if os.path.exists(model_path):
                self.drl_agent.load_model(model_path)
                print(f"Loaded pre-trained DRL model: {model_path}")
            
            print("DRL components initialized successfully")
            
        except Exception as e:
            print(f"Error initializing DRL components: {e}")
            self.drl_agent = None
            self.environment = None
    
    def update_spacecraft_status(self, fault_detection_results: Dict):
        """Update spacecraft health status based on fault detection"""
        
        self.healthy_spacecraft = []
        self.faulty_spacecraft = []
        
        print("\nUpdating spacecraft status from fault detection...")
        
        for spacecraft_name, detections in fault_detection_results.items():
            if detections and len(detections) > 0:
                # Spacecraft has detected faults
                self.faulty_spacecraft.append({
                    'name': spacecraft_name,
                    'faults': detections,
                    'severity': self._assess_fault_severity(detections),
                    'capabilities_lost': self._determine_lost_capabilities(detections)
                })
                print(f"   {spacecraft_name}: FAULTY - {len(detections)} faults detected")
            else:
                # Spacecraft is healthy
                self.healthy_spacecraft.append({
                    'name': spacecraft_name,
                    'capabilities': self._get_spacecraft_capabilities(spacecraft_name),
                    'load_capacity': self._get_load_capacity(spacecraft_name)
                })
                print(f"   {spacecraft_name}: HEALTHY")
        
        print(f"Status update: {len(self.healthy_spacecraft)} healthy, {len(self.faulty_spacecraft)} faulty")
    
    def _assess_fault_severity(self, detections: List[RealDetectionResult]) -> str:
        """Assess the severity of detected faults"""
        if not detections:
            return "none"
        
        max_confidence = max(d.confidence for d in detections)
        max_anomaly_score = max(d.anomaly_score for d in detections)
        
        if max_confidence > 0.8 or max_anomaly_score > 80:
            return "critical"
        elif max_confidence > 0.6 or max_anomaly_score > 50:
            return "moderate"
        else:
            return "minor"
    
    def _determine_lost_capabilities(self, detections: List[RealDetectionResult]) -> List[str]:
        """Determine what capabilities are lost due to faults"""
        lost_capabilities = []
        
        for detection in detections:
            fault_type = detection.fault_type.lower()
            
            # Check for any RW/attitude-related faults
            if any(keyword in fault_type for keyword in ['rw', 'wheel', 'reaction', 'speed', 'friction', 'torque']):
                lost_capabilities.extend(["attitude_control", "pointing_accuracy"])
                print(f"      Mapped '{detection.fault_type}' to attitude_control loss")
            
            elif "power" in fault_type:
                lost_capabilities.extend(["power_generation", "system_operations"])
                print(f"      Mapped '{detection.fault_type}' to power loss")
            
            elif "sensor" in fault_type or "sensing" in fault_type:
                lost_capabilities.extend(["sensing", "navigation"])
                print(f"      Mapped '{detection.fault_type}' to sensing loss")
            
            elif "communication" in fault_type or "comm" in fault_type:
                lost_capabilities.extend(["communication", "data_relay"])
                print(f"      Mapped '{detection.fault_type}' to communication loss")
            
            else:
                # Default: any detected fault causes some capability loss
                lost_capabilities.extend(["attitude_control"])
                print(f"      Mapped unknown fault '{detection.fault_type}' to default attitude_control loss")
        
        return list(set(lost_capabilities))  # Remove duplicates
    
    def _get_spacecraft_capabilities(self, spacecraft_name: str) -> List[str]:
        """Get the capabilities of a healthy spacecraft"""
        # Standard capabilities for a healthy spacecraft
        return [
            "attitude_control",
            "pointing_accuracy", 
            "power_generation",
            "system_operations",
            "sensing",
            "navigation",
            "communication",
            "data_relay"
        ]
    
    def _get_load_capacity(self, spacecraft_name: str) -> float:
        """Get the current load capacity of a spacecraft (0.0 to 1.0)"""
        # For now, assume healthy spacecraft have 70% available capacity
        return 0.7
    
    def trigger_drl_reassignment(self) -> Dict:
        """Trigger DRL-based task reassignment"""
        
        if not DRL_AVAILABLE or self.drl_agent is None:
            print("DRL not available - using fallback reassignment")
            return self._fallback_reassignment()
        
        print("\nTriggering DRL-based task reassignment...")
        
        try:
            # Prepare state vector for DRL agent
            state = self._prepare_drl_state()
            
            # Get DRL agent decision
            action = self.drl_agent.select_action(state, training=False)
            
            # Convert action to task reassignment
            reassignment_plan = self._convert_action_to_reassignment(action)
            
            print(f"DRL agent selected action: {action}")
            print(f"Reassignment plan generated: {len(reassignment_plan)} assignments")
            
            return reassignment_plan
            
        except Exception as e:
            print(f"Error in DRL reassignment: {e}")
            return self._fallback_reassignment()
    
    def _prepare_drl_state(self) -> np.ndarray:
        """Prepare state vector for DRL agent"""
        
        # State includes:
        # - Number of healthy/faulty spacecraft
        # - Severity of faults
        # - Available capabilities
        # - Current task loads
        
        state_vector = []
        
        # Basic constellation state
        total_spacecraft = len(self.healthy_spacecraft) + len(self.faulty_spacecraft)
        state_vector.extend([
            len(self.healthy_spacecraft) / total_spacecraft,  # Healthy ratio
            len(self.faulty_spacecraft) / total_spacecraft,   # Faulty ratio
        ])
        
        # Fault severity metrics
        critical_faults = sum(1 for sc in self.faulty_spacecraft if sc['severity'] == 'critical')
        moderate_faults = sum(1 for sc in self.faulty_spacecraft if sc['severity'] == 'moderate')
        minor_faults = sum(1 for sc in self.faulty_spacecraft if sc['severity'] == 'minor')
        
        state_vector.extend([
            critical_faults / total_spacecraft,
            moderate_faults / total_spacecraft,
            minor_faults / total_spacecraft
        ])
        
        # Available capacity
        total_capacity = sum(sc['load_capacity'] for sc in self.healthy_spacecraft)
        state_vector.append(total_capacity / len(self.healthy_spacecraft) if self.healthy_spacecraft else 0)
        
        # Pad or truncate to expected DRL input size
        expected_size = 10  # Adjust based on your DRL model
        if len(state_vector) < expected_size:
            state_vector.extend([0.0] * (expected_size - len(state_vector)))
        elif len(state_vector) > expected_size:
            state_vector = state_vector[:expected_size]
        
        return np.array(state_vector, dtype=np.float32)
    
    def _convert_action_to_reassignment(self, action) -> Dict:
        """Convert DRL action to task reassignment plan"""
        
        reassignment_plan = {
            "timestamp": datetime.now().isoformat(),
            "trigger": "drl_decision",
            "action_value": action,
            "assignments": [],
            "strategy": "unknown"
        }
        
        # Simple action mapping (extend based on your DRL action space)
        if action == 0:
            # Distribute tasks evenly among healthy spacecraft
            strategy = "even_distribution"
            reassignment_plan["strategy"] = strategy
            
            for i, healthy_sc in enumerate(self.healthy_spacecraft):
                assignment = {
                    "spacecraft": healthy_sc["name"],
                    "new_tasks": self._get_redistributed_tasks(strategy, i),
                    "priority": "normal",
                    "load_increase": 0.2
                }
                reassignment_plan["assignments"].append(assignment)
        
        elif action == 1:
            # Prioritize most capable spacecraft
            strategy = "capability_based"
            reassignment_plan["strategy"] = strategy
            
            # Sort by capability and assign critical tasks to best spacecraft
            sorted_healthy = sorted(self.healthy_spacecraft, 
                                  key=lambda x: x['load_capacity'], reverse=True)
            
            for i, healthy_sc in enumerate(sorted_healthy):
                priority = "high" if i == 0 else "normal"
                assignment = {
                    "spacecraft": healthy_sc["name"],
                    "new_tasks": self._get_redistributed_tasks(strategy, i),
                    "priority": priority,
                    "load_increase": 0.3 if i == 0 else 0.1
                }
                reassignment_plan["assignments"].append(assignment)
        
        elif action == 2:
            # Load balancing strategy
            strategy = "load_balanced"
            reassignment_plan["strategy"] = strategy
            
            for i, healthy_sc in enumerate(self.healthy_spacecraft):
                load_factor = 1.0 - healthy_sc['load_capacity']
                assignment = {
                    "spacecraft": healthy_sc["name"],
                    "new_tasks": self._get_redistributed_tasks(strategy, i),
                    "priority": "normal",
                    "load_increase": load_factor * 0.25
                }
                reassignment_plan["assignments"].append(assignment)
        
        return reassignment_plan
    
    def _get_redistributed_tasks(self, strategy: str, spacecraft_index: int) -> List[str]:
        """Get the list of redistributed tasks for a spacecraft"""
        
        # Tasks that need to be redistributed from faulty spacecraft
        orphaned_tasks = []
        
        for faulty_sc in self.faulty_spacecraft:
            lost_capabilities = faulty_sc['capabilities_lost']
            
            # Map lost capabilities to specific tasks
            for capability in lost_capabilities:
                if capability == "attitude_control":
                    orphaned_tasks.extend(["pointing_control", "stabilization"])
                elif capability == "sensing":
                    orphaned_tasks.extend(["target_observation", "environmental_monitoring"])
                elif capability == "communication":
                    orphaned_tasks.extend(["data_relay", "ground_communication"])
                elif capability == "navigation":
                    orphaned_tasks.extend(["position_determination", "orbit_maintenance"])
        
        # Remove duplicates
        orphaned_tasks = list(set(orphaned_tasks))
        
        # Distribute tasks based on strategy
        if strategy == "even_distribution":
            # Divide tasks evenly
            tasks_per_sc = len(orphaned_tasks) // len(self.healthy_spacecraft)
            start_idx = spacecraft_index * tasks_per_sc
            end_idx = start_idx + tasks_per_sc
            return orphaned_tasks[start_idx:end_idx]
        
        elif strategy == "capability_based":
            # Give more critical tasks to first (most capable) spacecraft
            if spacecraft_index == 0:
                return orphaned_tasks[:len(orphaned_tasks)//2]  # Give half to most capable
            else:
                remaining_tasks = orphaned_tasks[len(orphaned_tasks)//2:]
                tasks_per_sc = len(remaining_tasks) // (len(self.healthy_spacecraft) - 1)
                start_idx = (spacecraft_index - 1) * tasks_per_sc
                end_idx = start_idx + tasks_per_sc
                return remaining_tasks[start_idx:end_idx]
        
        elif strategy == "load_balanced":
            # Distribute based on current load capacity
            total_capacity = sum(sc['load_capacity'] for sc in self.healthy_spacecraft)
            sc_capacity = self.healthy_spacecraft[spacecraft_index]['load_capacity']
            proportion = sc_capacity / total_capacity
            num_tasks = int(len(orphaned_tasks) * proportion)
            start_idx = spacecraft_index * num_tasks
            return orphaned_tasks[start_idx:start_idx + num_tasks]
        
        return []
    
    def _fallback_reassignment(self) -> Dict:
        """Fallback reassignment when DRL is not available"""
        
        print("Using fallback task reassignment (rule-based)")
        
        reassignment_plan = {
            "timestamp": datetime.now().isoformat(),
            "trigger": "fallback_rules",
            "action_value": -1,
            "assignments": [],
            "strategy": "simple_redistribution"
        }
        
        # Simple rule: distribute tasks evenly among healthy spacecraft
        for i, healthy_sc in enumerate(self.healthy_spacecraft):
            assignment = {
                "spacecraft": healthy_sc["name"],
                "new_tasks": ["backup_attitude_control", "redundant_sensing"],
                "priority": "normal",
                "load_increase": 0.15
            }
            reassignment_plan["assignments"].append(assignment)
        
        return reassignment_plan
    
    def execute_reassignment(self, reassignment_plan: Dict) -> bool:
        """Execute the task reassignment plan"""
        
        print(f"\nExecuting task reassignment: {reassignment_plan['strategy']}")
        
        try:
            for assignment in reassignment_plan["assignments"]:
                spacecraft = assignment["spacecraft"]
                new_tasks = assignment["new_tasks"]
                priority = assignment["priority"]
                
                print(f"   {spacecraft}: +{len(new_tasks)} tasks ({priority} priority)")
                
                # Here you would interface with your Basilisk simulation
                # to actually reconfigure the spacecraft tasks
                self._apply_task_assignment(spacecraft, new_tasks, priority)
            
            print("Task reassignment executed successfully")
            return True
            
        except Exception as e:
            print(f"Error executing task reassignment: {e}")
            return False
    
    def _apply_task_assignment(self, spacecraft: str, tasks: List[str], priority: str):
        """Apply task assignment to specific spacecraft (interface with Basilisk)"""
        
        # This is where you would interface with your Basilisk simulation
        # to reconfigure the spacecraft's mission parameters
        
        print(f"      Applying {len(tasks)} tasks to {spacecraft}")
        
        # Example: Update spacecraft configuration
        # spacecraft_config = self.constellation_config.get_spacecraft_config(spacecraft)
        # spacecraft_config.update_tasks(tasks, priority)
        
        # For now, just store the assignment
        if spacecraft not in self.task_assignments:
            self.task_assignments[spacecraft] = []
        
        self.task_assignments[spacecraft].extend([
            {"task": task, "priority": priority, "assigned_at": datetime.now().isoformat()}
            for task in tasks
        ])


def integrate_fault_detection_with_drl(scenario, scenario_config, output_dir="."):
    """
    Main integration function that combines your fault detection with DRL reassignment
    """
    
    print("\n" + "="*80)
    print("INTEGRATED FAULT DETECTION + DRL TASK REASSIGNMENT")
    print("="*80)
    
    # Step 1: Run your existing fault detection
    print("\nSTEP 1: Running ML Fault Detection...")
    print("-" * 50)
    
    from real_ml_fault_detection import run_real_ml_detection_on_scenario
    
    fault_results = run_real_ml_detection_on_scenario(scenario, scenario_config, output_dir)
    
    if not fault_results or "error" in fault_results:
        print("Fault detection failed - cannot proceed with DRL reassignment")
        return None
    
    # Step 2: Initialize DRL task reassignment system
    print("\nSTEP 2: Initializing DRL Task Reassignment...")
    print("-" * 50)
    
    drl_system = DRLTaskReassignmentSystem(scenario_config)
    
    # Step 3: Update spacecraft status based on fault detection
    print("\nSTEP 3: Updating Spacecraft Status...")
    print("-" * 50)
    
    drl_system.update_spacecraft_status(fault_results["detections"])
    
    # Step 4: Trigger DRL-based task reassignment
    print("\nSTEP 4: DRL Decision Making...")
    print("-" * 50)
    
    reassignment_plan = drl_system.trigger_drl_reassignment()
    
    # Step 5: Execute task reassignment
    print("\nSTEP 5: Executing Task Reassignment...")
    print("-" * 50)
    
    success = drl_system.execute_reassignment(reassignment_plan)
    
    # Step 6: Generate integrated results
    integrated_results = {
        "fault_detection": fault_results,
        "drl_reassignment": {
            "plan": reassignment_plan,
            "executed": success,
            "healthy_spacecraft": drl_system.healthy_spacecraft,
            "faulty_spacecraft": drl_system.faulty_spacecraft,
            "task_assignments": drl_system.task_assignments
        },
        "summary": {
            "total_spacecraft": len(drl_system.healthy_spacecraft) + len(drl_system.faulty_spacecraft),
            "faults_detected": len([sc for sc in drl_system.faulty_spacecraft]),
            "tasks_reassigned": len(reassignment_plan.get("assignments", [])),
            "reassignment_strategy": reassignment_plan.get("strategy", "unknown"),
            "execution_success": success
        }
    }
    
    # Save integrated results
    _save_integrated_results(integrated_results, output_dir)
    
    print(f"\nINTEGRATED SYSTEM COMPLETE!")
    print(f"   Faults detected: {integrated_results['summary']['faults_detected']}")
    print(f"   Tasks reassigned: {integrated_results['summary']['tasks_reassigned']}")
    print(f"   Strategy: {integrated_results['summary']['reassignment_strategy']}")
    print(f"   Success: {integrated_results['summary']['execution_success']}")
    
    return integrated_results


def _save_integrated_results(results: Dict, output_dir: str):
    """Save integrated fault detection + DRL results"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save JSON results
    json_path = os.path.join(output_dir, "integrated_fault_drl_results.json")
    try:
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Integrated results saved: {json_path}")
    except Exception as e:
        print(f"Error saving JSON results: {e}")
    
    # Save text report
    report_path = os.path.join(output_dir, "integrated_fault_drl_report.txt")
    try:
        with open(report_path, 'w') as f:
            f.write("INTEGRATED FAULT DETECTION + DRL TASK REASSIGNMENT REPORT\n")
            f.write("=" * 65 + "\n\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Fault Detection Summary
            f.write("FAULT DETECTION SUMMARY:\n")
            fd_summary = results["fault_detection"]["summary"]
            f.write(f"  Total Spacecraft: {fd_summary['total_spacecraft']}\n")
            f.write(f"  Faults Detected: {fd_summary['total_detections']}\n")
            f.write(f"  Success Rate: {fd_summary['success_rate']:.1%}\n\n")
            
            # DRL Reassignment Summary
            f.write("DRL TASK REASSIGNMENT SUMMARY:\n")
            drl_summary = results["summary"]
            f.write(f"  Reassignment Strategy: {drl_summary['reassignment_strategy']}\n")
            f.write(f"  Tasks Reassigned: {drl_summary['tasks_reassigned']}\n")
            f.write(f"  Execution Success: {drl_summary['execution_success']}\n\n")
            
            # Detailed Assignments
            f.write("TASK ASSIGNMENTS:\n")
            for assignment in results["drl_reassignment"]["plan"].get("assignments", []):
                f.write(f"  {assignment['spacecraft']}:\n")
                f.write(f"    Tasks: {assignment['new_tasks']}\n")
                f.write(f"    Priority: {assignment['priority']}\n")
                f.write(f"    Load Increase: {assignment['load_increase']:.2f}\n")
        
        print(f"Integrated report saved: {report_path}")
        
    except Exception as e:
        print(f"Error saving integrated report: {e}")


if __name__ == "__main__":
    print("DRL Integration Bridge for Basilisk Fault Detection")
    print("This module connects ML fault detection with DRL task reassignment")
    print()
    print("Usage:")
    print("  1. Copy DRL files to bsk_rl/ folder")
    print("  2. Run: integrate_fault_detection_with_drl(scenario, config, output_dir)")
    print("  3. View integrated results")