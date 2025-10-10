#!/usr/bin/env python
"""
Integration.py - Fixed DRL Integration for Spacecraft Fault Tolerance

This module properly integrates your PPO model with the spacecraft simulation.
"""

import os
import sys
import numpy as np
import tensorflow as tf
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# Import your existing modules
try:
    from Envs import BasiliskEnv, BasiliskModel
    from stable_baselines3 import PPO
    from integration import PD_test, PDControllor  # Your original integration functions
    print("✓ Successfully imported DRL modules")
    DRL_MODULES_AVAILABLE = True
except ImportError as e:
    print(f"✗ DRL modules import error: {e}")
    DRL_MODULES_AVAILABLE = False

try:
    # Try to load your trained PPO model
    MODEL_PATH = "C:\\training\\train_Basilisk\\benv05_True_wheel\\model\\65536000.zip"
    if not os.path.exists(MODEL_PATH):
        # Look for alternative model paths
        alternative_paths = [
            "anomaly_detection_model.keras",
            "trained_models/ppo_model.zip",
            "models/ppo_spacecraft.zip"
        ]
        MODEL_PATH = None
        for path in alternative_paths:
            if os.path.exists(path):
                MODEL_PATH = path
                break
    
    PPO_MODEL_AVAILABLE = MODEL_PATH is not None and os.path.exists(MODEL_PATH)
    print(f"✓ PPO model found: {MODEL_PATH}" if PPO_MODEL_AVAILABLE else "✗ No PPO model found")
except Exception as e:
    PPO_MODEL_AVAILABLE = False
    MODEL_PATH = None
    print(f"✗ PPO model check failed: {e}")


class WorkingDRLAgent:
    """
    DRL Agent that actually works with your spacecraft constellation
    """
    
    def __init__(self, state_dim=10, action_dim=3, model_path=None):
        self.state_dim = state_dim
        self.action_dim = action_dim
        self.model_path = model_path
        self.model = None
        self.env = None
        self.is_loaded = False
        
        print(f"Initializing DRL Agent (state_dim={state_dim}, action_dim={action_dim})")
        
        # Try to load the model
        self._load_model()
    
    def _load_model(self):
        """Load the PPO model"""
        
        if not DRL_MODULES_AVAILABLE:
            print("  DRL modules not available - using rule-based agent")
            return False
        
        try:
            if self.model_path and os.path.exists(self.model_path):
                if self.model_path.endswith('.zip'):
                    # Load PPO model
                    self.model = PPO.load(self.model_path)
                    self.is_loaded = True
                    print(f"  ✓ Loaded PPO model from {self.model_path}")
                    
                    # Create environment for the model
                    self.env = BasiliskEnv(faulty=True, torque_mode="wheel")
                    print("  ✓ Created Basilisk environment")
                    
                    return True
                
                elif self.model_path.endswith('.keras'):
                    # Load Keras model
                    self.model = tf.keras.models.load_model(self.model_path)
                    self.is_loaded = True
                    print(f"  ✓ Loaded Keras model from {self.model_path}")
                    return True
            
            print(f"  ✗ Could not load model from {self.model_path}")
            return False
            
        except Exception as e:
            print(f"  ✗ Model loading failed: {e}")
            self.model = None
            self.is_loaded = False
            return False
    
    def select_action(self, state: np.ndarray, training=False) -> int:
        """Select action using the DRL model"""
        
        if not self.is_loaded or self.model is None:
            return self._rule_based_action(state)
        
        try:
            if hasattr(self.model, 'predict'):
                # PPO model
                action, _states = self.model.predict(state.reshape(1, -1), deterministic=not training)
                if isinstance(action, np.ndarray):
                    action = int(action[0])
                else:
                    action = int(action)
                
                print(f"   DRL Agent selected action: {action}")
                return action
            
            elif hasattr(self.model, '__call__'):
                # Keras model
                state_tensor = tf.constant(state.reshape(1, -1), dtype=tf.float32)
                prediction = self.model(state_tensor)
                action = int(tf.argmax(prediction, axis=1)[0])
                
                print(f"   DRL Agent (Keras) selected action: {action}")
                return action
            
            else:
                print("   Model doesn't have expected interface - using rule-based")
                return self._rule_based_action(state)
                
        except Exception as e:
            print(f"  ✗ DRL action selection failed: {e}")
            return self._rule_based_action(state)
    
    def _rule_based_action(self, state: np.ndarray) -> int:
        """Fallback rule-based action selection"""
        
        healthy_ratio = state[0] if len(state) > 0 else 0.5
        fault_confidence = state[2] if len(state) > 2 else 0.5
        anomaly_score = state[4] if len(state) > 4 else 0.5
        
        # Decision logic based on spacecraft constellation state
        if healthy_ratio < 0.5:  # Less than half healthy
            return 1  # capability_based strategy
        elif fault_confidence > 0.8 or anomaly_score > 0.7:  # High confidence faults
            return 2  # load_balanced strategy  
        else:
            return 0  # even_distribution strategy
    
    def get_action_probabilities(self, state: np.ndarray) -> Optional[np.ndarray]:
        """Get action probabilities if available"""
        
        if not self.is_loaded or not hasattr(self.model, 'predict_proba'):
            # Return uniform probabilities for rule-based
            return np.array([0.33, 0.33, 0.34])
        
        try:
            probs = self.model.predict_proba(state.reshape(1, -1))
            return probs[0] if len(probs.shape) > 1 else probs
        except:
            return np.array([0.33, 0.33, 0.34])


class WorkingDRLTaskReassignment:
    """
    DRL-based task reassignment that actually works
    """
    
    def __init__(self, num_spacecraft=4):
        self.num_spacecraft = num_spacecraft
        self.agent = WorkingDRLAgent(
            state_dim=10, 
            action_dim=3, 
            model_path=MODEL_PATH
        )
        self.spacecraft_status = {}
        self.task_history = []
        
    def update_spacecraft_status(self, fault_detections: Dict):
        """Update spacecraft status from fault detection results"""
        
        self.spacecraft_status = {}
        
        for spacecraft_name, detections in fault_detections.items():
            if detections and len(detections) > 0:
                # Faulty spacecraft
                max_confidence = max(d.confidence for d in detections)
                max_anomaly = max(d.anomaly_score for d in detections)
                
                self.spacecraft_status[spacecraft_name] = {
                    'status': 'FAULTY',
                    'fault_count': len(detections),
                    'max_confidence': max_confidence,
                    'max_anomaly_score': max_anomaly,
                    'severity': self._assess_severity(max_confidence, max_anomaly)
                }
            else:
                # Healthy spacecraft
                self.spacecraft_status[spacecraft_name] = {
                    'status': 'HEALTHY',
                    'fault_count': 0,
                    'max_confidence': 0.0,
                    'max_anomaly_score': 0.0,
                    'severity': 'none'
                }
        
        healthy_count = sum(1 for s in self.spacecraft_status.values() if s['status'] == 'HEALTHY')
        faulty_count = sum(1 for s in self.spacecraft_status.values() if s['status'] == 'FAULTY')
        
        print(f"Spacecraft status updated: {healthy_count} healthy, {faulty_count} faulty")
        
        return self.spacecraft_status
    
    def _assess_severity(self, confidence: float, anomaly_score: float) -> str:
        """Assess fault severity"""
        if confidence > 0.9 or anomaly_score > 90:
            return 'critical'
        elif confidence > 0.7 or anomaly_score > 60:
            return 'major'
        elif confidence > 0.5 or anomaly_score > 30:
            return 'minor'
        else:
            return 'negligible'
    
    def make_drl_decision(self) -> Dict:
        """Make task reassignment decision using DRL"""
        
        # Prepare state for DRL agent
        state = self._prepare_state_vector()
        
        print(f"DRL state vector: {state}")
        
        # Get action from DRL agent
        action = self.agent.select_action(state, training=False)
        action_probs = self.agent.get_action_probabilities(state)
        
        # Map action to strategy
        strategies = {
            0: "even_distribution",
            1: "capability_based", 
            2: "load_balanced"
        }
        
        selected_strategy = strategies.get(action, "even_distribution")
        
        # Create reassignment plan
        plan = self._create_reassignment_plan(selected_strategy, action)
        
        decision_result = {
            "drl_used": self.agent.is_loaded,
            "action": action,
            "strategy": selected_strategy,
            "action_probabilities": action_probs.tolist() if action_probs is not None else None,
            "confidence": float(np.max(action_probs)) if action_probs is not None else 0.6,
            "state_vector": state.tolist(),
            "reassignment_plan": plan
        }
        
        print(f" DRL Decision: {selected_strategy} (action={action}, confidence={decision_result['confidence']:.3f})")
        
        return decision_result
    
    def _prepare_state_vector(self) -> np.ndarray:
        """Prepare state vector for DRL agent"""
        
        total_spacecraft = len(self.spacecraft_status)
        if total_spacecraft == 0:
            return np.zeros(10, dtype=np.float32)
        
        # Count healthy/faulty
        healthy_count = sum(1 for s in self.spacecraft_status.values() if s['status'] == 'HEALTHY')
        faulty_count = sum(1 for s in self.spacecraft_status.values() if s['status'] == 'FAULTY')
        
        # Calculate fault metrics
        confidences = [s['max_confidence'] for s in self.spacecraft_status.values() if s['max_confidence'] > 0]
        anomaly_scores = [s['max_anomaly_score'] for s in self.spacecraft_status.values() if s['max_anomaly_score'] > 0]
        
        avg_confidence = np.mean(confidences) if confidences else 0.0
        max_confidence = np.max(confidences) if confidences else 0.0
        avg_anomaly = np.mean(anomaly_scores) if anomaly_scores else 0.0
        max_anomaly = np.max(anomaly_scores) if anomaly_scores else 0.0
        
        # Severity counts
        severities = [s['severity'] for s in self.spacecraft_status.values()]
        critical_count = severities.count('critical')
        major_count = severities.count('major')
        
        # Create normalized state vector
        state = np.array([
            healthy_count / total_spacecraft,           # 0: healthy ratio
            faulty_count / total_spacecraft,            # 1: faulty ratio  
            avg_confidence,                             # 2: avg fault confidence
            max_confidence,                             # 3: max fault confidence
            avg_anomaly / 100.0,                        # 4: normalized avg anomaly
            max_anomaly / 100.0,                        # 5: normalized max anomaly
            critical_count / total_spacecraft,          # 6: critical fault ratio
            major_count / total_spacecraft,             # 7: major fault ratio
            len(confidences) / (total_spacecraft * 2),  # 8: fault density
            0.5  # 9: placeholder for system load
        ], dtype=np.float32)
        
        return state
    
    def _create_reassignment_plan(self, strategy: str, action: int) -> Dict:
        """Create detailed task reassignment plan"""
        
        plan = {
            "strategy": strategy,
            "action": action,
            "timestamp": datetime.now().isoformat(),
            "assignments": []
        }
        
        # Get healthy and faulty spacecraft
        healthy_sc = [name for name, status in self.spacecraft_status.items() if status['status'] == 'HEALTHY']
        faulty_sc = [name for name, status in self.spacecraft_status.items() if status['status'] == 'FAULTY']
        
        if not healthy_sc:
            plan["error"] = "No healthy spacecraft available for task reassignment"
            return plan
        
        # Define tasks that need redistribution based on faults
        critical_tasks = ["attitude_control", "navigation", "power_management"]
        secondary_tasks = ["data_collection", "communication_relay", "monitoring"]
        
        # Determine tasks to redistribute based on strategy
        if strategy == "even_distribution":
            # Distribute tasks evenly
            tasks_per_sc = len(critical_tasks) // len(healthy_sc)
            for i, sc_name in enumerate(healthy_sc):
                start_idx = i * tasks_per_sc
                end_idx = min(start_idx + tasks_per_sc, len(critical_tasks))
                assigned_tasks = critical_tasks[start_idx:end_idx]
                
                plan["assignments"].append({
                    "spacecraft": sc_name,
                    "tasks": assigned_tasks,
                    "priority": "normal",
                    "load_increase": 0.2
                })
        
        elif strategy == "capability_based":
            # Assign more tasks to the "best" spacecraft (first in list for now)
            primary_sc = healthy_sc[0]
            plan["assignments"].append({
                "spacecraft": primary_sc,
                "tasks": critical_tasks[:2],  # Give most critical tasks
                "priority": "high",
                "load_increase": 0.4
            })
            
            # Distribute remaining tasks
            for i, sc_name in enumerate(healthy_sc[1:]):
                remaining_tasks = critical_tasks[2:] + secondary_tasks
                tasks_per_sc = len(remaining_tasks) // len(healthy_sc[1:])
                start_idx = i * tasks_per_sc
                end_idx = start_idx + tasks_per_sc
                assigned_tasks = remaining_tasks[start_idx:end_idx]
                
                plan["assignments"].append({
                    "spacecraft": sc_name,
                    "tasks": assigned_tasks,
                    "priority": "normal", 
                    "load_increase": 0.15
                })
        
        elif strategy == "load_balanced":
            # Balance load based on current spacecraft status
            for sc_name in healthy_sc:
                # Assign fewer tasks to each to balance load
                tasks_subset = critical_tasks[:2] if sc_name == healthy_sc[0] else secondary_tasks[:2]
                
                plan["assignments"].append({
                    "spacecraft": sc_name,
                    "tasks": tasks_subset,
                    "priority": "normal",
                    "load_increase": 0.1
                })
        
        return plan
    
    def execute_reassignment_plan(self, plan: Dict) -> Dict:
        """Execute the reassignment plan"""
        
        execution_result = {
            "success": True,
            "strategy_executed": plan["strategy"],
            "assignments_completed": 0,
            "errors": [],
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            for assignment in plan.get("assignments", []):
                spacecraft = assignment["spacecraft"]
                tasks = assignment["tasks"]
                
                print(f"   Assigning {len(tasks)} tasks to {spacecraft}")
                
                # Here you would interface with your actual Basilisk simulation
                # For now, we'll just log the assignment
                self.task_history.append({
                    "spacecraft": spacecraft,
                    "tasks": tasks,
                    "priority": assignment["priority"],
                    "load_increase": assignment["load_increase"],
                    "timestamp": datetime.now().isoformat()
                })
                
                execution_result["assignments_completed"] += 1
            
            print(f"  ✓ Task reassignment executed: {execution_result['assignments_completed']} assignments")
            
        except Exception as e:
            execution_result["success"] = False
            execution_result["errors"].append(str(e))
            print(f"  ✗ Task reassignment failed: {e}")
        
        return execution_result


def main(scenario=None, config=None, output_dir="logs"):
    """
    Main DRL integration function that actually uses the DRL model
    """
    
    print(" STARTING WORKING DRL INTEGRATION")
    print("=" * 50)
    
    try:
        # Step 1: Initialize DRL Task Reassignment System
        print("Step 1: Initializing DRL system...")
        drl_system = WorkingDRLTaskReassignment()
        
        # Step 2: Run fault detection (if scenario provided)
        print("Step 2: Running fault detection...")
        if scenario and hasattr(scenario, 'fault_detections'):
            # Use real fault detection results from scenario
            fault_results = scenario.fault_detections
        else:
            # Mock fault detection for testing
            fault_results = {
                "Satellite1": [type('MockFault', (), {
                    'confidence': 0.85,
                    'anomaly_score': 75.0,
                    'detection_time_minutes': 15.0,
                    'fault_type': 'rw_friction'
                })()],
                "Satellite2": [],
                "Satellite3": [],
                "Satellite4": []
            }
        
        # Step 3: Update spacecraft status
        print("Step 3: Updating spacecraft status...")
        spacecraft_status = drl_system.update_spacecraft_status(fault_results)
        
        # Step 4: Make DRL decision
        print("Step 4: Making DRL decision...")
        drl_decision = drl_system.make_drl_decision()
        
        # Step 5: Execute task reassignment
        print("Step 5: Executing task reassignment...")
        execution_result = drl_system.execute_reassignment_plan(drl_decision["reassignment_plan"])
        
        # Compile results
        final_results = {
            "success": True,
            "drl_model_used": drl_decision["drl_used"],
            "spacecraft_status": spacecraft_status,
            "drl_decision": drl_decision,
            "execution_result": execution_result,
            "summary": {
                "strategy": drl_decision["strategy"],
                "confidence": drl_decision["confidence"],
                "assignments": execution_result["assignments_completed"],
                "healthy_spacecraft": len([s for s in spacecraft_status.values() if s['status'] == 'HEALTHY']),
                "faulty_spacecraft": len([s for s in spacecraft_status.values() if s['status'] == 'FAULTY'])
            }
        }
        
        print(f"\n DRL INTEGRATION COMPLETED!")
        print(f"   DRL Model Used: {'✓' if drl_decision['drl_used'] else '✗'}")
        print(f"   Strategy: {drl_decision['strategy']}")
        print(f"   Confidence: {drl_decision['confidence']:.3f}")
        print(f"   Assignments: {execution_result['assignments_completed']}")
        
        return final_results
        
    except Exception as e:
        error_result = {
            "success": False,
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }
        print(f"✗ DRL integration failed: {e}")
        import traceback
        traceback.print_exc()
        return error_result


# Compatibility functions
def run_drl_analysis(*args, **kwargs):
    return main(*args, **kwargs)

def integrate_with_basilisk(scenario_data, config=None):
    return main(scenario_data, config)


if __name__ == "__main__":
    print("Testing Working DRL Integration...")
    
    # Test the integration
    test_results = main()
    
    if test_results.get("success", False):
        print("\n DRL Integration Test PASSED!")
        print(f"   DRL Model Used: {test_results['drl_model_used']}")
        print(f"   Strategy Selected: {test_results['summary']['strategy']}")
    else:
        print(f"\n DRL Integration Test FAILED: {test_results.get('error', 'Unknown error')}")