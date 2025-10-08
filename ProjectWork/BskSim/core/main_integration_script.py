#!/usr/bin/env python
"""
main_integration.py

Main integration script that combines:
1. Your existing Basilisk fault detection
2. DRL-based task reassignment 
3. Integration with cluster constellation system

This script should be called from your GUI after simulation completion.
"""

import os
import sys
import numpy as np
import json
from datetime import datetime
from typing import Dict, List, Optional

# Add paths for all components
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
bsk_rl_dir = os.path.join(current_dir, 'bsk_rl')

for path in [parent_dir, current_dir, bsk_rl_dir]:
    if path not in sys.path:
        sys.path.insert(0, path)

# Import your existing components
try:
    from real_ml_fault_detection import run_real_ml_detection_on_scenario
    from drl_integration_bridge import DRLTaskReassignmentSystem, integrate_fault_detection_with_drl
    from drl_config import ConfigManager, create_config_manager
    FAULT_DETECTION_AVAILABLE = True
    print("Fault detection components loaded")
except ImportError as e:
    FAULT_DETECTION_AVAILABLE = False
    print(f"Fault detection components not available: {e}")

# Import DRL task reassignment module
try:
    from simple_task_drl import ClusterTaskReassignmentDRL
    SIMPLE_DRL_AVAILABLE = True
    print("DRL task reassignment module loaded")
except ImportError as e:
    SIMPLE_DRL_AVAILABLE = False
    print(f"DRL task reassignment not available: {e}")

# Import TensorFlow for ML model
try:
    import tensorflow as tf
    ML_AVAILABLE = True
    print("TensorFlow available")
except ImportError:
    ML_AVAILABLE = False
    print("TensorFlow not available")


class IntegratedFaultDetectionDRLSystem:
    """
    Main class that integrates fault detection with DRL task reassignment
    """
    
    def __init__(self, config_profile: str = "development"):
        """Initialize the integrated system"""
        
        self.config_manager = create_config_manager(config_profile)
        self.config = self.config_manager.config
        
        # System components
        self.fault_detector = None
        self.drl_agent = None
        self.task_reassignment_system = None
        self.environment = None
        
        # Results storage
        self.simulation_results = {}
        self.fault_detection_results = {}
        self.drl_results = {}
        self.integrated_results = {}
        
        print(f"Integrated System initialized with '{config_profile}' profile")
        
        # Initialize components
        self._initialize_components()
    
    def _initialize_components(self):
        """Initialize all system components"""
        
        print("\nInitializing system components...")
        
        # Initialize fault detection
        if FAULT_DETECTION_AVAILABLE:
            try:
                from real_ml_fault_detection import RealMLFaultDetector
                self.fault_detector = RealMLFaultDetector()
                print("  Fault detector initialized")
            except Exception as e:
                print(f"  Fault detector initialization failed: {e}")
        
        # Initialize DRL agent for task reassignment
        if SIMPLE_DRL_AVAILABLE:
            try:
                # Check if a saved model exists
                model_path = "drl_task_reassignment.h5"
                
                self.drl_agent = ClusterTaskReassignmentDRL(
                    state_dim=self.config.drl.state_dim,
                    action_dim=self.config.drl.action_dim,
                    model_path=model_path if os.path.exists(model_path) else None
                )
                
                if self.drl_agent.is_loaded:
                    print("  DRL agent initialized successfully")
                    print("    Architecture: Neural network (64-64-32 layers)")
                    print("    Training: Expert rule-based supervision")
                    print("    Actions: 3 (even_distribution, capability_based, load_balanced)")
                    
                    # Save the trained model for future use
                    if not os.path.exists(model_path):
                        self.drl_agent.save_model(model_path)
                else:
                    print("  DRL agent initialization failed")
                    self.drl_agent = None
                    
            except Exception as e:
                print(f"  DRL initialization failed: {e}")
                import traceback
                traceback.print_exc()
                self.drl_agent = None
        else:
            print("  DRL not available - install TensorFlow to enable")
            self.drl_agent = None
        
        print("Component initialization complete\n")
    
    def run_complete_analysis(self, scenario, scenario_config=None, output_dir="logs") -> Dict:
        """
        Run the complete integrated analysis:
        1. Fault detection
        2. DRL decision making
        3. Task reassignment
        4. Results integration
        """
        
        print("="*80)
        print("STARTING COMPLETE INTEGRATED ANALYSIS")
        print("="*80)
        
        start_time = datetime.now()
        
        # Create output directory
        timestamp = start_time.strftime("%Y%m%d_%H%M%S")
        analysis_output_dir = os.path.join(output_dir, f"integrated_analysis_{timestamp}")
        os.makedirs(analysis_output_dir, exist_ok=True)
        
        try:
            # Step 1: Run fault detection
            print("\nSTEP 1: FAULT DETECTION")
            print("-" * 50)
            
            self.fault_detection_results = self._run_fault_detection(
                scenario, scenario_config, analysis_output_dir
            )
            
            if not self.fault_detection_results or "error" in self.fault_detection_results:
                raise Exception("Fault detection failed")
            
            # Step 2: Analyze results and prepare for DRL
            print("\nSTEP 2: DRL ANALYSIS PREPARATION")
            print("-" * 50)
            
            drl_state = self._prepare_drl_state()
            
            # Step 3: Run DRL decision making
            print("\nSTEP 3: DRL DECISION MAKING")
            print("-" * 50)
            
            self.drl_results = self._run_drl_decision_making(drl_state, analysis_output_dir)
            
            # Step 4: Execute task reassignment
            print("\nSTEP 4: TASK REASSIGNMENT")
            print("-" * 50)
            
            self.reassignment_results = self._execute_task_reassignment(analysis_output_dir)
            
            # Step 5: Integration and reporting
            print("\nSTEP 5: RESULTS INTEGRATION")
            print("-" * 50)
            
            self.integrated_results = self._integrate_results(
                analysis_output_dir, start_time
            )
            
            # Step 6: Generate comprehensive report
            print("\nSTEP 6: REPORT GENERATION")
            print("-" * 50)
            
            self._generate_comprehensive_report(analysis_output_dir)
            
            print("\nINTEGRATED ANALYSIS COMPLETED SUCCESSFULLY")
            print(f"Results saved to: {analysis_output_dir}")
            
            return self.integrated_results
            
        except Exception as e:
            print(f"\nINTEGRATED ANALYSIS FAILED: {e}")
            import traceback
            traceback.print_exc()
            
            # Save partial results
            self._save_partial_results(analysis_output_dir, str(e))
            
            return {"error": str(e), "partial_results": True}
    
    def _run_fault_detection(self, scenario, scenario_config, output_dir) -> Dict:
        """Run ML-based fault detection"""
        
        if not FAULT_DETECTION_AVAILABLE:
            print("Fault detection not available - using mock results")
            return self._generate_mock_fault_results()
        
        try:
            # Use your existing fault detection system
            results = run_real_ml_detection_on_scenario(
                scenario, scenario_config, output_dir
            )
            
            if results and "detections" in results:
                total_detections = sum(len(detections) for detections in results["detections"].values())
                print(f"  Fault detection completed: {total_detections} faults detected")
            
            return results
            
        except Exception as e:
            print(f"  Fault detection error: {e}")
            return {"error": str(e)}
    
    def _prepare_drl_state(self) -> np.ndarray:
        """Prepare state vector for DRL agent based on fault detection results"""
        
        # Extract information from fault detection results
        detections = self.fault_detection_results.get("detections", {})
        
        total_spacecraft = len(detections)
        faulty_spacecraft = sum(1 for det_list in detections.values() if det_list)
        healthy_spacecraft = total_spacecraft - faulty_spacecraft
        
        # Calculate fault severity metrics
        all_confidences = []
        all_anomaly_scores = []
        
        for spacecraft_detections in detections.values():
            for detection in spacecraft_detections:
                all_confidences.append(detection.confidence)
                all_anomaly_scores.append(detection.anomaly_score)
        
        avg_confidence = np.mean(all_confidences) if all_confidences else 0.0
        max_confidence = np.max(all_confidences) if all_confidences else 0.0
        avg_anomaly_score = np.mean(all_anomaly_scores) if all_anomaly_scores else 0.0
        max_anomaly_score = np.max(all_anomaly_scores) if all_anomaly_scores else 0.0
        
        # Create normalized state vector
        state = np.array([
            healthy_spacecraft / total_spacecraft if total_spacecraft > 0 else 1.0,
            faulty_spacecraft / total_spacecraft if total_spacecraft > 0 else 0.0,
            avg_confidence,
            max_confidence,
            avg_anomaly_score / 100.0,  # Normalize to [0,1]
            max_anomaly_score / 100.0,
            len(all_confidences) / 10.0,  # Number of faults normalized
            np.random.normal(0.5, 0.1),   # Environmental factor (placeholder)
            np.random.normal(0.5, 0.1),   # System load (placeholder)
            0.0  # Padding
        ], dtype=np.float32)
        
        # Ensure state is correct size
        if len(state) != self.config.drl.state_dim:
            if len(state) < self.config.drl.state_dim:
                state = np.pad(state, (0, self.config.drl.state_dim - len(state)))
            else:
                state = state[:self.config.drl.state_dim]
        
        print(f"  DRL state prepared: {state.shape}")
        print(f"    Healthy/Faulty ratio: {healthy_spacecraft}/{faulty_spacecraft}")
        print(f"    Avg confidence: {avg_confidence:.3f}")
        print(f"    Max anomaly score: {max_anomaly_score:.1f}")
        
        return state
    
    def _run_drl_decision_making(self, state: np.ndarray, output_dir: str) -> Dict:
        """Run DRL agent to make task reassignment decisions"""
        
        if not SIMPLE_DRL_AVAILABLE or self.drl_agent is None:
            print("DRL not available - using rule-based decisions")
            return self._generate_rule_based_decisions(state)
        
        try:
            # Get action from DRL agent
            action = self.drl_agent.select_action(state, training=False)
            
            # Get action probabilities
            action_probs = self.drl_agent.get_action_probabilities(state)
            
            # Map action to strategy
            strategy_mapping = {
                0: "even_distribution",
                1: "capability_based", 
                2: "load_balanced"
            }
            
            selected_strategy = strategy_mapping.get(action, "even_distribution")
            
            drl_results = {
                "action": int(action),
                "strategy": selected_strategy,
                "action_probabilities": action_probs.tolist(),
                "state_vector": state.tolist(),
                "decision_confidence": float(np.max(action_probs)),
                "drl_model_used": True,
                "model_type": "Neural Network (Expert-Supervised)"
            }
            
            print(f"  DRL decision: Action {action} -> {selected_strategy}")
            print(f"    Confidence: {drl_results['decision_confidence']:.3f}")
            print(f"    Probabilities: even={action_probs[0]:.3f}, capability={action_probs[1]:.3f}, load={action_probs[2]:.3f}")
            
            return drl_results
            
        except Exception as e:
            print(f"  DRL decision making error: {e}")
            import traceback
            traceback.print_exc()
            return self._generate_rule_based_decisions(state)
    
    def _generate_rule_based_decisions(self, state: np.ndarray) -> Dict:
        """Generate rule-based decisions when DRL is not available"""
        
        # Simple rule-based logic
        healthy_ratio = state[0]
        fault_severity = state[4]  # Normalized avg anomaly score
        
        if healthy_ratio < 0.5:  # More than half spacecraft faulty
            strategy = "capability_based"  # Focus on best remaining spacecraft
            action = 1
        elif fault_severity > 0.7:  # High severity faults
            strategy = "load_balanced"  # Distribute load carefully
            action = 2
        else:
            strategy = "even_distribution"  # Normal operation
            action = 0
        
        return {
            "action": action,
            "strategy": strategy,
            "action_probabilities": None,
            "state_vector": state.tolist(),
            "decision_confidence": 0.6,  # Lower confidence for rule-based
            "drl_model_used": False
        }
    
    def _execute_task_reassignment(self, output_dir: str) -> Dict:
        """Execute the task reassignment based on DRL decisions"""
        
        strategy = self.drl_results["strategy"]
        
        # Initialize task reassignment system if not done already
        if self.task_reassignment_system is None:
            from drl_integration_bridge import DRLTaskReassignmentSystem
            self.task_reassignment_system = DRLTaskReassignmentSystem(None)
        
        # Update spacecraft status
        detections = self.fault_detection_results.get("detections", {})
        self.task_reassignment_system.update_spacecraft_status(detections)
        
        # Create reassignment plan based on DRL strategy
        reassignment_plan = self._create_reassignment_plan_from_strategy(strategy)
        
        # Execute the plan
        success = self.task_reassignment_system.execute_reassignment(reassignment_plan)
        
        results = {
            "reassignment_plan": reassignment_plan,
            "execution_success": success,
            "healthy_spacecraft": self.task_reassignment_system.healthy_spacecraft,
            "faulty_spacecraft": self.task_reassignment_system.faulty_spacecraft,
            "task_assignments": self.task_reassignment_system.task_assignments
        }
        
        print(f"  Task reassignment executed: {strategy}")
        print(f"    Success: {success}")
        print(f"    Assignments: {len(reassignment_plan.get('assignments', []))}")
        
        return results
    
    def _create_reassignment_plan_from_strategy(self, strategy: str) -> Dict:
        """Create detailed reassignment plan from DRL strategy"""
        
        plan = {
            "timestamp": datetime.now().isoformat(),
            "strategy": strategy,
            "trigger": "drl_decision",
            "assignments": []
        }
        
        healthy_sc = self.task_reassignment_system.healthy_spacecraft
        faulty_sc = self.task_reassignment_system.faulty_spacecraft
        
        if not healthy_sc:
            print("    Warning: No healthy spacecraft available for reassignment")
            return plan
        
        # Collect tasks from faulty spacecraft
        # Collect tasks from faulty spacecraft
        orphaned_tasks = []
        for faulty in faulty_sc:
            capabilities_lost = faulty.get("capabilities_lost", [])
            print(f"      Processing faulty spacecraft: {faulty.get('name', 'Unknown')}")
            print(f"      Capabilities lost: {capabilities_lost}")
            
            for capability in capabilities_lost:
                if capability == "attitude_control" or capability == "pointing_accuracy":
                    tasks = ["attitude_stabilization", "pointing_control", "momentum_management"]
                    orphaned_tasks.extend(tasks)
                    print(f"        Adding attitude tasks: {tasks}")
                elif capability == "sensing" or capability == "navigation":
                    tasks = ["target_observation", "data_collection", "sensor_calibration"]
                    orphaned_tasks.extend(tasks)
                    print(f"        Adding sensing tasks: {tasks}")
                elif capability == "communication" or capability == "data_relay":
                    tasks = ["data_relay", "ground_link", "inter_satellite_comm"]
                    orphaned_tasks.extend(tasks)
                    print(f"        Adding communication tasks: {tasks}")
                elif capability == "power_generation" or capability == "system_operations":
                    tasks = ["power_management", "thermal_control", "backup_operations"]
                    orphaned_tasks.extend(tasks)
                    print(f"        Adding power tasks: {tasks}")
                else:
                    # Default tasks for unknown capability
                    tasks = ["backup_attitude_control", "emergency_operations"]
                    orphaned_tasks.extend(tasks)
                    print(f"        Adding default tasks for '{capability}': {tasks}")

        print(f"      Total orphaned tasks collected: {len(orphaned_tasks)}")
        
        orphaned_tasks = list(set(orphaned_tasks))  # Remove duplicates
        
        if not orphaned_tasks:
            print("    Info: No tasks need reassignment")
            return plan
        
        # Distribute tasks based on strategy
        if strategy == "even_distribution":
            tasks_per_sc = len(orphaned_tasks) // len(healthy_sc)
            for i, sc in enumerate(healthy_sc):
                start_idx = i * tasks_per_sc
                end_idx = start_idx + tasks_per_sc
                assigned_tasks = orphaned_tasks[start_idx:end_idx]
                
                plan["assignments"].append({
                    "spacecraft": sc["name"],
                    "new_tasks": assigned_tasks,
                    "priority": "normal",
                    "load_increase": 0.2
                })
        
        elif strategy == "capability_based":
            # Assign more tasks to most capable (highest load capacity)
            sorted_sc = sorted(healthy_sc, key=lambda x: x.get("load_capacity", 0.5), reverse=True)
            
            # Give 60% of tasks to most capable, distribute rest
            primary_tasks = orphaned_tasks[:int(len(orphaned_tasks) * 0.6)]
            remaining_tasks = orphaned_tasks[int(len(orphaned_tasks) * 0.6):]
            
            # Primary spacecraft gets most tasks
            plan["assignments"].append({
                "spacecraft": sorted_sc[0]["name"],
                "new_tasks": primary_tasks,
                "priority": "high",
                "load_increase": 0.4
            })
            
            # Distribute remaining tasks
            if len(sorted_sc) > 1 and remaining_tasks:
                tasks_per_backup = len(remaining_tasks) // (len(sorted_sc) - 1)
                for i, sc in enumerate(sorted_sc[1:]):
                    start_idx = i * tasks_per_backup
                    end_idx = start_idx + tasks_per_backup
                    assigned_tasks = remaining_tasks[start_idx:end_idx]
                    
                    plan["assignments"].append({
                        "spacecraft": sc["name"],
                        "new_tasks": assigned_tasks,
                        "priority": "normal",
                        "load_increase": 0.15
                    })
        
        elif strategy == "load_balanced":
            # Distribute based on available capacity
            total_capacity = sum(sc.get("load_capacity", 0.5) for sc in healthy_sc)
            
            for sc in healthy_sc:
                capacity_ratio = sc.get("load_capacity", 0.5) / total_capacity
                num_tasks = int(len(orphaned_tasks) * capacity_ratio)
                assigned_tasks = orphaned_tasks[:num_tasks]
                orphaned_tasks = orphaned_tasks[num_tasks:]  # Remove assigned tasks
                
                plan["assignments"].append({
                    "spacecraft": sc["name"],
                    "new_tasks": assigned_tasks,
                    "priority": "normal",
                    "load_increase": capacity_ratio * 0.3
                })
        
        return plan
    
    def _integrate_results(self, output_dir: str, start_time: datetime) -> Dict:
        """Integrate all results into comprehensive analysis"""
        
        end_time = datetime.now()
        analysis_duration = (end_time - start_time).total_seconds()
        
        # Calculate summary metrics
        fault_summary = self.fault_detection_results.get("summary", {})
        total_spacecraft = fault_summary.get("total_spacecraft", 0)
        total_faults = fault_summary.get("total_detections", 0)
        
        reassignment_summary = self.drl_results
        
        integrated_results = {
            "analysis_metadata": {
                "start_time": start_time.isoformat(),
                "end_time": end_time.isoformat(),
                "duration_seconds": analysis_duration,
                "config_profile": "development",
                "components_used": {
                    "fault_detection": FAULT_DETECTION_AVAILABLE,
                    "drl_available": SIMPLE_DRL_AVAILABLE,
                    "ml_model": ML_AVAILABLE
                }
            },
            "fault_detection_results": self.fault_detection_results,
            "drl_decision_results": self.drl_results,
            "task_reassignment_results": self.reassignment_results,
            "performance_metrics": {
                "total_spacecraft_analyzed": total_spacecraft,
                "total_faults_detected": total_faults,
                "fault_detection_success_rate": fault_summary.get("success_rate", 0.0),
                "drl_decision_confidence": reassignment_summary.get("decision_confidence", 0.0),
                "task_reassignment_success": self.reassignment_results.get("execution_success", False),
                "analysis_duration_minutes": analysis_duration / 60.0
            },
            "system_health": {
                "healthy_spacecraft_count": len(self.task_reassignment_system.healthy_spacecraft),
                "faulty_spacecraft_count": len(self.task_reassignment_system.faulty_spacecraft),
                "tasks_reassigned": len(self.reassignment_results.get("reassignment_plan", {}).get("assignments", [])),
                "overall_system_status": self._assess_overall_system_health()
            }
        }
        
        print(f"  Results integrated successfully")
        print(f"    Analysis duration: {analysis_duration:.1f} seconds")
        print(f"    System health: {integrated_results['system_health']['overall_system_status']}")
        
        return integrated_results
    
    def _assess_overall_system_health(self) -> str:
        """Assess overall system health based on all analysis results"""
        
        if not hasattr(self, 'task_reassignment_system') or self.task_reassignment_system is None:
            return "unknown"
        
        healthy_count = len(self.task_reassignment_system.healthy_spacecraft)
        faulty_count = len(self.task_reassignment_system.faulty_spacecraft)
        total_count = healthy_count + faulty_count
        
        if total_count == 0:
            return "unknown"
        
        healthy_ratio = healthy_count / total_count
        
        if healthy_ratio >= 0.8:
            return "excellent"
        elif healthy_ratio >= 0.6:
            return "good"
        elif healthy_ratio >= 0.4:
            return "degraded"
        elif healthy_ratio >= 0.2:
            return "critical"
        else:
            return "failure"
    
    def _generate_comprehensive_report(self, output_dir: str):
        """Generate comprehensive analysis report"""
        
        report_path = os.path.join(output_dir, "comprehensive_analysis_report.txt")
        
        try:
            with open(report_path, 'w') as f:
                f.write("COMPREHENSIVE FAULT DETECTION + DRL ANALYSIS REPORT\n")
                f.write("=" * 65 + "\n\n")
                
                # Metadata
                metadata = self.integrated_results["analysis_metadata"]
                f.write(f"Analysis Date: {metadata['start_time']}\n")
                f.write(f"Duration: {metadata['duration_seconds']:.1f} seconds\n")
                f.write(f"Components Used:\n")
                for comp, available in metadata["components_used"].items():
                    status = "AVAILABLE" if available else "UNAVAILABLE"
                    f.write(f"  {status}: {comp.replace('_', ' ').title()}\n")
                f.write("\n")
                
                # Fault Detection Summary
                f.write("FAULT DETECTION SUMMARY:\n")
                f.write("-" * 30 + "\n")
                fd_results = self.fault_detection_results
                if "summary" in fd_results:
                    summary = fd_results["summary"]
                    f.write(f"Spacecraft Analyzed: {summary.get('total_spacecraft', 0)}\n")
                    f.write(f"Faults Detected: {summary.get('total_detections', 0)}\n")
                    f.write(f"Success Rate: {summary.get('success_rate', 0.0):.1%}\n")
                    
                    if "confidence_scores" in summary and summary["confidence_scores"]:
                        f.write(f"Average Confidence: {np.mean(summary['confidence_scores']):.3f}\n")
                        f.write(f"Max Confidence: {np.max(summary['confidence_scores']):.3f}\n")
                f.write("\n")
                
                # DRL Decision Summary
                f.write("DRL DECISION MAKING SUMMARY:\n")
                f.write("-" * 30 + "\n")
                drl_results = self.drl_results
                f.write(f"Strategy Selected: {drl_results.get('strategy', 'unknown')}\n")
                f.write(f"Action Taken: {drl_results.get('action', 'unknown')}\n")
                f.write(f"Decision Confidence: {drl_results.get('decision_confidence', 0.0):.3f}\n")
                f.write(f"DRL Model Used: {drl_results.get('drl_model_used', False)}\n")
                f.write(f"Model Type: {drl_results.get('model_type', 'N/A')}\n")
                f.write("\n")
                
                # Task Reassignment Summary
                f.write("TASK REASSIGNMENT SUMMARY:\n")
                f.write("-" * 30 + "\n")
                if hasattr(self, 'reassignment_results'):
                    reassign_results = self.reassignment_results
                    plan = reassign_results.get("reassignment_plan", {})
                    f.write(f"Assignments Created: {len(plan.get('assignments', []))}\n")
                    f.write(f"Execution Success: {reassign_results.get('execution_success', False)}\n")
                    f.write(f"Healthy Spacecraft: {len(reassign_results.get('healthy_spacecraft', []))}\n")
                    f.write(f"Faulty Spacecraft: {len(reassign_results.get('faulty_spacecraft', []))}\n")
                f.write("\n")
                
                # System Health Assessment
                f.write("SYSTEM HEALTH ASSESSMENT:\n")
                f.write("-" * 30 + "\n")
                health = self.integrated_results["system_health"]
                f.write(f"Overall Status: {health['overall_system_status'].upper()}\n")
                f.write(f"Healthy Spacecraft: {health['healthy_spacecraft_count']}\n")
                f.write(f"Faulty Spacecraft: {health['faulty_spacecraft_count']}\n")
                f.write(f"Tasks Reassigned: {health['tasks_reassigned']}\n")
                f.write("\n")
                
                # Detailed Fault Analysis
                f.write("DETAILED FAULT ANALYSIS:\n")
                f.write("-" * 30 + "\n")
                if "detections" in self.fault_detection_results:
                    for spacecraft, detections in self.fault_detection_results["detections"].items():
                        f.write(f"\n{spacecraft}:\n")
                        if detections:
                            for i, detection in enumerate(detections):
                                f.write(f"  Fault {i+1}:\n")
                                f.write(f"    Type: {detection.fault_type}\n")
                                f.write(f"    Time: {detection.detection_time_minutes:.2f} min\n")
                                f.write(f"    Confidence: {detection.confidence:.3f}\n")
                                f.write(f"    Anomaly Score: {detection.anomaly_score:.2f}\n")
                        else:
                            f.write("  No faults detected\n")
                
                # Recommendations
                f.write("\nRECOMMENDATIONS:\n")
                f.write("-" * 30 + "\n")
                self._generate_recommendations(f)
                
        except Exception as e:
            print(f"Error generating comprehensive report: {e}")
    
    def _generate_recommendations(self, file_handle):
        """Generate recommendations based on analysis results"""
        
        health_status = self.integrated_results["system_health"]["overall_system_status"]
        fault_count = self.integrated_results["performance_metrics"]["total_faults_detected"]
        drl_confidence = self.integrated_results["performance_metrics"]["drl_decision_confidence"]
        
        if health_status in ["critical", "failure"]:
            file_handle.write("1. URGENT: Consider emergency procedures\n")
            file_handle.write("2. Evaluate manual override options\n")
            file_handle.write("3. Increase monitoring frequency\n")
        elif health_status == "degraded":
            file_handle.write("1. Monitor system closely\n")
            file_handle.write("2. Prepare backup procedures\n")
            file_handle.write("3. Consider preventive maintenance\n")
        else:
            file_handle.write("1. Continue normal operations\n")
            file_handle.write("2. Maintain routine monitoring\n")
        
        if fault_count > 0:
            file_handle.write(f"4. Investigate root causes of {fault_count} detected faults\n")
            file_handle.write("5. Update fault detection thresholds if needed\n")
        
        if drl_confidence < 0.7:
            file_handle.write("6. Review DRL model decisions\n")
            file_handle.write("7. Consider retraining with additional scenarios\n")
        
        file_handle.write("8. Archive analysis results for future reference\n")
    
    def _generate_mock_fault_results(self) -> Dict:
        """Generate mock fault detection results for testing"""
        
        mock_detections = {
            "Satellite1": [
                type('MockDetection', (), {
                    'fault_detected': True,
                    'fault_type': 'rw_speed_change',
                    'confidence': 0.85,
                    'detection_time_minutes': 15.0,
                    'anomaly_score': 75.5,
                    'affected_component': 'Satellite1'
                })()
            ],
            "Satellite2": [],
            "Satellite3": [],
            "Satellite4": []
        }
        
        return {
            "detections": mock_detections,
            "summary": {
                "total_spacecraft": 4,
                "total_detections": 1,
                "success_rate": 1.0,
                "confidence_scores": [0.85]
            }
        }
    
    def _save_partial_results(self, output_dir: str, error_message: str):
        """Save partial results when analysis fails"""
        
        partial_results = {
            "error": error_message,
            "timestamp": datetime.now().isoformat(),
            "partial_fault_detection": self.fault_detection_results,
            "partial_drl_results": self.drl_results,
            "components_status": {
                "fault_detection": FAULT_DETECTION_AVAILABLE,
                "drl_available": SIMPLE_DRL_AVAILABLE,
                "ml_available": ML_AVAILABLE
            }
        }
        
        try:
            with open(os.path.join(output_dir, "partial_results.json"), 'w') as f:
                json.dump(partial_results, f, indent=2, default=str)
            print(f"Partial results saved to {output_dir}")
        except Exception as e:
            print(f"Could not save partial results: {e}")


def run_integrated_analysis(scenario, scenario_config=None, output_dir="logs", config_profile="development"):
    """
    Main function to run the complete integrated fault detection + DRL analysis
    
    This function should be called from your GUI after Basilisk simulation completes.
    
    Args:
        scenario: Your Basilisk simulation scenario object
        scenario_config: Configuration for your constellation scenario  
        output_dir: Directory to save results
        config_profile: Configuration profile ("development", "production", "research")
    
    Returns:
        Dict containing complete analysis results
    """
    
    print("\nSTARTING INTEGRATED FAULT DETECTION + DRL ANALYSIS")
    print("=" * 60)
    
    # Initialize the integrated system
    system = IntegratedFaultDetectionDRLSystem(config_profile)
    
    # Print configuration summary
    system.config_manager.print_config_summary()
    
    # Run complete analysis
    results = system.run_complete_analysis(scenario, scenario_config, output_dir)
    
    if "error" not in results:
        print("\nANALYSIS COMPLETED SUCCESSFULLY!")
        print(f"Results available in: {output_dir}")
        
        # Print summary
        health = results["system_health"]
        performance = results["performance_metrics"]
        
        print(f"\nSUMMARY:")
        print(f"  System Health: {health['overall_system_status'].upper()}")
        print(f"  Faults Detected: {performance['total_faults_detected']}")
        print(f"  Tasks Reassigned: {health['tasks_reassigned']}")
        print(f"  DRL Confidence: {performance['drl_decision_confidence']:.3f}")
    else:
        print(f"\nANALYSIS FAILED: {results['error']}")
        if results.get("partial_results"):
            print("Partial results were saved")
    
    return results


# Additional utility functions for integration with existing GUI

def setup_drl_integration(bsk_sim_directory: str) -> bool:
    """
    Setup DRL integration by creating necessary directories and checking for files
    
    Args:
        bsk_sim_directory: Path to your BskSim directory
    
    Returns:
        bool: True if setup successful
    """
    
    print("Setting up DRL integration...")
    
    try:
        # Create necessary directories
        os.makedirs(os.path.join(bsk_sim_directory, "logs"), exist_ok=True)
        os.makedirs(os.path.join(bsk_sim_directory, "config"), exist_ok=True)
        
        # Check for simple_task_drl.py
        drl_module_path = os.path.join(bsk_sim_directory, "simple_task_drl.py")
        if not os.path.exists(drl_module_path):
            print(f"Warning: simple_task_drl.py not found at {drl_module_path}")
            print("Please ensure simple_task_drl.py is in the same directory as this file")
            return False
        
        print("DRL integration setup completed successfully")
        return True
        
    except Exception as e:
        print(f"Error setting up DRL integration: {e}")
        return False


def validate_integration_requirements() -> Dict[str, bool]:
    """
    Validate that all required components are available for integration
    
    Returns:
        Dict with component availability status
    """
    
    status = {
        "tensorflow": ML_AVAILABLE,
        "fault_detection": FAULT_DETECTION_AVAILABLE, 
        "drl_task_reassignment": SIMPLE_DRL_AVAILABLE,
        "basilisk": True  # Assume available since we're in Basilisk environment
    }
    
    print("Integration Requirements Status:")
    for component, available in status.items():
        symbol = "AVAILABLE" if available else "UNAVAILABLE"
        print(f"  {symbol}: {component.replace('_', ' ').title()}")
    
    return status


if __name__ == "__main__":
    print("Integrated Fault Detection + DRL Task Reassignment System")
    print("=" * 60)
    
    # Validate requirements
    requirements = validate_integration_requirements()
    
    if all(requirements.values()):
        print("\nAll requirements satisfied - ready for integration!")
    else:
        print("\nSome requirements missing:")
        for comp, avail in requirements.items():
            if not avail:
                print(f"  - {comp.replace('_', ' ').title()}")
        print("\nInstall missing components:")
        if not requirements["tensorflow"]:
            print("  pip install tensorflow")
    
    print("\nUsage:")
    print("1. Ensure simple_task_drl.py is in the same directory")
    print("2. Run your Basilisk constellation simulation")
    print("3. Call: run_integrated_analysis(scenario, config, output_dir)")
    print("4. View comprehensive results in output directory")