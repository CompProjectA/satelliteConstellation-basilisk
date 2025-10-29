#!/usr/bin/env python
"""
fault_detection_router.py

Routes fault detection to the appropriate detector (LSTM or Isolation Forest)
based on the model type selected in the GUI.

This ensures that when simulation runs, the correct detector is used.
"""

import os
import sys
import numpy as np
from typing import Dict, Optional


def run_fault_detection_on_scenario(scenario, fault_detection_tab, scenario_config=None, output_dir="."):
    """
    Main router function that selects and runs the appropriate fault detector.
    
    This should be called from your simulation code after the Basilisk scenario completes.
    
    Args:
        scenario: Basilisk scenario object
        fault_detection_tab: Reference to the FaultDetectionTab instance from GUI
        scenario_config: Optional scenario configuration
        output_dir: Directory to save results
        
    Returns:
        Detection results dictionary
    """
    
    # Get the selected model type from GUI
    model_type = fault_detection_tab.model_type_var.get()
    
    print(f"\n{'='*70}")
    print(f"FAULT DETECTION ROUTER")
    print(f"{'='*70}")
    print(f"Selected Model Type: {model_type}")
    
    # Route to appropriate detector
    if model_type == "Isolation":
        print("Routing to: ISOLATION FOREST")
        return run_isolation_forest_detection(scenario, fault_detection_tab, scenario_config, output_dir)
    
    elif model_type in ["Autoencoder", "LIME", "SHAP", "Ensemble"]:
        print(f"Routing to: LSTM/AUTOENCODER ({model_type})")
        return run_lstm_detection(scenario, fault_detection_tab, scenario_config, output_dir)
    
    else:
        print(f"WARNING: Unknown model type '{model_type}', defaulting to LSTM")
        return run_lstm_detection(scenario, fault_detection_tab, scenario_config, output_dir)


def run_isolation_forest_detection(scenario, fault_detection_tab, scenario_config=None, output_dir="."):
    """
    Run Isolation Forest detection on the scenario.
    """
    print("\n" + "="*70)
    print("RUNNING ISOLATION FOREST DETECTION")
    print("="*70)
    
    # Check if Isolation Forest detector is initialized
    if not hasattr(fault_detection_tab, 'isolation_forest_detector') or fault_detection_tab.isolation_forest_detector is None:
        error_msg = "ERROR: Isolation Forest detector not initialized!\n"
        error_msg += "Please load the Isolation Forest model in the GUI first."
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        return {"error": "Isolation Forest not initialized"}
    
    try:
        # Import the Isolation Forest detection function
        from isolation_forest_fault_detection import integrate_isolation_forest_with_basilisk
        
        fault_detection_tab.parent_app.add_log("Starting Isolation Forest detection...")
        
        # Run Isolation Forest detection
        results = integrate_isolation_forest_with_basilisk(scenario, scenario_config)
        
        if "error" in results:
            print(f"Isolation Forest detection failed: {results['error']}")
            fault_detection_tab.parent_app.add_log(f"Isolation Forest failed: {results['error']}")
            return results
        
        # Save results
        from isolation_forest_fault_detection import save_isolation_forest_results
        save_isolation_forest_results(results, output_dir)
        
        # Log success
        summary = results.get('summary', {})
        fault_detection_tab.parent_app.add_log(
            f"Isolation Forest complete: {summary.get('total_detections', 0)} anomalies detected"
        )
        
        print(f"\n{'='*70}")
        print(f"ISOLATION FOREST DETECTION COMPLETE")
        print(f"  Spacecraft: {summary.get('total_spacecraft', 0)}")
        print(f"  Anomalies: {summary.get('total_detections', 0)}")
        print(f"  Success Rate: {summary.get('success_rate', 0):.1%}")
        print(f"{'='*70}\n")
        
        # Display results in GUI
        fault_detection_tab.display_ml_results(results)
        
        return results
        
    except ImportError as e:
        error_msg = f"Could not import Isolation Forest module: {e}\n"
        error_msg += "Make sure isolation_forest_fault_detection.py is in the project directory."
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        return {"error": str(e)}
        
    except Exception as e:
        error_msg = f"Error during Isolation Forest detection: {e}"
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        import traceback
        traceback.print_exc()
        return {"error": str(e)}


def run_lstm_detection(scenario, fault_detection_tab, scenario_config=None, output_dir="."):
    """
    Run LSTM/Autoencoder detection on the scenario.
    """
    print("\n" + "="*70)
    print("RUNNING LSTM/AUTOENCODER DETECTION")
    print("="*70)
    
    # Check if LSTM detector is initialized
    if not hasattr(fault_detection_tab, 'ml_detector') or fault_detection_tab.ml_detector is None:
        error_msg = "ERROR: LSTM detector not initialized!\n"
        error_msg += "Please load the LSTM model in the GUI first."
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        return {"error": "LSTM not initialized"}
    
    if not fault_detection_tab.ml_detector.is_loaded:
        error_msg = "ERROR: LSTM model not loaded properly!"
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        return {"error": "LSTM model not loaded"}
    
    try:
        # Import the LSTM detection function
        from real_ml_fault_detection import (
            integrate_real_ml_with_basilisk, 
            save_real_detection_results
        )
        
        fault_detection_tab.parent_app.add_log("Starting LSTM detection...")
        
        # Run LSTM detection
        results = integrate_real_ml_with_basilisk(scenario, scenario_config)
        
        if "error" in results:
            print(f"LSTM detection failed: {results['error']}")
            fault_detection_tab.parent_app.add_log(f"LSTM failed: {results['error']}")
            return results
        
        # Save results
        save_real_detection_results(results, output_dir)
        
        # Log success
        summary = results.get('summary', {})
        fault_detection_tab.parent_app.add_log(
            f"LSTM detection complete: {summary.get('total_detections', 0)} faults detected"
        )
        
        print(f"\n{'='*70}")
        print(f"LSTM DETECTION COMPLETE")
        print(f"  Spacecraft: {summary.get('total_spacecraft', 0)}")
        print(f"  Detections: {summary.get('total_detections', 0)}")
        print(f"  Success Rate: {summary.get('success_rate', 0):.1%}")
        print(f"{'='*70}\n")
        
        # Display results in GUI
        fault_detection_tab.display_ml_results(results)
        
        return results
        
    except ImportError as e:
        error_msg = f"Could not import LSTM module: {e}\n"
        error_msg += "Make sure real_ml_fault_detection.py is in the project directory."
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        return {"error": str(e)}
        
    except Exception as e:
        error_msg = f"Error during LSTM detection: {e}"
        print(error_msg)
        fault_detection_tab.parent_app.add_log(error_msg)
        import traceback
        traceback.print_exc()
        return {"error": str(e)}


def get_active_detector_info(fault_detection_tab):
    """
    Get information about which detector is currently active.
    
    Returns:
        Dictionary with detector information
    """
    model_type = fault_detection_tab.model_type_var.get()
    
    info = {
        "model_type": model_type,
        "detector_ready": False,
        "detector_name": "None"
    }
    
    if model_type == "Isolation":
        if hasattr(fault_detection_tab, 'isolation_forest_detector') and fault_detection_tab.isolation_forest_detector is not None:
            info["detector_ready"] = True
            info["detector_name"] = "Isolation Forest"
            info["contamination"] = fault_detection_tab.if_contamination_var.get()
            info["n_estimators"] = fault_detection_tab.if_estimators_var.get()
    
    elif model_type in ["Autoencoder", "LIME", "SHAP", "Ensemble"]:
        if hasattr(fault_detection_tab, 'ml_detector') and fault_detection_tab.ml_detector is not None:
            if fault_detection_tab.ml_detector.is_loaded:
                info["detector_ready"] = True
                info["detector_name"] = f"LSTM {model_type}"
                info["model_path"] = fault_detection_tab.ml_detector.model_path
    
    return info


def print_detector_status(fault_detection_tab):
    """
    Print current detector status for debugging.
    """
    info = get_active_detector_info(fault_detection_tab)
    
    print("\n" + "="*70)
    print("DETECTOR STATUS")
    print("="*70)
    print(f"Model Type Selected: {info['model_type']}")
    print(f"Detector Name: {info['detector_name']}")
    print(f"Detector Ready: {info['detector_ready']}")
    
    if info['model_type'] == "Isolation":
        if 'contamination' in info:
            print(f"Contamination: {info['contamination']*100:.1f}%")
            print(f"N Estimators: {info['n_estimators']}")
    else:
        if 'model_path' in info:
            print(f"Model Path: {info['model_path']}")
    
    print("="*70 + "\n")


if __name__ == "__main__":
    print("Fault Detection Router")
    print("=" * 70)
    print("\nThis module routes fault detection to the correct detector.")
    print("\nUsage in your simulation code:")
    print("  from fault_detection_router import run_fault_detection_on_scenario")
    print("  results = run_fault_detection_on_scenario(scenario, fault_detection_tab)")
    print("\n" + "="*70)