#!/usr/bin/env python
"""
debug_detector_selection.py

Diagnostic script to help debug which fault detector is being used.
Run this to check your current setup.
"""

import sys
import os


def check_imports():
    """Check if required modules can be imported"""
    print("\n" + "="*70)
    print("CHECKING IMPORTS")
    print("="*70)
    
    # Check LSTM detector
    try:
        from real_ml_fault_detection import RealMLFaultDetector
        print("✓ LSTM Detector (real_ml_fault_detection.py): AVAILABLE")
        lstm_available = True
    except ImportError as e:
        print(f"✗ LSTM Detector: NOT AVAILABLE")
        print(f"  Error: {e}")
        lstm_available = False
    
    # Check Isolation Forest detector
    try:
        from isolation_forest_fault_detection import SatelliteIsolationForestDetector
        print("✓ Isolation Forest Detector (isolation_forest_fault_detection.py): AVAILABLE")
        if_available = True
    except ImportError as e:
        print(f"✗ Isolation Forest Detector: NOT AVAILABLE")
        print(f"  Error: {e}")
        if_available = False
    
    # Check Router
    try:
        from fault_detection_router import run_fault_detection_on_scenario
        print("✓ Fault Detection Router (fault_detection_router.py): AVAILABLE")
        router_available = True
    except ImportError as e:
        print(f"✗ Fault Detection Router: NOT AVAILABLE")
        print(f"  Error: {e}")
        router_available = False
    
    # Check scikit-learn
    try:
        from sklearn.ensemble import IsolationForest
        print("✓ scikit-learn (Isolation Forest dependency): AVAILABLE")
        sklearn_available = True
    except ImportError:
        print("✗ scikit-learn: NOT AVAILABLE (install with: pip install scikit-learn)")
        sklearn_available = False
    
    # Check TensorFlow
    try:
        import tensorflow as tf
        print(f"✓ TensorFlow (LSTM dependency): AVAILABLE (version {tf.__version__})")
        tf_available = True
    except ImportError:
        print("✗ TensorFlow: NOT AVAILABLE (install with: pip install tensorflow)")
        tf_available = False
    
    return {
        'lstm': lstm_available,
        'isolation_forest': if_available,
        'router': router_available,
        'sklearn': sklearn_available,
        'tensorflow': tf_available
    }


def check_file_locations():
    """Check if required files exist in the current directory"""
    print("\n" + "="*70)
    print("CHECKING FILE LOCATIONS")
    print("="*70)
    
    required_files = [
        'real_ml_fault_detection.py',
        'isolation_forest_fault_detection.py',
        'fault_detection_router.py',
        'fault_detection_tab.py'
    ]
    
    current_dir = os.getcwd()
    print(f"Current directory: {current_dir}\n")
    
    all_found = True
    for filename in required_files:
        if os.path.exists(filename):
            print(f"✓ {filename}: FOUND")
        else:
            print(f"✗ {filename}: NOT FOUND")
            all_found = False
    
    return all_found


def check_model_files():
    """Check if model files exist"""
    print("\n" + "="*70)
    print("CHECKING MODEL FILES")
    print("="*70)
    
    model_files = [
        'anomaly_detection_model.keras',
        'anomaly_detection_model.h5'
    ]
    
    found_model = None
    for filename in model_files:
        if os.path.exists(filename):
            print(f"✓ {filename}: FOUND")
            found_model = filename
        else:
            print(f"  {filename}: not found")
    
    if found_model:
        print(f"\nLSTM model available: {found_model}")
    else:
        print("\n⚠ No LSTM model file found (needed only for LSTM/Autoencoder mode)")
    
    return found_model


def analyze_code_usage():
    """Check where detection is being called in the code"""
    print("\n" + "="*70)
    print("CHECKING CODE INTEGRATION")
    print("="*70)
    
    print("\nSearching for detection calls...")
    
    # Common patterns to search for
    patterns = [
        'run_real_ml_detection_on_scenario',
        'integrate_real_ml_with_basilisk',
        'run_isolation_forest_detection',
        'integrate_isolation_forest_with_basilisk',
        'run_fault_detection_on_scenario'
    ]
    
    # Files to check
    check_files = []
    for root, dirs, files in os.walk('.'):
        # Skip hidden directories and common non-code directories
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ['__pycache__', 'venv', 'env']]
        
        for file in files:
            if file.endswith('.py'):
                filepath = os.path.join(root, file)
                check_files.append(filepath)
    
    findings = {}
    for filepath in check_files[:20]:  # Limit to first 20 files
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
                
            for pattern in patterns:
                if pattern in content:
                    if filepath not in findings:
                        findings[filepath] = []
                    findings[filepath].append(pattern)
        except:
            pass
    
    if findings:
        print("\nFound detection calls in:")
        for filepath, patterns_found in findings.items():
            print(f"\n  {filepath}:")
            for pattern in patterns_found:
                print(f"    - {pattern}")
    else:
        print("\n⚠ No detection function calls found in Python files")
        print("  This might be why detection isn't working!")


def provide_integration_advice(availability):
    """Provide advice on how to integrate based on what's available"""
    print("\n" + "="*70)
    print("INTEGRATION ADVICE")
    print("="*70)
    
    if not availability['router']:
        print("\n⚠ MISSING: fault_detection_router.py")
        print("  This is the key file that routes between LSTM and Isolation Forest!")
        print("\n  Solution:")
        print("  1. Make sure fault_detection_router.py is in your project directory")
        print("  2. In your simulation code, replace:")
        print("     from real_ml_fault_detection import run_real_ml_detection_on_scenario")
        print("     with:")
        print("     from fault_detection_router import run_fault_detection_on_scenario")
        print("\n  3. Then call:")
        print("     results = run_fault_detection_on_scenario(scenario, fault_detection_tab)")
    
    elif availability['router']:
        print("\n✓ Router is available!")
        print("\n  Make sure your simulation code uses:")
        print("  ```python")
        print("  from fault_detection_router import run_fault_detection_on_scenario")
        print("  results = run_fault_detection_on_scenario(")
        print("      scenario=scenario,")
        print("      fault_detection_tab=fault_detection_tab,")
        print("      scenario_config=scenario_config,")
        print("      output_dir='./results'")
        print("  )")
        print("  ```")
    
    if not availability['isolation_forest']:
        print("\n⚠ MISSING: isolation_forest_fault_detection.py")
        print("  Solution: Add this file to your project directory")
    
    if not availability['sklearn']:
        print("\n⚠ MISSING: scikit-learn")
        print("  Solution: pip install scikit-learn")
    
    if availability['isolation_forest'] and availability['sklearn'] and availability['router']:
        print("\n✓ All Isolation Forest components are available!")
        print("\n  To use Isolation Forest:")
        print("  1. In GUI: Select 'Isolation Forest (Unsupervised)' from model dropdown")
        print("  2. Click 'Load Selected Model'")
        print("  3. Click 'Start Detection'")
        print("  4. Run your simulation")
        print("  5. The router will automatically use Isolation Forest")


def main():
    """Run all diagnostic checks"""
    print("\n" + "="*70)
    print("FAULT DETECTION DIAGNOSTIC TOOL")
    print("="*70)
    
    # Check imports
    availability = check_imports()
    
    # Check files
    files_ok = check_file_locations()
    
    # Check models
    model_found = check_model_files()
    
    # Analyze code
    analyze_code_usage()
    
    # Provide advice
    provide_integration_advice(availability)
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    issues = []
    if not availability['router']:
        issues.append("Missing fault_detection_router.py - THIS IS CRITICAL!")
    if not availability['isolation_forest']:
        issues.append("Missing isolation_forest_fault_detection.py")
    if not availability['sklearn']:
        issues.append("Missing scikit-learn library")
    if not availability['lstm']:
        issues.append("Missing real_ml_fault_detection.py")
    if not availability['tensorflow']:
        issues.append("Missing TensorFlow library")
    
    if issues:
        print("\n⚠ Issues found:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
        
        print("\n🔧 MOST LIKELY PROBLEM:")
        print("  Your simulation code is probably calling the old LSTM function directly")
        print("  instead of using the router that selects between LSTM and Isolation Forest.")
        print("\n  Find where your simulation calls fault detection and update it to use:")
        print("  fault_detection_router.run_fault_detection_on_scenario(...)")
    else:
        print("\n✓ All components available!")
        print("  If Isolation Forest still isn't working, check your simulation code")
        print("  to ensure it's calling the router function.")
    
    print("\n" + "="*70 + "\n")


if __name__ == "__main__":
    main()