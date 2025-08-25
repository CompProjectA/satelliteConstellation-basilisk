#!/usr/bin/env python
"""
gui_ml_integration.py

Shows how to integrate real ML fault detection with your existing GUI system.
Add this to your spacecraft_simulation.py to get real ML detection after simulation.
"""

# Add this to your spacecraft_simulation.py imports:
try:
    from real_ml_fault_detection import run_real_ml_detection_on_scenario
    REAL_ML_AVAILABLE = True
    print("✅ Real ML fault detection available")
except ImportError as e:
    REAL_ML_AVAILABLE = False
    print(f"⚠️ Real ML fault detection not available: {e}")


def add_real_ml_detection_to_simulation():
    """
    This shows exactly where to add real ML detection to your existing simulation
    """
    
    # YOUR EXISTING CODE - run_custom_simulation function
    # ... all your existing simulation code ...
    
    # After this line in your run_custom_simulation():
    # scSim.ExecuteSimulation()
    
    # ADD THIS SECTION:
    
    print("\n" + "="*60)
    print("SPRINT 4: REAL ML FAULT DETECTION")
    print("="*60)
    
    if REAL_ML_AVAILABLE:
        try:
            print("🤖 Running client's ML model on REAL Basilisk data...")
            
            # Run real ML detection on the completed simulation
            ml_results = run_real_ml_detection_on_scenario(
                scenario=scenario,           # Your simulation scenario object
                scenario_config=config,      # Your SimulationConfig object  
                output_dir=output_dir        # Your output directory
            )
            
            if ml_results:
                print("✅ Real ML fault detection completed successfully!")
                
                # Add ML results to your figure list for the GUI
                if 'figureList' in locals():
                    figureList['Real_ML_Detection_Summary'] = create_ml_detection_plot(ml_results)
                
                # Update your return to include ML results
                return scenario, viz, figureList, output_dir, ml_results
            else:
                print("⚠️ Real ML detection had issues")
                return scenario, viz, figureList, output_dir, None
                
        except Exception as e:
            print(f"❌ Real ML detection error: {e}")
            return scenario, viz, figureList, output_dir, None
    else:
        print("⚠️ Real ML detection not available - check TensorFlow installation")
        return scenario, viz, figureList, output_dir, None


def create_ml_detection_plot(ml_results):
    """
    Create a plot showing real ML detection results
    """
    import matplotlib.pyplot as plt
    
    try:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot 1: Detection timeline
        detection_times = ml_results['summary']['detection_times']
        confidence_scores = ml_results['summary']['confidence_scores']
        
        if detection_times:
            ax1.scatter(detection_times, confidence_scores, c='red', s=100, alpha=0.7, label='ML Detections')
            ax1.axhline(y=0.5, color='orange', linestyle='--', alpha=0.7, label='Detection Threshold')
            ax1.set_xlabel('Time (minutes)')
            ax1.set_ylabel('ML Confidence Score')
            ax1.set_title('Real ML Fault Detection Results\n(Client\'s Model on Basilisk Data)')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
        else:
            ax1.text(0.5, 0.5, 'No ML Detections Found', ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title('Real ML Fault Detection Results - No Faults Detected')
        
        # Plot 2: Per-spacecraft summary
        spacecraft_names = list(ml_results['detections'].keys())
        detection_counts = [len(detections) for detections in ml_results['detections'].values()]
        
        if spacecraft_names:
            bars = ax2.bar(spacecraft_names, detection_counts, color=['red' if count > 0 else 'green' for count in detection_counts])
            ax2.set_xlabel('Spacecraft')
            ax2.set_ylabel('Number of ML Detections')
            ax2.set_title('ML Detections per Spacecraft')
            
            # Add count labels on bars
            for bar, count in zip(bars, detection_counts):
                if count > 0:
                    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                            str(count), ha='center', va='bottom')
        
        plt.tight_layout()
        return fig
        
    except Exception as e:
        print(f"❌ Error creating ML detection plot: {e}")
        return None


# EXACT INTEGRATION INSTRUCTIONS:

def integrate_with_spacecraft_simulation():
    """
    STEP-BY-STEP: How to add this to your spacecraft_simulation.py
    """
    
    instructions = """
    
    🔧 INTEGRATION STEPS FOR spacecraft_simulation.py:
    
    1. ADD IMPORT at top of spacecraft_simulation.py:
    
    ```python
    # Add this after your existing imports
    try:
        from real_ml_fault_detection import run_real_ml_detection_on_scenario
        REAL_ML_AVAILABLE = True
        print("✅ Real ML fault detection available")
    except ImportError as e:
        REAL_ML_AVAILABLE = False
        print(f"⚠️ Real ML fault detection not available: {e}")
    ```
    
    2. MODIFY run_custom_simulation() function:
    
    Find this section in your run_custom_simulation():
    ```python
    # Execute the simulation
    scSim.ExecuteSimulation()
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    ```
    
    3. ADD THIS CODE right after scSim.ExecuteSimulation():
    
    ```python
    # Execute the simulation
    scSim.ExecuteSimulation()
    
    # ===== ADD THIS SECTION =====
    print("\\n" + "="*60)
    print("SPRINT 4: REAL ML FAULT DETECTION")
    print("="*60)
    
    ml_results = None
    if REAL_ML_AVAILABLE:
        try:
            print("🤖 Running client's ML model on REAL Basilisk data...")
            
            ml_results = run_real_ml_detection_on_scenario(
                scenario=scenario,
                scenario_config=config,
                output_dir=output_dir
            )
            
            if ml_results:
                print("✅ REAL ML FAULT DETECTION COMPLETED!")
                summary = ml_results['summary']
                print(f"   Spacecraft: {summary['total_spacecraft']}")
                print(f"   ML Detections: {summary['total_detections']}")
                if summary['detection_times']:
                    print(f"   First Detection: {min(summary['detection_times']):.1f} min")
            else:
                print("⚠️ Real ML detection had issues")
                
        except Exception as e:
            print(f"❌ Real ML detection error: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("⚠️ Real ML detection not available")
        print("   Install TensorFlow: pip install tensorflow")
        print("   Copy client's model: anomaly_detection_model.keras")
    # ===== END ADDITION =====
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    ```
    
    4. UPDATE your return statement at the end of run_custom_simulation():
    
    Change from:
    ```python
    return scenario, viz, figureList, output_dir
    ```
    
    To:
    ```python
    return scenario, viz, figureList, output_dir, ml_results
    ```
    
    5. UPDATE your GUI code to handle ML results:
    
    In spacecraft_simulator_gui.py, modify the simulation call:
    ```python
    # Change from:
    scenario, viz, plots, output_dir = run_custom_simulation(config)
    
    # To:
    result = run_custom_simulation(config)
    if len(result) == 5:
        scenario, viz, plots, output_dir, ml_results = result
    else:
        scenario, viz, plots, output_dir = result
        ml_results = None
    
    # Add ML results display
    if ml_results:
        self.display_ml_results(ml_results)
    ```
    
    """
    
    return instructions


def display_ml_results_in_gui(ml_results):
    """
    Example of how to display ML results in your GUI
    """
    
    gui_integration_code = """
    
    def display_ml_results(self, ml_results):
        '''Add this method to your GUI class'''
        
        if not ml_results:
            return
        
        summary = ml_results['summary']
        
        # Update status in GUI
        status_text = f"ML Detection: {summary['total_detections']} faults found"
        if hasattr(self, 'status_label'):
            self.status_label.setText(status_text)
        
        # Add to results tab
        if hasattr(self, 'results_text'):
            ml_summary = f'''
    REAL ML FAULT DETECTION RESULTS:
    ================================
    Spacecraft Analyzed: {summary['total_spacecraft']}
    Total ML Detections: {summary['total_detections']}
    Success Rate: {summary['success_rate']:.1%}
    
    DETECTIONS BY SPACECRAFT:
    '''
            
            for spacecraft, detections in ml_results['detections'].items():
                ml_summary += f"\\n{spacecraft}: {len(detections)} detections"
                for detection in detections[:3]:  # Show first 3
                    ml_summary += f"\\n  - {detection.detection_time_minutes:.1f}min (conf: {detection.confidence:.3f})"
            
            self.results_text.append(ml_summary)
        
        print("✅ ML results displayed in GUI")
    
    """
    
    return gui_integration_code


if __name__ == "__main__":
    print("🔧 GUI Integration Helper for Real ML Detection")
    print("="*60)
    
    print("This file shows how to integrate real ML fault detection")
    print("with your existing Basilisk simulation GUI system.")
    print()
    
    print("Key Integration Points:")
    print("1. ✅ Extract REAL telemetry from Basilisk simulation")
    print("2. ✅ Feed to client's ML model (anomaly_detection_model.keras)")
    print("3. ✅ Get ML predictions on REAL spacecraft data")
    print("4. ✅ Display results in GUI")
    print()
    
    instructions = integrate_with_spacecraft_simulation()
    print(instructions)
    
    print("\n🎯 Sprint 4 Goal: REAL ML detection on REAL Basilisk data ✅")
    print("This replaces synthetic detection with actual client ML model!")