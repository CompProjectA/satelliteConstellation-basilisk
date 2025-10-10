"""
Satellite Fault Detection with LIME & SHAP XAI Integration

This adds explainability to your existing autoencoder-based fault detection system.
Works with your real_ml_fault_detection.py architecture.
"""

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import os  # Added for file operations

# XAI Libraries
import shap
from lime import lime_tabular

class SatelliteXAIExplainer:
    """
    XAI wrapper for your satellite fault detection model.
    Provides LIME and SHAP explanations for anomaly predictions.
    """
    
    def __init__(self, model, feature_names=None):
        """
        Initialize XAI explainer with your trained model.
        
        Args:
            model: Your trained TensorFlow/Keras autoencoder model
            feature_names: List of feature names (17 features in your case)
        """
        self.model = model
        self.feature_names = feature_names or self._default_feature_names()
        
        # Initialize LIME explainer
        self.lime_explainer = None
        
        # Initialize SHAP explainer
        self.shap_explainer = None
        self.shap_background = None
        
        print("XAI Explainer initialized")
        print(f"  Features: {len(self.feature_names)}")
    
    def _default_feature_names(self):
        """Default feature names for your 17-feature model"""
        return [
            'RW1_speed', 'RW2_speed', 'RW3_speed', 'RW4_speed',
            'RW1_speed_deriv', 'RW2_speed_deriv', 'RW3_speed_deriv', 'RW4_speed_deriv',
            'RW1_torque', 'RW2_torque', 'RW3_torque', 'RW4_torque',
            'attitude_error', 'attitude_deriv',
            'time_norm', 'time_sin', 'time_cos'
        ]
    
    def _prediction_wrapper(self, X):
        """
        Wrapper function that converts model output to anomaly scores.
        Both LIME and SHAP need a function that returns scalar scores.
        
        For autoencoder: returns reconstruction error as anomaly score
        """
        # Reshape if needed (LIME/SHAP flatten the input)
        if len(X.shape) == 2:
            # X is (batch_size, 170) - need to reshape to (batch_size, 10, 17)
            X = X.reshape(-1, 10, 17)
        
        # Get reconstruction from autoencoder
        reconstructions = self.model.predict(X, verbose=0)
        
        # Calculate reconstruction error (MSE per sample)
        errors = np.mean(np.square(X - reconstructions), axis=(1, 2))
        
        # Return as (n_samples, 1) for compatibility
        return errors.reshape(-1, 1)
    
    def prepare_lime_explainer(self, training_data, mode='tabular'):
        """
        Initialize LIME explainer with background training data.
        
        Args:
            training_data: numpy array of shape (n_samples, 10, 17)
            mode: 'tabular' (flatten time series) or 'timeseries' (preserve sequence)
        """
        print("Preparing LIME explainer...")
        
        # Flatten training data for LIME (it works on flat features)
        # (n_samples, 10, 17) -> (n_samples, 170)
        training_flat = training_data.reshape(training_data.shape[0], -1)
        
        # Create feature names for flattened data
        flat_feature_names = []
        for t in range(10):
            for f in self.feature_names:
                flat_feature_names.append(f"{f}_t{t}")
        
        # Initialize LIME tabular explainer
        self.lime_explainer = lime_tabular.LimeTabularExplainer(
            training_data=training_flat,
            feature_names=flat_feature_names,
            mode='regression',  # Predicting continuous anomaly score
            verbose=False
        )
        
        print(f"  LIME explainer ready with {len(flat_feature_names)} features")
        return self.lime_explainer
    
    def prepare_shap_explainer(self, background_data, explainer_type='kernel'):
        """
        Initialize SHAP explainer with background data.
        
        Args:
            background_data: numpy array of shape (n_samples, 10, 17)
            explainer_type: 'deep' (for neural networks) or 'kernel' (model-agnostic)
        """
        print(f"Preparing SHAP explainer (type: {explainer_type})...")
        
        # Store background data
        self.shap_background = background_data[:100]  # Use subset for efficiency
        
        if explainer_type == 'deep':
            # DeepExplainer for neural networks (faster)
            # This needs special handling for the model architecture
            print("  Warning: DeepExplainer requires specific model architecture")
            print("  Falling back to KernelExplainer for reliability")
            explainer_type = 'kernel'
        
        if explainer_type == 'kernel':
            # KernelExplainer (model-agnostic, slower but more general)
            # CRITICAL: Flatten background data for SHAP
            background_flat = self.shap_background.reshape(self.shap_background.shape[0], -1)
            print(f"  Background data shape: {self.shap_background.shape} -> {background_flat.shape}")
            
            def predict_flat(X_flat):
                """Wrapper that handles flattened input"""
                # Reshape back to (n_samples, 10, 17) for model prediction
                X_reshaped = X_flat.reshape(-1, 10, 17)
                # Get anomaly scores
                return self._prediction_wrapper(X_reshaped).flatten()
            
            self.shap_explainer = shap.KernelExplainer(
                predict_flat,
                background_flat
            )
        
        print(f"  SHAP explainer ready")
        return self.shap_explainer
    
    def explain_with_lime(self, instance, num_features=10):
        """
        Generate LIME explanation for a single instance.
        
        Args:
            instance: numpy array of shape (10, 17) - single time series
            num_features: number of top features to show
            
        Returns:
            explanation object with visualization methods
        """
        if self.lime_explainer is None:
            raise ValueError("LIME explainer not initialized. Call prepare_lime_explainer() first.")
        
        print(f"Generating LIME explanation...")
        
        # Flatten instance
        instance_flat = instance.reshape(-1)
        
        # Generate explanation
        explanation = self.lime_explainer.explain_instance(
            instance_flat,
            self._prediction_wrapper,
            num_features=num_features
        )
        
        return explanation
    
    def explain_with_shap(self, instances, max_samples=10):
        """
        Generate SHAP explanation for instances.
        
        Args:
            instances: numpy array of shape (n_samples, 10, 17)
            max_samples: maximum samples to explain (SHAP can be slow)
            
        Returns:
            shap_values: SHAP values for visualization
        """
        if self.shap_explainer is None:
            raise ValueError("SHAP explainer not initialized. Call prepare_shap_explainer() first.")
        
        print(f"Generating SHAP explanations for {min(len(instances), max_samples)} samples...")
        
        # Limit samples for efficiency
        instances_subset = instances[:max_samples]
        
        # CRITICAL FIX: Flatten instances for SHAP
        # SHAP expects 2D input: (n_samples, features)
        # We have: (n_samples, 10, 17) -> (n_samples, 170)
        print(f"  Original shape: {instances_subset.shape}")
        instances_flat = instances_subset.reshape(instances_subset.shape[0], -1)
        print(f"  Flattened shape for SHAP: {instances_flat.shape}")
        
        # Calculate SHAP values
        shap_values = self.shap_explainer.shap_values(instances_flat)
        
        print(f"  SHAP values shape: {shap_values.shape if hasattr(shap_values, 'shape') else type(shap_values)}")
        
        return shap_values
    
    def visualize_lime_explanation(self, explanation, save_path=None):
        """Visualize LIME explanation"""
        fig = explanation.as_pyplot_figure()
        plt.title("LIME: Feature Importance for Anomaly Detection")
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  LIME plot saved: {save_path}")
        
        # Don't call plt.show() - causes warnings in non-interactive environments
        # The figure will be saved to file
        plt.close(fig)
        return fig
    
    def visualize_shap_summary(self, shap_values, instances, save_path=None):
        """Visualize SHAP summary plot"""
        # Reshape for visualization
        instances_flat = instances.reshape(instances.shape[0], -1)
        shap_flat = shap_values.reshape(shap_values.shape[0], -1) if len(shap_values.shape) > 2 else shap_values
        
        plt.figure(figsize=(12, 8))
        shap.summary_plot(shap_flat, instances_flat, feature_names=self._get_flat_feature_names(), show=False)
        plt.title("SHAP: Feature Importance Summary")
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  SHAP summary saved: {save_path}")
        
        # Close to avoid warnings
        plt.close()
    
    def visualize_shap_waterfall(self, shap_values, instance_idx=0, save_path=None):
        """Visualize SHAP waterfall plot for a single instance"""
        instance_shap = shap_values[instance_idx]
        
        plt.figure(figsize=(10, 8))
        
        # Create SHAP Explanation object for waterfall plot
        try:
            shap.waterfall_plot(
                shap.Explanation(
                    values=instance_shap.flatten(),
                    base_values=np.mean(shap_values),
                    feature_names=self._get_flat_feature_names()
                ),
                show=False
            )
        except Exception as e:
            print(f"  Warning: Could not create waterfall plot: {e}")
            print(f"  Creating bar plot instead")
            # Fallback to bar plot
            feature_importance = instance_shap.flatten()
            top_indices = np.argsort(np.abs(feature_importance))[-10:]  # Top 10
            
            plt.barh(range(len(top_indices)), feature_importance[top_indices])
            plt.yticks(range(len(top_indices)), [self._get_flat_feature_names()[i] for i in top_indices])
            plt.xlabel('SHAP Value')
            plt.title(f'Top 10 Feature Importance (Instance {instance_idx})')
        
        plt.title(f"SHAP Feature Importance: Instance {instance_idx}")
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  SHAP waterfall saved: {save_path}")
        
        # Close to avoid warnings
        plt.close()
    
    def _get_flat_feature_names(self):
        """Get flattened feature names for all time steps"""
        flat_names = []
        for t in range(10):
            for f in self.feature_names:
                flat_names.append(f"{f}_t{t}")
        return flat_names


def integrate_xai_with_fault_detection(model_path, telemetry_data, anomaly_instances=None, output_dir=None):
    """
    Main integration function: Add XAI to your existing fault detection.
    
    Args:
        model_path: Path to your trained model (.keras file)
        telemetry_data: Training/background data (n_samples, 10, 17)
        anomaly_instances: Specific anomalous instances to explain (n_anomalies, 10, 17)
        output_dir: Directory to save XAI plots (default: current directory/xai_plots)
        
    Returns:
        Dictionary with XAI results and visualizations
    """
    print("="*70)
    print("INTEGRATING XAI (LIME & SHAP) WITH FAULT DETECTION")
    print("="*70)
    
    # CRITICAL FIX: Ensure output_dir is never None
    if output_dir is None:
        output_dir = os.getcwd()
        print(f"Warning: output_dir was None, using current directory: {output_dir}")
    
    # Create XAI output directory
    xai_output_dir = os.path.join(output_dir, "xai_plots")
    os.makedirs(xai_output_dir, exist_ok=True)
    
    print(f"\nXAI plots will be saved to: {xai_output_dir}")
    
    # Load your trained model
    print("\n1. Loading trained model...")
    model = tf.keras.models.load_model(model_path)
    print(f"   Model loaded: {model_path}")
    print(f"   Input shape: {model.input_shape}")
    print(f"   Output shape: {model.output_shape}")
    
    # Initialize XAI explainer
    print("\n2. Initializing XAI explainer...")
    xai_explainer = SatelliteXAIExplainer(model)
    
    # Prepare LIME
    print("\n3. Preparing LIME explainer...")
    xai_explainer.prepare_lime_explainer(telemetry_data)
    
    # Prepare SHAP
    print("\n4. Preparing SHAP explainer...")
    xai_explainer.prepare_shap_explainer(telemetry_data, explainer_type='kernel')
    
    results = {
        'explainer': xai_explainer,
        'lime_explanations': [],
        'shap_values': None,
        'output_dir': xai_output_dir  # Return the actual directory used
    }
    
    # Generate explanations for anomalous instances
    if anomaly_instances is not None and len(anomaly_instances) > 0:
        print(f"\n5. Generating explanations for {len(anomaly_instances)} anomalous instances...")
        
        # LIME explanations (individual instances)
        print("\n   LIME Explanations:")
        for idx in range(min(3, len(anomaly_instances))):  # Explain first 3 anomalies
            print(f"     Explaining anomaly {idx+1}...")
            lime_exp = xai_explainer.explain_with_lime(anomaly_instances[idx], num_features=15)
            results['lime_explanations'].append(lime_exp)
            
            # Visualize and save to XAI folder
            save_path = os.path.join(xai_output_dir, f"lime_anomaly_{idx+1}.png")
            print(f"     Saving LIME plot to: {save_path}")
            xai_explainer.visualize_lime_explanation(lime_exp, save_path=save_path)
        
        # SHAP explanations (batch)
        print("\n   SHAP Explanations:")
        shap_values = xai_explainer.explain_with_shap(anomaly_instances, max_samples=10)
        results['shap_values'] = shap_values
        
        # Visualize SHAP summary
        save_path = os.path.join(xai_output_dir, "shap_summary.png")
        print(f"     Saving SHAP summary to: {save_path}")
        xai_explainer.visualize_shap_summary(
            shap_values,
            anomaly_instances[:10],
            save_path=save_path
        )
        
        # Visualize SHAP waterfall for first anomaly
        save_path = os.path.join(xai_output_dir, "shap_waterfall_anomaly_1.png")
        print(f"     Saving SHAP waterfall to: {save_path}")
        xai_explainer.visualize_shap_waterfall(
            shap_values,
            instance_idx=0,
            save_path=save_path
        )
    
    print("\n" + "="*70)
    print("XAI INTEGRATION COMPLETE")
    print("="*70)
    print(f"LIME explanations generated: {len(results['lime_explanations'])}")
    print(f"SHAP values computed: {'Yes' if results['shap_values'] is not None else 'No'}")
    print(f"Plots saved to: {xai_output_dir}")
    
    return results


# Example usage function for your system
def add_xai_to_basilisk_detection(ml_detector, real_telemetry, detected_faults, output_dir=None):
    """
    Add this function to your real_ml_fault_detection.py
    Call it after detecting faults to explain why they were detected.
    
    Args:
        ml_detector: Your RealMLFaultDetector instance
        real_telemetry: Dictionary of telemetry data
        detected_faults: List of RealDetectionResult objects
        output_dir: Directory to save plots (default: current/xai_plots)
    """
    print("\n" + "="*70)
    print("ADDING XAI EXPLANATIONS TO FAULT DETECTIONS")
    print("="*70)
    
    # Collect anomalous instances
    anomaly_instances = []
    
    for sc_name, sc_telemetry in real_telemetry.items():
        if sc_telemetry.get('fault_injected', False):
            # Extract the window where fault occurred
            fault_idx = int(len(sc_telemetry['rw_speeds']) * 15.0 / 30.0)
            fault_window = ml_detector.prepare_ml_input_from_real_data(
                ml_detector._get_windowed_data(sc_telemetry, fault_idx-10, fault_idx+10)
            )
            anomaly_instances.append(fault_window[0])  # Remove batch dimension
    
    if len(anomaly_instances) == 0:
        print("No anomalous instances found for XAI analysis")
        return None
    
    anomaly_instances = np.array(anomaly_instances)
    
    # Prepare background data (normal operation)
    background_instances = []
    for sc_name, sc_telemetry in real_telemetry.items():
        # Extract pre-fault window (normal operation)
        normal_idx = int(len(sc_telemetry['rw_speeds']) * 5.0 / 30.0)  # 5 minutes
        normal_window = ml_detector.prepare_ml_input_from_real_data(
            ml_detector._get_windowed_data(sc_telemetry, normal_idx-10, normal_idx+10)
        )
        background_instances.append(normal_window[0])
    
    background_instances = np.array(background_instances)
    
    # Run XAI analysis with output directory
    xai_results = integrate_xai_with_fault_detection(
        model_path=ml_detector.model_path,
        telemetry_data=background_instances,
        anomaly_instances=anomaly_instances,
        output_dir=output_dir
    )
    
    return xai_results


if __name__ == "__main__":
    print("Satellite Fault Detection XAI Module")
    print("="*70)
    print("\nThis module adds LIME and SHAP explainability to your fault detection system.")
    print("\nIntegration steps:")
    print("1. Train your autoencoder model (already done)")
    print("2. Detect faults using your existing system")
    print("3. Call add_xai_to_basilisk_detection() to explain detections")
    print("\nXAI will show:")
    print("  - LIME: Which features contributed most to each anomaly detection")
    print("  - SHAP: Overall feature importance across all detections")