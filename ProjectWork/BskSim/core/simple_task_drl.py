#!/usr/bin/env python
"""
simple_task_drl.py

Lightweight DRL model for cluster-based task reassignment.
Trained on expert rules to provide neural network-based decision making.
"""

import numpy as np
import os

try:
    import tensorflow as tf
    TF_AVAILABLE = True
except ImportError:
    TF_AVAILABLE = False
    print("TensorFlow not available - DRL will not function")


class ClusterTaskReassignmentDRL:
    """
    Neural network-based DRL agent for task reassignment decisions.
    Maps constellation state (10D) to reassignment strategy (3 actions).
    """
    
    def __init__(self, state_dim=10, action_dim=3, model_path=None):
        """
        Args:
            state_dim: Dimension of state vector (default 10)
            action_dim: Number of actions (default 3: even, capability, load-balanced)
            model_path: Optional path to load pre-trained model
        """
        self.state_dim = state_dim
        self.action_dim = action_dim
        self.model = None
        self.is_loaded = False
        
        if not TF_AVAILABLE:
            print("Cannot initialize DRL - TensorFlow not installed")
            return
        
        # Try to load existing model
        if model_path and os.path.exists(model_path):
            try:
                self.model = tf.keras.models.load_model(model_path)
                self.is_loaded = True
                print(f"DRL model loaded from {model_path}")
                return
            except Exception as e:
                print(f"Could not load model from {model_path}: {e}")
        
        # Create and train new model
        self._build_network()
        self._train_from_expert_rules()
        self.is_loaded = True
    
    def _build_network(self):
        """Build the neural network architecture"""
        self.model = tf.keras.Sequential([
            tf.keras.layers.Dense(64, activation='relu', input_shape=(self.state_dim,)),
            tf.keras.layers.Dropout(0.2),
            tf.keras.layers.Dense(64, activation='relu'),
            tf.keras.layers.Dropout(0.2),
            tf.keras.layers.Dense(32, activation='relu'),
            tf.keras.layers.Dense(self.action_dim, activation='softmax')
        ])
        
        self.model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
            loss='categorical_crossentropy',
            metrics=['accuracy']
        )
        
        print("DRL network architecture created")
    
    def _train_from_expert_rules(self):
        """
        Train the network using expert rule-based logic.
        This creates a DRL agent that learns optimal decisions from domain expertise.
        """
        print("Training DRL model from expert rules...")
        
        X_train = []
        y_train = []
        
        # Generate diverse training scenarios (500 samples)
        np.random.seed(42)
        
        for _ in range(500):
            # Generate random state vector
            state = np.random.rand(self.state_dim).astype(np.float32)
            
            # Apply expert rules to determine optimal action
            healthy_ratio = state[0]
            faulty_ratio = state[1]
            avg_confidence = state[2]
            max_confidence = state[3]
            avg_anomaly = state[4]
            max_anomaly = state[5]
            
            # Expert Rule 1: Low healthy ratio -> capability-based
            # Focus on most capable remaining satellites
            if healthy_ratio < 0.4:
                label = [0.05, 0.90, 0.05]  # Strong preference for capability-based
            
            # Expert Rule 2: Moderate healthy ratio + high severity -> load-balanced
            # Distribute carefully to avoid overloading
            elif 0.4 <= healthy_ratio < 0.6 and (avg_anomaly > 0.6 or max_anomaly > 0.7):
                label = [0.05, 0.05, 0.90]  # Strong preference for load-balanced
            
            # Expert Rule 3: Very high severity regardless of ratio -> load-balanced
            # Critical situation requires careful distribution
            elif max_confidence > 0.85 or max_anomaly > 0.8:
                label = [0.1, 0.1, 0.8]  # Prefer load-balanced
            
            # Expert Rule 4: Low to moderate severity + good health -> even distribution
            # Normal operation, distribute evenly
            elif healthy_ratio >= 0.6 and avg_anomaly < 0.5:
                label = [0.85, 0.1, 0.05]  # Strong preference for even distribution
            
            # Expert Rule 5: Moderate healthy ratio + moderate severity -> mix
            # Balanced approach between capability and load
            elif 0.5 <= healthy_ratio < 0.7 and 0.4 <= avg_anomaly <= 0.6:
                label = [0.2, 0.5, 0.3]  # Prefer capability-based with load consideration
            
            # Default: Even distribution
            else:
                label = [0.7, 0.2, 0.1]
            
            X_train.append(state)
            y_train.append(label)
        
        X_train = np.array(X_train, dtype=np.float32)
        y_train = np.array(y_train, dtype=np.float32)
        
        # Train the model
        history = self.model.fit(
            X_train, y_train,
            epochs=100,
            batch_size=32,
            validation_split=0.2,
            verbose=0
        )
        
        final_accuracy = history.history['accuracy'][-1]
        final_val_accuracy = history.history['val_accuracy'][-1]
        
        print(f"DRL training complete:")
        print(f"  Training accuracy: {final_accuracy:.3f}")
        print(f"  Validation accuracy: {final_val_accuracy:.3f}")
        print(f"  Model learned from {len(X_train)} expert scenarios")
    
    def select_action(self, state, training=False):
        """
        Select action (task reassignment strategy) based on current state.
        
        Args:
            state: numpy array of shape (10,) representing constellation state
            training: whether in training mode (not used for inference)
        
        Returns:
            int: Action index (0=even_distribution, 1=capability_based, 2=load_balanced)
        """
        if not self.is_loaded or self.model is None:
            # Fallback to simple rule if model not available
            return self._fallback_rule_based(state)
        
        # Ensure state is correct shape
        state = np.array(state, dtype=np.float32)
        if len(state) != self.state_dim:
            if len(state) < self.state_dim:
                state = np.pad(state, (0, self.state_dim - len(state)))
            else:
                state = state[:self.state_dim]
        
        state = state.reshape(1, -1)
        
        # Get action probabilities from network
        action_probs = self.model.predict(state, verbose=0)[0]
        
        # Select action with highest probability
        action = int(np.argmax(action_probs))
        
        return action
    
    def get_action_probabilities(self, state):
        """
        Get probability distribution over actions.
        
        Args:
            state: numpy array of shape (10,)
        
        Returns:
            numpy array of shape (3,) with action probabilities
        """
        if not self.is_loaded or self.model is None:
            # Uniform distribution if model not available
            return np.array([0.33, 0.33, 0.34], dtype=np.float32)
        
        # Ensure state is correct shape
        state = np.array(state, dtype=np.float32)
        if len(state) != self.state_dim:
            if len(state) < self.state_dim:
                state = np.pad(state, (0, self.state_dim - len(state)))
            else:
                state = state[:self.state_dim]
        
        state = state.reshape(1, -1)
        
        # Get probabilities
        action_probs = self.model.predict(state, verbose=0)[0]
        
        return action_probs
    
    def _fallback_rule_based(self, state):
        """Fallback rule-based decision if model unavailable"""
        healthy_ratio = state[0] if len(state) > 0 else 0.5
        fault_severity = state[4] if len(state) > 4 else 0.5
        
        if healthy_ratio < 0.5:
            return 1  # capability_based
        elif fault_severity > 0.7:
            return 2  # load_balanced
        else:
            return 0  # even_distribution
    
    def save_model(self, path):
        """Save the trained model to disk"""
        if self.model is not None:
            self.model.save(path)
            print(f"DRL model saved to {path}")
    
    def evaluate(self, test_states):
        """
        Evaluate model on test states.
        
        Args:
            test_states: list of state vectors
        
        Returns:
            dict with evaluation metrics
        """
        if not self.is_loaded or self.model is None:
            return {"error": "Model not loaded"}
        
        actions = []
        confidences = []
        
        for state in test_states:
            action = self.select_action(state)
            probs = self.get_action_probabilities(state)
            confidence = float(np.max(probs))
            
            actions.append(action)
            confidences.append(confidence)
        
        return {
            "num_states": len(test_states),
            "actions": actions,
            "avg_confidence": np.mean(confidences),
            "action_distribution": {
                "even_distribution": actions.count(0),
                "capability_based": actions.count(1),
                "load_balanced": actions.count(2)
            }
        }


def test_drl_agent():
    """Test the DRL agent with sample scenarios"""
    print("\nTesting DRL Agent")
    print("=" * 60)
    
    agent = ClusterTaskReassignmentDRL()
    
    # Test scenario 1: Healthy constellation
    state1 = np.array([0.8, 0.2, 0.3, 0.4, 0.2, 0.3, 0.1, 0.1, 0.2, 0.0])
    action1 = agent.select_action(state1)
    probs1 = agent.get_action_probabilities(state1)
    print(f"\nScenario 1 (Healthy constellation):")
    print(f"  State: healthy_ratio=0.8, avg_anomaly=0.2")
    print(f"  Action: {action1} (0=even, 1=capability, 2=load)")
    print(f"  Probabilities: {probs1}")
    
    # Test scenario 2: Critical failures
    state2 = np.array([0.3, 0.7, 0.9, 0.95, 0.85, 0.9, 0.5, 0.3, 0.6, 0.0])
    action2 = agent.select_action(state2)
    probs2 = agent.get_action_probabilities(state2)
    print(f"\nScenario 2 (Critical failures):")
    print(f"  State: healthy_ratio=0.3, avg_anomaly=0.85")
    print(f"  Action: {action2}")
    print(f"  Probabilities: {probs2}")
    
    # Test scenario 3: Moderate situation
    state3 = np.array([0.5, 0.5, 0.6, 0.7, 0.5, 0.6, 0.3, 0.2, 0.4, 0.0])
    action3 = agent.select_action(state3)
    probs3 = agent.get_action_probabilities(state3)
    print(f"\nScenario 3 (Moderate situation):")
    print(f"  State: healthy_ratio=0.5, avg_anomaly=0.5")
    print(f"  Action: {action3}")
    print(f"  Probabilities: {probs3}")
    
    print("\nDRL Agent Test Complete")


if __name__ == "__main__":
    if TF_AVAILABLE:
        test_drl_agent()
    else:
        print("TensorFlow not available - cannot test DRL agent")
        print("Install TensorFlow: pip install tensorflow")