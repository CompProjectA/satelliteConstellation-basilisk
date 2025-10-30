#!/usr/bin/env python
"""
drl_model_manager.py

Centralized manager for all DRL models and checkpoints.
Place this file in: core/drl_model_manager.py
"""

import os
import re
import glob
import json
from typing import Dict, Optional, List, Any
from datetime import datetime

# Canonical project paths
CORE_DIR = os.path.dirname(os.path.abspath(__file__))           # .../BskSim/core
ROOT_DIR = os.path.dirname(CORE_DIR)                            # .../BskSim
DRL_DIR = os.path.join(ROOT_DIR, "DRL")
RESULT_DIR = os.path.join(DRL_DIR, "result")
MODELS_DIR = os.path.join(ROOT_DIR, "models")
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(MODELS_DIR, exist_ok=True)


class DRLModelManager:
    """Manages DRL model loading, saving, discovery and inventory."""

    SUPPORTED_EXTS = (".keras", ".h5", ".pth", ".pt", ".zip", ".pkl")

    def __init__(self, models_dir: Optional[str] = None):
        self.core_dir = CORE_DIR
        self.root_dir = ROOT_DIR
        self.drl_dir = DRL_DIR
        self.result_dir = RESULT_DIR
        self.models_path = os.path.abspath(models_dir or MODELS_DIR)
        os.makedirs(self.models_path, exist_ok=True)

        # Initial registry
        self.available_models: Dict[str, Dict[str, Any]] = {
            "anomaly_detection": {
                "path": os.path.join(self.core_dir, "anomaly_detection_model.keras"),
                "type": "keras",
                "purpose": "Fault detection (Keras)",
                "loaded": False,
                "model": None
            },
            "task_reassignment": {
                "path": os.path.join(self.core_dir, "drl_task_reassignment.h5"),
                "type": "h5",
                "purpose": "Baseline DRL (H5/Keras)",
                "loaded": False,
                "model": None
            }
        }

        self._scan_for_models()

    # Discovery
    def _scan_for_models(self) -> None:
        """Scan DRL/result and models/ for trained artifacts and add them to registry."""
        for directory in (self.result_dir, self.models_path):
            if not os.path.isdir(directory):
                continue
            for name in os.listdir(directory):
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    ext = os.path.splitext(name)[1].lower()
                    if ext in self.SUPPORTED_EXTS:
                        key = os.path.splitext(name)[0]
                        self.available_models.setdefault(key, {
                            "path": path,
                            "type": ext.lstrip("."),
                            "purpose": "Trained DRL policy (auto-discovered)",
                            "loaded": False,
                            "model": None
                        })
                else:
                    # RLlib checkpoints often appear as directories
                    if os.path.isdir(path) and any(
                        os.path.exists(os.path.join(path, f))
                        for f in ("params.json", "algorithm_state.pkl", ".is_checkpoint")
                    ):
                        key = os.path.basename(path)
                        self.available_models.setdefault(key, {
                            "path": path,
                            "type": "rllib_dir",
                            "purpose": "RLlib checkpoint (dir)",
                            "loaded": False,
                            "model": None
                        })

    def refresh(self) -> None:
        """Refresh the registry (non-destructive)."""
        self._scan_for_models()

    # Queries
    def list_models(self) -> List[Dict]:
        """List all available models with metadata."""
        return [{"name": name, **info} for name, info in self.available_models.items()]

    def get_model_info(self, model_name: str) -> Optional[Dict]:
        return self.available_models.get(model_name, None)

    def verify_all_models(self) -> Dict[str, bool]:
        """Check if files/dirs referenced by registry exist."""
        results = {}
        for name, info in self.available_models.items():
            exists = os.path.exists(info["path"])
            results[name] = exists
            if exists:
                print(f"{name}: Found at {info['path']}")
            else:
                print(f"{name}: NOT FOUND at {info['path']}")
        return results

    # Loading helpers
    def _load_keras(self, path: str):
        import tensorflow as tf
        return tf.keras.models.load_model(path, compile=False)

    def _load_pytorch(self, path: str):
        import torch
        map_location = "cpu" if not torch.cuda.is_available() else None
        return torch.load(path, map_location=map_location)

    def _load_stable_baselines3_zip(self, path: str):
        """Load Stable-Baselines3 .zip if SB3 is installed."""
        try:
            from stable_baselines3.common.save_util import load_from_zip_file
            data, params, pytorch_variables = load_from_zip_file(path)
            return {"sb3_data": bool(data), "sb3_params": params is not None, "raw": (data, params, pytorch_variables)}
        except Exception as e:
            raise RuntimeError(f"Stable-Baselines3 not available to load {os.path.basename(path)}") from e

    def _load_rllib_checkpoint(self, path: str):
        """Load an RLlib Algorithm from checkpoint directory."""
        try:
            from ray.rllib.algorithms.algorithm import Algorithm
            return Algorithm.from_checkpoint(path)
        except Exception as e:
            raise RuntimeError(f"Failed to load RLlib checkpoint at {path}") from e

    # Public loading API
    def load_model(self, model_name: str, force_reload: bool = False):
        """Load a model by registry name."""
        if model_name not in self.available_models:
            raise ValueError(f"Model '{model_name}' not found in registry")

        info = self.available_models[model_name]

        if info.get("loaded") and not force_reload:
            return info["model"]

        path = info["path"]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Model file not found: {path}")

        ext = info["type"].lower()
        obj = None

        try:
            if ext in ("keras", "h5"):
                obj = self._load_keras(path)
            elif ext in ("pt", "pth"):
                obj = self._load_pytorch(path)
            elif ext == "zip":
                obj = self._load_stable_baselines3_zip(path)
            elif ext == "pkl":
                import pickle
                with open(path, "rb") as f:
                    obj = pickle.load(f)
            elif ext == "rllib_dir":
                obj = self._load_rllib_checkpoint(path)
            else:
                raise ValueError(f"Unsupported model type: {ext}")

            info["model"] = obj
            info["loaded"] = True
            print(f"Loaded model: {model_name} ({ext})")
            return obj

        except Exception as e:
            print(f"Failed to load model '{model_name}': {e}")
            raise

    # Latest by prefix
    def load_latest_by_prefix(self, prefix: str):
        """
        Find and load the newest artifact in DRL/result or models/ whose base name
        starts with `prefix` (case-insensitive). Returns the loaded object.
        """
        candidates = []

        def add_files_from(d: str):
            if not os.path.isdir(d):
                return
            for ext in self.SUPPORTED_EXTS:
                for p in glob.glob(os.path.join(d, f"*{ext}")):
                    base = os.path.basename(p)
                    name_no_ext = os.path.splitext(base)[0]
                    if name_no_ext.lower().startswith(prefix.lower()):
                        try:
                            mtime = os.path.getmtime(p)
                        except Exception:
                            mtime = 0
                        candidates.append((mtime, name_no_ext, p, ext.lstrip(".")))

        add_files_from(self.result_dir)
        add_files_from(self.models_path)

        # RLlib checkpoint directories
        for d in (self.result_dir, self.models_path):
            if not os.path.isdir(d):
                continue
            for entry in os.listdir(d):
                full = os.path.join(d, entry)
                if os.path.isdir(full) and entry.lower().startswith(prefix.lower()):
                    if any(os.path.exists(os.path.join(full, f))
                           for f in ("params.json", "algorithm_state.pkl", ".is_checkpoint")):
                        try:
                            mtime = os.path.getmtime(full)
                        except Exception:
                            mtime = 0
                        candidates.append((mtime, entry, full, "rllib_dir"))

        if not candidates:
            raise FileNotFoundError(
                f"No artifact starting with '{prefix}' found in: {self.result_dir}, {self.models_path}"
            )

        candidates.sort(key=lambda t: t[0], reverse=True)
        _, base_name, path, type_str = candidates[0]

        if base_name not in self.available_models:
            self.available_models[base_name] = {
                "path": path,
                "type": type_str,
                "purpose": "Trained DRL policy (auto-discovered)",
                "loaded": False,
                "model": None
            }
        return self.load_model(base_name)

    def load_latest_family(self, family: str):
        """
        Convenience: load latest PPO/TDHD/TD3HD/DQN by common prefixes.
        family ∈ {'PPO','TDHD','TD3HD','DQN'}
        """
        prefix_map = {
            "PPO": ("PPO",),
            "TD3HD": ("TD3HD", "TDHD", "TD3_Hindsight", "TDHDYear2"),
            "TDHD": ("TDHD", "TD3HD"),
            "DQN": ("DQN", "DQNYear2")
        }
        candidates = prefix_map.get(family.upper(), (family,))
        last_err = None
        for pref in candidates:
            try:
                return self.load_latest_by_prefix(pref)
            except Exception as e:
                last_err = e
        raise last_err or RuntimeError(f"No models matched for family '{family}'")

    # Saving
    def save_model(self, model, model_name: str, model_type: str, purpose: str = "Custom"):
        """
        Save a model to the models directory; updates registry.
        model_type: 'keras' | 'h5' | 'pytorch' | 'rllib' | 'pkl'
        """
        ext_map = {
            "keras": ".keras",
            "h5": ".h5",
            "pytorch": ".pth",
            "rllib": "",        # RLlib saves to a directory
            "pkl": ".pkl"
        }

        if model_type == "rllib":
            save_dir = os.path.join(self.models_path, f"{model_name}_rllib")
            os.makedirs(save_dir, exist_ok=True)
            try:
                model.save(save_dir)  # RLlib Algorithm.save
            except Exception as e:
                raise RuntimeError(f"Failed to save RLlib model to {save_dir}") from e
            path = save_dir
            ftype = "rllib_dir"
        else:
            ext = ext_map.get(model_type, ".pkl")
            path = os.path.join(self.models_path, f"{model_name}{ext}")
            try:
                if model_type in ("keras", "h5"):
                    model.save(path)
                elif model_type == "pytorch":
                    import torch
                    torch.save(model, path)
                elif model_type == "pkl":
                    import pickle
                    with open(path, "wb") as f:
                        pickle.dump(model, f)
                else:
                    import pickle
                    with open(path, "wb") as f:
                        pickle.dump(model, f)
                ftype = model_type
            except Exception as e:
                raise RuntimeError(f"Failed to save model '{model_name}' to {path}") from e

        self.available_models[model_name] = {
            "path": path,
            "type": ftype,
            "purpose": purpose,
            "loaded": True,
            "model": model,
            "saved_at": datetime.now().isoformat()
        }

        print(f"Saved model: {model_name} -> {path}")
        return path

    # Export
    def export_model_registry(self, output_filename: str = "model_registry.json") -> str:
        """Export model registry to JSON (under models dir)."""
        out_path = os.path.join(self.models_path, output_filename)
        export = {}
        for name, info in self.available_models.items():
            export[name] = {
                "path": info["path"],
                "type": info["type"],
                "purpose": info.get("purpose", ""),
                "exists": os.path.exists(info["path"]),
                "loaded": info.get("loaded", False)
            }
        with open(out_path, "w") as f:
            json.dump(export, f, indent=2)
        print(f"Model registry exported to: {out_path}")
        return out_path


# Singleton helpers
def get_model_manager() -> DRLModelManager:
    if not hasattr(get_model_manager, "_instance"):
        get_model_manager._instance = DRLModelManager()
    return get_model_manager._instance


def verify_all_drl_models():
    mgr = get_model_manager()
    return mgr.verify_all_models()


if __name__ == "__main__":
    print("DRL Model Manager")
    print("=" * 60)
    manager = get_model_manager()

    print("\nAvailable Models:")
    print("-" * 60)
    for model in manager.list_models():
        status = "Loaded" if model.get("loaded") else "Not loaded"
        exists = "Yes" if os.path.exists(model["path"]) else "No"
        print(f"{model['name']}: {model.get('purpose','')}")
        print(f"   Path: {model['path']}")
        print(f"   Type: {model['type']}, Exists: {exists}, Status: {status}")
        print()

    print("Verifying All Models:")
    print("-" * 60)
    results = manager.verify_all_models()
    missing = [name for name, exists in results.items() if not exists]
    if missing:
        print(f"Missing models: {', '.join(missing)}")
    else:
        print("All models found.")

    manager.export_model_registry()
