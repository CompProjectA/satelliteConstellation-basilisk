#!/usr/bin/env python
"""
__init__.py

Initializes the gui_tab package and exports its components.
"""

# Import base class first
try:
    from .base_tab import BaseTab
except ImportError as e:
    print(f"Warning: Could not import BaseTab: {e}")
    
    # Create a minimal BaseTab if the import fails
    class BaseTab:
        def __init__(self, parent_app, parent_frame):
            self.parent_app = parent_app
            self.parent_frame = parent_frame
        
        def add_help_button(self, parent, title, topic=None, command=None):
            pass
        
        def show_help_content(self, title, topic=None, content=None):
            pass
        
        def show_message(self, title, message, message_type="info"):
            pass
        
        def add_log(self, message):
            if hasattr(self.parent_app, 'add_log'):
                self.parent_app.add_log(message)
        

# Import tab classes with error handling
try:
    from .constellation_tab import ConstellationTab
except ImportError as e:
    print(f"Warning: Could not import ConstellationTab: {e}")
    ConstellationTab = None

try:
    from .fault_tab import FaultTab
except ImportError as e:
    print(f"Warning: Could not import FaultTab: {e}")
    FaultTab = None

try:
    from .fault_detection_tab import FaultDetectionTab
except ImportError as e:
    print(f"Warning: Could not import FaultDetectionTab: {e}")
    FaultDetectionTab = None

try:
    from .target_tab import TargetTab
except ImportError as e:
    print(f"Warning: Could not import TargetTab: {e}")
    TargetTab = None

try:
    from .visualization_tab import VisualizationTab
except ImportError as e:
    print(f"Warning: Could not import VisualizationTab: {e}")
    VisualizationTab = None

try:
    from .output_tab import OutputTab
except ImportError as e:
    print(f"Warning: Could not import OutputTab: {e}")
    OutputTab = None

try:
    from .results_tab import ResultsTab
except ImportError as e:
    print(f"Warning: Could not import ResultsTab: {e}")
    ResultsTab = None

# NEW: Import CommunicationTab
try:
    from .communication_tab import CommunicationTab
except ImportError as e:
    print(f"Warning: Could not import CommunicationTab: {e}")
    CommunicationTab = None

# Export all available classes
__all__ = []

if BaseTab:
    __all__.append('BaseTab')
if ConstellationTab:
    __all__.append('ConstellationTab')
if FaultTab:
    __all__.append('FaultTab')
if FaultDetectionTab:
    __all__.append('FaultDetectionTab')
if TargetTab:
    __all__.append('TargetTab')
if VisualizationTab:
    __all__.append('VisualizationTab')
if OutputTab:
    __all__.append('OutputTab')
if ResultsTab:
    __all__.append('ResultsTab')
if CommunicationTab:
    __all__.append('CommunicationTab')