# faults/__init__.py
# Marks this directory as a package and provides friendly exports.

from .friction_fault import FrictionFaultScenario as _FFS  # noqa: F401
from .powerlimit_fault import PowerLimitFaultScenario as _PLFS  # noqa: F401
from .encoder_fault import EncoderFaultScenario as _EFS  # noqa: F401
# battery exposes only run() right now; scenarios are internal to that module.
