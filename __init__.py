"""
This module provides basic tools for analysis of dynamic system (described by time-depended functions)
"""

from .dynamics import (
    DampedHarmonicOscillator,
    HarmonicOscillator,
    LagrangesDynamicSystem,
    LinearDynamicSystem,
    WeakNonlinearOscillator,
)

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("Dyn4py")
except PackageNotFoundError:
    __version__ = "0.0.0"

# from .continuous import ContinuousSystem
