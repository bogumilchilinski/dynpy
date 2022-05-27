"""
This module provides basic tools for dynamic analysis of mechanical systems
"""

from .dynamics import (LagrangesDynamicSystem, LinearDynamicSystem,
                       HarmonicOscillator, WeakNonlinearOscillator, DampedHarmonicOscillator)

from .continuous import ContinuousSystem