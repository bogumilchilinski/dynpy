"""
This module provides basic tools for dynamic analysis of mechanical systems - it's change to generate synchronization conflict
"""

from .dynamics import (LagrangesDynamicSystem, LinearDynamicSystem,
                       HarmonicOscillator, WeakNonlinearOscillator, DampedHarmonicOscillator)

from .continuous import ContinuousSystem