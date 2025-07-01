import base64
import inspect
import random
from functools import cached_property, lru_cache

import IPython as IP
import numpy as np
from sympy import (
    Abs,
    Eq,
    Function,
    Matrix,
    N,
    Number,
    S,
    Subs,
    Symbol,
    Tuple,
    asin,
    cos,
    diag,
    diff,
    dsolve,
    factorial,
    flatten,
    fraction,
    hessian,
    im,
    latex,
    oo,
    pi,
    sin,
    solve,
    solveset,
    sqrt,
    symbols,
)
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint

from ..continuous import ContinuousSystem, PlaneStressProblem
from ..dynamics import HarmonicOscillator, LagrangesDynamicSystem, mech_comp
from .elements import (
    PID,
    Damper,
    Disk,
    Excitation,
    Force,
    GravitationalForce,
    MaterialPoint,
    RigidBody2D,
    Spring,
    base_frame,
    base_origin,
)
from .mechanics.bridge import BeamBridge
from .mechanics.engine import (
    BoxerEnginePerpendicularSprings,
    DampedEngine,
    DampedEngineVerticalSpringGravity,
    Engine,
    EngineVerticalSpringGravity,
    InlineEnginePerpendicularSprings,
    NonLinearBoxerEnginePerpendicularSprings,
    NonLinearInlineEnginePerpendicularSpringsGravity,
)
from .mechanics.pendulum import (
    DampedPendulum,
    ExcitedDampedPendulum,
    ExcitedPendulum,
    FreePendulum,
    Pendulum,
    PendulumKinematicExct,
    PulledPendulum,
)
from .mechanics.tensioner import BlowerToothedBelt, DampedBlowerToothedBelt
from .mechanics.trolley import (
    ForcedDampedTrolleyWithSpring,
    ForcedTrolleysWithSprings,
    ForcedTrolleyWithSpring,
    SpringDamperMassSystem,
    SpringMassSystem,
    TrolleyWithPendulum,
)
