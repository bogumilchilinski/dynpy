import base64
import random
from numbers import Number

import IPython as IP
import matplotlib.pyplot as plt
import numpy as np
from pylatex import TikZ, TikZNode
from sympy import (
    Derivative,
    Eq,
    Expr,
    Function,
    Heaviside,
    Integral,
    Matrix,
    S,
    Symbol,
    atan,
    cos,
    diag,
    diff,
    latex,
    pi,
    sin,
    sqrt,
    symbols,
)
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint

from ...dynamics import (
    GeometryScene,
    HarmonicOscillator,
    LagrangesDynamicSystem,
    base_frame,
    base_origin,
)
from ..elements import (
    # DerivativeElement,
    Force,
    # IntegralElement,
    # ProportionalElement,
)
from ...utilities.templates import tikz

# from dgeometry import GeometryScene,np
from ..mechanics.trolley import (
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)

#1194
class InertialThermalSystem(ComposedSystem):
    pass
