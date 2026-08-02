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
    MaterialPoint,
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
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    c = Symbol("c", positive=True)
    x = dynamicsymbols("x", positive=True)
    y = dynamicsymbols("y", positive=True)
    c0 = Symbol("c0", positive=True)
    ivar = Symbol("t")

    def __init__(self, m=None, g=None, c=None, x=None, y=None, ivar=None, **kwargs):

        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if c is not None:
            self.c = c
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.x, self.y]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        v_x = (self.x).diff(self.ivar)
        v_y = (self.y).diff(self.ivar)
        components = {}

        self._missile_airless = MaterialPoint(
            self.m, x=self.x, y=self.y, qs=self.qs
        )

        self._drag_x = Damper(self.c, self.x, qs=self.qs)(label="horizontal drag")

        self._drag_y = Damper(self.c, self.y, qs=self.qs)(label="drag")

        components["_missile_airless"] = self._missile_airless
        components["_drag_x"] = self._drag_x
        components["_drag_y"] = self._drag_y

        return components
