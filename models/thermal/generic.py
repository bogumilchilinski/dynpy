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
    Damper,
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

    m = Symbol('m', positive=True)
    cp = Symbol('cp', positive=True)

    UA = Symbol('UA', positive=True)

    T = dynamicsymbols('T')
    Qin = dynamicsymbols('Qin')

    ivar = Symbol('t')

    def __init__(self,
                 m=None,
                 cp=None,
                 UA=None,
                 T=None,
                 Qin=None,
                 ivar=None,
                 **kwargs):

        if m is not None:
            self.m = m

        if cp is not None:
            self.cp = cp

        if UA is not None:
            self.UA = UA

        if T is not None:
            self.T = T

        if Qin is not None:
            self.Qin = Qin

        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.T]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        Cth = self.m * self.cp

        self._thermal_capacity = MaterialPoint(
            Cth,
            pos1=self.T,
            qs=self.qs
        )

        self._heat_losses = Damper(
            self.UA,
            pos1=self.T,
            qs=self.qs
        )(label='heat losses')

        self._heat_source = Force(
            self.Qin,
            self.T,
            qs=self.qs
        )(label='heat input')

        components['_thermal_capacity'] = self._thermal_capacity
        components['_heat_losses'] = self._heat_losses
        components['_heat_source'] = self._heat_source

        return components
