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
from ...models.elements import (
    DerivativeElement,
    Force,
    IntegralElement,
    ProportionalElement,
)
from ...utilities.templates import tikz

# from dgeometry import GeometryScene,np


class PIController(Force):
    """
    A proportional–integral controller (PI controller) is a control loop mechanism
    employing feedback used in industrial control systems requiring continuously
    modulated control. A PI controller calculates an error value (e(t)) as the
    difference between a desired setpoint (SP) and a measured process variable (PV)
    and applies a correction based on proportional and integral terms. It consists of
    two independent elements - P and I that being combined return two-term controller.
    """

    k_P = Symbol("k_P", positive=True)
    k_I = Symbol("k_I", positive=True)
    error = Symbol("e", positive=True)
    target = Symbol("target", positive=True)
    reference = S.Zero

    def __init__(
        self,
        k_P,
        k_I,
        error=None,
        target=None,
        reference=None,
        qs=None,
        ivar=Symbol("t"),
        frame=base_frame,
        **kwargs
    ):

        if k_P is not None:
            self.k_P = k_P
        if k_I is not None:
            self.k_I = k_I

        if reference is not None:
            self.reference = reference

        if error is not None:
            self.error = error

        if target is None and isinstance(self.error, Function):
            target = self.error

        self.ivar = ivar
        self.qs = qs

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        k_P, k_I = self.k_P, self.k_I

        self._proporional_elem = ProportionalElement(
            k_P, self.error, self.reference, qs=self.qs, ivar=self.ivar
        )
        self._integral_elem = IntegralElement(
            k_I, self.error, self.reference, qs=self.qs, ivar=self.ivar
        )

        components["proportional_elem"] = self._proporional_elem
        components["integral_elem"] = self._integral_elem

        return components
