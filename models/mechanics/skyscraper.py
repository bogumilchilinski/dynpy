import base64
import inspect
import random

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
    atan,
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
from sympy.physics import units
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint

from ...continuous import ContinuousSystem, PlaneStressProblem
from ...dynamics import HarmonicOscillator, LagrangesDynamicSystem, mech_comp
from ...utilities.components.mech import en as mech_comp
from ..elements import (
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
from .pendulum import Pendulum, PendulumKinematicExct
from .principles import (
    REPORT_COMPONENTS_LIST,
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)
from .tmd import TunedMassDamperRelativeMotion

ureg = units

from functools import cached_property, lru_cache


class skyscraper(ComposedSystem):
    """Ready to use sample Single Degree of Freedom System with mass on spring
    Arguments:
    =========
        m = Mass
            -Mass of system on spring

        k = Spring coefficient
            -Spring carrying the system

        ivar = symbol object
            -Independant time variable

        qs = dynamicsymbol object
            -Generalized coordinates

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

    >>> t = symbols('t')
    >>> m, k = symbols('m, k')
    >>> qs = dynamicsymbols('z') # Generalized Coordinates
    >>> mass = SDoFHarmonicOscillator(m,k, qs=[z],) # Initialization of LagrangesDynamicSystem instance

    -We define the symbols and dynamicsymbols
    -Kinetic energy T and potential energy v are evaluated to calculate the lagrangian L
    -Reference frame was created with point P defining the position and the velocity determined on the z axis
    -external forces assigned
    -Next we determine the instance of the system using class LagrangeDynamicSystem
    -We call out the instance of the class
    -If necessary assign values for the default arguments

    """

    scheme_name = "engine.png"
    real_name = "engine_real.PNG"

    m = Symbol("m", positive=True)
    k = Symbol("k_b", positive=True)
    F = Symbol("F", positive=True)
    c = Symbol("c", positive=True)
    ivar = Symbol("t")
    g = Symbol("g", positive=True)
    h = Symbol("h", positive=True)

    x = dynamicsymbols("x")

    def __init__(
        self,
        g=None,
        h=None,
        m=None,
        k=None,
        F=None,
        x=None,
        c=None,
        ivar=None,
        **kwargs
    ):

        if g is not None:
            self.g = g
        if h is not None:
            self.h = h
        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if F is not None:
            self.F = F
        if x is not None:
            self.x = x
        if c is not None:
            self.c = c
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.x]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.x, qs=self.qs)
        self.stiffness = Spring(self.k, self.x, qs=self.qs)
        self.force = Force(self.F, self.x, qs=self.qs)
        self.damper = Damper(self.c, self.x, qs=self.qs)
        self._gravity = GravitationalForce(self.m, self.g, pos1=self.h, qs=self.qs)

        components["material_point"] = self.material_point
        components["stiffness"] = self.stiffness
        components["force"] = self.force
        components["damper"] = self.damper
        components["gravity"] = self._gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict
