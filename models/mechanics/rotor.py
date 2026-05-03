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

ureg = units

from .principles import (
    REPORT_COMPONENTS_LIST,
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)


class MasslessElasticShaft(ComposedSystem):
    """
    A simple elastic shaft model, consisting of a rigid body (the shaft) connected to a spring and damper.
    """
    

    def _components(self):
        return [

            ("spring", Spring("spring", stiffness=self.young_modulus)),
            ("damper", Damper("damper", damping_coefficient=self.density)),
        ]
        
class VibratingRotor(ComposedSystem):
    """
    A simple vibrating rotor model, consisting of a rigid body (the rotor) connected to a spring and damper.
    """
    
    def _components(self):
        return [
            ("rotor", RigidBody2D("rotor", base_frame, base_origin)),
            ("spring", MasslessElasticShaft("spring", stiffness=self.young_modulus)),
            ("damper", Damper("damper", damping_coefficient=self.density)),
        ]