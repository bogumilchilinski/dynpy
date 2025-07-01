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
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint

from ...dynamics import HarmonicOscillator, LagrangesDynamicSystem, mech_comp
from ..continuous import ContinuousSystem, PlaneStressProblem
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
from .principles import (
    REPORT_COMPONENTS_LIST,
    ComposedSystem,
    NonlinearComposedSystem,
    SpringMassSystem,
    base_frame,
    base_origin,
)


class TunedMassDamper(SpringMassSystem):
    """
    Ready to use sample Single Degree of Freedom System of Tuned Mass Damper


    Arguments:
    =========
        m = Mass
            -Mass of TMD

        k = Spring coefficient
            -Spring carrying the TMD

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

    scheme_name = "tmd_engine_vertical_spring_nogravity.png"
    real_name = "tmd_engine_real.jpg"

    m_E = Symbol("m_E", positive=True)
    k_E = Symbol("k_E", positive=True)
    z = dynamicsymbols("z")
    z_E = dynamicsymbols("z_E")

    def __init__(
        self, m_E=None, k_E=None, z_E=None, z=None, ivar=Symbol("t"), **kwargs
    ):

        if m_E is not None:
            self.m_E = m_E
        if k_E is not None:
            self.k_E = k_E
        if z is not None:
            self.z = z
        if z_E is not None:
            self.z_E = z_E
        self.ivar = ivar

        self.qs = [self.z_E]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self.m_E, self.z_E, qs=self.qs)(
            label="Material point - mass of the TMD"
        )
        self._spring = Spring(self.k_E, pos1=self.z_E, pos2=self.z, qs=self.qs)(
            label="Spring - stiffness of the spring"
        )

        components["_material_point"] = self._material_point
        components["_spring"] = self._spring

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict


class TunedMassDamperRelativeMotion(SpringMassSystem):
    """
    Ready to use sample Single Degree of Freedom System of Tuned Mass Damper


    Arguments:
    =========
        m = Mass
            -Mass of TMD

        k = Spring coefficient
            -Spring carrying the TMD

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

    scheme_name = "tmd_engine_vertical_spring_nogravity.png"
    real_name = "tmd_engine_real.jpg"

    m_E = Symbol("m_E", positive=True)
    k_E = Symbol("k_E", positive=True)
    z = dynamicsymbols("z")
    z_E = dynamicsymbols("z_E")

    def __init__(
        self, m_E=None, k_E=None, z_E=None, z=None, ivar=Symbol("t"), **kwargs
    ):

        if m_E is not None:
            self.m_E = m_E
        if k_E is not None:
            self.k_E = k_E
        if z is not None:
            self.z = z
        if z_E is not None:
            self.z_E = z_E
        self.ivar = ivar

        self.qs = [self.z_E, self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(
            self.m_E, pos1=self.z_E + self.z, pos2=self.z, qs=self.qs
        )(label="Material point - mass of the TMD")
        self._spring = Spring(
            self.k_E, pos1=self.z_E + self.z, pos2=self.z, qs=self.qs
        )(label="Spring - stiffness of the spring")

        components["_material_point"] = self._material_point
        components["_spring"] = self._spring

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict
