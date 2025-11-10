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
# ureg = units

from functools import cached_property, lru_cache


# Patryk
# dane domyślne i numeryczne
class SpringMassSystem(ComposedSystem):
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

    scheme_name = "harmonic_oscillator.png"
    real_name = "leaf_spring.png"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    ivar = Symbol("t")

    z = dynamicsymbols("z")

    def __init__(self, m=None, k=None, z=None, ivar=Symbol("t"), **kwargs):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if z is not None:
            self.z = z
        self.ivar = ivar

        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring = Spring(self.k, self.z, qs=self.qs)

        components["material_point"] = self.material_point
        components["spring"] = self.spring

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict


class ThreeSprings(ComposedSystem):

    scheme_name = "engine.png"
    real_name = "engine_real.PNG"

    k_l = Symbol("k_l", positive=True)
    k_cl = Symbol("k_cl", positive=True)
    ivar = Symbol("t")

    pos1 = dynamicsymbols("pos1")
    pos2 = dynamicsymbols("pos2")
    qs = dynamicsymbols("pos1 pos2")

    def __init__(
        self,
        k_l=None,
        k_cl=None,
        pos1=None,
        pos2=None,
        qs=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if k_l is not None:
            self.k_l = k_l
        if k_cl is not None:
            self.k_cl = k_cl
        if pos1 is not None:
            self.pos1 = pos1
        if pos2 is not None:
            self.pos2 = pos2

        self.ivar = ivar

        if qs is None:
            self.qs = [self.pos1, self.pos2]
        else:
            self.qs = qs

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._spring1 = Spring(self.k_l, pos1=self.pos1, pos2=self.pos2, qs=self.qs)
        self._spring2 = Spring(self.k_l, pos1=self.pos1, pos2=self.pos2, qs=self.qs)
        self._spring3 = Spring(self.k_cl, pos1=self.pos1, pos2=self.pos2, qs=self.qs)

        components["spring1"] = self._spring1
        components["spring2"] = self._spring2
        components["spring3"] = self._spring3

        return components


class ForcedTrolleyWithThreeSprings(ComposedSystem):

    scheme_name = "engine.png"
    real_name = "engine_real.PNG"

    m = Symbol("m", positive=True)
    k_l = Symbol("k_l", positive=True)
    k_cl = Symbol("k_cl", positive=True)
    F = Symbol("F", positive=True)
    ivar = Symbol("t")

    z = dynamicsymbols("z")

    def __init__(
        self,
        m=None,
        k_l=None,
        k_cl=None,
        F=None,
        z=None,
        qs=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if k_l is not None:
            self.k_l = k_l
        if k_cl is not None:
            self.k_cl = k_cl
        if F is not None:
            self.F = F
        if z is not None:
            self.z = z
        self.ivar = ivar

        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self._three_springs = ThreeSprings(
            self.k_l, self.k_cl, pos1=self.z, pos2=0, qs=self.qs
        )
        self._force = Force(self.F, pos1=self.z, qs=self.qs)

        components["material_point"] = self.material_point
        components["_three_springs"] = self._three_springs
        components["_force"] = self._force

        return components


# dane domyślne i numeryczne
class ForcedSpringMassSystem(SpringMassSystem):
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

    scheme_name = "forced_harmonic_oscillator_gravity.png"
    real_name = "leaf_spring.png"

    F = Symbol("F", positive=True)
    g = Symbol("g", positive=True)

    def __init__(self, m=None, k=None, F=None, g=None, z=None, ivar=None, **kwargs):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if F is not None:
            self.F = F
        if g is not None:
            self.g = g

        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}  # przerobić na kompozycję

        self._undamped_trolley = SpringMassSystem(self.m, self.k, self.z)
        self.force = Force(self.F, self.z, qs=self.qs)
        self.gravitational_force = GravitationalForce(
            self.m, self.g, self.z, qs=self.qs
        )

        components["_undamped_trolley"] = self._undamped_trolley
        components["force"] = self.force
        components["gravitational_force"] = self.gravitational_force

        return components

    def get_numerical_parameters(self):

        m, k, g, F = self.m, self.k, self.g, self.F

        get_numerical_parameters = {
            self.m: 0.5,
            self.k: 20,
            # self.F: 2,
            self.F: 2 * sin(10 * self.ivar),
            self.g: 9.81,
        }
        return get_numerical_parameters

    def symbols_description(self):

        self.sym_desc_dict = {
            self.m: r"Mass of system on the spring",
            self.k: r"Spring coefficient ",
            self.F: r"Force",
            self.g: r"Gravitational acceleration",
        }
        return {**super().symbols_description(), **self.sym_desc_dict}

    def unit_dict(self):

        from sympy.physics import units

        ureg = UnitRegistry()

        unit_dict = {
            self.k: ((ureg.kilogram * ureg.meter) / ureg.second / ureg.second)
            / ureg.meter,
            self.m: ureg.kilogram,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.F: (ureg.kilogram * ureg.meter) / ureg.second / ureg.second,
        }
        return unit_dict


# dziedziczenie po SpringMass
# dane domyślne i numeryczne
class SpringDamperMassSystem(ComposedSystem):
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

    scheme_name = "forced_damped_harmonic_oscillator.png"
    real_name = "leaf_spring.png"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    F = Symbol("F", positive=True)
    c = Symbol("c", positive=True)
    ivar = Symbol("t")

    z = dynamicsymbols("z")

    def __init__(self, m=None, k=None, F=None, z=None, c=None, ivar=None, **kwargs):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if F is not None:
            self.F = F
        if z is not None:
            self.z = z
        if c is not None:
            self.c = c
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring = Spring(self.k, self.z, qs=self.qs)
        self.force = Force(self.F, self.z, qs=self.qs)
        self.damper = Damper(self.c, self.z, qs=self.qs)

        components["material_point"] = self.material_point
        components["spring"] = self.spring
        components["force"] = self.force
        components["damper"] = self.damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict


# można spróbować dziedizczyć bo SPringMass
class SprungTrolley(ComposedSystem):
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
    k_l = Symbol("k_l", positive=True)
    k_r = Symbol("k_r", positive=True)
    ivar = Symbol("t")
    x = dynamicsymbols("x")

    def __init__(self, m=None, k_l=None, k_r=None, x=None, ivar=None, **kwargs):

        if m is not None:
            self.m = m
        if k_l is not None:
            self.k_l = k_l
        if k_r is not None:
            self.k_r = k_r
        if ivar is not None:
            self.ivar = ivar
        if x is not None:
            self.x = x

        self.qs = [self.x]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, pos1=self.x, qs=self.qs)("Trolley")
        self.spring_l = Spring(self.k_l, pos1=self.x, qs=self.qs)("Left joint")
        self.spring_r = Spring(self.k_r, pos1=self.x, qs=self.qs)("Right joint")

        components["material_point"] = self.material_point
        components["spring_l"] = self.spring_l
        components["spring_r"] = self.spring_r

        return components

    def get_default_data(self):

        m0, k0 = self.m0, self.k0

        default_data_dict = {
            self.m: [S.One * no * m0 / 100 for no in range(80, 120)],  # percentage
            self.k_r: [S.One * no * k0 / 100 for no in range(80, 120)],  # percentage
            self.k_l: [S.One * no * k0 / 100 for no in range(70, 110)],  # percentage
        }
        return default_data_dict

    def get_numerical_data(self):

        m0, k0 = 100, 1000

        default_data_dict = {
            self.m: [S.One * no * m0 / 100 for no in range(80, 120)],  # percentage
            self.k_r: [S.One * no * k0 / 100 for no in range(80, 120)],  # percentage
            self.k_l: [S.One * no * k0 / 100 for no in range(70, 110)],  # percentage
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k_l: r"Spring left coefficient ",  # ewentualnie zmienić na spring 1,2 coefficient
            self.k_r: r"Spring right coefficient ",
            self.x: r"displacement",
        }

        return self.sym_desc_dict


class NonLinearTrolley(NonlinearComposedSystem):

    scheme_name = "nonlin_trolley.PNG"
    real_name = "nonlin_trolley_real.PNG"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    d = Symbol("d", positive=True)
    l_0 = Symbol("l_0", positive=True)
    x = dynamicsymbols("x")

    def __init__(
        self, m=None, k=None, d=None, l_0=None, x=None, ivar=Symbol("t"), **kwargs
    ):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if d is not None:
            self.d = d
        if l_0 is not None:
            self.l_0 = l_0
        if x is not None:
            self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        # self._disk = Disk(S.One/2 * self.m1*self.R**2,pos1 = self.x / self.R, qs=[self.x]) ??? returns diffrend eoms than soultion above
        self._spring = Spring(
            self.k, pos1=sqrt(self.x**2 + self.d**2), pos2=-self.l_0, qs=[self.x]
        )

        components["trolley"] = self._trolley
        components["spring"] = self._spring

        return components

    #     def get_default_data(self):

    #         m0, k0, l0= symbols('m_0 k_0 l_0', positive=True)

    #         default_data_dict = {
    #             self.m: [S.One * m0 * no for no in range(10,20)],
    #             self.d: [S.One * l0 * no for no in range(1,10)],
    #             self.k: [S.One * k0 * no for no in range(20,30)],
    #             self.l_0: [S.One * l0 * no for no in range(10,15)],
    #         }

    #         return default_data_dict

    def get_numerical_data(
        self,
    ):  ### Brakowało numerical data. Przez to komenda get numerical parameters nie działa

        m0, k0, l0, c0 = symbols("m_0 k_0 l_0 c_0", positive=True)

        default_data_dict = {
            self.m1: [
                0.5 * m0,
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
            ],
            self.d: [
                5 * l0,
                2 * l0,
                3 * S.Half * l0,
                4 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],
            self.kl: [
                1 * k0,
                3 * k0,
                2 * k0,
                4 * k0,
                5 * k0,
                6 * k0,
                7 * k0,
                8 * k0,
                9 * k0,
            ],
            self.l_0: [
                1 * l0,
                3 * l0,
                2 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],  ### Brakowało nawiasu
            self.c: [
                1 * c0,
                3 * c0,
                2 * c0,
                4 * c0,
                5 * c0,
                6 * c0,
                7 * c0,
                8 * c0,
                9 * c0,
            ],
        }

        return default_data_dict


# Zrobione Amadi
class TrolleyWithPendulum(ComposedSystem):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "taipei101.png"

    l = Symbol("l", positive=True)
    m_t = Symbol("m_trolley", positive=True)
    m_p = Symbol("m_pendulum", positive=True)
    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("varphi")
    x = dynamicsymbols("x")
    Ek = Symbol("E_k", positive=True)
    Ep = Symbol("E_p", positive=True)

    def __init__(
        self,
        l=None,
        m_t=None,
        m_p=None,
        k=None,
        g=None,
        Omega=None,
        phi=None,
        x=None,
        F=None,
        Ek=None,
        Ep=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if l is not None:
            self.l = l
        if m_t is not None:
            self.m_t = m_t
        if m_p is not None:
            self.m_p = m_p
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if Ek is not None:
            self.Ek = Ek
        if Ep is not None:
            self.Ep = Ep
        self.ivar = ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        self._pendulum = PendulumKinematicExct(
            self.l, self.m_p, self.g, self.phi, self.x, self.ivar
        )(label="Pendulum")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.phi]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_pendulum"] = self._pendulum
        components["_force"] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.phi: r"Angle of a pendulum",
            self.m_p: r"Mass of pendulum",
            self.m_t: r"Mass of trolley",
            self.g: r"Gravity constant",
            self.F: r"Force",
            self.Ek: r"Kinetic Energy of the system",
            self.Ep: r"Potential Energy of the system",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols("m_0 l_0 F_0 Omega_0 k_0", positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * self.Omega],
            self.g: [S.One * self.g],
            self.k: [S.One * self.g * m0 / l0 * no for no in range(1, 30)],
            self.x: [self.x],
        }

        # default_data_dict.update({self.k: S.One * self.g * default_data_dict[self.m_t][0]+default_data_dict[self.l][0],
        #                         })

        return default_data_dict

    def get_numerical_data(self):

        m_taipei = 700 * 1e3
        omg_nat = 2 * pi / 10

        default_data_dict = {
            self.m_t: [m_taipei for no in range(20, 30)],
            self.m_p: [0.01 * m_taipei for no in range(1, 10)],
            self.l: [9.81 / omg_nat**2 for no in range(1, 10)],
            self.F: [300 * 1e3 for no in range(50, 100)],
            self.Omega: [omg_nat for no in range(1, 6)],
            self.k: [m_taipei * (omg_nat) ** 2 for no in range(50, 100)],
            self.g: [9.81],
        }
        return default_data_dict

    @lru_cache
    def max_static_cable_force(self):
        data = self._given_data
        ans = self.force_in_cable()
        free_coeff = ans.subs(
            {
                (cos(self.Omega * self.ivar)) ** 2: 0,
                (sin(self.Omega * self.ivar)) ** 2: 0,
            }
        ).subs(data)
        # display(free_coeff)
        return (abs(free_coeff)).subs(data)

    @lru_cache
    def max_dynamic_cable_force(self):
        ans = self.force_in_cable()
        data = self._given_data
        sin_coeff = (ans.coeff((sin(self.Omega * self.ivar)) ** 2)).subs(data)
        # display(sin_coeff)
        cos_coeff = (ans.coeff((cos(self.Omega * self.ivar)) ** 2)).subs(data)
        # display(cos_coeff)
        free_coeff = (
            ans.subs(
                {
                    (cos(self.Omega * self.ivar)) ** 2: 0,
                    (sin(self.Omega * self.ivar)) ** 2: 0,
                }
            )
        ).subs(data)
        # display(free_coeff)
        return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)

    @lru_cache
    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    @lru_cache
    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    @lru_cache
    def force_in_cable(self):
        data = self._given_data
        dyn_sys = self  # .subs(data)
        dyn_sys_lin = dyn_sys.linearized()
        phi = dyn_sys_lin._fodes_system.steady_solution[1]

        force_in_cable = (
            self.m_p * self.g * (1 - S.One / 2 * phi**2)
            + self.m_p * self.l * phi.diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(data)
        # force_subs=force_in_cable.subs({self.Omega:sqrt(self.g/self.l)}).subs(data)
        # display(force_subs.doit().expand())
        return force_subs.doit().expand()

    @property
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST,
            mech_comp.StaticAndDynamicCableForceComponent,
            mech_comp.StaticCableDiameter,
            mech_comp.DynamicCableDiameter,
            mech_comp.DynamicTensionForceComponent,
        ]

        return comp_list

    def unit_dict(self):

        from sympy.physics import units

        ureg = UnitRegistry()

        unit_dict = {
            self.m_t: ureg.kilogram,
            self.m_p: ureg.kilogram,
            self.l: ureg.meter,
            self.F: ureg.newton,
            self.Ek: ureg.newton * ureg.meter,
            self.Ep: ureg.newton * ureg.meter,
            self.Omega: ureg.radian,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.k: ureg.newton / ureg.meter,
            self.x: ureg.meter,
            self.x.diff(self.ivar): ureg.meter / ureg.second,
            self.x.diff(self.ivar, 2): ureg.meter / ureg.second / ureg.second,
            self.phi: ureg.radian,
            self.phi.diff(self.ivar): ureg.radian / ureg.second,
            self.phi.diff(self.ivar, 2): ureg.radian / ureg.second / ureg.second,
            self.ivar: ureg.second,
        }
        return unit_dict


class TrolleyWithElasticPendulum(TrolleyWithPendulum):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "taipei101.png"

    l = Symbol("l", positive=True)
    m_t = Symbol("m_t", positive=True)
    m_p = Symbol("m_p", positive=True)
    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")
    u = dynamicsymbols("u")
    y = dynamicsymbols("y")
    v = Symbol("v", positive=True)
    k_l = Symbol("k_l", positive=True)
    u_0 = Symbol("u_0", positive=True)

    def __init__(
        self,
        l=None,
        m_t=None,
        m_p=None,
        k=None,
        g=None,
        Omega=None,
        phi=None,
        x=None,
        F=None,
        u=None,
        v=None,
        y=None,
        k_l=None,
        u_0=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if l is not None:
            self.l = l
        if m_t is not None:
            self.m_t = m_t
        if m_p is not None:
            self.m_p = m_p
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if u is not None:
            self.u = u
        if v is not None:
            self.v = v
        if y is not None:
            self.y = y
        if k_l is not None:
            self.k_l = k_l
        if u_0 is not None:
            self.u_0 = self.m_p * self.g / self.k_l
        self.ivar = ivar
        self.qs = [self.x, self.phi, self.u]
        self.u_0 = self.g * self.l / self.k_l
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        # self._pendulum = PendulumKinematicExct(self.l+self.u, self.m_p, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.qs[0], qs=self.qs
        )(label="Force")
        self._material_point_1 = MaterialPoint(
            self.m_p, self.x + (self.l + self.u + self.u_0) * sin(self.phi), qs=self.qs
        )
        self._material_point_2 = MaterialPoint(
            self.m_p, (self.l + self.u + self.u_0) * cos(self.phi), qs=self.qs
        )
        self._gravity = GravitationalForce(
            self.m_p, self.g, -(self.l + self.u + self.u_0) * cos(self.phi), qs=self.qs
        )
        self._spring_2 = Spring(self.k_l, self.u, qs=self.qs)  # zostaw xD tu jest ok

        components["_trolley"] = self._trolley
        # components['_pendulum'] = self._pendulum
        components["_force"] = self._force
        components["_material_point_1"] = self._material_point_1
        components["_material_point_2"] = self._material_point_2
        # components['_spring_1'] = self._spring_1
        components["_spring_2"] = self._spring_2
        components["_gravity"] = self._gravity
        return components

    def symbols_description(self, lang="en"):
        if lang == "pl":
            self.sym_desc_dict = {
                self.m_t: r"masa początkowa wózka",
                self.m_p: r"masa początkowa wahadła",
                self.k: r"sztywność sprężyny mocującej",
                self.k_l: r"sztywność wahadła",
                self.l: r"długość wahadła",
                self.Omega: r"częstotliwość wymuszenia",
                self.F: r"wartość siły wymuszającej",
                self.g: r"przyspieszenie ziemskie",
                self.flow_coeff: r"współczynnik przepływu cieczy",
                self.t_init: r"czas aktywacji tłumienia",
            }
            return self.sym_desc_dict

        else:
            self.sym_desc_dict = {
                self.m_t: r"Mass of trolley",
                self.m_p: r"pendulum initial mass",
                self.k: r"mounting spring stiffness",
                self.k_l: r"pendulum stiffness",
                self.l: r"pendulum length",
                self.Omega: r"excitation frequency",
                self.F: r"excitation force",
                self.g: r"gravity constant",
                self.flow_coeff: r"flow coefficient",
                self.t_init: r"damping activation time",
            }
            return self.sym_desc_dict



class DampedTrolleyWithPendulumVariableInertia(TrolleyWithElasticPendulum):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol("l", positive=True)
    k_l = Symbol("k_l", positive=True)
    m_t = Symbol("m_t", positive=True)
    m_p = Symbol("m_p", positive=True)
    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")
    c = Symbol("c", positive=True)
    c_p = Symbol("c_p", positive=True)
    u = dynamicsymbols("u", positive=True)
    v = Symbol("v", positive=True)
    alpha = Symbol("\\alpha")
    beta = Symbol("\\beta")
    rayleigh_damping_matrix = Symbol("rayleigh_damping_matrix", positive=True)

    def __init__(
        self,
        l=None,
        k_l=None,
        m_t=None,
        m_p=None,
        k=None,
        g=None,
        alpha=None,
        beta=None,
        Omega=None,
        phi=None,
        x=None,
        rayleigh_damping_matrix=None,
        F=None,
        c=None,
        c_p=None,
        u=None,
        v=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if l is not None:
            self.l = l
        if k_l is not None:
            self.k_l = k_l
        if m_t is not None:
            self.m_t = m_t
        if m_p is not None:
            self.m_p = m_p
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if rayleigh_damping_matrix is not None:
            self.rayleigh_damping_matrix = rayleigh_damping_matrix
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if c is not None:
            self.c = c
        if c_p is not None:
            self.c_p = c_p
        if u is not None:
            self.u = u
        if v is not None:
            self.v = v
        self.ivar = ivar
        # self.qs = [self.x,self.phi,self.u]
        self.qs = [self.x, self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        # self.u=self.l+self.v*self.ivar
        # self._TrolleyWithPendulum = TrolleyWithElasticPendulum(self.l+self.qs[2], self.m_t, self.m_p, self.k, self.g, self.Omega, self.qs[2], self.qs[0], self.F,ivar=self.ivar)(label='Trolley')
        self._TrolleyWithPendulum = TrolleyWithElasticPendulum(
            self.l,
            self.m_t,
            self.m_p,
            self.k,
            self.g,
            self.Omega,
            F=self.F,
            ivar=self.ivar,
        )(label="Trolley")
        self._trolley_damper = Damper(
            self.c, pos1=self.x + (self.l + self.u) * sin(self.phi), qs=self.qs
        )(label="Trolley damper")
        self._pendulum_damper = Damper(self.c_p, pos1=self.phi, qs=self.qs)(
            label="Pendulum damper"
        )

        components["_TrolleyWithPendulum"] = self._TrolleyWithPendulum
        components["_trolley_damper"] = self._trolley_damper
        components["_pendulum_damper"] = self._pendulum_damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.k_l: r"Stiffness of a pendulum",
            self.x: r"Kinematic lateral excitation",
            self.phi: r"Angle of a pendulum",
            self.m_p: r"Mass of pendulum",
            self.m_t: r"Mass of trolley",
            self.alpha: r"inertia matrix Rayleigh damping coefficient",
            self.beta: r"stiffness matrix Rayleigh damping coefficient",
            self.g: r"Gravity constant",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    # def get_default_data(self):

    # return default_data_dict

    # def get_numerical_data(self):

    # return default_data_dict

    def symbols_description(self, lang="en"):

        if lang == "pl":

            self.sym_desc_dict = {
                self.m_t: r"masa początkowa wózka",
                self.m_p: r"masa początkowa wahadła",
                self.k: r"sztywność sprężyny mocującej",
                self.k_l: r"sztywność wahadła",
                self.l: r"długość wahadła",
                self.Omega: r"częstotliwość wymuszenia",
                self.alpha: r"współczynnik tłumienia Rayleigha przy macierzy bezwładności",
                self.beta: r"współczynnik tłumienia Rayleigha przy macierzy sztywności",
                self.F: r"wartość siły wymuszającej",
                self.g: r"przyspieszenie ziemskie",
                self.flow_coeff: r"współczynnik przepływu cieczy",
                self.t_init: r"czas aktywacji tłumienia",
            }
            return self.sym_desc_dict

        else:

            self.sym_desc_dict = {
                self.m_t: r"trolley initial mass",
                self.m_p: r"pendulum initial mass",
                self.k: r"mounting spring stiffness",
                self.k_l: r"pendulum stiffness",
                self.l: r"pendulum length",
                self.alpha: r"inertia matrix Rayleigh damping coefficient",
                self.beta: r"stiffness matrix Rayleigh damping coefficient",
                self.Omega: r"excitation frequency",
                self.F: r"excitation force",
                self.g: r"gravity constant",
                self.flow_coeff: r"flow coefficient",
                self.t_init: r"damping activation time",
            }
            return self.sym_desc_dict


# Zrobione Amadi
class DampedTrolleyWithPendulum(TrolleyWithPendulum):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol("l", positive=True)
    m_t = Symbol("m_trolley", positive=True)
    m_p = Symbol("m_pendulum", positive=True)
    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")
    c = Symbol("c", positive=True)

    def __init__(
        self,
        l=None,
        m_t=None,
        m_p=None,
        k=None,
        g=None,
        Omega=None,
        phi=None,
        x=None,
        F=None,
        c=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if l is not None:
            self.l = l
        if m_t is not None:
            self.m_t = m_t
        if m_p is not None:
            self.m_p = m_p
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if c is not None:
            self.c = c
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self._TrolleyWithPendulum = TrolleyWithPendulum(
            self.l,
            self.m_t,
            self.m_p,
            self.k,
            self.g,
            self.Omega,
            self.phi,
            self.x,
            self.F,
            self.ivar,
        )(label="Trolley")
        self._trolley_damper = Damper(self.c, pos1=self.x, qs=[self.x, self.phi])(
            label="Trolley damper"
        )
        self._pendulum_damper = Damper(self.c, pos1=self.phi, qs=[self.x, self.phi])(
            label="Pendulum damper"
        )

        components["_TrolleyWithPendulum"] = self._TrolleyWithPendulum
        components["_trolley_damper"] = self._trolley_damper
        components["_pendulum_damper"] = self._pendulum_damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.phi: r"Angle of a pendulum",
            self.m_p: r"Mass of pendulum",
            self.m_t: r"Mass of trolley",
            self.g: r"Gravity constant",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
            self.c: r"Damping coefficient",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0, lam, lam0 = symbols(
            "m_0 l_0 F_0 Omega_0 k_0 lambda lambda_0", positive=True
        )

        default_data_dict = {
            self.c: [S.One * self.k * lam],
            lam: [S.One * lam0 / 1000 * no for no in range(1, 10)],
        }
        return default_data_dict

    # def get_numerical_data(self):

    #     default_data_dict = {
    #         self.c: [no for no in range(5, 25)]
    #     }
    #     return default_data_dict

    def get_numerical_data(self):

        m_taipei = 700 * 1e3
        omg_nat = 2 * pi / 10

        default_data_dict = {
            self.m_t: [m_taipei for no in range(20, 30)],
            self.m_p: [0.01 * m_taipei for no in range(1, 10)],
            self.l: [9.81 / omg_nat**2 for no in range(1, 10)],
            self.F: [300 * 1e3 for no in range(50, 100)],
            self.Omega: [omg_nat for no in range(1, 6)],
            self.k: [m_taipei * (omg_nat) ** 2 for no in range(50, 100)],
            self.g: [9.81],
            self.c: [0.1 * m_taipei * (omg_nat) ** 2 for no in range(50, 100)],
        }
        return default_data_dict

    def unit_dict(self):

        from sympy.physics import units

        ureg = UnitRegistry()

        unit_dict = {
            self.l: ureg.meter,
            self.k: ((ureg.kilogram * ureg.meter) / ureg.second / ureg.second)
            / ureg.meter,
            self.x: ureg.meter,
            self.phi: ureg.radian,
            self.m_p: ureg.kilogram,
            self.m_t: ureg.kilogram,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.F: (ureg.kilogram * ureg.meter) / ureg.second / ureg.second,
            self.Omega: ureg.radian,
            self.c: ureg.kilogram / ureg.second,
        }
        return unit_dict


### Nieliniowe
class ForcedNonLinearTrolley(NonlinearComposedSystem):

    scheme_name = "nonlin_trolley.PNG"
    real_name = "nonlin_trolley_real.PNG"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    d = Symbol("d", positive=True)
    l_0 = Symbol("l_0", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    x = dynamicsymbols("x")

    def __init__(
        self,
        m=None,
        k=None,
        d=None,
        l_0=None,
        Omega=None,
        F=None,
        x=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if d is not None:
            self.d = d
        if l_0 is not None:
            self.l_0 = l_0
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if x is not None:
            self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(
            self.k, pos1=sqrt(self.x**2 + self.d**2), pos2=-self.l_0, qs=[self.x]
        )
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x]
        )

        components["trolley"] = self._trolley
        components["spring"] = self._spring
        components["force"] = self._force

        return components

    def get_default_data(self):

        m0, k0, l0 = symbols("m_0 k_0 l_0", positive=True)

        default_data_dict = {
            self.m: [
                0.5 * m0,
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
            ],
            self.d: [
                5 * l0,
                2 * l0,
                3 * S.Half * l0,
                4 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],
            self.k: [
                1 * k0,
                3 * k0,
                2 * k0,
                4 * k0,
                5 * k0,
                6 * k0,
                7 * k0,
                8 * k0,
                9 * k0,
            ],
            self.l_0: [
                1 * l0,
                3 * l0,
                2 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],
        }

        return default_data_dict

    def get_numerical_data(
        self,
    ):  ### Brakowało numerical data. Przez to komenda get numerical parameters nie działa

        m0, k0, l0, c0 = symbols("m_0 k_0 l_0 c_0", positive=True)

        default_data_dict = {
            self.m1: [
                0.5 * m0,
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
            ],
            self.d: [
                5 * l0,
                2 * l0,
                3 * S.Half * l0,
                4 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],
            self.kl: [
                1 * k0,
                3 * k0,
                2 * k0,
                4 * k0,
                5 * k0,
                6 * k0,
                7 * k0,
                8 * k0,
                9 * k0,
            ],
            self.l_0: [
                1 * l0,
                3 * l0,
                2 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
            ],  ### Brakowało nawiasu
            self.c: [
                1 * c0,
                3 * c0,
                2 * c0,
                4 * c0,
                5 * c0,
                6 * c0,
                7 * c0,
                8 * c0,
                9 * c0,
            ],
        }

        return default_data_dict


#     def subs(self, *args, **kwargs):

#         return super().subs(*args, **kwargs).subs(*args, **kwargs)


class ForcedTrolleyWithSpring(ComposedSystem):  ### 1 ODE

    scheme_name = "nonlin_trolley.PNG"
    real_name = "nonlin_trolley_real.PNG"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    G = Symbol("G", positive=True)
    x = dynamicsymbols("x")

    m0 = Symbol("m_0", positive=True)
    k0 = Symbol("k_0", positive=True)
    Omega0 = Symbol("Omega_0", positive=True)
    F0 = Symbol("F_0", positive=True)

    def __init__(
        self,
        m=None,
        k=None,
        Omega=None,
        F=None,
        m0=None,
        k0=None,
        Omega0=None,
        F0=None,
        G=None,
        x=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if G is not None:
            self.G = G
        if x is not None:
            self.x = x

        if m0 is not None:
            self.m0 = m0
        if k0 is not None:
            self.k0 = k0
        if Omega0 is not None:
            self.Omega0 = Omega0
        if F0 is not None:
            self.F0 = F0

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(self.k, pos1=self.x, qs=[self.x])
        self._force = Force(
            self.F * sin(self.Omega * self.ivar) + self.G, pos1=self.x, qs=[self.x]
        )

        components["trolley"] = self._trolley
        components["spring"] = self._spring
        components["force"] = self._force

        return components

    def get_default_data(self):

        m0, k0, Omega0, F0 = self.m0, self.k0, self.Omega0, self.F0

        default_data_dict = {
            self.m: [100 * m0],
            self.k: [50 * k0],
            self.Omega: [
                0.5 * 3.14 * Omega0,
                1 * 3.14 * Omega0,
                2 * 3.14 * Omega0,
                4 * 3.14 * Omega0,
            ],
            self.F: [0.5 * 100 * F0, 1 * 100 * F0, 2 * 100 * F0, 4 * 100 * F0],
        }
        default_data_dict.update(
            {
                self.G: [
                    4
                    * default_data_dict[self.F][0]
                    * cos(0.5 * default_data_dict[self.Omega][0] * self.ivar),
                    default_data_dict[self.F][0]
                    * cos(0.75 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                    1.5
                    * default_data_dict[self.F][0]
                    * 2
                    * cos(1.25 * default_data_dict[self.Omega][0] * self.ivar),
                    3
                    * default_data_dict[self.F][0]
                    * cos(2 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                ]
            }
        )

        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [100],
            self.k: [50],
            self.Omega: [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F: [0.5 * 100, 1 * 100, 2 * 100, 4 * 100],
        }
        default_data_dict.update(
            {
                self.G: [
                    4
                    * default_data_dict[self.F][0]
                    * cos(0.5 * default_data_dict[self.Omega][0] * self.ivar),
                    default_data_dict[self.F][0]
                    * cos(0.75 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                    1.5
                    * default_data_dict[self.F][0]
                    * 2
                    * cos(1.25 * default_data_dict[self.Omega][0] * self.ivar),
                    3
                    * default_data_dict[self.F][0]
                    * cos(2 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                ]
            }
        )

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"masa wózka",
            self.k: r"sztywność sprężyny",
            self.Omega: r"częstotliwość wymuszenia",
            self.F: r"wartość stałej siły wymuszającej",
            self.G: r"wartość zmiennej siły wymuszającej",
            self.x: "przemieszczenie wózka nr 1",
        }
        return self.sym_desc_dict

    def unit_dict(self):
        unit_dict = {
            self.m: ureg.kilogram,
            self.k: ureg.newton / ureg.meter,
            self.Omega: ureg.hertz,
            self.F: ureg.newton,
            self.G: ureg.newton,
            self.x: ureg.meter,
        }
        return unit_dict

    @property
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST
            # mech_comp.FundamentalMatrixComponent,
            # mech_comp.GeneralSolutionComponent,
            # mech_comp.SteadySolutionComponent,
        ]

        return comp_list


class ForcedFTrolleyWithSpring(ComposedSystem):  ### 2 ODE

    scheme_name = "nonlin_trolley.PNG"
    real_name = "nonlin_trolley_real.PNG"

    m = Symbol("m", positive=True)
    k = Symbol("k", positive=True)
    c = Symbol("c", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    G = Symbol("G", positive=True)
    x = dynamicsymbols("x")

    def __init__(
        self,
        m=None,
        k=None,
        c=None,
        Omega=None,
        F=None,
        G=None,
        x=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if c is not None:
            self.c = c
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if G is not None:
            self.G = G
        if x is not None:
            self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(self.k, pos1=self.x, qs=[self.x])
        self._damper = Damper(self.c, pos1=self.x, qs=[self.x])
        self._force = Force(
            self.F * sin(self.Omega * self.ivar) + self.G, pos1=self.x, qs=[self.x]
        )

        components["trolley"] = self._trolley
        components["spring"] = self._spring
        components["damper"] = self._damper
        components["force"] = self._force

        return components

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [100],
            self.k: [50],
            self.c: [100],
            self.Omega: [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F: [0.5 * 100, 1 * 100, 2 * 100, 4 * 100],
        }
        default_data_dict.update(
            {
                self.G: [
                    4
                    * default_data_dict[self.F][0]
                    * cos(0.5 * default_data_dict[self.Omega][0] * self.ivar),
                    default_data_dict[self.F][0]
                    * cos(0.75 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                    1.5
                    * default_data_dict[self.F][0]
                    * 2
                    * cos(1.25 * default_data_dict[self.Omega][0] * self.ivar),
                    3
                    * default_data_dict[self.F][0]
                    * cos(2 * default_data_dict[self.Omega][0] * self.ivar) ** 2,
                ]
            }
        )

        return default_data_dict


class VariableMassTrolleyWithPendulum(ComposedSystem):

    scheme_name = "kin_exct_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol("l", positive=True)

    m_f = Symbol("m_f", positive=True)

    m_t = Symbol("m_t", positive=True)
    m_t0 = Symbol("m_t0", positive=True)
    m_tf = Symbol("m_tf", positive=True)

    m_p = Symbol("m_p", positive=True)
    m_p0 = Symbol("m_p0", positive=True)
    m_pf = Symbol("m_pf", positive=True)

    flow_coeff = Symbol("\\lambda", positive=True)

    t_0 = Symbol("t_0", positive=True)

    m_0 = Symbol("m_0", positive=True)

    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)

    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")

    def __init__(
        self,
        l=None,
        m_f=None,
        m_t=None,
        m_t0=None,
        m_tf=None,
        m_p=None,
        m_p0=None,
        m_pf=None,
        m_0=None,
        flow_coeff=None,
        t_0=None,
        k=None,
        g=None,
        Omega=None,
        phi=None,
        x=None,
        F=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if l is not None:
            self.l = l
        if m_f is not None:
            self.m_f = m_f
        if m_t is not None:
            self.m_t = m_t
        if m_t0 is not None:
            self.m_t0 = m_t0
        if m_tf is not None:
            self.m_tf = m_tf
        if m_p is not None:
            self.m_p = m_p
        if m_p0 is not None:
            self.m_p0 = m_p0
        if m_pf is not None:
            self.m_pf = m_pf
        if m_0 is not None:
            self.m_0 = m_0
        if flow_coeff is not None:
            self.flow_coeff = flow_coeff
        if t_0 is not None:
            self.t_0 = t_0
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        self.ivar = ivar
        self._init_from_components(**kwargs)

        self.trans_expr = (
            S.One / 2 - atan(self.flow_coeff * (self.ivar - self.t_0)) / pi
        )

        # self.m_tf = self.m_f#((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        # self.m_pf = self.m_f#((S.One/2+atan(self.flow_coeff*(self.ivar-self.t0))/pi))

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(
            self.m_t + self.m_tf, self.k, self.x, self.ivar
        )(label="Trolley")
        self._pendulum = PendulumKinematicExct(
            self.l, self.m_p + self.m_pf, self.g, self.phi, self.x, self.ivar
        )(label="Pendulum")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.phi]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_pendulum"] = self._pendulum
        components["_force"] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.phi: r"Angle of a pendulum",
            self.m_p: r"Mass of pendulum",
            self.m_t: r"Mass of trolley",
            self.g: r"Gravity constant",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols("m_0 l_0 F_0 Omega_0 k_0", positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.g: [S.One * self.g],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m_p: [10],
            self.m_f: [50],
            #             self.m_tf: [50],
            #             self.m_pf: [50],
            self.F: [5],
            self.flow_coeff: [5],
            self.t_0: [10],
            self.Omega: [3.14 * 0.5],
            self.g: [9.81],
            # self.k: [2]
        }
        default_data_dict.update(
            {
                self.m_t: [10.0 * default_data_dict[self.m_p][0]],
                # self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
                # self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)]
            }
        )
        return default_data_dict

    def symbols_description(self, lang="en"):

        if lang == "pl":

            self.sym_desc_dict = {
                self.m_t: r"masa początkowa wózka",
                self.m_p: r"masa początkowa wahadła",
                self.m_tf: r"masa transportowanej cieczy wypływającej z wózka",
                self.m_pf: r"masa transportowanej cieczy wypływającej do wahadła",
                self.m_f: r"masa transferowanej cieczy",
                self.k: r"sztywność sprężyny mocującej",
                self.l: r"długość wahadła",
                self.Omega: r"częstotliwość wymuszenia",
                self.F: r"wartość siły wymuszającej",
                self.g: r"przyspieszenie ziemskie",
                self.flow_coeff: r"współczynnik przepływu cieczy",
                self.t_init: r"czas aktywacji tłumienia",
            }
            return self.sym_desc_dict

        else:

            self.sym_desc_dict = {
                self.m_t: r"trolley initial mass",
                self.m_p: r"pendulum initial mass",
                self.m_tf: r"fluid mass transported from the trolley",
                self.m_pf: r"fluid mass transported to the pendulum",
                self.m_f: r"transferred fluid mass",
                self.k: r"mounting spring stiffness",
                self.l: r"pendulum length",
                self.Omega: r"excitation frequency",
                self.F: r"excitation force",
                self.g: r"gravity constant",
                self.flow_coeff: r"flow coefficient",
                self.t_init: r"damping activation time",
            }
            return self.sym_desc_dict

    def _to_msm_form(self):

        M_p = Symbol("M_p", positive=True)
        M_t = Symbol("M_t", positive=True)
        Omega1 = Symbol("Omega1", positive=True)
        Omega2 = Symbol("Omega2", positive=True)

        subs_dict = {self.m_p + self.m_pf: M_p, self.m_t + self.m_tf: M_t}

        odes_lhs = self.linearized()._ode_system.lhs

        ode1 = (
            (odes_lhs.subs(subs_dict)[0] / (M_p + M_t))
            .expand()
            .doit()
            .ratsimp()
            .subs(self.k, Omega1**2 * (M_p + M_t))
            .ratsimp()
            .expand()
            .collect([self.phi.diff(self.ivar, 2) * self.l])
        )
        ode2 = (
            (odes_lhs.subs(subs_dict)[1] / M_p / self.l / self.l)
            .expand()
            .doit()
            .ratsimp()
            .subs(self.g, Omega2**2 * self.l)
            .simplify()
        )

        eps = Symbol("varepsilon")
        del1 = Symbol("delta1")
        f = Symbol("f", positive=True)

        from dynpy.solvers.linear import ODESystem

        return ODESystem(
            odes=Matrix(
                [
                    ode1.subs(
                        {
                            ode1.coeff(self.phi.diff(self.ivar, 2)): del1 / eps,
                            self.F: f * (M_p + M_t),
                        }
                    ),
                    ode2.subs(1 / self.l, eps),
                ]
            ),
            dvars=self.q,
        )


class VariableMassTrolleyWithPendulumRayleighDamping(ComposedSystem):

    scheme_name = "kin_exct_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol("l", positive=True)

    m_f = Symbol("m_f", positive=True)

    m_t = Symbol("m_t", positive=True)
    m_t0 = Symbol("m_t0", positive=True)
    m_tf = Symbol("m_tf", positive=True)
    m_tf_eq = Symbol("m_tf_eq", positive=True)

    m_p = Symbol("m_p", positive=True)
    m_p0 = Symbol("m_p0", positive=True)
    m_pf = Symbol("m_pf", positive=True)
    m_pf_eq = Symbol("m_pf_eq", positive=True)

    flow_coeff = Symbol("\\lambda", positive=True)

    t_0 = Symbol("t_0")

    m_0 = Symbol("m_0", positive=True)

    k = Symbol("k", positive=True)
    g = Symbol("g", positive=True)

    alpha = Symbol("\\alpha")
    beta = Symbol("\\beta")
    rayleigh_damping_matrix = Symbol("rayleigh_damping_matrix", positive=True)

    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")

    b = Symbol("b", positive=True)

    omega = Symbol("omega", positive=True)
    Ek = Symbol("E_k", positive=True)
    Ep = Symbol("E_p", positive=True)

    f = Symbol("f", positive=True)

    def __init__(
        self,
        l=None,
        m_f=None,
        m_t=None,
        m_t0=None,
        m_tf=None,
        m_tf_eq=None,
        m_p=None,
        m_p0=None,
        m_pf=None,
        m_pf_eq=None,
        m_0=None,
        flow_coeff=None,
        t_0=None,
        k=None,
        g=None,
        alpha=None,
        beta=None,
        b=None,
        Omega=None,
        omega=None,
        Ek=None,
        Ep=None,
        rayleigh_damping_matrix=None,
        F=None,
        f=None,
        phi=None,
        x=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if l is not None:
            self.l = l
        if m_f is not None:
            self.m_f = m_f
        if m_t is not None:
            self.m_t = m_t
        if m_t0 is not None:
            self.m_t0 = m_t0
        if m_tf is not None:
            self.m_tf = m_tf
        if m_pf_eq is not None:
            self.m_tf_eq = m_tf_eq
        if m_p is not None:
            self.m_p = m_p
        if m_p0 is not None:
            self.m_p0 = m_p0
        if m_pf is not None:
            self.m_pf = m_pf
        if m_pf_eq is not None:
            self.m_pf_eq = m_pf_eq
        if m_0 is not None:
            self.m_0 = m_0
        if flow_coeff is not None:
            self.flow_coeff = flow_coeff
        if rayleigh_damping_matrix is not None:
            self.rayleigh_damping_matrix = rayleigh_damping_matrix
        if t_0 is not None:
            self.t_0 = t_0
        if g is not None:
            self.g = g
        if b is not None:
            self.b = b
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if omega is not None:
            self.omega = omega
        if Ek is not None:
            self.Ek = Ek
        if Ep is not None:
            self.Ep = Ep
        if F is not None:
            self.F = F
        if f is not None:
            self.f = f
        self.ivar = ivar
        self.trans_expr = (
            S.One / 2 - atan(self.flow_coeff * (self.ivar - self.t_0)) / pi
        )
        self.alpha = self.b
        self.beta = self.b / 2
        # self.rayleigh_damping_matrix = self.alpha*VariableMassTrolleyWithPendulum().inertia_matrix() + self.beta*VariableMassTrolleyWithPendulum().stiffness_matrix()
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self.m_tf_eq = self.m_f * (
            (S.One / 2 - atan(self.flow_coeff * (self.ivar - self.t_0)) / pi)
        )
        self.m_pf_eq = self.m_f * (
            (S.One / 2 + atan(self.flow_coeff * (self.ivar - self.t_0)) / pi)
        )

        self._trolley = VariableMassTrolleyWithPendulum(m_tf=self.m_tf, m_pf=self.m_pf)(
            label="Trolley"
        )

        # damping_matrix = self.alpha*self._trolley.inertia_matrix() + self.beta*self._trolley.stiffness_matrix()

        #         self._trolley_damping = Damper((damping_matrix[0]+damping_matrix[2]), pos1 = self.x , qs = [self.x,self.phi])(label='Trolley dapming')
        #         self._pendulum_damping = Damper((damping_matrix[1]+damping_matrix[3]), pos1 = self.phi , qs = [self.x,self.phi])(label='Pendulum dapming')

        self._trolley_damping = Damper(self.b, pos1=self.x, qs=[self.x, self.phi])(
            label="Trolley dapming"
        )
        self._pendulum_damping = Damper(
            self.b / 2, pos1=self.phi, qs=[self.x, self.phi]
        )(label="Pendulum dapming")

        components["_trolley"] = self._trolley
        components["_trolley_damping"] = self._trolley_damping
        components["_pendulum_damping"] = self._pendulum_damping

        return components

    #     def get_default_data(self):

    #         m0, l0, F0, Omega0, k0 = symbols('m_0 l_0 F_0 Omega_0 k_0', positive=True)

    #         default_data_dict = {
    #             self.m_t: [S.One * m0 * no for no in range(20, 30)],
    #             self.m_p: [S.One * m0 * no for no in range(1, 10)],
    #             self.l: [S.Half * l0 * no for no in range(1, 10)],
    #             self.F: [S.One * F0 * no for no in range(50, 100)],
    #             self.Omega: [S.One * Omega0],
    #             self.g: [S.One * self.g],
    #             self.k: [S.One * k0 * no for no in range(50, 100)],
    #             self.x: [self.x]
    #         }
    #         return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m_p: [10],
            self.m_f: [50],
            # self.alpha: [0.8],
            # self.beta: [2.5],
            self.b: [2.5],
            self.F: [250],
            self.flow_coeff: [10],
            self.t_0: [30],
            self.Omega: [3.14 * 0.65],
            self.g: [9.81],
            # self.k: [2]
        }
        default_data_dict.update(
            {
                self.m_t: [10.0 * default_data_dict[self.m_p][0]],
                self.l: [
                    0.99
                    * default_data_dict[self.g][0]
                    / (
                        default_data_dict[self.Omega][0]
                        * default_data_dict[self.Omega][0]
                    )
                ],
                self.k: [
                    2.0
                    * default_data_dict[self.m_p][0]
                    * default_data_dict[self.g][0]
                    / (
                        default_data_dict[self.g][0]
                        / default_data_dict[self.Omega][0] ** 2
                    )
                ],
                self.m_tf: [
                    default_data_dict[self.m_f][0]
                    * (
                        (
                            S.One / 2
                            - atan(
                                default_data_dict[self.flow_coeff][0]
                                * (self.ivar - default_data_dict[self.t_0][0])
                            )
                            / pi
                        )
                    )
                ],
                self.m_pf: [
                    default_data_dict[self.m_f][0]
                    * (
                        (
                            S.One / 2
                            + atan(
                                default_data_dict[self.flow_coeff][0]
                                * (self.ivar - default_data_dict[self.t_0][0])
                            )
                            / pi
                        )
                    )
                ],
            }
        )

        default_data_dict.update(
            {
                self.m_t
                + self.m_tf: [
                    default_data_dict[self.m_t][0] + default_data_dict[self.m_tf][0]
                ],
                self.m_p
                + self.m_pf: [
                    default_data_dict[self.m_p][0] + default_data_dict[self.m_pf][0]
                ],
            }
        )

        return default_data_dict

    def symbols_description(self, lang="en"):

        if lang == "pl":

            self.sym_desc_dict = {
                self.m_t: r"masa początkowa wózka",
                self.m_p: r"masa początkowa wahadła",
                self.m_tf: r"masa transportowanej cieczy wypływającej z wózka",
                self.m_pf: r"masa transportowanej cieczy wypływającej do wahadła",
                self.m_f: r"masa transferowanej cieczy",
                self.k: r"sztywność sprężyny mocującej",
                self.l: r"długość wahadła",
                self.Omega: r"częstotliwość wymuszenia",
                self.F: r"wartość siły wymuszającej",
                self.g: r"przyspieszenie ziemskie",
                self.flow_coeff: r"współczynnik przepływu cieczy",
                self.t_0: r"czas aktywacji tłumienia",
                self.b: r"współczynnik tłumienia",
                self.alpha: r"współczynnik tłumienia Rayleigha przy macierzy bezwładności",
                self.beta: r"współczynnik tłumienia Rayleigha przy macierzy sztywności",
                self.x: r"przemieszczenie wózka",
                self.phi: r"kąt wychylenia wahadła",
                self.omega: r"częstość własna układu",
                self.Ek: r"energia kinetyczna układu",
                self.Ep: r"energia potencjalna układu",
                self.f: r"częstotliwość",
                self.ivar: r"czas",
                diff(self.x, self.ivar): "prędkość wózka",
                diff(self.x, self.ivar, self.ivar): r"przyspieszenie wózka",
                diff(self.phi, self.ivar): r"prędkość kątowa wahadła",
                diff(self.phi, self.ivar, self.ivar): r"przyspieszenie kątowe wahadła",
            }
            return self.sym_desc_dict

        else:
            self.sym_desc_dict = {
                self.m_t: r"trolley initial mass",
                self.m_p: r"pendulum initial mass",
                self.m_tf: r"fluid mass transported from trolley",
                self.m_pf: r"fluid mass transported to pendulum",
                self.m_f: r"transferred fluid mass",
                self.k: r"mounting spring stiffness",
                self.l: r"pendulum lenght",
                self.Omega: r"excitation frequency",
                self.F: r"excitation force",
                self.g: r"gravity constant",
                self.flow_coeff: r"flow coefficient",
                self.t_0: r"damping activation time",
                self.b: r"damping coefficient",
                self.alpha: r"inertia matrix Rayleigh damping coefficient",
                self.beta: r"stiffness matrix Rayleigh damping coefficient",
                self.x: r"trolley displacement",
                self.phi: r"pendulum angular displacement",
                self.omega: r"system natural frequency",
                self.Ek: r"system kinetic energy",
                self.Ep: r"system potential energy",
                self.f: r"frequency",
                self.ivar: r"time",
                diff(self.x, self.ivar): "trolley velocity",
                diff(self.x, self.ivar, self.ivar): r"trolley acceleration",
                diff(self.phi, self.ivar): r"pendulum angular velocity",
                diff(self.phi, self.ivar, self.ivar): r"pendulum angular acceleration",
            }
            return self.sym_desc_dict


# Amadi
# klasa zastępuje MDoFTMD
class TrolleyWithTMD(ComposedSystem):
    scheme_name = "MDoF_new.png"
    real_name = "mdof_tmd_real.png"

    m = Symbol("m_trolley", positive=True)
    m_TMD = Symbol("m_TMD", positive=True)
    k = Symbol("k", positive=True)
    k_TMD = Symbol("k_TMD", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    x_e = dynamicsymbols("x_e")
    x = dynamicsymbols("x")

    def __init__(
        self,
        m=None,
        m_TMD=None,
        k=None,
        k_TMD=None,
        Omega=None,
        x_e=None,
        x=None,
        F=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if m_TMD is not None:
            self.m_TMD = m_TMD
        if m is not None:
            self.m = m
        if x_e is not None:
            self.x_e = x_e
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if k_TMD is not None:
            self.k_TMD = k_TMD
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        self.ivar = ivar

        self.qs = [self.x, self.x_e]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        self._TMD = TunedMassDamperRelativeMotion(
            self.m_TMD, self.k_TMD, self.x_e, self.x, self.ivar
        )(label="Tuned Mass Damper")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.x_e]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_TMD"] = self._TMD
        components["_force"] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.x_e: r"Angle of a pendulum",
            self.m: r"Mass of trolley",
            self.m_TMD: r"Mass of TMD",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, F0, Omega0, k0 = symbols("m_0 F_0 Omega_0 k_0", positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(20, 30)],
            self.m_TMD: [S.One * m0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [no for no in range(20, 30)],
            self.m_TMD: [no for no in range(1, 10)],
            self.F: [no for no in range(50, 100)],
            self.Omega: [3.14 * no for no in range(1, 6)],
            self.k: [no for no in range(50, 100)],
        }
        return default_data_dict

    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,  # Mass of the trolley
            self.m_TMD: ureg.kilogram,  # Mass of TMD
            self.k: ureg.newton / ureg.meter,  # Stiffness of the spring for trolley
            self.k_TMD: ureg.newton / ureg.meter,  # Stiffness of the TMD spring
            self.F: ureg.newton,  # Force applied to the system
            self.Omega: ureg.radian / ureg.second,  # Frequency of excitation
            self.x: ureg.meter,  # Displacement of the trolley
            self.x.diff(self.ivar): ureg.meter / ureg.second,  # Velocity of the trolley
            self.x.diff(self.ivar, 2): ureg.meter
            / ureg.second**2,  # Acceleration of the trolley
            self.x_e: ureg.meter,  # Displacement of the TMD relative to trolley
            self.x_e.diff(self.ivar): ureg.meter
            / ureg.second,  # Velocity of the TMD relative to trolley
            self.x_e.diff(self.ivar, 2): ureg.meter
            / ureg.second**2,  # Acceleration of the TMD relative to trolley
            self.ivar: ureg.second,  # Time
        }
        return units_dict


# AF - system Amadiego, zmienne uspójnione zgodnie z obrazkiem (m_TMD-> m_e itd.)
class TrolleyWithTMD1(ComposedSystem):
    scheme_name = "MDoF_new.png"
    real_name = "mdof_tmd_real.png"

    m = Symbol("m", positive=True)
    m_e = Symbol("m_e", positive=True)
    k = Symbol("k", positive=True)
    k_e = Symbol("k_e", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    x_e = dynamicsymbols("x_e")
    x = dynamicsymbols("x")

    def __init__(
        self,
        m=None,
        m_e=None,
        k=None,
        k_e=None,
        Omega=None,
        x_e=None,
        x=None,
        F=None,
        ivar=Symbol("t"),
        **kwargs,
    ):
        if m_e is not None:
            self.m_e = m_e
        if m is not None:
            self.m = m
        if x_e is not None:
            self.x_e = x_e
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if k_e is not None:
            self.k_e = k_e
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        self.ivar = ivar

        self.qs = [self.x, self.x_e]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        self._TMD = TunedMassDamperRelativeMotion(
            self.m_e, self.k_e, self.x_e, self.x, self.ivar
        )(label="Tuned Mass Damper")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.x_e]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_TMD"] = self._TMD
        components["_force"] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.x_e: r"Angle of a pendulum",
            self.m: r"Mass of trolley",
            self.m_e: r"Mass of TMD",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, F0, Omega0, k0 = symbols("m_0 F_0 Omega_0 k_0", positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(20, 30)],
            self.m_e: [S.One * m0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [no for no in range(20, 30)],
            self.m_e: [no for no in range(1, 10)],
            self.F: [no for no in range(50, 100)],
            self.Omega: [3.14 * no for no in range(1, 6)],
            self.k: [no for no in range(50, 100)],
        }
        return default_data_dict

    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,  # Mass of the trolley
            self.m_e: ureg.kilogram,  # Mass of TMD
            self.k: ureg.newton / ureg.meter,  # Stiffness of the spring for trolley
            self.k_e: ureg.newton / ureg.meter,  # Stiffness of the TMD spring
            self.F: ureg.newton,  # Force applied to the system
            self.Omega: ureg.radian / ureg.second,  # Frequency of excitation
            self.x: ureg.meter,  # Displacement of the trolley
            self.x.diff(self.ivar): ureg.meter / ureg.second,  # Velocity of the trolley
            self.x.diff(self.ivar, 2): ureg.meter
            / ureg.second**2,  # Acceleration of the trolley
            self.x_e: ureg.meter,  # Displacement of the TMD relative to trolley
            self.x_e.diff(self.ivar): ureg.meter
            / ureg.second,  # Velocity of the TMD relative to trolley
            self.x_e.diff(self.ivar, 2): ureg.meter
            / ureg.second**2,  # Acceleration of the TMD relative to trolley
            self.ivar: ureg.second,  # Time
        }
        return units_dict


class VariableMassTrolleyWithPendulumFunction(ComposedSystem):

    scheme_name = "DampedTrolleyWithPendulumGravityOtherCoordinates.png"
    real_name = "taipei101.png"

    ivar = Symbol("t")
    l = Symbol("l", positive=True)
    rho = Symbol("rho", positive=True)
    m_f = Symbol("epsilon_f", positive=True)
    m_t = Function("m_t")(ivar)
    m_t0 = Symbol("m_t0", positive=True)
    m_p = Function("m_p")(ivar)
    m_p0 = Symbol("m_p0", positive=True)
    lam = Symbol("lambda", positive=True)
    t0 = Symbol("tau_0", real=True)
    k = Symbol("k", positive=True)
    omega0 = Symbol("omega_0", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("Phi")
    x = dynamicsymbols("X")
    eps = Symbol("epsilon", positive=True)
    m_s = Symbol("m_s", positive=True)
    chi = Symbol("chi", positive=True)
    psi = Symbol("psi", positive=True)
    delta = Symbol("delta", positive=True)
    Omega0 = Symbol("Omega_0", positive=True)
    tau = Symbol("tau", positive=True)

    def __init__(
        self,
        l=None,
        rho=None,
        m_f=None,
        m_t=None,
        m_t0=None,
        m_p=None,
        m_p0=None,
        lam=None,
        t0=None,
        k=None,
        g=None,
        Omega=None,
        omega0=None,
        phi=None,
        x=None,
        F=None,
        eps=None,
        m_s=None,
        chi=None,
        psi=None,
        delta=None,
        Omega0=None,
        tau=None,
        ivar=ivar,
        **kwargs,
    ):

        if l is not None:
            self.l = l
        if rho is not None:
            self.rho = rho
        if m_f is not None:
            self.m_f = m_f
        if m_t is not None:
            self.m_t = m_t
        if m_t0 is not None:
            self.m_t0 = m_t0
        if m_p is not None:
            self.m_p = m_p
        if m_p0 is not None:
            self.m_p0 = m_p0
        if lam is not None:
            self.lam = lam
        if t0 is not None:
            self.t0 = t0
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if omega0 is not None:
            self.omega0 = omega0
        if F is not None:
            self.F = F
        if eps is not None:
            self.eps = eps
        if m_s is not None:
            self.m_s = m_s
        if chi is not None:
            self.chi = chi
        if psi is not None:
            self.psi = psi
        if delta is not None:
            self.delta = delta
        if Omega0 is not None:
            self.Omega0 = Omega0
        if tau is not None:
            self.tau = tau
        self.ivar = ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        self._trolley_damper = Damper(
            c=self.k * self.lam, pos1=self.x, pos2=0, qs=[self.x, self.phi]
        )(label="Trolley damper")
        self._pendulum = PendulumKinematicExct(
            self.l, self.m_p, self.g, self.phi, self.x, self.ivar
        )(label="Pendulum")
        self._pendulum_damper = Damper(
            c=self.m_p * self.g * self.l * self.lam,
            pos1=self.phi,
            pos2=0,
            qs=[self.x, self.phi],
        )(label="Pendulum damper")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.phi]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_trolley_damper"] = self._trolley_damper
        components["_pendulum"] = self._pendulum
        components["_pendulum_damper"] = self._pendulum_damper
        components["_force"] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.k: r"Stiffness of a beam showed as a spring stiffness in trolley member",
            self.x: r"Kinematic lateral excitation",
            self.phi: r"Angle of a pendulum",
            self.m_p: r"Mass of pendulum",
            self.m_t: r"Mass of trolley",
            self.g: r"Gravity constant",
            self.F: r"Force",
            self.Omega: r"Excitation frequency",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols("m_0 l_0 F_0 Omega_0 k_0", positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.g: [S.One * self.g],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_t0: r"initial mass of the trolley",
            self.m_t: r"mass of the trolley",
            self.m_p0: r"initial mass of the pendulum",
            self.m_p: r"mass of the pendulum",
            self.m_f: r"mass of the transferred liquid to the mass of the system ratio",
            self.k: r"""trolley's equivalent stiffness""",
            self.lam: r"""Rayleigh's proportional to stiffness damping coefficient""",
            self.l: r"length of the pendulum",
            self.Omega: r"excitation frequency",
            self.omega0: r"natural frequency of the trolley",
            self.F: r"forcing amplitude",
            self.g: r"gravitational constant",
            self.rho: r"flow rate coefficient",
            self.t0: r"mass transfer activation time",
            self.eps: r"mass of the pendulum to the mass of the system ratio",
            self.m_s: r"total mass of the system",
            self.chi: r"members natural frequencies ratio",
            self.psi: r"excitation frequency to the natural frequency of the trolley ratio",
            self.delta: r"""forcing amplitude to structure's equivalent stiffness ratio""",
            self.Omega0: r"natural frequency of the pendulum",
            self.tau: "dimensionless time",
            self.x: """trolley's horizontal displacement""",
            self.x.diff(self.ivar): """trolley's horizontal velocity""",
            self.x.diff(self.ivar, 2): """trolley's horizontal acceleration""",
            self.phi.diff(self.ivar): """pendulum's angular velocity""",
            self.phi.diff(self.ivar, 2): """pendulum's angular acceleration""",
            self.phi: """pendulum's angular displacement""",
            self.ivar: "time",
            self.m_p.diff(self.ivar): "differentiated mass of the pendulum",
        }
        return self.sym_desc_dict

    def symbols_description_polish(self):
        self.sym_desc_dict = {
            self.m_t0: r"początkowa masa wózka",
            self.m_t: r"masa wózka",
            self.m_p0: r"początkowa masa wahadła",
            self.m_p: r"masa wahadła",
            self.m_f: r"stosunek masy transferowanej cieczy do masy systemu",
            self.k: r"""współczynnik równoważny sztywność konstrukcji""",
            self.lam: r"""współczynnik proporcjonalny do macierzy sztywności w tłumieniu Rayleigh'a""",
            self.l: r"długość wahadła",
            self.Omega: r"częstość wymuszenia",
            self.omega0: r"częstość własna wózka",
            self.F: r"amplituda siły wymuszającej",
            self.m_f: r"masa transferowanej cieczy",
            self.g: r"stała grawitacyjna",
            self.rho: r"współczynnik natężenia przepływu",
            self.t0: r"bezwymiarowy czas aktywacji transferu masy",
            self.eps: r"stosunek masy wahadła do masy systemu",
            self.m_s: r"całkowita masa systemu",
            self.chi: r"stosunek częstotliwości własnych elementów",
            self.psi: r"stosunek częstośći wymuszenia do częstści własnej wózka",
            self.delta: r"""stosunek amplitudy siły wymuszenia do współczynnika równoważnego sztywności konstrukcji""",
            self.Omega0: r"częstość własna wahadła",
            self.tau: "bezwymiarowy czas",
            self.x: """poziome przemieszczenie wózka""",
            self.x.diff(self.ivar): """prędkość pozioma wózka""",
            self.x.diff(self.ivar, 2): """przyspieszenie poziome wózka""",
            self.phi.diff(self.ivar): """prędkość kątowa wahadła""",
            self.phi.diff(self.ivar, 2): """przyspieszenie kątowe wahadła""",
            self.phi: """przemieszczenie kątowe wahadła""",
            self.ivar: "czas",
            self.m_p.diff(self.ivar): "pochodna masy wahadła",
        }
        return self.sym_desc_dict

    def unit_dict(self):
        from sympy.physics import units

        ureg = units  # UnitRegistry()

        unit_dict = {
            self.m_t0: ureg.kilogram,
            self.m_t: ureg.kilogram,
            self.m_p0: ureg.kilogram,
            self.m_p: ureg.kilogram,
            self.k: ureg.newton / ureg.meter,
            self.l: ureg.meter,
            self.Omega: ureg.rad,
            self.omega0: ureg.rad,
            self.F: ureg.newton,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.rho: (1 / ureg.second),
            self.m_s: ureg.kilogram,
            self.delta: ureg.meter,
            self.Omega0: ureg.rad,
            self.x: ureg.meter,
            self.x.diff(self.ivar): ureg.meter / ureg.second,
            self.x.diff(self.ivar, 2): ureg.meter / ureg.second / ureg.second,
            self.phi.diff(self.ivar): ureg.rad / ureg.second,
            self.phi.diff(self.ivar, 2): ureg.rad / ureg.second / ureg.second,
            self.phi: ureg.rad,
            self.ivar: ureg.second,
            self.m_p.diff(self.ivar): ureg.kilogram,
            #             self.m_f: ureg.dimensionless,
            #             self.lam: ureg.dimensionless,
            #             self.t0: ureg.dimensionless,
            #             self.eps: ureg.dimensionless,
            #             self.chi: ureg.dimensionless,
            #             self.psi: ureg.dimensionless,
            self.tau: "-",
        }

        return unit_dict

    def dimensionless(self):

        q = self.q
        m_p = self.m_p
        m_t = self.m_t
        ivar = self.ivar

        inertia_matrix = self.dimensionless_inertia_matrix()
        damping_matrix = self.dimensionless_damping_matrix()
        stiffness_matrix = self.dimensionless_stiffness_matrix()
        external_forces_matrix = self.dimensionless_external_forces_matrix()

        eoms = (
            inertia_matrix * q.diff(ivar, 2)
            + damping_matrix * q.diff(ivar)
            + stiffness_matrix * q
            - external_forces_matrix
        )
        eoms_left = (
            inertia_matrix * q.diff(ivar, 2)
            + damping_matrix * q.diff(ivar)
            + stiffness_matrix * q
        )

        return eoms

    def const_mass(self):

        subs_dict = {
            self.m_p: Symbol("m_p", positive=True),
            self.m_t: Symbol("m_t", positive=True),
        }

        return self.subs(subs_dict)

    def dimensionless_inertia_matrix(self):

        inertia_mat = self.inertia_matrix()
        first_inertia_row = inertia_mat.row(0) / inertia_mat[0]
        second_inertia_row = inertia_mat.row(1) / inertia_mat[3]
        dimensionless_inertia_mat = Matrix([[first_inertia_row], [second_inertia_row]])

        l = self.l
        m_p = self.m_p
        m_t = self.m_t
        inertia_matrix = Matrix([[1, l * m_p / (m_p + m_t)], [1 / l, 1]])

        return inertia_matrix

    def dimensionless_damping_matrix(self):

        inertia_mat = self.inertia_matrix()
        first_damping_row = self.damping_matrix().row(0) / inertia_mat[0]
        second_damping_row = self.damping_matrix().row(1) / inertia_mat[3]
        dimensionless_damping_mat = Matrix([[first_damping_row], [second_damping_row]])

        m_p = self.m_p
        m_t = self.m_t
        omega0 = self.omega0
        lam = self.lam
        l = self.l
        chi = self.chi

        ivar = self.ivar

        damping_matrix = Matrix(
            [
                [lam * omega0, l / omega0 * m_p.diff(ivar) / (m_p + m_t)],
                [
                    m_p.diff(ivar) / m_p / l / omega0,
                    (lam * omega0 * chi**2 + m_p.diff(ivar) / omega0 / m_p),
                ],
            ]
        )

        return damping_matrix

    def dimensionless_stiffness_matrix(self):

        inertia_mat = self.inertia_matrix()
        first_stiffness_row = self.stiffness_matrix().row(0) / inertia_mat[0]
        second_stiffness_row = self.stiffness_matrix().row(1) / inertia_mat[3]
        dimensionless_stiffness_mat = Matrix(
            [[first_stiffness_row], [second_stiffness_row]]
        )

        chi = self.chi

        stiffness_matrix = Matrix([[1, 0], [0, chi**2]])

        return stiffness_matrix

    def dimensionless_external_forces_matrix(self):

        inertia_mat = self.inertia_matrix()
        dimensionless_external_forces_mat = self.external_forces() / inertia_mat[0]
        dimensionless_external_forces_mat.doit()

        delta = self.delta
        psi = self.psi
        ivar = self.ivar

        external_forces_matrix = Matrix([[delta * sin(psi * ivar)], [0]])
        return external_forces_matrix

    def mass_transfer_function(self):

        subs_dict = {
            self.m_t: self.m_t0
            + self.m_f
            * (
                1 / 2
                - atan(
                    self.rho * self.ivar / self.omega0
                    - self.rho * self.t0 / self.omega0
                )
                / pi
            ),
            self.m_p: self.m_p0
            + self.m_f
            * (
                1 / 2
                + atan(
                    self.rho * self.ivar / self.omega0
                    - self.rho * self.t0 / self.omega0
                )
                / pi
            ),
        }
        return subs_dict

    def dimensionless_with_transfer(self):
        return self.dimensionless().subs(self.mass_transfer_function()).doit()

    def dimensionless_with_const_mass(self):

        m_p_const = Symbol("m_p", positive=True)
        m_t_const = Symbol("m_t", positive=True)

        masses_dict = {self.m_p: m_p_const, self.m_t: m_t_const}

        return (
            self.dimensionless()
            .subs(masses_dict)
            .doit()
            .subs({m_t_const: (-m_p_const + m_p_const / self.eps)})
        )  # it's a bit tricky, start with checking it if something doesn't work

    def system_mass_ratio(self):
        return Eq(Function("epsilon")(self.ivar), self.m_p / self.m_s)

    def total_system_mass(self):
        return Eq(self.m_s, self.m_p + self.m_t)

    def initial_system_mass_ratio(self):
        return Eq(Symbol("epsilon_0", positive=True), self.m_p0 / self.m_s)

    def excitation_freq_ratio(self):
        return Eq(self.psi, self.Omega / self.omega0)

    def frequency_ratio(self):
        return Eq(self.chi, self.Omega0 / self.omega0)

    def force_amplitude_ratio(self):
        return Eq(self.delta, self.F / self.k)

    def dimensionless_time(self):
        return Eq(self.tau, self.ivar * self.omega0)

    def get_numerical_data(self):

        params_dict = {
            self.m_s: 150000000,
            self.chi: 1,
            self.psi: 1,
            self.lam: 0.03,
            self.l: 25,
            self.omega0: sqrt(9.81 / 25),
            self.m_p0: 150000,
            self.m_f: 3 / 1000,
            self.rho: 10,
            pi: 3.1415,
            self.t0: 10 * 3.1415,
            self.delta: 300000 / (9.81 * (150000000) / 25),
        }

        return params_dict

    def get_analytical_data_before_transfer(self, psi=False, chi=None, eps=None):

        m_p_const = Symbol("m_p", positive=True)
        m_t_const = Symbol("m_t", positive=True)

        if psi == False:
            var = 1
        elif psi == True:
            var = self.psi
        elif isinstance(psi, float):
            var = psi

        if chi == None:
            var_chi = 1
        else:
            var_chi = chi

        if eps == None:
            var_eps = 1 / 1000
        else:
            var_eps = eps

        params_dict = {
            self.eps: var_eps,
            self.chi: var_chi,
            self.psi: var,
            self.lam: 0.03,
            self.l: 25,
            self.omega0: sqrt(9.81 / 25),
            m_p_const: 150000,
            m_t_const: 150000000 - 1 * 150000,
            pi: 3.1415,
            self.delta: 300000 / (9.81 * (150000000) / 25),
        }

        return params_dict

    def get_analytical_data_after_transfer(self, psi=False, chi=None, eps=None):

        m_p_const = Symbol("m_p", positive=True)
        m_t_const = Symbol("m_t", positive=True)

        if psi == False:
            var = 1
        elif psi == True:
            var = self.psi
        elif isinstance(psi, float):
            var = psi

        if eps == None:
            var_eps = 4 / 1000
        else:
            var_eps = eps

        if chi == None:
            var_chi = 1
        else:
            var_chi = chi

        params_dict = {
            self.eps: var_eps,
            self.chi: var_chi,
            self.psi: var,
            self.lam: 0.03,
            m_p_const: 4 * 150000,
            m_t_const: 150000000 - 4 * 150000,
            self.l: 25,
            self.omega0: sqrt(9.81 / 25),
            pi: 3.1415,
            self.delta: 300000 / (9.81 * (150000000) / 25),
        }

        return params_dict

    def get_analytical_data_without_pendulum(self, psi=False, omega0=None):

        m_p_const = Symbol("m_p", positive=True)
        m_t_const = Symbol("m_t", positive=True)

        if psi == False:
            var = 1
        else:
            var = self.psi

        if omega0 is not None:
            omega0_var = omega0
        else:
            omega0_var = sqrt(9.81 / 25)

        params_dict = {
            self.eps: 0,
            self.chi: 0,
            self.psi: var,
            self.lam: 0.03,
            m_p_const: 0,
            m_t_const: 150000000,
            self.l: 0,
            self.omega0: omega0_var,
            pi: 3.1415,
            self.delta: 300000 / (9.81 * (150000000) / 25),
        }

        return params_dict

    def get_numerical_data_for_solving_options(self):

        params_dict = {
            self.chi: 1,
            self.psi: 1,
            self.lam: 0.03,
            self.l: 25,
            self.omega0: sqrt(9.81 / 25),
            self.m_p0: 150000,
            self.m_t0: 150000000 - 4 * 150000,
            self.m_f: 3 * 150000,
            self.rho: 10,
            self.t0: 90,
            pi: 3.1415,
            self.delta: 300000 / (9.81 * (150000000) / 25),
        }

        return params_dict

    def _FRF_chart(self, chi=None):
        from sympy import lambdify
        from sympy.physics import units

        from ...solvers.linear import FirstOrderLinearODESystemWithHarmonics, ODESystem

        ureg = units
        import pandas as pd

        tau = Symbol("tau")
        sym = Symbol("A")
        sym_psi = Symbol("psi")
        sym_A = f"${latex(sym)}$[m]"
        ode_eoms = self.dimensionless_with_const_mass()

        if chi is not None:
            var_chi = chi
        else:
            var_chi = 1

        odes_przed = ode_eoms.subs(
            self.get_analytical_data_before_transfer(psi=True, chi=self.chi)
        ).n()
        odes_po = ode_eoms.subs(
            self.get_analytical_data_after_transfer(psi=True, chi=self.chi)
        ).n()
        odes_bez = Matrix(
            [ode_eoms[0].subs(self.get_analytical_data_without_pendulum(psi=True)).n()]
        )
        odes_bez_mniejsza_masa = Matrix(
            [
                (
                    ode_eoms[0].subs(
                        {
                            self.delta: self.delta * 1.00300903,
                            self.omega0: self.omega0**2 / 0.6264,
                        }
                    )
                    - 0.003 * self.x
                )
                .subs(
                    self.get_analytical_data_without_pendulum(psi=True, omega0=0.6273)
                )
                .n()
            ]
        )

        fode_system_przed = ODESystem(
            odes=odes_przed, dvars=self.q, ode_order=2
        )._as_fode()
        fode_system_po = ODESystem(odes=odes_po, dvars=self.q, ode_order=2)._as_fode()
        fode_system_bez = ODESystem(
            odes=odes_bez, dvars=Matrix([self.q[0]]), ode_order=2
        )._as_fode()
        fode_system_bez_mniejsza_masa = ODESystem(
            odes=odes_bez_mniejsza_masa, dvars=Matrix([self.q[0]]), ode_order=2
        )._as_fode()

        ana_sol_przed = fode_system_przed.steady_solution.subs(self.ivar, tau)[0].subs(
            self.chi, var_chi
        )

        ana_sol_po = fode_system_po.steady_solution.subs(self.ivar, tau)[0].subs(
            self.chi, var_chi
        )

        ana_sol_bez = fode_system_bez.steady_solution.subs(self.ivar, tau)[0]

        ana_sol_bez_mniejsza_masa = fode_system_bez_mniejsza_masa.steady_solution.subs(
            self.ivar, tau
        )[0]

        Amp_przed = sqrt(
            (ana_sol_przed.coeff(sin(self.psi * tau))) ** 2
            + (ana_sol_przed.coeff(cos(self.psi * tau))) ** 2
        )

        Amp_po = sqrt(
            (ana_sol_po.coeff(sin(self.psi * tau))) ** 2
            + (ana_sol_po.coeff(cos(self.psi * tau))) ** 2
        )

        Amp_bez = sqrt(
            (ana_sol_bez.coeff(sin(self.psi * tau))) ** 2
            + (ana_sol_bez.coeff(cos(self.psi * tau))) ** 2
        )

        Amp_bez_mniejsza_masa = sqrt(
            (ana_sol_bez_mniejsza_masa.coeff(sin(self.psi * tau))) ** 2
            + (ana_sol_bez.coeff(cos(self.psi * tau))) ** 2
        )

        psi_span = np.linspace(0.1, 1.9, 1801)

        data4plot_przed = lambdify(self.psi, Amp_przed, "numpy")(psi_span)
        data4plot_po = lambdify(self.psi, Amp_po, "numpy")(psi_span)
        data4plot_bez = lambdify(self.psi, Amp_bez, "numpy")(psi_span)
        data4plot_bez_mniejsza_masa = lambdify(
            self.psi, Amp_bez_mniejsza_masa, "numpy"
        )(psi_span)

        df = pd.DataFrame(
            data={
                f"bez TMD BT": data4plot_bez,
                f"bez TMD AT": data4plot_bez_mniejsza_masa,
                f"BT, ${latex(Eq(self.chi,var_chi,evaluate=False))}$": data4plot_przed,
                f"AT, ${latex(Eq(self.chi,var_chi,evaluate=False))}$": data4plot_po,
            },
            index=psi_span,
        )

        #         df = df.set_axis(pd.Index([ (val,sym_A)  for val  in df.columns]),axis=1)
        df = df.set_axis(pd.Index([(val, sym_A) for val in df.columns]), axis=1)
        df.index.name = self.psi

        return df

    def _simulate_FRFs(self, params_list=None):

        import pandas as pd

        sym_psi = Symbol("psi")
        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.9, 1, 1.1, 1.2, 1.3, 1.5]

        df = pd.DataFrame()
        psi_span = np.linspace(0.1, 1.9, 1801)
        df.index = psi_span

        for param in param_list:
            df_param = self._FRF_chart(
                chi=param
            )  # .reset_index(drop=True, inplace=True)

            df = pd.concat([df, df_param], axis=1)

        df = df.loc[:, ~df.T.duplicated(keep="first")]
        df.index.name = f"${latex(sym_psi)}$[$-$]"

        return df

    def excitation_frequency_sim(self, t_span=None, params_list=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.psi)
        params_dict[self.t0] = 150

        if t_span is None:
            t_span = np.linspace(0, 100, 1001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.8, 0.9, 1, 1.1, 1.2]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.psi,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def force_to_stiffness_ratio_sim(self, t_span=None, params_list=None, t0=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.delta)

        if t0 is None:
            t0 = 10 * 3.1415
        params_dict[self.t0] = t0

        if t_span is None:
            t_span = np.linspace(0, 100, 1001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.0011, 0.0031, 0.0051, 0.0071, 0.0091]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.delta,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def rayleigh_damping_coefficient_sim(self, t_span=None, params_list=None, t0=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.lam)

        if t0 is None:
            t0 = 31.415
        params_dict[self.t0] = t0

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.01, 0.02, 0.03, 0.04, 0.05]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.lam,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def flow_rate_coefficient_sim(self, t_span=None, params_list=None, t0=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.rho)

        if t0 is None:
            t0 = 10 * 3.1415
        params_dict[self.t0] = t0

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [1, 5, 10, 15, 20]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.rho,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def activation_time_sim(self, t_span=None, params_list=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.t0)

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0, 28.2735, 31.415, 91.1035, 330]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.t0,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def initial_mass_ratio_sim(self, t_span=None, params_list=None, t0=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}

        if t0 is None:
            t0 = 31.415
        params_dict[self.t0] = t0

        eps0 = Symbol("epsilon_0", positive=True)

        if t_span is None:
            t_span = np.linspace(0, 200, 2001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.0001, 0.0005, 0.001, 0.003, 0.005]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_p_dict = {self.m_p0: eps0 * self.m_s}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer()
            .subs(masses_dict)
            .subs(m_f_dict)
            .subs(m_p_dict)
            .doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            eps0,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def _activation_time_chart(
        self, transfer_func=None, data_dict=None, upper_func=True
    ):

        import pandas as pd
        from sympy import lambdify

        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        if upper_func == True:
            var = self.m_t
        else:
            var = self.m_p

        atan_expr = transfer_func.subs(data_dict)

        lambda_atan_expr = lambdify((self.tau, self.t0, self.rho), atan_expr.n())

        t_span = np.linspace(0, 100, 1001)

        t0_list = [-10, 25, 75, 110]

        fluid_flow_t0 = pd.DataFrame(data={})
        fluid_flow_t0.index.name = Symbol("tau")

        for t0num, t0val in enumerate(t0_list):
            fluid_flow_t0[t0val] = pd.DataFrame(
                (lambda_atan_expr(t_span, t0val, 10)), index=t_span
            )

        fluid_flow_t0 = TimeDataFrame(
            fluid_flow_t0.set_axis(
                pd.Index(
                    [
                        (Eq(self.t0, val, evaluate=False), var)
                        for val in fluid_flow_t0.columns
                    ]
                ),
                axis=1,
            )
        )
        fluid_flow_t0.index.name = Symbol("tau")
        fluid_flow_t0_tikz_data = TimeDataFrame(
            fluid_flow_t0[list((fluid_flow_t0.columns))]
        )

        return fluid_flow_t0_tikz_data

    def _flow_rate_chart(self, transfer_func=None, data_dict=None, upper_func=True):

        import pandas as pd
        from sympy import lambdify

        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        if upper_func == True:
            var = self.m_t
        else:
            var = self.m_p

        atan_expr = transfer_func.subs(data_dict)

        lambda_atan_expr = lambdify((self.tau, self.t0, self.rho), atan_expr.n())

        t_span = np.linspace(0, 100, 1001)

        flow_coeff_list = [1, 5, 10]

        fluid_flow_fc = pd.DataFrame(data={})
        fluid_flow_fc.index.name = Symbol("tau")

        for fcnum, fcval in enumerate(flow_coeff_list):
            fluid_flow_fc[fcval] = pd.DataFrame(
                (lambda_atan_expr(t_span, 50, fcval)), index=t_span
            )

        fluid_flow_fc = fluid_flow_fc.set_axis(
            pd.Index(
                [
                    (Eq(self.rho, val, evaluate=False), var)
                    for val in fluid_flow_fc.columns
                ]
            ),
            axis=1,
        )
        fluid_flow_fc.index.name = Symbol("tau")
        fluid_flow_tikz_data = TimeDataFrame(
            fluid_flow_fc[list((fluid_flow_fc.columns))]
        )

        return fluid_flow_tikz_data

    def _mass_transfer_influence(self, t_span=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.m_f)

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )

        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()

        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.m_f,
            [1 / 10000, 5 / 10000, 1 / 1000, 3 / 1000, 5 / 1000],
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def _skyscraper_stiffness_influence(self, t_span=None, params_list=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.chi)

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        if params_list is not None:
            param_list = params_list
        else:
            param_list = [0.8, 0.9, 1, 1.1, 1.2]

        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        m_f_dict = {self.m_f: self.m_f * self.m_s}
        eoms_num = (
            self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        )
        a_p = Symbol("a_p")
        dependencies = {a_p: self.x.diff().diff() + 25 * self.phi.diff().diff()}

        q = self.q

        dyn_sys = ODESystem(odes=eoms_num, dvars=q)  # .as_first_ode_linear_system()
        Y = self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(
            dyn_sys,
            self.chi,
            param_list,
            reference_data=params_dict,
            coordinates=Y,
            index=pd.Index(t_span, name=Symbol("tau")),
        )

        sym_wyn = num_df.perform_simulations(
            model_level_name=0, dependencies=dependencies
        )

        result = TimeDataFrame(sym_wyn)
        #         result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def _solving_options_analysis(self, t_span=None):
        from ...solvers.linear import FirstOrderLinearODESystemWithHarmonics, ODESystem
        from ...solvers.tools import TwoSystemsSimulation
        from ...utilities.adaptable import NumericalAnalysisDataFrame, TimeDataFrame, pd

        params_dict = {**self.get_numerical_data()}

        if t_span is None:
            t_span = np.linspace(0, 300, 3001)

        # ANALITYCZNIE
        eoms1 = (
            self.dimensionless_with_const_mass()
            .subs(self.get_analytical_data_before_transfer())
            .doit()
        )

        ode1 = ODESystem(odes=eoms1, dvars=self.q)

        eoms2 = (
            self.dimensionless_with_const_mass()
            .subs(self.get_analytical_data_after_transfer())
            .doit()
        )
        ode2 = ODESystem(odes=eoms2, dvars=self.q)

        odes = [ode1, ode2]

        symki = TwoSystemsSimulation(
            odes=(ode1, ode2),
            t_span=(np.linspace(0, 89.9, 900), np.linspace(90, 300, 2101)),
        )
        wyniki = symki.compute_solution(
            t_span=(np.linspace(0, 89.9, 900), np.linspace(90, 300, 2101))
        ).apply(np.real)

        analytical_sym_wyn = TimeDataFrame(data={("A-A ", "X"): wyniki[self.x]})
        analytical_sym_wyn.index.name = Symbol("tau")

        # NUMERYCZNIE
        masses_dict = {self.m_t0: self.m_s - self.m_p0 - self.m_f}
        eoms_num = (
            self.dimensionless_with_transfer()
            .subs(self.get_numerical_data_for_solving_options())
            .doit()
        )

        ode_2 = ODESystem(odes=eoms_num, dvars=self.q).numerized()
        num_sym_wyn = ode_2.compute_solution(
            t_span=t_span, ic_list=[0.0, 0.0, 0.0, 0.0]
        ).apply(np.real)

        numerical_sym_wyn = TimeDataFrame(data={("N-N-N", "X"): num_sym_wyn[self.x]})
        numerical_sym_wyn.index.name = Symbol("tau")

        # ANALITYCZNIE + NUMERYCZNIE
        fode_1 = ode1.as_first_ode_linear_system()
        fodes_1 = FirstOrderLinearODESystemWithHarmonics(fode_1, fode_1.dvars)
        ana_sol_1 = fodes_1.solution
        ics_list = [0.0, 0.0, 0.0, 0.0]
        ana_sol_sym_1 = ana_sol_1.with_ics(ics_list)
        t_span_1 = np.linspace(0, 84.9, 850)
        ana_sym_wyn_1 = (
            ana_sol_sym_1.numerized().compute_solution(t_span=t_span_1).apply(np.real)
        )
        ics_list_new = list(ana_sym_wyn_1.iloc[-1, 0 : len(ics_list)])

        ##Symulacja fazy 2.
        t_span_2 = np.linspace(85, 95, 101)
        num_sym_wyn_2 = ode_2.compute_solution(
            t_span=t_span_2, ic_list=ics_list_new
        ).apply(np.real)
        ics_list_the_newest = list(num_sym_wyn_2.iloc[-1, 0 : len(ics_list)])

        ##Symulacja fazy 3.
        fode_3 = ode2.as_first_ode_linear_system()
        fodes_3 = FirstOrderLinearODESystemWithHarmonics(fode_3, fode_3.dvars)
        ana_sol_3 = fodes_3.solution
        ana_sol_sym_3 = ana_sol_3.with_ics(ics_list_the_newest, ivar0=95.1)
        t_span_3 = np.linspace(95.1, 300, 2050)
        ana_sym_wyn_3 = (
            ana_sol_sym_3.numerized().compute_solution(t_span=t_span_3).apply(np.real)
        )

        results_ANA = pd.concat([ana_sym_wyn_1, num_sym_wyn_2, ana_sym_wyn_3])

        analytical_numerical_sym_wyn = TimeDataFrame(
            data={("A-N-A", "X"): results_ANA[self.x]}
        )
        analytical_numerical_sym_wyn.index.name = Symbol("tau")

        results = pd.concat(
            [analytical_numerical_sym_wyn, numerical_sym_wyn, analytical_sym_wyn]
        )
        return results


class InclinedSpringMassSystem(ComposedSystem):
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

    """

    scheme_name = "InclinedSpringMassSystem.png"

    m = Symbol("m", positive=True)
    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    k = Symbol("k", positive=True)
    ivar = Symbol("t")
    u = dynamicsymbols("u")
    alpha = Symbol("alpha", positive=True)
    x = dynamicsymbols("x")
    h0 = Symbol("h_0", positive=True)
    x_ramp = x
    x_mass = x + u * cos(alpha)
    y_mass = h0 - u * sin(alpha)
    p = 3 * cos(ivar)

    def __init__(
        self,
        m1=None,
        m2=None,
        k=None,
        u=None,
        x=None,
        h0=None,
        x_ramp=None,
        x_mass=None,
        y_mass=None,
        alpha=None,
        ivar=None,
        p=None,
        **kwargs,
    ):

        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if h0 is not None:
            self.h0 = h0
        if k is not None:
            self.k = k
        if alpha is not None:
            self.alpha = alpha
        if ivar is not None:
            self.ivar = ivar
        if x is not None:
            self.x = x
        if x_ramp is not None:
            self.x_ramp = x_ramp
        if x_mass is not None:
            self.x_mass = x_mass
        if y_mass is not None:
            self.y_mass = y_mass
        if p is not None:
            self.p = p

        self.qs = [self.x_ramp, self.u]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.mass_horizontal = MaterialPoint(self.m2, self.x_mass, qs=self.qs)
        self.mass_vertical = MaterialPoint(self.m2, self.y_mass, qs=self.qs)
        self.spring = Spring(self.k, self.u, qs=self.qs)
        self.incline = MaterialPoint(self.m1, self.x_ramp, qs=[self.qs])
        self.excitation = Force(self.p, self.x_ramp, qs=self.qs, ivar=self.ivar)

        components["spring"] = self.spring
        components["incline"] = self.incline
        components["excitation"] = self.excitation
        components["mass_horizontal"] = self.mass_horizontal
        components["mass_vertical"] = self.mass_vertical

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of system on the spring",
            self.k: r"Spring coefficient ",
        }

        return self.sym_desc_dict


class DDoFSliderPendulum(TrolleyWithPendulum):

    m_1 = Symbol("m_1", positive=True)
    m_2 = Symbol("m_2", positive=True)
    u_0 = Symbol("u_0", positive=True)
    omega = Symbol("omega", positive=True)

    m_p = m_2
    m_t = m_1
    F = TrolleyWithPendulum.k * u_0

    Omega = omega

    def symbols_description(self):
        self.sym_desc_dict = {
            self.g: r"Przyśpieszenie ziemskie",
            self.omega: r"Częstość wymuszenia",
            self.k: r"Sztywność",
            self.u_0: r"Wymuszenie początkowe",
            self.l: r"Długość belki",
            self.m_2: r"Masa obciążenia",
            self.m_1: r"Masa wózka",
            self.phi: r"Współrzędna uogólniona, kąt obrotu wahadła",
            self.x: r"Współrzędna uogólniona, przemieszczenie wózka",
            self.a: r"przyśpieszenie",
            self.ivar: r"Czas",
        }
        return self.sym_desc_dict

    def unit_dict(self):
        dyn_sys = self

        unit_dict = {
            dyn_sys.phi: ureg.radian,
            dyn_sys.x: ureg.meter,
            dyn_sys.g: (ureg.meter / ureg.second**2),
            dyn_sys.u_0: ureg.meter,
            dyn_sys.phi: ureg.radian,
            dyn_sys.k: ureg.newton / ureg.meter,
            dyn_sys.omega: (ureg.rad / ureg.second**2),
            dyn_sys.m_1: ureg.kilogram,
            dyn_sys.m_2: ureg.kilogram,
            dyn_sys.l_1: ureg.meter,
            dyn_sys.l_2: ureg.meter,
            dyn_sys.l: ureg.meter,
            Symbol("a"): (ureg.meter / ureg.second**2),
            t: ureg.second,
            f: ureg.hertz,
        }
        return unit_dict


class VariableMassTrolleyWithPendulumFunctionWithoutRayleigh(
    VariableMassTrolleyWithPendulumFunction
):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "taipei101.png"

    ivar = Symbol("t")
    l = Symbol("l", positive=True)
    rho = Symbol("rho", positive=True)
    m_f = Symbol("epsilon_f", positive=True)
    m_t = Function("m_t")(ivar)
    m_t0 = Symbol("m_t0", positive=True)
    m_p = Function("m_p")(ivar)
    m_p0 = Symbol("m_p0", positive=True)
    c = Symbol("c", positive=True)
    t0 = Symbol("tau_0", real=True)
    k = Symbol("k", positive=True)
    omega0 = Symbol("omega_0", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("Phi")
    x = dynamicsymbols("X")
    eps = Symbol("epsilon", positive=True)
    m_s = Symbol("m_s", positive=True)
    chi = Symbol("chi", positive=True)
    psi = Symbol("psi", positive=True)
    delta = Symbol("delta", positive=True)
    Omega0 = Symbol("Omega_0", positive=True)
    tau = Symbol("tau", positive=True)

    def __init__(
        self,
        l=None,
        rho=None,
        m_f=None,
        m_t=None,
        m_t0=None,
        m_p=None,
        m_p0=None,
        c=None,
        t0=None,
        k=None,
        g=None,
        Omega=None,
        omega0=None,
        phi=None,
        x=None,
        F=None,
        eps=None,
        m_s=None,
        chi=None,
        psi=None,
        delta=None,
        Omega0=None,
        tau=None,
        ivar=ivar,
        **kwargs,
    ):

        if l is not None:
            self.l = l
        if rho is not None:
            self.rho = rho
        if m_f is not None:
            self.m_f = m_f
        if m_t is not None:
            self.m_t = m_t
        if m_t0 is not None:
            self.m_t0 = m_t0
        if m_p is not None:
            self.m_p = m_p
        if m_p0 is not None:
            self.m_p0 = m_p0
        if c is not None:
            self.c = c
        if t0 is not None:
            self.t0 = t0
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if k is not None:
            self.k = k
        if Omega is not None:
            self.Omega = Omega
        if omega0 is not None:
            self.omega0 = omega0
        if F is not None:
            self.F = F
        if eps is not None:
            self.eps = eps
        if m_s is not None:
            self.m_s = m_s
        if chi is not None:
            self.chi = chi
        if psi is not None:
            self.psi = psi
        if delta is not None:
            self.delta = delta
        if Omega0 is not None:
            self.Omega0 = Omega0
        if tau is not None:
            self.tau = tau
        self.ivar = ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(
            label="Trolley"
        )
        self._trolley_damper = Damper(
            c=self.c, pos1=self.x, pos2=0, qs=[self.x, self.phi]
        )(label="Trolley damper")
        self._pendulum = PendulumKinematicExct(
            self.l, self.m_p, self.g, self.phi, self.x, self.ivar
        )(label="Pendulum")
        self._pendulum_damper = Damper(
            c=self.c, pos1=self.phi, pos2=0, qs=[self.x, self.phi]
        )(label="Pendulum damper")
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.phi]
        )(label="Force")

        components["_trolley"] = self._trolley
        components["_trolley_damper"] = self._trolley_damper
        components["_pendulum"] = self._pendulum
        components["_pendulum_damper"] = self._pendulum_damper
        components["_force"] = self._force

        return components


class DiskOnTrolleyWithFriction(ComposedSystem):

    scheme_name = "ForcedWinch.png"
    real_name = "winch_mechanism_real.PNG"

    r = Symbol("r", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    ivar = Symbol("t")
    phi = dynamicsymbols("\\varphi")
    x_1 = dynamicsymbols("x_1")
    x_2 = dynamicsymbols("x_2")
    y_w = dynamicsymbols("y_w")
    y_d = dynamicsymbols("y_d")
    F = Symbol("F", positive=True)
    Fx = Symbol("Fx", positive=True)
    Fy = Symbol("Fy", positive=True)
    Omega = Symbol("Omega", positive=True)
    K = Symbol("K", positive=True)
    M = Symbol("M", positive=True)
    N = Symbol("N", positive=True)
    m_t = Symbol("m_t", positive=True)
    mu = Symbol("mu", positive=True)
    I = Symbol("I", positive=True)

    def __init__(
        self,
        r=None,
        l=None,
        m=None,
        phi=None,
        F=None,
        Omega=None,
        g=None,
        Fx=None,
        Fy=None,
        M=None,
        N=None,
        K=None,
        x_1=None,
        x_2=None,
        x=None,
        y_w=None,
        y_d=None,
        m_t=None,
        mu=None,
        I=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if r is not None:
            self.r = r
        if l is not None:
            self.l = l
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if N is not None:
            self.N = N
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if Fx is not None:
            self.Fx = Fx
        if Fy is not None:
            self.Fy = Fy
        if K is not None:
            self.N = N
        if M is not None:
            self.M = M
        if x_1 is not None:
            self.x_1 = x_1
        if x_2 is not None:
            self.x_2 = x_2
        if y_w is not None:
            self.y_w = y_w
        if y_d is not None:
            self.y_w = y_d
        if m_t is not None:
            self.m_t = m_t
        if mu is not None:
            self.mu = mu
        if I is not None:
            self.I = I
        self.ivar = ivar

        self.qs = [self.x_1, self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.x_1 = self.x_1
        self.y_w = self.y_w
        self.y_d = self.y_d
        self.phi = self.phi  # self.x_1/self.r
        self.x_2 = self.x_1 + self.r * self.phi
        self.material_point_1 = MaterialPoint(self.m, self.x_1, qs=self.qs)
        self.material_point_2 = MaterialPoint(self.M, self.x_2, qs=self.qs)
        self.rolling_disk = Disk(self.m * self.r**2 / 2, self.phi, qs=self.qs)
        self.force_1 = Force(self.F - self.K, self.x_1, self.qs)
        self.force_2 = Force(self.K, self.x_2, self.qs)
        self.force_3 = Force(self.N, self.y_d, self.qs)
        self.force_4 = Force(-self.N, self.y_w, self.qs)
        self.force_5 = Force(self.K * self.r, self.phi, self.qs)
        self.gravity_1 = GravitationalForce(self.m, self.g, self.y_w)
        self.gravity_2 = GravitationalForce(self.M, self.g, self.y_d)

        components["_material_point_1"] = self.material_point_1
        components["_material_point_2"] = self.material_point_2
        components["_gravity_1"] = self.gravity_1
        components["_gravity_2"] = self.gravity_2
        components["_force_1"] = self.force_1
        components["_force_2"] = self.force_2
        components["_force_3"] = self.force_3
        components["_force_4"] = self.force_4
        components["_force_5"] = self.force_5
        components["_rolling_disk"] = self.rolling_disk

        return components

    @property
    def _report_components(self):

        comp_list = [
            WalecNaDesceForceComponent,
            WalecNaDesceProblemComponent,
            WalecNaDesceNewtonComponent,
            WalecNaDesceWiezyComponent,
        ]

        return comp_list

    def walecmoment(self):
        m = self.m
        r = self.r
        return (1 / 2) * m * r**2

    def accelerationxp(self):
        Mp = self.m
        Mw = self.M
        F = self.F
        x_1 = self.x_1
        K = self.K
        ax1 = x_1.diff(t, 2)
        ax1 = (F - K) / Mp
        return ax1

    def accelerationxw(self):
        Mw = self.M
        x_2 = self.x_2
        K = self.K
        ax2 = x_2.diff(t, 2)
        ax2 = K / Mw
        return ax2

    def accelerationphi(self):
        r = self.r
        K = self.K
        phi = self.phi
        I = self.walecmoment()
        aphi = phi.diff(t, 2)
        aphi = -K * r / I
        return aphi

    def accelerationphi2(self):
        r = self.r
        K = self.K
        phi = self.phi
        I = self.I
        aphi2 = phi.diff(t, 2)
        aphi2 = -K * r / I
        return aphi2

    def position(self):
        x1 = self.x_1
        x2 = self.x_2
        r = self.r
        phi = self.phi
        x2 = r * phi + x1
        return x2

    def posdiff(self):
        x2 = self.position()
        x2diff = x2.diff(t, 2)
        return x2diff

    def inertiamoment0(self):
        mw = self.M
        r = self.r
        I0 = mw * r
        return I0

    def kcoeff(self):
        F = self.F
        mw = self.M
        mp = self.m
        r = self.r
        K = F * mw / (mp * r + mp + mw)
        return K

    def tboundary(self):
        g = self.g
        m = self.m_t
        mu = self.mu
        T = g * m * mu
        return T

    def boundaryforce(self):
        T = self.tboundary()
        K = self.kcoeff()
        F = solve(Eq(T, K), self.F)
        return F


class DiskOnTrolley(ComposedSystem):

    scheme_name = "ForcedWinch.png"
    real_name = "winch_mechanism_real.PNG"

    r = Symbol("r", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    ivar = Symbol("t")
    phi = dynamicsymbols("\\varphi")
    x_1 = dynamicsymbols("x_1")
    x_2 = dynamicsymbols("x_2")
    y = dynamicsymbols("y")
    F = Symbol("F", positive=True)
    Fx = Symbol("Fx", positive=True)
    Fy = Symbol("Fy", positive=True)
    Omega = Symbol("Omega", positive=True)
    K = Symbol("K", positive=True)
    M = Symbol("M", positive=True)
    N = Symbol("N", positive=True)
    m_t = Symbol("m_t", positive=True)
    mu = Symbol("mu", positive=True)
    I = Symbol("I", positive=True)

    def __init__(
        self,
        r=None,
        l=None,
        m=None,
        phi=None,
        g=None,
        F=None,
        Omega=None,
        Fx=None,
        Fy=None,
        M=None,
        N=None,
        K=None,
        x_1=None,
        x_2=None,
        x=None,
        y=None,
        m_t=None,
        mu=None,
        I=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if m is not None:
            self.m = m
        if r is not None:
            self.r = r
        if l is not None:
            self.l = l
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if N is not None:
            self.N = N
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        if Fx is not None:
            self.Fx = Fx
        if Fy is not None:
            self.Fy = Fy
        if K is not None:
            self.N = N
        if M is not None:
            self.M = M
        if x_1 is not None:
            self.x_1 = x_1
        if x_2 is not None:
            self.x_2 = x_2
        if y is not None:
            self.y = y
        if m_t is not None:
            self.m_t = m_t
        if mu is not None:
            self.mu = mu
        if I is not None:
            self.I = I
        self.ivar = ivar

        self.qs = [self.x_1, self.x_2, self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.x_1 = self.x_1
        self.y = self.y
        self.x_2 = self.x_2
        self.phi = self.phi
        self.material_point_1 = MaterialPoint(self.m, self.x_1, qs=self.qs)
        self.material_point_2 = MaterialPoint(self.M, self.x_2, qs=self.qs)
        self.rolling_disk = Disk(self.m * self.r**2 / 2, self.phi, qs=self.qs)
        self.force_1 = Force(self.F - self.K, self.x_1, self.qs)
        self.force_2 = Force(self.K, self.x_2, self.qs)
        self.force_3 = Force(self.N, self.y, self.qs)
        self.force_4 = Force(self.K * self.r, self.phi, self.qs)
        self.gravity_1 = GravitationalForce(self.m, self.g, self.y)
        self.gravity_2 = GravitationalForce(self.M, self.g, self.y)

        components["_material_point_1"] = self.material_point_1
        components["_material_point_2"] = self.material_point_2
        components["_gravity_1"] = self.gravity_1
        components["_gravity_2"] = self.gravity_2
        components["_force_1"] = self.force_1
        components["_force_2"] = self.force_2
        components["_force_3"] = self.force_3
        components["_force_4"] = self.force_4
        components["_rolling_disk"] = self.rolling_disk

        return components

    @property
    def _report_components(self):

        comp_list = [
            WalecNaDesceForceComponent,
            WalecNaDesceProblemComponent,
            WalecNaDesceNewtonComponent,
            WalecNaDesceWiezyComponent,
        ]

        return comp_list

    def walecmoment(self):
        m = self.m
        r = self.r
        return (1 / 2) * m * r**2

    def accelerationxp(self):
        Mp = self.m
        Mw = self.M
        F = self.F
        x_1 = self.x_1
        K = self.K
        ax1 = x_1.diff(t, 2)
        ax1 = (F - K) / Mp
        return ax1

    def accelerationxw(self):
        Mw = self.M
        x_2 = self.x_2
        K = self.K
        ax2 = x_2.diff(t, 2)
        ax2 = K / Mw
        return ax2

    def accelerationphi(self):
        r = self.r
        K = self.K
        phi = self.phi
        I = self.walecmoment()
        aphi = phi.diff(t, 2)
        aphi = -K * r / I
        return aphi

    def accelerationphi2(self):
        r = self.r
        K = self.K
        phi = self.phi
        I = self.I
        aphi2 = phi.diff(t, 2)
        aphi2 = -K * r / I
        return aphi2

    def position(self):
        x1 = self.x_1
        x2 = self.x_2
        r = self.r
        phi = self.phi
        x2 = r * phi + x1
        return x2

    def posdiff(self):
        x2 = self.position()
        x2diff = x2.diff(t, 2)
        return x2diff

    def inertiamoment0(self):
        mw = self.M
        r = self.r
        I0 = mw * r
        return I0

    def kcoeff(self):
        F = self.F
        mw = self.M
        mp = self.m
        r = self.r
        K = F * mw / (mp * r + mp + mw)
        return K

    def tboundary(self):
        g = self.g
        m = self.m_t
        mu = self.mu
        T = g * m * mu
        return T

    def boundaryforce(self):
        T = self.tboundary()
        K = self.kcoeff()
        F = solve(Eq(T, K), self.F)
        return F


class TrolleyWithElasticPendulum(TrolleyWithPendulum):


    scheme_name = 'trolley_pendulum_tmd.png'
    real_name = 'taipei101.png'

    l = Symbol('l', positive=True)
    m_t = Symbol('m_t', positive=True)
    m_p = Symbol('m_p', positive=True)
    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')
    u_e = dynamicsymbols('u')
    y = dynamicsymbols('y')
    v = Symbol('v', positive=True)
    k_l = Symbol('k_l', positive=True)
    u_0 = Symbol('u_0', positive=True)

    def __init__(self,
                 l=None,
                 m_t=None,
                 m_p=None,
                 k=None,
                 g=None,
                 Omega=None,
                 phi=None,
                 x=None,
                 F=None,
                 u_e=None,
                 v=None,
                 y=None,
                 k_l=None,
                 u_0=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if l is not None: self.l = l
        if m_t is not None: self.m_t = m_t
        if m_p is not None: self.m_p = m_p
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if u_e is not None: self.u_e = u_e
        if v is not None: self.v = v
        if y is not None: self.y = y
        if k_l is not None: self.k_l = k_l
        if u_0 is not None: self.u_0 = self.m_p * self.g / self.k_l
      #  if u_0 is not None: u_0 = self.m_p * self.g / self.k_l
        self.ivar = ivar
        self.qs = [self.x,self.phi,self.u_e]
        self.u_0=self.g*self.m_p/self.k_l
        self._init_from_components(**kwargs)


    @cached_property
    def components(self):
        components = {}

        #self.u_0 = self.m_p * self.g / self.k_l
        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(label='Trolley')
        #self._pendulum = PendulumKinematicExct(self.l+self.u, self.m_p, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._force=Force(self.F*sin(self.Omega*self.ivar), pos1=self.qs[0], qs=self.qs)(label='Force')



        self._material_point_1=MaterialPoint(self.m_p,self.x+(self.l+self.u_e+self.u_0)*sin(self.phi),qs=self.qs)
        self._material_point_2=MaterialPoint(self.m_p,(self.l+self.u_e+self.u_0)*cos(self.phi),qs=self.qs)
        self._gravity = GravitationalForce(self.m_p, self.g, -(self.l+self.u_e+self.u_0)*cos(self.phi), qs=self.qs)



        self._spring_2=Spring(self.k_l,self.u_e,qs=self.qs)


        components['_trolley'] = self._trolley
        #components['_pendulum'] = self._pendulum
        components['_force'] = self._force
        components['_material_point_1'] = self._material_point_1
        components['_material_point_2'] = self._material_point_2
        #components['_spring_1'] = self._spring_1
        components['_spring_2'] = self._spring_2
        components['_gravity'] = self._gravity

        return components
    def get_default_data(self):
        """
        Zwraca słownik domyślnych wartości numerycznych 
        dla parametrów systemu.
        """

        
        # Mapowanie symboli (atrybutów klasy) na wartości numeryczne.
        # Kluczami są obiekty Symbol (np. self.l), a nie stringi ("l").
        data_dict = {
            self.l: 1.0,
            self.m_t: 2.0,
            self.m_p: 0.5,
            self.k: 10.0,
            self.g: 9.81,
            self.Omega: 1.5,
            self.F: 0.5,
            self.k_l: 20.0,
        }
        
        # Zauważ, że self.u_0 jest obliczane automatycznie w metodzie __init__
        # na podstawie m_p, g oraz k_l, więc nie ma potrzeby definiować go 
        # jako oddzielnej danej wejściowej.
        
        return data_dict

class DampedTrolleyWithPendulumVariableInertia(TrolleyWithElasticPendulum):

    scheme_name = "trolley_pendulum_tmd.png"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol('l', positive=True)
    k_l = Symbol('k_l', positive=True)
    m_t = Symbol('m_t', positive=True)#m_trolley
    m_p = Symbol('m_p', positive=True)#m_pendulum
    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')
    c = Symbol('c', positive=True)
    u_e = dynamicsymbols('u')
    v = Symbol('v', positive=True)
    c_phi=Symbol('\\hat c_{\\varphi}',positive=True)
    c_l=Symbol('c_l',positive=True)

    def __init__(self,
                 l=None,
                 k_l=None,
                 m_t=None,
                 m_p=None,
                 k=None,
                 g=None,
                 Omega=None,
                 phi=None,
                 x=None,
                 F=None,
                 c=None,
                 u_e=None,
                 v=None,
                 c_phi=None,
                 c_l=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if l is not None: self.l = l
        if k_l is not None: self.k_l = k_l
        if m_t is not None: self.m_t = m_t
        if m_p is not None: self.m_p = m_p
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if c is not None: self.c = c
        if u_e is not None: self.u_e = u_e
        if v is not None: self.v = v
        if c_phi is not None: self.c_phi = c_phi
        if c_l is not None: self.c_phi = c_l
        self.ivar = ivar
        self.qs = [self.x,self.phi,self.u_e]
        #self.qs = [self.x,self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        #self.u=self.l+self.v*self.ivar
        #self._TrolleyWithPendulum = TrolleyWithElasticPendulum(self.l+self.qs[2], self.m_t, self.m_p, self.k, self.g, self.Omega, self.qs[2], self.qs[0], self.F,ivar=self.ivar)(label='Trolley')
        self._TrolleyWithPendulum = TrolleyWithElasticPendulum(self.l, self.m_t, self.m_p, self.k, self.g, self.Omega, F=self.F,ivar=self.ivar)(label='Trolley')
        self._trolley_damper=Damper(self.c, pos1=self.x+(self.l+self.u_e)*sin(self.phi), qs=self.qs)(label='Trolley damper')
        self._pendulum_damper=Damper(self.c_phi, pos1=self.phi, qs=self.qs)(label='Pendulum damper')

        self._pendulum_line_damper_horizontal=Damper(self.c_l, pos1=self.x+self.u_e*sin(self.phi), qs=[self.u_e])(label='Pendulum horizontal line damper')
        self._pendulum_line_damper_vertical=Damper(self.c_l, pos1=self.u_e*cos(self.phi), qs=[self.u_e])(label='Pendulum vertical line damper')




        components['_TrolleyWithPendulum'] = self._TrolleyWithPendulum
        components['_trolley_damper'] = self._trolley_damper
        components['_pendulum_damper'] = self._pendulum_damper
        components['_pendulum_line_damper_horizontal'] = self._pendulum_line_damper_horizontal
        components['_pendulum_line_damper_vertical'] = self._pendulum_line_damper_vertical

        return components
