import base64
import inspect
import random

import IPython as IP
import numpy as np
import pandas as pd
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
    exp,
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

from ...dynamics import HarmonicOscillator, LagrangesDynamicSystem
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
    TorsionalSpring,
    base_frame,
    base_origin,
)
from ..mechanics.disk import RollingDisk
from .principles import (
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
    cached_property,
    lru_cache,
)

ureg = units


# DONE
class Pendulum(NonlinearComposedSystem):
    """
    Model of a sDoF mathematical Pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional's field acceleration

            l = length
                -Dimension of pendulum's strong

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """

    scheme_name = "undamped_pendulum.png"
    real_name = "pendulum_real.jpg"

    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    angle = dynamicsymbols("\phi")
    p_t = dynamicsymbols("p")
    qs = None

    m0 = Symbol("m_0", positive=True)
    l0 = Symbol("l_0", positive=True)

    def __init__(
        self, m=None, g=None, l=None, angle=None, p_t=None, ivar=None, **kwargs
    ):

        # if qs == None:
        #    qs = [angle]
        # else:
        #    qs = qs

        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if angle is not None:
            self.angle = angle
        if p_t is not None:
            self.p_t = p_t
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.angle]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self.gravitationalforce = GravitationalForce(
            self.m, self.g, self.l * (1 - cos(self.angle)), qs=self.qs
        )
        self.material_point = MaterialPoint(
            self.m, pos1=self.angle * self.l, qs=self.qs
        )

        components["Gravitational_Force"] = self.gravitationalforce
        components["Material_Point"] = self.material_point

        return components

    def get_default_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.angle: r"angle of the pendulum",
            self.angle.diff(self.ivar): r"velocity of the pendulum",
            self.angle.diff(self.ivar, 2): r"acceleration of the pendulum",
            self.ivar: r"time",
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [*REPORT_COMPONENTS_LIST]

        return comp_list

    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.l: ureg.meter,
            self.angle: ureg.radian,
            self.angle.diff(self.ivar): ureg.radian / ureg.second,
            self.angle.diff(self.ivar, 2): ureg.radian / ureg.second / ureg.second,
            self.ivar: ureg.second,
        }

        return units_dict

    def canonical_transformation(self):

        ham_approx = self.approximated(3).hamiltonian()

        p = Symbol(
            "p_" + "{" + str(self.q[0]).replace("(" + str(self.ivar) + ")", "") + "}"
        )

        dpdt = -ham_approx.rhs.subs(p, self.p_t).diff(self.q)[0]
        dqdt = ham_approx.rhs.subs(p, self.p_t).diff(self.p_t)

        return dpdt, dqdt


# DONE
class PulledPendulum(Pendulum):
    """
    Model of a sDoF mathematical Pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional's field acceleration

            l = lenght
                -Dimension of pendulum's strong

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') --> Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """

    scheme_name = "undamped_excited_pendulum.PNG"
    real_name = "pendulum_real.jpg"

    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    angle = dynamicsymbols("varphi")
    F = Symbol("F", positive=True)
    qs = None
    ivar = Symbol("t")
    Omega = Symbol("Omega", positive=True)

    def __init__(
        self,
        m=None,
        g=None,
        l=None,
        angle=None,
        ivar=None,
        F=None,
        Omega=None,
        **kwargs
    ):

        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if angle is not None:
            self.angle = angle
        if ivar is not None:
            self.ivar = ivar
        if F is not None:
            self.F = F
        if Omega is not None:
            self.Omega = Omega

        self.qs = [self.angle]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._pendulum = Pendulum(self.m, self.g, self.l, self.angle, self.ivar)(
            label="Pendulum"
        )
        self._force = Force(
            -self.F * self.l * sin(self.Omega * self.ivar),
            pos1=self.angle * self.l,
            qs=self.qs,
        )

        components["_pendulum"] = self._pendulum
        components["_force"] = self._force

        return components

    def get_default_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def get_numerical_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.l: [self.l],
            self.m: [self.m],
            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
            self.Omega: [50],
            self.F: [100],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
        }
        return self.sym_desc_dict

    def nonlinear_steady_solution(self):
        steady = self.eoms.steady_solution
        return steady

    def static_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        free_coeff = ans.subs(
            {cos(self.Omega * self.ivar): 0, sin(self.Omega * self.ivar): 0}
        ).subs(data)
        return free_coeff

    def max_static_cable_force(self, op_point=0):
        return abs(self.static_cable_force(op_point=op_point))

    def max_dynamic_cable_force(self, op_point=0):

        op_point = 0
        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        cos_amp = ans.subs(
            {cos(self.Omega * self.ivar): 1, sin(self.Omega * self.ivar): 0}
        ).subs(data)

        return abs(cos_amp)  # + self.max_static_cable_force()

    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def force_in_cable(self, op_point=0):

        op_point = 0

        data = self._given_data
        dyn_sys = self.subs(data)
        dyn_sys_lin = dyn_sys.linearized()
        phi = dyn_sys_lin.eoms.steady_solution[0]

        #         m=data[self.m]
        #         l=data[self.l]

        op_point = pi * op_point  # quick workaround - wrong implementation

        force_in_cable = (
            self.m * self.g * (1 - S.One / 2 * (phi - op_point) ** 2)
            + self.m * self.l * (phi - op_point).diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(
            data
        )  # .subs({self.Omega:0.999*dyn_sys_lin.natural_frequencies()[0]})

        return force_subs.doit().expand()

    def small_parameter(self):
        for n in self.approximated()._eoms:
            t = self.ivar
            coef = n.coeff((self.angle.diff(t, t)))
            smallparameter = n.coeff(self.angle**3) / coef
        return smallparameter

    def freq_amp_rel(self):
        A = Symbol("A", positive=True)
        free_vib_freq = self.natural_frequencies()
        return free_vib_freq**2 + S.Half * 3 / 2 * A**2 * self.small_parameter()

    @property
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST,
            mech_comp.SmallParameterComponent,
            mech_comp.StaticAndDynamicCableForceComponent,
            mech_comp.DynamicTensionForceComponent,
            mech_comp.DynamicCableDiameter,
            mech_comp.NonlinearSteadySolutionComponent,
            mech_comp.AmplitudeAndFrequencyRelationComponent,
        ]
        return comp_list


# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie


# DONE
class FreePendulum(Pendulum):
    """
    Model of a sDoF free pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> FreePendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SFoDFreePendulum()
    """

    scheme_name = "free_sdof_pendulum.png"
    real_name = "pendulum_real.jpg"


# 152
# ISSUE #152
class ExcitedPendulum(ComposedSystem):
    """
    Model of a sDoF Excited Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            F = Force
                -Pendulum's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l, F = symbols('m, g, l, F')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> SDoFExcitedPendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SDoFExcitedPendulum()
    """

    scheme_name = "damped_excited_pendulum.PNG"
    real_name = "pendulum2_real.jpg"

    def __init__(
        self,
        dummy=Symbol("dummy", positive=True),
        m=Symbol("m", positive=True),
        g=Symbol("g", positive=True),
        l=Symbol("l", positive=True),
        F=Symbol("F", positive=True),
        angle=dynamicsymbols("varphi"),
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):
        phi = angle
        self.phi = phi

        #         if qs == None:
        #             qs = [angle]
        #         else:
        #             qs = qs

        self.qs = [angle]

        self.dummy = dummy
        self.m = m
        self.g = g
        self.l = l
        self.F = F
        self.angle = angle
        self.ivar = ivar

        Omega = Symbol("Omega", positive=True)
        self.Omega = Omega

        self._init_from_components(**kwargs)

    @property
    def components(self):

        self.pendulum = Pendulum(self.m, self.g, self.l, angle=self.phi)
        self.force = Force(
            -self.F * self.l * sin(self.Omega * self.ivar), pos1=self.phi, qs=self.qs
        )
        # system = self.pendulum + self.force

        components = {}

        components["Pendulum"] = self.pendulum
        components["Force"] = self.force

        return components

    def get_default_data(self):

        m0, l0, g0, F0 = symbols("m_0 l_0 g_0 F_0", positive=True)

        default_data_dict = {
            self.m: [2 * m0, 1 * m0, S.Half * m0, S.Half**2 * m0, 3 * S.Half * m0],
            self.l: [2 * l0, 1 * l0, S.Half * l0, S.Half**2 * l0, 3 * S.Half * l0],
            self.g: [g0],
            self.F: [2 * F0, 1 * F0, S.Half * F0, S.Half**2 * F0, 3 * S.Half * F0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.F: r"Force",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0 = symbols("m_0 l_0 F_0", positive=True)

        default_data_dict = {
            self.m1: [
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
                10 * m0,
                11 * m0,
                12 * m0,
                13 * m0,
                14 * m0,
                15 * m0,
                16 * m0,
                17 * m0,
                18 * m0,
                19 * m0,
                20 * m0,
                21 * m0,
                22 * m0,
                23 * m0,
                24 * m0,
                25 * m0,
                26 * m0,
                27 * m0,
                28 * m0,
                29 * m0,
                30 * m0,
            ],
            self.l: [
                1 * l0,
                2 * l0,
                3 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
                10 * l0,
                11 * l0,
                12 * l0,
                13 * l0,
                14 * l0,
                15 * l0,
                16 * l0,
                17 * l0,
                18 * l0,
                19 * l0,
                20 * l0,
                21 * l0,
                22 * l0,
                23 * l0,
                24 * l0,
                25 * l0,
                26 * l0,
                27 * l0,
                28 * l0,
                29 * l0,
                30 * l0,
            ],
            self.F: [
                1 * F0,
                2 * F0,
                3 * F0,
                4 * F0,
                5 * F0,
                6 * F0,
                7 * F0,
                8 * F0,
                9 * F0,
                10 * F0,
                11 * F0,
                12 * F0,
                13 * F0,
                14 * F0,
                15 * F0,
                16 * F0,
                17 * F0,
                18 * F0,
                19 * F0,
                20 * F0,
                21 * F0,
                22 * F0,
                23 * F0,
                24 * F0,
                25 * F0,
                26 * F0,
                27 * F0,
                28 * F0,
                29 * F0,
                30 * F0,
            ],
        }
        return default_data_dict


# 151
class DampedPendulum(ComposedSystem):
    """
    Model of a sDoF damped Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            c = damper coefficient
                -value of damper coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l, c = symbols('m, g, l, c')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> SDoFDampedPendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFDampedPendulum()
    """

    scheme_name = "damped_pendulum.png"
    real_name = "pendulum2_real.jpg"

    def __init__(
        self,
        m=Symbol("m", positive=True),
        g=Symbol("g", positive=True),
        l=Symbol("l", positive=True),
        c=Symbol("c", positive=True),
        angle=dynamicsymbols("varphi"),
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):
        phi = angle
        self.phi = phi

        #         if qs == None:
        #             qs = [angle]
        #         else:
        #             qs = qs

        self.qs = [angle]

        self.m = m
        self.g = g
        self.l = l
        self.c = c

        self._init_from_components(**kwargs)

    @property
    def components(self):

        self.Pendulum = Pendulum(self.m, self.g, self.l, angle=self.phi)
        self.Damper = Damper(self.c, self.l * self.phi, qs=self.qs)
        # system = self.Pendulum + self.Damper

        components = {}

        components["Pendulum"] = self.Pendulum
        components["Damper"] = self.Damper

        return components

    def get_default_data(self):

        m0, l0, c0 = symbols("m_0 l_0 c_0", positive=True)

        default_data_dict = {
            self.m: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.l: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.c: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.c: r"Damping coefficient",
        }
        return self.sym_desc_dict


class ExcitedDampedPendulum(ComposedSystem):

    scheme_name = "damped_excited_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    def __init__(
        self,
        m=Symbol("m", positive=True),
        g=Symbol("g", positive=True),
        l=Symbol("l", positive=True),
        c=Symbol("c", positive=True),
        F=Symbol("F", positive=True),
        Omega=Symbol("Omega", positive=True),
        angle=dynamicsymbols("varphi"),
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):
        phi = angle
        self.phi = phi

        #         if qs == None:
        #             qs = [angle]
        #         else:
        #             qs = qs

        self.qs = [angle]

        self.m = m
        self.g = g
        self.l = l
        self.c = c
        self.F = F
        self.Omega = Omega

        # super().__init__(system, **kwargs)
        self._init_from_components(**kwargs)

    @property
    def components(self):

        self.Pendulum = Pendulum(self.m, self.g, self.l, angle=self.phi)
        self.Damper = Damper(self.c, self.l * self.phi, qs=self.qs)
        self.Force = Force(
            -self.F * sin(self.Omega * self.ivar), pos1=self.phi, qs=[self.phi]
        )
        # system = self.Pendulum + self.Damper + self.Force

        components = {}

        components["Pendulum"] = self.Pendulum
        components["Damper"] = self.Damper
        components["Force"] = self.Force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.c: r"Damping coefficient",
        }
        return self.sym_desc_dict


# DONE 160
# ISSUE #160
class PendulumKinematicExct(ComposedSystem):

    scheme_name = "kin_exct_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    phi = dynamicsymbols("\\varphi")
    x_e = dynamicsymbols("x_e")
    x_0 = Symbol("x_0", positive=True)

    def __init__(
        self, l=None, m=None, g=None, phi=None, x_e=None, ivar=Symbol("t"), **kwargs
    ):

        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x_e is not None:
            self.x_e = x_e
        self.ivar = ivar
        self.qs = [self.phi]

        self.x = self.l * sin(self.phi) + self.x_e
        self.y = self.l * cos(self.phi)

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._material_point_1 = MaterialPoint(self.m, self.x, qs=self.qs)
        self._material_point_2 = MaterialPoint(self.m, self.y, qs=self.qs)
        self._gravity = GravitationalForce(self.m, self.g, pos1=-self.y, qs=self.qs)

        components["_material_point_1"] = self._material_point_1
        components["_material_point_2"] = self._material_point_2
        components["_gravity"] = self._gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.x_e: r"Kinematic lateral excitation",
            self.m: r"Mass",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, x0, Omega = symbols("m_0 l_0 x_0 Omega", positive=True)

        default_data_dict = {
            self.m: [
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
                10 * m0,
                11 * m0,
                12 * m0,
                13 * m0,
                14 * m0,
                15 * m0,
                16 * m0,
                17 * m0,
                18 * m0,
                19 * m0,
                20 * m0,
                21 * m0,
                22 * m0,
                23 * m0,
                24 * m0,
                25 * m0,
                26 * m0,
                27 * m0,
                28 * m0,
                29 * m0,
                30 * m0,
            ],
            self.l: [
                1 * l0,
                2 * l0,
                3 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
                10 * l0,
                11 * l0,
                12 * l0,
                13 * l0,
                14 * l0,
                15 * l0,
                16 * l0,
                17 * l0,
                18 * l0,
                19 * l0,
                20 * l0,
                21 * l0,
                22 * l0,
                23 * l0,
                24 * l0,
                25 * l0,
                26 * l0,
                27 * l0,
                28 * l0,
                29 * l0,
                30 * l0,
            ],
            self.x_e: [x0 * sin(self.Omega * self.ivar)],
            x0: [l0 / (n * 20) for n in range(1, 11)],
        }
        return default_data_dict

    def static_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        free_coeff = ans.subs(
            {cos(self.Omega * self.ivar): 0, sin(self.Omega * self.ivar): 0}
        ).subs(data)
        return free_coeff

    def max_static_cable_force(self, op_point=0):
        return abs(self.static_cable_force(op_point=op_point))

    def max_dynamic_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        cos_amp = ans.subs(
            {cos(self.Omega * self.ivar): 1, sin(self.Omega * self.ivar): 0}
        ).subs(data)

        return abs(cos_amp)  # + self.max_static_cable_force())

    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def force_in_cable(self, op_point=0):

        data = self._given_data
        dyn_sys = self.subs(data)
        dyn_sys_lin = dyn_sys.linearized(op_point)
        phi = dyn_sys_lin._fodes_system.steady_solution[0]

        #         m=data[self.m]
        #         l=data[self.l]

        op_point = pi * op_point  # quick workaround - wrong implementation

        force_in_cable = (
            self.m * self.g * (1 - S.One / 2 * (phi - op_point) ** 2)
            + self.m * self.l * (phi - op_point).diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(
            data
        )  # .subs({self.Omega:0.999*dyn_sys_lin.natural_frequencies()[0]})

        return force_subs.doit().expand()

    def sin_coeff(self):

        phi = self.phi
        sin_coeff = (
            self._eoms[0].expand().coeff(sin(phi))
        )  ## eoms[0]: m * L**2 * phi" + m * g * L * sin(phi)  == 0

        return sin_coeff

    def cos_coeff(self):

        phi = self.phi
        cos_coeff = self._eoms[0].expand().coeff(cos(phi))

        return cos_coeff

    def equilibrium_position(self):

        x_e = self.x_e
        phi = self.phi
        equilibrium = (
            self.equilibrium_equation_demo(coordinates=True)
            .doit()
            .subs(x_e, 0)
            .doit()[0]
        )
        solution = solve(equilibrium, phi)
        solution += [i * solution[1] for i in range(2, 10)]

        return solution

    @property
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST,
            mech_comp.DynamicTensionForceComponent,
            mech_comp.StaticCableDiameter,
            mech_comp.DynamicCableDiameter,
        ]

        return comp_list


class KinematicallyExcitedInvertedPendulum(PendulumKinematicExct):

    scheme_name = "Inverse_pendulum.png"
    real_name = "elastic_pendulum_real.PNG"

    @cached_property
    def components(self):

        components = {}

        self._material_point_1 = MaterialPoint(self.m, self.x, qs=self.qs)
        self._material_point_2 = MaterialPoint(self.m, self.y, qs=self.qs)
        self._gravity = GravitationalForce(self.m, self.g, pos1=self.y, qs=self.qs)

        components["_material_point_1"] = self._material_point_1
        components["_material_point_2"] = self._material_point_2
        components["_gravity"] = self._gravity

        return components

    @property
    def _report_components(self):

        comp_list = [*REPORT_COMPONENTS_LIST, mech_comp.CriticalPointsComponent]
        return comp_list


class MDoFElasticPendulum(ComposedSystem):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    k = Symbol("k", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    z = dynamicsymbols("z")
    phi = dynamicsymbols("\\varphi")

    scheme_name = "elastic_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    def __init__(
        self,
        k=None,
        l=None,
        m=None,
        g=None,
        z=None,
        phi=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if k is not None:
            self.k = k
        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if z is not None:
            self.z = z
        self.ivar = ivar
        self.qs = [self.phi, self.z]

        self.x = (self.l + self.z) * sin(self.phi)
        self.y = (self.l + self.z) * cos(self.phi)

        #         self.frame = base_frame

        #         self.payload = Point('payload')

        #         self.payload.set_vel(self.frame,
        #         sqrt(diff(self.z, self.ivar)**2 + (diff(self.phi, self.ivar) * (l + self.z))**2) * self.frame(self.x))

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._spring = Spring(self.k, self.z, qs=[self.phi, self.z])
        self._material_point_1 = MaterialPoint(
            self.m, pos1=self.z, ivar=self.ivar, qs=[self.phi, self.z]
        )
        self._material_point_2 = MaterialPoint(
            self.m * (self.l + self.z) ** 2,
            pos1=self.phi,
            ivar=self.ivar,
            qs=[self.phi, self.z],
        )

        # self._material_point_1 = HarmonicOscillator(S.Half*self.m * (diff(self.z, self.ivar)**2 + (diff(self.phi, self.ivar) * (self.l + self.z))**2) , qs=[self.phi, self.z])

        # self._material_point_1 = HarmonicOscillator(S.Half*self.m * (diff(self.z, self.ivar)**2 + (diff(self.phi, self.ivar) * (self.l + self.z))**2) , qs=[self.phi, self.z])

        self._gravity = GravitationalForce(
            self.m, self.g, pos1=-self.y, qs=[self.phi, self.z]
        )

        components["_spring"] = self._spring
        components["_material_point_1"] = self._material_point_1
        components["_material_point_2"] = self._material_point_2
        components["_gravity"] = self._gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r"Spring stiffness",
            self.l: r"Pendulum length",
            self.m: r"Mass",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, k0 = symbols("m_0 l_0 k_0", positive=True)

        default_data_dict = {
            self.m: [1 * m0, S.Half**2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.l: [1 * l0, S.Half * l0, 2 * l0, 1 * l0, 2 * l0],
            self.k: [S.Half * k0, 1 * k0, 2 * k0, S.Half * k0, 2 * k0],
        }

        return default_data_dict

    def linearized(self, x0=None, op_point=True, hint=None, label=None):

        if hint is None:
            hint = [self.phi]

        temp_sys = HarmonicOscillator(Lagrangian=self.lagrangian, system=self)

        return temp_sys.linearized(x0=x0, op_point=True, hint=hint, label=label)

    @property
    def _report_components(self):

        comp_list = [*REPORT_COMPONENTS_LIST]
        return comp_list


class MDoFPendulumWithRod(ComposedSystem):
    """
    Model of a Double Degree of Freedom Pendulum

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    k = Symbol("k", positive=True)
    l = Symbol("l", positive=True)
    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    g = Symbol("g", positive=True)
    z = dynamicsymbols("z")
    phi = dynamicsymbols("\\varphi")

    scheme_name = "PendulumWithRod.png"
    # real_name = 'elastic_pendulum_real.PNG'

    def __init__(
        self,
        k=None,
        l=None,
        m1=None,
        m2=None,
        g=None,
        z=None,
        phi=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if k is not None:
            self.k = k
        if l is not None:
            self.l = l
        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if z is not None:
            self.z = z
        self.ivar = ivar
        self.qs = [self.phi, self.z]

        self.x = (self.l + self.z) * sin(self.phi)
        self.y = (self.l + self.z) * cos(self.phi)

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._spring = Spring(self.k, self.z, qs=[self.phi, self.z])
        self._material_point_1 = MaterialPoint(
            self.m2, pos1=self.z, ivar=self.ivar, qs=[self.phi, self.z]
        )
        self._material_point_2 = MaterialPoint(
            self.m2 * (self.l - self.z) ** 2,
            pos1=self.phi,
            ivar=self.ivar,
            qs=[self.phi, self.z],
        )
        self.rod = RigidBody2D(
            self.m1,
            self.m1 * self.l**2,
            pos_lin=self.z,
            pos_rot=self.phi,
            qs=[self.z, self.phi],
            frame=None,
            ivar=self.ivar,
        )
        self._gravity_rod = GravitationalForce(
            self.m1, self.g, pos1=-self.y, qs=[self.phi, self.z]
        )
        self._gravity_mass = GravitationalForce(
            self.m2, self.g, pos1=-self.y, qs=[self.phi, self.z]
        )

        components["_spring"] = self._spring
        components["_material_point_1"] = self._material_point_1
        components["_material_point_2"] = self._material_point_2
        components["_gravity_rod"] = self._gravity_rod
        components["_gravity_mass"] = self._gravity_mass
        components["rod"] = self.rod

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r"Spring stiffness",
            self.l: r"Pendulum length",
            self.m1: r"Mass of the rod",
            self.m2: r"Mass of the payload",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, k0 = symbols("m_0 l_0 k_0", positive=True)

        default_data_dict = {
            self.m1: [1 * m0, S.Half**2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m2: [1 * m0, S.Half**2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.l: [1 * l0, S.Half * l0, 2 * l0, 1 * l0, 2 * l0],
            self.k: [S.Half * k0, 1 * k0, 2 * k0, S.Half * k0, 2 * k0],
        }

        return default_data_dict

    def linearized(self, x0=None, op_point=True, hint=None, label=None):

        if hint is None:
            hint = [self.phi]

        temp_sys = HarmonicOscillator(Lagrangian=self.lagrangian, system=self)

        return temp_sys.linearized(x0=x0, op_point=True, hint=hint, label=label)


# DONE KrzysiekK
class PendulumWithRoller(NonlinearComposedSystem):
    """
    Model of a sDoF mathematical Pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.

        Arguments:
        =========
            m1 = Mass
                -Mass of the gear

            m2 = Mass
                -Mass of the carrier

            g = gravitional field
                -value of gravitional's field acceleration

            R = lenght
                -Length of the carrier

            r = length
                - Radius of the gear

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

    """

    scheme_name = "PlanetaryGearWithSpring.png"

    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    g = Symbol("g", positive=True)
    R = Symbol("R", positive=True)
    r = Symbol("r", positive=True)
    k = Symbol("k", positive=True)
    angle = dynamicsymbols("varphi")
    qs = None

    m0 = Symbol("m_0", positive=True)
    l0 = Symbol("l_0", positive=True)

    def __init__(
        self,
        m1=None,
        m2=None,
        k=None,
        R=None,
        r=None,
        g=None,
        angle=None,
        ivar=None,
        **kwargs
    ):

        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if k is not None:
            self.k = k
        if g is not None:
            self.g = g
        if r is not None:
            self.r = r
        if R is not None:
            self.R = R
        if angle is not None:
            self.angle = angle
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.angle]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.gravitationalforce_rod = GravitationalForce(
            self.m2, self.g, (self.R - self.r) * (1 - cos(self.angle)), qs=[self.angle]
        )
        self.gravitationalforce_gear = GravitationalForce(
            self.m1, self.g, self.r * (1 - cos(self.angle)), qs=[self.angle]
        )
        self.Rod = RigidBody2D(
            self.m2,
            self.m2 * (self.R - self.r) ** 2,
            pos_lin=self.angle * (self.R - self.r),
            pos_rot=self.angle,
            qs=[self.angle],
            frame=None,
            ivar=self.ivar,
        )
        self.Gear = RollingDisk(
            self.m1,
            self.r,
            F=None,
            x=self.angle * (self.R - self.r),
            qs=[self.angle],
            ivar=self.ivar,
        )
        self.Spring = Spring(self.k, self.angle, qs=[self.angle])

        components["Gravitational_Force_Rod"] = self.gravitationalforce_rod
        components["Gravitational_Force_Gear"] = self.gravitationalforce_gear
        components["Rod"] = self.Rod
        components["Gear"] = self.Gear
        components["Spring"] = self.Spring

        return components

    def get_default_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
        }
        return self.sym_desc_dict


# DONE #128
# ISSUE #128
class MDoFLinearizedThreePendulumsWithSprings(ComposedSystem):
    scheme_name = "three_pendulums_forced.PNG"
    real_name = "lifting_tandem.png"

    def __init__(
        self,
        R=Symbol("R", positive=True),
        phi_1=dynamicsymbols("\\varphi_1"),
        phi_2=dynamicsymbols("\\varphi_2"),
        phi_3=dynamicsymbols("\\varphi_3"),
        phi_l=dynamicsymbols("\\varphi_l"),
        phi_c=dynamicsymbols("\\varphi_c"),
        phi_r=dynamicsymbols("\\varphi_r"),
        m1=Symbol("m_1", positive=True),
        m2=Symbol("m_2", positive=True),
        m3=Symbol("m_3", positive=True),
        k_1=Symbol("k_1", positive=True),
        k_2=Symbol("k_2", positive=True),
        k_3=Symbol("k_3", positive=True),
        k_4=Symbol("k_4", positive=True),
        l=Symbol("l", positive=True),
        g=Symbol("g", positive=True),
        Omega=Symbol("Omega", positive=True),
        F=Symbol("F", positive=True),
        qs=dynamicsymbols("\\varphi_l, \\varphi_c, \\varphi_r"),
        ivar=Symbol("t"),
        **kwargs
    ):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_1 = k_1
        self.k_2 = k_2
        self.k_3 = k_3
        self.k_4 = k_4
        self.l = l
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.phi_l = phi_l
        self.phi_c = phi_c
        self.phi_r = phi_r
        self.g = g
        self.Omega = Omega
        self.F = F
        self.qs = qs

        # super().__init__(system(qs),**kwargs)
        self._init_from_components(**kwargs)

    @property
    def components(self):

        self.Pendulum1 = Pendulum(
            self.m1, self.g, self.l, angle=self.phi_l, qs=[self.phi_l]
        ).linearized() + Force(
            self.F * self.l * cos(self.Omega * self.ivar),
            pos1=self.phi_l,
            qs=[self.phi_l],
        )
        self.Pendulum2 = Pendulum(
            self.m2, self.g, self.l, angle=self.phi_c, qs=[self.phi_c]
        ).linearized()
        self.Pendulum3 = Pendulum(
            self.m3, self.g, self.l, angle=self.phi_r, qs=[self.phi_r]
        ).linearized() + Force(
            2 * self.F * self.l * cos(self.Omega * self.ivar),
            pos1=self.phi_r,
            qs=[self.phi_r],
        )
        self.Spring1 = Spring(
            self.k_1,
            pos1=(self.phi_l * (self.l / 2)),
            pos2=(self.phi_c * (self.l / 2)),
            qs=[self.phi_l, self.phi_c],
        )
        self.Spring2 = Spring(
            self.k_2,
            pos1=(self.phi_c * (self.l / 2)),
            pos2=(self.phi_r * (self.l / 2)),
            qs=[self.phi_c, self.phi_r],
        )
        self.Spring3 = Spring(
            self.k_3,
            pos1=(self.phi_l * self.l),
            pos2=(self.phi_c * self.l),
            qs=[self.phi_l, self.phi_c],
        )
        self.Spring4 = Spring(
            self.k_4,
            pos1=(self.phi_c * self.l),
            pos2=(self.phi_r * self.l),
            qs=[self.phi_c, self.phi_r],
        )

        # system = self.Pendulum1 + self.Pendulum2 + self.Pendulum3 + self.Spring1 + self.Spring2 + self.Spring3 + self.Spring4

        components = {}

        components["Pendulum1"] = self.Pendulum1
        components["Pendulum2"] = self.Pendulum2
        components["Pendulum3"] = self.Pendulum3
        components["Spring1"] = self.Spring1
        components["Spring2"] = self.Spring2
        components["Spring3"] = self.Spring3
        components["Spring4"] = self.Spring4

        return components

    def get_default_data(self):

        m0, k0, l0, F0 = symbols("m_0 k_0 l_0 F_0", positive=True)

        default_data_dict = {
            #             self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            #             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            #             self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m1: [m0 * no for no in range(1, 9)],
            self.m2: [m0 * no for no in range(1, 9)],
            self.m3: [m0 * no for no in range(1, 9)],
            self.l: [l0 * no for no in range(1, 9)],
            self.F: [F0 * no for no in range(1, 9)],
            self.k_1: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.phi_l: [self.phi_1, 0],
            self.phi_c: [self.phi_1, self.phi_2],
            self.phi_r: [self.phi_2, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if (
            parameters_dict[self.phi_l] != self.phi_1
            or parameters_dict[self.phi_c] != self.phi_1
        ):

            parameters_dict[self.phi_l] = self.phi_1

        if (
            parameters_dict[self.phi_c] != self.phi_2
            or parameters_dict[self.phi_r] != self.phi_2
        ):

            parameters_dict[self.phi_r] = self.phi_2

        return parameters_dict

    def max_static_cable_force(self):
        return (self.m1 * self.g).subs(self._given_data)

    def max_dynamic_cable_force(self):

        omg_amp = ComposedSystem(self.linearized())._frf()[0] * self.Omega

        return self.m1 * self.l * (omg_amp) ** 2 + self.max_static_cable_force()


class ForcedTriplePendulum(ComposedSystem):

    scheme_name = "forced_triple_pendulum.png"
    real_name = "TriplePendulum_real.jpg"

    m = Symbol("m", positive=True)
    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    m3 = Symbol("m_3", positive=True)
    l1 = Symbol("l_1", positive=True)
    l2 = Symbol("l_2", positive=True)
    l3 = Symbol("l_3", positive=True)
    F1 = Symbol("F_1", positive=True)
    F2 = Symbol("F_2", positive=True)
    F3 = Symbol("F_3", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    ivar = Symbol("t")

    phi1 = dynamicsymbols("\\varphi_1")
    phi2 = dynamicsymbols("\\varphi_2")
    phi3 = dynamicsymbols("\\varphi_3")
    phi_u = dynamicsymbols("\\varphi_u")
    phi_l = dynamicsymbols("\\varphi_l")
    phi = dynamicsymbols("\\varphi")

    def __init__(
        self,
        m=None,
        m1=None,
        m2=None,
        m3=None,
        l1=None,
        l2=None,
        l3=None,
        F1=None,
        F2=None,
        F3=None,
        g=None,
        phi1=None,
        phi2=None,
        phi3=None,
        phi_u=None,
        phi_l=None,
        Omega=None,
        ivar=None,
        **kwargs
    ):

        if m is not None:
            self.m = m
        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if m3 is not None:
            self.m3 = m3
        if F1 is not None:
            self.F1 = F1
        if F2 is not None:
            self.F2 = F2
        if F3 is not None:
            self.F3 = F3
        if l1 is not None:
            self.l1 = l1
        if l2 is not None:
            self.l2 = l2
        if l3 is not None:
            self.l3 = l3
        if g is not None:
            self.g = g
        if Omega is not None:
            self.Omega = Omega
        if phi1 is not None:
            self.phi1 = phi1
        if phi2 is not None:
            self.phi2 = phi2
        if phi3 is not None:
            self.phi3 = phi3
        if phi_u is not None:
            self.phi_u = phi_u
        if phi_l is not None:
            self.phi_l = phi_l
        self.qs = [self.phi1, self.phi2, self.phi3]

        self.x_2 = sin(self.phi1) * self.l1 + sin(self.phi2) * self.l2
        self.y_2 = cos(self.phi1) * self.l1 + cos(self.phi2) * self.l2
        self.x_3 = self.x_2 + sin(self.phi3) * self.l3
        self.y_3 = self.y_2 + cos(self.phi3) * self.l3

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._Pendulum1 = Pendulum(self.m1, self.g, self.l1, self.phi1, self.ivar)
        self._material_point_11 = MaterialPoint(
            self.m2, self.x_2, qs=[self.phi1, self.phi2]
        )
        self._material_point_21 = MaterialPoint(
            self.m2, self.y_2, qs=[self.phi1, self.phi2]
        )
        self._gravity_1 = GravitationalForce(
            self.m2, self.g, pos1=-self.y_2, qs=[self.phi2]
        )
        self._material_point_12 = MaterialPoint(
            self.m3, self.x_3, qs=[self.phi1, self.phi2, self.phi3]
        )
        self._material_point_22 = MaterialPoint(
            self.m3, self.y_3, qs=[self.phi1, self.phi2, self.phi3]
        )
        self._gravity_2 = GravitationalForce(
            self.m3, self.g, pos1=-self.y_3, qs=[self.phi3]
        )
        self._Force1 = Force(
            self.F1 * self.l1 * sin(self.Omega * self.ivar),
            pos1=self.phi1,
            qs=[self.phi1],
        )
        self._Force2 = Force(
            self.F2 * self.l2 * sin(self.Omega * self.ivar),
            pos1=self.phi2,
            qs=[self.phi1, self.phi2],
        )
        self._Force3 = Force(
            self.F3 * self.l3 * sin(self.Omega * self.ivar),
            pos1=self.phi3,
            qs=[self.phi1, self.phi2, self.phi3],
        )

        components["_Pendulum1"] = self._Pendulum1
        components["_material_point_11"] = self._material_point_11
        components["_material_point_21"] = self._material_point_21
        components["_gravity_1"] = self._gravity_1
        components["_material_point_12"] = self._material_point_12
        components["_material_point_22"] = self._material_point_22
        components["_gravity_2"] = self._gravity_2
        components["_Force1"] = self._Force1
        components["_Force2"] = self._Force2
        components["_Force3"] = self._Force3

        return components

    def get_default_data(self):

        m0, l0, F0 = symbols("m_0 l_0 F", positive=True)

        default_data_dict = {
            self.m1: [m0 * no for no in range(1, 9)],
            self.m2: [m0 * no for no in range(1, 9)],
            self.m3: [m0 * no for no in range(1, 9)],
            self.l1: [l0 * no for no in range(1, 9)],
            self.l2: [l0 * no for no in range(1, 9)],
            self.l3: [l0 * no for no in range(1, 9)],
            self.phi1: [self.phi_u, 0],
            self.phi2: [self.phi_u, self.phi_l],
            self.phi3: [self.phi_l],
            self.F1: [F0 * no for no in range(1, 2)],
            self.F2: [F0 * no for no in range(1, 2)],
            self.F3: [F0 * no for no in range(1, 2)],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi2] == parameters_dict[self.phi3]:

            parameters_dict[self.phi2] = self.phi_u

        #         display(parameters_dict)
        return parameters_dict

    @property
    def l(self):

        if self._given_data[self.phi1] == 0:
            l1 = 0
        else:
            l1 = self.l1

        if self._given_data[self.phi2] == 0:
            l2 = 0
        else:
            l2 = self.l2

        return (
            l1 * self.m1 + (l1 + l2) * self.m2 + (l1 + l2 + self.l3) * self.m3
        ) / self.m

    @property
    def m(self):

        return self.m1 + self.m2 + self.m3

    def nonlinear_steady_solution(self):
        steady = self._fodes_system.steady_solution
        return steady

    def static_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        free_coeff = ans.subs(
            {cos(self.Omega * self.ivar): 0, sin(self.Omega * self.ivar): 0}
        ).subs(data)
        return free_coeff

    def max_static_cable_force(self, op_point=0):
        return abs(self.static_cable_force(op_point=op_point))

    def max_dynamic_cable_force(self, op_point=0):

        op_point = 0
        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        cos_amp = ans.subs(
            {cos(self.Omega * self.ivar): 1, sin(self.Omega * self.ivar): 0}
        ).subs(data)

        return abs(cos_amp)  # + self.max_static_cable_force()

    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def force_in_cable(self, op_point=0):

        op_point = 0

        data = self._given_data
        dyn_sys = self.subs(data)
        #         display(type(dyn_sys))
        dyn_sys_lin = dyn_sys.linearized()
        #         display(type(dyn_sys_lin))
        phi = dyn_sys_lin._fodes_system.steady_solution[0]

        #         m=data[self.m]
        #         l=data[self.l]

        op_point = pi * op_point  # quick workaround - wrong implementation

        force_in_cable = (
            self.m * self.g * (1 - S.One / 2 * (phi - op_point) ** 2)
            + self.m * self.l * (phi - op_point).diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(
            data
        )  # .subs({self.Omega:0.999*dyn_sys_lin.natural_frequencies()[0]})

        return force_subs.doit().expand()


class SDoFForcedTriplePendulum(ForcedTriplePendulum):

    def get_default_data(self):

        m0, l0, F0 = symbols("m_0 l_0 F", positive=True)

        default_data_dict = {
            self.m1: [m0 * no for no in range(1, 9)],
            self.m2: [m0 * no for no in range(1, 9)],
            self.m3: [m0 * no for no in range(1, 9)],
            self.l1: [l0 * no for no in range(1, 4)],
            self.l2: [l0 * no for no in range(1, 4)],
            self.l3: [l0 * no for no in range(1, 4)],
            # self.phi1: [self.phi,S.Zero],
            self.phi1: [S.Zero],
            self.phi2: [self.phi, S.Zero],
            self.phi3: [self.phi],
            self.F1: [F0 * no for no in range(1, 3)],
            self.F2: [F0 * no for no in range(1, 3)],
            self.F3: [F0 * no for no in range(1, 3)],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi1] != parameters_dict[self.phi2]:

            parameters_dict[self.phi1] = self.phi
            parameters_dict[self.phi2] = self.phi

        #         display(parameters_dict)
        return parameters_dict


class Winch(ComposedSystem):

    scheme_name = "ForcedWinch.png"
    real_name = "winch_mechanism_real.PNG"

    r = Symbol("r", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    ivar = Symbol("t")
    phi = dynamicsymbols("\\varphi")
    F = Symbol("F", positive=True)
    Omega = Symbol("Omega", positive=True)

    def __init__(
        self,
        r=None,
        l=None,
        m=None,
        g=None,
        phi=None,
        F=None,
        Omega=None,
        ivar=Symbol("t"),
        **kwargs
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
        if Omega is not None:
            self.Omega = Omega
        if F is not None:
            self.F = F
        self.ivar = ivar

        self.qs = [self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.x = self.r * cos(self.phi) + (self.l + self.r * self.phi) * sin(self.phi)
        self.y = -(self.r) * sin(self.phi) + (self.l + self.r * self.phi) * cos(
            self.phi
        )

        self.material_point_1 = MaterialPoint(self.m, self.x, qs=self.qs)
        self.material_point_2 = MaterialPoint(self.m, self.y, qs=self.qs)
        self.gravity = GravitationalForce(self.m, self.g, pos1=-self.y, qs=self.qs)
        self.force = Force(
            self.F * self.l * sin(self.Omega * self.ivar), self.phi, self.qs
        )

        components["_material_point_1"] = self.material_point_1
        components["_material_point_2"] = self.material_point_2
        components["_gravity"] = self.gravity
        components["_force"] = self.force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r"Winch radius",
            self.l: r"Winch length",
            self.m: r"Mass",
            self.g: "Gravity constant",
            self.F: "Force",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0 = symbols("m_0 l_0 F_0", positive=True)

        default_data_dict = {
            self.r: [S.One * no * 2 * l0 for no in range(1, 5)],
            self.m: [S.One * no * m0 for no in range(2, 10)],
            self.l: [S.One * no * l0 for no in range(1, 5)],
            self.F: [S.One * no * F0 for no in range(5, 10)],
        }
        return default_data_dict

    def force_in_cable(self, op_point=0):

        data = self._given_data
        dyn_sys = self.subs(data)
        dyn_sys_lin = dyn_sys.linearized()
        phi = dyn_sys_lin._fodes_system.steady_solution[0]

        #         m=data[self.m]
        #         l=data[self.l]

        op_point = pi * op_point  # quick workaround - wrong implementation

        force_in_cable = (
            self.m * self.g * (1 - S.Half * (phi - op_point) ** 2)
            + self.m * self.l * (phi - op_point).diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(
            data
        )  # .subs({self.Omega:0.999*dyn_sys_lin.natural_frequencies()[0]})

        return force_subs.doit().expand()

    def sin_coeff(self):

        phi = self.phi
        sin_coeff = (
            self._eoms[0].expand().coeff(sin(phi))
        )  ## eoms[0]: m * L**2 * phi" + m * g * L * sin(phi)  == 0

        return sin_coeff

    def cos_coeff(self):

        phi = self.phi
        cos_coeff = self._eoms[0].expand().coeff(cos(phi))

        return cos_coeff


# Franek
class DDOFCoupledPendulum(ComposedSystem):
    """
    Model of a DDoF Coupled Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            k = spring coefficient
                -value of spring coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l, k = symbols('m, g, l, k')
        >>> qs = dynamicsymbols('varphi_1, varphi_2') # Generalized Coordinates
        >>> DDoFCouplePendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = "ddof_coupled_pendulum.png"
    real_name = "lifting_tandem.png"
    # phi_1,phi_2=dynamicsymbols('\\varphi_1, \\varphi_2')
    # _hint = [phi_1,phi_2]

    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    l = Symbol("l", positive=True)
    k = Symbol("k", positive=True)
    phi1 = dynamicsymbols("\\varphi_1")
    phi2 = dynamicsymbols("\\varphi_2")
    g = Symbol("g", positive=True)
    F1 = Symbol("F_1", positive=True)
    F2 = Symbol("F_2", positive=True)
    Omega = Symbol("Omega", positive=True)
    l0 = Symbol("l_0", positive=True)
    m0 = Symbol("m_0", positive=True)
    k0 = Symbol("k_0", positive=True)
    F0 = Symbol("F_0", positive=True)
    Omega0 = Symbol("Omega_0", positive=True)

    def __init__(
        self,
        m=None,
        l=None,
        k=None,
        phi1=None,
        phi2=None,
        g=None,
        F1=None,
        F2=None,
        Omega=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if m is not None:
            self.m = m
        if l is not None:
            self.l = l
        if k is not None:
            self.k = k
        if phi1 is not None:
            self.phi1 = phi1
        if phi2 is not None:
            self.phi2 = phi2
        if g is not None:
            self.g = g
        if F1 is not None:
            self.F1 = F1
        if F2 is not None:
            self.F2 = F2
        if Omega is not None:
            self.Omega = Omega

        self.qs = [self.phi1, self.phi2]
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._spring = Spring(
            self.k,
            pos1=(self.phi1 * self.l / 2),
            pos2=(self.phi2 * self.l / 2),
            qs=[self.qs],
        )
        self._pendulum_1 = Pendulum(
            self.m1, self.g, self.l, angle=self.phi1, qs=[self.qs]
        )
        self._pendulum_2 = Pendulum(
            self.m2, self.g, self.l, angle=self.phi2, qs=[self.qs]
        )
        self._mass_1 = MaterialPoint(self.m1 * self.l**2, self.phi1, qs=[self.qs])
        self._mass_2 = MaterialPoint(self.m2 * self.l**2, self.phi2, qs=[self.qs])
        self._force_1 = Force(
            self.F1 * self.l * sin(self.Omega * self.ivar), pos1=self.phi1
        )
        self._force_2 = Force(
            self.F2 * self.l * sin(self.Omega * self.ivar), pos1=self.phi2
        )

        components["_spring"] = self._spring
        components["_pendulum_1"] = self._pendulum_1
        components["_pendulum_2"] = self._pendulum_2
        components["_mass_1"] = self._mass_1
        components["_mass_2"] = self._mass_2
        components["_force_1"] = self._force_1
        components["_force_2"] = self._force_2

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.k: r"Stifness coefficient",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, l0, F0, Omega0 = self.m0, self.k0, self.l0, self.F0, self.Omega0

        default_data_dict = {
            self.m1: [m0 * no for no in range(1, 8)],
            self.m2: [m0 * no for no in range(1, 8)],
            #             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            #             self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.k: [k0 * no for no in range(1, 8)],
            #             self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            #             self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            #             self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.l: [l0 * no for no in range(1, 8)],
            #             self.phi_1: [self.phi_l, 0],
            #             self.phi_2: [self.phi_l, self.phi_r, 0],
            #             self.phi_3: [self.phi_r, 0],
            self.F1: [F0 * no for no in range(1, 8)],
            self.F2: [F0 * no for no in range(1, 8)],
            self.Omega: [Omega0 * no for no in range(1, 2)],
            self.l0: [self.m0 / self.k0 * self.g],
        }

        return default_data_dict

    def max_static_g_force(self):
        return (self.m1 * self.g).subs(self._given_data)

    def max_dynamic_cable_force(self):

        omg_amp = ComposedSystem(self.linearized())._frf()[0] * self.Omega

        return self.m1 * self.l * (omg_amp) ** 2 + self.max_static_cable_force()

    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)


# Franek
class DDOFLinearizedCoupledPendulum(ComposedSystem):
    """
    Model of a DDoF Linearized Coupled Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            k = spring coefficient
                -value of spring coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l, k = symbols('m, g, l, k')
        >>> qs = dynamicsymbols('varphi_1, varphi_2') # Generalized Coordinates
        >>> DDoFCouplePendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = "mdof_dpendulum.png"
    real_name = "lifting_tandem.png"

    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    k = Symbol("k", positive=True)
    phi1 = dynamicsymbols("\\varphi_1")
    phi2 = dynamicsymbols("\\varphi_2")

    def __init__(
        self,
        m=None,
        g=None,
        l=None,
        k=None,
        phi1=None,
        phi2=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if k is not None:
            self.k = k
        if phi1 is not None:
            self.phi1 = phi1
        if phi2 is not None:
            self.phi2 = phi2

        self.ivar = ivar
        self.qs = [self.phi1, self.phi2]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._spring = Spring(
            self.k,
            pos1=(self.phi1 * (self.l / 2)),
            pos2=(self.phi2 * (self.l / 2)),
            qs=[self.qs],
        )
        self._pendulum_1 = Pendulum(
            self.m, self.g, self.l, angle=self.phi1, qs=[self.qs]
        ).linearized()
        self._pendulum_2 = Pendulum(
            self.m, self.g, self.l, angle=self.phi2, qs=[self.qs]
        ).linearized()

        components["_spring"] = self._spring
        components["_pendulum_1"] = self._pendulum_1
        components["_pendulum_2"] = self._pendulum_2

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.k: r"Stifness coefficient",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, l0 = symbols("m_0 k_0 l_0", positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            #             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            #             self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.k: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            #             self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            #             self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            #             self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.l: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            #             self.phi_1: [self.phi_l, 0],
            #             self.phi_2: [self.phi_l, self.phi_r, 0],
            #             self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict


# Marcel_k(zmienione na nowy sposób, ale nie wiem czy sprawdzone dokładnie)
class DoublePendulum(ComposedSystem):

    scheme_name = "MDOFTriplePendulum.PNG"
    real_name = "TriplePendulum_real.jpg"
    detali_scheme_name = "MDOFTriplePendulum.PNG"
    detali_real_name = "TriplePendulum_real.jpg"

    m = Symbol("m", positive=True)
    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    l_1 = Symbol("l_1", positive=True)
    l_2 = Symbol("l_2", positive=True)
    g = Symbol("g", positive=True)
    omega = Symbol("Omega", positive=True)
    phi_1 = dynamicsymbols("varphi_1")
    phi_2 = dynamicsymbols("varphi_2")
    phi_u = dynamicsymbols("varphi_u")
    phi_l = dynamicsymbols("varphi_l")
    qs = dynamicsymbols("varphi_1 varphi_2")
    ivar = Symbol("t")

    def __init__(
        self,
        m=None,
        m1=None,
        m2=None,
        l_1=None,
        l_2=None,
        g=None,
        omega=None,
        phi_1=None,
        phi_2=None,
        phi_u=None,
        phi_l=None,
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if m is not None:
            self.m = m
        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if l_1 is not None:
            self.l_1 = l_1
        if l_2 is not None:
            self.l_2 = l_2
        if g is not None:
            self.g = g
        if omega is not None:
            self.omega = omega
        if phi_1 is not None:
            self.phi_1 = phi_1
        if phi_2 is not None:
            self.phi_2 = phi_2
        if phi_u is not None:
            self.phi_u = phi_u
        if phi_l is not None:
            self.phi_l = phi_l
        self.qs = [self.phi_1, self.phi_2]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self.x_2 = sin(self.phi_1) * self.l_1 + sin(self.phi_2) * self.l_2
        self.y_2 = cos(self.phi_1) * self.l_1 + cos(self.phi_2) * self.l_2

        self.Pendulum1 = Pendulum(
            self.m1, self.g, self.l_1, angle=self.phi_1, qs=[self.phi_1]
        )
        self.material_point_11 = MaterialPoint(
            self.m2, self.x_2, qs=[self.phi_1, self.phi_2]
        )
        self.material_point_21 = MaterialPoint(
            self.m2, self.y_2, qs=[self.phi_1, self.phi_2]
        )
        self.gravity_1 = GravitationalForce(
            self.m2, self.g, pos1=-self.y_2, qs=[self.phi_2]
        )

        components["Pendulum 1"] = self.Pendulum1
        components["material_point_11"] = self.material_point_11
        components["material_point_21"] = self.material_point_21
        components["gravity_1"] = self.gravity_1

        return components


# Grzes
class InvertedPendulumDDoF(ComposedSystem):

    scheme_name = "InvertedPendulumDDoF.PNG"
    real_name = "elastic_pendulum_real.PNG"

    M = Symbol("M", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    phi = dynamicsymbols("\\varphi")
    x = dynamicsymbols("x")
    x_0 = Symbol("x_0", positive=True)
    F = Symbol("F", positive=True)

    def __init__(
        self,
        M=None,
        l=None,
        m=None,
        g=None,
        phi=None,
        x=None,
        F=None,
        ivar=Symbol("t"),
        **kwargs
    ):
        if M is not None:
            self.M = M
        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if phi is not None:
            self.phi = phi
        if x is not None:
            self.x = x
        if F is not None:
            self.F = F
        self.ivar = ivar
        self.qs = [self.phi, self.x]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.M, self.x, qs=[self.x])
        self._pendulum = KinematicallyExcitedInvertedPendulum(
            self.l, self.m, self.g, self.phi, self.x, self.ivar
        )
        self._force = Force(
            self.F * sin(self.Omega * self.ivar), pos1=self.x, qs=[self.x, self.phi]
        )

        #         self._gravity = GravitationalForce(self.m, self.g, pos1=self.y, qs=self.qs)

        components["_trolley"] = self._trolley
        components["_pendulum"] = self._pendulum
        components["_force"] = self._force
        #         components['_gravity'] = self._gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r"Pendulum length",
            self.x: r"Kinematic lateral excitation",
            self.m: r"Mass",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, x0, Omega = symbols("m_0 l_0 x_0 Omega", positive=True)

        default_data_dict = {
            self.m: [
                1 * m0,
                2 * m0,
                3 * m0,
                4 * m0,
                5 * m0,
                6 * m0,
                7 * m0,
                8 * m0,
                9 * m0,
                10 * m0,
                11 * m0,
                12 * m0,
                13 * m0,
                14 * m0,
                15 * m0,
                16 * m0,
                17 * m0,
                18 * m0,
                19 * m0,
                20 * m0,
                21 * m0,
                22 * m0,
                23 * m0,
                24 * m0,
                25 * m0,
                26 * m0,
                27 * m0,
                28 * m0,
                29 * m0,
                30 * m0,
            ],
            self.l: [
                1 * l0,
                2 * l0,
                3 * l0,
                4 * l0,
                5 * l0,
                6 * l0,
                7 * l0,
                8 * l0,
                9 * l0,
                10 * l0,
                11 * l0,
                12 * l0,
                13 * l0,
                14 * l0,
                15 * l0,
                16 * l0,
                17 * l0,
                18 * l0,
                19 * l0,
                20 * l0,
                21 * l0,
                22 * l0,
                23 * l0,
                24 * l0,
                25 * l0,
                26 * l0,
                27 * l0,
                28 * l0,
                29 * l0,
                30 * l0,
            ],
            self.x_e: [x0 * sin(self.Omega * self.ivar)],
            x0: [l0 / (n * 20) for n in range(1, 11)],
        }
        return default_data_dict

    def static_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        free_coeff = ans.subs(
            {cos(self.Omega * self.ivar): 0, sin(self.Omega * self.ivar): 0}
        ).subs(data)
        return free_coeff

    def max_static_cable_force(self, op_point=0):
        return abs(self.static_cable_force(op_point=op_point))

    def max_dynamic_cable_force(self, op_point=0):

        data = self._given_data
        ans = self.force_in_cable(op_point=op_point)
        cos_amp = ans.subs(
            {cos(self.Omega * self.ivar): 1, sin(self.Omega * self.ivar): 0}
        ).subs(data)

        return abs(cos_amp)  # + self.max_static_cable_force())

    def static_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol("k_r", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re)) ** (1 / 2)

    def force_in_cable(self, op_point=0):

        data = self._given_data
        dyn_sys = self.subs(data)
        dyn_sys_lin = dyn_sys.linearized(op_point)
        phi = dyn_sys_lin._fodes_system.steady_solution[0]

        #         m=data[self.m]
        #         l=data[self.l]

        op_point = pi * op_point  # quick workaround - wrong implementation

        force_in_cable = (
            self.m * self.g * (1 - S.One / 2 * (phi - op_point) ** 2)
            + self.m * self.l * (phi - op_point).diff(self.ivar) ** 2
        )
        force_subs = force_in_cable.subs(
            data
        )  # .subs({self.Omega:0.999*dyn_sys_lin.natural_frequencies()[0]})

        return force_subs.doit().expand()

    def sin_coeff(self):

        phi = self.phi
        sin_coeff = (
            self._eoms[0].expand().coeff(sin(phi))
        )  ## eoms[0]: m * L**2 * phi" + m * g * L * sin(phi)  == 0

        return sin_coeff

    def cos_coeff(self):

        phi = self.phi
        cos_coeff = self._eoms[0].expand().coeff(cos(phi))

        return cos_coeff

    def equilibrium_position(self):

        x_e = self.x_e
        phi = self.phi
        equilibrium = (
            self.equilibrium_equation_demo(coordinates=True)
            .doit()
            .subs(x_e, 0)
            .doit()[0]
        )
        solution = solve(equilibrium, phi)
        solution += [i * solution[1] for i in range(2, 10)]

        return solution

    @property
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST,
            mech_comp.PendulumLongitudinalForce,
        ]

        return comp_list


# Kuba - DONE
class LinearizedTriplePendulum(ForcedTriplePendulum):
    scheme_name = "MDOFTriplePendulum.PNG"
    real_name = "TriplePendulum_real.jpg"

    @cached_property
    def components(self):
        components = {}

        components["Pendulum1"] = ForcedTriplePendulum(
            m1=self.m1,
            m2=self.m2,
            m3=self.m3,
            l1=self.l1,
            l2=self.l2,
            l3=self.l3,
            F1=0,
            F2=0,
            F3=0,
            g=self.g,
            phi1=self.phi1,
            phi2=self.phi2,
            phi3=self.phi3,
            phi_u=self.phi_u,
            phi_l=self.phi_l,
            Omega=0,
            ivar=self.ivar,
        ).linearized()

        return components

    def get_default_data(self):

        m0, l0 = symbols("m_0 l_0", positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.l1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l3: [2 * l0, 4 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.phi1: [self.phi_u, 0],
            self.phi2: [self.phi_u, self.phi_l],
            self.phi3: [self.phi_l],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi2] == parameters_dict[self.phi3]:

            parameters_dict[self.phi2] = self.phi_u

        return parameters_dict


# Done Marcelka
class ThreePendulumsWithSprings(ComposedSystem):
    scheme_name = "MDOF_Three_Pendulum_Spring.png"
    real_name = "three_carriages.PNG"

    phi_1 = dynamicsymbols("varphi_1")
    phi_2 = dynamicsymbols("varphi_2")
    phi_3 = dynamicsymbols("varphi_3")
    phi_l = dynamicsymbols("varphi_l")
    phi_r = dynamicsymbols("varphi_r")
    m1 = Symbol("m_1", positive=True)
    m2 = Symbol("m_2", positive=True)
    m3 = Symbol("m_3", positive=True)
    k_1 = Symbol("k_1", positive=True)
    k_2 = Symbol("k_2", positive=True)
    k_3 = Symbol("k_3", positive=True)
    k_4 = Symbol("k_4", positive=True)
    l = Symbol("l", positive=True)
    g = Symbol("g", positive=True)

    qs = dynamicsymbols("varphi_1, varphi_2, varphi_3")

    def __init__(
        self,
        phi_1=None,
        phi_2=None,
        phi_3=None,
        phi_l=None,
        phi_r=None,
        m1=None,
        m2=None,
        m3=None,
        k_1=None,
        k_2=None,
        k_3=None,
        k_4=None,
        l=None,
        g=None,
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if phi_1 is not None:
            self.phi_1 = phi_1
        if phi_2 is not None:
            self.phi_2 = phi_2
        if phi_3 is not None:
            self.phi_3 = phi_3
        if phi_l is not None:
            self.phi_l = phi_l
        if phi_r is not None:
            self.phi_r = phi_r
        if m1 is not None:
            self.m1 = m1
        if m2 is not None:
            self.m2 = m2
        if m3 is not None:
            self.m3 = m3
        if k_1 is not None:
            self.k_1 = k_1
        if k_2 is not None:
            self.k_2 = k_2
        if k_3 is not None:
            self.k_3 = k_3
        if k_4 is not None:
            self.k_4 = k_4
        if l is not None:
            self.l = l
        if g is not None:
            self.g = g
        self.qs = [self.phi_1, self.phi_2, self.phi_3]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self.Pendulum1 = Pendulum(
            self.m1, self.g, self.l, angle=self.phi_1, qs=[self.phi_1]
        )

        self.Pendulum2 = Pendulum(
            self.m2, self.g, self.l, angle=self.phi_2, qs=[self.phi_2]
        )

        self.Pendulum3 = Pendulum(
            self.m3, self.g, self.l, angle=self.phi_3, qs=[self.phi_3]
        )

        self.Spring1 = Spring(
            self.k_1,
            pos1=(
                sqrt(
                    (self.l / 2 * (sin(self.phi_1) - sin(self.phi_2))) ** 2
                    + (self.l / 2 * (cos(self.phi_1) - cos(self.phi_2))) ** 2
                )
            ),
            qs=[self.phi_1, self.phi_2],
        )

        self.Spring2 = Spring(
            self.k_2,
            pos1=(
                sqrt(
                    (self.l / 2 * (sin(self.phi_2) - sin(self.phi_3))) ** 2
                    + (self.l / 2 * (cos(self.phi_2) - cos(self.phi_3))) ** 2
                )
            ),
            qs=[self.phi_2, self.phi_3],
        )

        self.Spring3 = Spring(
            self.k_3,
            pos1=(
                sqrt(
                    (self.l * (sin(self.phi_1) - sin(self.phi_2))) ** 2
                    + (self.l * (cos(self.phi_1) - cos(self.phi_2))) ** 2
                )
            ),
            qs=[self.phi_1, self.phi_2],
        )

        self.Spring4 = Spring(
            self.k_4,
            pos1=(
                sqrt(
                    (self.l * (sin(self.phi_2) - sin(self.phi_3))) ** 2
                    + (self.l * (cos(self.phi_2) - cos(self.phi_3))) ** 2
                )
            ),
            qs=[self.phi_2, self.phi_3],
        )

        components["Pendulum1"] = self.Pendulum1
        components["Pendulum2"] = self.Pendulum2
        components["Pendulum3"] = self.Pendulum3
        components["Spring1"] = self.Spring1
        components["Spring2"] = self.Spring2
        components["Spring3"] = self.Spring3
        components["Spring4"] = self.Spring4

        return components

    def get_default_data(self):

        m0, k0, l0 = symbols("m_0 k_0 l_0", positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.k_1: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.phi_1: [self.phi_l, 0],
            self.phi_2: [self.phi_l, self.phi_r, 0],
            self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if (
            parameters_dict[self.phi_1] != self.phi_l
            or parameters_dict[self.phi_2] != self.phi_l
        ):

            parameters_dict[self.phi_1] = self.phi_l

        if (
            parameters_dict[self.phi_2] != self.phi_r
            or parameters_dict[self.phi_3] != self.phi_r
        ):

            parameters_dict[self.phi_2] = self.phi_r

        return parameters_dict


# DONE - Kuba & Madi
class DampedElasticPendulum(MDoFElasticPendulum):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = "damped_elastic_pendulum.PNG"
    real_name = "elastic_pendulum_real.PNG"

    # class DampedElasticPendulum(MDoFElasticPendulum):
    c = Symbol("c", positive=True)
    k = Symbol("k", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    ivar = Symbol("t")
    z = dynamicsymbols("z")
    phi = dynamicsymbols("\\varphi")

    def __init__(
        self,
        c=None,
        k=None,
        l=None,
        m=None,
        g=None,
        ivar=Symbol("t"),
        z=dynamicsymbols("z"),
        phi=dynamicsymbols("\\varphi"),
        **kwargs
    ):

        if c is not None:
            self.c = c
        if k is not None:
            self.k = k
        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if z is not None:
            self.z = z

        self.qs = [self.phi, self.z]
        self.ivar = ivar
        # self.undamped = undamped_system
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self.pendulum = MDoFElasticPendulum(
            k=self.k,
            l=self.l,
            m=self.m,
            g=self.g,
            z=self.z,
            phi=self.phi,
            ivar=self.ivar,
        )
        self.damper = Damper(c=self.c, pos1=self.z, qs=self.qs)

        components["pendulum"] = self.pendulum
        components["damper"] = self.damper
        # components['mdofelasticpendulum'] = self.mdofelasticpendulum

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r"Spring stiffness",
            self.l: r"Winch length",
            self.m: r"Mass",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict


# DONE - Kuba & Karolina
class DDoFWinch(ComposedSystem):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = "mdof_winch.png"
    real_name = "mdof_winch_real.png"
    I = Symbol("I", positive=True)
    k = Symbol("kappa", positive=True)
    r = Symbol("r", positive=True)
    l = Symbol("l", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    theta = dynamicsymbols("theta")
    phi = dynamicsymbols("varphi")

    def __init__(
        self,
        I=None,
        k=None,
        r=None,
        l=None,
        m=None,
        g=None,
        ivar=Symbol("t"),
        theta=None,
        phi=None,
        **kwargs
    ):

        if I is not None:
            self.I = I
        if k is not None:
            self.k = k
        if r is not None:
            self.r = r
        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if theta is not None:
            self.theta = theta
        if phi is not None:
            self.phi = phi
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        x = self.r * cos(self.phi) + (self.l + self.r * (self.phi - self.theta)) * sin(
            self.phi
        )
        y = self.r * sin(self.phi) + (self.l + self.r * (self.phi - self.theta)) * cos(
            self.phi
        )
        engine_torque = self.m * self.g * self.r

        self.disc_1 = Disk(I=self.I, pos1=self.theta, qs=[self.phi, self.theta])
        self.spring = Spring(self.k, self.theta, qs=[self.phi, self.theta])
        self.material_point_1 = MaterialPoint(self.m, x, qs=[self.phi, self.theta])
        self.material_point_2 = MaterialPoint(self.m, y, qs=[self.phi, self.theta])
        self.M_engine = Force(engine_torque, self.theta, qs=[self.phi, self.theta])
        self.gravity = GravitationalForce(self.m, self.g, pos1=-y, qs=[self.phi])

        components["disc_1"] = self.disc_1
        components["spring"] = self.spring
        components["material_point_1"] = self.material_point_1
        components["material_point_2"] = self.material_point_2
        components["M_engine"] = self.M_engine
        components["gravity"] = self.gravity

        return components

    def get_default_data(self):

        m0, l0, I0, k0, r0 = symbols("m_0 l_0 I_0 k_0 r_0", positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(2, 10)],
            self.l: [S.One * l0 * no for no in range(2, 10)],
            self.I: [2 * I0, 3 * I0, 4 * I0, 5 * I0, 6 * I0],
            self.k: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r"Moment of Inertia of the winch",
            self.k: r"Spring stiffness",
            self.r: r"Winch radius",
            self.l: r"Winch length",
            self.m: r"Mass",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict


class DampedMeasuringTool(ComposedSystem):

    scheme_name = "measure_tool.PNG"
    real_name = "measure_tool_real.PNG"
    # detail_scheme_name =
    # detail_real_name =

    m = Symbol("m", positive=True)
    I = Symbol("I", positive=True)
    l = Symbol("l", positive=True)
    k = Symbol("k", positive=True)
    k_t = Symbol("k_t", positive=True)
    Omega = Symbol("Omega", positive=True)
    F = Symbol("F", positive=True)
    phi = dynamicsymbols("\\varphi")
    c = Symbol("c", positive=True)
    c_t = Symbol("c_t", positive=True)
    lam = Symbol("lambda", positive=True)
    l0 = Symbol("l_0", positive=True)
    lam0 = Symbol("lambda_0", positive=True)

    def __init__(
        self,
        m=None,
        l=None,
        k=None,
        k_t=None,
        ivar=Symbol("t"),
        Omega=None,
        F=None,
        phi=None,
        qs=None,
        c=None,
        c_t=None,
        lam=None,
        l0=None,
        lam0=None,
        **kwargs
    ):
        if l is not None:
            self.l = l
        if m is not None:
            self.m = m
        if k is not None:
            self.k = k
        if k_t is not None:
            self.k_t = k_t
        if F is not None:
            self.F = F
        if Omega is not None:
            self.Omega = Omega
        if phi is not None:
            self.phi = phi
        if c is not None:
            self.c = c
        if c_t is not None:
            self.ct = ct
        if lam is not None:
            self.lam = lam
        if l0 is not None:
            self.l0 = l0
        if lam0 is not None:
            self.lam0 = lam0

        self.qs = [self.phi]
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        # self.bar = Pendulum(self.m, g=0, l=self.l, angle=self.phi, ivar=self.ivar)
        self.bar = MaterialPoint(self.m * self.l**2 / 3, pos1=self.phi, qs=[self.phi])

        self._upper_spring = Spring(self.k, pos1=self.l * self.phi, qs=[self.phi])
        self._lower_spring = Spring(self.k, pos1=self.l * self.phi, qs=[self.phi])
        self._spiral_spring = TorsionalSpring(self.k_t, self.phi, qs=[self.phi])
        self._force = Force(self.F * self.l, pos1=self.phi)
        self._spring_1_damping = Damper(self.c, pos1=self.l * self.phi, qs=[self.phi])
        self._spring_2_damping = Damper(self.c, pos1=self.l * self.phi, qs=[self.phi])
        self._spiral_spring_damping = Damper(self.c_t, pos1=self.phi, qs=[self.phi])

        components["_bar"] = self.bar
        components["_upper_spring"] = self._upper_spring
        components["_lower_spring"] = self._lower_spring
        components["_spiral_spring"] = self._spiral_spring
        components["_force"] = self._force
        components["_spring_1_damping"] = self._spring_1_damping
        components["_spring_2_damping"] = self._spring_2_damping
        components["_spiral_spring_damping"] = self._spiral_spring_damping

        return components

    def get_default_data(self):

        m0, l0, I0, k0, F0, Omg = symbols("m_0 l_0 I_0 k_0 F_0 Omega", positive=True)

        default_data_dict = {
            self.c: [self.lam * (self.k)],
            self.c_t: [self.lam * (self.k_t)],
            self.F: [F0 * cos(Omg * self.ivar)],
            self.m
            * self.l**2
            / 3: [
                S.Half * m0 * self.l**2 / 3,
                m0 * self.l**2 / 3,
                2 * self.m0 * self.l**2 / 3,
                3 * self.m0 * self.l**2 / 3,
                4 * self.m0 * self.l**2 / 3,
            ],
            # self.l : [l0, 2*l0],
        }

        return default_data_dict

    def dynamic_force(self):
        data = self._given_data
        amps = self._fodes_system.steady_solution.as_dict()
        dyn_force = (
            (self.components["_spiral_spring"].force().subs(amps))
            .subs(data)
            .expand()
            .doit()
        )

        return dyn_force

    def static_force(self):
        data = self._given_data
        ans = self.dynamic_force()
        free_coeff = ans.subs(
            {cos(self.Omega * self.ivar): 0, sin(self.Omega * self.ivar): 0}
        ).subs(data)
        return free_coeff

    def steady_state(self):
        return 3 * (S.One / 2 * self.damping_coefficient()) ** (-1)

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force_pin(self):
        lin_sys = self.linearized()

        dyn_comp = (lin_sys.frequency_response_function() * self.l * self.k).subs(
            self._given_data
        )

        total_force = dyn_comp + self.max_static_force_pin()

        return total_force

    def static_force_pin_diameter(self):
        kt = Symbol("k_t", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re)) ** (1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol("k_t", positive=True)
        Re = Symbol("R_e", positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re)) ** (1 / 2)


# issue 589 karolina ziabrowska piotr loch AM
class PendulumWithVaryingMassDamper(ComposedSystem):

    # do zmiany rysunek juz sie wczytuje
    scheme_name = "MDOFDoublePendulum.png"
    real_name = "TriplePendulum_real.jpg"
    detali_scheme_name = "MDOFDoublePendulum.png"
    detali_real_name = "TriplePendulum_real.jpg"
    # do zmiany

    # m=Symbol('m', positive=True)
    ivar = Symbol("t")
    t0 = Symbol("t0")
    m1 = Function("m_alpha")(ivar)
    m2 = Function("m_beta")(ivar)
    l_1 = Symbol("l_alpha", positive=True)
    l_2 = Symbol("l_beta", positive=True)
    g = Symbol("g", positive=True)
    omega = Symbol("Omega", positive=True)
    phi_1 = dynamicsymbols("alpha")
    phi_2 = dynamicsymbols("beta")
    phi_u = dynamicsymbols("varphi_u")
    phi_l = dynamicsymbols("varphi_l")
    u0 = Symbol("u_0", positive=True)
    k = Symbol("k", positive=True)
    b = Symbol("b", positive=True)
    flow_coeff = Symbol("lambda")
    F = Symbol("F", positive=True)
    qs = dynamicsymbols("alpha beta")
    m0 = Symbol("m_0")

    def __init__(
        self,
        # m=None,
        m1=None,
        m2=None,
        l_1=None,
        l_2=None,
        g=None,
        omega=None,
        phi_1=None,
        phi_2=None,
        phi_u=None,
        phi_l=None,
        u0=None,
        k=None,
        b=None,
        F=None,
        qs=None,
        ivar=Symbol("t"),
        t0=t0,
        flow_coeff=flow_coeff,
        m0=None,
        **kwargs
    ):

        if m1 is not None:
            self.m1 = m1
        # if m is not None: self.m=m
        if m2 is not None:
            self.m2 = m2
        if l_1 is not None:
            self.l_1 = l_1
        if l_2 is not None:
            self.l_2 = l_2
        if g is not None:
            self.g = g
        if omega is not None:
            self.omega = omega
        if phi_1 is not None:
            self.phi_1 = phi_1
        if phi_2 is not None:
            self.phi_2 = phi_2
        if phi_u is not None:
            self.phi_u = phi_u
        if phi_l is not None:
            self.phi_l = phi_l
        if u0 is not None:
            self.u0 = u0
        if k is not None:
            self.k = k
        if b is not None:
            self.b = b
        if flow_coeff is not None:
            self.flow_coeff = flow_coeff
        if F is not None:
            self.F = F
        if t0 is not None:
            self.t0 = t0
        if m0 is not None:
            self.m0 = m0

        self.qs = [self.phi_1, self.phi_2]
        self._init_from_components(**kwargs)

        from sympy import atan, pi

        #         display(type(ivar),type(t0))
        self.trans_expr = S.One / 2 - atan(flow_coeff * (ivar - t0)) / pi

    @cached_property
    def components(self):
        components = {}

        self.x_2 = sin(self.phi_1) * self.l_1 + sin(self.phi_2) * self.l_2
        self.y_2 = cos(self.phi_1) * self.l_1 + cos(self.phi_2) * self.l_2

        self.Pendulum1 = Pendulum(
            self.m1, self.g, self.l_1, angle=self.phi_1, qs=[self.phi_1]
        )
        self.material_point_11 = MaterialPoint(
            self.m2, self.x_2, qs=[self.phi_1, self.phi_2]
        )
        self.material_point_21 = MaterialPoint(
            self.m2, self.y_2, qs=[self.phi_1, self.phi_2]
        )
        self.gravity_1 = GravitationalForce(
            self.m2, self.g, pos1=-self.y_2, qs=[self.phi_2]
        )
        self.damper_1 = Damper(self.b, self.phi_1, qs=[self.phi_1])
        self.damper_2 = Damper(self.b, self.phi_2, qs=[self.phi_2])
        self.force = Force(
            self.k * self.u0 * sin(self.omega * self.ivar) * self.l_1,
            self.phi_1,
            qs=[self.phi_1],
        )

        components["Pendulum 1"] = self.Pendulum1
        components["material_point_11"] = self.material_point_11
        components["material_point_21"] = self.material_point_21
        components["gravity_1"] = self.gravity_1
        components["damper_1"] = self.damper_1
        components["damper_2"] = self.damper_2
        components["force"] = self.force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            #             self.m: r'mass of the system',
            self.m1: r"initial mass of the upper pendulum member",
            self.m2: r"initial mass of the lower pendulum member",
            #             self.m_f1: r'mass of the upper pendulum',
            #             self.m_f2: r'mass of the lower pendulum',
            self.m1.diff(
                self.ivar
            ): r"mass of the upper pendulum and fl and fluid flow",
            self.m2.diff(
                self.ivar
            ): r"mass of the of the lower pendulum and fluid flow",
            #             self.m_t: r'total mass of the system',
            #             self.m_t2: r'total mass of the lower pendulum',
            self.m0: r"transferred damping mass",
            self.l_1: r"length of the upper pendulum",
            self.l_2: r"length of the lower pendulum",
            #             self.l_3: r'length of the damper pendulum',
            #             self.w: r'angular velocity',
            self.k: r"spring stiffness",
            #             self.ls: r'length of the spring',
            # self.u: r'difference of spring length and initial spring length',
            #             self.x0: r'initial position of the system',
            #             self.l0: r'initial spring length',
            self.u0: r"initial difference of spring length and starting spring length",
            #             self.x_k: r'forcing force',
            self.b: r"damping coefficient value",
            self.omega: r"excitation frequency",
            # self.f: r'????',
            self.g: r"acceleration of gravity",
            #             self.rho_f: r'fluid density',
            self.phi_1: r"angle of the upper pendulum swing",
            self.phi_2: r"angle of the lower pendulum swing",
            self.phi_1.diff(
                self.ivar
            ): r"angular (generalised) velocity of the upper pendulum member",
            self.phi_2.diff(
                self.ivar
            ): r"angular (generalised) velocity of the lower pendulum member",
            self.phi_1.diff(
                self.ivar, self.ivar
            ): r"angular (generalised) acceleration of the upper pendulum member",
            self.phi_2.diff(
                self.ivar, self.ivar
            ): r"angular (generalised) acceleration of the lower pendulum member",
            self.t0: r"activation time",
            self.flow_coeff: r"fluid flow rate",
            self.ivar: r"time",
            #             self.h_fa:r'centre of gravity of the higher member',
            #             self.h_fb:r'centre of gravity of the lower member',
            #             self.A_pipe1:r'cross-section area of the upper member',
            #             self.A_pipe2:r'cross-section area of the lower member',
        }
        return self.sym_desc_dict

    def flow_equation(self):
        return_eq = Eq(self.m2, (self.m0 - self.m0 * self.trans_expr).expand())

        return return_eq

    def upper_mass_equation(self):
        return_eq = Eq(self.m1, (self.m0 * self.trans_expr).expand())

        return return_eq

    def amplitude_envelope(self, signal, var, frame_size=128, hop_length=56):

        max_vals = [
            max(signal.iloc[i : i + frame_size][var])
            for i in range(0, len(signal), hop_length)
        ]

        envelope_averaged = pd.DataFrame(max_vals).rolling(10, min_periods=2).mean()

        return envelope_averaged

    def instantaneous_frequency(
        self, signal, var, t_param, frame_size=128, hop_length=56
    ):

        time_id = [
            signal.iloc[i : i + frame_size][var].idxmax()
            for i in range(0, len(signal), hop_length)
        ]

        insta_freq = pd.DataFrame(time_id).apply(lambda row: t_param * row)

        return insta_freq


class DampedPendulumWithVaryingStiffness(ComposedSystem):

    # do zmiany rysunek juz sie wczytuje
    scheme_name = "MDOFDoublePendulum.png"
    real_name = "TriplePendulum_real.jpg"
    detali_scheme_name = "MDOFDoublePendulum.png"
    detali_real_name = "TriplePendulum_real.jpg"
    # do zmiany

    # m=Symbol('m', positive=True)
    ivar = Symbol("t")

    m1 = Symbol("m_alpha")  # (ivar)
    m2 = Symbol("m_beta")  # (ivar)
    l_1 = Symbol("l_alpha", positive=True)
    l_2 = Symbol("l_beta", positive=True)
    g = Symbol("g", positive=True)
    omega = Symbol("Omega", positive=True)
    phi_1 = dynamicsymbols("alpha")
    phi_2 = dynamicsymbols("beta")
    phi_u = dynamicsymbols("varphi_u")
    phi_l = dynamicsymbols("varphi_l")
    F = Symbol("F", positive=True)
    k = Function("k")(ivar)
    b = Symbol("b", positive=True)
    tau0 = Symbol("tau_0", positive=True)
    flow_coeff = Symbol("lambda")
    qs = dynamicsymbols("alpha beta")

    def __init__(
        self,
        # m=None,
        m1=None,
        m2=None,
        l_1=None,
        l_2=None,
        g=None,
        omega=None,
        phi_1=None,
        phi_2=None,
        phi_u=None,
        phi_l=None,
        u0=None,
        k=None,
        b=None,
        tau0=None,
        flow_coeff=None,
        qs=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if m1 is not None:
            self.m1 = m1
        # if m is not None: self.m=m
        if m2 is not None:
            self.m2 = m2
        if l_1 is not None:
            self.l_1 = l_1
        if l_2 is not None:
            self.l_2 = l_2
        if g is not None:
            self.g = g
        if omega is not None:
            self.omega = omega
        if phi_1 is not None:
            self.phi_1 = phi_1
        if phi_2 is not None:
            self.phi_2 = phi_2
        if phi_u is not None:
            self.phi_u = phi_u
        if phi_l is not None:
            self.phi_l = phi_l
        if u0 is not None:
            self.u0 = u0
        if k is not None:
            self.k = k
        if b is not None:
            self.b = b
        if flow_coeff is not None:
            self.flow_coeff = flow_coeff

        self.qs = [self.phi_1, self.phi_2]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self.x_2 = sin(self.phi_1) * self.l_1 + sin(self.phi_2) * self.l_2
        self.y_2 = cos(self.phi_1) * self.l_1 + cos(self.phi_2) * self.l_2

        self.Pendulum1 = Pendulum(
            self.m1, self.g, self.l_1, angle=self.phi_1, qs=[self.phi_1]
        )
        self.material_point_11 = MaterialPoint(
            self.m2, self.x_2, qs=[self.phi_1, self.phi_2]
        )
        self.material_point_21 = MaterialPoint(
            self.m2, self.y_2, qs=[self.phi_1, self.phi_2]
        )
        self.gravity_1 = GravitationalForce(
            self.m2, self.g, pos1=-self.y_2, qs=[self.phi_2]
        )
        self.damper_1 = Damper(self.b, self.phi_1, qs=[self.phi_1])
        self.damper_2 = Damper(self.b, self.phi_2, qs=[self.phi_2])
        # self.force = Force( (1-exp(-self.ivar/self.tau0))*self.F*sin(self.omega*self.ivar)*self.l_1,self.phi_1,qs=[self.phi_1])

        self.force_thrust = Force(
            self.F * (1 + sin(self.omega * self.ivar)),
            self.x_2,
            qs=[self.phi_1, self.phi_2],
        )
        self.force_lateral = Force(
            self.F * cos(self.omega * self.ivar), self.y_2, qs=[self.phi_1, self.phi_2]
        )

        self.spring = Spring(
            self.k, self.phi_1, self.phi_2, qs=[self.phi_1, self.phi_2]
        )
        self.damper_c = Damper(
            self.b, self.phi_1, self.phi_2, qs=[self.phi_1, self.phi_2]
        )

        components["Pendulum 1"] = self.Pendulum1
        components["material_point_11"] = self.material_point_11
        components["material_point_21"] = self.material_point_21
        components["gravity_1"] = self.gravity_1
        components["damper_1"] = self.damper_1
        components["damper_2"] = self.damper_2
        components["force_thrust"] = self.force_thrust
        components["force_lateral"] = self.force_lateral
        components["spring"] = self.spring
        components["damper_c"] = self.damper_c

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"mass of the system",
            self.m1: r"initial mass of the upper pendulum member",
            self.m2: r"initial mass of the lower pendulum member",
            self.m_f1: r"mass of the upper pendulum",
            self.m_f2: r"mass of the lower pendulum",
            self.m_f1.diff(t): r"mass fluid flow of the upper pendulum",
            self.m_f2.diff(t): r"mass fluid flow of the lower pendulum",
            self.m_t: r"total mass of the system",
            self.m_t2: r"total mass of the lower pendulum",
            self.m_0: r"transferred damping mass",
            self.l_1: r"length of the upper pendulum",
            self.l_2: r"length of the lower pendulum",
            self.l_3: r"length of the damper pendulum",
            self.w: r"angular velocity",
            self.k: r"spring stiffness",
            self.ls: r"length of the spring",
            # self.u: r'difference of spring length and initial spring length',
            self.x0: r"initial position of the system",
            self.l0: r"initial spring length",
            self.u0: r"initial difference of spring length and starting spring length",
            self.x_k: r"forcing force",
            self.b: r"damping coefficient value",
            self.omega: r"excitation frequency",
            # self.f: r'????',
            self.g: r"acceleration of gravity",
            self.rho_f: r"fluid density",
            self.phi_1: r"angle of the upper pendulum swing",
            self.phi_2: r"angle of the lower pendulum swing",
            self.phi_1.diff(t): r"generalized velocity of the upper pendulum member",
            self.phi_2.diff(t): r"generalized velocity of the lower pendulum member",
            self.phi_1.diff(
                t, t
            ): r"generalized acceleration of the upper pendulum member",
            self.phi_2.diff(
                t, t
            ): r"generalized acceleration of the lower pendulum member",
            self.t0: r"activation time",
            self.flow_coeff: r"fluid flow rate",
            self.ivar: r"time",
            self.h_fa: r"centre of gravity of the higher member",
            self.h_fb: r"centre of gravity of the lower member",
            self.A_pipe1: r"cross-section area of the upper member",
            self.A_pipe2: r"cross-section area of the lower member",
        }
        return self.sym_desc_dict


### Classes of RollingBar
class RollingBar(ComposedSystem):
    """
    Model of a Single Degree of Freedom - RollingBar

        Arguments:
        =========
            m = Mass
                -Mass of the Bar

            l = length
                -Length of th Bar

            h = thickness
                -thickness of the Bar

            g = gravitional field
                -value of gravitional field acceleration

            r = radius
                -radius of the Half Disc

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -angle of the Bar

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A rolling bar on half disc with radius R.

        >>> t = symbols('t')
        >>> R = symbols(' R ',positive=True)
        >>> RollingBar(r=R)

        -We define the symbols and dynamicsymbols
        -determine the instance of the RollingBar by using class RollingBar()
    """

    scheme_name = "rolling_bar.png"
    real_name = "rolling_bar.png"
    r = Symbol("r", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    h = Symbol("h", positive=True)
    phi = dynamicsymbols("varphi")

    def __init__(
        self,
        r=None,
        m=None,
        l=None,
        h=None,
        g=None,
        ivar=Symbol("t"),
        phi=None,
        **kwargs
    ):

        if r is not None:
            self.r = r
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if h is not None:
            self.h = h

        if phi is not None:
            self.phi = phi
        self.ivar = ivar
        self.qs = [self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        #         self.x = (self.r+self.h/2) / cos(self.phi)
        #         self.y =  (self.r+self.h/2) * cos(self.phi) - self.r * self.phi * sin(self.phi)
        self.x = (self.r + self.h / 2) * sin(self.phi) - (
            self.r * self.phi * cos(self.phi)
        )
        self.y = (self.r + self.h / 2) * cos(self.phi) + self.r * self.phi * sin(
            self.phi
        )
        I = self.m / 12 * (self.l**2 + self.h**2)

        self.material_point_1 = MaterialPoint(
            self.m, self.x, ivar=self.ivar, qs=[self.phi]
        )
        self.material_point_2 = MaterialPoint(
            self.m, self.y, ivar=self.ivar, qs=[self.phi]
        )
        self.material_point_3 = MaterialPoint(
            I, self.phi, ivar=self.ivar, qs=[self.phi]
        )
        self.gravity = GravitationalForce(
            self.m, self.g, pos1=-self.y, ivar=self.ivar, qs=[self.phi]
        )

        components["material_point_1"] = self.material_point_1
        components["material_point_2"] = self.material_point_2
        components["material_point_3"] = self.material_point_3
        components["gravity"] = self.gravity

        return components

    def get_default_data(self):

        m0, r0, l0, h0 = symbols("m_0 r_0 l_0 h_0", positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(2, 7)],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
            self.l: [1.4 * l0, 1.5 * l0, 1.6 * l0, 1.7 * l0, 1.8 * l0],
            self.h: [S.One * h0 * no for no in range(1, 7)],
        }
        return default_data_dict

    def get_numerical_data(self):
        m0, l0, r0 = symbols("m_0 l_0 r_0", positive=True)

        default_data_dict = {
            self.m: [1],
            self.r: [1],
            self.h: [0.05],
            self.l: [1.8],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r"Half disc radius",
            self.m: r"Mass",
            self.l: r"Length of the Bar",
            self.h: r"Thickness of the Bar",
            self.g: "Gravity constant",
        }
        return self.sym_desc_dict

    def unit_dict(self):

        unit_dict = {
            self.r: ureg.meter,
            self.l: ureg.meter,
            self.h: ureg.meter,
            self.m: ureg.kilogram,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.ivar: ureg.second,
            self.phi: ureg.radian,
            self.phi.diff(self.ivar): ureg.radian / ureg.second,
            self.phi.diff(self.ivar, 2): ureg.radian / ureg.second / ureg.second,
        }
        return unit_dict


class ForcedRollingBar(RollingBar):
    """
    Model of a Single Degree of Freedom - ForcedRollingBar

        Arguments:
        =========
            m = Mass
                -Mass of the Bar

            l = length
                -Length of th Bar

            h = thickness
                -thickness of the Bar

            g = gravitional field
                -value of gravitional field acceleration

            r = radius
                -radius of the Half Disc

            M_0 = torque amplitude
                -amplitude value of horisontal force

            Omega = Frequency of force
                -value of frequency force

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -angle of the Bar

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A rolling bar on half disc with radius R.

        >>> t = symbols('t')
        >>> R = symbols(' R ',positive=True)
        >>> ForcedRollingBar(r=R)

        -We define the symbols and dynamicsymbols
        -determine the instance of the ForcedRollingBar by using class ForcedRollingBar()
    """

    scheme_name = "forced_rolling_bar.png"
    real_name = "forced_rolling_bar.png"
    r = Symbol("r", positive=True)
    m = Symbol("m", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    h = Symbol("h", positive=True)
    M_0 = Symbol("M_0", positive=True)
    Omega = Symbol("Omega", positive=True)
    phi = dynamicsymbols("varphi")

    def __init__(
        self,
        r=None,
        m=None,
        l=None,
        h=None,
        g=None,
        M_0=None,
        Omega=None,
        ivar=Symbol("t"),
        phi=None,
        **kwargs
    ):

        if r is not None:
            self.r = r
        if m is not None:
            self.m = m
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if h is not None:
            self.h = h
        if M_0 is not None:
            self.M_0 = M_0
        if Omega is not None:
            self.Omega = Omega

        if phi is not None:
            self.phi = phi
        self.ivar = ivar
        self.qs = [self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        #         self.x = (self.r+self.h/2) / cos(self.phi)
        #         self.y =  (self.r+self.h/2) * cos(self.phi) - self.r * self.phi * sin(self.phi)
        self.x = (self.r + self.h / 2) * sin(self.phi) - (
            self.r * self.phi * cos(self.phi)
        )
        self.y = (self.r + self.h / 2) * cos(self.phi) + self.r * self.phi * sin(
            self.phi
        )
        I = self.m / 12 * (self.l**2 + self.h**2)

        self.roling_bar = RollingBar(
            r=self.r, m=self.m, l=self.l, h=self.h, g=self.g, phi=self.phi
        )
        self.torque = Force(
            self.M_0 * cos(self.Omega * self.ivar), self.phi, qs=self.qs
        )

        components["roling_bar"] = self.roling_bar
        components["torque"] = self.torque

        return components

    def get_default_data(self):

        m0, r0, l0, h0 = symbols("m_0 r_0 l_0 h_0", positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(2, 7)],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
            self.l: [1.4 * l0, 1.5 * l0, 1.6 * l0, 1.7 * l0, 1.8 * l0],
            self.h: [S.One * h0 * no for no in range(1, 7)],
        }
        return default_data_dict

    def get_numerical_data(self):
        m0, l0, r0 = symbols("m_0 l_0 r_0", positive=True)

        default_data_dict = {
            self.m: [1],
            self.r: [1],
            self.h: [0.05],
            self.l: [1.8],
            self.M_0: [1],
            self.Omega: [50],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r"Half disc radius",
            self.m: r"Mass",
            self.l: r"Length of the Bar",
            self.h: r"Thickness of the Bar",
            self.g: r"Gravity constant",
            self.M_0: r"Torque amplitude",
            self.Omega: r"Angular frequency of torque",
        }
        return self.sym_desc_dict

    def unit_dict(self):

        unit_dict = {
            self.r: ureg.meter,
            self.l: ureg.meter,
            self.h: ureg.meter,
            self.m: ureg.kilogram,
            self.g: ureg.meter / ureg.second / ureg.second,
            self.ivar: ureg.second,
            self.phi: ureg.radian,
            self.phi.diff(self.ivar): ureg.radian / ureg.second,
            self.phi.diff(self.ivar, 2): ureg.radian / ureg.second / ureg.second,
            self.M_0: ureg.newton * ureg.meter,
            self.Omega: ureg.radian / ureg.second,
        }
        return unit_dict


class CompoundPendulum(ComposedSystem):
    """
    The CompoundPendulum class models the motion of a compound pendulum by describing its dynamics based on fundamental physical parameters such as mass, moment of inertia, length, gravity, and deflection angle. These parameters are then used to calculate the position of the mass in space, as well as to model the forces acting on the system.


        Arguments:
        =========
            m = mass
                -The mass of the entire pendulum system

            I = moment of inertia
                -The moment of inertia associated with the rotation of the pendulum around its axis of rotation

            g = gravitional field
                -value of gravitional's field acceleration

            l = lenght
                -Distance of the center of mass of a Compound pendulum from the axis of rotation

            phi = angle
                -A dynamic variable that represents the angle of rotation of the pendulum as a function of time

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Represents time as an independent variable, upon which all dynamic changes in the system are based

        Example
        =======
        A compound pendulum swinging with a rotation angle `phi` around its pivot.

        >>> t = symbols('t') # Time variable
        >>> m, I, g, l = symbols('m, I, g, l') # Physical parameters
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates (angle)
        >>> CompoundPendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the Compound Pendulum by using class Pendulum()
    """

    scheme_name = "compound_pendulum.png"
    real_name = "compound_pendulum.png"

    m = Symbol("m", positive=True)
    I = Symbol("I", positive=True)
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    phi = dynamicsymbols("\phi")
    qs = None

    m0 = Symbol("m_0", positive=True)
    l0 = Symbol("l_0", positive=True)

    def __init__(
        self, m=None, I=None, g=None, l=None, phi=None, ivar=Symbol("t"), **kwargs
    ):

        if m is not None:
            self.m = m
        if I is not None:
            self.I = I
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if phi is not None:
            self.phi = phi
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        """
        System Components:
        =========
            x = The position of the pendulum's mass along the horizontal axis, dependent on the angle phi
            y = The position of the pendulum's mass along the vertical axis, dependent on the angle phi

            material_point_1 = Material point associated with motion along the x-axis
            material_point_2 = Material point associated with motion along the y-axis
            material_point_3 = Material point associated with rotational motion
            gravity = Gravitational force
        """
        components = {}
        self.x = self.l * sin(self.phi)
        self.y = self.l * (1 - cos(self.phi))

        self.material_point_1 = MaterialPoint(
            self.m, self.x, ivar=self.ivar, qs=[self.phi]
        )
        self.material_point_2 = MaterialPoint(
            self.m, self.y, ivar=self.ivar, qs=[self.phi]
        )
        self.material_point_3 = MaterialPoint(
            self.I, self.phi, ivar=self.ivar, qs=[self.phi]
        )
        self.gravity = GravitationalForce(
            self.m, self.g, pos1=self.y, ivar=self.ivar, qs=[self.phi]
        )

        components["material_point_1"] = self.material_point_1
        components["material_point_2"] = self.material_point_2
        components["material_point_3"] = self.material_point_3
        components["gravity"] = self.gravity

        return components

    def get_default_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            self.I: [S.One * no for no in range(3, 9)],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.phi: r"angle of the pendulum",
            self.phi.diff(self.ivar): r"velocity of the pendulum",
            self.phi.diff(self.ivar, 2): r"acceleration of the pendulum",
            self.ivar: r"time",
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [*REPORT_COMPONENTS_LIST]

        return comp_list

    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,  # Mass of the pendulum
            self.I: ureg.kilogram * ureg.meter**2,  # Moment of inertia
            self.g: ureg.meter / ureg.second**2,  # Gravitational acceleration
            self.l: ureg.meter,  # Length of the pendulum
            self.phi: ureg.radian,  # Angular displacement (radians)
            self.phi.diff(self.ivar): ureg.radian
            / ureg.second,  # Angular velocity (rad/s)
            self.phi.diff(self.ivar, 2): ureg.radian
            / ureg.second**2,  # Angular acceleration (rad/s^2)
            self.ivar: ureg.second,  # Time (seconds)
        }

        return units_dict


class ForcedCompoundPendulum(ComposedSystem):
    """
    Model of a Forced Compound Pendulum with an external periodic force applied horizontally.


            Arguments:
            =========
                m = mass
                    -The mass of the entire pendulum system

                I = moment of inertia
                    -The moment of inertia associated with the rotation of the pendulum around its axis of rotation

                g = gravitional field
                    -value of gravitional's field acceleration

                l = lenght
                    -Distance of the center of mass of a Compound pendulum from the axis of rotation

                phi = angle
                    -A dynamic variable that represents the angle of rotation of the pendulum as a function of time

                M_0 = force amplitude
                    -Amplitude of the horizontal force applied to the pendulum

                Omega = Frequency of force
                    -Frequency of the applied horizontal force.

                ivar = symbol object
                    -Independant time variable

                qs = dynamicsymbol object
                    -Represents time as an independent variable, upon which all dynamic changes in the system are based

            Example
            =======
            A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

            >>> t = symbols('t')  # Time variable
            >>> m, I, g, l, M_0, Omega = symbols('m I g l M_0 Omega')  # Physical parameters

            >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
            >>> ForcedCompoundPendulum()

            -We define the symbols and dynamicsymbols
            -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
            -determine the instance of the Compound Pendulum by using class Pendulum()
    """

    scheme_name = "forced_compound_pendulum.png"
    real_name = "forced_compound_pendulum.png"

    m = Symbol("m", positive=True)
    I = Symbol("I", positive=True)
    g = Symbol("g", positive=True)
    Omega = Symbol("Omega", positive=True)
    M_0 = Symbol("M_0", positive=True)
    l = Symbol("l", positive=True)
    phi = dynamicsymbols("\phi")
    qs = None
    m0 = Symbol("m_0", positive=True)
    l0 = Symbol("l_0", positive=True)

    def __init__(
        self,
        m=None,
        I=None,
        g=None,
        l=None,
        M_0=None,
        Omega=None,
        phi=None,
        ivar=None,
        **kwargs
    ):

        if m is not None:
            self.m = m
        if I is not None:
            self.I = I
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if M_0 is not None:
            self.M_0 = M_0
        if Omega is not None:
            self.Omega = Omega
        if phi is not None:
            self.phi = phi
        if ivar is not None:
            self.ivar = ivar

        self.qs = [self.phi]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        self.compound_pendulum = CompoundPendulum(
            m=self.m, l=self.l, phi=self.phi, I=self.I, g=self.g, qs=self.qs
        )
        self.force = Force(self.M_0 * cos(self.Omega * self.ivar), self.phi, qs=self.qs)
        components["compound_pendulum"] = self.compound_pendulum
        components["force"] = self.force

        return components

    def get_default_data(self):

        m0, l0 = self.m0, self.l0

        default_data_dict = {
            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            self.m * self.l**2: [self.m * self.l**2],
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            self.I: [S.One * no for no in range(3, 9)],
            self.M_0: [S.One * no for no in range(3, 9)],
            self.Omega: [S.One * no for no in range(60, 63)],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r"Mass of pendulum",
            self.g: r"Gravity constant",
            self.l: r"Pendulum length",
            self.phi: r"angle of the pendulum",
            self.phi.diff(self.ivar): r"velocity of the pendulum",
            self.phi.diff(self.ivar, 2): r"acceleration of the pendulum",
            self.ivar: r"time",
            self.M_0: r"Torque amplitude",
            self.Omega: r"Angular frequency of torque",
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [*REPORT_COMPONENTS_LIST]

        return comp_list

    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,  # Mass of the pendulum
            self.I: ureg.kilogram * ureg.meter**2,  # Moment of inertia
            self.g: ureg.meter / ureg.second**2,  # Gravitational acceleration
            self.l: ureg.meter,  # Length of the pendulum
            self.M_0: ureg.newton
            * ureg.meter,  # Moment of force - torque (newton-metre)
            self.Omega: ureg.radian / ureg.second,  # Frequency of the force (rad/s)
            self.phi: ureg.radian,  # Angular displacement (radians)
            self.phi.diff(self.ivar): ureg.radian
            / ureg.second,  # Angular velocity (rad/s)
            self.phi.diff(self.ivar, 2): ureg.radian
            / ureg.second**2,  # Angular acceleration (rad/s^2)
            self.ivar: ureg.second,  # Time (seconds)
        }

        return units_dict
