from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from ..continuous import ContinuousSystem, PlaneStressProblem

from .mechanics.pendulum import Pendulum, FreePendulum, PulledPendulum, ExcitedPendulum, DampedPendulum, ExcitedDampedPendulum, PendulumKinematicExct
from .mechanics.engine import EngineVerticalSpringGravity

import base64
import random
import IPython as IP
import numpy as np
import inspect


class ComposedSystem(HarmonicOscillator):
    """Base class for all systems

    """
    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'
    _default_args = ()
    _default_folder_path = "./dynpy/models/images/"

    z = dynamicsymbols('z')

    m0 = Symbol('m_0', positive=True)
    k0 = Symbol('k_0', positive=True)
    F0 = Symbol('F_0', positive=True)
    Omega0 = Symbol('Omega_0', positive=True)

    @classmethod
    def _scheme(cls):

        path = cls._default_folder_path + cls.scheme_name

        return path

    @classmethod
    def _real_example(cls):
        path = cls._default_folder_path + cls.real_name

        return path

    @classmethod
    def _detail_real(cls):
        path = cls._default_folder_path + cls.detail_real_name

        return path

    @classmethod
    def _detail_scheme(cls):
        path = cls._default_folder_path + cls.detail_scheme_name

        return path

    def _init_from_components(self, *args, system=None, **kwargs):

        if system is None:
            composed_system = self._elements_sum
        else:
            composed_system = system

        #print('CS',composed_system._components)
        super().__init__(None, system=composed_system)

        #print('self',self._components)
        if self._components is None:
            comps = {}
        else:
            comps = self._components

        self._components = {**comps, **self.components}

    def __init__(self,
                 Lagrangian=None,
                 m0=None,
                 qs=None,
                 forcelist=None,
                 bodies=None,
                 frame=None,
                 hol_coneqs=None,
                 nonhol_coneqs=None,
                 label=None,
                 ivar=None,
                 evaluate=True,
                 system=None):

        if ivar is not None: self.ivar = ivar
        if m0 is not None: self.m0 = m0

        if qs is not None:
            self.qs = qs
        else:
            self.qs = [self.z]

        self._init_from_components(system=system)

    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self.m0, self.qs[0],
                                             self.qs)('Material Point')
        components['_material_point'] = self._material_point

        return components

    @property
    def elements(self):

        return {**super().components, **self.components}

    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        elif example == 'detail_scheme_name':
            path = cls._detail_scheme()
        elif example == 'detail_real_name':
            path = cls._detail_real()
        else:
            path = cls._scheme()
        print(path)
        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))

    def get_default_data(self):
        return None

    def get_numerical_data(self):
        return None

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    def get_numerical_parameters(self):

        default_data_dict = self.get_numerical_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            mech_comp.GoverningEquationComponent,
            mech_comp.FundamentalMatrixComponent,
            mech_comp.GeneralSolutionComponent,
            mech_comp.SteadySolutionComponent,
        ]

        return comp_list

    def linearized(self):

        return type(self).from_system(super().linearized())

    def tensioner_belt_force(self):
        return self.k_tensioner * self.steady_solution()

    def left_belt_force(self):
        return self.k_belt * self.steady_solution()

    def right_belt_force(self):
        return self.k_belt * self.steady_solution()

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force_pin(self):
        return self.frequency_response_function() * self.stiffness_matrix(
        )[0] + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


class CompoundSystem(ComposedSystem):

    z = dynamicsymbols('z')
    _p = Symbol('p')

    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self._p, self.qs[0],
                                             self.qs)('Material Point')
        components['_material_point'] = self._material_point

        return components


class NonlinearComposedSystem(ComposedSystem):

    def frequency_response_function(self,
                                    frequency=Symbol('Omega', positive=True),
                                    amplitude=Symbol('a', positive=True)):

        omega = ComposedSystem(self.linearized()).natural_frequencies()[0]
        eps = self.small_parameter()

        exciting_force = self.external_forces()[0]

        comps = exciting_force.atoms(sin, cos)
        exciting_amp = sum([exciting_force.coeff(comp) for comp in comps])
        inertia = self.inertia_matrix()[0]

        return amplitude * (-frequency**2 + omega**2) * inertia + S(
            3) / 4 * eps * amplitude**3 - exciting_amp

    def amplitude_from_frf(self, amplitude=Symbol('a', positive=True)):

        return solveset(self.frequency_response_function(), amplitude)

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            mech_comp.LinearizationComponent,
            mech_comp.GoverningEquationComponent,
            mech_comp.FundamentalMatrixComponent,
            mech_comp.GeneralSolutionComponent,
            mech_comp.SteadySolutionComponent,
        ]

        return comp_list


#pati
class MaterialPointMovement(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('\\varphi')

    c0 = Symbol('C0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar

        self.qs = [self.phi]

        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity_ = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)

        composed_system = self._mass_x + self._mass_y + self._gravity_

        super().__init__(composed_system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [C0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict

    def max_static_force(self):
        return S.Zero

    def max_dynamic_force(self):
        return S.Zero


from .mechanics.tensioner import BlowerToothedBelt, DampedBlowerToothedBelt



from .mechanics.trolley import SpringMassSystem, SpringDamperMassSystem

from .mechanics.bridge import BeamBridge


#Amadi
class BeamBridgeDamped(ComposedSystem):
    """Ready to use model of damped bridge represented by the mass supported by elastic beam.
        Arguments:
        =========
            m = Symbol object
                -Mass embedded on beam.

            k = Symbol object
                -Bending stiffness of the beam

            g = Symbol object
                -Gravitational field acceleration

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass hanged on the elastic beam with the stiffness k in the gravitational field

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
    scheme_name = 'bridge_dmp.png'
    real_name = 'beam_bridge_real.PNG'

    m = Symbol('m', positive=True)
    k_beam = Symbol('k_beam', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    l = Symbol('l', positive=True)
    module = Symbol('E', positive=True)
    inertia = Symbol('I', positive=True)
    z = dynamicsymbols('z')
    c = Symbol('c', positive=True)
    E0 = Symbol('E_0', positive=True)
    I0 = Symbol('I_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    m0 = Symbol('m_0', positive=True)
    lam0 = Symbol('lambda_0', positive=True)
    F0 = Symbol('F_0', positive=True)
    lam = Symbol('lambda', positive=True)

    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 l=None,
                 module=None,
                 inertia=None,
                 z=None,
                 c=None,
                 lam=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if l is not None: self.l = l
        if z is not None: self.z = z
        if module is not None: self.module = module
        if inertia is not None: self.inertia = inertia
        if c is not None: self.c = c
        if lam is not None: self.lam = lam
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z,
                                   qs=[self.z])(label='Material point')
        self._spring = Spring(self.k_beam, self.z,
                              qs=[self.z])(label='Beam stiffness')
        self._gravity_force = GravitationalForce(self.m, self.g,
                                                 self.z)(label='Gravity field')
        self._force = Force(-self.F * sin(self.Omega * self.ivar),
                            pos1=self.z)(label='External force')
        self._damper = Damper(self.c, pos1=self.z,
                              qs=[self.z])(label='Damping of a beam')

        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravity_force'] = self._gravity_force
        components['_force'] = self._force
        components['_damper'] = self._damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'Gravitational field acceleration',
            self.Omega: r'Excitation frequency',
            self.F: r'Force acting on a bridge',
            self.l: r'Lenght of a beam',
            self.module: r'Youngs modulus',
            self.inertia: r'Interia of a beam',
            self.c: r'Damping coefficient'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        #E0, I0, l0, m0, k0,c0, lam0= symbols('E_0 I_0 l_0 m_0 k_0 c_0 lambda_0', positive=True)
        E0, I0, l0, m0, lam0, F0 = self.E0, self.I0, self.l0, self.m0, self.lam0, self.F0
        default_data_dict = {

            #self.lam:[10,
            self.m: [S.One * no * m0 * 10 for no in range(2, 8)],
            self.module: [S.One * E0 * no for no in range(10, 20)],
            self.inertia: [S.One * no * I0 for no in range(10, 20)],
            self.l: [S.One * no * l0 for no in range(10, 25)],
            self.lam: [S.One * no * lam0 for no in range(1, 5)],
            self.F: [S.One * no * F0 for no in range(10, 25)],
            self.k_beam: [S.One * 48 * self.module * self.inertia / self.l**3],
            self.c: [self.k_beam * self.lam]
        }
        return default_data_dict


class SDOFSemiSubmersiblePlatform(ComposedSystem):

    #scheme_name = '.PNG'
    #real_name = '.jpg'

    z = dynamicsymbols('z')
    m = Symbol('m', positive=True)
    m_p = Symbol('m_p', positive=True)
    k_1 = Symbol('k', positive=True)
    c_1 = Symbol('c', positive=True)
    c_2 = Symbol('c_water', positive=True)
    F_load = Symbol('F_{load}', positive=True)
    g = Symbol('g', positive=True)
    p_p = Symbol('P_p', positive=True)
    rho = Symbol('\\rho', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)

    def __init__(self,
                 m=None,
                 m_p=None,
                 k_1=None,
                 c_1=None,
                 c_2=None,
                 F_load=None,
                 ivar=Symbol('t'),
                 z=None,
                 g=None,
                 p_p=None,
                 rho=None,
                 Omega=None,
                 F=None,
                 **kwargs):

        if z is not None: self.z = z
        if m is not None: self.m = m
        if k_1 is not None: self.k_1 = k_1
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if F_load is not None: self.F_load = F_load
        if g is not None: self.g = g
        if rho is not None: self.rho = rho
        if p_p is not None: self.p_p = p_p
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._platform = MaterialPoint(
            self.m, pos1=self.z,
            qs=[self.z])(label='Material point - mass of the platform')
        self._spring_1 = Spring(
            self.k_1,
            pos1=self.z,
            pos2=(self.k_1 * self.z) / (self.p_p * self.rho),
            qs=[self.z])(label='Spring - stiffness of the platform')
        self._spring_2 = Spring(self.p_p * self.rho,
                                pos1=(self.k_1 * self.z) /
                                (self.p_p * self.rho),
                                qs=[self.z])(label='Spring - buoyancy')
        self._damper = Damper(self.c_1,
                              pos1=self.z,
                              pos2=(self.k_1 * self.z) / (self.p_p * self.rho),
                              qs=[self.z
                                  ])(label='Damper - damping of the platform')
        self._water = Damper(self.c_2,
                             pos1=(self.k_1 * self.z) / (self.p_p * self.rho),
                             qs=[self.z
                                 ])(label='Damper - damping of the water')
        self._force = Force(-self.F_load, pos1=self.z,
                            qs=[self.z])(label='Force - load on the platform')
        self._MassOfPlatform = Force(
            -self.m * self.g, pos1=self.z, qs=[
                self.z
            ])(label='Force - weight of the platform in gravitational field')
        self._excitation = Force(
            self.F * sin(self.Omega * self.ivar),
            pos1=(self.k_1 * self.z) / (self.p_p * self.rho),
            qs=[self.z])(label='Force - Force caused by the waves')

        components['_platform'] = self._platform
        components['_spring_1'] = self._spring_1
        components['_spring_2'] = self._spring_2
        components['_damper'] = self._damper
        components['_water'] = self._water
        components['_force'] = self._force
        components['_MassOfPlatform'] = self._MassOfPlatform
        components['_excitation'] = self._excitation

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of a platform',
            self.g: r'Gravity field acceleration',
            self.Omega: r'Excitation frequency',
            self.c_1: r'Platform damping coefficient',
            self.c_2: r'Water damping coefficient',
            self.rho: r'Water density',
            self.k_1: r'Platform stiffness coefficient',
            self.p_p: r'Area submerged in water',
            self.F_load: r'Force caused by a load on a rig',
            self.F: r'Force caused by a waves',
            self.ivar: r'Time'
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, Omega0, k0, p_p0, F0, lam0 = Symbol('m_0', positive=True), Symbol(
            'Omega0', positive=True), Symbol('k0', positive=True), Symbol(
                'p_{p0}',
                positive=True), Symbol('F0',
                                       positive=True), Symbol('lambda0',
                                                              positive=True)
        c_2 = self.c_2
        g = self.g
        rho = self.rho

        default_data_dict = {
            self.m: [S.One * m0 * no * 1000 for no in range(20, 30)],
            self.g: [S.One * g * no for no in range(1, 2)],
            self.Omega: [S.One * Omega0 * no for no in range(1, 2)],
            self.c_2: [S.One * c_2 * no for no in range(1, 2)],
            self.rho: [S.One * rho * no for no in range(1, 2)],
            self.k_1: [S.One * k0 * no * 10000 for no in range(5, 10)],
            self.c_1: [S.One * k0 * lam0 * no * 10000 for no in range(5, 10)],
            self.p_p: [S.One * p_p0 * no for no in range(9, 16)],
            self.F_load: [S.One * F0 * no * 1000 for no in range(5, 10)],
            self.F: [S.One * F0 * no * 1000 for no in range(1, 10)]
        }

        return default_data_dict


class DampedHarmonicOscillator(SpringDamperMassSystem):
    pass






class Winch(ComposedSystem):

    scheme_name = 'sdof_winch.PNG'
    real_name = 'winch_mechanism_real.PNG'
    r = Symbol('r', positive=True)
    l = Symbol('l', positive=True)
    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    ivar = Symbol('t')
    phi = dynamicsymbols('\\varphi')
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)

    def __init__(self,
                 r=None,
                 l=None,
                 m=None,
                 g=None,
                 ivar=None,
                 phi=None,
                 Omega=None,
                 F=None,
                 **kwargs):

        if r is not None: self.r = r
        if l is not None: self.l = l
        if m is not None: self.m = m
        if g is not None: self.g = g
        if ivar is not None: self.ivar = ivar
        if phi is not None: self.phi = phi
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        x = self.r * cos(self.phi) + (self.l + self.r * self.phi) * sin(
            self.phi)
        y = -self.r * sin(self.phi) + (self.l + self.r * self.phi) * cos(
            self.phi)

        self.material_point_1 = MaterialPoint(self.m, x, qs=[self.phi])
        self.material_point_2 = MaterialPoint(self.m, y, qs=[self.phi])
        self.gravity = GravitationalForce(self.m,
                                          self.g,
                                          pos1=-y,
                                          qs=[self.phi])
        self.force = Force(-self.F * (self.l + self.r * self.phi) *
                           sin(self.Omega * self.ivar),
                           pos1=self.phi,
                           qs=[self.phi])
        system = self.material_point_1 + self.material_point_2 + self.gravity + self.force

        super().__init__(system, **kwargs)

    def linearized(self):

        return type(self).from_system(super().linearized())

    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r'Winch radius',
            self.l: r'Winch length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.r: [l0, 2 * l0, 4 * l0, l0, 8 * l0],
            self.m:
            [2 * m0, S.Half * m0, 4 * m0, m0, S.Half**2 * m0, 8 * m0, 16 * m0],
            self.l: [
                2 * l0, S.Half * l0, 4 * l0, S.Half**2 * l0, 3 * l0,
                3 * S.Half * l0, 9 * l0, 3 * S.Half**2 * l0
            ],
        }
        return default_data_dict

    def max_static_force_pin(self):
        return (self.m * self.g).subs(self._given_data)

    def max_dynamic_force_pin(self):

        omg_amp = ComposedSystem(
            self.linearized()).frequency_response_function() * self.Omega

        return (self.m * self.l * (omg_amp)**2 + self.max_static_force_pin())

    def static_force_pin_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kr * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kr * Re))**(1 / 2)


class NonlinearEngine(ComposedSystem):
    scheme_name = 'nonline_engine_angled_springs.png'
    real_name = 'engine_real.PNG'
    """
    Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            m_e = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            d = distance
                -distance between mounting point and engine

            k_m = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            beta = angle
                -angle of springs that hold the system

            l0 = length
                -Initial length of non-linear springs

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        >>> t = symbols('t')
        >>> M, me, e, km, beta, l0 = symbols('M, m_e, e, k_m, beta, l_0')
        >>> qs = dynamicsymbols('z') 
        >>> SDoFNonlinearEngine()
    """
    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    d = Symbol('d', positive=True)
    l_0 = Symbol('l_0', positive=True)
    z = dynamicsymbols('z')
    phi = dynamicsymbols('phi')
    Omega = Symbol('Omega', positive=True)
    g = Symbol('g', positive=True)

    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 e=None,
                 d=None,
                 l_0=None,
                 z=None,
                 phi=None,
                 ivar=Symbol('t', positive=True),
                 Omega=None,
                 g=None,
                 **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if d is not None: self.d = d
        if phi is not None: self.phi = phi
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if l_0 is not None: self.l_0 = l_0
        if z is not None: self.z = z
        if Omega is not None: self.Omega = Omega
        if g is not None: self.g = g

        self.materialPoint_1 = MaterialPoint(self.M, self.z, qs=[self.z])
        self.materialPoint_2 = MaterialPoint(self.m_e,
                                             self.z + self.e * cos(self.phi),
                                             qs=[self.z])
        self.spring = Spring(2 * self.k_m,
                             pos1=(self.z**2 + self.d**2)**0.5 - self.l_0,
                             qs=[self.z])
        self.gravity_force1 = GravitationalForce(self.M,
                                                 self.g,
                                                 self.z,
                                                 qs=[self.z])
        self.gravity_force2 = GravitationalForce(self.m_e,
                                                 self.g,
                                                 self.z +
                                                 self.e * cos(self.phi),
                                                 qs=[self.z])
        system = self.spring + self.materialPoint_1 + self.materialPoint_2 + self.gravity_force1 + self.gravity_force2
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
            self.l_0: r'',
            self.beta: r'',
        }

        return self.sym_desc_dict

    def linearized(self):

        return type(self).from_system(super().linearized())

    def max_static_force(self):
        return abs(self.static_load().doit()[0] / 2)

    def max_dynamic_force(self):
        return self.frequency_response_function() * self.stiffness_matrix(
        )[0] + self.max_static_force()

    def get_default_data(self):

        m0, k0, e0, l0 = symbols('m_0 k_0 e_0 l_0', positive=True)

        default_data_dict = {
            self.M: [
                100 * m0, 300 * m0, 500 * m0, 700 * m0, 900 * m0, 200 * m0,
                400 * m0, 600 * m0, 800 * m0
            ],
            self.m_e: [
                m0, 3 * m0, 5 * m0, 7 * m0, 9 * m0, 2 * m0, 4 * m0, 6 * m0,
                8 * m0
            ],
            self.k_m: [
                k0, 2 * k0, 4 * k0, 6 * k0, 8 * k0, 3 * k0, 5 * k0, 7 * k0,
                9 * k0
            ],
            self.e: [
                2 * e0, S.Half * e0, 4 * e0, S.Half**2 * e0, 3 * e0,
                3 * S.Half * e0, 9 * e0, 3 * S.Half**2 * e0
            ],
            self.phi: [self.Omega * self.ivar],
            self.d:
            [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0, 9 * l0]
        }
        return default_data_dict


class StraightNonlinearEngine(NonlinearEngine):
    """
    Model of an exemplary Engine with nonlinear suspension aligned horizontally.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            me = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            km = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            beta = angle
                -angle of springs that hold the system

            l0 = length
                -Initial length of non-linear springs

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        >>> t = symbols('t')
        >>> M, me, e, km, beta, l0 = symbols('M, m_e, e, k_m, beta, l_0')
        >>> qs = dynamicsymbols('z') 
        >>> SDoFNonlinearEngine()
    """
    scheme_name = 'non_linear_engine.png'

    def get_default_data(self):

        m0, k0, e0, l, omega = symbols('m_0 k_0 e_0 l Omega', positive=True)

        default_data_dict = {
            self.M: [
                100 * m0, 300 * m0, 500 * m0, 700 * m0, 900 * m0, 200 * m0,
                400 * m0, 600 * m0, 800 * m0
            ],
            self.m_e: [
                m0, 3 * m0, 5 * m0, 7 * m0, 9 * m0, 2 * m0, 4 * m0, 6 * m0,
                8 * m0
            ],
            self.k_m: [
                k0, 2 * k0, 4 * k0, 6 * k0, 8 * k0, 3 * k0, 5 * k0, 7 * k0,
                9 * k0
            ],
            self.e: [
                2 * e0, S.Half * e0, 4 * e0, S.Half**2 * e0, 3 * e0,
                3 * S.Half * e0, 9 * e0, 3 * S.Half**2 * e0
            ],
            self.l_0: [S.Half * l, l, 2 * l],
            self.d: [4 * l, 8 * l],
            self.beta: [S.Half * pi],
            self.phi: [omega * self.ivar]
        }
        return default_data_dict

    
from .mechanics.engine import (
    Engine, DampedEngine, BoxerEnginePerpendicularSprings,
    InlineEnginePerpendicularSprings, EngineVerticalSpringGravity,
    DampedEngineVerticalSpringGravity,
    NonLinearInlineEnginePerpendicularSpringsGravity,
    NonLinearBoxerEnginePerpendicularSprings)
    

class ForcedNonLinearTrolley(ComposedSystem):
    scheme_name = 'sdof_nonlin_trolley.PNG'
    real_name = 'tension_leg_platform.png'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t', positive=True),
                 F=Symbol('F_0', positive=True),
                 x=dynamicsymbols('x'),
                 Omega=Symbol('Omega', positive=True),
                 **kwargs):
        """
        Model of Single Degree of Freedom Trolley with nonlinear spring (type of inverted pendulum)

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            l = length
                -actual length of the non-linear spring

            l_0 = lenght
                -Initial length of non-linear spring

            k = spring coefficient
                -value of spring coefficient

            F = Force
                -Trolley's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, l, l_0, k, F = symbols('m, l, l_0, k, F')
        >>> qs = dynamicsymbols('x') # Generalized Coordinate
        >>> Omega = symbols('Omega')
        >>> SDoFTrolleyWithNonlinearSpring()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFTrolleyWithNonlinearSpring()
    """

        self.m = m
        self.k = k
        self.d = d
        self.l_0 = l_0
        self.F = F

        self.MaterialPoint = MaterialPoint(m, x, qs=[x])
        self.Spring = Spring(k, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])
        self.Force = Force(F * cos(Omega * ivar), pos1=x, qs=[x])

        system = self.MaterialPoint + self.Spring + self.Force
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.k: 'Spring Stiffness',
            self.l: r'length',
            self.l_0: r'length',
            self.F: r'Force',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
            self.d: [
                1 * l0, 2 * l0, S.Half * l0, 3 * S.Half * l0, 4 * l0, 5 * l0,
                6 * l0, 7 * l0, 8 * l0, 9 * l0
            ],
            self.k: [
                S.Half * k0, 2 * k0, 1 * k0, 3 * S.Half * k0, 4 * k0, 5 * k0,
                6 * k0, 7 * k0, 8 * k0, 9 * k0, 3 * k0
            ],
            self.l_0: [l0]
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if parameters_dict[self.k] - parameters_dict[self.k] * parameters_dict[
                self.l_0] / parameters_dict[self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]
        return parameters_dict


class NonLinearTrolley(NonlinearComposedSystem):

    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x'),
                 **kwargs):

        self.m = m
        self.k = k
        self.d = d
        self.l_0 = l_0
        self.x = x

        self.trolley = MaterialPoint(m, x, qs=[x]) + Spring(
            k, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])

        super().__init__(self.trolley, **kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
            self.d: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0
            ],
            self.k: [
                S.Half * k0, 2 * k0, 1 * k0, 3 * S.Half * k0, 4 * k0, 5 * k0,
                6 * k0, 7 * k0, 8 * k0, 9 * k0, 3 * k0
            ],
            self.l_0: [l0]
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if parameters_dict[self.k] - parameters_dict[self.k] * parameters_dict[
                self.l_0] / parameters_dict[self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]
        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Trolley Mass',
            self.k: 'Spring Stiffness',
            self.d: r'length',
            self.l_0: r'length',
        }
        return self.sym_desc_dict


class NonLinearDisc(NonlinearComposedSystem):
    scheme_name = 'nonlinear_disc.png'
    real_name = 'roller_tightener.png'

    def __init__(self,
                 m1=Symbol('m', positive=True),
                 kl=Symbol('k', positive=True),
                 R=Symbol('R', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x'),
                 **kwargs):

        self.m1 = m1
        self.kl = kl
        self.R = R
        self.l_0 = l_0
        self.d = d
        self.x = x

        self.Disk1 = MaterialPoint(m1, x, qs=[x]) + MaterialPoint(
            m1 / 2 * R**2, x / R, qs=[x]) + Spring(
                kl, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])

        system = self.Disk1
        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [
                0.5 * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
            self.d: [
                5 * l0, 2 * l0, 3 * S.Half * l0, 4 * l0, 6 * l0, 7 * l0,
                8 * l0, 9 * l0
            ],
            self.kl: [
                1 * k0, 3 * k0, 2 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0,
                9 * k0
            ],
            self.l_0: [l0],
        }

        return default_data_dict




#Sav
class Shaft(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents the Kinematicly excited shaft with two disks.
    =========
            I = Moment of Inertia
                -Moment of Inertia in case of both disc

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

    >>> t = symbols('t')
    >>> I, k1, k2 = symbols('I, k_1, k_2')
    >>> qs = dynamicsymbols('phi_1, phi_2') # Generalized Coordinates
    >>> DDoFShaft()

    -defines the symbols and dynamicsymbols
    -finally determines the instance of the system using class DDoFShaft
    """

    scheme_name = 'shaft_sdof_scheme.png'
    real_name = 'ddof_shaft_real.png'
    detail_scheme_name = 'parallel_key_load.png'
    detail_real_name = 'shaft_with_key.png'

    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    I = Symbol('I', positive=True)
    l_1 = Symbol('l_1', positive=True)
    l_2 = Symbol('l_2', positive=True)
    I_1 = Symbol('I_1', positive=True)
    I_2 = Symbol('I_2', positive=True)
    Ms = Symbol('M_s', positive=True)
    Omega = Symbol('Omega', positive=True)

    theta = dynamicsymbols('theta')
    phi = dynamicsymbols('\\varphi')

    def __init__(self,
                 l0=None,
                 G=None,
                 I=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 Ms=None,
                 phi=None,
                 theta=None,
                 ivar=Symbol('t'),
                 qs=None,
                 **kwargs):
        if G is not None: self.G = G

        if I is not None: self.I = I
        if Ms is not None: self.Ms = Ms
        #if Omega is not None: self.Omega = Omega
        if l_1 is not None: self.l_1 = l_1
        if l_2 is not None: self.l_2 = l_2
        if I_1 is not None: self.I_1 = I_1
        if I_2 is not None: self.I_2 = I_2
        if phi is not None: self.phi = phi
        if theta is not None: self.theta = theta

        theta = self.theta
        self.qs = [self.phi]
        self.ivar = ivar

        self.k_1 = (self.G * self.I_1) / self.l_1
        self.k_2 = (self.G * self.I_2) / self.l_2

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_2 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1),
                               pos1=self.phi,
                               pos2=theta,
                               qs=self.qs)  # right spring
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)
        system = self.disc_1 + self.spring_2 + self.moment
        self.system = system

        super(HarmonicOscillator,self).__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, G, l = symbols('m_0 l_0 G l', positive=True)
        theta0, Omega = symbols('theta_0, Omega', positive=True)

        default_data_dict = {
            self.I: [S.Half * m0 * (l0**2) * no for no in range(1, 3)],
            self.I_1: [S.Half**(no) * (l0**4) for no in range(1, 8)],
            self.I_2: [S.Half**no * (l0**4) for no in range(1, 8)],
            self.l_1: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.l_2: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.theta: [theta0 * cos(Omega * self.ivar)],
        }

        return default_data_dict

    def disc_force(self):
        t = self.ivar
        return self.I * self.steady_solution().diff(t, t)

    def max_static_force_pin(self):
        d = Symbol('d', positive=True)
        return 2 * self.Ms / d

    def max_dynamic_force_pin(self):
        d = Symbol('d', positive=True)
        return self.frequency_response_function(
            self.natural_frequencies()[0]) * self.stiffness_matrix()[0]

    def max_static_bearing_force(self):
        d = Symbol('d', positive=True)
        return abs(2 * self.static_load()[0] / d)

    def max_dynamic_bearing_force(self):
        d = Symbol('d', positive=True)
        acc_amp = self.frequency_response_function() * self.Omega**2

        return abs(
            2 * (self.I * acc_amp) /
            d) + self.max_static_bearing_force()  #.subs(self._given_data)

    def static_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_static_bearing_force()) / (kd * h)

    def dynamic_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_dynamic_bearing_force()) / (kd * h)


#############################################################################################
class SDOFDampedShaft(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents the Kinematicly excited shaft with two disks.
    =========
            I = Moment of Inertia
                -Moment of Inertia in case of both disc

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

    >>> t = symbols('t')
    >>> I, k1, k2 = symbols('I, k_1, k_2')
    >>> qs = dynamicsymbols('phi_1, phi_2') # Generalized Coordinates
    >>> DDoFShaft()

    -defines the symbols and dynamicsymbols
    -finally determines the instance of the system using class DDoFShaft
    """

    scheme_name = 'damped_shaft_phi.png'
    real_name = 'ddof_shaft_real.png'
    detail_scheme_name = 'parallel_key_load.png'
    detail_real_name = 'shaft_with_key.png'

    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    I = Symbol('I', positive=True)
    l_1 = Symbol('l_1', positive=True)
    l_2 = Symbol('l_2', positive=True)
    I_1 = Symbol('I_1', positive=True)
    I_2 = Symbol('I_2', positive=True)
    c_1 = Symbol('c_1', positive=True)
    c_2 = Symbol('c_2', positive=True)
    k_1 = Symbol('k_1', positive=True)
    k_2 = Symbol('k_2', positive=True)
    Ms = Symbol('M_s', positive=True)
    Omega = Symbol('Omega', positive=True)

    theta = dynamicsymbols('theta')
    phi = dynamicsymbols('\\varphi')

    def __init__(self,
                 l0=None,
                 G=None,
                 I=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 c_1=None,
                 c_2=None,
                 k_1=None,
                 k_2=None,
                 Ms=None,
                 phi=None,
                 theta=None,
                 ivar=Symbol('t'),
                 qs=None,
                 **kwargs):
        if G is not None: self.G = G

        if I is not None: self.I = I
        if Ms is not None: self.Ms = Ms
        #if Omega is not None: self.Omega = Omega
        if l_1 is not None: self.l_1 = l_1
        if l_2 is not None: self.l_2 = l_2
        if I_1 is not None: self.I_1 = I_1
        if I_2 is not None: self.I_2 = I_2
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if k_1 is not None: self.k_1 = k_1
        if k_2 is not None: self.k_2 = k_2
        if phi is not None: self.phi = phi
        if theta is not None: self.theta = theta

        theta = self.theta
        self.qs = [self.phi]
        self.ivar = ivar

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_1 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1),
                               pos1=self.phi,
                               pos2=self.theta,
                               qs=self.qs)
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)
        self.damper_1 = Damper(self.c_1 * self.c_2 / (self.c_2 + self.c_1),
                               pos1=self.phi,
                               pos2=self.theta,
                               qs=self.qs)

        system = self.disc_1 + self.spring_1 + self.moment + self.damper_1
        self.system = system

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, G, l, lamb = symbols('m_0 l_0 G l lambda', positive=True)
        theta0, Omega = symbols('theta_0, Omega', positive=True)

        default_data_dict = {
            self.c_1: [lamb * self.k_1],
            self.c_2: [lamb * self.k_2],
            self.k_1: [(self.G * self.I_1) / self.l_1],
            self.k_2: [(self.G * self.I_2) / self.l_2],
            self.I: [S.Half * m0 * (l0**2) * no for no in range(1, 3)],
            self.I_1: [S.Half**(no) * (l0**4) for no in range(1, 8)],
            self.I_2: [S.Half**no * (l0**4) for no in range(1, 8)],
            self.l_1: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.l_2: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.theta: [theta0 * cos(Omega * self.ivar)],
        }

        return default_data_dict

    def disc_force(self):
        t = self.ivar
        return self.I * self.steady_solution().diff(t, t)

    def max_static_force_pin(self):
        d = Symbol('d', positive=True)
        return 2 * self.Ms / d

    def max_dynamic_force_pin(self):
        d = Symbol('d', positive=True)
        return self.frequency_response_function(
            self.natural_frequencies()[0]) * self.stiffness_matrix()[0]

    def max_static_bearing_force(self):
        d = Symbol('d', positive=True)
        return abs(2 * self.static_load()[0] / d)

    def max_dynamic_bearing_force(self):
        d = Symbol('d', positive=True)
        acc_amp = self.frequency_response_function() * self.Omega**2

        return abs(
            2 * (self.I * acc_amp) /
            d) + self.max_static_bearing_force()  #.subs(self._given_data)

    def static_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_static_bearing_force()) / (kd * h)

    def dynamic_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_dynamic_bearing_force()) / (kd * h)


#Dominik
class KinematicClutchWithSprings(ComposedSystem):
    """NotReady to use sample Double Degree of Freedom System represents the Kinematicly excited clutch with spring between input shaft and disc.
    =========
            I = Moment of Inertia
                -Moment of Inertia in case of both disc

            k = Right spring coefficient
                -Right spring carrying the system
            n = Number of springs
            
            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

    Example
    =======
    A disc oscillating while being held up by a springs with number n and spring constant k

    >>> t = symbols('t')
    >>> I, k1, k2 = symbols('I, k')
    >>> qs = dynamicsymbols('phi_1, phi_2') # Generalized Coordinates
    >>> SDOFClutch()

    -defines the symbols and dynamicsymbols
    -finally determines the instance of the system using class DDoFShaft
    """

    scheme_name = ''
    real_name = ''
    detail_scheme_name = ''
    detail_real_name = ''

    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    I = Symbol('I', positive=True)
    l_1 = Symbol('l_1', positive=True)
    l_2 = Symbol('l_2', positive=True)
    I_1 = Symbol('I_1', positive=True)
    I_2 = Symbol('I_2', positive=True)
    Ms = Symbol('M_s', positive=True)
    Omega = Symbol('Omega', positive=True)

    theta = dynamicsymbols('theta')
    phi = dynamicsymbols('\\varphi')

    def __init__(self,
                 l0=None,
                 G=None,
                 I=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 Ms=None,
                 phi=None,
                 theta=None,
                 ivar=Symbol('t'),
                 qs=None,
                 **kwargs):
        if G is not None: self.G = G
        if I is not None: self.I = I
        if Ms is not None: self.Ms = Ms
        if l_1 is not None: self.l_1 = l_1
        if l_2 is not None: self.l_2 = l_2
        if I_1 is not None: self.I_1 = I_1
        if I_2 is not None: self.I_2 = I_2
        if phi is not None: self.phi = phi
        if theta is not None: self.theta = theta

        theta = self.theta
        self.qs = [self.phi]
        self.ivar = ivar

        self.k_1 = (self.G * self.I_1) / self.l_1
        self.k_2 = (self.G * self.I_2) / self.l_2

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_2 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1),
                               pos1=self.phi,
                               pos2=theta,
                               qs=self.qs)  # right spring
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)
        system = self.disc_1 + self.spring_2 + self.moment
        self.system = system

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, G, l = symbols('m_0 l_0 G l', positive=True)
        theta0, Omega = symbols('theta_0, Omega', positive=True)

        default_data_dict = {
            self.I: [S.Half * m0 * (l0**2) * no for no in range(1, 3)],
            self.I_1: [S.Half**(no) * (l0**4) for no in range(1, 8)],
            self.I_2: [S.Half**no * (l0**4) for no in range(1, 8)],
            self.l_1: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.l_2: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.theta: [theta0 * cos(Omega * self.ivar)],
        }

        return default_data_dict

    def disc_force(self):
        t = self.ivar
        return self.I * self.steady_solution().diff(t, t)

    def max_static_force_pin(self):
        d = Symbol('d', positive=True)
        return 2 * self.Ms / d

    def max_dynamic_force_pin(self):
        d = Symbol('d', positive=True)
        return self.frequency_response_function(
            self.natural_frequencies()[0]) * self.stiffness_matrix()[0]

    def max_static_bearing_force(self):
        d = Symbol('d', positive=True)
        return abs(2 * self.static_load()[0] / d)

    def max_dynamic_bearing_force(self):
        d = Symbol('d', positive=True)
        acc_amp = self.frequency_response_function() * self.Omega**2

        return abs(
            2 * (self.I * acc_amp) /
            d) + self.max_static_bearing_force()  #.subs(self._given_data)

    def static_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_static_bearing_force()) / (kd * h)

    def dynamic_key_length(self):

        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_dynamic_bearing_force()) / (kd * h)


class DampedMeasuringTool(ComposedSystem):

    scheme_name = 'measure_tool.PNG'
    real_name = 'measure_tool_real.PNG'
    #detail_scheme_name =
    #detail_real_name =

    m = Symbol('m', positive=True)
    l = Symbol('l', positive=True)
    k = Symbol('k', positive=True)
    k_t = Symbol('k_t', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    c = Symbol('c', positive=True)
    c_t = Symbol('c_t', positive=True)
    lam = Symbol('lambda', positive=True)
    l0 = Symbol('l_0', positive=True)
    lam0 = Symbol('lambda_0', positive=True)

    def __init__(self,
                 m=None,
                 l=None,
                 k=None,
                 k_t=None,
                 ivar=Symbol('t'),
                 Omega=None,
                 F=None,
                 phi=None,
                 qs=None,
                 c=None,
                 c_t=None,
                 lam=None,
                 l0=None,
                 lam0=None,
                 **kwargs):
        if l is not None: self.l = l
        if m is not None: self.m = m
        if k is not None: self.k = k
        if k_t is not None: self.k_t = k_t
        if F is not None: self.F = F
        if Omega is not None: self.Omega = Omega
        if phi is not None: self.phi = phi
        if c is not None: self.c = c
        if c_t is not None: self.ct = ct
        if lam is not None: self.lam = lam
        if l0 is not None: self.l0 = l0
        if lam0 is not None: self.lam0 = lam0

        self.qs = [self.phi]
        self.ivar = ivar

        self._moment_of_inertia = MaterialPoint(
            (S.One / 3) * self.m * self.l**2, self.phi, qs=[self.phi])
        self._upper_spring = Spring(self.k,
                                    pos1=self.l * self.phi,
                                    qs=[self.phi])
        self._lower_spring = Spring(self.k,
                                    pos1=self.l * self.phi,
                                    qs=[self.phi])
        self._spiral_spring = Spring(self.k_t, self.phi, qs=[self.phi])
        self._force = Force(self.F * self.l, pos1=self.phi)
        self._springs_damping = Damper(2 * self.c,
                                       pos1=self.l * self.phi,
                                       qs=[self.phi])
        self._spiral_spring_damping = Damper(self.c_t,
                                             pos1=self.phi,
                                             qs=[self.phi])
        composed_system = self._moment_of_inertia + self._upper_spring + self._lower_spring + self._spiral_spring + self._force + self._springs_damping + self._spiral_spring_damping

        super().__init__(composed_system, **kwargs)

    def get_default_data(self):

        m0, k0, F0, Omega0, lam0, l0 = self.m0, self.k0, self.F0, self.Omega0, self.lam0, self.l0

        default_data_dict = {
            self.c: [self.lam * (self.k)],
            self.c_t: [self.lam * (self.k_t)],
            self.m: [m0 * S.One * no for no in range(1, 8)],
            self.k: [k0 * S.One * no for no in range(1, 8)],
            self.k_t: [k0 * l0**2 * S.One * no for no in range(1, 8)],
            self.F: [
                F0 * S.One * no * cos(self.Omega * self.ivar)
                for no in range(1, 8)
            ],
            self.Omega: [self.Omega],
            self.lam: [self.lam],
            self.l: [l0 * S.One * no for no in range(1, 8)],
        }

        return default_data_dict

    def steady_state(self):
        return 3 * (S.One / 2 * self.damping_coefficient())**(-1)

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force_pin(self):
        lin_sys = self.linearized()

        dyn_comp = (lin_sys.frequency_response_function() * self.l *
                    self.k).subs(self._given_data)

        total_force = (dyn_comp + self.max_static_force_pin())

        return total_force

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


# class DampedShaft(ComposedSystem):

#     def __init__(self,
#                  I=Symbol('I', positive=True),
#                  k_2=Symbol('k_2', positive=True),
#                  k_1=Symbol('k_1', positive=True),
#                  c_1=Symbol('c_1', positive=True),
#                  c_2=Symbol('c_1', positive=True),
#                  input_displacement=dynamicsymbols('theta'),
#                  ivar=Symbol('t'),
#                  qs=dynamicsymbols('\\varphi_1, \\varphi_2'),
#                  **kwargs):

#         phi1, phi2 = qs
#         theta = input_displacement

#         self.k_2 = k_2  # left spring
#         self.k_1 = k_1  # right spring
#         self.c_1 = c_1  # right spring
#         self.c_2 = c_2  # right spring
#         self.I = I  # moment of inertia of a rod
#         self.input_displacement = input_displacement
#         self.qs = qs

#         self.disc_1 = Disk(I, pos1=phi1, qs=qs)
#         self.spring_1 = Spring(k_2, phi1, phi2, qs=qs)  # left spring
#         self.disc_2 = Disk(I, pos1=phi2, qs=qs)
#         self.spring_2 = Spring(k_1, pos1=phi2, pos2=theta,
#                                qs=qs)  # right spring
#         self.damper_1 = Damper(c_2, phi1, phi2, qs=qs)  # left spring
#         self.damper_2 = Damper(c_1, pos1=phi2, pos2=theta,
#                                qs=qs)  # right spring
#         system = self.disc_1 + self.disc_2 + self.spring_1 + self.spring_2 + self.damper_1 + self.damper_2

#         super().__init__(system,**kwargs)

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.I: r'Moment of Inertia',
#             self.k_1: r'',
#             self.k_2: r'',
#         }
#         return self.sym_desc_dict
#     def get_default_data(self):

#         I0, k0, lamb = symbols('I_0 k_0 lambda', positive=True)

#         default_data_dict = {
#             self.k_2: [2 * k0, 4 * k0,6*k0,8*k0,10*k0],
#             self.k_1: [k0, 3 * k0,5*k0,7*k0,9*k0],
#             self.I: [2 * I0, S.Half * I0, 4 * I0, S.Half**2 * I0,3 * I0,3* S.Half * I0, 9 * I0, 3*S.Half**2 * I0],
#             self.c_1: [lamb * self.k_1],
#             self.c_2: [lamb * self.k_2],
#         }
#         return default_data_dict

#     def get_random_parameters(self):

#         default_data_dict = self.get_default_data()

#         parameters_dict = {
#             key: random.choice(items_list)
#             for key, items_list in default_data_dict.items()
#             }

#         return parameters_dict
