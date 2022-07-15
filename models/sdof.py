from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from ..continuous import ContinuousSystem, PlaneStressProblem

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
    detail_scheme_name = 'damped_car_new.PNG'
    detail_real_name = 'car_real.jpg'
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


    def _init_from_components(self,*args,system=None,**kwargs):
        
        if system is None:
            composed_system = self._elements_sum
        else:
            composed_system = system
            
        print('CS',composed_system._components)
        super().__init__(None,system = composed_system)
        
        print('self',self._components)
        if self._components is None:
            comps = {}
        else:
            comps=self._components

        self._components = {**comps,**self.components}

    def __init__(self,
                 Lagrangian=None,
                 m0 = None,
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
        
        
        return {**super().components,**self.components}



    
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


#Mateusz
class BlowerToothedBelt(ComposedSystem):

    scheme_name = 'blower_toothed_belt.png'
    real_name = 'blown_440_big_block.jpg'
    detail_scheme_name = 'blower_roller_bolt.png'
    detail_real_name = 'tensioner_pulley.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k_belt=Symbol('k_b', positive=True),
                 k_tensioner=Symbol('k_t', positive=True),
                 ivar=Symbol('t'),
                 Omega=Symbol('Omega', positive=True),
                 F=Symbol('F', positive=True),
                 Q=Symbol('Q', positive=True),
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.z = z
        self.m = m
        self.k_belt = k_belt
        self.k_tensioner = k_tensioner
        self.F = F
        self.Q = Q
        self.P0 = Symbol('P_0', positive=True)
        self.Omega = Omega
        self.mass = MaterialPoint(m, z, qs=[z])
        self.upper_belt = Spring(k_belt, z, qs=[z])
        self.lower_belt = Spring(k_belt, z, qs=[z])
        self.tensioner = Spring(k_tensioner, z, qs=[z])
        self.force = Force(F * sin(Omega * ivar), pos1=z) + Force(-Q, pos1=z)
        
        
        composed_system = self.mass + self.upper_belt + self.lower_belt + self.tensioner + self.force

        super().__init__(composed_system, **kwargs)

    def get_default_data(self):

        m0, k0, F0, Omega0 = symbols('m_0 k_0 F_0 Omega_0', positive=True)

        default_data_dict = {
            self.m: [0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0],
            self.k_belt: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0],
            self.k_tensioner: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0],
            self.F: [F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0],
            #self.Omega: [Omega0, 2 * Omega0, 3 * Omega0, 4 * Omega0, 5 * Omega0, 6 * Omega0],
            self.Q: [2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0, F0],
        }

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_belt: r'Belt stiffnes',
            self.k_tensioner: r'Tensioner spring stiffnes',
        }

        return self.sym_desc_dict

#Mateusz
class DampedBlowerToothedBelt(ComposedSystem):

    scheme_name = 'damped_blower_toothed_belt.png'
    real_name = 'blown_440_big_block.jpg'
    detail_scheme_name = 'blower_roller_bolt.png'
    detail_real_name = 'tensioner_pulley.jpg'

    m = Symbol('m', positive=True)
    k_belt = Symbol('k_b', positive=True)
    k_tensioner = Symbol('k_t', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    #Q=Symbol('Q', positive=True)
    F0 = Symbol('F_0', positive=True)
    z = dynamicsymbols('z')
    c_belt = Symbol('c_b', positive=True)
    c_tensioner = Symbol('c_t', positive=True)
    lam = Symbol('lambda', positive=True)
    z0 = Symbol('z0', positive=True)

    lam0 = Symbol('lambda_0', positive=True)

    def __init__(
            self,
            m=None,
            k_belt=None,
            k_tensioner=None,
            ivar=Symbol('t'),
            Omega=None,
            F=None,
            #Q=None,
            c_belt=None,
            c_tensioner=None,
            lam=None,
            z=None,
            qs=None,
            P0=None,
            z0=None,
            **kwargs):
        if z is not None: self.z = z
        if m is not None: self.m = m
        if k_belt is not None: self.k_belt = k_belt
        if k_tensioner is not None: self.k_tensioner = k_tensioner
        if F is not None: self.F = F
        if P0 is not None: self.P0 = P0
        if Omega is not None: self.Omega = Omega
        if lam is not None: self.lam = lam
        if c_belt is not None: self.c_belt = c_belt
        if c_tensioner is not None: self.c_tensioner = c_tensioner
        if z0 is not None: self.z0 = z0

        self.qs = [self.z]
        self.ivar = ivar

        self._mass = MaterialPoint(self.m, self.z, qs=[self.z])
        self._upper_belt = Spring(self.k_belt,
                                  self.z,
                                  pos2=self.z0,
                                  qs=[self.z])
        self._lower_belt = Spring(self.k_belt,
                                  self.z,
                                  pos2=self.z0,
                                  qs=[self.z])
        self._tensioner = Spring(self.k_tensioner,
                                 self.z,
                                 pos2=self.z0,
                                 qs=[self.z])
        self._force = Force(self.F * sin(self.Omega * self.ivar),
                            pos1=self.z)  # + Force(-self.Q, pos1=self.z)
        self._damping = Damper(2 * self.c_belt + self.c_tensioner,
                               pos1=self.z,
                               qs=[self.z])
        composed_system = self._mass + self._upper_belt + self._lower_belt + self._tensioner + self._force + self._damping

        super().__init__(composed_system, **kwargs)

    def get_default_data(self):

        m0, k0, F0, Omega0, lam0, z0 = self.m0, self.k0, self.F0, self.Omega0, self.lam0, self.z0

        default_data_dict = {
            self.c_belt: [self.lam * (self.k_belt)],
            self.c_tensioner: [self.lam * (self.k_tensioner)],
            self.m: [m0 / 10 * S.One * no for no in range(1, 8)],
            self.k_belt: [k0 * S.One * no for no in range(1, 8)],
            self.k_tensioner: [k0 * S.One * no for no in range(4, 8)],
            self.F: [F0 * S.One * no for no in range(1, 8)],
            #self.Omega: [Omega0*S.One*no for no in range(1,2)],
            #self.Q: [Q*S.One*no for no in range(1,8)],
            #self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.z0: [F0 / k0 * S.One / 4 * no for no in range(1, 2)],
        }

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_belt: r'Belt stiffnes',
            self.k_tensioner: r'Tensioner spring stiffnes',
        }

        return self.sym_desc_dict

    def tensioner_belt_force(self):
        return self.k_tensioner * self.steady_solution()

    def left_belt_force(self):
        return self.k_belt * self.steady_solution()

    def right_belt_force(self):
        return self.k_belt * self.steady_solution()

    def tensioner_damper_force(self):
        return self.c_tensioner * self.steady_solution().diff(self.ivar)

    def left_belt_damper_force(self):
        return self.c_belt * self.steady_solution().diff(self.ivar)

    def right_belt_damper_force(self):
        return self.c_belt * self.steady_solution().diff(self.ivar)

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


#Szymek #Grześ
class EngineVerticalSpringGravity(ComposedSystem):
    scheme_name = 'engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'

    M = Symbol('M', positive=True)
    m_e = Symbol('m_e', positive=True)
    phi = dynamicsymbols('varphi')
    g = Symbol('g', positive=True)
    k_m = Symbol('k_m', positive=True)
    c_m = Symbol('c_m', positive=True)
    e = Symbol('e', positive=True)
    z = dynamicsymbols('z')

    def __init__(self,
                 M=None,
                 m_e=None,
                 phi=None,
                 g=None,
                 q=None,
                 z=None,
                 k_m=None,
                 c_m=None,
                 e=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if phi is not None: self.phi = phi
        if g is not None: self.g = g
        if k_m is not None: self.k_m = k_m
        if c_m is not None: self.c_m = c_m
        if e is not None: self.e = e
        if z is not None: self.z = z

        self.qs = [self.z]
        self.ivar = ivar

        self.mass = MaterialPoint(self.M, pos1=self.z, qs=[self.z])
        self.crank = MaterialPoint(self.m_e,
                                   pos1=self.z + self.e * cos(self.phi),
                                   qs=[self.z])

        self.left_junction = Spring(self.k_m, self.z, qs=[self.z])

        #self.left_damper = Damper(self.c_m,self.z,qs=[self.z])

        self.right_junction = Spring(self.k_m, self.z, qs=[self.z])

        #self.right_damper = Damper(self.c_m,self.z,qs=[self.z])

        composed_system = self.mass + self.left_junction + self.right_junction + self.crank  # + self.left_damper + self.right_damper
        super().__init__(composed_system, **kwargs)

    def get_default_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)

        default_data_dict = {
            self.c_m: [lam * self.k_m],
            self.M: [m0 * no * 10 for no in range(5, 8)],
            self.m_e: [m0 * no for no in range(1, 8)],
            self.k_m: [k0 * no for no in range(1, 8)],
            self.e: [S.One / 10 * e0 * no for no in range(1, 8)],
        }

        return default_data_dict

    def right_spring_force(self):
        return self.k_m * self.steady_solution()

    def right_damper_force(self):
        return self.c_m * self.steady_solution().diff(self.ivar)

    def left_spring_force(self, ):
        return self.k_m * self.steady_solution()

    def left_damper_force(self):
        return self.c_m * self.steady_solution().diff(self.ivar)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_belt: r'Belt stiffnes',
            self.k_tensioner: r'Tensioner spring stiffnes',
        }

        return self.sym_desc_dict


class DampedEngineVerticalSpringGravity(ComposedSystem):
    scheme_name = 'damped_engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l', positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 Omega=Symbol('Omega', positive=True),
                 phi=dynamicsymbols('varphi'),
                 ivar=Symbol('t'),
                 c_m=Symbol('c_m', positive=True),
                 g=Symbol('g', positive=True),
                 lam=Symbol('lambda', positive=True),
                 **kwargs):
        self.z = z
        self.Omega = Omega
        self.t = ivar
        self.M = M
        self.k_m = k_m
        self.c_m = c_m
        self.m_e = m_e
        self.e = e
        self.g = g
        self.phi = phi
        self.lam = lam
        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.SpringVer = Spring(2 * k_m, pos1=z, qs=[z])
        self.gravity_force1 = GravitationalForce(self.M, self.g, z, qs=[z])
        self.gravity_force2 = GravitationalForce(self.m_e,
                                                 self.g,
                                                 z + e * cos(phi),
                                                 qs=[z])
        self.damping = Damper((2 * c_m), pos1=z, qs=[z])
        system = self.SpringVer + self.MaterialPoint_1 + self.MaterialPoint_2 + self.gravity_force1 + self.gravity_force2 + self.damping
        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)

        default_data_dict = {
            self.c_m: [self.lam * (self.k_m)],
            self.M: [
                200 * m0, 350 * m0, 400 * m0, 550 * m0, 650 * m0, 700 * m0,
                800 * m0
            ],
            self.k_m: [
                2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,
                10 * k0
            ],
            self.m_e: [
                0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0, 0.7 * m0,
                0.8 * m0, 0.9 * m0
            ],
            self.e: [2 * e0, 3 * e0, 4 * e0, 5 * e0, 6 * e0],
            self.g: [g],

            #             self.phi:[self.Omega*self.t],
            self.phi: [self.Omega * self.t]
        }

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class InlineEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'inline_engine_perpendicular_springs.png'
    real_name = 'paccar.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l', positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.SpringVer = Spring(2 * k_m, pos1=z, qs=[z])
        self.SpringHor = Spring(2 * k_m, pos1=x, qs=[z])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class BoxerEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'boxer_engine_perpendicular_springs.png'
    real_name = 'f6c_valkyrie.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l', positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1 = MaterialPoint(M, pos1=x, qs=[x])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=x + e * sin(phi),
                                             qs=[x])
        self.SpringVer = Spring(2 * k_m, pos1=x, qs=[x])
        self.SpringHor = Spring(2 * k_m, pos1=z, qs=[x])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class NonLinearInlineEnginePerpendicularSpringsGravity(NonlinearComposedSystem
                                                       ):
    scheme_name = 'nonlinear_inline_engine_perpendicular_springs.png'
    real_name = 'paccar.jpg'
    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    l = Symbol('l', positive=True)
    g = Symbol('g', positive=True)
    z = dynamicsymbols('z')
    phi = dynamicsymbols('phi')
    Omega = Symbol('Omega', positive=True)
    ivar = Symbol('t')
    d = Symbol('d', positive=True)

    e0 = Symbol('e_0', positive=True)

    def __init__(
            self,
            d1=None,
            #d2=None,
            M=None,
            k_m=None,
            m_e=None,
            e=None,
            l=None,
            g=None,
            z=None,
            phi=None,
            Omega=None,
            ivar=None,
            d=None,
            **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if phi is not None: self.phi = phi
        if Omega is not None: self.Omega = Omega
        if d is not None: self.d = d
        if l is not None: self.l = l
        if g is not None: self.g = g
        if ivar is not None: self.ivar = ivar

        self.MaterialPoint_1 = MaterialPoint(self.M, pos1=self.z, qs=[self.z])
        self.MaterialPoint_2 = MaterialPoint(self.m_e,
                                             pos1=self.z +
                                             self.e * cos(self.phi),
                                             qs=[self.z])
        self.SpringVer = Spring(2 * self.k_m, pos1=self.z, qs=[self.z])
        self.SpringHor = Spring(2 * self.k_m,
                                pos1=(sqrt(self.d**2 + self.z**2) - self.l),
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
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1 + self.MaterialPoint_2 + self.gravity_force1 + self.gravity_force2
        super().__init__(system, **kwargs)

    def linearized(self):

        return type(self).from_system(super().linearized())

    def get_default_data(self):

        m0, k0, e0, l0 = symbols('m_0 k_0 e_0 l_0', positive=True)

        default_data_dict = {
            self.M: [
                200 * m0, 350 * m0, 400 * m0, 550 * m0, 650 * m0, 700 * m0,
                800 * m0
            ],
            self.k_m: [
                2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,
                10 * k0
            ],
            self.m_e: [
                0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0, 0.7 * m0,
                0.8 * m0, 0.9 * m0
            ],
            self.e: [2 * e0, 3 * e0, 4 * e0, 5 * e0, 6 * e0],
            self.d: [l0 * (S.One + no / 20) for no in range(5, 10)],
            self.l: [l0],
            #             self.g:[g],
            #             self.phi:[self.Omega*self.t],
            self.phi: [self.Omega * self.ivar]
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if 4 * parameters_dict[self.k_m] - 2 * parameters_dict[
                self.k_m] * parameters_dict[self.l] / parameters_dict[
                    self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0]) / 2

    def max_dynamic_force_pin(self):
        lin_sys = ComposedSystem(self.linearized())
        #k_m = self._given_data[self.k_m]
        k_m = self.k_m
        #         display(lin_sys.stiffness_matrix()[0])

        return lin_sys.frequency_response_function() * (
            lin_sys.stiffness_matrix()[0]) / 2 + self.max_static_force_pin()

    def max_dynamic_nonlinear_force_pin(self):
        lin_sys = ComposedSystem(self.linearized())

        amp = list(self.amplitude_from_frf())
        display(amp)
        #k_m = self._given_data[self.k_m]
        k_m = self.k_m

        return amp[0] * k_m + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


class NonLinearBoxerEnginePerpendicularSprings(NonlinearComposedSystem):
    scheme_name = 'nonlin_boxer_engine_perpendicular_springs.png'
    real_name = 'f6c_valkyrie.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l', positive=True),
                 x=dynamicsymbols('x'),
                 phi=dynamicsymbols('phi'),
                 Omega=Symbol('\Omega', positive=True),
                 ivar=Symbol('t'),
                 d=Symbol('d', positive=True),
                 **kwargs):
        self.t = ivar
        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e
        self.phi = phi
        self.Omega = Omega
        self.d = d
        self.l = l
        self.MaterialPoint_1 = MaterialPoint(M, pos1=x, qs=[x])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=x + e * sin(phi),
                                             qs=[x])
        self.SpringHor = Spring(2 * k_m, pos1=x, qs=[x])
        self.SpringVer = Spring(2 * k_m, pos1=(sqrt(d**2 + x**2) - l), qs=[x])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, k0, e0, l0 = symbols('m_0 k_0 e_0 l_0', positive=True)

        default_data_dict = {
            self.M: [
                200 * m0, 350 * m0, 400 * m0, 550 * m0, 650 * m0, 700 * m0,
                800 * m0
            ],
            self.k_m: [
                2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,
                10 * k0
            ],
            self.m_e: [
                0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0, 0.7 * m0,
                0.8 * m0, 0.9 * m0
            ],
            self.e: [2 * e0, 3 * e0, 4 * e0, 5 * e0, 6 * e0],
            self.d: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.l: [l0],
            #             self.g:[g],
            #             self.phi:[self.Omega*self.t],
            self.phi: [self.Omega * self.t]
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if 4 * parameters_dict[self.k_m] - 2 * parameters_dict[
                self.k_m] * parameters_dict[self.l] / parameters_dict[
                    self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force_pin(self):
        lin_sys = self.linearized()

        return lin_sys.frequency_response_function(
        ) * lin_sys.stiffness_matrix()[0] + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)

#Patryk
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
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    ivar=Symbol('t', positive=True)
    
    z=dynamicsymbols('z')
    
    def __init__(self,
                 m=None,
                 k=None,
                 ivar=None,
                 z=None,
                 **kwargs):

        
        
        if m is not None: self.m = m
        if k is not None: self.k = k
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        
   
        self.qs = [self.z]


    @property
    def components(self):

        components = {}
        
        self.material_point = MaterialPoint(self.m, self.z, self.qs)
        self.spring = Spring(self.k, self.z, self.qs)
        
        components['material_point'] = self.material_point
        components['spring'] = self.spring
        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict

#Amadi
class BeamBridge(ComposedSystem):
    """Ready to use model of bridge represented by the mass supported by elastic beam.
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
    scheme_name = 'beam_bridge.PNG'
    real_name = 'beam_bridge_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F=Symbol('F', positive=True),
                 l=Symbol('l', positive=True),
                 module=Symbol('E', positive=True),
                 inertia=Symbol('I', positive=True),
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.m = m

        self.k_beam = k_beam

        self.g = g
        self.Omega = Omega
        self.F = F
        self.l = l
        self.z = z
        self.module = module
        self.inertia = inertia

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(m, g, z)
        self.force = Force(-F * sin(Omega * ivar), pos1=z)

        composed_system = (self.mass + self.spring + self.gravity_force +
                           self.force)

        super().__init__(composed_system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        #         E0, I0, l0, m0, k0,c0, lam0= symbols('E_0 I_0 l_0 m_0 k_0 c_0 lambda_0', positive=True)
        E0, I0, l0, m0, lam0, F0 = symbols('E_0 I_0 l_0 m_0 lambda_0 F_0',
                                           positive=True)
        default_data_dict = {

            #             self.lam:[10],
            #             self.c:[self.k_beam*self.lam],
            self.k_beam: [S.One * 48 * self.module * self.inertia / self.l**3],
            self.m: [
                10 * m0, 20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0, 70 * m0,
                80 * m0, 90 * m0
            ],
            #             self.l:[l0,2*l0,3*l0,4*l0,5*l0,6*l0,7*l0,8*l0,9*l0],
            #             self.E:[E0,2*E0,3*E0,4*E0,5*E0,6*E0,7*E0,8*E0,9*E0],
            #             self.I:[I0,2*I0,3*I0,4*I0,5*I0,6*I0,7*I0,8*I0,9*I0],
            #,100* m0, 200 * m0, 300 * m0, 400 * m0, 500 * m0, 600 * m0, 700 * m0, 800 * m0, 900 * m0
            self.module: [
                E0,
                2 * E0,
                3 * E0,
                4 * E0,
                5 * E0,
                6 * E0,
                7 * E0,
                8 * E0,
                9 * E0,
                10 * E0,
                11 * E0,
                12 * E0,
                13 * E0,
                14 * E0,
                15 * E0,
                16 * E0,
                17 * E0,
                18 * E0,
                19 * E0,
            ],
            self.inertia: [
                I0,
                2 * I0,
                3 * I0,
                4 * I0,
                5 * I0,
                6 * I0,
                7 * I0,
                8 * I0,
                9 * I0,
                10 * I0,
                11 * I0,
                12 * I0,
                13 * I0,
                14 * I0,
                15 * I0,
                16 * I0,
                17 * I0,
                18 * I0,
                19 * I0,
            ],
            self.l: [
                l0,
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
            ],
            #             self.lam:[lam0,2*lam0,3*lam0,4*lam0,5*lam0,6*lam0,7*lam0,8*lam0,9*lam0],
            self.F: [
                F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0, 7 * F0, 8 * F0,
                9 * F0
            ]
        }
        return default_data_dict

#Amadi
class BeamBridgeDamped(ComposedSystem):

    scheme_name = 'bridge_dmp.png'
    real_name = 'beam_bridge_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F=Symbol('F', positive=True),
                 c=Symbol('c', positive=True),
                 l=Symbol('l', positive=True),
                 module=Symbol('E', positive=True),
                 inertia=Symbol('I', positive=True),
                 lam=Symbol('lambda', positive=True),
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.m = m
        self.c = c
        self.k_beam = k_beam
        self.lam = lam
        self.g = g
        self.Omega = Omega
        self.F = F
        self.l = l
        self.z = z
        self.module = module
        self.inertia = inertia
        #         c=self.lam*k_beamhaft

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(m, g, z)
        self.force = Force(-F * sin(Omega * ivar), pos1=z)
        self.damper = Damper(c, pos1=z, qs=[z])
        composed_system = (self.mass + self.spring + self.gravity_force +
                           self.force + self.damper)

        super().__init__(composed_system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        #         E0, I0, l0, m0, k0,c0, lam0= symbols('E_0 I_0 l_0 m_0 k_0 c_0 lambda_0', positive=True)
        E0, I0, l0, m0, lam0, F0 = symbols('E_0 I_0 l_0 m_0 lambda_0 F_0',
                                           positive=True)
        default_data_dict = {

            #             self.lam:[10],
            self.c: [self.k_beam * self.lam],
            self.k_beam: [S.One * 48 * self.module * self.inertia / self.l**3],
            self.m: [
                10 * m0, 20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0, 70 * m0,
                80 * m0, 90 * m0
            ],
            #             self.l:[l0,2*l0,3*l0,4*l0,5*l0,6*l0,7*l0,8*l0,9*l0],
            #             self.E:[E0,2*E0,3*E0,4*E0,5*E0,6*E0,7*E0,8*E0,9*E0],
            #             self.I:[I0,2*I0,3*I0,4*I0,5*I0,6*I0,7*I0,8*I0,9*I0],
            #,100* m0, 200 * m0, 300 * m0, 400 * m0, 500 * m0, 600 * m0, 700 * m0, 800 * m0, 900 * m0
            self.module: [
                E0,
                2 * E0,
                3 * E0,
                4 * E0,
                5 * E0,
                6 * E0,
                7 * E0,
                8 * E0,
                9 * E0,
                10 * E0,
                11 * E0,
                12 * E0,
                13 * E0,
                14 * E0,
                15 * E0,
                16 * E0,
                17 * E0,
                18 * E0,
                19 * E0,
            ],
            self.inertia: [
                I0,
                2 * I0,
                3 * I0,
                4 * I0,
                5 * I0,
                6 * I0,
                7 * I0,
                8 * I0,
                9 * I0,
                10 * I0,
                11 * I0,
                12 * I0,
                13 * I0,
                14 * I0,
                15 * I0,
                16 * I0,
                17 * I0,
                18 * I0,
                19 * I0,
            ],
            self.l: [
                l0,
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
            ],
            self.lam: [
                lam0, 2 * lam0, 3 * lam0, 4 * lam0, 5 * lam0, 6 * lam0,
                7 * lam0, 8 * lam0, 9 * lam0
            ],
            self.F: [
                F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0, 7 * F0, 8 * F0,
                9 * F0
            ]
        }

        return default_data_dict

#Laura
class DampedSpringMassSystem(ComposedSystem):

    scheme_name = '???'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 c=Symbol('c', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.m = m
        self.k = k
        self.c = c

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k, z, qs=[z])
        self.damper = Damper(c, pos1=z, qs=[z])
        system = self.mass + self.spring + self.damper

        super().__init__(system, **kwargs)


class DampedHarmonicOscillator(DampedSpringMassSystem):
    pass

#Patryk
class Pendulum(NonlinearComposedSystem):
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
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """
    scheme_name = 'undamped_pendulum.png'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('\\varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l

        self.potential = GravitationalForce(self.m,
                                            self.g,
                                            l * (1 - cos(angle)),
                                            qs=qs)
        self.kinen = MaterialPoint(self.m * self.l**2, pos1=qs[0], qs=[angle])
        print(self.kinen)
        system = self.potential + self.kinen
        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            #         mech_comp.LinearizationComponent,
            #         mech_comp.GoverningEquationComponent,
            #         mech_comp.FundamentalMatrixComponent,
            #         mech_comp.GeneralSolutionComponent,
            #         mech_comp.SteadySolutionComponent,
        ]

        return comp_list


class PulledPendulum(ComposedSystem):
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
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """
    scheme_name = 'undamped_pendulum.png'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('\\varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l

        self.potential = GravitationalForce(self.m,
                                            self.g,
                                            l * (1 - cos(angle)),
                                            qs=qs)
        self.kinen = MaterialPoint(self.m * self.l**2, pos1=qs[0], qs=[angle])
        self.force = Force(-2 * self.m * self.l * (self.g / self.l * cos(pi)),
                           angle)
        print(self.kinen)
        system = self.potential + self.kinen + self.force
        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            #         mech_comp.LinearizationComponent,
            mech_comp.GoverningEquationComponent,
            mech_comp.FundamentalMatrixComponent,
            mech_comp.GeneralSolutionComponent,
            mech_comp.SteadySolutionComponent,
        ]
        return comp_list


# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie

# Sav
class FreePendulum(ComposedSystem):
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
        >>> SDOFFreePendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SFoDFreePendulum()
    """
    scheme_name = 'free_sdof_pendulum.png'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.g = g
        self.l = l

        self.pendulum = Pendulum(m, g, l, angle=angle)
        system = self.pendulum

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m:
            [2 * m0, 1 * m0, S.Half * m0, S.Half**2 * m0, 3 * S.Half * m0],
            self.l:
            [2 * l0, 1 * l0, S.Half * l0, S.Half**2 * l0, 3 * S.Half * l0],
        }
        return default_data_dict


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
    scheme_name = 'damped_excited_pendulum.PNG'
    real_name = 'pendulum2_real.jpg'

    def __init__(self,
                 dummy=Symbol('dummy', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 F=Symbol('F', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        phi = angle
        self.phi = phi

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.F = F

        Omega = Symbol('Omega', positive=True)
        self.Omega = Omega
        self.pendulum = Pendulum(m, g, l, angle=phi)
        self.force = Force(-F * l * sin(Omega * ivar), pos1=phi, qs=qs)
        system = self.pendulum + self.force

        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, l0, g0, F0 = symbols('m_0 l_0 g_0 F_0', positive=True)

        default_data_dict = {
            self.m:
            [2 * m0, 1 * m0, S.Half * m0, S.Half**2 * m0, 3 * S.Half * m0],
            self.l:
            [2 * l0, 1 * l0, S.Half * l0, S.Half**2 * l0, 3 * S.Half * l0],
            self.g: [g0],
            self.F:
            [2 * F0, 1 * F0, S.Half * F0, S.Half**2 * F0, 3 * S.Half * F0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.F: r'Force',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0 = symbols('m_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.m1: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
            self.F: [
                1 * F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0, 7 * F0, 8 * F0,
                9 * F0, 10 * F0, 11 * F0, 12 * F0, 13 * F0, 14 * F0, 15 * F0,
                16 * F0, 17 * F0, 18 * F0, 19 * F0, 20 * F0, 21 * F0, 22 * F0,
                23 * F0, 24 * F0, 25 * F0, 26 * F0, 27 * F0, 28 * F0, 29 * F0,
                30 * F0
            ],
        }
        return default_data_dict


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
    scheme_name = 'damped_pendulum.png'
    real_name = 'pendulum2_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 c=Symbol('c', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        phi = angle
        self.phi = phi

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.c = c

        self.Pendulum = Pendulum(m, g, l, angle=phi)
        self.Damper = Damper(c, l * phi, qs=qs)
        system = self.Pendulum + self.Damper

        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, l0, c0 = symbols('m_0 l_0 c_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.l: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.c: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0]
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict


class ExcitedDampedPendulum(ComposedSystem):

    scheme_name = 'damped_excited_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 c=Symbol('c', positive=True),
                 F=Symbol('F', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        phi = angle

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.c = c
        self.F = F
        self.Omega = Omega

        self.Pendulum = Pendulum(m, g, l, angle=phi)
        self.Damper = Damper(c, l * phi, qs=qs)
        self.Force = Force(-F * sin(Omega * ivar), pos1=phi, qs=[phi])
        system = self.Pendulum + self.Damper + self.Force

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict


class PendulumKinematicExct(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    l = Symbol('l', positive=True)
    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    phi = dynamicsymbols('\\varphi')
    x_e = dynamicsymbols('x_e')

    def __init__(self,
                 l0=None,
                 l=None,
                 m=None,
                 g=None,
                 phi=None,
                 x_e=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if l is not None: self.l = l
        if m is not None: self.m = m
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x_e is not None: self.x_e = x_e

        self.ivar = ivar
        self.qs = [self.phi]

        x = self.l * sin(self.phi) + self.x_e
        y = self.l * cos(self.phi)

        self.material_point_1 = MaterialPoint(self.m, x, qs=self.qs)
        self.material_point_2 = MaterialPoint(self.m, y, qs=self.qs)
        self.gravity = GravitationalForce(self.m, self.g, pos1=-y, qs=self.qs)

        system = self.material_point_1 + self.material_point_2 + self.gravity

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.x_e: r'Kinematic lateral excitation',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, x0, Omega = symbols('m_0 l_0 x_0 Omega', positive=True)

        default_data_dict = {
            self.m: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
            self.x_e: [x0 * sin(self.Omega * self.ivar)]
        }
        return default_data_dict

    def max_static_cable_force(self):
        return (self.m * self.g).subs(self._given_data)

    def max_dynamic_cable_force(self):

        omg_amp = ComposedSystem(
            self.linearized()).frequency_response_function() * self.Omega

        return (self.m * self.l * (omg_amp)**2 + self.max_static_cable_force())

    def static_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re))**(1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re))**(1 / 2)


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


class Engine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'
    """
    Model of a SDoF engine.

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

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A SDoF engine oscillating up and down while being held up by a two springs with the given stiffness

        >>> t = symbols('t')
        >>> M, me, e, km = symbols('M, m_e, e, k_m')
        >>  phi = dynamicsymbosl('phi')
        >>> qs = dynamicsymbols('z') #Generalized coordinate
        >>> Engine()

    """

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.Spring = Spring(2 * k_m, pos1=z, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

#Pati
class DampedEngine(ComposedSystem):
    scheme_name = 'engine_with_damper.png'
    real_name = 'engine_real.PNG'
    """
    Model of a SDoF engine.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            me = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            k_m = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A SDoF engine oscillating up and down while being held up by a two springs with the given stiffness

        >>> t = symbols('t')
        >>> M, me, e, km = symbols('M, m_e, e, k_m')
        >>  phi = dynamicsymbosl('phi')
        >>> qs = dynamicsymbols('z') #Generalized coordinate
        >>> SDoFEngine()

    """

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 c_m=Symbol('c_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.c_m = c_m
        self.m_e = m_e
        self.e = e
        self.phi = phi

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.Spring = Spring(2 * k_m, pos1=z, qs=[z])
        self.damper = Damper(2 * c_m, pos1=z, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2 + self.damper
        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'unbalanced rotating mass',
            self.e: r'radius of rotation',
        }
        return self.sym_desc_dict


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

#Tomek
class ForcedNonLinearDisc(NonlinearComposedSystem):
    scheme_name = 'nonlinear_disc.png'
    real_name = 'roller_tightener.png'

    def __init__(self,
                 m1=Symbol('m', positive=True),
                 kl=Symbol('k', positive=True),
                 R=Symbol('R', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 F=Symbol('F', positive=True),
                 Omega=Symbol('Omega', positive=True),
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
        self.F = F

        self.disk1 = MaterialPoint(m1, x, qs=[x]) + MaterialPoint(
            m1 / 2 * R**2, x / R, qs=[x]) + Spring(
                kl, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])
        self.force = Force(F * cos(Omega * ivar), pos1=x, qs=[x])

        system = self.disk1 + self.force
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
