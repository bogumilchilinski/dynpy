from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from  ..continuous import ContinuousSystem, PlaneStressProblem

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


    def _init_from_components(self,*args,system=None,**kwargs):
        
        if system is None:
            composed_system = self._elements_sum
        else:
            composed_system = system
            
        #print('CS',composed_system._components)
        super().__init__(None,system = composed_system)
        
        #print('self',self._components)
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

#     def tensioner_belt_force(self):
#         return self.k_tensioner * self.steady_solution()

#     def left_belt_force(self):
#         return self.k_belt * self.steady_solution()

#     def right_belt_force(self):
#         return self.k_belt * self.steady_solution()

#     def max_static_force_pin(self):
#         return abs(self.static_load().doit()[0])

#     def max_dynamic_force_pin(self):
#         return self.frequency_response_function() * self.stiffness_matrix(
#         )[0] + self.max_static_force_pin()

#     def static_force_pin_diameter(self):
#         kt = Symbol('k_t', positive=True)
#         Re = Symbol('R_e', positive=True)
#         return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

#     def dynamic_force_pin_diameter(self):
#         kt = Symbol('k_t', positive=True)
#         Re = Symbol('R_e', positive=True)
#         return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


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

    
    
    
class SDOFSemiSubmersiblePlatform(ComposedSystem):

    #scheme_name = '.PNG'
    #real_name = '.jpg'
    
    z=dynamicsymbols('z')
    m=Symbol('m', positive=True)
    m_p=Symbol('m_p', positive=True)
    k_1=Symbol('k', positive=True)
    c_1=Symbol('c', positive=True)
    c_2=Symbol('c_water', positive=True)
    F_load=Symbol('F_{load}', positive=True)
    g=Symbol('g', positive=True)
    p_p=Symbol('P_p', positive=True)
    rho=Symbol('\\rho', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    
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
        
        if z is not None: self.z=z
        if m is not None: self.m = m
        if k_1 is not None: self.k_1 = k_1
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if F_load is not None: self.F_load = F_load
        if g is not None: self.g=g
        if rho is not None: self.rho=rho
        if p_p is not None: self.p_p=p_p
        if Omega is not None: self.Omega=Omega
        if F is not None: self.F = F
        self.ivar=ivar
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}
        
        self._platform = MaterialPoint(self.m, pos1=self.z, qs=[self.z])(label='Material point - mass of the platform')
        self._spring_1 = Spring(self.k_1, pos1=self.z, pos2=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Spring - stiffness of the platform')
        self._spring_2 = Spring(self.p_p*self.rho, pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Spring - buoyancy')
        self._damper = Damper(self.c_1, pos1=self.z, pos2=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Damper - damping of the platform')
        self._water = Damper(self.c_2, pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Damper - damping of the water')
        self._force = Force(-self.F_load, pos1=self.z, qs=[self.z])(label='Force - load on the platform')
        self._MassOfPlatform = Force(-self.m*self.g, pos1=self.z, qs=[self.z])(label='Force - weight of the platform in gravitational field')
        self._excitation = Force(self.F*sin(self.Omega*self.ivar), pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Force - Force caused by the waves')
        
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

        m0, Omega0, k0, p_p0, F0, lam0 = Symbol('m_0', positive=True), Symbol('Omega0', positive=True), Symbol('k0', positive=True), Symbol('p_{p0}', positive=True), Symbol('F0', positive=True), Symbol('lambda0', positive=True)
        c_2 = self.c_2
        g=self.g
        rho=self.rho

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

    
    
class DDOFSemiSubmersiblePlatform(ComposedSystem):

    #scheme_name = '.PNG'
    #real_name = '.jpg'
    
    z_p=dynamicsymbols('z_p')
    z=dynamicsymbols('z')
    m=Symbol('m', positive=True)
    m_p=Symbol('m_p', positive=True)
    k_1=Symbol('k', positive=True)
    c_1=Symbol('c', positive=True)
    c_2=Symbol('c_p', positive=True)
    F_load=Symbol('F_{load}', positive=True)
    g=Symbol('g', positive=True)
    p_p=Symbol('P_p', positive=True)
    rho=Symbol('rho', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    
    def __init__(self,
                 m=None,
                 m_p=None,
                 k_1=None,
                 c_1=None,
                 c_2=None,
                 F_load=None,
                 ivar=Symbol('t'),
                 z=None,
                 z_p=None,
                 g=None,
                 p_p=None,
                 rho=None,
                 Omega=None,
                 F=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if z_p is not None: self.z_p=z_p
        if m is not None: self.m = m
        if m_p is not None: self.m_p = m_p
        if k_1 is not None: self.k_1 = k_1
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if F_load is not None: self.F_load = F_load
        if g is not None: self.g=g
        if rho is not None: self.rho=rho
        if p_p is not None: self.p_p=p_p
        if Omega is not None: self.Omega=Omega
        if F is not None: self.F=F
        self.ivar=ivar
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}
        
        self._platform = MaterialPoint(self.m, pos1=self.z, qs=[self.z])(label='Material point - mass of the platform')
        self._pontoon = MaterialPoint(self.m_p, pos1=self.z_p, qs=[self.z_p])(label='Material point - mass of the pontoon')
        self._spring_1 = Spring(self.k_1, pos1=self.z, pos2=self.z_p, qs=[self.z, self.z_p])(label='Spring - stiffness of the platform')
        self._spring_2 = Spring(self.p_p*self.rho, pos1=self.z_p, qs=[self.z_p, self.z])(label='Spring - buoyancy')
        self._damper = Damper(self.c_1, pos1=self.z, pos2=self.z_p, qs=[self.z, self.z_p])(label='Damper - damping of the platform')
        self._water = Damper(self.c_2, pos1=self.z_p, qs=[self.z_p, self.z])(label='Damper - damping of the water')
        self._force = Force(-self.F_load, pos1=self.z, qs=[self.z, self.z_p])(label='Force - load on the platform')
        self._MassOfPlatform = Force(-self.m*self.g, pos1=self.z, qs=[self.z, self.z_p])(label='Force - weight of the platform in gravitational field')
        self._MassOfPontoon = Force(-self.m_p*self.g, pos1=self.z_p, qs=[self.z_p, self.z])(label='Force - weight of the pontoon in gravitational field')
        self._excitation = Force(self.F*sin(self.Omega*self.ivar), pos1=self.z_p, qs=[self.z_p, self.z])(label='Force - Force caused by the waves')
        
        
        components['_platform'] = self._platform
        components['_spring_1'] = self._spring_1
        components['_spring_2'] = self._spring_2
        components['_damper'] = self._damper
        components['_water'] = self._water
        components['_force'] = self._force
        components['_MassOfPlatform'] = self._MassOfPlatform
        components['_MassOfPontoon'] = self._MassOfPontoon
        components['_pontoon'] = self._pontoon
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
            self.ivar: r'Time',
            self.m_p: r'Mass of a pontoon'
        }
        return self.sym_desc_dict
    
    
    def get_default_data(self):

        m0, Omega0, k0, p_p0, F0, lam0 = Symbol('m_0', positive=True), Symbol('Omega0', positive=True), Symbol('k0', positive=True), Symbol('p_{p0}', positive=True), Symbol('F0', positive=True), Symbol('lambda_0', positive=True)
        c_2 = self.c_2
        g=self.g
        rho=self.rho

        default_data_dict = {
            self.m: [S.One * m0 * no * 1000 for no in range(20, 30)],
            self.m_p: [S.One * m0 * no * 1000 for no in range(10, 20)],
            self.g: [S.One * g * no for no in range(1, 2)],
            self.Omega: [S.One * Omega0 * no for no in range(1, 2)],
            self.c_2: [S.One * c_2 * no for no in range(1, 2)],
            self.rho: [S.One * rho * no for no in range(1, 2)],
            self.k_1: [S.One * k0 * no * 10000 for no in range(5, 10)],
            self.c_1: [S.One * k0 * lam0 * no * 10000 for no in range(5, 10)],
            self.p_p: [S.One * p_p0 * no for no in range(9, 16)],
            self.F_load: [S.One * F0 * 1000 * no  for no in range(5, 10)],
            self.F: [S.One * F0 * 1000 * no for no in range(1, 10)]
        }

        return default_data_dict
