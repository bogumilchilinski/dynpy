import importlib
from . import dynamics as dyn
importlib.reload(dyn)
from . import dynamics as dyn

from sympy.physics.vector import dynamicsymbols
from sympy import *
from sympy.physics.mechanics import *
import sympy as sym
import numpy as np
import pylatex
from pylatex import Command, NoEscape
import pint
from dynpy.dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp
from dynpy.models.elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from dynpy.models.mechanics.trolley import SpringMassSystem, SpringDamperMassSystem
# from dynpy.models.sdof import ComposedSystem, SpringMassSystem, DampedSpringMassSystem, BeamBridgeDamped
ureg = pint.UnitRegistry()

mechanics_printing()

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

t=Symbol('t') #independent variable - time

class LEWEL16_simple(SpringMassSystem):
    scheme_name = 'train_sdof.png'
    real_name = 'lewel16.jpg'

class LEWEL16_simple_damped(SpringDamperMassSystem):
    scheme_name = 'train_sdof_damped.png'
    real_name = 'lewel16.jpg'

class LeafSpring(ComposedSystem):

    m=Symbol('m', positive=True),
    k_beam=Symbol('k_beam', positive=True)
    ivar=Symbol('t')
    g=Symbol('g', positive=True)
    Omega=Symbol('Omega', positive=True)
    F_0=Symbol('F_0', positive=True)
    c=Symbol('c', positive=True)
    l=Symbol('l', positive=True)
    module=Symbol('E', positive=True)
    inertia=Symbol('I', positive=True)
    lam=Symbol('lambda',positive=True)
    z=dynamicsymbols('z')
    
    scheme_name = 'bridge_dmp.png'
    real_name = 'leaf_spring.jpg'
    
    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=None,
                 g=None,
                 Omega=None,
                 F_0=None,
                 c=None,
                 l=None,
                 module=None,
                 inertia=None,
                 lam=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if c is not None: self.c=c
        if k_beam is not None: self.k_beam = k_beam
        if lam is not None: self.lam=lam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F_0 is not None: self.F_0 = F_0
        if l is not None: self.l=l
        if z is not None: self.z = z
        if module is not None: self.module=module
        if inertia is not None: self.inertia=inertia
#         c=self.lam*k_beamhaft
        
        self.qs=[self.z]
        self._init_from_components()


        
#         self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)

        composed_system = (self.mass + self.spring + self.damper)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z, qs=[self.z])(label = 'Material point')
        self._spring = Spring(self.k_beam, self.z, qs=[self.z])(label = 'Structural stiffness')
        self._gravitational_force = GravitationalForce(self.m, self.g, self.z)(label = 'Gravitational force')
        self._damper = Damper(self.c, pos1=self.z, qs=[self.z])(label = 'Internal Rayleigh damping')
        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravitational_force'] = self._gravitational_force
        components['_damper'] = self._damper
        return components


    def set_equivalent_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        self.b=b
        self.h=h
        self.rho=rho
        default_data_dict = {

#             self.lam:[10],
            self.c:self.k_beam*self.lam,
            self.k_beam: S.One*48*self.module * self.inertia / self.l**3,
            self.inertia: b*h**3/12,
#             self.m:b*h*self.l*rho
        }

        return default_data_dict
    
    def get_default_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        default_data_dict = {self.m:3000,self.module:2.1*10**11,self.rho:7800,self.b:0.08,self.h:9*0.012,self.l:1.5,self.lam:0.1
        }

        return default_data_dict
    def get_real_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        default_data_dict = {self.m:3000,self.module:2.1*10**11,self.rho:7800,self.b:0.08,self.h:8*0.0135,self.l:0.9,self.lam:0.1
        }

        return default_data_dict