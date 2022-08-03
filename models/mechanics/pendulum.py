from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin


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
    #scheme_name = 'undamped_pendulum.png'
    #real_name = 'pendulum_real.jpg'

    
    m=Symbol('m', positive=True)
    g=Symbol('g', positive=True)
    l=Symbol('l', positive=True)
    angle=dynamicsymbols('\\varphi')
    qs=None
    ivar=Symbol('t', positive=True)
    
    def __init__(self,
                 m=None,
                 g=None,
                 l=None,
                 angle=None,
                 ivar=None,
                 **kwargs):

        #if qs == None:
        #    qs = [angle]
        #else:
        #    qs = qs

        if m is not None: self.m = m
        if g is not None: self.g = g
        if l is not None: self.l = g
        if angle is not None: self.angle = angle
        if ivar is not None: self.ivar = ivar
        if t is not None: self.t = t
        
        self.qs = [self.angle]

        self._init_from_components(**kwargs)
        
        
        
    @property
    def components(self):

        components = {}
        
        self.gravitationalforce = GravitationalForce(self.m, self.g, l * (1 - cos(angle)), qs = self.qs)
        self.material_point = MaterialPoint(self.m * self.l**2, pos1=self.angle, qs=self.qs)
        
        components['material_point'] = self.material_point
        components['gravitationalforce'] = self.gravitationalforce
        
        return components
        
    def get_default_data(self):

       
        m0, l0 = self.m0, self.l0
        
        
        default_data_dict = {

            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            
        }
        return default_data_dict

    def get_numerical_data(self):

       
        m0, l0 = self.m0, self.l0
        
        
        default_data_dict = {

            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            
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
#Sav
class FreePendulum(ComposedSystem):
    
    scheme_name = 'free_sdof_pendulum.png'
    real_name = 'pendulum_real.jpg'
    
    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    l=Symbol('l',positive=True)
    angle=dynamicsymbols('\\varphi') # musi być zgodność
    ivar=Symbol('t')


    def __init__(self,
                 m=None,
                 g=None,
                 l=None,
                 angle=dynamicsymbols('\\varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if l is not None: self.l = l
        if angle is not None: self.angle=angle
        if qs is not None: self.qs=qs
        if ivar is not None: self.ivar=ivar
    
        self.qs = [angle]

        self._init_from_components(**kwargs)
        

    @property
    def components(self):

        components = {}

        x = self.l * sin(self.angle)
        y = self.l * cos(self.angle)
        

        self._mass_x = MaterialPoint(self.m,
                                     pos1=x,
                                     qs=self.qs)(label='Payload - x component')
        
        self._mass_y = MaterialPoint(self.m,
                                     pos1=y,
                                     qs=self.qs)(label='Payload - y component')

        self._gravity = GravitationalForce(self.m,
                                            self.g,
                                            pos1=-y,
                                            qs=self.qs)(label='Gravity field')


        components['_mass_x'] = self._mass_x #ok
        components['_mass_y'] = self._mass_y #ok
        components['_gravity'] = self._gravity #ok
        return components
        

    def symbols_description(self):
        
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.angle: r'Angle',
            self.ivar: r'Time'
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.l: [l0 * no for no in range(1, 8)],
        }

        return default_data_dict

