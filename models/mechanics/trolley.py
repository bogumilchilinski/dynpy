from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from ...continuous import ContinuousSystem, PlaneStressProblem

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
    ivar=Symbol('t')

    
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
        super(HarmonicOscillator,self).__init__(None, system=composed_system)

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
                 system=None,
                 **kwargs):

        if ivar is not None: self.ivar = ivar
        if m0 is not None: self.m0 = m0

        if qs is not None:
            self.qs = qs
        else:
            self.qs = [self.z]

        self._init_from_components(system=system, **kwargs)

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

    def _components_default_data(self):
        
        data=[elem._all_default_data()   for elem in self.elements.values()]

        
        return {key:value for elem in data for key, value in elem.items()}    
    
    def _components_numerical_data(self):
        
        data=[elem._all_numerical_data()   for elem in self.elements.values()]
        
        
        return {key:value for elem in data for key, value in elem.items()}    
    
    def _all_default_data(self):
        
        

        
        return {**self._components_default_data(),**self.get_default_data()}    
    
    def _all_numerical_data(self):
        
        return {**self._components_numerical_data(),**self.get_numerical_data()}  
    
    
    def get_default_data(self):
        return {}

    def get_numerical_data(self):
        return {}

    
    
    
    def get_random_parameters(self):

        
        #print('preview for',self)
        #display(self._all_default_data())
        #display(self.get_default_data())
        
        default_data_dict = {**self._components_default_data(),**self.get_default_data()}

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    def get_numerical_parameters(self):

        default_data_dict = {**self._components_numerical_data(),**self.get_numerical_data()}

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



    def tensioner_belt_force(self):
        return self.k_tensioner * self.steady_solution()

    def left_belt_force(self):
        return self.k_belt * self.steady_solution()

    def right_belt_force(self):
        return self.k_belt * self.steady_solution()


#     def max_static_force_pin(self):
#         return abs(self.static_load().doit()[0])

#     def max_dynamic_force_pin(self):
#         return self.frequency_response_function() * self.stiffness_matrix(
#         )[0] + self.max_static_force_pin()

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0]) / 2

    def max_dynamic_force_pin(self):
        return self._frf()[0] * self.k_m + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)
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
    

    

#Patryk
#dane domyślne i numeryczne
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
                 z=None,
                 ivar=None,
                 **kwargs):

        
        
        if m is not None: self.m = m
        if k is not None: self.k = k
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        
   
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring = Spring(self.k, self.z, qs=self.qs)
        
        components['material_point'] = self.material_point
        components['spring'] = self.spring
        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict

#dane domyślne i numeryczne
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
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'


    F=Symbol('F', positive=True)
    g=Symbol('g', positive=True)

    
    def __init__(self,
                 m=None,
                 k=None,
                 F=None,
                 g=None,
                 z=None,
                 ivar=None,
                 **kwargs):


        if m is not None: self.m = m
        if k is not None: self.k = k
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if F is not None: self.F = F
        if g is not None: self.g = g


        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}  # przerobić na kompozycję
       
        
        self._undamped_trolley = SpringMassSystem(self.m, self.k, self.z)
        self.force = Force(self.F, self.z, qs=self.qs)
        self.gravitational_force = GravitationalForce(self.m, self.g, self.z, qs = self.qs)
        

        components['_undamped_trolley'] = self._undamped_trolley
        components['force'] = self.force
        components['gravitational_force'] = self.gravitational_force
        
        return components


    def symbols_description(self):

        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return {**super().symbols_description(),**self.sym_desc_dict}

#dziedziczenie po SpringMass
#dane domyślne i numeryczne
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
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    F=Symbol('F', positive=True)
    c=Symbol('c',positive=True)
    ivar=Symbol('t', positive=True)
    
    
    z=dynamicsymbols('z')
    
    def __init__(self,
                 m=None,
                 k=None,
                 F=None,
                 z=None,
                 c=None,
                 ivar=None,
                 **kwargs):

        
        
        if m is not None: self.m = m
        if k is not None: self.k = k
        if F is not None: self.F = F
        if z is not None: self.z = z
        if c is not None: self.c= c
        if ivar is not None: self.ivar = ivar

        
   
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring = Spring(self.k, self.z, qs=self.qs)
        self.force = Force(self.F, self.z, qs=self.qs)
        self.damper= Damper(self.c, self.z, qs=self.qs)
        
        components['material_point'] = self.material_point
        components['spring'] = self.spring
        components['force'] = self.force
        components['damper']=self.damper
        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict

    #można spróbować dziedizczyć bo SPringMass
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
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    k_l=Symbol('k_l', positive=True)
    k_r=Symbol('k_r', positive=True)
    ivar=Symbol('t', positive=True)
    x=dynamicsymbols('x')

    def __init__(self,
                 m=None,
                 k_l=None,
                 k_r=None,
                 x=None,
                 ivar=None,
                 **kwargs):



        if m is not None: self.m = m
        if k_l is not None: self.k_l = k_l
        if k_r is not None: self.k_r = k_r
        if ivar is not None: self.ivar = ivar
        if x is not None: self.x = x


        self.qs = [self.x]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, pos1=self.x, qs=self.qs)('Trolley')
        self.spring_l = Spring(self.k_l, pos1=self.x, qs=self.qs)('Left joint')
        self.spring_r = Spring(self.k_r, pos1=self.x, qs=self.qs)('Right joint')

        components['material_point'] = self.material_point
        components['spring_l'] = self.spring_l
        components['spring_r'] = self.spring_r
        
        return components

    
    def get_default_data(self):

        m0, k0 = self.m0, self.k0
        
        default_data_dict = {

            self.m: [S.One * no * m0 /100 for no in range(80, 120)], # percentage
            self.k_r: [S.One * no * k0/100 for no in range(80, 120)], # percentage
            self.k_l: [S.One * no * k0/100 for no in range(70, 110)], # percentage
        }
        return default_data_dict

    def get_numerical_data(self):

       
        m0, k0 = 100, 1000
        
        
        default_data_dict = {
            self.m: [S.One * no * m0 /100 for no in range(80, 120)], # percentage
            self.k_r: [S.One * no * k0/100 for no in range(80, 120)], # percentage
            self.k_l: [S.One * no * k0/100 for no in range(70, 110)], # percentage

        }
        return default_data_dict

    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_l: r'Spring left coefficient ', # ewentualnie zmienić na spring 1,2 coefficient
            self.k_r: r'Spring right coefficient ',
            self.x: r'displacement',
        }

        return self.sym_desc_dict

