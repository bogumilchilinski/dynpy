from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset, atan)

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

from .pendulum import Pendulum, PendulumKinematicExct
from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin
from .tmd import TunedMassDamperRelativeMotion
from ...utilities.components.mech import en as mech_comp

from pint import UnitRegistry
ureg = UnitRegistry()

from functools import cached_property, lru_cache



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
    ivar=Symbol('t')
    
    z=dynamicsymbols('z')
    
    def __init__(self,
                 m=None,
                 k=None,
                 z=None,
                 ivar=Symbol('t'),
                 **kwargs):

        
        
        if m is not None: self.m = m
        if k is not None: self.k = k
        if z is not None: self.z = z
        self.ivar = ivar
   
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
    ivar=Symbol('t')
    
    
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
    ivar=Symbol('t')
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

    

class NonLinearTrolley(NonlinearComposedSystem):
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    d=Symbol('d', positive=True)
    l_0=Symbol('l_0', positive=True)
    x=dynamicsymbols('x')
   
    def __init__(self,
                 m=None,
                 k=None,
                 d=None,
                 l_0=None,
                 x=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k= k
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if x is not None: self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        #self._disk = Disk(S.One/2 * self.m1*self.R**2,pos1 = self.x / self.R, qs=[self.x]) ??? returns diffrend eoms than soultion above
            
        self._spring = Spring(self.k, pos1 = sqrt(self.x**2 + self.d**2), pos2 = - self.l_0 , qs=[self.x])


        components['trolley'] = self._trolley
        components['spring'] = self._spring
        
        return components

    def get_default_data(self):

        m0, k0, l0= symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                0.5 * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
            self.d: [
                5 * l0, 2 * l0, 3 * S.Half * l0, 4 * l0, 6 * l0, 7 * l0,
                8 * l0, 9 * l0
            ],
            self.k: [
                1 * k0, 3 * k0, 2 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0,
                9 * k0
            ],
            self.l_0: [
                1 * l0, 3 * l0, 2 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0,
            ],
        }

        return default_data_dict

    def get_numerical_data(self): ### Brakowało numerical data. Przez to komenda get numerical parameters nie działa

        m0, k0, l0, c0= symbols('m_0 k_0 l_0 c_0', positive=True)

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
            self.l_0: [
                1 * l0, 3 * l0, 2 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0,
            ], ### Brakowało nawiasu
            self.c:  [
                1 * c0, 3 * c0, 2 * c0, 4 * c0, 5 * c0, 6 * c0, 7 * c0, 8 * c0,
                9 * c0
            ],
        }

        return default_data_dict



#Zrobione Amadi
class TrolleyWithPendulum(ComposedSystem):

    scheme_name = 'trolley_pendulum_tmd.png'
    real_name = 'taipei101.png'

    l = Symbol('l', positive=True)
    m_t = Symbol('m_trolley', positive=True)
    m_p = Symbol('m_pendulum', positive=True)
    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')

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
        self.ivar = ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(label='Trolley')
        self._pendulum = PendulumKinematicExct(self.l, self.m_p, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._force=Force(self.F*sin(self.Omega*self.ivar), pos1=self.x, qs=[self.x, self.phi])(label='Force')

        components['_trolley'] = self._trolley
        components['_pendulum'] = self._pendulum
        components['_force'] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.k: r'Stiffness of a beam showed as a spring stiffness in trolley member',
            self.x: r'Kinematic lateral excitation',
            self.phi: r'Angle of a pendulum',
            self.m_p: r'Mass of pendulum',
            self.m_t: r'Mass of trolley',
            self.g: r'Gravity constant',
            self.F: r'Force',
            self.Omega: r'Excitation frequency',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols('m_0 l_0 F_0 Omega_0 k_0', positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * self.Omega],
            self.g: [S.One * self.g],
            self.k: [S.One * self.g * m0 / l0 * no for no in range(1, 30)],
            self.x: [self.x]
        }
        
        #default_data_dict.update({self.k: S.One * self.g * default_data_dict[self.m_t][0]+default_data_dict[self.l][0],
        #                         })
        
        return default_data_dict
    
    def get_numerical_data(self):

        m_taipei = 700 * 1e3
        omg_nat = 2*pi/10
        
        default_data_dict = {
            self.m_t: [m_taipei for no in range(20, 30)],
            self.m_p: [0.01*m_taipei for no in range(1, 10)],
            self.l: [9.81/omg_nat**2 for no in range(1, 10)],           
            self.F: [300 * 1e3  for no in range(50, 100)],
            self.Omega: [omg_nat   for no in range(1,6)],
            self.k: [m_taipei*(omg_nat)**2 for no in range(50, 100)],
            self.g: [9.81]
        }        
        return default_data_dict

    @lru_cache
    def max_static_cable_force(self):
        data=self._given_data
        ans=self.force_in_cable()
        free_coeff=ans.subs({(cos(self.Omega*self.ivar))**2:0, (sin(self.Omega*self.ivar))**2:0}).subs(data)
        #display(free_coeff)
        return (abs(free_coeff)).subs(data)

    @lru_cache
    def max_dynamic_cable_force(self):
        ans = self.force_in_cable()
        data=self._given_data
        sin_coeff=(ans.coeff((sin(self.Omega*self.ivar))**2)).subs(data)
        #display(sin_coeff)
        cos_coeff=(ans.coeff((cos(self.Omega*self.ivar))**2)).subs(data)
        #display(cos_coeff)
        free_coeff=(ans.subs({(cos(self.Omega*self.ivar))**2:0, (sin(self.Omega*self.ivar))**2:0})).subs(data)
        #display(free_coeff)
        return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)

    @lru_cache
    def static_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re))**(1 / 2)

    @lru_cache
    def dynamic_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re))**(1 / 2)

    @lru_cache  
    def force_in_cable(self):
        data=self._given_data
        dyn_sys=self#.subs(data)
        dyn_sys_lin=dyn_sys.linearized()
        phi=dyn_sys_lin._fodes_system.steady_solution[1]
        
        force_in_cable = self.m_p*self.g*(1-S.One/2*phi**2) + self.m_p * self.l * phi.diff(self.ivar)**2
        force_subs=force_in_cable.subs(data)
        #force_subs=force_in_cable.subs({self.Omega:sqrt(self.g/self.l)}).subs(data)
        #display(force_subs.doit().expand())
        return force_subs.doit().expand()
    


# class DampedTrolleysWithSprings(ComposedSystem):
#     scheme_name = 'MDOF_Damped_Trolleys_With_Springs.PNG'
#     real_name = 'two_damped_trolleys.png'

#     def __init__(self,
#                  R=Symbol('R', positive=True),
#                  m=Symbol('m', positive=True),
#                  m1=Symbol('m_1', positive=True),
#                  m2=Symbol('m_2', positive=True),
#                  k_l=Symbol('k_l', positive=True),
#                  c_cl=Symbol('c_cl', positive=True),
#                  k_c=Symbol('k_c', positive=True),
#                  c_cc=Symbol('c_cc', positive=True),
#                  k_r=Symbol('k_r', positive=True),
#                  c_cr=Symbol('c_cr', positive=True),
#                  x_l=dynamicsymbols('x_l'),
#                  x_r=dynamicsymbols('x_r'),
#                  x=dynamicsymbols('x'),
#                  qs=dynamicsymbols('x_l x_r'),
#                  ivar=Symbol('t'),
#                  **kwargs):

#         self.m1 = m1
#         self.m2 = m2
#         self.m = m
#         self.k_l = k_l
#         self.c_cl = c_cl
#         self.k_c = k_c
#         self.c_cc = c_cc
#         self.k_r = k_r
#         self.c_cr = c_cr
#         self.x_l = x_l
#         self.x_r = x_r
#         self.x = x
#         self.R = R

#         self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
#                           + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
#                           + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
#                           + MaterialPoint(m/2, x_l/2, qs = [x_l]))

#         self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
#                           + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
#                           + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
#                           + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) 
#                           + MaterialPoint(m/2, x_r/2, qs = [x_r]))

#         system = self.Trolley_1 + self.Trolley_2
#         super().__init__(system(qs),**kwargs)

#     def get_default_data(self):

#         m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

#         default_data_dict = {
#             self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

#             self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
#             self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
#             self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

#             self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
#             self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
#             self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

#             self.x_l: [self.x, 0],
#             self.x_r: [self.x, 0],
#         }

#         return default_data_dict

#     def get_random_parameters(self):

#         default_data_dict = self.get_default_data()

#         parameters_dict = {
#             key: random.choice(items_list)
#             for key, items_list in default_data_dict.items()
#         }

#         if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

#             parameters_dict[self.x_l] = self.x


#         return parameters_dict

    
#Zrobione Amadi
class DampedTrolleyWithPendulum(TrolleyWithPendulum):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    l = Symbol('l', positive=True)
    m_t = Symbol('m_trolley', positive=True)
    m_p = Symbol('m_pendulum', positive=True)
    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')
    c = Symbol('c', positive=True)

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
                 c=None,
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
        if c is not None: self.c = c
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self._TrolleyWithPendulum = TrolleyWithPendulum(self.l, self.m_t, self.m_p, self.k, self.g, self.Omega, self.phi, self.x, self.F, self.ivar)(label='Trolley')
        self._trolley_damper=Damper(self.c, pos1=self.x, qs=[self.x, self.phi])(label='Trolley damper')
        self._pendulum_damper=Damper(self.c, pos1=self.phi, qs=[self.x, self.phi])(label='Pendulum damper')

        components['_TrolleyWithPendulum'] = self._TrolleyWithPendulum
        components['_trolley_damper'] = self._trolley_damper
        components['_pendulum_damper'] = self._pendulum_damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.k: r'Stiffness of a beam showed as a spring stiffness in trolley member',
            self.x: r'Kinematic lateral excitation',
            self.phi: r'Angle of a pendulum',
            self.m_p: r'Mass of pendulum',
            self.m_t: r'Mass of trolley',
            self.g: r'Gravity constant',
            self.F: r'Force',
            self.Omega: r'Excitation frequency',
            self.c: r'Damping coefficient'
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0, lam, lam0 = symbols('m_0 l_0 F_0 Omega_0 k_0 lambda lambda_0', positive=True)

        default_data_dict = {
            self.c: [S.One * self.k * lam],
            lam: [S.One * lam0 / 1000 * no for no in range(1, 10)]
        }
        return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.c: [no for no in range(5, 25)]
        }
        return default_data_dict
### Nieliniowe
class ForcedNonLinearTrolley(NonlinearComposedSystem):
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    d=Symbol('d', positive=True)
    l_0=Symbol('l_0', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F',positive=True)
    x=dynamicsymbols('x')
   
    def __init__(self,
                 m=None,
                 k=None,
                 d=None,
                 l_0=None,
                 Omega=None,
                 F=None,
                 x=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k= k
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if x is not None: self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(self.k, pos1 = sqrt(self.x**2 + self.d**2), pos2 = - self.l_0 , qs=[self.x])
        self._force = Force(self.F*sin(self.Omega*self.ivar), pos1 = self.x, qs=[self.x])


        components['trolley'] = self._trolley
        components['spring'] = self._spring
        components['force'] = self._force
        
        return components

    def get_default_data(self):

        m0, k0, l0= symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                0.5 * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
            self.d: [
                5 * l0, 2 * l0, 3 * S.Half * l0, 4 * l0, 6 * l0, 7 * l0,
                8 * l0, 9 * l0
            ],
            self.k: [
                1 * k0, 3 * k0, 2 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0,
                9 * k0
            ],
            self.l_0: [
                1 * l0, 3 * l0, 2 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0,
            ],
        }

        return default_data_dict

    def get_numerical_data(self): ### Brakowało numerical data. Przez to komenda get numerical parameters nie działa

        m0, k0, l0, c0= symbols('m_0 k_0 l_0 c_0', positive=True)

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
            self.l_0: [
                1 * l0, 3 * l0, 2 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0,
            ], ### Brakowało nawiasu
            self.c:  [
                1 * c0, 3 * c0, 2 * c0, 4 * c0, 5 * c0, 6 * c0, 7 * c0, 8 * c0,
                9 * c0
            ],
        }

        return default_data_dict
### NotODE
class TrolleysWithSprings(NonlinearComposedSystem):
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    k=Symbol('k', positive=True)
    #Omega=Symbol('Omega', positive=True)
    #F=Symbol('F', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 k=None,
                 #Omega=None,
                 #F=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if k is not None: self.k= k
        #if Omega is not None: self.Omega = Omega
        #if F is not None: self.F = F
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2

        self.qs = [self.x_1, self.x_2]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        #self._force = Force(self.F*sin(self.Omega*self.ivar), pos1 = self.x_1, qs=[self.x_1])

        self._spring12 = Spring(self.k, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(self.k, pos1 = self.x_2, qs=[self.x_2])


        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['spring_12'] = self._spring12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        #components['force'] = self._force
        
        return components



class ForcedTrolleysWithSprings(NonlinearComposedSystem): ### 3 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    k=Symbol('k', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    G=Symbol('G', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
    
    m0=Symbol('m_0', positive=True)
    k0=Symbol('k_0', positive=True)
    Omega0=Symbol('Omega_0', positive=True)
    F0=Symbol('F_0', positive=True)
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 k=None,
                 Omega=None,
                 F=None,
                 m0=None,
                 k0=None,
                 Omega0=None,
                 F0=None,
                 G=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if k is not None: self.k= k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2
            
        if m0 is not None: self.m0 = m0
        if k0 is not None: self.k0= k0
        if Omega0 is not None: self.Omega0 = Omega0
        if F0 is not None: self.F0 = F0

        self.qs = [self.x_1, self.x_2]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x_1, qs=[self.x_1])

        self._spring12 = Spring(self.k, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(self.k, pos1 = self.x_2, qs=[self.x_2])


        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['spring_12'] = self._spring12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        components['force'] = self._force
        
        return components

    def get_default_data(self):

        m0, k0, Omega0, F0 =  self.m0, self.k0, self.Omega0, self.F0
        
        default_data_dict = {
            self.m_1 : [100*m0],
            self.m_2 : [200*m0],
            self.k : [50*k0],
            self.Omega : [0.5 * 3.14*Omega0, 1 * 3.14*Omega0, 2 * 3.14*Omega0, 4 * 3.14*Omega0],
            self.F : [0.5 * 100*F0, 1 * 100*F0, 2 * 100*F0, 4 * 100*F0]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_1 : [100],
            self.m_2 : [200],
            self.k : [50],
            #self.c : [100],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_1: r'masa wózka nr 1',
            self.m_2: r'masa wózka nr 2',
            self.k: r'sztywność sprężyny',
            self.Omega: r'częstotliwość wymuszenia',
            self.F: r'wartość stałej siły wymuszającej',
            self.G: r'wartość zmiennej siły wymuszającej',
            self.x_1: 'przemieszczenie wózka nr 1',
            self.x_2: 'przemieszczenie wózka nr 2',
        }
        return self.sym_desc_dict
    
    def unit_dict(self):
        self.sym_desc_dict = {
            self.m_1: r'masa wózka nr 1',
            self.m_2: r'masa wózka nr 2',
            self.k: r'sztywność sprężyny',
            self.Omega: r'częstotliwość wymuszenia',
            self.F: r'wartość stałej siły wymuszającej',
            self.G: r'wartość zmiennej siły wymuszającej',
            self.x_1: 'przemieszczenie wózka nr 1',
            self.x_2: 'przemieszczenie wózka nr 2',
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        #mech_comp.SteadySolutionComponent,
        ]
        return comp_list
    
    
class ForcedDampedTrolleysWithSprings(ComposedSystem):
    
    scheme_name = 'MDOF_Damped_Trolleys_With_Springs.png'
    real_name = 'two_trolleys_damped_with_springs_real.png'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    m=Symbol('m', positive=True)
    R=Symbol('R', positive=True)
    k_l=Symbol('k_l', positive=True)
    k_r=Symbol('k_r', positive=True)
    k_c=Symbol('k_c', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_r=Symbol('c_r', positive=True)
    c_c=Symbol('c_c', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    x_1=dynamicsymbols('x_l')
    x_2=dynamicsymbols('x_r')
   
    def __init__(self,
                 R=None,
                 m_1=None,
                 m_2=None,
                 m=None,
                 k_l=None,
                 k_r=None,
                 k_c=None,
                 c_l=None,
                 c_r=None,
                 c_c=None,
                 Omega=None,
                 F=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if m is not None: self.m = m
        if k_l is not None: self.k_l = k_l
        if k_r is not None: self.k_r = k_r
        if k_c is not None: self.k_c = k_c
        if c_l is not None: self.c_l = c_l
        if c_r is not None: self.c_r = c_r
        if c_c is not None: self.c_c = c_c
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2

        self.qs = [self.x_1, self.x_2]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._wheel11 = MaterialPoint(self.m, self.x_1, qs = [self.x_1])
        self._wheel11_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_1/self.R, qs=[self.x_1])
        self._wheel12= MaterialPoint(self.m, self.x_1, qs = [self.x_1])
        self._wheel12_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_1/self.R, qs=[self.x_1])
        
        self._spring1_top = Spring(self.k_l, pos1 = self.x_1, qs=[self.x_1])
        self._spring1_bottom = Spring(self.k_l, pos1 = self.x_1, qs=[self.x_1])
        
        self._damper1 = Damper(self.c_l, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.F, pos1 = self.x_1, qs=[self.x_1])

        self._spring12 = Spring(2*self.k_c, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._damper12 = Damper(self.c_c, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._wheel21 = MaterialPoint(self.m, self.x_2, qs = [self.x_2])
        self._wheel21_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_2/self.R, qs=[self.x_2])
        self._wheel22 = MaterialPoint(self.m, self.x_2, qs = [self.x_2])
        self._wheel22_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_2/self.R, qs=[self.x_2])
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(2*self.k_r, pos1 = self.x_2, qs=[self.x_2])
        self._damper2 = Damper(self.c_r, pos1 = self.x_2, qs=[self.x_2])


        components['trolley_1'] = self._trolley1
        components['spring_1_top'] = self._spring1_top
        components['spring_1_bottom'] = self._spring1_bottom
        components['damper_1'] = self._damper1
        components['spring_12'] = self._spring12
        components['wheel12'] = self._wheel12
        components['wheel11'] = self._wheel11
        components['wheel11_disk'] = self._wheel11_disk
        components['wheel12_disk'] = self._wheel12_disk
        components['damper_12'] = self._damper12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        components['damper_2'] = self._damper2
        components['wheel22'] = self._wheel22
        components['wheel21'] = self._wheel21
        components['wheel22_disk'] = self._wheel22_disk
        components['wheel21_disk'] = self._wheel21_disk
        components['force'] = self._force
        
        return components

    def get_default_data(self):

        m0, k0, lam, Omega0, F0 = symbols('m_0 k_0 lambda Omega_0 F_0', positive=True)

        default_data_dict = {
            self.k_l: [S.One * self.k_c],
            self.k_r: [S.One * self.k_c],
            self.c_r: [S.One * lam * self.k_r],
            self.c_l: [S.One * lam * self.k_l],
            self.c_c: [S.One * lam * self.k_c],

            self.m_2: [S.One * self.m],
            self.m_1: [S.One* self.m],
            
#             self.Omega: [S.One * Omega0],
            self.F: [S.One * F0 * no for no in range(5,25)],
            self.m: [S.One * m0 * no for no in range(1,5)],
            self.m_2: [S.One * m0 * no for no in range(1,5)],
            self.m_1: [S.One * m0 * no for no in range(1,5)],
            
            

            self.k_c: [S.One * k0 * no for no in range(60, 90)],
            self.k_l: [S.One * k0 * no for no in range(60, 90)],
            self.k_r: [S.One * k0 * no for no in range(60, 90)],            
#             self.c_r: [S.One * lam * self.k_c],
#             self.c_l: [S.One * lam * self.k_c],
#             self.c_c: [S.One * lam * self.k_c],
            
        }

        return default_data_dict

    
#     def get_default_data(self):

#         m0, k0, lam, Omega0, F0 = symbols('m_0 k_0 lambda Omega_0 F_0', positive=True)

#         default_data_dict = {
#             self.k_l: [S.One * self.k_c],
#             self.k_r: [S.One * self.k_c],
#             self.c_r: [S.One * lam * self.k_r],
#             self.c_l: [S.One * lam * self.k_l],
#             self.c_c: [S.One * lam * self.k_c],

#             self.m_2: [S.One * self.m],
#             self.m_1: [S.One* self.m],
            
# #             self.Omega: [S.One * Omega0],
#             self.F: [S.One * F0 * 20],
#             self.m: [S.One * m0 * 2],
#             self.m_2: [S.One * m0 * no for no in range(1,5)],
#             self.m_1: [S.One * m0 * no for no in range(1,5)],
            
            

#             self.k_c: [S.One * k0 * no for no in  [74,74.01,74.02]],
#             self.k_l: [S.One * k0 * 65],
#             self.k_r: [S.One * k0 * 65],            
# #             self.c_r: [S.One * lam * self.k_c],
# #             self.c_l: [S.One * lam * self.k_c],
# #             self.c_c: [S.One * lam * self.k_c],
            
#         }

#         return default_data_dict    
    
    
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_1 : [100],
            self.m_2 : [200],
            self.k_l : [50],
            self.k_r : [50],
            self.k_c : [50],
            self.c_l : [100],
            self.c_r : [100],
            self.c_c : [100],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
#         default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        parameters_dict[self.k_r] = parameters_dict[self.k_l]
        parameters_dict[self.m_1] = parameters_dict[self.m]
        parameters_dict[self.m_2] = parameters_dict[self.m]

        return parameters_dict
#     def subs(self, *args, **kwargs):
        
#         return super().subs(*args, **kwargs).subs(*args, **kwargs)



    def max_static_force_pin(self):

        ans=self.static_force()
        
        return abs(ans)


    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def static_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_static_force_pin())/(pi*kt*Re))**(1/2)


    
    def dynamic_force(self):
        
        amps = self._fodes_system.steady_solution.as_dict()
        #force = self.components['_spring_l'].force().doit().expand()#.doit()
        data=self._given_data

        
        return (self.components['spring_1_top'].force().subs(amps)).subs(data).expand().doit()

    def max_dynamic_force_pin(self):
        
        ans = self.dynamic_force()
        
        #display(abs(self.components['_spring_l'].force()))
        data=self._given_data
        sin_coeff=ans.coeff(sin(self.Omega*self.ivar)).subs(data)
        #display(sin_coeff)
        cos_coeff=ans.coeff(cos(self.Omega*self.ivar)).subs(data)
        #display(cos_coeff)
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        #display(free_coeff)

        return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)

    
    def dynamic_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        
        
        load = abs(self.max_dynamic_force_pin())
        
        return ((4*load)/(pi*kt*Re))**(1/2)
    
     
    def right_spring_force(self):
        k_r=self.k_r
        x_r=self.x_2
        sol_dict=self._fodes_system.steady_solution.as_dict()
        Fr=(-k_r*x_r).subs(sol_dict).subs(self._given_data)
        
        return Fr
    
    def centre_spring_force(self):
        k_c=self.k_c
        x_r=self.x_2
        x_l=self.x_1
        sol_dict=self._fodes_system.steady_solution.as_dict()
        Fc=(-k_c*(x_r-x_l)).subs(sol_dict).subs(self._given_data)
        
        return Fc
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.RightSpringForceComponent,
        mech_comp.CentralSpringForceComponent

        ]

        return comp_list



### Nieliniowe
class ForcedTrolleysWithNonLinearSprings(NonlinearComposedSystem):
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    k=Symbol('k', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    G=Symbol('G', positive=True)
    mu = Symbol('\\mu', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 k=None,
                 Omega=None,
                 F=None,
                 G=None,
                 mu=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if k is not None: self.k= k
        if mu is not None: self.mu= mu
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2

        self.qs = [self.x_1, self.x_2]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k*(self.mu*self.x_1**2), pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar), pos1 = self.x_1, qs=[self.x_1])

        self._spring12 = Spring(self.k*(self.mu*(self.x_1-self.x_2)**2), pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(self.k*(self.mu*self.x_2**2), pos1 = self.x_2, qs=[self.x_2])


        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['spring_12'] = self._spring12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        components['force'] = self._force
        
        return components
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def dynamic_force(self):
        amps = self._fodes_system.steady_solution.as_dict()
        data=self._given_data

        return (self.components['spring_1'].force().subs(amps)).subs(data).expand().doit()
    


class ForcedTrolleyWithSpring(ComposedSystem): ### 1 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F',positive=True)
    G=Symbol('G',positive=True)
    x=dynamicsymbols('x')
    
    m0=Symbol('m_0', positive=True)
    k0=Symbol('k_0', positive=True)
    Omega0=Symbol('Omega_0', positive=True)
    F0=Symbol('F_0', positive=True)
   
    def __init__(self,
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
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k= k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x is not None: self.x = x
            
        if m0 is not None: self.m0 = m0
        if k0 is not None: self.k0= k0
        if Omega0 is not None: self.Omega0 = Omega0
        if F0 is not None: self.F0 = F0

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(self.k, pos1 = self.x , qs=[self.x])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x, qs=[self.x])



        components['trolley'] = self._trolley
        components['spring'] = self._spring
        components['force'] = self._force
        
        return components

    def get_default_data(self):

        m0, k0, Omega0, F0 =  self.m0, self.k0, self.Omega0, self.F0
        
        default_data_dict = {
            self.m : [100*m0],
            self.k : [50*k0],
            self.Omega : [0.5 * 3.14*Omega0, 1 * 3.14*Omega0, 2 * 3.14*Omega0, 4 * 3.14*Omega0],
            self.F : [0.5 * 100*F0, 1 * 100*F0, 2 * 100*F0, 4 * 100*F0]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m : [100],
            self.k : [50],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict


    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'masa wózka',
            self.k: r'sztywność sprężyny',
            self.Omega: r'częstotliwość wymuszenia',
            self.F: r'wartość stałej siły wymuszającej',
            self.G: r'wartość zmiennej siły wymuszającej',
            self.x: 'przemieszczenie wózka nr 1',
        }
        return self.sym_desc_dict
    
    def unit_dict(self):
        unit_dict = {
            self.m: ureg.kilogram,
            self.k: ureg.newton/ureg.meter,
            self.Omega: ureg.hertz,
            self.F: ureg.newton,
            self.G: ureg.newton,
            self.x: ureg.meter,
        }
        return unit_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        #mech_comp.FundamentalMatrixComponent,
        #mech_comp.GeneralSolutionComponent,
        #mech_comp.SteadySolutionComponent,
            
            
        ]
        
        return comp_list
    


class ForcedDampedTrolleyWithSpring(ComposedSystem): ### 2 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    c=Symbol('c', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F',positive=True)
    G=Symbol('G',positive=True)
    x=dynamicsymbols('x')
   
    def __init__(self,
                 m=None,
                 k=None,
                 c=None,
                 Omega=None,
                 F=None,
                 G=None,
                 x=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k= k
        if c is not None: self.c= c
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x is not None: self.x = x

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley = MaterialPoint(self.m, self.x, qs=[self.x])
        self._spring = Spring(self.k, pos1 = self.x , qs=[self.x])
        self._damper = Damper(self.c, pos1 = self.x , qs=[self.x])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x, qs=[self.x])


        components['trolley'] = self._trolley
        components['spring'] = self._spring
        components['damper'] = self._damper
        components['force'] = self._force
        
        return components

    def get_numerical_data(self):

        default_data_dict = {
            self.m : [100],
            self.k : [50],
            self.c : [100],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
class ForcedThreeTrolleysWithSprings(ComposedSystem): ### 7 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    m_3=Symbol('m_3', positive=True)
    k=Symbol('k', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    G=Symbol('G', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
    x_3=dynamicsymbols('x_3')
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 m_3=None,
                 k=None,
                 Omega=None,
                 F=None,
                 G=None,
                 x_1=None,
                 x_2=None,
                 x_3=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if m_3 is not None: self.m_3 = m_3
        if k is not None: self.k= k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2
        if x_3 is not None: self.x_3 = x_3

        self.qs = [self.x_1, self.x_2, self.x_3]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x_1, qs=[self.x_1])
        
        self._spring12 = Spring(self.k, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(self.k, pos1 = self.x_2, qs=[self.x_2])
        
        self._spring23 = Spring(self.k, pos1 = self.x_2, pos2 = self.x_3, qs=[self.x_2, self.x_3])
        self._trolley3 = MaterialPoint(self.m_3, self.x_3, qs=[self.x_3])
        self._spring3 = Spring(self.k, pos1 = self.x_3, qs=[self.x_3])

        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['spring_12'] = self._spring12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        components['spring_23'] = self._spring23
        components['trolley_3'] = self._trolley3
        components['spring_3'] = self._spring3
        components['force'] = self._force
        
        return components

    def get_numerical_data(self):

        default_data_dict = {
            self.m_1 : [100],
            self.m_2 : [200],
            self.m_3 : [300],
            self.k : [50],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    


class ForcedDisconnectedTrolleysWithSprings(ForcedThreeTrolleysWithSprings): ### 5 ODE

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar), pos1 = self.x_1, qs=[self.x_1])
        
        self._spring12 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        #self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring23 = Spring(self.k, pos1= -self.x_3, qs=[self.x_3])
        
        self._trolley3 = MaterialPoint(self.m_3, self.x_3, qs=[self.x_3])
        self._spring3 = Spring(self.k, pos1 = self.x_3, qs=[self.x_3])

        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['spring_12'] = self._spring12
        #components['trolley_2'] = self._trolley2
        components['spring_23'] = self._spring23
        components['trolley_3'] = self._trolley3
        components['spring_3'] = self._spring3
        components['force'] = self._force
        
        return components    
    


class ForcedDampedThreeTrolleysWithSprings(ComposedSystem): ### 8 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    m_3=Symbol('m_3', positive=True)
    k=Symbol('k', positive=True)
    c=Symbol('c', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    G=Symbol('G', positive=True) ## Siła
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
    x_3=dynamicsymbols('x_3')
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 m_3=None,
                 k=None,
                 c=None,
                 Omega=None,
                 F=None,
                 G=None,
                 x_1=None,
                 x_2=None,
                 x_3=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if m_3 is not None: self.m_3 = m_3
        if k is not None: self.k= k
        if c is not None: self.c= c
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if G is not None: self.G = G
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None: self.x_2 = x_2
        if x_3 is not None: self.x_3 = x_3

        self.qs = [self.x_1, self.x_2, self.x_3]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._damper1 = Damper(self.c, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x_1, qs=[self.x_1])
        
        self._spring12 = Spring(self.k, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._damper12 = Damper(self.c, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        
        self._spring23 = Spring(self.k, pos1 = self.x_2, pos2 = self.x_3, qs=[self.x_2, self.x_3])
        self._damper23 = Damper(self.c, pos1 = self.x_2, pos2 = self.x_3, qs=[self.x_2, self.x_3])
        
        self._trolley3 = MaterialPoint(self.m_3, self.x_3, qs=[self.x_3])
        self._spring3 = Spring(self.k, pos1 = self.x_3, qs=[self.x_3])
        self._damper3 = Damper(self.c, pos1 = self.x_3, qs=[self.x_3])

        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['damper_1'] = self._damper1
        components['spring_12'] = self._spring12
        components['damper_12'] = self._damper12
        components['trolley_2'] = self._trolley2
        components['spring_23'] = self._spring23
        components['damper_23'] = self._damper23
        components['trolley_3'] = self._trolley3
        components['spring_3'] = self._spring3
        components['damper_3'] = self._damper3
        components['force'] = self._force
        
        return components
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_1 : [100],
            self.m_2 : [200],
            self.m_3 : [300],
            self.k : [50],
            self.c : [100],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict



class ForcedDampedDisconnectedTrolleysWithSprings(ForcedDampedThreeTrolleysWithSprings): ### 6 ODE

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_1, qs=[self.x_1])
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._damper1 = Damper(self.c, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x_1, qs=[self.x_1])
        
        self._spring12 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._damper12 = Damper(self.c, pos1 = self.x_1, qs=[self.x_1])
        
        #self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        
        self._spring23 = Spring(self.k, pos1 = -self.x_3, qs=[self.x_3])
        self._damper23 = Damper(self.c, pos1 = -self.x_3, qs=[self.x_3])
        
        self._trolley3 = MaterialPoint(self.m_3, self.x_3, qs=[self.x_3])
        self._spring3 = Spring(self.k, pos1 = self.x_3, qs=[self.x_3])
        self._damper3 = Damper(self.c, pos1 = self.x_3, qs=[self.x_3])

        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['damper_1'] = self._damper1
        components['spring_12'] = self._spring12
        components['damper_12'] = self._damper12
        #components['trolley_2'] = self._trolley2
        components['spring_23'] = self._spring23
        components['damper_23'] = self._damper23
        components['trolley_3'] = self._trolley3
        components['spring_3'] = self._spring3
        components['damper_3'] = self._damper3
        components['force'] = self._force
        
        return components



class DDoFTwoNonLinearTrolleys(NonlinearComposedSystem):
    
    scheme_name = 'ForcedTwoNonLinearTrolleys.png'
    real_name = 'tension_leg_platform.png'

    g=Symbol('g', positive=True)
    m1=Symbol('m_1', positive=True)
    m2=Symbol('m_2', positive=True)
    k1=Symbol('k_1', positive=True)
    k2=Symbol('k_2', positive=True)
    k3=Symbol('k_3', positive=True)
    d=Symbol('d', positive=True)
    l_0=Symbol('l_0', positive=True)
    ivar=Symbol('t')
    x1=dynamicsymbols('x1')
    x2=dynamicsymbols('x2')
    x=dynamicsymbols('x')
    qs=dynamicsymbols('x1, x2')
    F1=Symbol('F_1', positive=True)
    F2=Symbol('F_2', positive=True)
    Omega=Symbol('Omega', positive=True)
   
    def __init__(self,
                 g=None,
                 m1=None,
                 m2=None,
                 k1=None,
                 k2=None,
                 k3=None,
                 d=None,
                 l_0=None,
                 ivar=Symbol('t'),
                 x1=None,
                 x2=None,
                 x=None,
                 qs=None,
                 F1=None,
                 F2=None,
                 Omega=None,
                 **kwargs):

        if m1 is not None: self.m1 = m1
        if m2 is not None: self.m2 = m2
        if k1 is not None: self.k1= k1
        if k2 is not None: self.k2= k2
        if k3 is not None: self.k3= k3
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if Omega is not None: self.Omega = Omega
        if F1 is not None: self.F1 = F1
        if F2 is not None: self.F2 = F2
        if x is not None: self.x = x
        if x1 is not None: self.x1 = x1
        if x2 is not None: self.x2 = x2
        if qs is not None: self.qs = qs
        if g is not None: self.g = g
        self.ivar=ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m1, self.x1, qs=self.qs)
        self._spring1 = Spring(self.k1, pos1=(sqrt(self.x1**2 + self.d**2) - self.l_0), qs=self.qs)
        self._trolley2 = MaterialPoint(self.m2, self.x2, qs=self.qs)
        self._spring2 = Spring(self.k2, pos1=(sqrt(self.x2**2 + self.d**2) - self.l_0), qs=self.qs)
        self._spring12 = Spring(self.k3, self.x1, self.x2,qs=self.qs)
        self._force1 = Force(self.F1*sin(self.Omega*self.ivar), self.x1, qs=self.qs)
        self._force2 = Force(self.F2*sin(self.Omega*self.ivar), self.x2, qs=self.qs)

        components['_trolley1'] = self._trolley1
        components['_trolley2'] = self._trolley2
        components['_spring1'] = self._spring1
        components['_spring2'] = self._spring2
        components['_spring12'] = self._spring12
        components['_force1'] = self._force1
        components['_force2'] = self._force2
        
        return components

    def get_default_data(self):

        m0, k0, l0, F0 = symbols('m_0 k_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.m1: [S.One * m0 * no for no in range(20,50)],
            self.m2: [S.One * m0 * no for no in range(20,50)],
            self.d: [S.One * l0 * no for no in range(2,5)],
            self.k1: [S.One * k0 * no for no in range(10,20)],
            self.k2: [S.One * k0 * no for no in range(10,20)],
            self.k3: [S.One * k0 * no for no in range(10,20)],
            self.F1: [S.One * F0 * no for no in range(5,15)],
            self.F2: [S.One * F0 * no for no in range(5,15)],
            self.x1: [self.x, 0],
            self.x2: [self.x, S.Zero],

        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x1] == S.Zero:
            parameters_dict[self.x2] = self.x

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley Mass',
            self.m2: r'Trolley Mass',
            self.k1: 'Spring Stiffness',
            self.k2: 'Spring Stiffness',
            self.k3: 'Spring Stiffness',
            self.d: r'length',
            self.l_0: r'length',
        }
        return self.sym_desc_dict

    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)

#     def dynamic_force(self):

#         data=self._given_data
#         amps = self._fodes_system.steady_solution.as_dict()

#         return (self.components['_spring12'].force().subs(amps)).subs(data).expand().doit()

    def dynamic_force(self):

        data=self._given_data
        amps = self.linearized().subs(data)._fodes_system.steady_solution.as_dict()


        return (self.components['_spring12'].force().subs(data).subs(amps)).expand().doit()



class VariableMassTrolleyWithPendulum(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    l = Symbol('l', positive=True)
    
    m_f = Symbol('m_f', positive=True)

    m_t = Symbol('m_t', positive=True)
    m_t0 = Symbol('m_t0', positive=True)
    m_tf = Symbol('m_tf', positive=True)

    m_p = Symbol('m_p', positive=True)
    m_p0 = Symbol('m_p0', positive=True)
    m_pf = Symbol('m_pf', positive=True)
    
    flow_coeff = Symbol('\\lambda', positive=True)
    
    t0 = Symbol('t_0', positive=True)
    
    m_0 = Symbol('m_0', positive=True)

    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')

    def __init__(self,
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
                 t0=None,
                 k=None,
                 g=None,
                 Omega=None,
                 phi=None,
                 x=None,
                 F=None,
                 ivar=Symbol('t'),
                 **kwargs):
        
        if l is not None: self.l = l
        if m_f is not None: self.m_f = m_f
        if m_t is not None: self.m_t = m_t
        if m_t0 is not None: self.m_t0 = m_t0
        if m_tf is not None: self.m_tf = m_tf
        if m_p is not None: self.m_p = m_p
        if m_p0 is not None: self.m_p0 = m_p0
        if m_pf is not None: self.m_pf = m_pf
        if m_0 is not None: self.m_0 = m_0
        if flow_coeff is not None: self.flow_coeff = flow_coeff
        if t0 is not None: self.t0 = t0
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        self.ivar = ivar
        self._init_from_components(**kwargs)

        self.trans_expr = ((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        
        #self.m_tf = self.m_f#((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        #self.m_pf = self.m_f#((S.One/2+atan(self.flow_coeff*(self.ivar-self.t0))/pi))

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t + self.m_tf, self.k, self.x, self.ivar)(label='Trolley')
        self._pendulum = PendulumKinematicExct(self.l, self.m_p + self.m_pf, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._force=Force(self.F*sin(self.Omega*self.ivar), pos1=self.x, qs=[self.x, self.phi])(label='Force')

        components['_trolley'] = self._trolley
        components['_pendulum'] = self._pendulum
        components['_force'] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.k: r'Stiffness of a beam showed as a spring stiffness in trolley member',
            self.x: r'Kinematic lateral excitation',
            self.phi: r'Angle of a pendulum',
            self.m_p: r'Mass of pendulum',
            self.m_t: r'Mass of trolley',
            self.g: r'Gravity constant',
            self.F: r'Force',
            self.Omega: r'Excitation frequency',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols('m_0 l_0 F_0 Omega_0 k_0', positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.g: [S.One * self.g],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x]
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
            self.t0: [10],
            self.Omega: [3.14 * 0.5],
            self.g: [9.81],
            #self.k: [2]
        }
        default_data_dict.update({self.m_t: [10.0*default_data_dict[self.m_p][0]],
                                  #self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
                                  #self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)]
                                 })
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_t: r'masa początkowa wózka',
            self.m_p: r'masa początkowa wahadła',
            self.m_tf: r'masa transportowanej cieczy wypływającej z wózka',
            self.m_pf: r'masa transportowanej cieczy wypływającej do wahadła',
            self.m_f: r'masa transferowanej cieczy',
            self.k: r'sztywność sprężyny mocującej',
            self.l: r'długość wahadła',
            self.Omega: r'częstotliwość wymuszenia',
            self.F: r'wartość siły wymuszającej',
            self.g: r'przyspieszenie ziemskie',
            self.flow_coeff: r'współczynnik przepływu cieczy',
            self.t0: r'czas aktywacji tłumienia',
        }
        return self.sym_desc_dict
    


class VariableMassTrolleyWithPendulumRayleighDamping(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    l = Symbol('l', positive=True)
    
    m_f = Symbol('m_f', positive=True)

    m_t = Symbol('m_t', positive=True)
    m_t0 = Symbol('m_t0', positive=True)
    m_tf = Symbol('m_tf', positive=True)
    m_tf_eq = Symbol('m_tf_eq', positive=True)

    m_p = Symbol('m_p', positive=True)
    m_p0 = Symbol('m_p0', positive=True)
    m_pf = Symbol('m_pf', positive=True)
    m_pf_eq = Symbol('m_pf_eq', positive=True)
    
    flow_coeff = Symbol('\\lambda', positive=True)
    
    t0 = Symbol('t_0', positive=True)
    
    m_0 = Symbol('m_0', positive=True)

    k = Symbol('k', positive=True)
    g = Symbol('g', positive=True)
    
    alpha = Symbol('\\alpha')
    beta = Symbol('\\beta')
    rayleigh_damping_matrix = Symbol('rayleigh_damping_matrix', positive=True)
    
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('\\varphi')
    x = dynamicsymbols('x')
    
    b = Symbol('b', positive=True)
    
    omega = Symbol('omega', positive=True)
    Ek = Symbol('E_k', positive=True)
    Ep = Symbol('E_p', positive=True)
    
    f = Symbol('f', positive=True)

    def __init__(self,
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
                 t0=None,
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
                 ivar=Symbol('t'),
                 **kwargs):
        
        if l is not None: self.l = l
        if m_f is not None: self.m_f = m_f
        if m_t is not None: self.m_t = m_t
        if m_t0 is not None: self.m_t0 = m_t0
        if m_tf is not None: self.m_tf = m_tf
        if m_pf_eq is not None: self.m_tf_eq = m_tf_eq
        if m_p is not None: self.m_p = m_p
        if m_p0 is not None: self.m_p0 = m_p0
        if m_pf is not None: self.m_pf = m_pf
        if m_pf_eq is not None: self.m_pf_eq = m_pf_eq
        if m_0 is not None: self.m_0 = m_0
        if flow_coeff is not None: self.flow_coeff = flow_coeff
        if rayleigh_damping_matrix is not None: self.rayleigh_damping_matrix = rayleigh_damping_matrix
        if t0 is not None: self.t0 = t0
        if g is not None: self.g = g
        if b is not None: self.b = b
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if omega is not None: self.omega = omega
        if Ek is not None: self.Ek = Ek
        if Ep is not None: self.Ep = Ep
        if F is not None: self.F = F
        if f is not None: self.f = f
        self.ivar = ivar
        self.trans_expr = ((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        self.alpha = self.b
        self.beta = self.b/2
        #self.rayleigh_damping_matrix = self.alpha*VariableMassTrolleyWithPendulum().inertia_matrix() + self.beta*VariableMassTrolleyWithPendulum().stiffness_matrix()
        self._init_from_components(**kwargs)



    @cached_property
    def components(self):
        components = {}
        
        self.m_tf_eq = self.m_f*((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        self.m_pf_eq = self.m_f*((S.One/2+atan(self.flow_coeff*(self.ivar-self.t0))/pi))

        self._trolley = VariableMassTrolleyWithPendulum(m_tf = self.m_tf, m_pf = self.m_pf)(label='Trolley')
        
        #damping_matrix = self.alpha*self._trolley.inertia_matrix() + self.beta*self._trolley.stiffness_matrix()
        
#         self._trolley_damping = Damper((damping_matrix[0]+damping_matrix[2]), pos1 = self.x , qs = [self.x,self.phi])(label='Trolley dapming')
#         self._pendulum_damping = Damper((damping_matrix[1]+damping_matrix[3]), pos1 = self.phi , qs = [self.x,self.phi])(label='Pendulum dapming')

        self._trolley_damping = Damper(self.b, pos1 = self.x , qs = [self.x,self.phi])(label='Trolley dapming')
        self._pendulum_damping = Damper(self.b/2, pos1 = self.phi , qs = [self.x,self.phi])(label='Pendulum dapming')
        
        components['_trolley'] = self._trolley
        components['_trolley_damping'] = self._trolley_damping
        components['_pendulum_damping'] = self._pendulum_damping

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
            #self.alpha: [0.8],
            #self.beta: [2.5],
            self.b: [2.5],
            self.F: [250],
            self.flow_coeff: [10],
            self.t0: [30],
            self.Omega: [3.14 * 0.65],
            self.g: [9.81],
            #self.k: [2]
        }
        default_data_dict.update({self.m_t: [10.0*default_data_dict[self.m_p][0]],
                                  self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
                                  self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)],
                                  self.m_tf: [default_data_dict[self.m_f][0]*((S.One/2-atan(default_data_dict[self.flow_coeff][0]*(self.ivar-default_data_dict[self.t0][0]))/pi))],
                                  self.m_pf: [default_data_dict[self.m_f][0]*((S.One/2+atan(default_data_dict[self.flow_coeff][0]*(self.ivar-default_data_dict[self.t0][0]))/pi))],
                                          })
        
        default_data_dict.update({self.m_t+self.m_tf: [default_data_dict[self.m_t][0]+default_data_dict[self.m_tf][0]],
                                  self.m_p+self.m_pf: [default_data_dict[self.m_p][0]+default_data_dict[self.m_pf][0]],
                                          })
        
        return default_data_dict
    

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_t: r'masa początkowa wózka',
            self.m_p: r'masa początkowa wahadła',
            self.m_tf: r'masa transportowanej cieczy wypływającej z wózka',
            self.m_pf: r'masa transportowanej cieczy wypływającej do wahadła',
            self.m_f: r'masa transferowanej cieczy',
            self.k: r'sztywność sprężyny mocującej',
            self.l: r'długość wahadła',
            self.Omega: r'częstotliwość wymuszenia',
            self.F: r'wartość siły wymuszającej',
            self.g: r'przyspieszenie ziemskie',
            self.flow_coeff: r'współczynnik przepływu cieczy',
            self.t0: r'czas aktywacji tłumienia',
            self.b: r'współczynnik tłumienia',
            self.alpha: r'współczynnik tłumienia Rayleigha przy macierzy bezwładności',
            self.beta: r'współczynnik tłumienia Rayleigha przy macierzy sztywności',
            self.x: r'przemieszczenie wózka',
            self.phi: r'kąt wychylenia wahadła',
            self.omega: r'częstość własna układu',
            self.Ek: r'energia kinetyczna układu',
            self.Ep: r'energia potencjalna układu',
            self.f: r'częstotliwość',
            self.ivar: r'czas',
            diff(self.x,self.ivar): 'prędkość wózka',
            diff(self.x,self.ivar,self.ivar): r'przyspieszenie wózka',
            diff(self.phi,self.ivar): r'prędkość kątowa wahadła',
            diff(self.phi,self.ivar,self.ivar): r'przyspieszenie kątowe wahadła',
        }
        return self.sym_desc_dict


#Amadi
class TrolleyWithTMD(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    
    m = Symbol('m_trolley', positive=True)
    m_TMD = Symbol('m_TMD', positive=True)
    k = Symbol('k', positive=True)
    k_TMD = Symbol('k_TMD', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    x_e = dynamicsymbols('x_e')
    x = dynamicsymbols('x')

    def __init__(self,
                 m=None,
                 m_TMD=None,
                 k=None,
                 k_TMD=None,
                 Omega=None,
                 x_e=None,
                 x=None,
                 F=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if m_TMD is not None: self.m_TMD = m_TMD
        if m is not None: self.m = m
        if x_e is not None: self.x_e = x_e
        if x is not None: self.x = x
        if k is not None: self.k = k
        if k_TMD is not None: self.k_TMD = k_TMD
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        self.ivar = ivar
        self._init_from_components(**kwargs)
        

    @property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m, self.k, self.x, self.ivar)(label='Trolley')
        self._TMD = TunedMassDamperRelativeMotion(self.m_TMD, self.k_TMD, self.x_e, self.x, self.ivar)(label='Tuned Mass Damper')
        self._force=Force(self.F*sin(self.Omega*self.ivar), pos1=self.x, qs=[self.x, self.x_e])(label='Force')

        components['_trolley'] = self._trolley
        components['_TMD'] = self._TMD
        components['_force'] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r'Stiffness of a beam showed as a spring stiffness in trolley member',
            self.x: r'Kinematic lateral excitation',
            self.x_e: r'Angle of a pendulum',
            self.m: r'Mass of trolley',
            self.m_TMD: r'Mass of TMD',
            self.F: r'Force',
            self.Omega: r'Excitation frequency',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, F0, Omega0, k0 = symbols('m_0 F_0 Omega_0 k_0', positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(20, 30)],
            self.m_TMD: [S.One * m0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x]
        }
        return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m: [no for no in range(20, 30)],
            self.m_TMD: [no for no in range(1, 10)],
            self.F: [no for no in range(50, 100)],
            self.Omega: [3.14 * no for no in range(1,6)],
            self.k: [no for no in range(50, 100)],
        }
        return default_data_dict



class VariableMassTrolleyWithPendulumFunction(ComposedSystem):

    scheme_name = 'trolley_pendulum_tmd.png'
    real_name = 'taipei101.png'

    ivar=Symbol('t')
    l = Symbol('l', positive=True)
    rho = Symbol('rho', positive=True)
    m_f = Symbol('epsilon_f', positive=True)
    m_t = Function('m_t')(ivar)
    m_t0 = Symbol('m_t0', positive=True)
    m_p = Function('m_p')(ivar)
    m_p0 = Symbol('m_p0', positive=True)
    lam = Symbol('lambda', positive=True)
    t0 = Symbol('tau_0', real=True)
    k = Symbol('k', positive=True)
    omega0=Symbol('omega_0', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    phi = dynamicsymbols('Phi')
    x = dynamicsymbols('X')
    eps = Symbol('epsilon', positive=True)
    m_s = Symbol('m_s', positive=True)
    chi = Symbol('chi', positive=True)
    psi = Symbol('psi', positive=True)
    delta = Symbol('delta', positive=True)
    Omega0=Symbol('Omega_0', positive=True)
    tau=Symbol('tau', positive=True)

    def __init__(self,
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
                 **kwargs):
        
        if l is not None: self.l = l
        if rho is not None: self.rho = rho
        if m_f is not None: self.m_f = m_f
        if m_t is not None: self.m_t = m_t
        if m_t0 is not None: self.m_t0 = m_t0
        if m_p is not None: self.m_p = m_p
        if m_p0 is not None: self.m_p0 = m_p0
        if lam is not None: self.lam = lam
        if t0 is not None: self.t0 = t0
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if omega0 is not None: self.omega0 = omega0
        if F is not None: self.F = F
        if eps is not None: self.eps = eps
        if m_s is not None: self.m_s = m_s
        if chi is not None: self.chi = chi
        if psi is not None: self.psi = psi
        if delta is not None: self.delta = delta
        if Omega0 is not None: self.Omega0 = Omega0
        if tau is not None: self.tau = tau
        self.ivar = ivar
        self._init_from_components(**kwargs)

        #self.trans_expr = ((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        #self.m_tf = self.m_f#((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        #self.m_pf = self.m_f#((S.One/2+atan(self.flow_coeff*(self.ivar-self.t0))/pi))

    @cached_property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(label='Trolley')
        self._trolley_damper = Damper(c=self.k*self.lam, pos1=self.x, pos2=0, qs=[self.x, self.phi])(label='Trolley damper')
        self._pendulum = PendulumKinematicExct(self.l, self.m_p, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._pendulum_damper = Damper(c=self.m_p*self.g*self.l*self.lam, pos1=self.phi, pos2=0, qs=[self.x, self.phi])(label='Pendulum damper')
        self._force=Force(self.F*sin(self.Omega*self.ivar), pos1=self.x, qs=[self.x, self.phi])(label='Force')

        components['_trolley'] = self._trolley
        components['_trolley_damper'] = self._trolley_damper
        components['_pendulum'] = self._pendulum
        components['_pendulum_damper'] = self._pendulum_damper
        components['_force'] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.k: r'Stiffness of a beam showed as a spring stiffness in trolley member',
            self.x: r'Kinematic lateral excitation',
            self.phi: r'Angle of a pendulum',
            self.m_p: r'Mass of pendulum',
            self.m_t: r'Mass of trolley',
            self.g: r'Gravity constant',
            self.F: r'Force',
            self.Omega: r'Excitation frequency',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0, Omega0, k0 = symbols('m_0 l_0 F_0 Omega_0 k_0', positive=True)

        default_data_dict = {
            self.m_t: [S.One * m0 * no for no in range(20, 30)],
            self.m_p: [S.One * m0 * no for no in range(1, 10)],
            self.l: [S.Half * l0 * no for no in range(1, 10)],
            self.F: [S.One * F0 * no for no in range(50, 100)],
            self.Omega: [S.One * Omega0],
            self.g: [S.One * self.g],
            self.k: [S.One * k0 * no for no in range(50, 100)],
            self.x: [self.x]
        }
        return default_data_dict

#     def get_numerical_data(self):

#         default_data_dict = {
#             self.m_p: [10],
#             self.m_f: [50],
# #             self.m_tf: [50],
# #             self.m_pf: [50],
#             self.F: [5],
#             self.flow_coeff: [5],
#             self.t0: [10],
#             self.Omega: [3.14 * 0.5],
#             self.g: [9.81],
#             #self.k: [2]
#         }
#         default_data_dict.update({self.m_t: [10.0*default_data_dict[self.m_p][0]],
#                                   #self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
#                                   #self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)]
#                                  })
#         return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_t0: r'initial mass of the trolley',
            self.m_t: r'mass of the trolley',
            self.m_p0: r'initial mass of the pendulum',
            self.m_p: r'mass of the pendulum',
            self.m_f: r'mass of the transferred liquid to the mass of the system ratio',
            self.k: r'''trolley's equivalent stiffness''',
            self.lam: r'''Rayleigh's proportional damping coefficient''',
            self.l: r'length of the pendulum',
            self.Omega: r'excitation frequency',
            self.omega0: r'natural frequency of the trolley',
            self.F: r'forcing amplitude',
            self.g: r'gravitational constant',
            self.rho: r'flow rate coefficient',
            self.t0: r'mass transfer activation time',
            self.eps: r'mass of the pendulum to the mass of the system ratio',
            self.m_s: r'total mass of the system',
            self.chi: r'members natural frequencies ratio',
            self.psi: r'excitation frequency to the natural frequency of the trolley ratio',
            self.delta: r'''forcing amplitude to structure's equivalent stiffness ratio''',
            self.Omega0: r'natural frequency of the pendulum',
            self.tau: 'dimensionless time',
            self.x: '''trolley's horizontal displacement''',
            self.x.diff(self.ivar): '''trolley's horizontal velocity''',
            self.phi.diff(self.ivar): '''pendulum's angular velocity''',
            self.phi: '''pendulum's angular displacement''',
        }
        return self.sym_desc_dict


    def dimensionless(self):
        
        q=self.q
        m_p = self.m_p
        m_t = self.m_t
        ivar = self.ivar

        inertia_matrix=self.dimensionless_inertia_matrix()
        damping_matrix=self.dimensionless_damping_matrix()
        stiffness_matrix=self.dimensionless_stiffness_matrix()
        external_forces_matrix=self.dimensionless_external_forces_matrix()

        eoms=inertia_matrix*q.diff(ivar,2)+damping_matrix*q.diff(ivar)+stiffness_matrix*q-external_forces_matrix
        eoms_left=inertia_matrix*q.diff(ivar,2)+damping_matrix*q.diff(ivar)+stiffness_matrix*q
        
        
        return eoms


    
    def dimensionless_inertia_matrix(self):
        
        #### nie wiem co z tym, więc zostawię
        inertia_mat=self.inertia_matrix()
        first_inertia_row=inertia_mat.row(0)/inertia_mat[0]
        second_inertia_row=inertia_mat.row(1)/inertia_mat[3]
        dimensionless_inertia_mat=Matrix([[first_inertia_row],[second_inertia_row]])
        #### nie wiem co z tym, więc zostawię
    
        l = self.l
        m_p = self.m_p
        m_t = self.m_t
        inertia_matrix=Matrix([
        [1, l*m_p/(m_p+m_t)],
        [1/l,1]
        ])
        
        return inertia_matrix
        
    def dimensionless_damping_matrix(self):
        
        #### nie wiem co z tym, więc zostawię
        inertia_mat=self.inertia_matrix()
        first_damping_row=self.damping_matrix().row(0)/inertia_mat[0]
        second_damping_row=self.damping_matrix().row(1)/inertia_mat[3]
        dimensionless_damping_mat=Matrix([[first_damping_row],[second_damping_row]])
        #### nie wiem co z tym, więc zostawię
        
        m_p = self.m_p
        m_t = self.m_t
        omega0 = self.omega0
        lam = self.lam
        l = self.l
        chi = self.chi
        
        ivar = self.ivar
        
        damping_matrix=Matrix([
                [lam*omega0, l/omega0*m_p.diff(ivar)/(m_p+m_t)],
                [m_p.diff(ivar)/m_p/l/omega0, (lam*omega0*chi**2+m_p.diff(ivar)/omega0/m_p)]
            ])
        
        return damping_matrix
    
    def dimensionless_stiffness_matrix(self):
        
        #### nie wiem co z tym, więc zostawię
        inertia_mat=self.inertia_matrix()
        first_stiffness_row=self.stiffness_matrix().row(0)/inertia_mat[0]
        second_stiffness_row=self.stiffness_matrix().row(1)/inertia_mat[3]
        dimensionless_stiffness_mat=Matrix([[first_stiffness_row],[second_stiffness_row]])
        #return dimensionless_stiffness_mat.doit()
        #### nie wiem co z tym, więc zostawię

        chi = self.chi
        
        stiffness_matrix=Matrix([
        [1, 0],
        [0, chi**2]
        ])
        
        return stiffness_matrix
    
    def dimensionless_external_forces_matrix(self):
        
        #### nie wiem co z tym, więc zostawię
        inertia_mat=self.inertia_matrix()
        dimensionless_external_forces_mat=self.external_forces()/inertia_mat[0]
        dimensionless_external_forces_mat.doit()
        #### nie wiem co z tym, więc zostawię
        
        delta = self.delta
        psi = self.psi
        ivar = self.ivar
        
        external_forces_matrix=Matrix([
        [delta*sin(psi*ivar)],
        [0]
    ])
        return external_forces_matrix
        
    def mass_transfer_function(self):

        subs_dict={
            self.m_t:self.m_t0 + self.m_f*(1/2-atan(self.rho*self.ivar/self.omega0 - self.rho*self.t0/self.omega0)/pi),
            self.m_p:self.m_p0 + self.m_f*(1/2+atan(self.rho*self.ivar/self.omega0 - self.rho*self.t0/self.omega0)/pi),
        }
        return subs_dict

    def dimensionless_with_transfer(self):
        return self.dimensionless().subs(self.mass_transfer_function()).doit()

    def dimensionless_with_const_mass(self):
        
        m_p_const = Symbol('m_p', positive=True)
        m_t_const = Symbol('m_t', positive=True)
        
        masses_dict = ({self.m_p:m_p_const,self.m_t:m_t_const })
        
        return self.dimensionless().subs(masses_dict).doit().subs({m_t_const:(-m_p_const+m_p_const/self.eps)}) # it's a bit tricky, start with checking it if something doesn't work
    
    
    def system_mass_ratio(self):
        return Eq(Function('epsilon')(self.ivar), self.m_p/self.m_s)
    
    def total_system_mass(self):
        return Eq(self.m_s, self.m_p+self.m_t)
    
    def initial_system_mass_ratio(self):
        return Eq(Symbol('epsilon_0', positive=True), self.m_p0/self.m_s)
    
    def excitation_freq_ratio(self):
        return Eq(self.psi, self.Omega/self.omega0)
    
    def frequency_ratio(self):
        return Eq(self.chi, self.Omega0/self.omega0)
    
    def force_amplitude_ratio(self):
        return Eq(self.delta, self.F/self.k)
    
    def dimensionless_time(self):
        return Eq(self.tau, self.ivar*self.omega0)
    
    
    def get_numerical_data(self):
        
        params_dict = {
                self.m_s:150000000,
                self.chi:1,
                self.psi:1,
                self.lam:0.03,
                self.l:25,
                self.omega0:sqrt(9.81/25),
                self.m_p0:150000,
                #m_f:15000,
                self.m_f:3/1000,
                self.rho:10,
                pi:3.1415,
                self.t0:30,
                self.delta: 300000/(9.81*(150000000)/25)}
        
        return params_dict
    
    def get_analytical_data_before_transfer(self):

        m_p_const = Symbol('m_p', positive=True)
        m_t_const = Symbol('m_t', positive=True)

        params_dict = {
            self.eps:1/1000,
            self.chi:1,
            self.psi:1,
            self.lam:0.03,
            self.l:25,
            self.omega0:sqrt(9.81/25),
            m_p_const:150000,
            m_t_const:150000000-1*150000,
            pi:3.1415,
            self.delta: 300000/(9.81*(150000000)/25)}
        
        return params_dict

    def get_analytical_data_after_transfer(self):

        m_p_const = Symbol('m_p', positive=True)
        m_t_const = Symbol('m_t', positive=True)

        params_dict = {
            self.eps:4/1000,
            self.chi:1,
            self.psi:1,
            self.lam:0.03,
            m_p_const:4*150000,
            m_t_const:150000000-4*150000,
            self.l:25,
            self.omega0:sqrt(9.81/25),
            pi:3.1415,
            self.delta: 300000/(9.81*(150000000)/25)}
        
        return params_dict
    
    def get_numerical_data_for_solving_options(self):

        params_dict = {
            self.chi:1,
            self.psi:1,
            self.lam:0.03,
            self.l:25,
            self.omega0:sqrt(9.81/25),
            self.m_p0:150000,
            self.m_t0:150000000-4*150000,
            self.m_f:3*150000,
            self.rho:10,
            self.t0:90,
            pi:3.1415,
            self.delta: 300000/(9.81*(150000000)/25)}
        
        return params_dict
    
    def _mass_transfer_influence(self, t_span=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame,pd,TimeDataFrame
        
        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.m_f)

        if t_span is None: t_span=np.linspace(0,300,3001)

        masses_dict = {self.m_t0:self.m_s-self.m_p0-self.m_f}
        m_f_dict = {self.m_f:self.m_f*self.m_s}
        eoms_num = self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
#         display(eoms_num)
        a_p=Symbol('a_p')
        dependencies={a_p:self.x.diff().diff() + 25*self.phi.diff().diff()}

        q=self.q
        
        dyn_sys = ODESystem(odes=eoms_num, dvars=q).as_first_ode_linear_system()
        Y=self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(dyn_sys,
                                                      self.m_f,
                                                      [1/10000, 5/10000, 1/1000, 3/1000, 5/1000],
                                                      reference_data = params_dict,
                                                      coordinates=Y,
                                                      index=pd.Index(t_span,name=Symbol('t'))
                                                     )

        #sym_wyn = num_df.reindex(pd.Index(t_span,name=t))#.perform_simulations(model_level_name=0)
        sym_wyn = num_df.perform_simulations(model_level_name=0, dependencies=dependencies)


        result = TimeDataFrame(sym_wyn)
        result.index.name = Symbol('tau')
        
        return result.droplevel(0, axis=1)
    
    def _skyscraper_stiffness_influence(self, t_span=None):
        from ...solvers.linear import ODESystem
        from ...utilities.adaptable import NumericalAnalysisDataFrame,pd,TimeDataFrame
        
        params_dict = {**self.get_numerical_data()}
        params_dict.pop(self.chi)

        if t_span is None: t_span=np.linspace(0,300,3001)

        masses_dict = {self.m_t0:self.m_s-self.m_p0-self.m_f}
        m_f_dict = {self.m_f:self.m_f*self.m_s}
        eoms_num = self.dimensionless_with_transfer().subs(masses_dict).subs(m_f_dict).doit()
        a_p=Symbol('a_p')
        dependencies={a_p:self.x.diff().diff() + 25*self.phi.diff().diff()}

        q=self.q
        
        dyn_sys = ODESystem(odes=eoms_num, dvars=q).as_first_ode_linear_system()
        Y=self._coords_with_acceleration_list + [a_p]
        num_df = NumericalAnalysisDataFrame.from_model(dyn_sys,
                                                      self.chi,
                                                      [0.8, 0.9, 1, 1.1, 1.2],
                                                      reference_data = params_dict,
                                                      coordinates=Y,
                                                      index=pd.Index(t_span,name=Symbol('t'))
                                                     )

        #sym_wyn = num_df.reindex(pd.Index(t_span,name=t))#.perform_simulations(model_level_name=0)
        sym_wyn = num_df.perform_simulations(model_level_name=0, dependencies=dependencies)


        result = TimeDataFrame(sym_wyn)
        result.index.name = Symbol('tau')

        return result.droplevel(0, axis=1)

    def _solving_options_analysis(self, t_span=None):
        from ...solvers.linear import ODESystem, FirstOrderLinearODESystemWithHarmonics
        from ...solvers.tools import TwoSystemsSimulation
        from ...utilities.adaptable import NumericalAnalysisDataFrame,pd,TimeDataFrame

        params_dict = {**self.get_numerical_data()}

        if t_span is None: t_span=np.linspace(0,300,3001)

        #ANALITYCZNIE
        eoms1=self.dimensionless_with_const_mass().subs(self.get_analytical_data_before_transfer()).doit()

        ode1=ODESystem(odes=eoms1, dvars=self.q)

        eoms2=self.dimensionless_with_const_mass().subs(self.get_analytical_data_after_transfer()).doit()
        ode2=ODESystem(odes=eoms2, dvars=self.q)

        odes=[ode1, ode2]

        symki=TwoSystemsSimulation(odes=(ode1,ode2),t_span=(np.linspace(0,89.9,900),np.linspace(90,300,2101)))
        wyniki=symki.compute_solution(t_span=(np.linspace(0,89.9,900),np.linspace(90,300,2101))).apply(np.real)

        analytical_sym_wyn=TimeDataFrame(data={('A-A ','X'): wyniki[self.x]})
        analytical_sym_wyn.index.name = Symbol('tau')

        #NUMERYCZNIE
        masses_dict = {self.m_t0:self.m_s-self.m_p0-self.m_f}
        eoms_num = self.dimensionless_with_transfer().subs(self.get_numerical_data_for_solving_options()).doit()
#         display(eoms_num)
        ode_2=ODESystem(odes=eoms_num, dvars=self.q).numerized()
        num_sym_wyn=ode_2.compute_solution(t_span=t_span, ic_list=[0.0, 0.0, 0.0, 0.0]).apply(np.real)
        
        numerical_sym_wyn=TimeDataFrame(data={('N-N-N','X'): num_sym_wyn[self.x]})
        numerical_sym_wyn.index.name = Symbol('tau')

        #ANALITYCZNIE + NUMERYCZNIE
        fode_1=ode1.as_first_ode_linear_system()
        fodes_1=FirstOrderLinearODESystemWithHarmonics(fode_1,fode_1.dvars)
        ana_sol_1=fodes_1.solution
        ics_list=[0.0,0.0,0.0,0.0]
        ana_sol_sym_1=ana_sol_1.with_ics(ics_list)
        t_span_1=np.linspace(0,84.9,850)
        ana_sym_wyn_1=ana_sol_sym_1.numerized().compute_solution(t_span=t_span_1).apply(np.real)
        ics_list_new=list(ana_sym_wyn_1.iloc[-1,0:len(ics_list)])

        ##Symulacja fazy 2.
        t_span_2=np.linspace(85,95,101)
        num_sym_wyn_2=ode_2.compute_solution(t_span=t_span_2, ic_list=ics_list_new).apply(np.real)
        ics_list_the_newest=list(num_sym_wyn_2.iloc[-1,0:len(ics_list)])

        ##Symulacja fazy 3.
        fode_3=ode2.as_first_ode_linear_system()
        fodes_3=FirstOrderLinearODESystemWithHarmonics(fode_3,fode_3.dvars)
        ana_sol_3=fodes_3.solution
        ana_sol_sym_3=ana_sol_3.with_ics(ics_list_the_newest, ivar0=95.1)
        t_span_3=np.linspace(95.1,300,2050)
        ana_sym_wyn_3=ana_sol_sym_3.numerized().compute_solution(t_span=t_span_3).apply(np.real)
        
        results_ANA=pd.concat([ana_sym_wyn_1, num_sym_wyn_2, ana_sym_wyn_3])
        
        analytical_numerical_sym_wyn=TimeDataFrame(data={('A-N-A','X'): results_ANA[self.x]})
        analytical_numerical_sym_wyn.index.name = Symbol('tau')

        results=pd.concat([analytical_numerical_sym_wyn, numerical_sym_wyn, analytical_sym_wyn])
        return results
        


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU mdof.py
class DoubleTrolleyWithNonlinearSprings(ComposedSystem):
    scheme_name = 'MDOF_Double_Trolley_With_Springs.PNG'
    real_name = 'sdof_nonlinspring_trolleys_real.PNG'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 k_cc=Symbol('k_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 mu=Symbol('mu', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_c = k_c
        self.k_cc = k_cc
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R
        self.mu = mu

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_cl*(1-mu*x_l**2), pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))


        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_cc*(1-mu*(x_r - x_l)**2), pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_cr*(1-mu*x_r**2), pos1=x_r, qs=[x_r]) 
                          + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cc: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU mdof.py
class DampedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'MDOF_Damped_Trolleys_With_Springs.PNG'
    real_name = 'two_damped_trolleys.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m = m
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))

        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) 
                          + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU mdof.py
class DoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'DoubleTrolley_real.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m_w1=Symbol('m_w1', positive=True),
                 m_w2=Symbol('m_w2', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 lam=Symbol('lambda', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.x_l = x_l
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x=x
        self.m1 = m1
        self.m2 = m2
        self.m_w1 = m_w1
        self.m_w2 = m_w2
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.R = R
        self.lam = lam

        self.trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m_w1, x_l/2, qs = [x_l]) + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]) + MaterialPoint(m_w1, x_l/2, qs = [x_l])
                          + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]))

        self.trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r])
                          + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]))

        system = self.trolley_1 + self.trolley_2
        super().__init__(system(qs),**kwargs)
        
        

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            
            
            self.lam:[0],


        }

        return default_data_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU mdof.py
class ForcedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'mdof_three_trolleys.PNG'
    real_name = 'three_carriages.PNG'

    def __init__(self,
                 dum=Symbol('dumy'),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_12=Symbol('k_12', positive=True),
                 k_c12=Symbol('k_c12', positive=True),
                 k_23=Symbol('k_23', positive=True),
                 k_c23=Symbol('k_c23', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 F=Symbol('F_0', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 
                 x_l=dynamicsymbols('x_l'),
                 x_c=dynamicsymbols('x_c'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x_3=dynamicsymbols('x_3'),
                 
                 qs=dynamicsymbols('x_l x_c x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_12 = k_12
        self.k_c12 = k_c12
        self.k_23 = k_23
        self.k_c23 = k_c23
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_c = x_c
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x_3 = x_3
        self.Omega = Omega

        self.Trolley_1 = MaterialPoint(m1, x_l, qs=[x_l]) + Spring(
            k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[
                x_l
            ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(
                -F * cos(Omega * ivar), pos1=x_l, qs=[x_l])
        
        
        
        self.Trolley_2 = MaterialPoint(m2, x_c, qs=qs) + Spring(
            k_12, x_l, x_c,qs=qs) + Spring(k_12, x_l, x_c,qs=qs) + Spring(
                k_23, x_c, x_r,qs=qs) + Spring(k_23, x_c, x_r,qs=[x_c,x_l]) + Spring(
                    k_c12, x_l, x_c,qs=qs) + Spring(k_c23, x_c, x_r,qs=qs)
        self.Trolley_3 = MaterialPoint(m3, x_r, qs=qs) + Spring(
            k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[
                x_r
            ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(
                -2 * F * cos(Omega * ivar), pos1=x_r, qs=[x_r])

        system = self.Trolley_1 + self.Trolley_2 + self.Trolley_3
        super().__init__(system(x_l,x_c),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.k_l: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_cl: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_12: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2, 0],
            self.x_r: [self.x_2, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        print(parameters_dict)
        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict
    

#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU systems.py
class MDoFDoubleTrolleyWithNonlinearSprings(ComposedSystem):
    scheme_name = 'MDOF_Double_Trolley_With_Springs.PNG'
    real_name = 'sdof_nonlinspring_trolleys_real.PNG'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 k_cc=Symbol('k_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 mu=Symbol('mu', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_c = k_c
        self.k_cc = k_cc
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R
        self.mu = mu

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_cl*(1-mu*x_l**2), pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))
        
        
       
        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_cc*(1-mu*(x_r - x_l)**2), pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_cr*(1-mu*x_r**2), pos1=x_r, qs=[x_r]) 
                          + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cc: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU systems.py
class MDoFDampedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'MDOF_Damped_Trolleys_With_Springs.PNG'
    real_name = 'two_damped_trolleys.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m = m
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))

        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) 
                          + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU systems.py
class SDoFDoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'

    def __init__(self,
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 m_w1=Symbol('m_w1', positive=True),
                 m_w2=Symbol('m_w2', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 R=Symbol('R', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 lam=Symbol('lambda', positive=True),

                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m_w1 = m_w1
        self.m_w2 = m_w2
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.R = R
        self.lam = lam

        self.trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m_w1, x_l/2, qs = [x_l]) + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]) + MaterialPoint(m_w1, x_l/2, qs = [x_l])
                          + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]))

        self.trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r])
                          + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]))

        system = self.trolley_1 + self.trolley_2
        super().__init__(system(qs),**kwargs)

        
    def get_default_data(self):

        m0, k0, l0, lam = symbols('m_0 k_0 l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


#TODO
#DO SPRAWDZENIA - - PRZENIESIONO Z MODUŁU systems.py
class MDoFDoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'DoubleTrolley_real.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m_w1=Symbol('m_w1', positive=True),
                 m_w2=Symbol('m_w2', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 lam=Symbol('lambda', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.x_l = x_l
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x=x
        self.m1 = m1
        self.m2 = m2
        self.m_w1 = m_w1
        self.m_w2 = m_w2
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.R = R
        self.lam = lam

        self.trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m_w1, x_l/2, qs = [x_l]) + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]) + MaterialPoint(m_w1, x_l/2, qs = [x_l])
                          + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]))

        self.trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r])
                          + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]))

        system = self.trolley_1 + self.trolley_2
        super().__init__(system(qs),**kwargs)
        
        

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            
            
            self.lam:[0],


        }

        return default_data_dict


#TODO
#DO SPRAWDZENIA - PRZENIESIONO Z MODUŁU systems.py
class MDoFForcedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'mdof_three_trolleys.PNG'
    real_name = 'three_carriages.PNG'

    def __init__(self,

                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_12=Symbol('k_12', positive=True),
                 k_c12=Symbol('k_c12', positive=True),
                 k_23=Symbol('k_23', positive=True),
                 k_c23=Symbol('k_c23', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 F=Symbol('F', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 
                 x_l=dynamicsymbols('x_l'),
                 x_c=dynamicsymbols('x_c'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x_3=dynamicsymbols('x_3'),
                 
                 qs=dynamicsymbols('x_l x_c x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_12 = k_12
        self.k_c12 = k_c12
        self.k_23 = k_23
        self.k_c23 = k_c23
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_c = x_c
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x_3 = x_3
        self.Omega = Omega

        self.Trolley_1 = MaterialPoint(m1, x_l, qs=[x_l]) + Spring(
            k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[
                x_l
            ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(
                -F * cos(Omega * ivar), pos1=x_l, qs=[x_l])
        
        
        
        self.Trolley_2 = MaterialPoint(m2, x_c, qs=qs) + Spring(
            k_12, x_l, x_c,qs=qs) + Spring(k_12, x_l, x_c,qs=qs) + Spring(
                k_23, x_c, x_r,qs=qs) + Spring(k_23, x_c, x_r,qs=[x_c,x_l]) + Spring(
                    k_c12, x_l, x_c,qs=qs) + Spring(k_c23, x_c, x_r,qs=qs)
        self.Trolley_3 = MaterialPoint(m3, x_r, qs=qs) + Spring(
            k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[
                x_r
            ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(
                -2 * F * cos(Omega * ivar), pos1=x_r, qs=[x_r])

        system = self.Trolley_1 + self.Trolley_2 + self.Trolley_3
        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2, 0],
            self.x_r: [self.x_2, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict


#TODO
#ZAPYTAC SIE BOGUSIA CO Z TYM - PRZENIESIONO Z MODUŁU systems.py
class MDoFTMD(ComposedSystem):
    scheme_name = 'mdof_tmd.png'
    real_name = 'mdof_tmd_real.png'
    """
    Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            me = Mass
                -Mass of tuned mass damper

            k = spring coefficient
                -value of spring coefficient

            ke = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            F = Force
                -Trolley's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        >>> t = symbols('t')
        >>> m, me, k, ke, F = symbols('m, m_e, k, k_e, F')

    """

    def __init__(self,
                 m=Symbol('m', positive=True),
                 me=Symbol('m_e', positive=True),
                 k=Symbol('k', positive=True),
                 ke=Symbol('k_e', positive=True),
                 F=Symbol('F', positive=True),
                 xe=dynamicsymbols('x_e'),
                 xb=dynamicsymbols('x_b'),
                 angle=dynamicsymbols('Omega'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.m = m
        self.me = me
        self.k = k
        self.ke = ke
        self.F = F

        self.MaterialPoint_1 = MaterialPoint(m, pos1=xe, qs=[xe])
        self.MaterialPoint_2 = MaterialPoint(me, pos1=xb, qs=[xb])
        self.Spring_1 = Spring(k, pos1=xe, qs=[xe])
        self.Spring_2 = Spring(ke, pos1=xe, pos2=xb, qs=[xe, xb])
        self.Force = Force(F * sin(angle * ivar), pos1=xe, qs=[xe])

        system = self.Spring_1 + self.Spring_2 + \
            self.MaterialPoint_1 + self.MaterialPoint_2 + self.Force
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of main object',
            self.me: r'Mass of the TMD',
            self.k: r'Stiffness coefficient',
            self.ke: r'Stiffness coefficient',
            self.F: r'Force',
        }
        return self.sym_desc_dict

