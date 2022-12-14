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

    @property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem(self.m_t, self.k, self.x, self.ivar)(label='Trolley')
        self._pendulum = PendulumKinematicExct(self.l, self.m_p, self.g, self.phi, self.x, self.ivar)(label='Pendulum')
        self._force=Force(2*self.F+self.F*sin(self.Omega*self.ivar)+3*self.F*cos(self.Omega*self.ivar), pos1=self.x, qs=[self.x, self.phi])(label='Force')

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
            self.m_t: [no for no in range(20, 30)],
            self.m_p: [no for no in range(1, 10)],
            self.l: [1/2 * no for no in range(1, 10)],
            self.F: [no for no in range(50, 100)],
            self.Omega: [3.14 * no for no in range(1,6)],
            self.k: [no for no in range(50, 100)],
            self.g: [9.81]
        }
        return default_data_dict

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
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 k=None,
                 Omega=None,
                 F=None,
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
    
class ForcedDampedTrolleysWithSprings(NonlinearComposedSystem): ### 4 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    k=Symbol('k', positive=True)
    c=Symbol('c', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    G=Symbol('G', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
   
    def __init__(self,
                 m_1=None,
                 m_2=None,
                 k=None,
                 c=None,
                 Omega=None,
                 F=None,
                 G=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_1 is not None: self.m_1 = m_1
        if m_2 is not None: self.m_2 = m_2
        if k is not None: self.k= k
        if c is not None: self.c= c
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
        self._spring1 = Spring(self.k, pos1 = self.x_1, qs=[self.x_1])
        self._damper1 = Damper(self.c, pos1 = self.x_1, qs=[self.x_1])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.G, pos1 = self.x_1, qs=[self.x_1])

        self._spring12 = Spring(self.k, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        self._damper12 = Damper(self.c, pos1 = self.x_1, pos2 = self.x_2, qs=[self.x_1, self.x_2])
        
        self._trolley2 = MaterialPoint(self.m_2, self.x_2, qs=[self.x_2])
        self._spring2 = Spring(self.k, pos1 = self.x_2, qs=[self.x_2])
        self._damper2 = Damper(self.c, pos1 = self.x_2, qs=[self.x_2])


        components['trolley_1'] = self._trolley1
        components['spring_1'] = self._spring1
        components['damper_1'] = self._damper1
        components['spring_12'] = self._spring12
        components['damper_12'] = self._damper12
        components['trolley_2'] = self._trolley2
        components['spring_2'] = self._spring2
        components['damper_2'] = self._damper2
        components['force'] = self._force
        
        return components
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_1 : [100],
            self.m_2 : [200],
            self.k : [50],
            self.c : [100],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
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
    
class ForcedTrolleyWithSpring(ComposedSystem): ### 1 ODE
    
    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    m=Symbol('m', positive=True)
    k=Symbol('k', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F',positive=True)
    G=Symbol('G',positive=True)
    x=dynamicsymbols('x')
   
    def __init__(self,
                 m=None,
                 k=None,
                 Omega=None,
                 F=None,
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

    def get_numerical_data(self):

        default_data_dict = {
            self.m : [100],
            self.k : [50],
            self.Omega : [0.5 * 3.14, 1 * 3.14, 2 * 3.14, 4 * 3.14],
            self.F : [0.5 * 100, 1 * 100, 2 * 100, 4 * 100]
        }
        default_data_dict.update({self.G: [4*default_data_dict[self.F][0]*cos(0.5*default_data_dict[self.Omega][0]*self.ivar) , default_data_dict[self.F][0]*cos(0.75*default_data_dict[self.Omega][0]*self.ivar)**2 , 1.5*default_data_dict[self.F][0]*2*cos(1.25*default_data_dict[self.Omega][0]*self.ivar) , 3*default_data_dict[self.F][0]*cos(2*default_data_dict[self.Omega][0]*self.ivar)**2]})

        return default_data_dict
    
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

    @property
    def components(self):
        components = {}

        self._trolley = SpringMassSystem((self.m_t + self.m_tf), self.k, self.x, self.ivar)(label='Trolley')
        self._pendulum = PendulumKinematicExct(self.l, (self.m_p + self.m_pf), self.g, self.phi, self.x, self.ivar)(label='Pendulum')
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
    
#     def get_numerical_data(self):

#         default_data_dict = {
#             self.m_t: [no for no in range(20, 30)],
#             self.m_p: [no for no in range(1, 10)],
#             self.l: [1/2 * no for no in range(1, 10)],
#             self.F: [0.1*no for no in range(50, 100)],
#             self.Omega: [3.14 * no for no in range(1,2)],
#             self.k: [no for no in range(50, 100)],
#             self.g: [9.81]
#         }
#         return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_p: [10],
            self.m_f: [50],
            self.m_tf: [50],
            self.m_pf: [50],
            self.F: [5],
            self.flow_coeff: [5],
            self.t0: [10],
            self.Omega: [3.14 * 0.5],
            self.g: [9.81],
            self.k: [2]
        }
        default_data_dict.update({self.m_t: [10.0*default_data_dict[self.m_p][0]],
                                  self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
                                  self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)]})
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

    m_p = Symbol('m_p', positive=True)
    m_p0 = Symbol('m_p0', positive=True)
    m_pf = Symbol('m_pf', positive=True)
    
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
                 alpha=None,
                 beta=None,
                 b=None,
                 Omega=None,
                 rayleigh_damping_matrix=None,
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
        if rayleigh_damping_matrix is not None: self.rayleigh_damping_matrix = rayleigh_damping_matrix
        if t0 is not None: self.t0 = t0
        if g is not None: self.g = g
        if b is not None: self.b = b
        if phi is not None: self.phi = phi
        if x is not None: self.x = x
        if k is not None: self.k = k
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        self.ivar = ivar
        self.trans_expr = ((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        self.alpha = self.b
        self.beta = self.b/2
        #self.rayleigh_damping_matrix = self.alpha*VariableMassTrolleyWithPendulum().inertia_matrix() + self.beta*VariableMassTrolleyWithPendulum().stiffness_matrix()
        self._init_from_components(**kwargs)



    @property
    def components(self):
        components = {}
        
        m_tf_eq = self.m_f*((S.One/2-atan(self.flow_coeff*(self.ivar-self.t0))/pi))
        m_pf_eq = self.m_f*((S.One/2+atan(self.flow_coeff*(self.ivar-self.t0))/pi))

        self._trolley = VariableMassTrolleyWithPendulum(m_tf = m_tf_eq, m_pf = m_pf_eq)(label='Trolley')
        
        damping_matrix = self.alpha*self._trolley.inertia_matrix() + self.beta*self._trolley.stiffness_matrix()
        
        self._trolley_damping = Damper((damping_matrix[0]+damping_matrix[2]), pos1 = self.x , qs = [self.x,self.phi])(label='Trolley dapming')
        self._pendulum_damping = Damper((damping_matrix[1]+damping_matrix[3]), pos1 = self.phi , qs = [self.x,self.phi])(label='Pendulum dapming')

        components['_trolley'] = self._trolley
        components['_trolley_damping'] = self._trolley_damping
        components['_pendulum_damping'] = self._pendulum_damping

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
#             self.m_t: [no for no in range(20, 30)],
#             self.m_p: [no for no in range(1, 10)],
#             self.l: [1/2 * no for no in range(1, 10)],
#             self.F: [0.1*no for no in range(50, 100)],
#             self.Omega: [3.14 * no for no in range(1,2)],
#             self.k: [no for no in range(50, 100)],
#             self.g: [9.81]
#         }
#         return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_p: [10],
            self.m_f: [50],
            self.m_tf: [50],
            self.m_pf: [50],
            self.alpha: [0.8],
            self.beta: [0.6],
            self.b: [0.02],
            self.F: [250],
            self.flow_coeff: [10],
            self.t0: [30],
            self.Omega: [3.14 * 0.5],
            self.g: [9.81],
            self.k: [2]
        }
        default_data_dict.update({self.m_t: [10.0*default_data_dict[self.m_p][0]],
                                  self.l: [0.99*default_data_dict[self.g][0]/(default_data_dict[self.Omega][0]*default_data_dict[self.Omega][0])],
                                  self.k: [2.0*default_data_dict[self.m_p][0]*default_data_dict[self.g][0]/(default_data_dict[self.g][0]/default_data_dict[self.Omega][0]**2)]})
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
            self.x: 'przemieszczenie wózka',
            self.phi: 'kąt wychylenia wahadła',
        }
        return self.sym_desc_dict

dyn_sys_RD = VariableMassTrolleyWithPendulumRayleighDamping()


    
    