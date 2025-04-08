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
from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST
from .tmd import TunedMassDamperRelativeMotion
from ...utilities.components.mech import en as mech_comp
from .trolley import ThreeSprings, ForcedTrolleyWithThreeSprings
from .disk import RollingDisk

from sympy.physics import units
ureg = units

from functools import cached_property, lru_cache


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



class ForcedNonLinTrolleysWithSprings(ComposedSystem): ### 3 ODE
    
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

        m0, k0, F0 =  self.m0, self.k0, self.F0
        
        default_data_dict = {
            self.m_1 : [S.One * m0 * no for no in range(10,20)],
            self.m_2 : [S.One * m0 * no for no in range(20,30)],
            self.k : [S.One * k0 * no for no in range(30,50)],
            self.F : [S.One * F0 * no for no in range(1,10)]
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
        *REPORT_COMPONENTS_LIST,
#         mech_comp.TitlePageComponent,
#         mech_comp.SchemeComponent,
#         mech_comp.ExemplaryPictureComponent,
#         mech_comp.KineticEnergyComponent,
#         mech_comp.PotentialEnergyComponent,
#         mech_comp.LagrangianComponent,
#         mech_comp.GoverningEquationComponent,
#         mech_comp.FundamentalMatrixComponent,
#         mech_comp.GeneralSolutionComponent,
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
    x_l=dynamicsymbols('x_l')
    x_r=dynamicsymbols('x_r')

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
                 x_l=None,
                 x_r=None,
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
        if x_l is not None: self.x_l = x_l
        if x_r is not None: self.x_r = x_r
        self.ivar=ivar

        self.qs = [self.x_l, self.x_r]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._trolley1 = MaterialPoint(self.m_1, self.x_l, qs=[self.x_l])
        self._wheel11 = MaterialPoint(self.m, self.x_l, qs = [self.x_l])
        self._wheel11_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_l/self.R, qs=[self.x_l])
        self._wheel12= MaterialPoint(self.m, self.x_l, qs = [self.x_l])
        self._wheel12_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_l/self.R, qs=[self.x_l])
        
        self._spring1_top = Spring(self.k_l, pos1 = self.x_l, qs=[self.x_l])
        self._spring1_bottom = Spring(self.k_l, pos1 = self.x_l, qs=[self.x_l])
        
        self._damper1 = Damper(self.c_l, pos1 = self.x_l, qs=[self.x_l])
        self._force = Force(self.F*sin(self.Omega*self.ivar) + self.F, pos1 = self.x_l, qs=[self.x_l])

        self._spring12 = Spring(2*self.k_c, pos1 = self.x_l, pos2 = self.x_r, qs=[self.x_l, self.x_r])
        self._damper12 = Damper(self.c_c, pos1 = self.x_l, pos2 = self.x_r, qs=[self.x_l, self.x_r])
        self._wheel21 = MaterialPoint(self.m, self.x_r, qs = [self.x_r])
        self._wheel21_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_r/self.R, qs=[self.x_r])
        self._wheel22 = MaterialPoint(self.m, self.x_r, qs = [self.x_r])
        self._wheel22_disk = Disk(I=S.Half*self.m*self.R**2, pos1=self.x_r/self.R, qs=[self.x_r])
        self._trolley2 = MaterialPoint(self.m_2, self.x_r, qs=[self.x_r])
        self._spring2 = Spring(2*self.k_r, pos1 = self.x_r, qs=[self.x_r])
        self._damper2 = Damper(self.c_r, pos1 = self.x_r, qs=[self.x_r])

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

    def unit_dict(self):
        units_dict = {
        self.m: ureg.kilogram,  # Mass of the trolley
        self.Omega: ureg.radian / ureg.second,  # Frequency of the force (rad/s)
        self.ivar: ureg.second,  # Time (seconds)
        self.k_l: ureg.newton / ureg.meter,  # Spring constant (units depend on spring type)
        self.k_r: ureg.newton / ureg.meter,  # Spring constant (units depend on spring type)
        self.k_c: ureg.newton / ureg.meter,  # Spring constant (units depend on spring type)
        self.x_l: ureg.meter,
        self.x_r: ureg.meter,
        self.x_l.diff(self.ivar): ureg.meter/ureg.second,
        self.x_r.diff(self.ivar): ureg.meter/ureg.second,
        self.x_l.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
        self.x_r.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
        self.c_l: ureg.newton*ureg.second / ureg.meter,  # Damping coefficient
        self.c_r: ureg.newton*ureg.second / ureg.meter,  # Damping coefficient
        self.c_c: ureg.newton*ureg.second / ureg.meter,  # Damping coefficient
        self.F: ureg.newton,  # Force amplitude
        }

        return units_dict


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



#issue 204 zrobione karolinka + bogus
class DoubleTrolleyWithNonlinearSprings(ComposedSystem):
    scheme_name = 'MDOF_Double_Trolley_With_Springs.PNG'
    real_name = 'sdof_nonlinspring_trolleys_real.PNG'

    R=Symbol('R', positive=True)
    m=Symbol('m', positive=True)
    m1=Symbol('m_1', positive=True)
    m2=Symbol('m_2', positive=True)
    k_l=Symbol('k_l', positive=True)
    k_cl=Symbol('k_cl', positive=True)
    k_c=Symbol('k_c', positive=True)
    k_cc=Symbol('k_cc', positive=True)
    k_r=Symbol('k_r', positive=True)
    k_cr=Symbol('k_cr', positive=True)
    mu=Symbol('mu', positive=True)
    x_l=dynamicsymbols('x_l')
    x_r=dynamicsymbols('x_r')
    x=dynamicsymbols('x')
    qs=dynamicsymbols('x_l x_r')
    ivar=Symbol('t')



    def __init__(self,
                 R=None,
                 m=None,
                 m1=None,
                 m2=None,
                 k_l=None,
                 k_cl=None,
                 k_c=None,
                 k_cc=None,
                 k_r=None,
                 k_cr=None,
                 mu=None,
                 x_l=None,
                 x_r=None,
                 x=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if m1 is not None: self.m1 = m1
        if m2 is not None: self.m2 = m2
        if k_l is not None: self.k_l = k_l
        if k_cl is not None: self.k_cl = k_cl
        if k_c is not None: self.k_c = k_c
        if k_cc is not None: self.k_cc = k_cc
        if k_r is not None: self.k_r = k_r
        if k_cr is not None: self.k_cr = k_cr
        if x_l is not None: self.x_l = x_l
        if x_r is not None: self.x_r = x_r
        if x is not None: self.x = x
        if R is not None: self.R = R
        if mu is not None: self.mu = mu
        if qs is not None: self.qs=qs

        self.ivar=ivar
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._trolley_1 = ForcedTrolleyWithThreeSprings(self.m1, self.k_l, self.k_cl*(1-self.mu*self.x_l**2), 0, z=self.x_l, qs=[self.x_l])
        self._three_springs=ThreeSprings(self.k_c, self.k_cc*(1-self.mu*(self.x_r-self.x_l)**2),pos1=self.x_l,pos2=self.x_r, qs=self.qs)
        self._rolling_disk1=RollingDisk(self.m,self.R,self.x_l/2,qs=[self.x_l])
        self._rolling_disk2=RollingDisk(self.m,self.R,self.x_l/2,qs=[self.x_l])
        self._trolley_2= ForcedTrolleyWithThreeSprings(self.m2, self.k_r, self.k_cr*(1-self.mu*self.x_r**2), 0, z=self.x_r, qs=[self.x_r])
        self._rolling_disk3=RollingDisk(self.m,self.R,self.x_r/2,qs=[self.x_r])
        self._rolling_disk4=RollingDisk(self.m,self.R,self.x_r/2,qs=[self.x_r])

        components['_trolley_1'] = self._trolley_1
        components['_three_springs'] = self._three_springs
        components['_trolley_2'] = self._trolley_2
        components['_rolling_disk1'] = self._rolling_disk1
        components['_rolling_disk2'] = self._rolling_disk2
        components['_rolling_disk3'] = self._rolling_disk3
        components['_rolling_disk4'] = self._rolling_disk4

        return components



#DONE(#206)
class DoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'DoubleTrolley_real.png'
    
    R=Symbol('R', positive=True)
    m_w1=Symbol('m_w1', positive=True)
    m_w2=Symbol('m_w2', positive=True)
    m1=Symbol('m_1', positive=True)
    m2=Symbol('m_2', positive=True)
    k_l=Symbol('k_l', positive=True)
    c_cl=Symbol('c_cl', positive=True)
    k_c=Symbol('k_c', positive=True)
    c_cc=Symbol('c_cc', positive=True)
    k_r=Symbol('k_r', positive=True)
    c_cr=Symbol('c_cr', positive=True)
    lam=Symbol('lambda', positive=True)
    x_l=dynamicsymbols('x_l')
    x_r=dynamicsymbols('x_r')

    qs=[x_l, x_r]
    
    def __init__(self,
                 R=None,
                 m_w1=None,
                 m_w2=None,
                 m1=None,
                 m2=None,
                 k_l=None,
                 c_cl=None,
                 k_c=None,
                 c_cc=None,
                 k_r=None,
                 c_cr=None,
                 lam=None,
                 x_l=None,
                 x_r=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if x_l is not None: self.x_l = x_l
        if x_r is not None: self.x_r = x_r
        if m1 is not None: self.m1 = m1
        if m2 is not None: self.m2 = m2
        if m_w1 is not None: self.m_w1 = m_w1
        if m_w2 is not None: self.m_w2 = m_w2
        if k_l is not None: self.k_l = k_l
        if c_cl is not None: self.c_cl = c_cl
        if k_c is not None: self.k_c = k_c
        if c_cc is not None: self.c_cc = c_cc
        if k_r is not None: self.k_r = k_r
        if c_cr is not None: self.c_cr = c_cr
        if R is not None: self.R = R
        if lam is not None: self.lam = lam
        if qs is not None: self.qs=qs
        
        self.ivar=ivar
        self._init_from_components(**kwargs)

        
    @property
    def components(self):

        components = {}


        self._trolley_1 = MaterialPoint(self.m1, self.x_l, qs=self.qs)
        self._rolling_disk1=RollingDisk(self.m_w1,self.R,self.x_l/2, qs=self.qs)
        self._rolling_disk2=RollingDisk(self.m_w1,self.R,self.x_l/2, qs=self.qs)
        self._k_l_spring_1 = Spring(self.k_l, pos1=self.x_l,  qs=self.qs)
        self._k_l_spring_2 = Spring(self.k_l, pos1=self.x_l,  qs=self.qs)
        self._c_cl_damper = Damper(self.c_cl, pos1=self.x_l,  qs=self.qs)
        
        
        self._trolley_2 = MaterialPoint(self.m1, self.x_r, qs=self.qs)
        self._rolling_disk3=RollingDisk(self.m_w2,self.R,self.x_r/2, qs=self.qs)
        self._rolling_disk4=RollingDisk(self.m_w2,self.R,self.x_r/2, qs=self.qs)
        self._k_r_spring_1 = Spring(self.k_r, pos1=self.x_r,  qs=self.qs)
        self._k_r_spring_2 = Spring(self.k_r, pos1=self.x_r,  qs=self.qs)
        self._c_cr_damper = Damper(self.c_cr, pos1=self.x_r,  qs=self.qs)
        
        self._c_cc_damper = Damper(self.c_cc, pos1=self.x_l, pos2=self.x_r, qs=self.qs)
        self._k_c_spring_1 = Spring(self.k_c, pos1=self.x_l, pos2=self.x_r, qs=self.qs)
        self._k_c_spring_2 = Spring(self.k_c, pos1=self.x_l, pos2=self.x_r, qs=self.qs)


        components['_trolley_1'] = self._trolley_1
        
        components['_trolley_2'] = self._trolley_2
        components['_rolling_disk1'] = self._rolling_disk1
        components['_rolling_disk2'] = self._rolling_disk2
        
        components['_k_l_spring_1'] = self._k_l_spring_1
        components['_k_l_spring_2'] = self._k_l_spring_2
        components['_c_cl_damper'] = self._c_cl_damper
        

        components['_rolling_disk3'] = self._rolling_disk3
        components['_rolling_disk4'] = self._rolling_disk4

        components['_k_r_spring_1'] = self._k_r_spring_1
        components['_k_r_spring_2'] = self._k_r_spring_2
        components['_c_cr_damper'] = self._c_cr_damper
        
        components['_k_c_spring_1'] = self._k_c_spring_1
        components['_k_c_spring_2'] = self._k_c_spring_2
        components['_c_cc_damper'] = self._c_cc_damper
        
        
        return components

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

    def get_numerical_data(self):
        
        numerical_data_dict = {
            self.m1: [1],
            self.m2: [1],
            self.m_w1: [0.1],
            self.m_w2: [0.5],
            self.k_l: [50],
            self.k_c: [70],
            self.k_r: [50],

            self.c_cr: [0.3],
            self.c_cc: [0.1],
            self.c_cl: [0.3],
            self.R: [0.05],
        }

        return numerical_data_dict

    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley Mass',
            self.m2: r'Trolley Mass',
            self.m_w1: r'Wheels 1 Mass',
            self.m_w2: r'Wheels 2 Mass',
            self.k_l: r'Spring Stiffness',
            self.k_c: r'Spring Stiffness',
            self.k_r: r'Spring Stiffness',

            self.c_cr: r'Damper damping',
            self.c_cc: r'Damper damping',
            self.c_cl: r'Damper damping',
            self.R: r'Wheels radius',

        }
        return self.sym_desc_dict
    
    def unit_dict(self):

        unit_dict = {

            self.ivar: ureg.second,

            self.m1: ureg.kilogram,
            self.m2: ureg.kilogram,
            self.m_w1: ureg.kilogram,
            self.m_w2: ureg.kilogram,
            self.k_l: ureg.newton/ureg.meter,
            self.k_c: ureg.newton/ureg.meter,
            self.k_r: ureg.newton/ureg.meter,

            self.c_cr: ureg.newton/(ureg.meter/ureg.second),
            self.c_cc: ureg.newton/(ureg.meter/ureg.second),
            self.c_cl: ureg.newton/(ureg.meter/ureg.second),
            self.x_l: ureg.meter,
            self.x_r: ureg.meter,
            self.x_l.diff(self.ivar): ureg.meter/ureg.second,
            self.x_r.diff(self.ivar): ureg.meter/ureg.second,
            self.x_l.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.x_r.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.R: ureg.meter,
        }
        return unit_dict

    
#DONE(#205)
class DampedTrolleysWithSprings(DoubleTrolleyDifferentWheels):
#    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.png'
    scheme_name = 'MDOFDoubleTrolley.png'
    real_name = 'two_damped_trolleys.png'

    R=Symbol('R', positive=True)
    m=Symbol('m', positive=True)
    m1=Symbol('m_1', positive=True)
    m2=Symbol('m_2', positive=True)
    k_l=Symbol('k_l', positive=True)
    c_cl=Symbol('c_cl', positive=True)
    k_c=Symbol('k_c', positive=True)
    c_cc=Symbol('c_cc', positive=True)
    k_r=Symbol('k_r', positive=True)
    c_cr=Symbol('c_cr', positive=True)
    lam=Symbol('lambda', positive=True)
    x_l=dynamicsymbols('x_l')
    x_r=dynamicsymbols('x_r')

    qs=[x_l, x_r]

    def __init__(self,
                 R=None,
                 m1=None,
                 m2=None,
                 m=None,
                 k_l=None,
                 k_r=None,
                 k_c=None,
                 c_cl=None,
                 c_cr=None,
                 c_cc=None,
                 x_l=None,
                 x_r=None,
                 lam=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if x_l is not None: self.x_l = x_l
        if x_r is not None: self.x_r = x_r
        if m1 is not None: self.m1 = m1
        if m2 is not None: self.m2 = m2
        if m is not None: self.m = m
        if k_l is not None: self.k_l = k_l
        if c_cl is not None: self.c_cl = c_cl
        if k_c is not None: self.k_c = k_c
        if c_cc is not None: self.c_cc = c_cc
        if k_r is not None: self.k_r = k_r
        if c_cr is not None: self.c_cr = c_cr
        if R is not None: self.R = R
        if lam is not None: self.lam = lam
        if qs is not None: self.qs=qs

        self.ivar=ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._system = DoubleTrolleyDifferentWheels(R=self.R,
                 m_w1=self.m,
                 m_w2=self.m,
                 m1=self.m1,
                 m2=self.m2,
                 k_l=self.k_l,
                 c_cl=self.c_cl,
                 k_c=self.k_c,
                 c_cc=self.c_cc,
                 k_r=self.k_r,
                 c_cr=self.c_cr,
                 lam=self.lam,
                 x_l=self.x_l,
                 x_r=self.x_r,
                 qs=self.qs,
                 ivar=self.ivar
                 )

        components['_system'] = self._system
        return components

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

        }

        return default_data_dict


    def get_numerical_data(self):

        numerical_data_dict = {
            self.m1: [1],
            self.m2: [1],
            self.m: [0.1],
            self.k_l: [50],
            self.k_c: [70],
            self.k_r: [100],

            self.c_cr: [3.5],
            self.c_cc: [4.0],
            self.c_cl: [3.0],
            self.R: [0.05],
        }

        return numerical_data_dict

    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley Mass',
            self.m2: r'Trolley Mass',
            self.m: r'Wheels Mass',
            self.k_l: r'Spring Stiffness',
            self.k_c: r'Spring Stiffness',
            self.k_r: r'Spring Stiffness',

            self.c_cr: r'Damper damping',
            self.c_cc: r'Damper damping',
            self.c_cl: r'Damper damping',
            self.R: r'Wheels radius',
        }
        return self.sym_desc_dict
    
    def unit_dict(self):

        unit_dict = {

            self.ivar: ureg.second,

            
            self.m1: ureg.kilogram,
            self.m2: ureg.kilogram,
            self.m: ureg.kilogram,
            self.k_l: ureg.newton/ureg.meter,
            self.k_c: ureg.newton/ureg.meter,
            self.k_r: ureg.newton/ureg.meter,

            self.c_cr: ureg.newton/(ureg.meter/ureg.second),
            self.c_cc: ureg.newton/(ureg.meter/ureg.second),
            self.c_cl: ureg.newton/(ureg.meter/ureg.second),
            self.x_l: ureg.meter,
            self.x_r: ureg.meter,
            self.x_l.diff(self.ivar): ureg.meter/ureg.second,
            self.x_r.diff(self.ivar): ureg.meter/ureg.second,
            self.x_l.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.x_r.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.R: ureg.meter,
        }
        return unit_dict


# ISSUE 211 DONE karolina i bogus
# SAME AS ISUUE 207
class ForcedTrolleysWithSprings(ComposedSystem):

    scheme_name = 'mdof_three_trolleys.PNG'
    real_name = 'three_carriages.PNG'

    m1=Symbol('m_1', positive=True)
    m2=Symbol('m_2', positive=True)
    m3=Symbol('m_3', positive=True)
    k_l=Symbol('k_l', positive=True)
    k_cl=Symbol('k_cl', positive=True)
    k_12=Symbol('k_12', positive=True)
    k_c12=Symbol('k_c12', positive=True)
    k_23=Symbol('k_23', positive=True)
    k_c23=Symbol('k_c23', positive=True)
    k_r=Symbol('k_r', positive=True)
    k_cr=Symbol('k_cr', positive=True)
    F=Symbol('F', positive=True)
    F_1=Symbol('F_1', positive=True)
    F_2=Symbol('F_2', positive=True)
    Omega=Symbol('Omega', positive=True)

    x_l=dynamicsymbols('x_l')
    x_c=dynamicsymbols('x_c')
    x_r=dynamicsymbols('x_r')

    qs=dynamicsymbols('x_l x_c x_r')
    ivar=Symbol('t')

    def __init__(self,
                 m1=None,
                 m2=None,
                 m3=None,
                 k_l=None,
                 k_cl=None,
                 k_12=None,
                 k_c12=None,
                 k_23=None,
                 k_c23=None,
                 k_r=None,
                 k_cr=None,
                 x_l=None,
                 x_c=None,
                 x_r=None,
                 F=None,
                 F_1=None,
                 F_2=None,
                 Omega=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m1 is not None: self.m1 = m1
        if m2 is not None: self.m2 = m2
        if m3 is not None: self.m3 = m3
        if k_l is not None: self.k_l = k_l
        if k_cl is not None: self.k_cl = k_cl
        if k_12 is not None: self.k_12 = k_12
        if k_c12 is not None: self.k_c12 = k_c12
        if k_23 is not None: self.k_23 = k_23
        if k_c23 is not None: self.k_c23 = k_c23
        if k_r is not None: self.k_r = k_r
        if k_cr is not None: self.k_cr = k_cr
        if x_l is not None: self.x_l = x_l
        if x_c is not None: self.x_c = x_c
        if x_r is not None: self.x_r = x_r
        if F is not None: self.F = F
        if F_1 is not None: self.F_1 = F_1
        if F_2 is not None: self.F_2 = F_2
        if Omega is not None: self.Omega = Omega


        if qs is not None: self.qs =qs

        self.ivar=ivar

        self._init_from_components(**kwargs)



    @cached_property
    def components(self):

            components = {}
            F_1=-self.F * cos(self.Omega * self.ivar)
            F_2=-2 * self.F * cos(self.Omega * self.ivar)

            self._trolley_1 = ForcedTrolleyWithThreeSprings(self.m1, self.k_l, self.k_cl, F_1, z=self.x_l, qs=[self.x_l] )

            self._trolley_2 = MaterialPoint(self.m2, pos1=self.x_c, qs=self.qs)
            self._springs1=ThreeSprings(self.k_12, self.k_c12, pos1=self.x_l, pos2=self.x_c, qs=self.qs)
            self._springs2=ThreeSprings(self.k_23, self.k_c23, pos1=self.x_c, pos2=self.x_r, qs=self.qs)

            self._trolley_3 = ForcedTrolleyWithThreeSprings(self.m3, self.k_r, self.k_cr, F_2, z=self.x_r, qs=[self.x_r])


            components['_trolley_1'] = self._trolley_1
        
            components['_trolley_2'] = self._trolley_2
            components['_springs1'] = self._springs1
            components['_springs2'] = self._springs2
        
            components['_trolley_3'] = self._trolley_3
        
            return components
    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            
            self.k_12: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
        }

        return default_data_dict


    def get_numerical_data(self):
        
        numerical_data_dict = {
            self.m1: [1],
            self.m2: [2],
            self.m3: [3],
            
            self.k_l: [5],
            self.k_cl: [8],
            
            self.k_23: [5],
            self.k_c23: [8],
            
            self.k_12: [5],
            self.k_c12: [8],
            
            self.k_r: [5],
            self.k_cr: [8],
            self.F: [5],
            
            self.Omega: [60],

        }

        return numerical_data_dict

    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley 1 Mass',
            self.m2: r'Trolley 2 Mass',
            self.m3: r'Trolley 3 Mass',

            self.k_l: r'Spring Stiffness',
            self.k_cl: r'Spring Stiffness',
            
            self.k_23: r'Spring Stiffness',
            self.k_c23: r'Spring Stiffness',
            
            self.k_12: r'Spring Stiffness',
            self.k_c12: r'Spring Stiffness',
            
            self.k_r: r'Spring Stiffness',
            self.k_cr: r'Spring Stiffness',


        }
        return self.sym_desc_dict
    
    def unit_dict(self):

        unit_dict = {

            self.ivar: ureg.second,

            
            self.m1: ureg.kilogram,
            self.m2: ureg.kilogram,
            self.m3: ureg.kilogram,

            self.k_l: ureg.newton/ureg.meter,
            self.k_cl: ureg.newton/ureg.meter,
            
            self.k_23: ureg.newton/ureg.meter,
            self.k_c23: ureg.newton/ureg.meter,
            
            self.k_12: ureg.newton/ureg.meter,
            self.k_c12: ureg.newton/ureg.meter,
            
            self.k_r: ureg.newton/ureg.meter,
            self.k_cr: ureg.newton/ureg.meter,
            

            self.x_l.diff(self.ivar): ureg.meter/ureg.second,
            self.x_r.diff(self.ivar): ureg.meter/ureg.second,
            self.x_c.diff(self.ivar): ureg.meter/ureg.second,
            self.x_l.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.x_r.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,
            self.x_c.diff(self.ivar,2): ureg.meter/ureg.second/ureg.second,

            self.F: ureg.newton,
            self.F_1: ureg.newton,
            self.F_2: ureg.newton,
            self.Omega: ureg.radian/ureg.second,
        }
        return unit_dict
    

#poprawione MDofTMD- ready to check
#klasa MDoFTMD została usunięta- zastępuje ją teraz TrolleyWithTMD
#TODO(#212 - check)




