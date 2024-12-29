from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset, atan,Heaviside,sign)

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
                 ivar=None,
                 **kwargs):

        
        
        if m is not None: self.m = m
        if k is not None: self.k = k
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        
   
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @cached_property
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

        Exampled
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

    @cached_property
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
#owner: KrzychT
#supervision: Boogi
class EquivalentGearModel(ComposedSystem):
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

    m=Symbol('m_eq', positive=True)
    k=Symbol('k_g', positive=True)
    F=Symbol('F', positive=True)
    c=Symbol('c_g',positive=True)
    T=Symbol('T',positive=True)
    omega=Symbol('omega',positive=True)
    ivar=Symbol('t')
    k_var=Symbol('kappa_mesh', positive=True)
    eps=Symbol('varepsilon', positive=True)
    c_var=Symbol('c_var', positive=True)
    f=Symbol('f', positive=True)
    
    z=dynamicsymbols('z')
    
    def __init__(self,
                 m=None,
                 k=None,
                 F=None,
                 z=None,
                 c=None,
                 T=None,
                 ivar=None,
                 k_var=None,
                 eps=None,
                 c_var=None,
                 f=None,
                 **kwargs):

        if k_var is not None: self.k_var = k_var
        if eps is not None: self.eps = eps
        if m is not None: self.m = m
        if k is not None: self.k = k
        if F is not None: self.F = F
        if z is not None: self.z = z
        if c is not None: self.c= c
        if T is not None: self.T= T
        if c_var is not None: self.c_var = c_var
        if ivar is not None: self.ivar = ivar
        if f is not None: self.f=f

        
   
        self.qs = [self.z]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}
        self.c_var = self.c*(1+0*self.eps*self.k_var )
        stiffness = self.k*(1+self.eps*self.k_var  ) 
        
        
        self.gear_inertia = MaterialPoint(self.m, self.z, qs=self.qs)
        self.gear_stiffness = Spring(stiffness, self.z, qs=self.qs)
        self.force = Force(self.F, self.z, qs=self.qs)
        self.gear_damping = Damper(self.c_var, self.z, qs=self.qs)
        
        components['gear_inertia'] = self.gear_inertia
        components['gear_stiffness'] = self.gear_stiffness
        components['force'] = self.force
        components['gear_damping']=self.gear_damping
        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of the system on spring',
            self.F: r'excitation force',
            self.c: r'damping constant',
            self.k: r'stiffness',
            self.T: r'torque [not sure]',
        }

        return self.sym_desc_dict
    def _report_components(self):

        comp_list = [
            *REPORT_COMPONENTS_LIST
        ]

        return comp_list

    def units(self):
        f=Symbol('f')
        units_dict={
            self.k:ureg.newton/ureg.meter,
            self.m:ureg.kilogram,
            self.F:ureg.newton,
            self.ivar:ureg.second,
            self.T:ureg.newton*ureg.meter,
            self.z.diff(self.ivar,2):ureg.meter/ureg.second**2,
            self.c:ureg.newton*ureg.second/ureg.meter,
            self.f:ureg.hertz
           }
        return units_dict
    
    def get_numerical_data(self):

        k,m,c,F,eps = symbols('k m c F epsilon', positive=True)

        default_data_dict = {
            self.k : [1e6],
            self.m : [1],
            self.c : [2e-4*1e6],
            self.F : [100],
            self.eps : [0.1],
        }
        return default_data_dict
    
    def trig_stiff(self, angle=2*pi):
        trig = sin(self.omega * self.ivar)
        new_eq = self.subs(self.k_var,trig)
        
        return new_eq
    
    def wave_stiff(self):
        wave1=(100*(sin(self.ivar/self.T))+1000)*Heaviside(sin(self.ivar/self.T))
        wave2=(100*(-sin(self.ivar/self.T)))*Heaviside(-sin(self.ivar/self.T))
        waves=wave1+wave2
        new_eq = self.subs(self.k_var,waves)

        return new_eq
        
    def rect_stiff(self,no=6,numerical=False):
        t=self.ivar
        omg = self.omega
        
        trig=sum([Heaviside(omg*t-1) + 0* Heaviside(omg*t-2)  for ind in range(no)])
        new_eq=self.subs(self.k_var,trig)
        
        return new_eq
    
    def approx_rect(self,no=6,numerical=False):
        
        if numerical is True:
            amps_list = [2.27348466531425,0, 0.757408805199249,0, 0.453942816897038,0, 0.323708002807428,0, 0.25121830779797,0, 0.204977919963796,0, 0.172873394602606,0, 0.149252079729775,0, 0.131121653619234,0, 0.116749954968057,0]
        else:
            amps_list = symbols(f'a_0:{no}')
        
        rectangular_approx = sum([N(amp,3)*sin(((ind)+1)*self.omega*self.ivar) for ind,amp in enumerate(amps_list[0:no])])
        new_eq=self.subs({self.k_var:(rectangular_approx)})
        
        return new_eq
    
    def ode_with_delta(self):
        delta=Symbol('delta', positive=True)
        eps =Symbol('varepsilon', positive=True)
        
        with_delta = self._eoms[0]+self.k*eps*delta*self.z
        delta_sys = type(self._ode_system)(with_delta, Matrix([self.z]), ivar=self.ivar, ode_order=2)

        return delta_sys


    
    
    def _stiffness_models(self):
    
        #parabola
        #para_series = 80*(t_values%(T_value/2))*((t_values%(T_value/2))-(T_value/2))*(-0.4)
        #para_series_heave =  10*np.heaviside(np.sin(2*np.pi/(T_value)*t_values),0.5)
        #para_values=(0.5*(para_series+para_series_heave-6))
        

        #wave
        t = self.ivar
        T = self.T
        
        wave1=(1*(sin(2*pi*t/T))+2.0)*(1/2 + 1/2*sign(sin(2*pi*t/T)))
        wave2=(1*(-sin(2*pi*t/T))-3.0)*(1/2+ 1/2*sign(-sin(2*pi*t/T)))
        waves=wave1+wave2



        #rectangular
        rectangular=5*((1/2+1/2*sign(sin(2*pi*t/T)))-S.Half)



        #rectangular_approx
        amps_list = [2.27348466531425, 0.757408805199249, 0.453942816897038, 0.323708002807428, 0.25121830779797, 0.204977919963796, 0.172873394602606, 0.149252079729775, 0.131121653619234, 0.116749954968057]
        rectangular_approx = sum([amp*2/sqrt(2)*sin((2*(no)+1)*2*pi*t/T) for no,amp in enumerate(amps_list[0:])])


        
        return {'wave':waves, 'rect':rectangular, 'approx':rectangular_approx}

    def _stiffness_waveforms(self):
        
        from sympy import lambdify
        
        t = self.ivar
        T = self.T
        
        return {label:lambdify((t, T), waveform)   for label,waveform in self._stiffness_models().items() }
    
    
    
class DDOFGearMechanism(ComposedSystem):
    scheme_name = 'MDOF_Forced_Disks_With_Serial_Springs.PNG'
    real_name = 'three_carriages.PNG'

    r1 = Symbol('r_1', positive=True)
    r2 = Symbol('r_2', positive=True)
    J1 = Symbol('J_1', positive=True)
    J2 = Symbol('J_2', positive=True)
    k = Symbol('k', positive=True)
    c = Symbol('c', positive=True)
    T1 = Symbol('T_1', positive=True)
    T2 = Symbol('T_2', positive=True)
    Omega = Symbol('Omega', positive=True)
    phi1 = dynamicsymbols('phi_1')
    phi2 = dynamicsymbols('phi_2')
    qs = dynamicsymbols('phi_1 phi_2')
    ivar = Symbol('t')

    def __init__(self,
                 r1=None,
                 r2=None,
                 J1=None,
                 J2=None,
                 k=None,
                 c=None,
                 T1=None,
                 T2=None,
                 Omega=None,
                 phi1=None,
                 phi2=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if r1 is not None: self.r1 = r1
        if r2 is not None: self.r2 = r2
        if J1 is not None: self.J1 = J1
        if J2 is not None: self.J2 = J2
        if k is not None: self.k = k
        if c is not None: self.c = c
        if T1 is not None: self.T1 = T1
        if T2 is not None: self.T2 = T2
        if Omega is not None: self.Omega = Omega
        if phi1 is not None: self.phi1 = phi1
        if phi2 is not None: self.phi2 = phi2
        if qs is not None: self.qs = qs

        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        self.gear1 = Disk(self.J1, pos1=self.phi1, qs=self.qs, ivar = self.ivar)(label='Drive gear')
        self.gear2 = Disk(self.J2, pos1=self.phi2, qs=self.qs, ivar = self.ivar)(label='Driven gear')


        self.gears_equivalent_stiffness = Spring(self.k, pos1=self.phi1*self.r1, pos2=self.phi2*self.r2, qs=self.qs)(label='Gears equivalent stiffness of teeth')
        self.gears_equivalent_damping = Damper(self.c, pos1=self.phi1*self.r1, pos2=self.phi2*self.r2, qs=self.qs)(label='Gears equivalent damping coefficient of teeth')


        self.torque1 = Force(self.T1 * cos(self.Omega * self.ivar), pos1 = self.phi1, qs = self.qs)(label='Torque of the driving gear')
        self.torque2 = Force(self.T2 * cos(self.Omega * self.ivar), pos1 = self.phi2, qs = self.qs)(label='Torque of the driven gear')


        components['_gear1'] = self.gear1
        components['_gear2'] = self.gear2
        components['_gears_equivalent_stiffness'] = self.gears_equivalent_stiffness
        components['_gears_equivalent_damping'] = self.gears_equivalent_damping
        components['_torque1'] = self.torque1
        components['_torque2'] = self.torque2

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.J1: r'moment of inertia of the driving gear',
            self.J2: r'moment of inertia of the driven gear',
            self.r1: r'rolling radius of the driving gear',
            self.r2: r'rolling radius of the driven gear',
            self.T1: r'torque of the driving gear',
            self.T2: r'torque of the driven gear',
            self.k: r'meshing stiffness',
            self.c: r'meshing damping',
            self.Omega: r'rotation frequency',
        }

        return self.sym_desc_dict

    def units(self):
        units_dict={
            self.J1:ureg.kilogram*ureg.meter**2,
            self.J2:ureg.kilogram*ureg.meter**2,
            self.Omega:ureg.radian/ureg.second,
            self.r1:ureg.meter,
            self.r2:ureg.meter,
            self.ivar:ureg.second,
            self.T1:ureg.newton*ureg.meter,
            self.T2:ureg.newton*ureg.meter,
           }
        return units_dict
    
    def as_equivalent_sdof(self):
        return EquivalentSDOFGearModel()

class EquivalentSDOFGearModel(EquivalentGearModel):
    pass
<<<<<<< HEAD
=======

>>>>>>> 425a651 ( On branch master)
