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

from pint import UnitRegistry
ureg = UnitRegistry()

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST

    
class SDoFShaft(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents the Kinematicly excited shaft with two disks.
    =========
            I = Moment of Inertiafrom sympy import *
from sympy.plotting import plot
from pylatex import Document, Section, Subsection, Subsubsection, Command, Itemize, Package, HorizontalSpace, Description, Marker
from pylatex.section import Paragraph, Chapter
from pylatex.utils import italic, NoEscape

from dynpy.utilities.report import * #to do tworzenia raportu pozniej

from dynpy.utilities.adaptable import *
import pypandoc as ppd
import importlib


from pylatex import Package, Document, NoEscape, Command, Section
from sympy.physics.mechanics import *
from dynpy.solvers.linear import * #ODESystem here
from dynpy.models.mechanics.pendulum import *

from dynpy.models.mechanics.pendulum import *  #skad zaciagasz modele
 
mechanics_printing() #mechaniczny zapis pochodnej to kropki, ta linijka zmienia z d/dt na kropkis
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
        self.qs = [self.phi]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.k_1 = (self.G * self.I_1) / self.l_1
        self.k_2 = (self.G * self.I_2) / self.l_2

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_2 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1),
                               pos1=self.phi,
                               pos2=self.theta,
                               qs=self.qs)  # right spring
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)


        components['disc_1'] = self.disc_1
        components['spring_2'] = self.spring_2
        components['moment'] = self.moment

        return components


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
    
    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticKeyLength,
        mech_comp.DynamicKeyLength,
        mech_comp.DynamicTorqueComponent

        ]

        return comp_list

    def disc_force(self):
        t = self.ivar
        return self.I * self.steady_solution()[0].diff(t, t)

    def max_static_force(self):
        d = Symbol('d', positive=True)
        return 2 * self.Ms / d

    def max_dynamic_force(self):
        d, theta0, Omega = symbols('d theta_0, Omega', positive=True)
#         return self.frequency_response_function(
#             self.natural_frequencies()[0]) * self.stiffness_matrix()[0]
        return self.subs({self.theta:theta0*cos(Omega*self.ivar)}).frequency_response_function() * self.spring_2.stiffness *2/d + self.max_static_force()

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
        return (2 * self.max_static_force()) / (kd * h)

    def dynamic_key_length(self):
        kd = Symbol('k_d', positive=True)
        h = Symbol('h', positive=True)
        return (2 * self.max_dynamic_force()) / (kd * h)
    
    def dynamic_torque(self):
        t=self.ivar
        steady=self._fodes_system.steady_solution[0]
        torque=steady.diff(t,t)*self.I+self.Ms
        return torque
        
class SDoFShaftHarmonicExcitation(SDoFShaft):
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
    theta0=Symbol('theta_0', positive=True)
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
                 theta0=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 qs=None,
                 **kwargs):
        if G is not None: self.G = G

        if I is not None: self.I = I
        if Ms is not None: self.Ms = Ms
        if Omega is not None: self.Omega = Omega
        if l_1 is not None: self.l_1 = l_1
        if l_2 is not None: self.l_2 = l_2
        if I_1 is not None: self.I_1 = I_1
        if I_2 is not None: self.I_2 = I_2
        if phi is not None: self.phi = phi
        if theta0 is not None: self.theta0 = theta0
        self.qs = [self.phi]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.k_1 = (self.G * self.I_1) / self.l_1
        self.k_2 = (self.G * self.I_2) / self.l_2

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_2 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1),
                               pos1=self.phi,
                               pos2=self.theta0*cos(self.Omega*self.ivar),
                               qs=self.qs)  # right spring
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)


        components['disc_1'] = self.disc_1
        components['spring_2'] = self.spring_2
        components['moment'] = self.moment

        return components
    
    def get_default_data(self):

        m0, l0, G, l, theta0, Omega0 = symbols('m_0 l_0 G l theta_0 Omega_0', positive=True)

        default_data_dict = {
            self.I: [S.Half * m0 * (l0**2) * no for no in range(1, 3)],
            self.I_1: [S.Half**(no) * (l0**4) for no in range(1, 8)],
            self.I_2: [S.Half**no * (l0**4) for no in range(1, 8)],
            self.l_1: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.l_2: [S.Half**(no - 6) * l0 for no in range(1, 8)],
            self.theta0: [S.One * self.theta0 * no for no in range(1,5)],
            self.Omega: [S.One * Omega0],
        }

        return default_data_dict

    def max_dynamic_force(self):
        d, theta0= symbols('d theta_0', positive=True)
#         return self.frequency_response_function(
#             self.natural_frequencies()[0]) * self.stiffness_matrix()[0]
        return self.subs({self.theta:theta0}).frequency_response_function() * self.spring_2.stiffness *2/d + self.max_static_force()

class DoubleDiskShaft(ComposedSystem):
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

    scheme_name = 'ddof_undamped_shaft.png'
    real_name = 'ddof_shaft_real.png'
    
    d=Symbol('d', positive=True)
    G=Symbol('G', positive=True)
    
    l=Symbol('l', positive=True)
    l_1=Symbol('l_1', positive=True)
    l_2=Symbol('l_2', positive=True)
    
    I_1=Symbol('I_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    I_m1=Symbol('I_m1', positive=True)
    I_m2=Symbol('I_m2', positive=True)
    
    k_2=Symbol('k_2', positive=True)
    k_1=Symbol('k_1', positive=True)
    
    T_1=Symbol('T_1',positive=True)
    T_2=Symbol('T_2',positive=True)
    
    Omega=Symbol('Omega',positive=True)
 
    theta=dynamicsymbols('theta')
    phi_1=dynamicsymbols('\\varphi_1')                 
    phi_2=dynamicsymbols('\\varphi_2')                 
    phi=dynamicsymbols('\\varphi')

    m0 = Symbol('m_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    T0 = Symbol('T_0', positive=True)

    theta0 = Symbol('theta_0', positive=True)
    Omega =Symbol('Omega', positive=True)

    def __init__(self,
                 l=None,
                 I_m1=None,
                 I_m2=None,
                 k_2=None,
                 k_1=None,
                 T_1=None,
                 T_2=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 phi_1=None,                 
                 phi_2=None,                 
                 phi=None,
                 d=None,
                 theta=None,
                 ivar=Symbol('t'),
                 qs=None,
                 
                 **kwargs):
        
        if  l_1 is not None: self.l_1 = l_1
        if  l_2 is not None: self.l_2 = l_2
        if  I_1 is not None: self.I_1 = I_1
        if  I_2 is not None: self.I_2 = I_2
        if  I_m1 is not None: self.I_m1 = I_m1
        if  I_m2 is not None: self.I_m2 = I_m2
        if T_1 is not None: self.T_1 = T_1
        if T_2 is not None: self.T_2 = T_2
        #if Omega is not None: self.Omega = Omega
        if  k_1 is not None: self.k_1 = k_1
        if  k_2 is not None: self.k_2 = k_2
        if  phi_1 is not None: self.phi_1 = phi_1
        if  phi_2 is not None: self.phi_2 = phi_2 
        if  phi is not None: self.phi = phi 
        if  theta is not None: self.theta = theta 
        #if d is not None: self.d=d    
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.qs = [self.phi_1,self.phi_2]
        self.k_1 = (self.G*self.I_1)/self.l_1
        self.k_2 = (self.G*self.I_2)/self.l_2

        self.disc_1 = Disk(self.I_m1, pos1=self.phi_1, qs=self.qs)
        self.spring_1 = Spring(self.k_1, self.phi_1, self.phi_2, qs=self.qs)  # left spring
        
        self.disc_2 = Disk(self.I_m2, pos1=self.phi_2, qs=self.qs)
        self.spring_2 = Spring(self.k_2, pos1=self.phi_2, pos2=self.theta,
                               qs=self.qs)  # right spring
        self.moment_disc1=Force(self.T_1, pos1=self.phi_1)
        self.moment_disc2=Force(self.T_2, pos1=self.phi_2)

        components['disc_1'] = self.disc_1
        components['disc_2'] = self.disc_2
        components['spring_1'] = self.spring_1
        components['spring_2'] = self.spring_2
        components['moment_disc1'] = self.moment_disc1
        components['moment_disc2'] = self.moment_disc2

        return components


    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 , G, l,  d, T0 = self.m0, self.l0 , self.G, self.l,  self.d, self.T0
        theta0, Omega = self.theta0, self.Omega

        default_data_dict = {
            self.I_m1: [S.Half*m0*(l0**2)*no for no in range(1,9)],
            self.I_m2: [S.Half*m0*(l0**2)*no for no in range(1,9)],
            self.I_1: [S.Half**(no)*(l0**4) for no in range(1,9)],
            self.I_2: [S.Half**no*(l0**4) for no in range(1,9)],
            self.l_1: [S.Half**(no-6)*l0 for no in range(1,9)],
            self.l_2: [S.Half**(no-6)*l0 for no in range(1,9)],
            self.T_1: [T0 * (no) for no in range(1,12)],
            self.T_2: [T0 * (no) for no in range(1,12)],
            self.theta:[theta0* cos(self.Omega * self.ivar) ],
            
        }

        return default_data_dict


    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST
        ]

        return comp_list
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def dynamic_force(self):
        data=self._given_data
        amps = self._fodes_system.steady_solution.as_dict()
        dyn_force=(self.components['spring_2'].force().subs(amps)).subs(data).expand().doit()
        
        return dyn_force
    
    def disc_1_force(self):
        t=self.ivar
        return self.I_m1 * self.steady_solution().diff(t,t)

    def max_static_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[0]/d)
    
    def max_dynamic_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[0]*self.Omega**2

        return  abs(2*(self.I_m1*acc_amp)/d) + self.max_static_disk_1_bearing_force()#.subs(self._given_data)

    def max_static_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[1]/d)
    
    def max_dynamic_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[1]*self.Omega**2

        return  abs(2*(self.I_m2*acc_amp)/d) + self.max_static_disk_2_bearing_force()#.subs(self._given_data)
    
    
    def static_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_1_bearing_force())/(kd*h)
    
    def dynamic_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_1_bearing_force())/(kd*h)
    
    
    def static_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_2_bearing_force())/(kd*h)
    
    def dynamic_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_2_bearing_force())/(kd*h)
    

    
class DoubleDiskShaftHarmonicExcitation(DoubleDiskShaft):
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

    scheme_name = 'ddof_undamped_shaft.png'
    real_name = 'ddof_shaft_real.png'
    
    d=Symbol('d', positive=True)
    G=Symbol('G', positive=True)
    l=Symbol('l', positive=True)
    l_1=Symbol('l_1', positive=True)
    l_2=Symbol('l_2', positive=True)
    I_1=Symbol('I_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    I_m1=Symbol('I_m1', positive=True)
    I_m2=Symbol('I_m2', positive=True)
    k_2=Symbol('k_2', positive=True)
    k_1=Symbol('k_1', positive=True)
    T_1=Symbol('T_1',positive=True)
    T_2=Symbol('T_2',positive=True)
    Omega=Symbol('Omega',positive=True)
#     theta=dynamicsymbols('theta')
    phi_1=dynamicsymbols('\\varphi_1')                 
    phi_2=dynamicsymbols('\\varphi_2')                 
    phi=dynamicsymbols('\\varphi')
    m0 = Symbol('m_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    T0 = Symbol('T_0', positive=True)
    theta = Symbol('theta', positive=True)
    Omega =Symbol('Omega', positive=True)

    def __init__(self,
                 l=None,
                 I_m1=None,
                 I_m2=None,
                 k_2=None,
                 k_1=None,
                 T_1=None,
                 T_2=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 phi_1=None,                 
                 phi_2=None,                 
                 phi=None,
                 d=None,
                 theta=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 qs=None,
                 
                 **kwargs):
        
        if  l_1 is not None: self.l_1 = l_1
        if  l_2 is not None: self.l_2 = l_2
        if  I_1 is not None: self.I_1 = I_1
        if  I_2 is not None: self.I_2 = I_2
        if  I_m1 is not None: self.I_m1 = I_m1
        if  I_m2 is not None: self.I_m2 = I_m2
        if T_1 is not None: self.T_1 = T_1
        if T_2 is not None: self.T_2 = T_2
        if Omega is not None: self.Omega = Omega
        if  k_1 is not None: self.k_1 = k_1
        if  k_2 is not None: self.k_2 = k_2
        if  phi_1 is not None: self.phi_1 = phi_1
        if  phi_2 is not None: self.phi_2 = phi_2 
        if  phi is not None: self.phi = phi 
        if  theta is not None: self.theta = theta
        #if d is not None: self.d=d
        self._init_from_components(**kwargs)
            
    @property
    def components(self):
        components = {}

        self.qs = [self.phi_1,self.phi_2]
        self.k_1 = (self.G*self.I_1)/self.l_1
        self.k_2 = (self.G*self.I_2)/self.l_2

        self.disc_1 = Disk(self.I_m1, pos1=self.phi_1, qs=self.qs)
        self.spring_1 = Spring(self.k_1, self.phi_1, self.phi_2, qs=self.qs)  # left spring
        self.disc_2 = Disk(self.I_m2, pos1=self.phi_2, qs=self.qs)
        self.spring_2 = Spring(self.k_2, pos1=self.phi_2, pos2=self.theta*cos(self.Omega*self.ivar),qs=self.qs)  # right spring
        self.moment_disc1=Force(self.T_1, pos1=self.phi_1)
        self.moment_disc2=Force(self.T_2, pos1=self.phi_2)

        components['disc_1'] = self.disc_1
        components['disc_2'] = self.disc_2
        components['spring_1'] = self.spring_1
        components['spring_2'] = self.spring_2
        components['moment_disc1'] = self.moment_disc1
        components['moment_disc2'] = self.moment_disc2

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

#     def get_default_data(self):

#         m0, l0 , G, l,  d, T0 = self.m0, self.l0 , self.G, self.l,  self.d, self.T0
#         theta0, Omega0 = symbols('theta_0 Omega_0', positive=True)

#         default_data_dict = {
#             self.I_m1: [S.Half*m0*(l0**2)*no for no in range(1,9)],
#             self.I_m2: [S.Half*m0*(l0**2)*no for no in range(1,9)],
#             self.I_1: [S.Half**(no)*(l0**4) for no in range(1,9)],
#             self.I_2: [S.Half**no*(l0**4) for no in range(1,9)],
#             self.l_1: [S.Half**(no-6)*l0 for no in range(1,9)],
#             self.l_2: [S.Half**(no-6)*l0 for no in range(1,9)],
#             self.T_1: [T0 * (no) for no in range(1,12)],
#             self.T_2: [T0 * (no) for no in range(1,12)],
#             self.theta:[S.One * theta0 * (no) for no in range(1,5)],
# #             self.Omega:[S.One * Omega0],
            
#         }

#         return default_data_dict

    def get_default_data(self):

        m0, l0 , G, l,  d, T0 = self.m0, self.l0 , self.G, self.l,  self.d, self.T0
        theta0, Omega = self.theta0, self.Omega

        default_data_dict = {
            self.k_1:[self.k_1],
            self.k_2:[self.k_2],
            self.I_m1: [S.Half*m0*(l0**2)],
            self.I_m2: [S.Half*m0*(l0**2)],
            self.I_1: [S.Half*(l0**4)],
            self.I_2: [S.Half*(l0**4)],
            self.l_1: [S.Half*l0,2*S.Half*l0,4*S.Half*l0],
            self.l_2: [S.Half*l0,2*S.Half*l0,4*S.Half*l0],
            self.T_1: [T0],
            self.T_2: [T0],
        }

        return default_data_dict

    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST
        ]

        return comp_list
    
    def disc_1_force(self):
        t=self.ivar
        return self.I_m1 * self.steady_solution().diff(t,t)

    def max_static_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[0]/d)
    
    def max_dynamic_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[0]*self.Omega**2

        return  abs(2*(self.I_m1*acc_amp)/d) + self.max_static_disk_1_bearing_force()#.subs(self._given_data)

    def max_static_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[1]/d)
    
    def max_dynamic_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[1]*self.Omega**2

        return  abs(2*(self.I_m2*acc_amp)/d) + self.max_static_disk_2_bearing_force()#.subs(self._given_data)
    
    
    def static_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_1_bearing_force())/(kd*h)
    
    def dynamic_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_1_bearing_force())/(kd*h)
    
    
    def static_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_2_bearing_force())/(kd*h)
    
    def dynamic_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_2_bearing_force())/(kd*h)
    
class DampedDoubleDiskShaftHarmonicExcitation(DoubleDiskShaftHarmonicExcitation):
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

    scheme_name = 'ddof_undamped_shaft.png'
    real_name = 'ddof_shaft_real.png'
    
    d=Symbol('d', positive=True)
    G=Symbol('G', positive=True)
    l=Symbol('l', positive=True)
    l_1=Symbol('l_1', positive=True)
    l_2=Symbol('l_2', positive=True)
    I_1=Symbol('I_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    I_m1=Symbol('I_m1', positive=True)
    I_m2=Symbol('I_m2', positive=True)
    k_2=Symbol('k_2', positive=True)
    k_1=Symbol('k_1', positive=True)
    T_1=Symbol('T_1',positive=True)
    T_2=Symbol('T_2',positive=True)
    Omega=Symbol('Omega',positive=True)
#     theta=dynamicsymbols('theta')
    phi_1=dynamicsymbols('\\varphi_1')                 
    phi_2=dynamicsymbols('\\varphi_2')                 
    phi=dynamicsymbols('\\varphi')
    m0 = Symbol('m_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    G = Symbol('G', positive=True)
    T0 = Symbol('T_0', positive=True)
    theta = Symbol('theta', positive=True)
    Omega =Symbol('Omega', positive=True)
    c1=Symbol('c_1',positive=True)
    c2=Symbol('c_2',positive=True)

    def __init__(self,
                 l=None,
                 I_m1=None,
                 I_m2=None,
                 k_2=None,
                 k_1=None,
                 T_1=None,
                 T_2=None,
                 l_1=None,
                 l_2=None,
                 I_1=None,
                 I_2=None,
                 phi_1=None,                 
                 phi_2=None,                 
                 phi=None,
                 d=None,
                 theta=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 qs=None,
                 c1=None,
                 c2=None,
                 
                 **kwargs):
        
        if  l_1 is not None: self.l_1 = l_1
        if  l_2 is not None: self.l_2 = l_2
        if  I_1 is not None: self.I_1 = I_1
        if  I_2 is not None: self.I_2 = I_2
        if  I_m1 is not None: self.I_m1 = I_m1
        if  I_m2 is not None: self.I_m2 = I_m2
        if T_1 is not None: self.T_1 = T_1
        if T_2 is not None: self.T_2 = T_2
        if Omega is not None: self.Omega = Omega
        if  k_1 is not None: self.k_1 = k_1
        if  k_2 is not None: self.k_2 = k_2
        if  phi_1 is not None: self.phi_1 = phi_1
        if  phi_2 is not None: self.phi_2 = phi_2 
        if  phi is not None: self.phi = phi 
        if  theta is not None: self.theta = theta
        if  c2 is not None: self.c2 = c2
        if  c1 is not None: self.c1 = c1
        #if d is not None: self.d=d
        self._init_from_components(**kwargs)
            
    @property
    def components(self):
        components = {}

        self.qs = [self.phi_1,self.phi_2]
        self.k_1 = (self.G*self.I_1)/self.l_1
        self.k_2 = (self.G*self.I_2)/self.l_2

        self.disc_1 = Disk(self.I_m1, pos1=self.phi_1, qs=self.qs)
        self.spring_1 = Spring(self.k_1, self.phi_1, self.phi_2, qs=self.qs)  # left spring
        self.disc_2 = Disk(self.I_m2, pos1=self.phi_2, qs=self.qs)
        self.spring_2 = Spring(self.k_2, pos1=self.phi_2, pos2=self.theta*cos(self.Omega*self.ivar),qs=self.qs)  # right spring
        self.moment_disc1=Force(self.T_1, pos1=self.phi_1)
        self.moment_disc2=Force(self.T_2, pos1=self.phi_2)
        self.damper_1 = Damper(self.c1, pos1=self.phi_1, pos2=self.phi_2, qs=self.qs)
        self.damper_2 = Damper(self.c2, pos1=self.phi_2, pos2=self.phi_2, qs=self.qs)
        

        components['disc_1'] = self.disc_1
        components['disc_2'] = self.disc_2
        components['spring_1'] = self.spring_1
        components['spring_2'] = self.spring_2
        components['moment_disc1'] = self.moment_disc1
        components['moment_disc2'] = self.moment_disc2
        components['damper_1'] = self.damper_1
        components['damper_2'] = self.damper_2

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict


    def get_default_data(self):

        m0, l0 , G, l,  d, T0, lam = self.m0, self.l0 , self.G, self.l,  self.d, self.T0 ,Symbol('lambda')
        theta0, Omega = self.theta0, self.Omega

        default_data_dict = {
            self.I_m1: [S.Half*m0*(l0**2)*no for no in range(1,9)],
            self.I_m2: [S.Half*m0*(l0**2)*no for no in range(1,9)],
            self.I_1: [S.Half**(no)*(l0**4) for no in range(1,9)],
            self.I_2: [S.Half**no*(l0**4) for no in range(1,9)],
            self.l_1: [S.Half**(no-6)*l0 for no in range(1,9)],
            self.l_2: [S.Half**(no-6)*l0 for no in range(1,9)],
            self.T_1: [T0 * (no) for no in range(1,12)],
            self.T_2: [T0 * (no) for no in range(1,12)],
            #self.theta:[theta0* cos(self.Omega * self.ivar) ],
            self.c1: [self.k_1 * (no) * lam for no in range (1,9)],
            self.c2: [self.k_2 * (no) * lam for no in range (1,9)]
            
        }

        return default_data_dict

    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST
        ]

        return comp_list
    
    def disc_1_force(self):
        t=self.ivar
        return self.I_m1 * self.steady_solution().diff(t,t)

    def max_static_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[0]/d)
    
    def max_dynamic_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[0]*self.Omega**2

        return  abs(2*(self.I_m1*acc_amp)/d) + self.max_static_disk_1_bearing_force()#.subs(self._given_data)

    def max_static_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[1]/d)
    
    def max_dynamic_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        acc_amp = self._frf()[1]*self.Omega**2

        return  abs(2*(self.I_m2*acc_amp)/d) + self.max_static_disk_2_bearing_force()#.subs(self._given_data)
    
    
    def static_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_1_bearing_force())/(kd*h)
    
    def dynamic_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_1_bearing_force())/(kd*h)
    
    
    def static_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_2_bearing_force())/(kd*h)
    
    def dynamic_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_2_bearing_force())/(kd*h)
    
    
#TO_DO 158 AF
class TripleShaft(ComposedSystem):
    """
     Triple Shaft system representing a kinematically excited shaft with three disks.
     Describes a physical model of a triple elastic shaft system with three rotating discs subjected to external moments. The class simulates the dynamics of such a system by defining the relevant parameters, dynamics and components.
    Arguments:
    =========
            m = Mass
            - The total mass of the shaft. (Całkowita masa wału.)

        m_1, m_2, m_3 = Masses of the Disks
            - The masses of the first, second, and third disks, respectively. 
              (Masy pierwszego, drugiego i trzeciego dysku.)

        l_1, l_2, l_3 = Lengths of the Shaft Segments
            - The lengths of the shaft segments between the disks. 
              (Długości segmentów wału pomiędzy dyskami.)

        d, d_1, d_2, d_3 = Diameters of the Shaft and Disks
            - The diameters of the shaft and each disk, respectively.
              (Średnica wału i każdego z trzech dysków.)

        G = Kirchhoff module

        M_1, M_2, M_3 = Applied Moments
            - The moments applied to the first, second, and third disks, respectively.
              (Moment siły przyłożony do pierwszego, drugiego i trzeciego dysku.)

        Omega = Angular Frequency
            - The angular frequency of the system's motion. (Częstotliwość kątowa ruchu układu.)

        delta = Phase Shift
            - The phase shift of the oscillatory input. (Przesunięcie fazowe wejściowego sygnału oscylacyjnego.)

        input_displacement = Angular Displacement
            - The angular displacement of the shaft input. (Przemieszczenie kątowe wejścia wału.)

        ivar = Symbol object
            - The independent time variable. (Zmienna niezależna czasu.)

        qs = Dynamicsymbol object
            - Generalized coordinates of the system. (Uogólnione współrzędne układu.)
            


    Example
    =======
   A shaft with three disks oscillating due to an applied moment and connected by springs.
 
    >>> t = symbols('t')
    >>> m, m_1, m_2, m_3 = symbols('m m_1 m_2 m_3')
    >>> l_1, l_2, l_3 = symbols('l_1 l_2 l_3')
    >>> G, Omega, delta = symbols('G Omega delta')
    >>> qs = dynamicsymbols('phi_1, phi_2, phi_3') # Generalized Coordinates
    >>> TripleShaft()

    -defines the symbols and dynamicsymbols
    -finally determines the instance of the Triple Shaft system using the TripleShaft class.
    """

    scheme_name = 'MDoFTripleShaft.PNG'
    real_name = 'MDoFTripleShaft_real.jpg'

    # Symbol definitions with default positive values
    m=Symbol('m', positive=True)
    m_1=Symbol('m_1', positive=True)
    m_2=Symbol('m_2', positive=True)
    m_3=Symbol('m_3', positive=True)
    l_1=Symbol('l_1', positive=True)
    l_2=Symbol('l_2', positive=True)
    l_3=Symbol('l_3', positive=True)
    d=Symbol('d', positive=True)
    d_1=Symbol('d_1', positive=True)
    d_2=Symbol('d_2', positive=True)
    d_3=Symbol('d_3', positive=True)
    G = Symbol('G', positive=True)
    M_1=Symbol('M_1',positive=True)
    M_2=Symbol('M_2',positive=True)
    M_3=Symbol('M_3',positive=True)
    Omega=Symbol('Omega',positive=True)
    delta =Symbol('delta', positive=True)
    
    # Dynamicsymbols for generalized coordinates and input displacement
    input_displacement=dynamicsymbols('theta')
    phi_1=dynamicsymbols('\\varphi_1')
    phi_2=dynamicsymbols('\\varphi_2')
    phi_3=dynamicsymbols('\\varphi_3')
    
    phi_r=dynamicsymbols('\\varphi_r')
    phi_l=dynamicsymbols('\\varphi_l')

    def __init__(self,
                 m = symbols('m', positive=True),
                 m_1=Symbol('m_1', positive=True),
                 m_2=Symbol('m_2', positive=True),
                 m_3=Symbol('m_3', positive=True),
                 l_1=Symbol('l_1', positive=True),
                 l_2=Symbol('l_2', positive=True),
                 l_3=Symbol('l_3', positive=True),
                 d=Symbol('d', positive=True),
                 d_1=Symbol('d_1', positive=True),
                 d_2=Symbol('d_2', positive=True),
                 d_3=Symbol('d_3', positive=True),
                 G=Symbol('G', positive=True),
                 M_1=Symbol('M_1', positive=True),
                 M_2=Symbol('M_2', positive=True),
                 M_3=Symbol('M_3', positive=True),
                 Omega=Symbol('\\Omega', positive=True),
                 delta=Symbol('\\delta', positive=True),

                 input_displacement=dynamicsymbols('theta'),
                 phi_1=dynamicsymbols('\\varphi_1'),
                 phi_2=dynamicsymbols('\\varphi_2'),
                 phi_3=dynamicsymbols('\\varphi_3'),

                 phi_l=dynamicsymbols('\\varphi_L'),
                 phi_r=dynamicsymbols('\\varphi_R'),

                 ivar=Symbol('t'),
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2, \\varphi_3,'),
                 **kwargs):

        # Initializing class attributes with provided arguments or default values
        if  m is not None: self.m = m
        if  m_1 is not None: self.m_1 = m_1
        if  m_2 is not None: self.m_2 = m_2
        if  m_3 is not None: self.m_3 = m_3
        if  l_1 is not None: self.l_1 = l_1
        if  l_2 is not None: self.l_2 = l_2
        if  l_3 is not None: self.l_3 = l_3
        if  d is not None: self.d = d
        if  d_1 is not None: self.d_1 = d_1
        if  d_2 is not None: self.d_2 = d_2
        if  d_3 is not None: self.d_3 = d_3
        if  G is not None: self.G = G
        if  M_1 is not None: self.M_1 = M_1
        if  M_2 is not None: self.M_2 = M_2
        if  M_3 is not None: self.M_3 = M_3
        if Omega is not None: self.Omega = Omega
        if delta is not None: self.delta = delta
        if  phi_1 is not None: self.phi_1 = phi_1
        if  phi_2 is not None: self.phi_2 = phi_2
        if  phi_3 is not None: self.phi_3 = phi_3
        if  phi_l is not None: self.phi_l = phi_l
        if  phi_r is not None: self.phi_r = phi_r
        if  input_displacement is not None: self.input_displacement = input_displacement

        self.qs = [self.phi_1,self.phi_2,self.phi_3]
               
        self.theta = input_displacement
        self._init_from_components(**kwargs)

        
        


    @property
    def components(self):
        components = {}

        # Moment of inertia calculations
        I_0=S.Half*pi*(self.d/2)**4
        I_1=S.Half*self.m_1*(self.d_1/2)**2
        I_2=S.Half*self.m_2*(self.d_2/2)**2
        I_3=S.Half*self.m_3*(self.d_3/2)**2
        
        # Assigning calculated moments of inertia to class attributes
        self.I_m0 = I_0
        self.I_m1 = I_1
        self.I_m2 = I_2
        self.I_m3 = I_3
        
        # Spring coefficients
        self.k_1=(self.G*I_0)/self.l_1
        self.k_2=(self.G*I_0)/self.l_2
        self.k_3=(self.G*I_0)/self.l_3

        # Defining system components
        self.disc_1 = Disk(self.I_m1, pos1=self.phi_1, qs=self.qs)
        self.spring_1 = Spring(self.k_1, self.phi_1, self.phi_2, qs=self.qs)  # left spring
        
        self.disc_2 = Disk(self.I_m2, pos1=self.phi_2, qs=self.qs)
        self.spring_2 = Spring(self.k_2, pos1=self.phi_2, pos2=self.theta,qs=self.qs)  # right spring
               
        self.disc_3 = Disk(self.I_m3, pos1=self.phi_3, qs=self.qs)
        self.spring_3 = Spring(self.k_3, pos1=self.phi_3, pos2=self.theta,qs=self.qs)
               
        self.moment_disc1=Force(self.M_1, pos1=self.phi_1)
        self.moment_disc2=Force(self.M_2, pos1=self.phi_2)
        self.moment_disc3=Force(self.M_3, pos1=self.phi_3)

        components['disc_1'] = self.disc_1
        components['disc_2'] = self.disc_2
        components['disc_3'] = self.disc_3
        components['spring_1'] = self.spring_1
        components['spring_2'] = self.spring_2
        components['spring_3'] = self.spring_3
        components['moment_disc1'] = self.moment_disc1
        components['moment_disc2'] = self.moment_disc2
        components['moment_disc3'] = self.moment_disc3

        return components
#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.I: r'Moment of Inertia',
#             self.k_1: r'',
#             self.k_2: r'',
#         }
#         return self.sym_desc_dict
    def get_numerical_data(self):

        default_data_dict = {

            
            
            self.m_1: [1],
            self.m_2: [1],
            self.m_3: [1],

            self.l_1: [1],
            self.l_2: [1],
            self.l_3: [2],

            self.d_1: [1],
            self.d_2: [1],
            self.d_3: [1],
            self.d: [1],

            self.G:[1],
            
            self.delta:[0],

            self.input_displacement:[0],
            
            self.theta:[0],
            
            self.M_1:[1],
            self.M_2:[1],
            self.M_3:[1],
        }
        return default_data_dict
    def get_default_data(self):


        m0, M0, l0 , d0 ,M0= symbols('m_0 M_0 l_0 d_0 M_0', positive=True)
        theta0, Omega = symbols('theta_0 Omega', positive=True)

        default_data_dict = {
            
            Symbol('I_m'):[S.Half*Symbol('m')*Symbol('R')**2],
            Symbol('I_S'):[S.Half*pi*Symbol('R_S')**4],
            Symbol('k_S'):[Symbol('G')*Symbol('I_S')/Symbol('L_S')],
            
            
            self.m_1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m_2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m_3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_3: [2 * l0, 4 * l0, S.Half * l0, 2 * l0, S.Half * l0],

            self.d_1: [1 * d0, 2 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d_2: [1 * d0, 2 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d_3: [2 * d0, 4 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d: [2 * d0, 4 * d0, S.Half * d0, 2 * d0, S.Half * d0],

#             self.phi_1:[self.phi_l,0],
#             self.phi_2:[self.phi_l,self.phi_r],
#             self.phi_3:[self.phi_r,0],
            
            self.delta:[0],

            self.input_displacement:[0],
            
            self.theta:[0],
            
            self.M_1:[M0],
            self.M_2:[M0],
            self.M_3:[M0],
            

        }

        return default_data_dict


    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi_1] != self.phi_l or parameters_dict[self.phi_2] != self.phi_l:

            parameters_dict[self.phi_1] = self.phi_l

        if parameters_dict[self.phi_2] != self.phi_r or parameters_dict[self.phi_3] != self.phi_r:

            parameters_dict[self.phi_3] = self.phi_r
            

        return parameters_dict
    
    def max_static_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[0]/d)
    
    def max_dynamic_disk_1_bearing_force(self):
        d=Symbol('d',positive=True)

        
        frf = self._frf()
        acc_amp = (self.phi_1).subs(self._given_data).subs({self.phi_l:frf[0],self.phi_r:frf[1]})
        

        return  abs(2*(self.I_m1*acc_amp)/d) + self.max_static_disk_1_bearing_force()#.subs(self._given_data)

    def max_static_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)
        return abs(2*self.static_load().doit()[1]/d)
    
    def max_dynamic_disk_2_bearing_force(self):
        d=Symbol('d',positive=True)

        
        frf = self._frf()
        acc_amp = (self.phi_2).subs(self._given_data).subs({self.phi_l:frf[0],self.phi_r:frf[1]})

        return  abs(2*(self.I_m2*acc_amp)/d) + self.max_static_disk_2_bearing_force()#.subs(self._given_data)
    
    
    def static_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_1_bearing_force())/(kd*h)
    
    def dynamic_disk_1_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_1_bearing_force())/(kd*h)
    
    
    def static_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_static_disk_2_bearing_force())/(kd*h)
    
    def dynamic_disk_2_key_length(self):
        kd=Symbol('k_d', positive=True)
        h=Symbol('h', positive=True)
        return (2*self.max_dynamic_disk_2_bearing_force())/(kd*h)
    
    def unit_dict(self):
        units_dict = {
            self.m: ureg.kilogram,  # Mass of the shaft
            self.m_1: ureg.kilogram,  # Mass of the first disk
            self.m_2: ureg.kilogram,  # Mass of the second disk
            self.m_3: ureg.kilogram,  # Mass of the third disk
            self.l_1: ureg.meter,  # Length of the first shaft segment
            self.l_2: ureg.meter,  # Length of the second shaft segment
            self.l_3: ureg.meter,  # Length of the third shaft segment
            self.d: ureg.meter,  # Diameter of the shaft
            self.d_1: ureg.meter,  # Diameter of the first disk
            self.d_2: ureg.meter,  # Diameter of the second disk
            self.d_3: ureg.meter,  # Diameter of the third disk
            self.G: ureg.pascal,  # Shear modulus (Moduł Kirchhoffa)
            self.M_1: ureg.newton * ureg.meter,  # Moment applied on the first disk (torque)
            self.M_2: ureg.newton * ureg.meter,  # Moment applied on the second disk (torque)
            self.M_3: ureg.newton * ureg.meter,  # Moment applied on the third disk (torque)
            self.Omega: ureg.radian / ureg.second,  # Angular velocity (rad/s)
            self.delta: ureg.radian,  # Phase difference (radians)
            self.input_displacement: ureg.radian,  # Input angular displacement (radians)
            self.phi_1: ureg.radian,  # Angular displacement of the first disk (radians)
            self.phi_2: ureg.radian,  # Angular displacement of the second disk (radians)
            self.phi_3: ureg.radian,  # Angular displacement of the third disk (radians)
            self.phi_r: ureg.radian,  # Angular displacement on the right side (radians)
            self.phi_l: ureg.radian,  # Angular displacement on the left side (radians)
            self.ivar: ureg.second,  # Time (seconds)
            }
        
        return units_dict

#Grześ
class SDoFDampedShaft(ComposedSystem):
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

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components={}

        self.disc_1 = Disk(self.I, pos1=self.phi, qs=self.qs)
        self.spring_1 = Spring(self.k_1 * self.k_2 / (self.k_2 + self.k_1), pos1=self.phi, pos2=self.theta, qs=self.qs)
        self.moment = Force(self.Ms, pos1=self.phi, qs=self.qs)
        self.damper_1 = Damper(self.c_1 * self.c_2/(self.c_2 + self.c_1), pos1=self.phi, pos2=self.theta, qs=self.qs)

        components['_disc_1'] = self.disc_1
        components['_spring_1'] = self.spring_1
        components['_moment'] = self.moment
        components['_damper_1']=self.damper_1

        return components

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
