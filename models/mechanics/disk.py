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

from functools import cached_property, lru_cache


from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST

from pint import UnitRegistry
ureg = UnitRegistry()



class RollingDisk(ComposedSystem):
    m1 = Symbol('m_1', positive=True)
    R = Symbol('R', positive=True)
    ivar = Symbol('t')
    x = dynamicsymbols('x')
    qs = dynamicsymbols('x')

    def __init__(self,
                 m1=None,
                 R=None,
                 x=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m1 is not None: self.m1 = m1
        if R is not None: self.R = R
        if x is not None: self.x = x

        if qs is None:
            self.qs = [self.x]
        else:
            self.qs = qs

        self.ivar=ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._disk_lin = MaterialPoint(self.m1, self.x, qs=self.qs)
        self._disk_rot = MaterialPoint(self.m1 / 2 * self.R**2, self.x / self.R, qs=self.qs)


        components['_disk_lin'] = self._disk_lin
        components['_disk_rot'] = self._disk_rot

        return components

#     def get_default_data(self):

#         m0, l0 = symbols('m_0 l_0', positive=True)

#         default_data_dict = {
#             self.m1: [S.One * m0 * no for no in range(1,20)],
#             self.R: [S.Half * l0 * no for no in range(1,20)],
#         }

#         return default_data_dict


class DiskMountingBlock(ComposedSystem):

    m = Symbol('m', positive=True)
    m1 = Symbol('m_1', positive=True)
    R = Symbol('R', positive=True)
    ivar = Symbol('t')
    x = dynamicsymbols('x')
    qs = dynamicsymbols('x')

    def __init__(self,
                 m=None,
                 m1=None,
                 R=None,
                 x=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):


        if m is not None: self.m = m
        if m1 is not None: self.m1 = m1
        if R is not None: self.R = R
        if x is not None: self.x = x
        if qs is None:
            self.qs = [self.x]
        else:
            self.qs = qs

        self.ivar=ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._disk = RollingDisk(m1=self.m, R=self.R, x=self.x, qs=[self.x], ivar=self.ivar)
        self._disk_inner = MaterialPoint(self.m1, self.x, qs = [self.x])


        components['_disk'] = self._disk
        components['_disk_inner'] = self._disk_inner

        return components

#     def get_default_data(self):

#         m0 = symbols('m_0', positive=True)

#         default_data_dict = {
#             self.m1: [S.One * m0 * no for no in range(1,20)],
#         }

#         return default_data_dict

#     def get_numerical_data(self):

#         m0 = symbols('m_0', positive=True)


#         default_data_dict = {
#             self.m1: [
#                 0.5 * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
#                 7 * m0, 8 * m0, 9 * m0
#             ],
#         }

#         return default_data_dict



class DiskWithNonlinearSpring(NonlinearComposedSystem):
    scheme_name = 'nonlinear_disc.png'
    real_name = 'roller_tightener.png'

    m1 = Symbol('m1', positive=True)
    kl = Symbol('k', positive=True)
    R = Symbol('R', positive=True)
    d = Symbol('d', positive=True)
    l_0 = Symbol('l_0', positive=True)
    ivar = Symbol('t')
    x = dynamicsymbols('x')

    def __init__(self,
                 m1=None,
                 kl=None,
                 R=None,
                 d=None,
                 l_0=None,
                 x=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m1 is not None: self.m1 = m1
        if kl is not None: self.kl = kl
        if R is not None: self.R = R
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if x is not None: self.x = x
        self.ivar=ivar
        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._disk = RollingDisk(self.m1, self.R, self.x, self.qs, self.ivar)
        self._spring = Spring(self.kl, pos1=(sqrt(self.x**2 + self.d**2) - self.l_0), qs=self.qs)


        components['_disk'] = self._disk
        components['_spring'] = self._spring
        return components


    def get_default_data(self):

        m0, k0, l0, F0 = symbols('m_0 k_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.m1: [m0 * no for no in range(5, 15)],
            self.d: [self.l_0*S.Half*no for no in range(4,16)]
        }

        return default_data_dict



class ForcedDiskWithNonlinearSpring(NonlinearComposedSystem):
    scheme_name = 'nonlinear_disc.png'
    real_name = 'roller_tightener.png'

    m1 = Symbol('m1', positive=True)
    kl = Symbol('k', positive=True)
    R = Symbol('R', positive=True)
    d = Symbol('d', positive=True)
    l_0 = Symbol('l_0', positive=True)
    F = Symbol('F', positive=True)
    Omega = Symbol('Omega', positive=True)
    ivar = Symbol('t')
    x = dynamicsymbols('x')

    def __init__(self,
                 m1=None,
                 kl=None,
                 R=None,
                 d=None,
                 l_0=None,
                 F=None,
                 Omega=None,
                 x=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m1 is not None: self.m1 = m1
        if kl is not None: self.kl = kl
        if R is not None: self.R = R
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if F is not None: self.F = F
        if Omega is not None: self.Omega = Omega
        if x is not None: self.x = x
        self.ivar=ivar
        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._disk = DiskWithNonlinearSpring(self.m1, self.kl, self.R, self.d, self.l_0, self.x, self.ivar)
        self._force = Force(self.F * cos(self.Omega * self.ivar), pos1=self.x, qs=self.qs)


        components['_disk'] = self._disk
        components['_force'] = self._force

        return components


class DampedDiskWithNonlinearSpring(NonlinearComposedSystem): #trzeba zupdateowac scheme - dodac tlumienie liniowe
    scheme_name = 'nonlinear_disc.png'
    real_name = 'roller_tightener.png'

    c = Symbol('c', positive=True)
    m1 = Symbol('m1', positive=True)
    kl = Symbol('k', positive=True)
    R = Symbol('R', positive=True)
    d = Symbol('d', positive=True)
    l_0 = Symbol('l_0', positive=True)
    F = Symbol('F', positive=True)
    Omega = Symbol('Omega', positive=True)
    ivar = Symbol('t')
    x = dynamicsymbols('x')
    qs = dynamicsymbols('x')

    def __init__(self,
                 m1=None,
                 c=None,
                 kl=None,
                 R=None,
                 d=None,
                 l_0=None,
                 F=None,
                 x=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if c is not None: self.c = c
        if m1 is not None: self.m1 = m1
        if kl is not None: self.kl = kl
        if R is not None: self.R = R
        if d is not None: self.d = d
        if l_0 is not None: self.l_0 = l_0
        if F is not None: self.F = F
        if x is not None: self.x = x

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._disk = ForcedDiskWithNonlinearSpring(self.m1, self.kl, self.R, self.d, self.l_0, self.F, self.Omega, self.x, self.ivar)
        self._damper = Damper(self.c, pos1 = self.x, pos2 = 0 , qs=[self.x])
        

        components['_disk'] = self._disk
        components['damper'] = self._damper
        
        return components



class ForcedNonLinearDisk(ForcedDiskWithNonlinearSpring):
    pass


class ForcedNonLinearDiskSpring(DampedDiskWithNonlinearSpring):
    pass



### Half Disk classes implementation


class RollingHalfDisk(ComposedSystem):
    """
    Model of a Single Degree of Freedom - RollingHalfDisc

        Arguments:
        =========
            m = Mass
                -Mass of the Half Disc

            g = gravitional field
                -value of gravitional field acceleration

            r = radius
                -radius of the Half Disc

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A half disc of radius R with optional horisontal force.

        >>> t = symbols('t')
        >>> R = symbols(' R ',positive=True)
        >>> RollingHalfDisc(r=R)

        -We define the symbols and dynamicsymbols
        -determine the instance of the RollingHalfDisc by using class RollingHalfDisc()
    """

    scheme_name = 'rolling_half_disk.png'
    real_name = 'rolling_half_disk.png'
    r=Symbol('r', positive=True)
    m=Symbol('m', positive=True)
    g=Symbol('g', positive=True)
    phi=dynamicsymbols('varphi')
    qs=[phi]

    def __init__(self,
                 r=None,
                 m=None,
                 g=None,
                 ivar=Symbol('t'),
                 phi=None,
                 qs=None,
                 **kwargs):


        if r is not None: self.r = r
        if m is not None: self.m = m
        if g is not None: self.g = g

        if phi is not None: self.phi = phi
        self.ivar=ivar
        
        if qs is not None: self.qs = qs
        
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        d = 4*self.r/3/pi
        self.y = self.r - d*cos(self.phi)
        self.x =  self.r*self.phi - d*sin(self.phi)

        self.material_point_1 = MaterialPoint(self.m, self.x,ivar=self.ivar, qs=[self.phi])
        self.material_point_2 = MaterialPoint(self.m, self.y,ivar=self.ivar, qs=[self.phi])
        self.material_point_3 = MaterialPoint(((self.m*(self.r**2))/2)-self.m*d**2, self.phi,ivar=self.ivar, qs=[self.phi])
        self.gravity = GravitationalForce(self.m, self.g, pos1=self.y,ivar=self.ivar, qs=[self.phi])

        components['material_point_1'] = self.material_point_1
        components['material_point_2'] = self.material_point_2
        components['material_point_3'] = self.material_point_3
        components['gravity'] = self.gravity

        return components

    def get_default_data(self):

        m0, r0 = symbols('m_0 r_0', positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(2,10)],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
        }
        return default_data_dict

    def get_numerical_data(self):
        m0, l0, r0 = symbols('m_0 l_0 r_0', positive=True)

        default_data_dict = {
            self.m: [1],
            self.r: [1],
            #self.pi: [3.14],
        }
        return default_data_dict
    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r'Half disc radius',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict
    
    def unit_dict(self):

        unit_dict = {

            self.r: ureg.meter,
            self.m: ureg.kilogram,
            self.g: ureg.meter/ureg.second/ureg.second,
            self.ivar: ureg.second,
            self.phi: ureg.radian,
            self.phi.diff(self.ivar): ureg.radian/ureg.second,
            self.phi.diff(self.ivar,2): ureg.radian/ureg.second/ureg.second,
        }
        return unit_dict

class ForcedRollingHalfDisk(RollingHalfDisk):
    """
    Model of a Forced Single Degree of Freedom - ForcedRollingHalfDisk

        Arguments:
        =========
            m = Mass
                -Mass of the Half Disc

            g = gravitional field
                -value of gravitional field acceleration

            r = radius
                -radius of the Half Disc
            
            M_0 = torque amplitude
                -amplitude value of horisontal force
            
            Omega = Frequency of force
                -value of frequency force
            
            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A half disc of radius R with optional horisontal force.

        >>> t = symbols('t')
        >>> R, Omega, M_0 = symbols(' R Omega M_0 ',positive=True)
        >>> ForcedRollingHalfDisk(r=R, Omega=Omega, M_0=M_0 )

        -We define the symbols and dynamicsymbols
        -determine the instance of the ForcedRollingHalfDisk by using class ForcedRollingHalfDisk()
    """

    scheme_name = 'forced_rolling_half_disk.png'
    real_name = 'forced_rolling_half_disk.png'
    r=Symbol('r', positive=True)
    l=Symbol('l', positive=True)
    m=Symbol('m', positive=True)
    g=Symbol('g', positive=True)
    Omega=Symbol('Omega', positive=True)
    M_0=Symbol('M_0', positiv=True)
    phi=dynamicsymbols('varphi')

    def __init__(self,
                 r=None,
                 m=None,
                 g=None,
                 M_0=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 phi=None,
                 **kwargs):

        if r is not None: self.r = r
        if m is not None: self.m = m
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if M_0 is not None: self.M_0 = M_0
        if Omega is not None: self.Omega = Omega
        self.ivar=ivar
        self.qs=[self.phi]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        #self.rolling_half_disk.x
        d = 4*self.r/3/pi
        #self.y = self.r - d*cos(self.phi)
        x =  self.r*self.phi - d*sin(self.phi)
        self.rolling_half_disk = RollingHalfDisk(m=self.m, r=self.r, phi=self.phi, qs=self.qs)
        self.torque = Force(self.M_0*cos(self.Omega*self.ivar), self.phi, qs=self.qs)(label='force')

        components['rolling_half_disk'] = self.rolling_half_disk
        components['torque'] = self.torque
        
        return components

    def get_default_data(self):

        m0, r0 = symbols('m_0 r_0', positive=True)

        default_data_dict = {
            self.m: [S.One * m0 * no for no in range(2,10)],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
        }
        return default_data_dict

    def get_numerical_data(self):
        m0, l0, r0 = symbols('m_0 l_0 r_0', positive=True)

        default_data_dict = {
            self.m: [1],
            self.r: [1],
            self.M_0: [1],
            self.Omega: [1],
            #self.pi: [3.14],
        }
        return default_data_dict
    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r'Half disc radius',
            self.m: r'Mass',
            self.g: 'Gravity constant',
            self.M_0: 'Torque amplitude value',
            self.Omega: 'Frequency of torque'
        }
        return self.sym_desc_dict

    def unit_dict(self):

        unit_dict = {

            self.r: ureg.meter,
            self.m: ureg.kilogram,
            self.g: ureg.meter/ureg.second/ureg.second,
            self.ivar: ureg.second,
            self.phi: ureg.radian,
            self.phi.diff(self.ivar): ureg.radian/ureg.second,
            self.phi.diff(self.ivar,2): ureg.radian/ureg.second/ureg.second,
            self.M_0: ureg.newton*ureg.meter,
            self.Omega: ureg.radian/ureg.second,
        }
        return unit_dict
