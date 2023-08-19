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


from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin

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

        self.DiskLin = MaterialPoint(self.m1, self.x, qs=self.qs)
        self.DiskRot = MaterialPoint(self.m1 / 2 * self.R**2, self.x / self.R, qs=self.qs)


        components['_DiskLin'] = self.DiskLin
        components['_DiskRot'] = self.DiskRot

        return components

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.One * m0 * no for no in range(1,20)],
            self.R: [S.Half * l0 * no for no in range(1,20)],
        }

        return default_data_dict


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

#         self.DiskLin = MaterialPoint(self.m, self.x, qs=[self.x])
#         self.DiskRot = MaterialPoint(self.m / 2 * self.R**2, self.x / self.R, qs=[self.x])
        self.Disk = RollingDisk(m1=self.m, R=self.R, x=self.x, qs=[self.x], ivar=self.ivar)
        self.DiskInner = MaterialPoint(self.m1, self.x, qs = [self.x])

#         components['_DiskLin'] = self.DiskLin
#         components['_DiskRot'] = self.DiskRot
        components['_Disk'] = self.Disk
        components['_DiskInner'] = self.DiskInner

        return components

    def get_default_data(self):

        m0 = symbols('m_0', positive=True)

        default_data_dict = {
            self.m1: [S.One * m0 * no for no in range(1,20)],
        }

        return default_data_dict

    def get_numerical_data(self):

        m0 = symbols('m_0', positive=True)


        default_data_dict = {
            self.m1: [
                0.5 * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0,
                7 * m0, 8 * m0, 9 * m0
            ],
        }

        return default_data_dict

#dane numeryczne
class ForcedNonLinearDisk(NonlinearComposedSystem):
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

        self.qs = [self.x]
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self.disk1_lin = MaterialPoint(self.m1, pos1=self.x, qs=self.qs)
        self.disk1_rot= MaterialPoint(self.m1 / 2 * self.R**2, pos1=self.x/self.R, qs=self.qs)
        self.spring = Spring(self.kl, pos1=(sqrt(self.x**2 + self.d**2) - self.l_0), qs=self.qs)
        self.force = Force(self.F * cos(self.Omega * self.ivar), pos1=self.x, qs=self.qs)



        components['_disk1_lin'] = self.disk1_lin
        components['_disk1_rot'] = self.disk1_rot
        components['_spring'] = self.spring
        components['_force'] = self.force

        return components

    def get_default_data(self):

        m0, k0, l0, F0 = symbols('m_0 k_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.m1: [S.One * no * m0 for no in range(5,15)],
            self.d: [S.Half * no * self.l_0 for no in range (4,16)],
            self.kl: [S.One * k0 * no for no in range (20,30)],
            self.F: [S.One * F0 * no for no in range(10,20)]
        }

        return default_data_dict

    def get_numerical_data(self): ### Brakowało numerical data. Przez to komenda get numerical parameters nie działa

        #m0, k0, l0, c0= symbols('m_0 k_0 l_0 c_0', positive=True)
        m0, k0, l0, c0 = 100, 10e3, 0.5, 5e3

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
        }

        return default_data_dict
    
class ForcedNonLinearDiskSpring(NonlinearComposedSystem): ### Miałeś bład z tabulacją - pilnuj wcięć 
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

        self._disk = MaterialPoint(self.m1, self.x, qs=[self.x]) +  MaterialPoint(self.m1 / 2 * self.R**2, self.x / self.R, qs=[self.x])
        #self._disk = Disk(S.One/2 * self.m1*self.R**2,pos1 = self.x / self.R, qs=[self.x]) ??? returns diffrend eoms than soultion above
            
        self._spring = Spring(self.kl, pos1 = sqrt(self.x**2 + self.d**2), pos2 = - self.l_0 , qs=[self.x])
        
        self._damper = Damper(self.c, pos1 = self.l_0, pos2 = - sqrt(self.x**2 + self.d**2) , qs=[self.x]) #not sure about pos2
        
        self._force = Force(self.F * cos(self.Omega * self.ivar), pos1=self.x, qs=[self.x])


        components['disk'] = self._disk
        components['spring'] = self._spring
        components['force'] = self._force
        components['damper'] = self._damper
        
        return components

    def get_default_data(self):

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
#             self.c:  [
#                 1 * c0, 3 * c0, 2 * c0, 4 * c0, 5 * c0, 6 * c0, 7 * c0, 8 * c0,
#                 9 * c0
#             ],
            self.c:  [
                1 * l0, 3 * l0, 2 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0
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
