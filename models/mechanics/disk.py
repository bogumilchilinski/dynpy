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