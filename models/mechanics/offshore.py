from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from  ..continuous import ContinuousSystem, PlaneStressProblem

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin

class CompoundSystem(ComposedSystem):

    z = dynamicsymbols('z')
    _p = Symbol('p')


    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self._p, self.qs[0],
                                             self.qs)('Material Point')
        components['_material_point'] = self._material_point

        return components

    
class SDOFSemiSubmersiblePlatform(ComposedSystem):

    #scheme_name = '.PNG'
    #real_name = '.jpg'
    
    z=dynamicsymbols('z')
    m=Symbol('m', positive=True)
    m_p=Symbol('m_p', positive=True)
    k_1=Symbol('k', positive=True)
    c_1=Symbol('c', positive=True)
    c_2=Symbol('c_water', positive=True)
    F_load=Symbol('F_{load}', positive=True)
    g=Symbol('g', positive=True)
    p_p=Symbol('P_p', positive=True)
    rho=Symbol('\\rho', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    
    def __init__(self,
                 m=None,
                 m_p=None,
                 k_1=None,
                 c_1=None,
                 c_2=None,
                 F_load=None,
                 ivar=Symbol('t'),
                 z=None,
                 g=None,
                 p_p=None,
                 rho=None,
                 Omega=None,
                 F=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if m is not None: self.m = m
        if k_1 is not None: self.k_1 = k_1
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if F_load is not None: self.F_load = F_load
        if g is not None: self.g=g
        if rho is not None: self.rho=rho
        if p_p is not None: self.p_p=p_p
        if Omega is not None: self.Omega=Omega
        if F is not None: self.F = F
        self.ivar=ivar
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}
        
        self._platform = MaterialPoint(self.m, pos1=self.z, qs=[self.z])(label='Material point - mass of the platform')
        self._spring_1 = Spring(self.k_1, pos1=self.z, pos2=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Spring - stiffness of the platform')
        self._spring_2 = Spring(self.p_p*self.rho, pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Spring - buoyancy')
        self._damper = Damper(self.c_1, pos1=self.z, pos2=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Damper - damping of the platform')
        self._water = Damper(self.c_2, pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Damper - damping of the water')
        self._force = Force(-self.F_load, pos1=self.z, qs=[self.z])(label='Force - load on the platform')
        self._MassOfPlatform = Force(-self.m*self.g, pos1=self.z, qs=[self.z])(label='Force - weight of the platform in gravitational field')
        self._excitation = Force(self.F*sin(self.Omega*self.ivar), pos1=(self.k_1*self.z)/(self.p_p*self.rho), qs=[self.z])(label='Force - Force caused by the waves')
        
        components['_platform'] = self._platform
        components['_spring_1'] = self._spring_1
        components['_spring_2'] = self._spring_2
        components['_damper'] = self._damper
        components['_water'] = self._water
        components['_force'] = self._force
        components['_MassOfPlatform'] = self._MassOfPlatform
        components['_excitation'] = self._excitation
        
        return components


    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of a platform',
            self.g: r'Gravity field acceleration',
            self.Omega: r'Excitation frequency',
            self.c_1: r'Platform damping coefficient',
            self.c_2: r'Water damping coefficient',
            self.rho: r'Water density',
            self.k_1: r'Platform stiffness coefficient',
            self.p_p: r'Area submerged in water',
            self.F_load: r'Force caused by a load on a rig',
            self.F: r'Force caused by a waves',
            self.ivar: r'Time'
        }
        return self.sym_desc_dict
    
    
    def get_default_data(self):

        m0, Omega0, k0, p_p0, F0, lam0 = Symbol('m_0', positive=True), Symbol('Omega0', positive=True), Symbol('k0', positive=True), Symbol('p_{p0}', positive=True), Symbol('F0', positive=True), Symbol('lambda0', positive=True)
        c_2 = self.c_2
        g=self.g
        rho=self.rho

        default_data_dict = {
            self.m: [S.One * m0 * no * 1000 for no in range(20, 30)],
            self.g: [S.One * g * no for no in range(1, 2)],
            self.Omega: [S.One * Omega0 * no for no in range(1, 2)],
            self.c_2: [S.One * c_2 * no for no in range(1, 2)],
            self.rho: [S.One * rho * no for no in range(1, 2)],
            self.k_1: [S.One * k0 * no * 10000 for no in range(5, 10)],
            self.c_1: [S.One * k0 * lam0 * no * 10000 for no in range(5, 10)],
            self.p_p: [S.One * p_p0 * no for no in range(9, 16)],
            self.F_load: [S.One * F0 * no * 1000 for no in range(5, 10)],
            self.F: [S.One * F0 * no * 1000 for no in range(1, 10)]
        }

        return default_data_dict

    
    
class DDOFSemiSubmersiblePlatform(ComposedSystem):

    #scheme_name = '.PNG'
    #real_name = '.jpg'
    
    z_p=dynamicsymbols('z_p')
    z=dynamicsymbols('z')
    m=Symbol('m', positive=True)
    m_p=Symbol('m_p', positive=True)
    k_1=Symbol('k', positive=True)
    c_1=Symbol('c', positive=True)
    c_2=Symbol('c_p', positive=True)
    F_load=Symbol('F_{load}', positive=True)
    g=Symbol('g', positive=True)
    p_p=Symbol('P_p', positive=True)
    rho=Symbol('rho', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    
    def __init__(self,
                 m=None,
                 m_p=None,
                 k_1=None,
                 c_1=None,
                 c_2=None,
                 F_load=None,
                 ivar=Symbol('t'),
                 z=None,
                 z_p=None,
                 g=None,
                 p_p=None,
                 rho=None,
                 Omega=None,
                 F=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if z_p is not None: self.z_p=z_p
        if m is not None: self.m = m
        if m_p is not None: self.m_p = m_p
        if k_1 is not None: self.k_1 = k_1
        if c_1 is not None: self.c_1 = c_1
        if c_2 is not None: self.c_2 = c_2
        if F_load is not None: self.F_load = F_load
        if g is not None: self.g=g
        if rho is not None: self.rho=rho
        if p_p is not None: self.p_p=p_p
        if Omega is not None: self.Omega=Omega
        if F is not None: self.F=F
        self.ivar=ivar
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}
        
        self._platform = MaterialPoint(self.m, pos1=self.z, qs=[self.z])(label='Material point - mass of the platform')
        self._pontoon = MaterialPoint(self.m_p, pos1=self.z_p, qs=[self.z_p])(label='Material point - mass of the pontoon')
        self._spring_1 = Spring(self.k_1, pos1=self.z, pos2=self.z_p, qs=[self.z, self.z_p])(label='Spring - stiffness of the platform')
        self._spring_2 = Spring(self.p_p*self.rho, pos1=self.z_p, qs=[self.z_p, self.z])(label='Spring - buoyancy')
        self._damper = Damper(self.c_1, pos1=self.z, pos2=self.z_p, qs=[self.z, self.z_p])(label='Damper - damping of the platform')
        self._water = Damper(self.c_2, pos1=self.z_p, qs=[self.z_p, self.z])(label='Damper - damping of the water')
        self._force = Force(-self.F_load, pos1=self.z, qs=[self.z, self.z_p])(label='Force - load on the platform')
        self._MassOfPlatform = Force(-self.m*self.g, pos1=self.z, qs=[self.z, self.z_p])(label='Force - weight of the platform in gravitational field')
        self._MassOfPontoon = Force(-self.m_p*self.g, pos1=self.z_p, qs=[self.z_p, self.z])(label='Force - weight of the pontoon in gravitational field')
        self._excitation = Force(self.F*sin(self.Omega*self.ivar), pos1=self.z_p, qs=[self.z_p, self.z])(label='Force - Force caused by the waves')
        
        
        components['_platform'] = self._platform
        components['_spring_1'] = self._spring_1
        components['_spring_2'] = self._spring_2
        components['_damper'] = self._damper
        components['_water'] = self._water
        components['_force'] = self._force
        components['_MassOfPlatform'] = self._MassOfPlatform
        components['_MassOfPontoon'] = self._MassOfPontoon
        components['_pontoon'] = self._pontoon
        components['_excitation'] = self._excitation
        
        return components
       
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of a platform',
            self.g: r'Gravity field acceleration',
            self.Omega: r'Excitation frequency',
            self.c_1: r'Platform damping coefficient',
            self.c_2: r'Water damping coefficient',
            self.rho: r'Water density',
            self.k_1: r'Platform stiffness coefficient',
            self.p_p: r'Area submerged in water',
            self.F_load: r'Force caused by a load on a rig',
            self.F: r'Force caused by a waves',
            self.ivar: r'Time',
            self.m_p: r'Mass of a pontoon'
        }
        return self.sym_desc_dict
    
    
    def get_default_data(self):

        m0, Omega0, k0, p_p0, F0, lam0 = Symbol('m_0', positive=True), Symbol('Omega0', positive=True), Symbol('k0', positive=True), Symbol('p_{p0}', positive=True), Symbol('F0', positive=True), Symbol('lambda_0', positive=True)
        c_2 = self.c_2
        g=self.g
        rho=self.rho

        default_data_dict = {
            self.m: [S.One * m0 * no * 1000 for no in range(20, 30)],
            self.m_p: [S.One * m0 * no * 1000 for no in range(10, 20)],
            self.g: [S.One * g * no for no in range(1, 2)],
            self.Omega: [S.One * Omega0 * no for no in range(1, 2)],
            self.c_2: [S.One * c_2 * no for no in range(1, 2)],
            self.rho: [S.One * rho * no for no in range(1, 2)],
            self.k_1: [S.One * k0 * no * 10000 for no in range(5, 10)],
            self.c_1: [S.One * k0 * lam0 * no * 10000 for no in range(5, 10)],
            self.p_p: [S.One * p_p0 * no for no in range(9, 16)],
            self.F_load: [S.One * F0 * 1000 * no  for no in range(5, 10)],
            self.F: [S.One * F0 * 1000 * no for no in range(1, 10)]
        }

        return default_data_dict
