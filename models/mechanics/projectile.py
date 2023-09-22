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

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST

class MissileTrajectoryAirless(ComposedSystem):
    
    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    c=Symbol('c',positive=True)
    x=dynamicsymbols('x',positive=True)
    y=dynamicsymbols('y',positivie=True)
    c0=Symbol('C0',positive=True)
    ivar=Symbol('t')
    
    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 x=None,
                 y=None,
                 ivar=None,
                 **kwargs):
        
        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if x is not None: self.x = x
        if y is not None: self.y = y
        if ivar is not None: self.ivar= ivar

        self.qs = [self.x,self.y]
        
        self._init_from_components(**kwargs)
       
    @property
    def components(self):

        v_x=(self.x).diff(self.ivar)
        v_y=(self.y).diff(self.ivar)
        components = {}

        self._mass_x = MaterialPoint(self.m,pos1=self.x,qs=self.qs)
        self._mass_y = MaterialPoint(self.m,pos1=self.y,qs=self.qs)

        self._gravity = GravitationalForce(self.m,self.g,pos1=self.y,qs=self.qs)
        
        

        
        components['_mass_x'] = self._mass_x
        components['_mass_y'] = self._mass_y
        components['_gravity'] = self._gravity
        #components['_drag_x'] = self._drag_x
        #components['_drag_y'] = self._drag_y
        return components

#     def get_default_data(self):

#         m0, c0,= self.m0,self.c0

#         default_data_dict = {
#             self.m: [m0*no for no in range (1,8)],

#             self.c: [c0*no for no in range (0.3,1)],

#         }

#         return default_data_dict

    @property
    def _report_components(self):
        
        comp_list=[
       *REPORT_COMPONENTS_LIST
        ]
        
        return comp_list
    
class MissileTrajectory(MissileTrajectoryAirless):

    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    c=Symbol('c',positive=True)
    x=dynamicsymbols('x',positive=True)
    y=dynamicsymbols('y',positivie=True)
    c0=Symbol('c0',positive=True)
    ivar=Symbol('t')

    
    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 x=None,
                 y=None,
                 ivar=None,
                 **kwargs):
        
        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if x is not None: self.x = x
        if y is not None: self.y = y
        if ivar is not None: self.ivar= ivar

        self.qs = [self.x,self.y]
        
        self._init_from_components(**kwargs)


        
    @property
    def components(self):

        v_x=(self.x).diff(self.ivar)
        v_y=(self.y).diff(self.ivar)
        components = {}

        self._missile_airless = MissileTrajectoryAirless(self.m, x=self.x, y=self.y, qs=self.qs)

        self._drag_x = Damper(self.c*v_x,self.x,qs=self.qs)(label = 'horizontal drag')
         
        self._drag_y = Damper(self.c* v_y,self.y,qs=self.qs)(label = 'drag')
        
        

        
        components['_missile_airless'] = self._missile_airless
        components['_drag_x'] = self._drag_x
        components['_drag_y'] = self._drag_y
        
        return components

    def get_default_data(self):

        m0, c0,= self.m0,self.c0

        default_data_dict = {
            self.m: [m0*no for no in range (1,8)],

            self.c: [c0*no for no in range (1,8)],

        }
        
        return default_data_dict
        
    def get_numerical_data(self):

        default_data_dict = {
            self.m: [140],
            self.g: [10],
            self.c: [0.95],
        }

        return default_data_dict

    
class MissileTrajectoryAerodynamic(MissileTrajectory):
    
    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    c=Symbol('c',positive=True)
    x=dynamicsymbols('x',positive=True)
    y=dynamicsymbols('y',positivie=True)
    c=Symbol('c',positive=True)
    ro=Symbol('\\rho',positive=True)
    area=Symbol('A',positive=True)
    c_x=Symbol('c_x',positive=True)
    c0=Symbol('c0',positive=True)
    
    
    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 x=None,
                 y=None,
                 ca=None,
                 ro=None,
                 area=None,
                 c_x=None,
                 ivar=None,
                 **kwargs):
        
        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if x is not None: self.x = x
        if y is not None: self.y = y
        if ca is not None: self.ca = ca
        if ro is not None: self.ro = ro
        if area is not None: self.area = area
        if c_x is not None: self.c_x = c_x
        if ivar is not None: self.ivar = ivar

        self.qs = [self.x,self.y]
        
        self._init_from_components(**kwargs)
       
    @property
    def components(self):
        area = self.area
        v_x=(self.x).diff(self.ivar)
        v_y=(self.y).diff(self.ivar)
        ca = 1/2*self.ro*area*self.c_x*((((self.x).diff(self.ivar))**2+((self.y).diff(self.ivar))**2))
        components = {}

        self.missile = MissileTrajectory(m = self.m, g = self.g, c = ca, x=self.x, y=self.y, qs=self.qs)
 
        components['missile'] = self.missile

        return components

    def get_default_data(self):

        m0, c0,= self.m0,self.c0

        default_data_dict = {
            self.m: [m0*no for no in range (1,8)],

            self.c: [c0*no for no in range (1,8)],

        }
        
        return default_data_dict
        
#     def get_numerical_data(self):

#         default_data_dict = {
#             self.m: [140],
#             self.g: [10],
#             self.c: [0.95],
#             self.area:[1]
#             self.ro:[1]
#             selfc_x:[1]
#             }
        

        return default_data_dict