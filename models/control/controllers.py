from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Derivative,Integral, Expr,Function,latex, Heaviside, atan, pi)
from numbers import Number
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

import numpy as np

from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, GeometryScene,base_frame,base_origin
from ...utilities.templates import tikz

#from dgeometry import GeometryScene,np

from pylatex import TikZ,TikZNode
import matplotlib.pyplot as plt

import base64
import random
import IPython as IP

from ...models.elements import Force,IntegralElement,ProportionalElement,DerivativeElement

class PIController(Force):
    
    '''
    A proportional–integral controller (PI controller) is a control loop mechanism
    employing feedback used in industrial control systems requiring continuously
    modulated control. A PI controller calculates an error value (e(t)) as the
    difference between a desired setpoint (SP) and a measured process variable (PV)
    and applies a correction based on proportional and integral terms. It consists of
    two independent elements - P and I that being combined return two-term controller.
    '''
    
    k_P = Symbol('k_P', positive=True)
    k_I = Symbol('k_I', positive=True)
    error = Symbol('e', positive=True)
    target = Symbol('target', positive=True)
    reference = S.Zero

    def __init__(self,
                 k_P,
                 k_I,
                 error = None,
                 target=None,
                 reference = None,
                 qs=None,
                 ivar=Symbol('t'),
                 frame = base_frame,
                 **kwargs):

        if k_P is not None: self.k_P = k_P
        if k_I is not None: self.k_I = k_I
        
        if reference is not None: self.reference = reference

        if error is not None: self.error = error

        if target is None and isinstance(self.error,Function):
               target=self.error

        self.ivar = ivar
        self.qs = qs

        self._init_from_components(**kwargs)

    @property
    def components(self):
        
        components = {}
        k_P, k_I = self.k_P, self.k_I
        
        self._proporional_elem = ProportionalElement(k_P, self.error, self.reference, qs = self.qs, ivar = self.ivar)
        self._integral_elem = IntegralElement(k_I, self.error, self.reference, qs = self.qs, ivar = self.ivar)
        
        components['proportional_elem'] = self._proporional_elem
        components['integral_elem'] = self._integral_elem
        
        return components