import importlib
from . import dynamics as dyn
importlib.reload(dyn)
from . import dynamics as dyn

from sympy.physics.vector import dynamicsymbols
from sympy import *
from sympy.physics.mechanics import *
import sympy as sym
import numpy as np
import pylatex
from pylatex import Command, NoEscape
import pint
from dynpy.dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp
from dynpy.models.elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from dynpy.models.mechanics.trolley import SpringMassSystem, SpringDamperMassSystem
# from dynpy.models.sdof import ComposedSystem, SpringMassSystem, DampedSpringMassSystem, BeamBridgeDamped
ureg = pint.UnitRegistry()
from dynpy.models.mechanics.principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin
mechanics_printing()



t=Symbol('t') #independent variable - time

class LEWEL16_simple(SpringMassSystem):
    scheme_name = 'train_sdof.png'
    real_name = 'lewel16.jpg'

class LEWEL16_simple_damped(SpringDamperMassSystem):
    scheme_name = 'train_sdof_damped.png'
    real_name = 'lewel16.jpg'

class LeafSpring(ComposedSystem):

    m=Symbol('m', positive=True)
    k_beam=Symbol('k_beam', positive=True)
    ivar=Symbol('t')
    g=Symbol('g', positive=True)
    l=Symbol('l', positive=True)
    n=Symbol('n',positive=True)
    t=Symbol('s',positive=True)
    w=Symbol('w',positive=True)
    module=Symbol('E', positive=True)
    z=dynamicsymbols('z')
    
    scheme_name = 'bridge_dmp.png'
    real_name = 'leaf_spring.jpg'
    
    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=None,
                 g=None,
                 l=None,
                 n=None,
                 t=None,
                 w=None,
                 module=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if l is not None: self.l=l
        if n is not None: self.n=n
        if t is not None: self.t=t
        if w is not None: self.w=w
        if z is not None: self.z = z
        if module is not None: self.module=module


        
        self.qs=[self.z]
        self._init_from_components(**kwargs)


        
#         self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)



    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z, qs=self.qs)
        self._spring = Spring(self.k_beam, self.z, qs=self.qs)
        self._gravitational_force = GravitationalForce(self.m, self.g, self.z, qs=self.qs)
        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravitational_force'] = self._gravitational_force

        return components


    def set_equivalent_data(self):
        #E - Young's modulus, n - number of leaves, w - width of leaves, t - thickness of leaves, l - span; https://www.piping-designer.com/index.php/disciplines/mechanical/stationary-equipment/fastener/2486-leaf-spring-stiffness

        equivalent_data_dict = {
            self.k_beam: (8*self.module*self.n*self.w*self.t**3)/(3*self.l**3),
        }

        return equivalent_data_dict
    
    def get_default_data(self):
        default_data_dict = {self.m:3000,self.module:2.1*10**11,self.l:0.9,self.g:9.81,self.n:8,self.t:0.0135,self.w:0.08
        }

        return default_data_dict
    def get_real_data(self):

        real_data_dict = {self.m:3000,self.module:2.1*10**11,self.l:0.9,self.g:9.81,self.n:8,self.t:0.0135,self.w:0.08
        }

        return real_data_dict
    
