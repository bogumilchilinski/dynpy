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

    m=Symbol('m', positive=True),
    k_beam=Symbol('k_beam', positive=True)
    ivar=Symbol('t')
    g=Symbol('g', positive=True)
    Omega=Symbol('Omega', positive=True)
    F_0=Symbol('F_0', positive=True)
    c=Symbol('c', positive=True)
    l=Symbol('l', positive=True)
    module=Symbol('E', positive=True)
    inertia=Symbol('I', positive=True)
    lam=Symbol('lambda',positive=True)
    z=dynamicsymbols('z')
    
    scheme_name = 'bridge_dmp.png'
    real_name = 'leaf_spring.jpg'
    
    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=None,
                 g=None,
                 Omega=None,
                 F_0=None,
                 c=None,
                 l=None,
                 module=None,
                 inertia=None,
                 lam=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if c is not None: self.c=c
        if k_beam is not None: self.k_beam = k_beam
        if lam is not None: self.lam=lam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F_0 is not None: self.F_0 = F_0
        if l is not None: self.l=l
        if z is not None: self.z = z
        if module is not None: self.module=module
        if inertia is not None: self.inertia=inertia
#         c=self.lam*k_beamhaft
        
        self.qs=[self.z]
        self._init_from_components(**kwargs)


        
#         self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)



    @property
    def components(self):
        z=self.z
        m=self.m
        components = {}

        self._mass = MaterialPoint(m, z, z)(label = 'Material point')
        self._spring = Spring(self.k_beam, self.z, qs=[self.z])(label = 'Structural stiffness')
        self._gravitational_force = GravitationalForce(self.m, self.g, self.z)(label = 'Gravitational force')
        self._damper = Damper(self.c, pos1=self.z, qs=[self.z])(label = 'Internal Rayleigh damping')
        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravitational_force'] = self._gravitational_force
        components['_damper'] = self._damper
        return components


    def set_equivalent_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        self.b=b
        self.h=h
        self.rho=rho
        default_data_dict = {

#             self.lam:[10],
            self.c:self.k_beam*self.lam,
            self.k_beam: S.One*48*self.module * self.inertia / self.l**3,
            self.inertia: b*h**3/12,
#             self.m:b*h*self.l*rho
        }

        return default_data_dict
    
    def get_default_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        default_data_dict = {self.m:3000,self.module:2.1*10**11,self.rho:7800,self.b:0.08,self.h:9*0.012,self.l:1.5,self.lam:0.1
        }

        return default_data_dict
    def get_real_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        default_data_dict = {self.m:3000,self.module:2.1*10**11,self.rho:7800,self.b:0.08,self.h:8*0.0135,self.l:0.9,self.lam:0.1
        }

        return default_data_dict