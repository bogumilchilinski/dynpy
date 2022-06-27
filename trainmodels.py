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
from dynpy.models.elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from dynpy.models.sdof import ComposedSystem, SpringMassSystem, DampedSpringMassSystem, BeamBridgeDamped
ureg = pint.UnitRegistry()

mechanics_printing()



t=Symbol('t') #independent variable - time

class LEWEL16_simple(SpringMassSystem):
    scheme_name = 'train_sdof.png'
    real_name = 'lewel16.jpg'

class LEWEL16_simple_damped(DampedSpringMassSystem):
    scheme_name = 'train_sdof_damped.png'
    real_name = 'lewel16.jpg'

class LeafSpring(ComposedSystem):

    scheme_name = 'bridge_dmp.png'
    real_name = 'leaf_spring.jpg'
    def __init__(self,
                 m=Symbol('m', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 c=Symbol('c', positive=True),
                 l=Symbol('l', positive=True),
                 module=Symbol('E', positive=True),
                 inertia=Symbol('I', positive=True),
                 lam=Symbol('lambda',positive=True),
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.m = m
        self.c=c
        self.k_beam = k_beam
        self.lam=lam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0
        self.l=l
        self.z = z
        self.module=module
        self.inertia=inertia
#         c=self.lam*k_beamhaft

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(m, g, z)
#         self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)
        self.damper = Damper(c, pos1=z, qs=[z])
        composed_system = (self.mass + self.spring + self.damper)

        super().__init__(composed_system,**kwargs)


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
            self.m:b*h*self.l*rho
        }

        return default_data_dict
    
    def get_numerical_data(self):
        b,h,rho=symbols('b h rho', positive=True)
        default_data_dict = {self.module:2.1*10**11,self.rho:7800,self.b:0.08,self.h:9*0.012,self.l:1.5,self.lam:0.1
        }

        return default_data_dict