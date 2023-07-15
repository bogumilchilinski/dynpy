from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, pi, latex, dsolve,
                   solve, fraction, factorial,Subs, Number, oo, Abs, N)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from .sdof import Pendulum, EngineVerticalSpringGravity
from ..continuous import ContinuousSystem, PlaneStressProblem

import os

import base64
import random
import IPython as IP
import numpy as np

import inspect

from functools import cached_property, lru_cache

class ComposedSystem(HarmonicOscillator):
    """Base class for all systems

    """
    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'
    detail_scheme_name = 'damped_car_new.PNG'
    detail_real_name = 'car_real.jpg'
    _default_args = ()
    _default_folder_path = "./dynpy/models/images/"
    
    @classmethod
    def _scheme(cls):

        path = cls._default_folder_path + cls.scheme_name
        
        return path

    @classmethod
    def _real_example(cls):
        path = cls._default_folder_path + cls.real_name

        return path
    
    @classmethod
    def _detail_real(cls):
        path = cls._default_folder_path + cls.detail_real_name

        return path
    
    @classmethod
    def _detail_scheme(cls):
        path = cls._default_folder_path + cls.detail_scheme_name

        return path

    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))


    def get_default_data(self):
        return None

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()


        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict=None

        return parameters_dict

    @lru_cache   
    def linearized(self, x0=None, op_point=False, hint=[], label=None):

        return type(self).from_system(super().linearized(x0=x0,op_point=op_point,hint=hint,label=label))



#DO ZROBIENIA
class Inverted_Pendulum(HarmonicOscillator):
    def __init__(self,
                 M=symbols('M', positive=True),
                 m=symbols('m', positive=True),
                 I=symbols('I', positive=True),
                 g=symbols('g', positive=True),
                 b=symbols('b', positive=True),
                 l=symbols('l', positive=True),
                 F=symbols('F', positive=True),
                 var=dynamicsymbols('x, phi'),
                 ivar=Symbol('t'),
                 **kwargs):

        x, phi = var

        self.rod = (
            RigidBody2D(m,
                        I,
                        pos_lin=0,
                        pos_rot=0,
                        pos_lin_c=(x + l * sin(phi)),
                        pos_rot_c=phi,
                        qs=[x, phi]) +
            MaterialPoint(m, l * cos(phi), qs=[phi]) +
            GravitationalForce(m, g, pos1=0, pos_c=l * cos(phi), qs=[phi]))

        self.cart = MaterialPoint(M, x, qs=[x])

        self.force = Force(F, x)

        self.friction = Damper(b, x)

        system = self.rod + self.cart + self.friction + self.force

        super().__init__(system,**kwargs)

#DO ZROBIENIA
class TwoDisksWithThreeSprings(ComposedSystem):
    scheme_name = 'ddof_disks_3_springs_scheme.png'
    real_name = 'nonlin_trolley_real.PNG'

    def __init__(self,
                 d=Symbol('d', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 kl=Symbol('k_l', positive=True),
                 kc=Symbol('k_c', positive=True),
                 kr=Symbol('k_r', positive=True),
                 R=Symbol('R', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 xl=dynamicsymbols('x_l'),
                 xr=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l, x_r'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.kl = kl
        self.kc = kc
        self.kr = kr
        self.R = R
        self.l_0 = l_0
        self.d = d
        self.xl = xl
        self.xr = xr
        self.x = x

        self.Disk1 = MaterialPoint(m1, xl, qs=[xl]) + MaterialPoint(m1/2*R**2, xl/R, qs=[xl]) + Spring(kl, xl, qs=[xl])
        self.Disk2 = MaterialPoint(m2, xr, qs=[xr]) + MaterialPoint(m2/2*R**2, xr/R, qs=[xr]) + Spring(kr,xr, qs=[xr])
        self.Spring = Spring(kc, xl, xr, qs=[xl, xr])

        system = self.Disk1 + self.Spring + self.Disk2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

#             self.d: [1 * l0, 2 * l0, S.Half * l0, 4 * S.Half * l0,  S.Half**2 * l0],

            self.kl: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
            self.kr: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
            self.kc: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }



        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Disk Mass',
            self.m2: r'Disk Mass',
            self.kl: 'Left Spring Stiffness',
            self.kr: 'Right Spring Stiffness',
            self.kc: 'Central Spring Stiffness',
            self.l: r'Length',
            self.l_0: r'initial Spring Length',
        }
        return self.sym_desc_dict


#DO ZROBIENIA
class VeeEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'vee_engine_perpendicular_springs.png'
    real_name = '440_magnum_v8.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l',positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1x = MaterialPoint(M, pos1=x, qs=[x])
        self.MaterialPoint_1z = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2x = MaterialPoint(m_e,
                                             pos1=x + e * sin(phi),
                                             qs=[x])
        self.MaterialPoint_2z = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.SpringVer = Spring(2 * k_m, pos1=z, qs=[z])
        self.SpringHor = Spring(2 * k_m, pos1=x, qs=[x])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1x + self.MaterialPoint_2x + self.MaterialPoint_1z + self.MaterialPoint_2z
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


#DO ZROBIENIA
class QuarterOfVehicle(ComposedSystem):
    
    def __init__(self,
                 k_d=None,
                 k_u=None,
                 m_d=None,
                 m_u=None,
                 c=None,
                 omega=None,
                 A_0=None,
                 x_1 =None,
                 x_2=None,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
    

        self.qs = [self.x_1,self.x_2]
    
        if m_d is not None: self.m_d = m_d
        if k_d is not None: self.k_d = k_d
        if m_u is not None: self.m_u = m_u
        if k_u is not None:self.k_u = k_u
        if c is not None: self.c=c
        if omega is not None: self.omega=omega
        if A_0 is not None: self.A_0=A_0
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None : self.x_2 = x_2
        if qs is not None: self.qs = qs


        
        self.material_point_1 = MaterialPoint(self.m_d, self.x_1, qs=self.qs)
        self.spring_1 = Spring(self.k_d,self.x_1,self.A_0*sin(ivar*self.omega),qs=self.qs)
        self.material_point_2 = MaterialPoint(self.m_u, self.x_2, qs=self.qs)
        self.spring_2 = Spring(self.k_u,self.x_1,self.x_2,qs=self.qs)
        self.damper_1=Damper(self.c,self.x_1,self.x_2, qs=self.qs)


        system =  self.material_point_1 + self.spring_1+ self.material_point_2+self.spring_2+self.damper_1
        super().__init__(system,**kwargs)



    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_d: r'unsprung mass',
            self.m_u: r'sprung mass',
            self.k_d: r'wheel spring stiffness coefficient',
            self.k_u: r'body spring stiffness coefficient',
            self.c: r'współczynnik tłumienia zawieszenia',
            self.A_0: r'amplitude',
            self.x_1: r'wheel displacement',
            self.x_2: r'body displacement',
            self.x_1.diff(t): r'wheel velocity',
            self.x_2.diff(t): r'body velocity',
            self.omega: r'circular frequency',
            self.ivar: r'time',
            
        }
        return self.sym_desc_dict 

#DO ZROBIENIA
class DDoFTwoNonLinearTrolleys(ComposedSystem):
    scheme_name = 'ddof_nonlin_trolleys.PNG'
    real_name = 'tension_leg_platform.png'

    def __init__(self,
                 g=Symbol('g', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k1=Symbol('k_1', positive=True),
                 k2=Symbol('k_2', positive=True),
                 k3=Symbol('k_3', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 x1=dynamicsymbols('x1'),
                 x2=dynamicsymbols('x2'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x1, x2'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.d = d
        self.l_0 = l_0
        self.x1 = x1
        self.x2 = x2
        self.x = x

        self.Trolley1 = MaterialPoint(m1, x1, qs=[x1]) + Spring(
            k1, pos1=(sqrt(x1**2 + d**2) - l_0), qs=[x1])
        self.Trolley2 = MaterialPoint(m2, x2, qs=[x2]) + Spring(
            k2, pos1=(sqrt(x2**2 + d**2) - l_0), qs=[x2])
        self.Spring = Spring(k3, x1, x2,qs=[x1,x2])

        system = self.Trolley1 + self.Spring + self.Trolley2
        super().__init__(system(qs=qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.d: [1 * l0, 2 * l0, S.Half * l0, 3 * S.Half * l0, 1 * l0],
            self.k1: [S.Half * k0, S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0],
            self.k2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k3: [S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0, 5 * S.Half * k0],
            self.x1: [self.x, 0],
            self.x2: [self.x, S.Zero],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x1] == S.Zero:
            parameters_dict[self.x2] = self.x

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley Mass',
            self.m2: r'Trolley Mass',
            self.k1: 'Spring Stiffness',
            self.k2: 'Spring Stiffness',
            self.k3: 'Spring Stiffness',
            self.d: r'length',
            self.l_0: r'length',
        }
        return self.sym_desc_dict

