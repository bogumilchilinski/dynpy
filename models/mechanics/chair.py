from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from ..mechanics.disk import RollingDisk


import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem,  base_frame, base_origin, REPORT_COMPONENTS_LIST
from functools import cached_property

t,f= symbols('t, f')


from sympy import * # provides mathematical interface for symbolic calculations

from pint import UnitRegistry
ureg = UnitRegistry()

# from sympy.physics.mechanics import *

# import sympy as sym
# import numpy as np
# import numpy.fft as fft

# import matplotlib.pyplot as plt
# import math

# from dynpy.utilities.report import SystemDynamicsAnalyzer

# import pandas as pd

# from dynpy.utilities.adaptable import *




# ix = pd.IndexSlice



#TODO - road profiles

#u0=Symbol('u_0',positive=True)
#l_ramp=Symbol('l_ramp',positive=True)
#t_0=Symbol('t_0',positive=True)
#delta_t=Symbol('Delta_t',positive=True)

# # u_rr=u0/2*(1+(Heaviside(x-t0-2*R,0.5) - Heaviside(x-t0-delta_t-2*R,0.5) ))#road profile given in time domain
# # u_rf=u0/2*(1+(Heaviside(x-t0,0.5) - Heaviside(x-t0-delta_t,0.5) ))#road profile given in time domain

# u_rf_bump_on=u0/2*(1+2/pi*atan(5*(t-t0)))-u0/2*(1+2/pi*atan(5*(t-t0-t_l)))
# u_rr_bump_on=u_rf_bump.subs(t,t-1) 
#  #road profile given in time domain
    
# u_rf_mul_bumps_on=sum(u_rf_bump_on.subs(x,x-ii*l_bumps) for ii in range(5))
# u_rr_mul_bumps_on=u_rf_mul_bumps_on.subs(x,x-2*R)
#  #road profile given in time domain
    
# # u_rf_ramp=u0/l_ramp*x
# # u_rr_ramp=u_rf_ramp.subs(x,x-2*R)

# u_road=Function('u_r')
# # u_rf=u_road(x)
# # u_rr=u_road(x-2*R)

# # u_rr=u0*(sin(x-2*R))  #road profile given in time domain
# # u_rf=u_rr.subs(x,x+2*R) #road profile given in time domain




###################+++++++++++++++++++++++++++++++++++++++++++ NEWEST SYSTEMS

class UndampedChair5DOF(ComposedSystem):

    scheme_name = 'chair5dof.jpg'
    real_name = 'gtm.png'

    x=dynamicsymbols('x')
    z=dynamicsymbols('z')
    phi=dynamicsymbols('varphi')
    z_rear=dynamicsymbols('z_r')
    z_front=dynamicsymbols('z_f')

    M=Symbol('M', positive=True)
    I_ch=Symbol('I_chair',positive=True)
    m_rear=Symbol('m_rw',positive=True)
    I_w=Symbol('I_rw', positive=True)
    R=Symbol('R',positive=True)
    m_fr=Symbol('m_fw',positive=True)
    g=Symbol('g',positive=True)
    z_c=Symbol('z_c',positive=True)
    k_rt=Symbol('k_rt',positive=True)
    u_rrp=Symbol('u_rrp',positive=True)
    k_r=Symbol('k_r',positive=True)
    u_rsd=Symbol('u_rsd',positive=True)
    k_ft=Symbol('k_ft',positive=True)
    u_frp=Symbol('u_frp',positive=True)
    k_f=Symbol('k_f',positive=True)
    u_fsd=Symbol('u_fsd',positive=True)
    F=Symbol('F',positive=True)
    Omega=Symbol('Omega',positive=True)
    pm=Symbol('PM',positive=True)
    
    ivar=Symbol('t')

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components={}

        self._body_horizontal = MaterialPoint(self.M, self.x+self.z_c*self.phi, qs=self.qs,frame=base_frame)(label='body_horizontal')
        self._body_vertical_rotational = RigidBody2D(self.M, self.I_ch, pos_lin=self.z, pos_rot=self.phi, qs=self.qs,frame=base_frame)(label='body_vertical_rotational')
        self._rear_wheel_vertical_rotational = RigidBody2D(self.m_rear, self.I_w, pos_lin=self.z_rear, pos_rot=self.x/self.R, qs=self.qs,frame=base_frame)(label='rear_wheel_vertical_rotational')
        self._rear_wheel_horizontal = MaterialPoint(self.m_rear, self.x, qs=self.qs,frame=base_frame)(label='rear_wheel_horizontal')
        self._front_wheel_vertical = MaterialPoint(self.m_fr, self.z_front, qs=self.qs,frame=base_frame)(label='front_wheel_vertical')
        self._front_wheel_horizontal = MaterialPoint(self.m_fr, self.x, qs=self.qs,frame=base_frame)(label='front_wheel_horizontal')
        self._front_wheel_gravity = GravitationalForce(self.m_fr, self.g, self.z_front)(label='front_wheel_gravity')
        self._rear_wheel_gravity = GravitationalForce(self.m_rear, self.g, self.z_rear)(label='rear_wheel_gravity')
        self._body_gravity = GravitationalForce(self.M, self.g, self.z+self.z_c*cos(self.phi), qs=self.qs)(label='body_gravity')
        self._rear_tire_spring = Spring(self.k_rt, self.z_rear-self.u_rrp, qs=self.qs,frame=base_frame)('rear_tire_spring')
        self._rear_wheel_spring = Spring(self.k_r, self.u_rsd, qs=self.qs,frame=base_frame)('rear_wheel_spring')
        self._front_tire_spring = Spring(self.k_ft, self.z_front-self.u_frp, qs=self.qs,frame=base_frame)('front_tire_spring')
        self._front_wheel_spring = Spring(self.k_f, self.u_fsd, qs=self.qs,frame=base_frame)('front_wheel_spring')
        self._hand_drive_force = Force(self.F*(sign(cos(self.Omega*self.ivar)-self.pm)+1), pos1=self.x , qs=self.qs,frame=base_frame)('hand_drive_force')
        self._hand_drive_force_phi = Force(self.F*self.Omega*self.R*((sign(cos(self.Omega*self.ivar)-self.pm)+1)*cos(self.Omega*self.ivar)), pos1=self.phi , qs=self.qs,frame=base_frame)('hand_drive_force_phi')


        components['_body_horizontal'] = self._body_horizontal
        components['_body_vertical_rotational'] = self._body_vertical_rotational
        components['_rear_wheel_vertical_rotational'] = self._rear_wheel_vertical_rotational
        components['_rear_wheel_horizontal'] = self._rear_wheel_horizontal
        components['_front_wheel_vertical'] = self._front_wheel_vertical
        components['_front_wheel_horizontal'] = self._front_wheel_horizontal
        components['_front_wheel_gravity'] = self._front_wheel_gravity
        components['_rear_wheel_gravity'] = self._rear_wheel_gravity
        components['_body_gravity'] = self._body_gravity
        components['_rear_tire_spring'] = self._rear_tire_spring
        components['_rear_wheel_spring'] = self._rear_wheel_spring
        components['_front_tire_spring'] = self._front_tire_spring
        components['_front_wheel_spring'] = self._front_wheel_spring
#         components['_hand_drive_force'] = self._hand_drive_force
        components['_hand_drive_force_phi'] = self._hand_drive_force_phi
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
                                self.x:'main body center of mass horizontal displacement',
                                self.z:'main body center of mass vertical displacement',
                                self.phi:'angular displacement of the main mass',
                                self.z_rear:'rear wheel vertical displacement',
                                self.M:'''mass of the main body including user''',
                                self.I_ch:'''wheelchair's moment of inertia''',
                                self.m_rear:'rear wheel mass',
                                self.R:'rear wheel radius',
                                self.I_w:'''rear wheel moment of inertia''',
                                self.m_fr:'front wheel mass',
                                self.g:'gravitational acceleration',
                                self.z_c:'''vertical distance between wheelchair's center of gravity and center of rotation''',
                                self.k_rt:'rear wheel tire stiffness',
                                self.u_rrp:'rear road profile',
                                self.k_r:'rear wheel stiffness',
                                self.k_ft:'front wheel tire stiffness',
                                self.u_frp:'front road profile',
                                self.k_r:'front wheel stiffness',
                                self.F:'disabled driver arm force',
                                self.Omega:'driving force frequency',
                                self.pm:'force duty cycle',
        }

        return self.sym_desc_dict

    def get_param_values(self):

        self.default_data_dict={
                            self.M:75,
                            self.I_ch:9.479342+self.M*0.39*0.39,
                            self.m_rear:1.5,
                            self.R:0.3,
                            self.I_w:self.m_rear*self.R**2,
                            self.m_fr:0.6,
                            self.g:9.81,
                            self.z_c:0.4,
                            self.k_rt:109000,
                            self.k_r:750000,
                            self.k_ft:300000,
                            self.k_f:750000,
                            self.F:150,
                            self.Omega:0.3*np.pi*2,
                            self.pm:0.1,
                          }

        return self.default_data_dict
    
    def get_table_values(self):
        table_data_dict={self.F:150,
                   
                   self.c_mu:0.0001,
                   self.c_lam:0.0001,
                   self.l_l:0.2,
                   self.l_r:0.4,
                   
                   self.k_f:607500,
                   self.k_ft:475000,
                   self.k_r:580000,
                   self.k_rt:400000,
                   self.m_3:75,
                   self.I_ch:20.8868,
                   self.m_rear:1.5,
                   self.m_fr:0.6,
                   self.pm:0.1,
                   self.Omega:0.3454,
                   self.R:0.3,
                   self.z_c3:0.4,
                   self.g:9.81,
                   self.I_w:0.135,
                   self.l_fr:0.2,
                   self.l_rear:0.01,
                   self.u0:0.005,

                   self.l_bumps:0.15,
                   self.amplitude:0.0165,
                   self.length:0.19,
                   self.speed:1.7,
                   self.axw:0.47}
        return table_data_dict  
    
    
class DampedChair5DOF(UndampedChair5DOF):

    c_mu=Symbol('c_mu',positive=True)
    c_lam=Symbol('c_lambda',positive=True)

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                c_mu=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z-rear
        if z_front is not None: self.z_front=_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if c_mu is not None: self.c_mu=c_mu
        if c_lam is not None: self.c_lam=c_lam
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components={}

#         self._undamped_chair_5DOF = UndampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,self.u_rrp,self.k_r,self.u_rsd,self.k_ft,self.u_frp,self.k_f,self.ivar,self.u_fsd,self.F,self.Omega,self.pm)('undamped_chair_5DOF')
        
        self._undamped_chair_5DOF = UndampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,0,self.k_r,((self.z+self.R*self.phi+Symbol('l_r',positive=True) -self.z_rear)**2 )-Symbol('l_r',positive=True),self.k_ft,0,self.k_f,self.ivar,((self.z-self.R*self.phi+Symbol('l_fr',positive=True)   -self.z_front)**2  )-Symbol('l_fr',positive=True),self.F,self.Omega,self.pm)('undamped_chair_5DOF')
        self._rear_wheel_damping = Damper(self.I_w*self.c_mu,self.x/self.R,qs=self.qs)('rear_wheel_damping')
#         self._rayleigh_damping_inertia_x = Damper(self.rayleigh_damping_inertia_x()*self.c_mu,self.x,qs=self.qs)
        self._damping_body_horizontal = Damper(self.M*self.c_mu, self.x+self.z_c*self.phi, qs=self.qs)(label='_damping_body_horizontal')
        self._damping_body_vertical = Damper(self.M*self.c_mu, self.z, qs=self.qs)(label='damping_body_vertical')
        self._damping_body_rotational = Damper(self.I_ch*self.c_mu, self.phi, qs=self.qs)(label='damping_body_rotational')
        self._damping_rear_wheel_vertical = Damper(self.m_rear*self.c_mu, self.z_rear, qs=self.qs)(label='damping_rear_wheel_vertical')
        self._damping_rear_wheel_rotational = Damper(self.I_w*self.c_mu, self.x/self.R, qs=self.qs)(label='damping_rear_wheel_rotational')
        self._damping_rear_wheel_horizontal = Damper(self.m_rear*self.c_mu, self.x, qs=self.qs)(label='damping_rear_wheel_horizontal')
        self._damping_front_wheel_vertical = Damper(self.m_fr*self.c_mu, self.z_front, qs=self.qs)(label='damping_front_wheel_vertical')
        self._damping_front_wheel_horizontal = Damper(self.m_fr*self.c_mu, self.x, qs=self.qs)(label='damping_front_wheel_horizontal')
            
        components['_undamped_chair_5DOF'] = self._undamped_chair_5DOF
        components['_damping_body_horizontal'] = self._damping_body_horizontal
        components['_damping_body_vertical'] = self._damping_body_vertical
        components['_damping_body_rotational'] = self._damping_body_rotational
        components['_damping_rear_wheel_vertical'] = self._damping_rear_wheel_vertical
        components['_damping_rear_wheel_rotational'] = self._damping_rear_wheel_rotational
        components['_damping_rear_wheel_horizontal'] = self._damping_rear_wheel_horizontal
        components['_damping_front_wheel_vertical'] = self._damping_front_wheel_vertical
        components['_damping_front_wheel_horizontal'] = self._damping_front_wheel_horizontal


        return components
    
    def rayleigh_damping_inertia_x(self):
        
        return (UndampedChair5DOF().mass_matrix*Matrix(self.qs))#[0]#*self.c_mu#[0].coeff(self.c_mu)


    def get_param_values(self):

        default_data_dict={
                            self.c_mu:0.0001,
                          }

        return {**super().get_param_values(),**self.default_data_dict}
    
    def symbols_description(self):
        self.sym_desc_dict = {
                                self.c_mu:"Rayileigh's damping mass matrix coefficient"
        }

        return {**super().symbols_description(),**self.sym_desc_dict}


class BasicUndampedChair5DOF(ComposedSystem):

    scheme_name = 'chair5dof.jpg'
    real_name = 'gtm.png'

    x=dynamicsymbols('x')
    z=dynamicsymbols('z')
    phi=dynamicsymbols('varphi')
    z_rear=dynamicsymbols('z_r')
    z_front=dynamicsymbols('z_f')

    M=Symbol('M', positive=True)
    I_ch=Symbol('I_chair',positive=True)
    m_rear=Symbol('m_rw',positive=True)
    I_w=Symbol('I_rw', positive=True)
    R=Symbol('R',positive=True)
    m_fr=Symbol('m_fw',positive=True)
    g=Symbol('g',positive=True)
    z_c=Symbol('z_c',positive=True)
    k_rt=Symbol('k_rt',positive=True)
    u_rrp=Symbol('u_rrp',positive=True)
    k_r=Symbol('k_r',positive=True)
    u_rsd=Symbol('u_rsd',positive=True)
    k_ft=Symbol('k_ft',positive=True)
    u_frp=Symbol('u_frp',positive=True)
    k_f=Symbol('k_f',positive=True)
    u_fsd=Symbol('u_fsd',positive=True)
    F=Symbol('F',positive=True)
    Omega=Symbol('Omega',positive=True)
    pm=Symbol('PM',positive=True)
    
    ivar=Symbol('t')

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components={}

        self._body_horizontal = MaterialPoint(self.M, self.x+self.z_c*self.phi, qs=self.qs)(label='body_horizontal')
        self._body_vertical_rotational = RigidBody2D(self.M, self.I_ch, pos_lin=self.z, pos_rot=self.phi, qs=self.qs)(label='body_vertical_rotational')
        self._rear_wheel_vertical_rotational = RigidBody2D(self.m_rear, self.I_w, pos_lin=self.z_rear, pos_rot=self.x/self.R, qs=self.qs)(label='rear_wheel_vertical_rotational')
        self._rear_wheel_horizontal = MaterialPoint(self.m_rear, self.x, qs=self.qs)(label='rear_wheel_horizontal')
        self._front_wheel_vertical = MaterialPoint(self.m_fr, self.z_front, qs=self.qs)(label='front_wheel_vertical')
        self._front_wheel_horizontal = MaterialPoint(self.m_fr, self.x, qs=self.qs)(label='front_wheel_horizontal')
        self._front_wheel_gravity = GravitationalForce(self.m_fr, self.g, self.z_front)(label='front_wheel_gravity')
        self._rear_wheel_gravity = GravitationalForce(self.m_rear, self.g, self.z_rear)(label='rear_wheel_gravity')
        self._body_gravity = GravitationalForce(self.M, self.g, self.z+self.z_c*cos(self.phi), qs=self.qs)(label='body_gravity')
        self._rear_tire_spring = Spring(self.k_rt, self.z_rear-self.u_rrp, qs=self.qs)('rear_tire_spring')
        self._rear_wheel_spring = Spring(self.k_r, self.u_rsd, qs=self.qs)('rear_wheel_spring')
        self._front_tire_spring = Spring(self.k_ft, self.z_front-self.u_frp, qs=self.qs)('front_tire_spring')
        self._front_wheel_spring = Spring(self.k_f, self.u_fsd, qs=self.qs)('front_wheel_spring')
        self._hand_drive_force = Force(self.F*(sign(cos(self.Omega*self.ivar)-self.pm)+1), pos1=self.x , qs=self.qs)('hand_drive_force')
        self._hand_drive_force_phi = Force(self.F*self.Omega*self.R*((sign(cos(self.Omega*self.ivar)-self.pm)+1)*cos(self.Omega*self.ivar)), pos1=self.phi , qs=self.qs)('hand_drive_force_phi')

        components['_body_horizontal'] = self._body_horizontal
        components['_body_vertical_rotational'] = self._body_vertical_rotational
        components['_rear_wheel_vertical_rotational'] = self._rear_wheel_vertical_rotational
        components['_rear_wheel_horizontal'] = self._rear_wheel_horizontal
        components['_front_wheel_vertical'] = self._front_wheel_vertical
        components['_front_wheel_horizontal'] = self._front_wheel_horizontal
        components['_front_wheel_gravity'] = self._front_wheel_gravity
        components['_rear_wheel_gravity'] = self._rear_wheel_gravity
        components['_body_gravity'] = self._body_gravity
        components['_rear_tire_spring'] = self._rear_tire_spring
        components['_rear_wheel_spring'] = self._rear_wheel_spring
        components['_front_tire_spring'] = self._front_tire_spring
        components['_front_wheel_spring'] = self._front_wheel_spring
        components['_hand_drive_force'] = self._hand_drive_force
        components['_hand_drive_force_phi'] = self._hand_drive_force_phi
        
        return components

    def get_param_values(self):

        default_data_dict={
                            self.M:75,
                            self.I_ch:9.479342+self.M*0.39*0.39,
                            self.m_rear:1.5,
                            self.R:0.3,
                            self.I_w:self.m_rear*self.R**2,
                            self.m_fr:0.6,
                            self.g:9.81,
                            self.z_c:0.4,
                            self.k_rt:109000,
                            self.k_r:750000,
                            self.k_ft:300000,
                            self.k_f:750000,
                            self.F:150,
                            self.Omega:0.3*np.pi*2,
                            self.pm:0.1,
                          }

        return default_data_dict
    def symbols_description(self):
        self.sym_desc_dict = {
                                self.x:'main body center of mass horizontal displacement',
                                self.z:'main body center of mass vertical displacement',
                                self.phi:'angular displacement of the main mass',
                                self.z_rear:'rear wheel vertical displacement',
                                self.M:'''mass of the main body including user''',
                                self.I_ch:'''wheelchair's moment of inertia''',
                                self.m_rear:'rear wheel mass',
                                self.R:'rear wheel radius',
                                self.I_w:'''rear wheel moment of inertia''',
                                self.m_fr:'front wheel mass',
                                self.g:'gravitational acceleration',
                                self.z_c:'''vertical distance between wheelchair's center of gravity and center of rotation''',
                                self.k_rt:'rear wheel tire stiffness',
                                self.u_rrp:'rear road profile',
                                self.k_r:'rear wheel stiffness',
                                self.k_ft:'front wheel tire stiffness',
                                self.u_frp:'front road profile',
                                self.k_r:'front wheel stiffness',
                                self.F:'disabled driver arm force',
                                self.Omega:'driving force frequency',
                                self.pm:'force duty cycle',
        }

        return self.sym_desc_dict
    
class BasicDampedChair5DOF(BasicUndampedChair5DOF):

    c_mu=Symbol('c_mu',positive=True)

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                c_mu=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if c_mu is not None: self.c_mu=c_mu
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components={}

        self._undamped_chair_5DOF = BasicUndampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,self.u_rrp,self.k_r,self.u_rsd,self.k_ft,self.u_frp,self.k_f,self.ivar,self.u_fsd,self.F,self.Omega,self.pm)('undamped_chair_5DOF')

        self._damping_body_horizontal = Damper(self.M*self.c_mu, self.x+self.z_c*self.phi, qs=self.qs)(label='_damping_body_horizontal')
        self._damping_body_vertical = Damper(self.M*self.c_mu, self.z, qs=self.qs)(label='damping_body_vertical')
        self._damping_body_rotational = Damper(self.I_ch*self.c_mu, self.phi, qs=self.qs)(label='damping_body_rotational')
        self._damping_rear_wheel_vertical = Damper(self.m_rear*self.c_mu, self.z_rear, qs=self.qs)(label='damping_rear_wheel_vertical')
        self._damping_rear_wheel_rotational = Damper(self.I_w*self.c_mu, self.x/self.R, qs=self.qs)(label='damping_rear_wheel_rotational')
        self._damping_rear_wheel_horizontal = Damper(self.m_rear*self.c_mu, self.x, qs=self.qs)(label='damping_rear_wheel_horizontal')
        self._damping_front_wheel_vertical = Damper(self.m_fr*self.c_mu, self.z_front, qs=self.qs)(label='damping_front_wheel_vertical')
        self._damping_front_wheel_horizontal = Damper(self.m_fr*self.c_mu, self.x, qs=self.qs)(label='damping_front_wheel_horizontal')

        components['_undamped_chair_5DOF'] = self._undamped_chair_5DOF
        components['_damping_body_horizontal'] = self._damping_body_horizontal
        components['_damping_body_vertical'] = self._damping_body_vertical
        components['_damping_body_rotational'] = self._damping_body_rotational
        components['_damping_rear_wheel_vertical'] = self._damping_rear_wheel_vertical
        components['_damping_rear_wheel_rotational'] = self._damping_rear_wheel_rotational
        components['_damping_rear_wheel_horizontal'] = self._damping_rear_wheel_horizontal
        components['_damping_front_wheel_vertical'] = self._damping_front_wheel_vertical
        components['_damping_front_wheel_horizontal'] = self._damping_front_wheel_horizontal


        return components

    def get_param_values(self):

        self.default_data_dict={
                            self.c_mu:0.0001,
                          }

        return {**self.default_data_dict,**super().get_param_values()}
    
    def symbols_description(self):
        self.sym_desc_dict = {
                                self.c_mu:"Rayileigh's damping mass matrix coefficient"
        }

        return {**self.sym_desc_dict,**super().symbols_description(),}
    
class LinearSpringDeflectionUndampedChair5DOF(BasicUndampedChair5DOF):
    
    l_rear=Symbol('l_r',positive=True)
    l_front=Symbol('l_f',positive=True)

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                l_rear=None,
                l_front=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if l_rear is not None: self.l_rear=l_rear
        if l_front is not None: self.l_front=l_front
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):
        rear_spring_deflection=((self.z+self.R*self.phi+self.l_rear-self.z_rear)**2 )-self.l_rear
        front_spring_deflection=((self.z-self.R*self.phi+self.l_front-self.z_front)**2  )-self.l_front
        rear_road_profile=0
        front_road_profile=0
        components={}

        self._basic_undamped_chair_5DOF = BasicUndampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,rear_road_profile,self.k_r,rear_spring_deflection,self.k_ft,front_road_profile,self.k_f,self.ivar,front_spring_deflection,self.F,self.Omega,self.pm)('basic_undamped_chair_5DOF')
        
        components['_basic_undamped_chair_5DOF'] = self._basic_undamped_chair_5DOF
        
        return components

    def get_param_values(self):

        self.default_data_dict={
                            self.l_rear:0.1,
                            self.l_front:0.1,
                          }

        return {**self.default_data_dict,**super().get_param_values()}
    
    def symbols_description(self):
        self.sym_desc_dict = {
                                self.l_rear:'???',
                                self.l_front:'???'
        }

        return {**self.sym_desc_dict,**super().symbols_description()}
    
class LinearSpringDeflectionDampedChair5DOF(LinearSpringDeflectionUndampedChair5DOF,BasicDampedChair5DOF):

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                l_rear=None,
                l_front=None,
                c_mu=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if l_rear is not None: self.l_rear=l_rear
        if l_front is not None: self.l_front=l_front
        if c_mu is not None: self.c_mu=c_mu
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):
        rear_spring_deflection=((self.z+self.R*self.phi+self.l_rear-self.z_rear)**2 )-self.l_rear
        front_spring_deflection=((self.z-self.R*self.phi+self.l_front-self.z_front)**2  )-self.l_front
        rear_road_profile=0
        front_road_profile=0
        components={}

        self._basic_damped_chair_5DOF = BasicDampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,rear_road_profile,self.k_r,rear_spring_deflection,self.k_ft,front_road_profile,self.k_f,self.ivar,front_spring_deflection,self.F,self.Omega,self.pm,self.c_mu)('basic_damped_chair_5DOF')
        
        components['_basic_damped_chair_5DOF'] = self._basic_damped_chair_5DOF
        
        return components

class NonLinearSpringDeflectionUndampedChair5DOF(LinearSpringDeflectionUndampedChair5DOF):

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                l_rear=None,
                l_front=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if l_rear is not None: self.l_rear=l_rear
        if l_front is not None: self.l_front=l_front
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):
        rear_spring_deflection=sqrt((self.z+self.R*sin(self.phi)+self.l_rear-self.z_rear)**2+(self.R-self.R*cos(self.phi))**2)-self.l_rear
        front_spring_deflection=sqrt((self.z-self.R*sin(self.phi)+self.l_front-self.z_front)**2+(self.R-self.R*cos(self.phi))**2)-self.l_front
        rear_road_profile=0
        front_road_profile=0
        components={}

        self._basic_undamped_chair_5DOF = BasicUndampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,rear_road_profile,self.k_r,rear_spring_deflection,self.k_ft,front_road_profile,self.k_f,self.ivar,front_spring_deflection,self.F,self.Omega,self.pm)('basic_undamped_chair_5DOF')
        
        components['_basic_undamped_chair_5DOF'] = self._basic_undamped_chair_5DOF
        
        return components
    
class NonLinearSpringDeflectionDampedChair5DOF(LinearSpringDeflectionDampedChair5DOF):

    def __init__(self,
                x=None,
                z=None,
                phi=None,
                z_rear=None,
                z_front=None,
                M=None,
                I_ch=None,
                m_rear=None,
                I_w=None,
                R=None,
                m_fr=None,
                g=None,
                z_c=None,
                k_rt=None,
                u_rrp=None,
                k_r=None,
                u_rsd=None,
                k_ft=None,
                u_frp=None,
                k_f=None,
                ivar=None,
                u_fsd=None,
                F=None,
                Omega=None,
                pm=None,
                l_rear=None,
                l_front=None,
                c_mu=None,
                **kwargs):

        if x is not None: self.x=x
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        if z_rear is not None: self.z_rear=z_rear
        if z_front is not None: self.z_front=z_front

        if M is not None: self.M=M
        if I_ch is not None: self.I_ch=I_ch
        if m_rear is not None: self.m_rear=m_rear
        if I_w is not None: self.I_w=I_w
        if R is not None: self.R=R
        if m_fr is not None: self.m_fr=m_fr
        if g is not None: self.g=g
        if z_c is not None: self.z_c=z_c
        if k_rt is not None: self.k_rt=k_rt
        if u_rrp is not None: self.u_rrp=u_rrp
        if k_r is not None: self.k_r=k_r
        if u_rsd is not None: self.u_rsd=u_rsd
        if k_ft is not None: self.k_ft=k_ft
        if u_frp is not None: self.u_frp=u_frp
        if k_f is not None: self.k_f=k_f
        if u_fsd is not None: self.u_fsd=u_fsd
        if F is not None: self.F=F
        if Omega is not None: self.Omega=Omega
        if pm is not None: self.pm=pm
        if l_rear is not None: self.l_rear=l_rear
        if l_front is not None: self.l_front=l_front
        if c_mu is not None: self.c_mu=c_mu
        
        if ivar is not None: self.ivar=ivar

        self.qs = [self.x,self.z,self.phi,self.z_rear,self.z_front]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):
        rear_spring_deflection=sqrt((self.z+self.R*sin(self.phi)+self.l_rear-self.z_rear)**2+(self.R-self.R*cos(self.phi))**2)-self.l_rear
        front_spring_deflection=sqrt((self.z-self.R*sin(self.phi)+self.l_front-self.z_front)**2+(self.R-self.R*cos(self.phi))**2)-self.l_front
        rear_road_profile=0
        front_road_profile=0
        components={}

        self._basic_damped_chair_5DOF = BasicDampedChair5DOF(self.x,self.z,self.phi,self.z_rear,self.z_front,self.M,self.I_ch,self.m_rear,self.I_w,self.R,self.m_fr,self.g,self.z_c,self.k_rt,rear_road_profile,self.k_r,rear_spring_deflection,self.k_ft,front_road_profile,self.k_f,self.ivar,front_spring_deflection,self.F,self.Omega,self.pm,self.c_mu)('basic_damped_chair_5DOF')
        
        components['_basic_damped_chair_5DOF'] = self._basic_damped_chair_5DOF
        
        return components

###################+++++++++++++++++++++++++++++++++++++++++++ OLDER SYSTEMS
class DampedChair4DOF(ComposedSystem):

    scheme_name = 'chair5dof.jpg'
    real_name = 'gtm.png'
    
    z,phi,z_fr,z_rear,x=dynamicsymbols('z varphi z_f z_r x')
    M=Symbol('M', positive=True)
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_k=Symbol('m_k', positive=True)
    m_fr=Symbol('m_fw', positive=True)
    m_rear=Symbol('m_rw', positive=True)
    I_ch=Symbol('I_ch', positive=True)
    I_w=Symbol('I_rw', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_k=Symbol('l_k', positive=True)
    l_r=Symbol('l_r', positive=True)
    l_rear=Symbol('l_r', positive=True)
    l_bumps=Symbol('l_bumps', positive=True)
    l_fr=Symbol('l_fr', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_r=Symbol('k_r', positive=True)
    k_rt=Symbol('k_rt', positive=True)
    k_f=Symbol('k_f', positive=True)
    k_ft=Symbol('k_ft', positive=True)
    c_rs=Symbol('c_rs', positive=True)
    c_fs=Symbol('c_fs', positive=True)
    c_rt=Symbol('c_rt', positive=True)
    c_ft=Symbol('c_ft', positive=True)
    c=Symbol('c', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    F=Symbol('F_0', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    Omega=Symbol('Omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    ivar=Symbol('t')
    f=Symbol('f', positive=True)
    k=Symbol('k', positive=True)
    g=Symbol('g', positive=True)
    R=Symbol('R_rw', positive=True)
    amplitude=Symbol('amplitude',positive=True)
    length=Symbol('length',positive=True)
    speed=Symbol('speed',positive=True)
    Leng=Symbol('Leng',positive=True)
    k_rot=Symbol('k_rot',positive=True)
    F_1=Symbol('F_1',positive=True)
    F_2=Symbol('F_2',positive=True)
    v0=Symbol('v_0',positive=True)
    u0=Symbol('u_0',positive=True)
    z_c3=Symbol('z_c3',positive=True)
    xm_3=Symbol('xm_3',positive=True)
    t_l=Symbol('t_l',positive=True)
    delta_t=Symbol('delta_t',positive=True)
    l_ramp=Symbol('l_ramp',positive=True)
    axw=Symbol('axw',positive=True)
    pm=Symbol('PM',positive=True)
    M_m=Symbol('M_m',positive=True)
    K_m=Symbol('K_m',positive=True)
    C_m=Symbol('C_m',positive=True)
    T=Symbol('T')
    V=Symbol('V')
    D=Symbol('D')
#     xt=dynamicsymbols('qwe') #dummy variable for deleting t symbol from descriptions
    def __init__(self,
                 m=None,
                 M=None,
                 m_k=None,
                 m_fr=None,
                 m_rear=None,
                 m_b=None,
                 I_ch=None,
                 I_w=None,
                 l_rod=None,
                 l_l=None,
                 l_r=None,
                 l_fr=None,
                 l_b=None,
                 l_k=None,
                 k_r=None,
                 k_rt=None,
                 k_f=None,
                 k_ft=None,
                 c_l=None,
                 c_p=None,
                 c_rt=None,
                 c_ft=None,
                 c=None,
                 f=None,
                 g=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 x=None,
                 s=None,
                 phi=None,
                 z_fr=None,
                 z_rear=None,
                 F=None,
                 A=None,
                 R=None,
                 Omega=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
                 l_bumps=None,
                 amplitude=None,
                 length=None,
                 speed=None,
                 Leng=None,
                 k=None,
                 k_rot=None,
                 F_1=None,
                 F_2=None,
                 v0=None,
                 u0=None,
                 z_c3=None,
                 xm_3=None,
                 t_l=None,
                 delta_t=None,
                 l_ramp=None,
                 axw=None,
                 pm=None,
#                  xt=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if x is not None: self.x=x
        if z_rear is not None: self.z_rear=z_rear
        if z_fr is not None: self.z_pw=z_pw
        if z_l is not None: self.z_lw=z_lw
        if z_p is not None: self.z_pw=z_pw
        if z_b is not None: self.z_b=z_b
        if phi is not None: self.phi=phi

        if M is not None: self.M = M  # mass of a rod
        if m is not None: self.m = m
        if m_fr is not None: self.m_fr = m_fr
        if m_rear is not None: self.m_rear = m_rear
        if m_b is not None: self.m_b = m_b
        if m_k is not None: self.m_k = m_k
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_fr is not None: self.l_fr = l_fr
        if l_r is not None: self.l_r = l_r #offset of right spring
        if l_b is not None: self.l_b = l_b  
        if l_k is not None: self.l_k = l_k  
        if l_rod is not None: self.l_rod = l_rod
        if k_r is not None: self.k_r = k_r
        if k_rt is not None: self.k_rt = k_rt
        if k_f is not None: self.k_f = k_f
        if k_ft is not None: self.k_ft = k_ft
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_rt is not None: self.c_rt=c_rt
        if c_ft is not None: self.c_ft=c_ft
        if c is not None: self.c=c
        if F is not None: self.F=F
        if A is not None: self.A=A
        if R is not None: self.R=R
        if f is not None: self.f=f
        if g is not None: self.g=g
        if Omega is not None: self.Omega=Omega
        if omega is not None: self.omega=omega
        if I_ch is not None: self.I_ch = I_ch  # moment of inertia of a rod
        if I_w is not None: self.I_w = I_w
        if F_engine is not None: self.F_engine = F_engine
        if l_bumps is not None: self.l_bumps = l_bumps
        if amplitude is not None: self.amplitude = amplitude
        if length is not None: self.length = length
        if speed is not None: self.speed = speed
        if Leng is not None: self.Leng = Leng
        if k is not None: self.k=k
        if k_rot is not None: self.k_rot=k_rot
        if F_1 is not None: self.F_1=F_1
        if F_2 is not None: self.F_2=F_2
        if v0 is not None: self.v0=v0
        if u0 is not None: self.u0=u0
        if z_c3 is not None: self.z_c3=z_c3
        if xm_3 is not None: self.xm_3=xm_3
        if t_l is not None: self.t_l=t_l
        if delta_t is not None: self.delta_t=delta_t
        if l_ramp is not None: self.l_ramp=l_ramp
        if axw is not None: self.axw=axw
        if pm is not None: self.pm=pm
#         if xt is not None: self.xt=xt
        self.s=self.A*sin(ivar*self.omega)
        self.qs=[self.z_fr,self.z_rear,self.x,self.z,self.phi]
       
        
        self.right_wheel = MaterialPoint(self.m_fr, pos1=self.z_fr, qs=[self.z_fr])
        self.left_wheel = MaterialPoint(self.m_rear, pos1=self.z_rear, qs=[self.z_rear])
        self.spring_pw = Spring(self.k_ft, pos1=self.z_fr, pos2=self.s, qs=[self.z_fr])
        self.spring_lw = Spring(self.k_rt, pos1=self.z_rear, pos2=self.s, qs=[self.z_rear])
        self.damper_pw = Damper(self.c_ft, pos1=self.z_fr, pos2=self.s, qs=[self.z_fr])
        self.damper_lw = Damper(self.c_rt, pos1=self.z_rear, pos2=self.s, qs=[self.z_rear])
        
   
        
#         self.body = RigidBody2D(self.m, self.I, pos_lin=self.z+self.x, pos_rot=self.phi, qs=[self.z, self.phi,self.x])  # rod
        self.body_hor= MaterialPoint(self.M, pos1=self.x, qs=[self.x])
        self.body_ver= MaterialPoint(self.M, pos1=self.z, qs=[self.z])
        self.body_rot= MaterialPoint(self.I_ch, pos1=self.phi, qs=[self.phi])
        self.front_wheel_hor=MaterialPoint(self.m_fr,pos1=self.x,qs=[self.x])
        self.rear_wheel_hor=MaterialPoint(self.m_rear,pos1=self.x,qs=[self.x])
        self.front_wheel_rot=MaterialPoint(self.I_w,pos1=self.x/self.R,qs=[self.x])
#         self.pasazer= MaterialPoint(self.m_k, pos1=self.z+self.phi*self.l_k, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_r, pos1=self.z+self.phi*self.l_l , pos2 = self.z_rear , qs=[self.z, self.phi, self.z_rear])  # left spring
        self.spring_2 = Spring(self.k_f, pos1=self.z-self.phi*self.l_r , pos2 = self.z_fr , qs=[self.z, self.phi, self.z_fr])
        self.damper_1= Damper(self.c_rs, pos1=self.z+self.phi*self.l_l , pos2 = self.z_rear , qs=[self.z, self.phi, self.z_rear])
        self.damper_2 = Damper(self.c_fs, pos1=self.z-self.phi*self.l_r , pos2 = self.z_fr , qs=[self.z, self.phi, self.z_fr])
        self.damper_x=Damper(self.c,pos1=self.x,qs=[self.x])
        self.damper_z=Damper(self.c,pos1=self.z,qs=[self.z])
        self.damper_phi=Damper(self.c,pos1=self.phi,qs=[self.phi])
        self.force_1 = Force(self.F*(sign(cos(self.Omega*self.ivar)-self.pm)+1), pos1=self.x , qs=[self.x])
        self.gravitational_force = GravitationalForce(self.M, self.g, (self.z+self.z_c3)*cos(self.phi), qs = self.qs)
#         self.force_2 = Force(self.A*sin(ivar*self.omega), pos1=self.z_pw, qs=[self.z_pw])
        
        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body_hor + self.body_ver + self.body_rot  +self.spring_1 + self.spring_2 +self.damper_1+self.damper_2 + self.right_wheel + self.left_wheel + self.spring_pw+ self.spring_lw+self.damper_lw+self.damper_pw+self.damper_x +self.force_1+self.front_wheel_hor+self.front_wheel_rot + self.rear_wheel_hor + self.gravitational_force + self.damper_z + self.damper_phi

        super().__init__(**{'system':system,**kwargs})
        
    def get_param_values(self):
        default_data_dict={self.F:150,
                   
#                    c_mu:0.0001,
#                    c_lam:0.0001,
                   self.l_l:0.2,
                   self.l_r:0.4,
                   
                   self.k_f:750000,
                   self.k_ft:300000,
                   self.k_r:750000,
                   self.k_rt:109000,
                   self.M:75,
                   self.I_ch:9.479342+self.M*0.39*0.39,
                   self.m_rear:1.5,
                   self.m_fr:0.6,
                   self.pm:0.1,
                   self.Omega:0.3*np.pi*2,
                   self.R:0.3,

                   self.g:9.81,
                   self.I_w:self.m_rear*self.R**2,
                   self.l_fr:0.2,
                   self.l_rear:0.01,
                   self.u0:0.005,
                   self.z_c3:0.01,
                   self.l_bumps:0.15,
                   self.amplitude:0.0165,
                   self.length:0.19,
                   self.speed:1.7,
                   self.axw:0.47,
                   self.c_ft:107,
                   self.c_rt:5500,
                   self.c_fs:107,
                   self.c_rs:107,
                   self.c:100,
                   self.A:0.001*10,
                   self.omega:4*np.pi,
                          }
        
        return default_data_dict  
    def get_table_values(self):
        default_data_dict={
                   
                   self.c:100,
                   self.c_rt:5500,
                   self.m_rear:1.5,
                   self.M:75,
                   self.k_ft:300000,
                   self.k_r:750000,
                          }
        
        return default_data_dict  
    def symbols_description(self):
        self.sym_desc_dict = {
              self.m_fr:'front wheel mass',
              self.m_rear:'rear wheel mass',
              self.M:'''mass of the main body including user''',
              self.k_r:'rear wheel stiffness',
              self.k_rt:'rear wheel tire stiffness',
              self.k_f:'front wheel stiffness',
              self.k_ft:'front wheel tire stiffness',
              self.k_rot:'rotational stiffness',
              self.m:'wheelchair mass',
              self.k:'wheelchair stiffness',
              self.g:'gravitational acceleration',
              self.F_1:'initial force',
              self.F_2:'initial force',
              self.Omega:'driving force frequency',
              self.F:'disabled driver arm force',
              self.R:'rear wheel radius',
              self.v0:'wheelchair initial velocity ',
              self.u0:'road profil amplitude',
              self.I_ch:'''wheelchair's moment of inertia''',
              self.I_w:'''rear wheel moment of inertia''',
              self.z_c3:'''vertical distance between wheelchair's center of mass and constant reference height''',
              self.l_fr:'front wheelchair spring initial length',
              self.l_rear:'rear wheelchair spring initial length',
#               self.m_RC:'RapidChair drive rocker arm mass',
#               self.I_RC:'RapidChair drive moment of inertia',
#               self.l_RC:'Rapidchair drive rocker arm length',
#               self.k_RC:'RapidChair drive stiffness',
#               self.phi0:'initial rocker arm axis angle relative to horizontal plane',
#               self.m_w:'RapidChair drive wheel mass',
#               self.I_wrc:'RapidChair drive wheel moment of inertia',
#               self.r_w:'RapidChair drive wheel radius', #do sprawdzenia, nie jestem pewien
#               self.k_w:'RapidChair drive wheel striffness',
#               self.k_fix:'fastening stiffness',
#               self.k_tire:'RapidChair drive wheel tire striffness',
              self.c:'general resistance to motion coefficient',
              self.c_fs:'vertical damping coefficient of front suspension',
              self.c_rs:'vertical damping coefficient of rear suspension',
              self.c_ft:'vertical damping coefficient of caster band',
              self.c_rt:'vertical damping coefficient of rear tire',
#               self.c_mu:'inertia damping decrement',
#               self.c_lam:'stiffness damping decrement',
              self.xm_3:'center of mass location in x axis',
              self.l_l: r'offset of left spring',
              self.l_r: r'offset of right spring',
#               self.a_ox:'acceleration in x axis - longitudal axis', # acc osi
#               self.a_oz:'acceleration in z axis - vertical axis',
#               self.a_rz:'acceleration of wishbone',
#               self.a_rcz:'????',
              self.t_l:'phase shift of the wheels',
              self.delta_t:'the moment of occurrence of the obstacle',
              self.l_ramp:'ramp length',
              self.l_bumps:'bumps length',
              self.amplitude:'????',
              self.length:'????',
              self.speed:'steady state velocity',
              self.axw:'wheelbase',
              self.pm:'force duty cycle',
              self.A:'ground forcing amplitude',
              self.Leng:'????',
            
            
               self.x:'main body center of mass horizontal displacement',
               self.z_fr:'front wheel vertical displacement',
               self.z_rear:'rear wheel vertical displacement',
               self.z:'main body center of mass vertical displacement',
#                self.z_wrc:'vertical displacement of the RapidChair driven wheel',
               self.phi:'angular displacement of the main mass',
#                self.phi_rc:'angular displacement of the RapidChair drive',
#                self.theta:'angular displacement of the RapidChair driven wheel',
            
               self.x.diff(t):'main body center of mass horizontal velocity',
               self.z_fr.diff(t):'front wheel vertical velocity',
               self.z_rear.diff(t):'rear wheel vertical velocity',
               self.z.diff(t):'main body center of mass vertical velocity',
#                self.z_wrc.diff(t):'vertical velocity of the RapidChair driven wheel',
               self.phi.diff(t):'angular velocity of the main mass',
#                self.phi_rc.diff(t):'angular velocity of the RapidChair drive',
#                self.theta.diff(t):'angular velocity of the RapidChair driven wheel',
            

            
               self.x.diff(t,t):'main body center of mass horizontal acceleration',
               self.z_fr.diff(t,t):'front wheel vertical acceleration',
               self.z_rear.diff(t,t):'rear wheel vertical acceleration',
               self.z.diff(t,t):'main body center of mass vertical acceleration',
#                self.z_wrc.diff(t,t):'vertical acceleration of the RapidChair driven wheel',
               self.phi.diff(t,t):'angular acceleration of the main mass',
#                self.phi_rc.diff(t,t):'angular acceleration of the RapidChair drivetrain',
#                self.theta.diff(t,t):'angular acceleration of the RapidChair driven wheel',
#                self.T:'overall kinetic energy',
#                self.V:'overall potential energy',
               Symbol('L',positive=True):'''Lagrange's function''',
               Symbol('v',positive=True):'''steady state velocity''',
               self.ivar:'time',
               self.T:'total kinetic energy',
               self.V:'total potential energy',
               self.D:'total dissipative potential',
               self.omega:'frequency of road profile'
        }
        return self.sym_desc_dict
#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m: r'Masa resorowana',
#             self.m_p: r'Masa przedniej osi',
#             self.m_l: r'Masa tylnej osi',
#             self.m_b: r'Masa baterii trakcyjnej',
#             self.m_b: r'Masa pasazera',
#             self.I: r'Moment bezwładności',
#             self.l_rod: r'Długość osi',
#             self.l_p: r'Odległość przedniej osi do środka ciężkości autobusu',
#             self.l_l: r'Odległość tylnej osi do środka ciężkości autobusu',
#             self.l_b: r'Odległość środka cieżkości baterii do środka ciężkości autobusu',
#             self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
#             self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
#             self.k_pw: r'Współczynnik sztywności dla przedniej opony',
#             self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
#             self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
#             self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
#             self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
#             self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
#             self.phi:  r'Kąt obrotu masy resorowanej',
#             self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
#             self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
#             self.z_lw: r'Przemieszczenie pionowe tylnego koła',
#             self.z_b: r'Przemieszczenie pionowe baterii',
#             self.A: r'Amplituda siły wymuszającej',
#             self.omega: r'Częstość siły wymuszającej',
#             self.ivar: r'Czas',
#             self.f: r'Częstotliwość wymuszenia',
            
            
#         }
#         return self.sym_desc_dict
    
sys4=DampedChair4DOF()
y=[sys4.q,sys4.z, sys4.z.diff(t,t)]
ic_list = [0.1,0.1,0,0,0.1,0.1,0,0]
units_dict = {
               sys4.F:ureg.kilogram*ureg.meter/ureg.second**2,
                   #dyn_sys.c_mu:0.0001,
                   #dyn_sys.c_lam:0.0001,
                   sys4.l_l:ureg.meter,
                   sys4.l_r:ureg.meter,
                   sys4.k_f:ureg.newton/ureg.meter,
                   sys4.k_ft:ureg.kilogram/ureg.second**2,
                   sys4.k_r:ureg.newton/ureg.meter,
                   sys4.k_rt:ureg.kilogram/ureg.second**2,
                   sys4.c:ureg.newton*ureg.second/ureg.meter,
                    sys4.c_fs:ureg.newton*ureg.second/ureg.meter,
                    sys4.c_rs:ureg.newton*ureg.second/ureg.meter,
                    sys4.c_ft:ureg.newton*ureg.second/ureg.meter,
                    sys4.c_rt:ureg.newton*ureg.second/ureg.meter,
                   sys4.I_ch:ureg.meter**2*ureg.kilogram,
                    sys4.M:ureg.kilogram,
                   sys4.m_rear:ureg.kilogram,
                   sys4.m_fr:ureg.kilogram,
                   sys4.Omega:ureg.radian/ureg.second,
                   sys4.R:ureg.meter,
                   #dyn_sys.z_c3:0.4,
                   sys4.g:ureg.meter/ureg.second**2,
                   sys4.I_w:ureg.meter**2*ureg.kilogram,
                   sys4.l_fr:ureg.meter,
                   sys4.l_rear:ureg.meter,
                   #dyn_sys.u0:,
                   sys4.l_bumps:ureg.meter,
                   sys4.amplitude:ureg.meter,
                   sys4.length:ureg.meter,
                   sys4.speed:ureg.meter/ureg.second,
                   #self.axw:0.47,
                   sys4.Leng:ureg.meter,
               sys4.m:ureg.kilogram,
               sys4.k:(ureg.newton / ureg.meter),
               t:ureg.second,
                f:ureg.hertz,
                sys4.z:ureg.meter,
                sys4.phi: ureg.radian,
                sys4.z.diff(t,2):ureg.meter/ureg.second**2,
                sys4.phi.diff(t,2): ureg.radian/ureg.second**2,
                sys4.pm:ureg.meter/ureg.meter,
                sys4.u0:ureg.meter,
              }

# units_dict = {sys4.m:ureg.kilogram,
#               sys4.m_p:ureg.kilogram,
#               sys4.m_l:ureg.kilogram,
#               sys4.m_b:ureg.kilogram,
#               sys4.m_k:ureg.kilogram,
#               sys4.A:ureg.meter,
#               sys4.k_l:(ureg.newton / ureg.meter),
#               sys4.k_lw:(ureg.newton / ureg.meter),
#               sys4.k_p:(ureg.newton / ureg.meter),
#               sys4.k_pw:(ureg.newton / ureg.meter),
#               sys4.c_l:(ureg.newton*ureg.second / ureg.meter),
#               sys4.c_lw:(ureg.newton*ureg.second / ureg.meter),
#               sys4.c_p:(ureg.newton*ureg.second / ureg.meter),
#               sys4.c_pw:(ureg.newton*ureg.second / ureg.meter),
#               sys4.omega:ureg.radian,
#               sys4.l_l:ureg.meter,
#               sys4.l_p:ureg.meter,
#               sys4.l_b:ureg.meter,
#               sys4.l_rod:ureg.meter,
#               sys4.I:(ureg.kilogram*ureg.meter*ureg.meter),
#               sys4.phi:ureg.radian,
#               sys4.z:ureg.meter,
#               sys4.z.diff(t,t):ureg.meter/ureg.second**2,
#               sys4.z_lw:ureg.meter,
#               sys4.z_pw:ureg.meter,
#               sys4.z_b:ureg.meter,
#               t:ureg.second,
#               f:ureg.hertz,
              
#              }

unit=units_dict
        
class DampedChairDDOF(ComposedSystem):

   
    
    z,phi,z_pw,z_lw=dynamicsymbols('z, \\varphi z_pw z_lw')
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_p=Symbol('m_p', positive=True)
    m_l=Symbol('m_l', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_p=Symbol('l_p', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_l=Symbol('k_l1', positive=True)
    k_lw=Symbol('k_lw', positive=True)
    k_p=Symbol('k_p1', positive=True)
    k_pw=Symbol('k_pw', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_p=Symbol('c_p', positive=True)
    c_lw=Symbol('c_lw', positive=True)
    c_pw=Symbol('c_pw', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    ivar=Symbol('t')
    
    def __init__(self,
                  m=None,
                 m_p=None,
                 m_l=None,
                 m_b=None,
                 I=None,
                 l_rod=None,
                 l_l=None,
                 l_p=None,
                 l_b=None,
                 k_l=None,
                 k_lw=None,
                 k_p=None,
                 k_pw=None,
                 c_l=None,
                 c_p=None,
                 c_lw=None,
                 c_pw=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 s=None,
                 phi=None,
                 z_pw=None,
                 z_lw=None,
                 A=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        

        if m is not None: self.m = m  # mass of a rod
        if m_p is not None: self.m_w = m_w
        if m_l is not None: self.m_w = m_w
        if m_b is not None: self.m_b = m_b
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_p is not None: self.l_p = l_p #offset of right spring
        if l_b is not None: self.l_p = l_b  
        if l_rod is not None: self.l_rod = l_rod
        if k_l is not None: self.k_l1 = k_l1
        if k_lw is not None: self.k_lw = k_lw
        if k_p is not None: self.k_p1 = k_p1
        if k_pw is not None: self.k_pw = k_pw
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_lw is not None: self.c_lw=c_lw
        if c_pw is not None: self.c_pw=c_pw
        if A is not None: self.A=A
        if omega is not None: self.omega=omega
        if I is not None: self.I = I  # moment of inertia of a rod
        if F_engine is not None: self.F_engine = F_engine
            
        self.s=self.A*sin(ivar*self.omega)
        
        self.z_l=self.z+self.phi*self.l_l
        self.z_p=self.z-self.phi*self.l_p
        self.z_lw=(self.k_l*self.z_l)/(self.k_lw+self.k_l)
        self.z_pw=(self.k_p*self.z_p)/(self.k_pw+self.k_p)
        
        
        self.right_wheel = MaterialPoint(self.m_p, pos1=self.z_pw, qs=[self.z, self.phi])
        self.left_wheel = MaterialPoint(self.m_l, pos1=self.z_lw, qs=[self.z, self.phi])
        self.spring_pw = Spring(self.k_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        self.spring_lw = Spring(self.k_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        self.damper_pw = Damper(self.c_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        self.damper_lw = Damper(self.c_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring=Spring(k,pos1,po2,qs)
        
        self.body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=[self.z, self.phi])  # rod
        self.battery= MaterialPoint(self.m_b, pos1=self.z+self.phi*self.l_b, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_l, pos1=self.z_l , pos2 = self.z_lw , qs=[self.z, self.phi])  # left spring
        self.spring_2 = Spring(self.k_p, pos1=self.z_p, pos2=self.z_pw, qs=[self.z, self.phi])
        self.damper_1=Damper(self.c_l, pos1=self.z_p , pos2 = self.z_lw , qs=[self.z, self.phi])
        self.damper_2 = Damper(self.c_p, pos1=self.z_l , pos2 = self.z_pw , qs=[self.z, self.phi])

        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body +self.battery+ self.spring_1 + self.spring_2 + self.damper_1+ self.damper_2 + self.right_wheel + self.left_wheel + self.spring_pw + self.spring_lw+self.damper_lw+self.damper_pw
        
#         display(type(system))
        
#         system_new = system.subs({self.z_lw:self.z_l2,self.z_pw:self.z_p2}) 
        
#         display(type(system_new))
#         display(type(self.body))
        
        #super().__init__(system_new,**kwargs)

        super().__init__(**{'system':system,**kwargs})
        


        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Masa resorowana',
            self.m_p: r'Masa przedniej osi',
            self.m_l: r'Masa tylnej osi',
            self.m_b: r'Masa baterii trakcyjnej',
            self.I: r'Moment bezwładności',
            self.l_rod: r'Długość osi',
            self.l_p: r'Odległość od przedniej osi do środka ciężkości autobusu',
            self.l_l: r'Odległość od tylnej osi do środka ciężkości autobusu',
            self.l_b: r'Odległość od środka cieżkości baterii do środka ciężkości autobusu',
            self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
            self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
            self.k_pw: r'Współczynnik sztywności dla przedniej opony',
            self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
            self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
            self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
            self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
            self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
            self.phi:  r'Kąt obrotu masy resorowanej',
            self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
            self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
            self.z_lw: r'Przemieszczenie pionowe tylnego koła',
            self.z_b: r'Przemieszczenie pionowe baterii',
            self.A: r'Amplituda siły wymuszającej',
            self.omega: r'Częstość siły wymuszającej',
            self.ivar: r'Czas',
            
        }
        return self.sym_desc_dict
    
sys2=DampedChairDDOF()



units_dict1 = {sys2.m:ureg.kilogram,
              sys2.m_p:ureg.kilogram,
              sys2.m_l:ureg.kilogram,
              sys2.m_b:ureg.kilogram,
              sys2.A:ureg.meter,
              sys2.k_l:(ureg.newton / ureg.meter),
              sys2.k_lw:(ureg.newton / ureg.meter),
              sys2.k_p:(ureg.newton / ureg.meter),
              sys2.k_pw:(ureg.newton / ureg.meter),
              sys2.c_l:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_lw:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_p:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_pw:(ureg.newton*ureg.second / ureg.meter),
              sys2.omega:ureg.radian,
              sys2.l_l:ureg.meter,
              sys2.l_p:ureg.meter,
              sys2.l_b:ureg.meter,
              sys2.l_rod:ureg.meter,
              sys2.I:(ureg.kilogram*ureg.meter*ureg.meter),
              sys2.phi:ureg.radian,
              sys2.z:ureg.meter,
              sys2.z.diff(t,t):ureg.meter/ureg.second**2,
              sys2.z_lw:ureg.meter,
              sys2.z_pw:ureg.meter,
              sys2.z_b:ureg.meter,
              t:ureg.second,
              f:ureg.hertz,
             }

unit1=units_dict1
        
class DampedChairSimplifiedDDOF2(ComposedSystem):

   
    
    z,phi,z_pw,z_lw=dynamicsymbols('z, \\varphi z_pw z_lw')
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_p=Symbol('m_p', positive=True)
    m_l=Symbol('m_l', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_p=Symbol('l_p', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_l=Symbol('k_l1', positive=True)
    k_lw=Symbol('k_lw', positive=True)
    k_p=Symbol('k_p1', positive=True)
    k_pw=Symbol('k_pw', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_p=Symbol('c_p', positive=True)
    c_lw=Symbol('c_lw', positive=True)
    c_pw=Symbol('c_pw', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    k_l_zas=Symbol('k_l_zas', positive=True)
    k_p_zas=Symbol('k_p_zas', positive=True)
    c_l_zas=Symbol('c_l_zas', positive=True)
    c_p_zas=Symbol('c_lp_zas', positive=True)
    ivar=Symbol('t')
    
    def __init__(self,
                 m=None,
                 m_p=None,
                 m_l=None,
                 m_b=None,
                 I=None,
                 l_rod=None,
                 l_l=None,
                 l_p=None,
                 l_b=None,
                 k_l=None,
                 k_lw=None,
                 k_p=None,
                 k_pw=None,
                 c_l=None,
                 c_p=None,
                 c_lw=None,
                 c_pw=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 phi=None,
                 z_pw=None,
                 z_lw=None,
                 A=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
                 s=None,
                 k_l_zas=None,
                 k_p_zas=None,
                 c_l_zas=None,
                 c_p_zas=None,
                 
                 **kwargs):
        if z is not None: self.z=z
        if z_lw is not None: self.z_lw=z_lw
        if z_pw is not None: self.z_pw=z_pw
        if z_l is not None: self.z_lw=z_lw
        if z_p is not None: self.z_pw=z_pw
        if z_b is not None: self.z_b=z_b
        if phi is not None: self.phi=phi

        if m is not None: self.m = m  # mass of a rod
        if m_p is not None: self.m_w = m_w
        if m_l is not None: self.m_w = m_w
        if m_b is not None: self.m_b = m_b
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_p is not None: self.l_p = l_p #offset of right spring
        if l_b is not None: self.l_p = l_b  
        if l_rod is not None: self.l_rod = l_rod
        if k_l is not None: self.k_l1 = k_l1
        if k_lw is not None: self.k_lw = k_lw
        if k_p is not None: self.k_p1 = k_p1
        if k_pw is not None: self.k_pw = k_pw
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_lw is not None: self.c_lw=c_lw
        if c_pw is not None: self.c_pw=c_pw
        if A is not None: self.A=A
        if omega is not None: self.omega=omega
        if I is not None: self.I = I  # moment of inertia of a rod
        if F_engine is not None: self.F_engine = F_engine
        
                 
        self.s=self.A*sin(ivar*self.omega)
        
        self.k_l_zas=(self.k_l*self.k_lw)/(self.k_l+self.k_lw)
        self.k_p_zas=(self.k_p*self.k_pw)/(self.k_p+self.k_pw)
        self.c_l_zas=(self.c_l*self.c_lw)/(self.c_l+self.c_lw)
        self.c_p_zas=(self.c_p*self.c_pw)/(self.c_p+self.c_pw)
        
        self.z_l=self.z+self.phi*self.l_l
        self.z_p=self.z-self.phi*self.l_p
        self.z_lw=(self.k_l*self.z_l)/(self.k_lw+self.k_l)
        self.z_pw=(self.k_p*self.z_p)/(self.k_pw+self.k_p)
        
        #self.right_wheel = MaterialPoint(self.m_w, pos1=self.z_pw, qs=[self.z, self.phi])
        #self.left_wheel = MaterialPoint(self.m_w, pos1=self.z_lw, qs=[self.z, self.phi])
        #self.spring_pw = Spring(self.k_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring_lw = Spring(self.k_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring=Spring(k,pos1,po2,qs)
        
        self.body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=[self.z, self.phi])  # rod
        self.battery= MaterialPoint(self.m_b, pos1=self.z+self.phi*self.l_b, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_l_zas, pos1=self.z_l , pos2 = self.s , qs=[self.z, self.phi])  # left spring
        self.spring_2 = Spring(self.k_p_zas, pos1=self.z_p, pos2=self.s, qs=[self.z, self.phi])
        self.damper_1=Damper(self.c_l_zas, pos1=self.z_p , pos2 = self.s, qs=[self.z, self.phi])
        self.damper_2 = Damper(self.c_p_zas, pos1=self.z_l , pos2 = self.s , qs=[self.z, self.phi])

        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body +self.battery+ self.spring_1 + self.spring_2 + self.damper_1+ self.damper_2 
        
#         display(type(system))
        
#         system_new = system.subs({self.z_lw:self.z_l2,self.z_pw:self.z_p2}) 
        
#         display(type(system_new))
#         display(type(self.body))
        
        #super().__init__(system_new,**kwargs)

        super().__init__(**{'system':system,**kwargs})        

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Masa resorowana',
            self.m_p: r'Masa przedniej osi',
            self.m_l: r'Masa tylnej osi',
            self.m_b: r'Masa baterii trakcyjnej',
            self.I: r'Moment bezwładności',
            self.l_rod: r'Długość osi',
            self.l_p: r'Odległość od przedniej osi do środka ciężkości autobusu',
            self.l_l: r'Odległość od tylnej osi do środka ciężkości autobusu',
            self.l_b: r'Odległość od środka cieżkości baterii do środka ciężkości autobusu',
            self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
            self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
            self.k_pw: r'Współczynnik sztywności dla przedniej opony',
            self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
            self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
            self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
            self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
            self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
            self.phi:  r'Kąt obrotu masy resorowanej',
            self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
            self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
            self.z_lw: r'Przemieszczenie pionowe tylnego koła',
            self.z_b: r'Przemieszczenie pionowe baterii',
            self.A: r'Amplituda siły wymuszającej',
            self.omega: r'Częstość siły wymuszającej',
            self.ivar: r'Czas',
            self.k_l_zas: r'Zastępcza wartość współczynnika sztywności zawieszenia tylnej osi',
            self.k_p_zas: r'Zastępcza wartość współczynnika sztywności zawieszenia przedniej osi',
            self.c_l_zas: r'Zastępcza wartość współczynnika tłumienia zawieszenia tylnej osi',
            self.c_p_zas: r'Zastępcza wartość współczynnika tłumienia zawieszenia przedniej osi',
        }
        return self.sym_desc_dict
                 
sys22=DampedChairSimplifiedDDOF2()      
                 
units_dict2 = {sys22.m:ureg.kilogram,
              sys22.m_p:ureg.kilogram,
              sys22.m_l:ureg.kilogram,
              sys22.m_b:ureg.kilogram,
              sys22.A:ureg.meter,
              sys22.k_l:(ureg.newton / ureg.meter),
              sys22.k_lw:(ureg.newton / ureg.meter),
              sys22.k_p:(ureg.newton / ureg.meter),
              sys22.k_pw:(ureg.newton / ureg.meter),
              sys22.c_l:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_lw:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_p:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_pw:(ureg.newton*ureg.second / ureg.meter),
              sys22.omega:ureg.radian,
              sys22.l_l:ureg.meter,
              sys22.l_p:ureg.meter,
              sys22.l_b:ureg.meter,
              sys22.l_rod:ureg.meter,
              sys22.I:(ureg.kilogram*ureg.meter*ureg.meter),
              sys22.phi:ureg.radian,
              sys22.z:ureg.meter,
              sys22.z.diff(t,t):ureg.meter/ureg.second**2,
              sys22.z_lw:ureg.meter,
              sys22.z_pw:ureg.meter,
              sys22.z_b:ureg.meter,
              t:ureg.second,
              f:ureg.hertz,
              sys22.k_p_zas:(ureg.newton / ureg.meter),
              sys22.k_l_zas:(ureg.newton / ureg.meter),
              sys22.c_l_zas:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_p_zas:(ureg.newton*ureg.second / ureg.meter),
             }

unit2=units_dict2
   

