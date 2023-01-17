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
from functools import cached_property




        


    
class UndampedVehicleSuspension(ComposedSystem):

    scheme_name = 'car_undamped.png'
    real_name = 'car_real.jpg'

    m=Symbol('m', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_r=Symbol('l_r', positive=True)
                
    r=Symbol('k_r', positive=True)
                
    k_l=Symbol('k_l', positive=True)
    k_r=Symbol('k_r', positive=True)
                
    F_engine=Symbol('F_{engine}', positive=True)
    Omega=Symbol('Omega',positive=True)
    ivar=Symbol('t', positive=True)


    
          
    #qs=dynamicsymbols('z, \\varphi')
    phi=dynamicsymbols('\\varphi')
    z=dynamicsymbols('z')

    #default symbols for substitution
    m0=Symbol('m_0', positive=True)
    l0=Symbol('l_0', positive=True)
    c0=Symbol('c_0', positive=True)
    k0=Symbol('k_0', positive=True)
    lam0=Symbol('\lambda_0', positive=True)
    Omega0=Symbol('Omega_0', positive=True)
    F0=Symbol('F_0', positive=True)
    

    def __init__(self,
                 m=None,
                 I=None,
                 l_rod=None,
                 l_r=None,
                 k_r=None,
                 k_l=None,
                 F_engine=None,
                 l_l=None,
                 phi=None,
                 z=None,
                 ivar=Symbol('t'),
                 qs=None,
                 
                 **kwargs):
  
        if  m is not None: self.m =m # mass of a rod
        if  I is not None: self.I = I # moment of inertia of a rod
        if  l_rod is not None: self.l_rod =l_rod# length of a rod
        if  l_r is not None: self.l_r = l_r 
        if  l_l is not None: self.l_l = l_l 
        if  k_r is not None: self.k_r =k_r# left spring
        if  k_l is not None: self.k_l = k_l# right spring
        if  F_engine is not None: self.F_engine = F_engine
        if  phi is not None: self.phi = phi
        if z is not None: self.z=z    
          
            # self.z, self.phi = self.qs
        self.qs = [self.z,self.phi]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=self.qs)(label='rod')
        self._spring_l = Spring(self.k_l, pos1=self.z + self.phi * self.l_l, qs=self.qs)(label='left spring')
        self._spring_r = Spring(self.k_r, pos1=self.z - self.phi * self.l_r, qs=self.qs)(label='right spring')
        self._force = Force(self.F_engine, pos1=self.z - self.l_r * self.phi, qs=self.qs)(label='force')


        components['_body'] = self._body
        components['_spring_l'] = self._spring_l
        components['_spring_r'] = self._spring_r
        components['_force'] = self._force
        
        
        
        return components

    
#    def max_dynamic_force_pin(self):
#        return self.frequency_response_function()*self.stiffness_matrix()[0]
    
#    def dynamic_bearing_force(self):
       
#       return self.max_dynamic_force_pin() * self.sqrt(L)
    def max_static_force_pin(self):
        return self.components['_spring_l'].subs({coord:0 for coord in self.q})

    

    
    
    
    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])/2
    
    def static_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_static_force_pin())/(pi*kt*Re))**(1/2)

#     def max_dynamic_force_pin(self):
#         return (self._frf()[0] + self._frf()[1]*self.l_l)*self.k_l

    
    def max_dynamic_force_pin(self):
        
        amps = self._frf()
        force = self.components['_spring_l'].force().doit().expand()#.doit()
        #print('in class method')
        #display(force)
        #print('in class method')
        
        return self.components['_spring_l'].force().subs({coord:amp for amp,coord in zip(amps,self.q)})



    def dynamic_bearing_force(self):
        L=Symbol('L')#wymagana trwałość
        return self.max_dynamic_force_pin() * L**(S.One/3)
    
    def dynamic_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_dynamic_force_pin())/(pi*kt*Re))**(1/2)
    

    def get_default_data(self):

        #m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0', positive=True)
        c0, k_0, l_l0, omega, F_0 = symbols('c_0 k_0 l_l0 Omega F_0', positive=True)
        m0, l0  = self.m0,self.l0
        c0, k0= self.c0,self.k0
        lam0=self.lam0
        Omega0=self.Omega0
        F0=self.F0

        default_data_dict = {
            self.l_l: [self.l_r],
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            self.I: [S.One/12*self.m *self.l_rod**2],
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            #self.m: [m0,2*m0,3*m0,4*m0,5*m0,6*m0,7*m0,8*m0,9*m0],
            self.m: [m0*S.One*no for no in range(1,8)],
            #self.c_r: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
#             self.c_r: [self.lam*(self.k_r)],
            #self.k: [k0*S.One*no for no in range(1,8)],
#             self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.k_l: [k0*S.One*no for no in range(1,8)],
            self.k_r: [k0*S.One*no for no in range(1,8)],
            
            
#             self.c_l: [self.lam*(self.k_l)],
            self.l_r: [l0*S.One*no for no in range(1, 8)],
            #self.l_l: [l0*S.One*no for no in range(1, 8)],
           
            self.Omega: [self.Omega],
            self.F_engine: [F0*cos(self.Omega*self.ivar)*S.One*no for no in range(1,8)],
            
          #  self.F_engine: [
           #     2 * F_0 * sin(omega * self.ivar),
            #    3 * F_0 * sin(omega * self.ivar),
             #   4 * F_0 * sin(omega * self.ivar),
              #  5 * F_0 * sin(omega * self.ivar),
                #6 * F_0 * sin(omega * self.ivar)
          #  ]
        }

        return default_data_dict
    
class DampedVehicleSuspension(UndampedVehicleSuspension):

    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'

    c_l=Symbol('c_l', positive=True)
    c_r=Symbol('c_r', positive=True)
                
    l_cl=Symbol('l_{cl}', positive=True)
    l_cr=Symbol('l_{cr}', positive=True)
    lam=Symbol('lambda',positive=True)
    phi=dynamicsymbols('\\varphi')
    z=dynamicsymbols('z')


    def __init__(self,
                 m=None,
                 I=None,
                 l_rod=None,
                 l_r=None,
                 k_r=None,
                 k_l=None,
                 F_engine=None,
                 c_l=None,
                 c_r=None,
                 l_cl=None,
                 l_l=None,
                 l_cr=None,
                 phi=None,
                 z=None,
                 ivar=Symbol('t'),
                 qs=None,
                 **kwargs):

        if  l_cr is not None: self.l_cr = l_cr
        if  l_cl is not None: self.l_cl = l_cl

        if  c_l is not None: self.c_l =c_l
        if  c_r is not None: self.c_r =  c_r


        self.qs = [self.phi,self.z]

        super().__init__(m=m,I=I, l_rod= l_rod,l_r=l_r,k_r=k_r,k_l=k_l, F_engine= F_engine,z=z,l_l=l_l,ivar=ivar,**kwargs)


    @property
    def components(self):

        components =super().components

        
        self._damper_l = Damper(self.c_l, pos1=self.z + self.phi * self.l_l,
                               qs=self.qs)(label='left damper')
        self._damper_r = Damper(c=self.c_r, pos1=self.z - self.phi * self.l_r,
                               qs=self.qs)(label='right damper')



        components['_damper_l'] = self._damper_l
        components['_damper_r'] = self._damper_r
       

        
        return components

        
#    def max_dynamic_force_pin(self):
#        return self.frequency_response_function()*self.stiffness_matrix()[0]
    
#    def dynamic_bearing_force(self):
       
#       return self.max_dynamic_force_pin() * self.sqrt(L)
    
    def max_dynamic_force_pin(self):
        return abs(self.frequency_response_function()*self.stiffness_matrix()[0])

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])
    
    def static_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_static_force_pin())/(pi*kt*Re))**(1/2)

#     def max_dynamic_force_pin(self):
#         return (self._frf()[0].subs(self.lam0,0) + self._frf()[1].subs(self.lam0,0)*self.l_rod/2)*self.stiffness_matrix()[0]

    def dynamic_bearing_force(self):
        L=Symbol('L')#wymagana trwałość
        return self.max_dynamic_force_pin() * L**(S.One/3)
    
    def dynamic_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_dynamic_force_pin())/(pi*kt*Re))**(1/2)
    
    
    def dynamic_bearing_force(self):
        L=Symbol('L')
        return self.max_dynamic_force_pin() * self.sqrt(L)
    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.I: r'',
            self.l_rod: r'',
            self.l_r: r'',
            self.k_r: r'',
            self.F_engine: r'',
            self.k_l: r'',
            self.c_l: r'',
            self.c_r: r'',
            self.l_cl: r'',
            self.l_l: r'',
            self.c_r: r'',
            self.phi: r'',
            self.Z: r'',
            self.ivar: r'',
            self.qs: r'',
        }
    
    def get_default_data(self):

        #m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0', positive=True)
        c0, k_0, l_l0, Omega0, F_0 = symbols('c_0 k_0 l_l0 Omega0 F_0', positive=True)
        m0, l0  = self.m0,self.l0
        c0, k0= self.c0,self.k0
        lam0=self.lam0
        Omega0=self.Omega0
        F0=self.F0

        default_data_dict = {
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            self.I: [S.One/12*self.m *self.l_rod],
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            #self.m: [m0,2*m0,3*m0,4*m0,5*m0,6*m0,7*m0,8*m0,9*m0],
            self.m: [m0*S.One*no for no in range(1,8)],
            #self.c_r: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.c_r: [self.lam*(self.k_r)],
            #self.k: [k0*S.One*no for no in range(1,8)],
            self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.k_l: [k0*S.One*no for no in range(1,8)],
            self.k_r: [k0*S.One*no for no in range(1,8)],
            
            
            self.c_l: [self.lam*(self.k_l)],
            self.l_r: [l0*S.One*no for no in range(1, 8)],
            self.l_l: [l0*S.One*no for no in range(1, 8)],
           
            self.Omega: [Omega0*S.One*no for no in range(1,2)],
            self.F_engine: [F0*sin(self.Omega*self.ivar)*S.One*no for no in range(1,8)],}
            


        return default_data_dict

class DampedSymmetricalVehicleSuspension(DampedVehicleSuspension):
    def get_default_data(self):

        #m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0', positive=True)
        c0, k_0, l_l0, Omega0, F_0 = symbols('c_0 k_0 l_l0 Omega0 F_0', positive=True)
        m0, l0  = self.m0,self.l0
        c0, k0= self.c0,self.k0
        lam0=self.lam0
        Omega0=self.Omega0
        F0=self.F0

        default_data_dict = {
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            self.I: [S.One/12*self.m *self.l_rod**2],
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            #self.m: [m0,2*m0,3*m0,4*m0,5*m0,6*m0,7*m0,8*m0,9*m0],
            self.m: [m0*S.One*no for no in range(1,8)],
            #self.c_r: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.c_r: [self.lam*(self.k_r)],
            #self.k: [k0*S.One*no for no in range(1,8)],
            self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.k_l: [self.k_r],
            self.k_r: [k0*S.One*no for no in range(1,8)],
            
            
            self.c_l: [self.lam*(self.k_l)],
            self.l_r: [self.l_l],
            self.l_l: [l0*S.One*no for no in range(1, 8)],
           
            self.Omega: [Omega0*S.One*no for no in range(1,2)],
            self.F_engine: [F0*cos(self.Omega*self.ivar)*S.One*no for no in range(1,8)],
            
          #  self.F_engine: [
           #     2 * F_0 * sin(omega * self.ivar),
            #    3 * F_0 * sin(omega * self.ivar),
             #   4 * F_0 * sin(omega * self.ivar),
              #  5 * F_0 * sin(omega * self.ivar),
                #6 * F_0 * sin(omega * self.ivar)
          #  ]
        }

        return default_data_dict

