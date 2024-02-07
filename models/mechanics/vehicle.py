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
    ivar=Symbol('t')
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
                 l_l=None,
                 F_engine=None,
                 Omega=None,
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
        if  Omega is not None: self.Omega = Omega
        if  phi is not None: self.phi = phi
        if z is not None: self.z=z
        # self.z, self.phi = self.qs
        self.qs = [self.z,self.phi]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=self.qs)(label='rod')
        self._left_mount = Spring(self.k_l, pos1=self.z + self.phi * self.l_l, qs=self.qs)(label='left spring')
        self._right_mount = Spring(self.k_r, pos1=self.z - self.phi * self.l_r, qs=self.qs)(label='right spring')
        self._force = Force(self.F_engine*cos(self.Omega*self.ivar), pos1=self.z - self.l_r * self.phi, qs=self.qs)(label='force')
        self._static_force = Force(self.F_engine, pos1=self.z - self.l_r * self.phi, qs=self.qs)(label='force')

        components['_body'] = self._body
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        components['_force'] = self._force
        components['_static_force'] = self._static_force

        return components

    
#    def max_dynamic_force_pin(self):
#        return self.frequency_response_function()*self.stiffness_matrix()[0]
    
#    def dynamic_bearing_force(self):
       
#       return self.max_dynamic_force_pin() * self.sqrt(L)

#     def max_static_force_pin(self):
#         return self.components['_spring_l'].subs({coord:0 for coord in self.q})

    
    def max_static_force_pin(self):

        ans=self.static_force()
        
        return abs(ans)


    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def static_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        return ((4*self.max_static_force_pin())/(pi*kt*Re))**(1/2)


#     def max_dynamic_force_pin(self):
#         frf=self._frf
#         return (frf[0] + frf[1]*self.l_l)*self.k_l

    
#     def max_dynamic_force_pin(self):
        
#         amps = self._frf()
#         force = self.components['_spring_l'].force().doit().expand()#.doit()
#         data=self._given_data
#         display(abs(self.components['_spring_l'].force()))
#         return abs((self.components['_spring_l'].force().subs({coord:amp for amp,coord in zip(amps,self.q)})).subs(data)).doit()
    
    
    def dynamic_force(self):
        
        amps = self._fodes_system.steady_solution.as_dict()
        force = self.components['_left_mount'].force().doit().expand()#.doit()
        data=self._given_data
        
        #display(abs(self.components['_spring_l'].force()))
        return (self.components['_left_mount'].force().subs(amps)).subs(data).expand().doit()
    
    def max_static_force(self):
        return self.static_force()

    def max_dynamic_force_pin(self):
        
        ans = self.dynamic_force()
        
        #display(abs(self.components['_spring_l'].force()))
        data=self._given_data
        sin_coeff=ans.coeff(sin(self.Omega*self.ivar)).subs(data)
        #display(sin_coeff)
        cos_coeff=ans.coeff(cos(self.Omega*self.ivar)).subs(data)
        #display(cos_coeff)
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        #display(free_coeff)

        return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)
    
    def max_dynamic_force(self):
        return self.max_dynamic_force_pin()

    def dynamic_bearing_force(self):
        L=Symbol('L')#wymagana trwałość
        return self.max_dynamic_force_pin() * L**(S.One/3)
    
    def dynamic_force_pin_diameter(self):
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)
        
        
        load = abs(self.max_dynamic_force_pin())
        
        return ((4*load)/(pi*kt*Re))**(1/2)
    

    def get_default_data(self):

        #m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0', positive=True)
        c0, k_0, l_l0, omega, F_0 = symbols('c_0 k_0 l_l0 Omega F_0', positive=True)
        m0, l0  = self.m0,self.l0
        c0, k0= self.c0,self.k0
        lam0=self.lam0
        Omega0=self.Omega0
        F0=self.F0
        Omega = self.Omega

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
           
#             self.Omega: [self.Omega],
            #self.F_engine: [4*no*F0*S.One*abs(cos(pi/2*no))  + F0*S.One*no*cos(self.Omega*self.ivar) for no in range(1,8)],
            self.F_engine: [4 * no * F0 * S.One for no in range(1,8)],
        }

        return default_data_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST,
        mech_comp.MaxDynamicForce,
        mech_comp.MaxStaticForce,
        mech_comp.DynamicPinDiameter,
        mech_comp.StaticPinDiameter,
        mech_comp.SpringForce
        ]
        
        return comp_list
    
class UndampedSymmetricalVehicleSuspension(UndampedVehicleSuspension):

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
    ivar=Symbol('t')
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
                 l_l=None,
                 F_engine=None,
                 Omega=None,
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
        if  Omega is not None: self.Omega = Omega
        if  phi is not None: self.phi = phi
        if z is not None: self.z=z
        # self.z, self.phi = self.qs
        self.qs = [self.z,self.phi]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=self.qs)(label='rod')
        self._left_mount = Spring(self.k_l, pos1=self.z + self.phi * self.l_l, qs=self.qs)(label='left spring')
        self._right_mount = Spring(self.k_r, pos1=self.z - self.phi * self.l_r, qs=self.qs)(label='right spring')
        self._force = Force(self.F_engine*cos(self.Omega*self.ivar), pos1=self.z - self.l_r * self.phi, qs=self.qs)(label='force')
        self._static_force = Force(self.F_engine, pos1=self.z - self.l_r * self.phi, qs=self.qs)(label='force')

        components['_body'] = self._body
        components['_left_mount'] = self._left_mount
        components['_right_monut'] = self._right_mount
        components['_force'] = self._force
        components['_static_force'] = self._static_force

        return components

    def get_default_data(self):

        #m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0', positive=True)
        c0, k_0, l_l0, omega, F_0 = symbols('c_0 k_0 l_l0 Omega F_0', positive=True)
        m0, l0  = self.m0,self.l0
        c0, k0= self.c0,self.k0
        lam0=self.lam0
        Omega0=self.Omega0
        F0=self.F0
        Omega = self.Omega

        default_data_dict = {
            self.l_l: [self.l_r],
            self.l_rod:[l0*S.One*no for no in range(1, 8)],
            self.I: [S.One/12*self.m *self.l_rod**2],
            self.m: [m0*S.One*no for no in range(1,50)],
            self.k_l: [self.k_r],
            self.k_r: [k0*S.One*no for no in range(1,40)],
            self.l_r: [l0*S.One*no for no in range(1, 8)],
#           self.Omega: [self.Omega],
            self.F_engine: [4 * no * F0 * S.One for no in range(1,8)],
        }

        return default_data_dict
    
    def reference_data_dict(self):

        default_data_dict = {
            self.l_l: [self.l_r],
#             self.I: [S.One/12*self.m *self.l_rod**2],
            self.k_l: [self.k_r],
        }

        return default_data_dict

    
class DampedVehicleSuspension(UndampedVehicleSuspension):

    scheme_name = 'damped_car.png'
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

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components =super().components

        
        self._damper_l = Damper(self.c_l, pos1=self.z + self.phi * self.l_l,
                               qs=self.qs)(label='left damper')
        self._damper_r = Damper(c=self.c_r, pos1=self.z - self.phi * self.l_r,
                               qs=self.qs)(label='right damper')



        components['_damper_l'] = self._damper_l
        components['_damper_r'] = self._damper_r
       

        
        return components

#     def max_static_force_pin(self):
#         data=self._given_data
#         ans=self.dynamic_force()
#         free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
#         return abs(free_coeff)
    
#     def static_force_pin_diameter(self):
#         kt=Symbol('k_t', positive=True)
#         Re=Symbol('R_e', positive=True)
#         return ((4*self.max_static_force_pin())/(pi*kt*Re))**(1/2)
    
#     def dynamic_force(self):
        
#         amps = self._fodes_system.steady_solution.as_dict()
#         force = self.components['_spring_l'].force().doit().expand()#.doit()
#         data=self._given_data
        
#         #display(abs(self.components['_spring_l'].force()))
#         return (self.components['_spring_l'].force().subs(amps)).subs(data).expand().doit()

#     def max_dynamic_force_pin(self):
        
#         ans = self.dynamic_force()
        
#         #display(abs(self.components['_spring_l'].force()))
#         data=self._given_data
#         sin_coeff=ans.coeff(sin(self.Omega*self.ivar)).subs(data)
#         #display(sincoeff)
#         cos_coeff=ans.coeff(cos(self.Omega*self.ivar)).subs(data)
#         #display(coscoeff)
#         free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
#         #display(freecoeff)

#         return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)

    
#    def dynamic_bearing_force(self):
       
#       return self.max_dynamic_force_pin() * self.sqrt(L)
    

#     def max_dynamic_force_pin(self):
#         return (self._frf()[0].subs(self.lam0,0) + self._frf()[1].subs(self.lam0,0)*self.l_rod/2)*self.stiffness_matrix()[0]

    def dynamic_bearing_force(self):
        L=Symbol('L')#wymagana trwałość
        return self.max_dynamic_force_pin() * L**(S.One/3)

    
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
            #self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.k_l: [k0*S.One*no for no in range(1,8)],
            self.k_r: [k0*S.One*no for no in range(1,8)],
            
            
            self.c_l: [self.lam*(self.k_l)],
            self.l_r: [l0*S.One*no for no in range(1, 8)],
            self.l_l: [l0*S.One*no for no in range(1, 8)],
           
#             self.Omega: [Omega0*S.One*no for no in range(1,2)],
            self.F_engine: [F0*S.One*no for no in range(1,20)],}
            


        return default_data_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST,
        mech_comp.DampedVibrationFrequency,
        mech_comp.MaxDynamicForce,
        mech_comp.DynamicForceComponent,
        mech_comp.MaxStaticForce,
        mech_comp.DynamicPinDiameter,
        mech_comp.DampingMatrixComponent
        ]
        
        return comp_list

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
            #self.lam: [lam0/10*S.One*no for no in range(1,8)],
            self.k_l: [self.k_r],
            self.k_r: [k0*S.One*no for no in range(1,8)],
            
            
            self.c_l: [self.lam*(self.k_l)],
            self.l_r: [self.l_l],
            self.l_l: [l0*S.One*no for no in range(1, 8)],
           
#             self.Omega: [Omega0*S.One*no for no in range(1,2)],
            self.F_engine: [F0 * S.One * no for no in range(1,20)],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        parameters_dict[self.k_l] = parameters_dict[self.k_r]
        parameters_dict[self.l_r] = parameters_dict[self.l_l]
#        parameters_dict[self.m_2] = parameters_dict[self.m]

        return parameters_dict

    def reference_data_dict(self):

        lam, m0, k0, l0, F0 = symbols('lambda m_0 k_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.c_r: [self.lam*(self.k_r)],
            self.c_l: [self.lam*(self.k_l)],
            self.l_l: [self.l_r],
            self.k_l: [self.k_r],
        }

        return default_data_dict

class SimplifySuspension(UndampedVehicleSuspension):
    """
    Ready to use sample Double Degree of Freedom System represents symmetrical kinematically excited beam with two springs.
        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            I = Moment of inertia
                -Moment of inertia of a rod

            l_rod = lenght of a rod
                -Dimension of a rod imitating vehicle's floor plate

            l_r = offset of right spring
                -Dimension of a right spring's offset

            l_l = offset of left spring
                -Dimension of a left spring's offset

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            F_engine = thrust of an engine
                -Force made by engine, which allows vehicle to move forward

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, k_1, k_2, l_rod, l_l, l_r, F_engine = symbols('m, k_1, k_2, l_rod, l_l, l_r, F_{engine}')
        >>> qs = dynamicsymbols('z, varphi') 
        >>> DDoFSymplifyVehicleSuspension()

        -We define the symbols and dynamicsymbols
        -and then we determine the instance of the system using class DDoFSymplifyVehicleSuspension()
    """
    scheme_name = 'vehicle_suspension.PNG'
    real_name = 'vehicle_suspension_real.PNG'

    m=Symbol('m', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_r=Symbol('l_r', positive=True)
    k_1=Symbol('k_1', positive=True)
    k_2=Symbol('k_2', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    Omega=Symbol('Omega', positive=True)
    ivar=Symbol('t')
    phi=dynamicsymbols('\\varphi')
    z=dynamicsymbols('z')

    def __init__(self,
                 m=None,
                 I=None,
                 l_rod=None,
                 l_r=None,
                 k_2=None,
                 k_1=None,
                 l_l=None,
                 F_engine=None,
                 phi=None,
                 z=None,
                 Omega=None,
                 ivar=Symbol('t'),
#                  qs=None,
                 **kwargs):

        if m is not None: self.m =m # mass of a rod
        if I is not None: self.I = I # moment of inertia of a rod
        if l_rod is not None: self.l_rod =l_rod# length of a rod
        if l_r is not None: self.l_r = l_r 
        if l_l is not None: self.l_l = l_l 
        if k_2 is not None: self.k_2 =k_2# left spring
        if k_1 is not None: self.k_1 = k_1# right spring
        if F_engine is not None: self.F_engine = F_engine
        if phi is not None: self.phi = phi
        if Omega is not None: self.Omega = Omega
        if z is not None: self.z=z
        self.ivar=ivar
        self.qs = [self.z,self.phi]

        self._init_from_components(**kwargs)


    @cached_property
    def components(self):

        components = {}

#         self._body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=self.qs)(label='rod')
#         self._spring_l = Spring(self.k_1, pos1=self.z + self.phi * self.l_l, qs=self.qs)
#         self._spring_2 = Spring(self.k_2, pos1=self.z - self.phi * self.l_r, qs=self.qs)(label='right spring')
#         self._force = Force(self.F_engine*cos(self.Omega*self.ivar), pos1=self.z - self.phi * self.l_r, qs=self.qs)(label='force')

#         components['_body'] = self._body
#         components['_spring_l'] = self._spring_l
#         components['_spring_2'] = self._spring_2
#         components['_force'] = self._force

        self._vehicle = UndampedVehicleSuspension(m=self.m, I=self.I, l_rod=self.l_rod, l_r=self.l_r, k_r=self.k_2, k_l=self.k_1, l_l=self.l_l, F_engine=self.F_engine, Omega=self.Omega, phi=self.phi, z=self.z, ivar=Symbol('t'), qs=None)

        components['_vehicle'] = self._vehicle

        return components

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)
        e,m_e = symbols('e m_e',positive=True)

        default_data_dict = {
            self.I: [(self.m / 12) * (2 * self.l_rod)**2],
            self.m: [m0 * S.One * no for no in range(10,21)],
            self.k_1: [k0 * S.One * no for no in range(100,201)],
            self.k_2: [k0 * S.One * no for no in range(100,201)],
            self.l_l: [l0 * S.One * no for no in range(1,11)],
            self.l_r: [l0 * S.One * no for no in range(1,11)],
            self.l_rod: [l0 * S.One * no for no in range(3,7)],
            self.F_engine: [e*m_e*self.Omega**2],
            m_e: [m0 * S.One * no for no in range(10,11)],
            e: [l0 * S.One * no for no in range(1,6)],
        }

        return default_data_dict
    
class QuarterOfVehicle(ComposedSystem):
    
    k_d=Symbol('k_d', positive=True)
    k_u=Symbol('k_u', positive=True)
    m_d=Symbol('m_d', positive=True)
    m_u=Symbol('m_u', positive=True)
    c=Symbol('c', positive=True)
    omega=Symbol('omega', positive=True)
    A_0=Symbol('A_0', positive=True)
    x_1=dynamicsymbols('x_1')
    x_2=dynamicsymbols('x_2')
    ivar=Symbol('t')
    
    def __init__(self,
                 k_d=None,
                 k_u=None,
                 m_d=None,
                 m_u=None,
                 c=None,
                 omega=None,
                 A_0=None,
                 x_1=None,
                 x_2=None,
                 ivar=Symbol('t'),
                 **kwargs):
    
    
        if m_d is not None: self.m_d = m_d
        if k_d is not None: self.k_d = k_d
        if m_u is not None: self.m_u = m_u
        if k_u is not None:self.k_u = k_u
        if c is not None: self.c=c
        if omega is not None: self.omega=omega
        if A_0 is not None: self.A_0=A_0
        if x_1 is not None: self.x_1 = x_1
        if x_2 is not None : self.x_2 = x_2

        self.ivar = ivar
        self.qs = [self.x_1, self.x_2]

        self._init_from_components(**kwargs)
        
    @cached_property
    def components(self):
        
        components = {}
        
        self._material_point_1 = MaterialPoint(self.m_d, self.x_1, qs=self.qs)
        self._spring_1 = Spring(self.k_d,self.x_1,self.A_0*sin(self.ivar*self.omega),qs=self.qs)
        self._material_point_2 = MaterialPoint(self.m_u, self.x_2, qs=self.qs)
        self._spring_2 = Spring(self.k_u,self.x_1,self.x_2,qs=self.qs)
        self._damper_1=Damper(self.c,self.x_1,self.x_2, qs=self.qs)

        components['material_point_1'] = self._material_point_1
        components['spring_1'] = self._spring_1
        components['material_point_2'] = self._material_point_2
        components['spring_2'] = self._spring_2
        components['damper'] = self._damper_1
        
        return components
       
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

    