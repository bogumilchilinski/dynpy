from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, Element
from  ..continuous import ContinuousSystem, PlaneStressProblem

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST


from sympy.physics import units
ureg = units

#TODO
#ISSUE #149
class VesselWithHeaveAndRotation(ComposedSystem):

    scheme_name = 'vessel.jpg'
    
    m_vessel=Symbol('M_vessel', positive=True)
    
    I_5=Symbol('I_5', positive=True)
    
    wave_level=dynamicsymbols('W')
    wave_slope=dynamicsymbols('S')
    rho=Symbol('rho', positive=True)
    g=Symbol('g', positive=True)
    A_wl=Symbol('A_wl',positive=True)
    V=Symbol('V', positive=True)
    GM_L=Symbol('GM_L', positive=True)
    CoB=Symbol('CoB', positive=True)
    CoF=Symbol('CoF', positive=True)
    A_h=Symbol('A_h',positive=True)
    Phi_h=Symbol('Phi_h',positive=True)
    omega=Symbol('omega',positive=True)
    
    H = dynamicsymbols('H')
    Phi = dynamicsymbols('Phi')
    qs=dynamicsymbols('H, Phi')
    ivar=Symbol('t')
    
    def __init__(self,
                 m_vessel=None,
                 I_5=None,

                 wave_level=None,
                 wave_slope=None,
                 rho=None,
                 g=None,
                 A_wl=None,
                 V=None,
                 GM_L=None,
                 CoB=None,
                 CoF=None,
                 A_h=None,
                 Phi_h=None,
                 omega=None,
                 qs=None,
                 H=None,
                 Phi=None,
                 
                 ivar=Symbol('t'),
                 **kwargs):

        if m_vessel is not None:self.m_vessel = m_vessel
        if I_5 is not None:self.I_5 = I_5
        if wave_level is not None:self.wave_level = wave_level
        if wave_slope is not None:self.wave_slope = wave_slope
        if rho is not None:self.rho = rho
        if g is not None:self.g = g
        if A_wl is not None:self.A_wl = A_wl
        if V is not None:self.V = V
        if GM_L is not None:self.GM_L = GM_L
        if CoB is not None:self.CoB = CoB
        if CoF is not None:self.CoF = CoF
        if A_h is not None:self.A_h = A_h
        if Phi_h is not None:self.Phi_h = Phi_h
        if omega is not None:self.omega = omega
        if H is not None:self.H = H
        if Phi is not None:self.H = Phi
        
        self.ivar = ivar

        
        
        self.M_matrix = Matrix([[self.m_vessel, 0], [0, self.I_5]]) #wszedzie teraz pobiera None # dopisz wszędzie self tut
        self.K_matrix = Matrix([[self.rho * self.g * self.A_wl, -self.rho * self.g * self.A_wl * (self.CoF - self.CoB)],#trzeba ify popsać? #ify jużsą tylko na to H i Phi brakue chyba
                           [-self.rho * self.g * self.A_wl * (self.CoF - self.CoB),
                            self.rho * self.g * self.V * self.GM_L]])
        
        if qs is not None: self.qs = qs
        self.qs = [self.H, self.Phi]
        

        
        # vessel mass and stiffness matrix
        #M_matrix = Matrix([[m_vessel, 0], [0, I_5]])

        #K_matrix = Matrix([[rho * g * A_wl, -rho * g * A_wl * (CoF - CoB)],
                           #[-rho * g * A_wl * (CoF - CoB),
                            #rho * g * V * GM_L]])

        # generalized displacements and velocities
        
        #wave_level = A_h*cos(omega*self.ivar + Phi_h)
        #wave_slope = diff(wave_level,self.ivar)*(omega/g)
        
        #dh_wave_slope = diff(wave_level,t,t)*(omega/g)   # TU JEST TEN MAGICZNY FORMUŁ, KTÓRY CHYBA JEDNAK POWINNIŚMY WPROWADZIĆ BO dq go tak nie policzy XD #
        
        q = Matrix(self.qs) + Matrix([self.wave_level, self.wave_slope])
        dq = q.diff(self.ivar)
        
        

        # lagrangian components definition
        
        self.kinetic_energy = S.Half * sum(dq.T * self.M_matrix * dq)
        self.potential_energy = S.Half * sum(Matrix(self.qs).T * self.K_matrix * Matrix(self.qs))
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}

        self._vessel_inertia = Element(self.kinetic_energy,
                                 qs=self.qs)
        self._buoyancy_stiffness = Element(-self.potential_energy,
                                 qs=self.qs)

        components['_vessel_inertia']=self._vessel_inertia
        components['_buoyancy_stiffnes']=self._buoyancy_stiffness

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_vessel: r'mass of vessel,',
            self.I_5: r'moment of inertia of 5-th degree (with respect to y axis, determined by the radius of gyration),',
            self.q: r'vessel generalized coordinates,',
            self.wave_level: r'wave level,',
            self.wave_slope: r'wave slope,',
            self.rho: r'fluid density,',
            self.g: r'acceleration of gravity,',
            self.A_wl: r'wetted area,',
            self.V: r'submerged volume of the vessel,',
            self.GM_L: r'longitudinal metacentric height,',
            self.CoB: r'centre of buoyancy,',
            self.CoF: r'centre of floatation.',

        }

        return self.sym_desc_dict
    
    def numerical_data(self):
        numbers_dict={
            self.I_5:21*self.m_vessel,
            self.m_vessel: 1e7,
            self.wave_level: 0.5 * cos(2/7 * pi * self.ivar),
            self.wave_slope: 0.5 * cos(2/7 * pi * self.ivar),
            self.rho:1025,
            self.g: 9.81,
            self.A_wl:4025,
            self.V:25000,
            self.GM_L:290,
            self.CoB:3.6,
            self.CoF:13,
            self.A_h:0.45,
            self.Phi_h:-2.95,
        }
        
        return numbers_dict
    
    def get_numerical_parameters(self):

        m_vessel, I_5, g, rho, A_wl, V, GM_L, CoB, CoF, A_h, Phi_h, wave_level, wave_slope = self.m_vessel, self.I_5, self.g, self.rho, self.A_wl, self.V, self.GM_L, self.CoB, self.CoF, self.A_h, self.Phi_h, self.wave_level, self.wave_slope
        
        get_numerical_parameters = {

            self.I_5:21*self.m_vessel,
            self.m_vessel: 1e7,
            self.wave_level: 0.5 * cos(2/7 * pi * self.ivar),
            self.wave_slope: 0.5 * cos(2/7 * pi * self.ivar),
            self.rho:1025,
            self.g: 9.81,
            self.A_wl:4025,
            self.V:25000,
            self.GM_L:290,
            self.CoB:3.6,
            self.CoF:13,
            self.A_h:0.45,
            self.Phi_h:-2.95, 
        }
        return get_numerical_parameters
    
    def units(self):
        units_dict={

            self.q[1]:ureg.radian,
            self.q[0]:ureg.meter,
            self.q[1].diff(self.ivar):ureg.radian/ureg.second,
            self.q[0].diff(self.ivar):ureg.meter/ureg.second,

           }
        
        return units_dict    
    
#TODO
#ISSUE #150
class TDoFCompensatedPayload(ComposedSystem):

    scheme_name = '3dofs_new.PNG'

    def __init__(self,
                 phi=dynamicsymbols('varphi'),
                 h=dynamicsymbols('h'),
                 h_c=dynamicsymbols('h_c'),
                 m_p=Symbol('m_p', positive=True),
                 k_w=Symbol('k_w', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 qs=dynamicsymbols('varphi h h_c'),
                 y_e=dynamicsymbols('y_e'),
                 z_e=dynamicsymbols('z_e'),
                 m_c=Symbol('m_c', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 l_c=Symbol('l_c', positive=True),
                 l_w=Symbol('l_w', positive=True),
                 g=Symbol('g', positive=True),
                 h_eq=Symbol('h_eq', positive=True),
                 h_ceq=Symbol('h_ceq', positive=True),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m_p = m_p
        self.k_w = k_w
        self.l_0 = l_0
        self.y_e = y_e
        self.z_e = z_e
        self.m_c = m_c
        self.k_c = k_c
        self.l_c = l_c
        self.l_w = l_w
        self.g = g
        self.h_eq = h_eq
        self.h_ceq = h_ceq
        self.ivar = ivar
        self.phi = phi
        self.h = h
        self.h_c = h_c
        self.qs = qs

        phi, h, h_c = qs

        y = (h + h_eq + l_0 + l_c) * sin(phi) + y_e
        z = (h + h_eq + l_0 + l_c) * cos(phi) + z_e

        y_c = (h_c + h_ceq + l_0) * sin(phi) + y_e
        z_c = (h_c + h_ceq + l_0) * cos(phi) + z_e

        y, z, y_c, z_c = dynamicsymbols('y_{p},z_{p},y_c,z_c')

        v_c = sqrt(diff(y_c, ivar)**2 + diff(z_c, ivar)**2)
        v = sqrt(diff(y, ivar)**2 + diff(z, ivar)**2)

        positions_dict = {
            y: (h + h_eq + l_0 + l_c) * sin(phi) + y_e,
            z: (h + h_eq + l_0 + l_c) * cos(phi) + z_e,
            y_c: (h_c + h_ceq + l_0) * sin(phi) + y_e,
            z_c: (h_c + h_ceq + l_0) * cos(phi) + z_e,
        }

        self.kinetic_energy = S.Half * m_p * v**2 + S.Half * m_c * v_c**2

        self.potential_energy = (S.Half * k_w * (h_c + h_ceq)**2 +
                                 S.Half* k_c * (h + h_eq - (h_c + h_ceq))**2 -
                                 m_p * g * z - m_c * g * z_c)

        super(HarmonicOscillator,self).__init__((self.kinetic_energy - self.potential_energy).subs(positions_dict).doit(),
                         qs=qs,
                         ivar=ivar,**kwargs)


    def single_dof_transformation(self,p=None):
        
        if p is None:
            p = S.One
        
        
        omega, delta, eps = symbols('omega, delta, varepsilon',positive=True)
        A_y, A_z, Phi_y, Phi_z = symbols('A_y, A_z, Phi_y, Phi_z,',positive=True)
        
        sdof_payload = self(self.phi)
        
        nonlin_mathieu_model = sdof_payload.subs({self.h:0,self.h_c:0,self.m_c:0}).subs({self.h_eq:self.l_c,self.h_ceq:self.l_c,self.l_0:self.l_w-self.l_c})

        mathieu_nonlin_approx = nonlin_mathieu_model.approximated(n=3)

        mathieu_acc_form = mathieu_nonlin_approx._to_acc()

        mathieu_basis = mathieu_acc_form.subs({self.l_c:0,
                                               self.g:omega**2*self.l_w})

        mathieu_instability_eq = mathieu_basis.subs(S.Half*omega**2*self.phi**2,(S.Half)**3*omega**2*self.phi**2*(1+eps*delta),simultenouos=True)

        mathieu_transformed_eq = mathieu_instability_eq.subs({self.l_w:1/eps,
                                                              self.phi.diff(self.ivar)*self.phi**2*self.y_e.diff(self.ivar):0,
                                                              self.z_e:A_z*cos(omega/p*self.ivar+Phi_z),
                                                              self.y_e:A_y*cos(omega/p*self.ivar+Phi_y),
                                                              A_y*self.phi**2*eps:0
                                                             })

        return mathieu_transformed_eq._eoms

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m_p: r'mass of payload \si{[\kilogram]}',
#             self.k_w: r'wire stiffness \si{[\newton\per\meter]}',
#             self.l_0: r'length of the lifting cable \si{[\metre]}',
#             tuple(self.q): r'generalized coordinates',
#             self.y_e:
#             r'lateral displacement at crane tip obtained from RAOs (a regular wave excitation) \si{[\metre]}',
#             self.z_e:
#             r'vertical displacement at crane tip obtained from RAOs (a regular wave excitation) \si{[\metre]}',
#             self.m_c: r'mass of compensator \si{[\kilogram]}',
#             self.k_c:
#             r'stiffness of heave compensator \si{[\newton\per\meter]}',
#             self.l_c:
#             r'length of the attached compensating element \si{[\metre]}',
#             self.g: r'acceleration of gravity \si{[\metre/\second\squared]}',
#             self.h_eq: r'equilibrium point of payload \si{[\metre]}',
#             self.h_ceq: r'equilibrium point of compensator \si{[\metre]}',
#             self.ivar: r'independent time variable',
#         }

#         return self.sym_desc_dict

    def symbols_description(self):

        parent_symbols_dict=super().symbols_description()

        self.sym_desc_dict = {** parent_symbols_dict , **{
            self.m_p: r'mass of payload,',
            self.k_w: r'wire stiffness,',
            self.l_0: r'length of the lifting cable,',
            self.l_w: r'length of the lifting cable-payload system,',
            tuple(self.q): r'payload generalized coordinates,',
            self.y_e:
            r'lateral displacement at crane tip obtained from RAOs (a regular wave excitation),',
            self.z_e:
            r'vertical displacement at crane tip obtained from RAOs (a regular wave excitation),',
            self.m_c: r'mass of compensator,',
            self.k_c:
            r'stiffness of heave compensator,',
            self.l_c:
            r'length of the attached compensating element,',
            self.g: r'acceleration of gravity,',
            self.h_eq: r'equilibrium point of payload,',
            self.h_ceq: r'equilibrium point of compensator,'
        }}

        return self.sym_desc_dict

    def numerical_data(self):
        numbers_dict={
            self.h_eq: ((self.m_c+self.m_p)/self.k_w+(self.m_p)/self.k_c)*self.g,
            self.h_ceq: self.g*(self.m_c+self.m_p)/self.k_w,
            self.m_p: 1e5,
            self.k_w: 4.09e6,
            self.l_0: self.l_w - self.l_c,
            self.y_e: 0.5 * cos(2/7 * pi * self.ivar -2.396),
            self.z_e: 1.5 * cos(2/7 * pi * self.ivar + 0.315),
            self.m_c: 1e4,
            self.k_c: 2.5e6,
            self.l_c:5,
            self.l_w:48,
            self.g: 9.81,
        }
        
        return numbers_dict    
    

    
    def units(self):
        units_dict={
            self.q[2]:ureg.meter,
            self.q[1]:ureg.meter,
            self.q[0]:ureg.radian,
            self.q[2].diff(self.ivar):ureg.meter/ureg.second,            
            self.q[1].diff(self.ivar):ureg.meter/ureg.second,
            self.q[0].diff(self.ivar):ureg.radian/ureg.second,

           }
        
        return units_dict 
# class PayloadVesselSystem(ComposedSystem):
    
#     def __init__(self,
#                  y_e=dynamicsymbols('y_e'),
#                  z_e=dynamicsymbols('z_e'),
#                  wave_level=dynamicsymbols('W'),
#                  wave_slope=dynamicsymbols('S'),
# #                  payload=TDoFCompensatedPayload(),
# #                  vessel=DDoFVessel(),
#                  system=None
#                 ):
        
#         self.payload = TDoFCompensatedPayload(y_e,z_e)
#         self.vessel = DDofVessel(wave_level,wave_slope)
        
#         system = self.payload + self.vessel
        
#         super().__init__(system)





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
    
class ROVImpactDamped(ComposedSystem):

    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    m_n=Symbol('m_n', positive=True)
    c=Symbol('c',positive=True)
    k=Symbol('k', positive=True)
    ivar=Symbol('t')
    
    z, z_s=dynamicsymbols('z z_s')
    
    def __init__(self,
                 m=None,
                 c=None,
                 z=None,
                 z_s=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if c is not None: self.c = c
        if z is not None: self.z = z
        if z_s is not None: self.z_s = z_s

        self.ivar = ivar
        self.qs = [self.z,self.z_s]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.damper = Damper(self.c, pos1=self.z_s, qs=self.qs)

        components['material_point'] = self.material_point
        components['damper'] = self.damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict

class ROVImpactElasticDamped(ROVImpactDamped):

    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    m_n=Symbol('m_n', positive=True)
    c=Symbol('c',positive=True)
    k=Symbol('k', positive=True)
    ivar=Symbol('t')
    
    z, z_s=dynamicsymbols('z z_s')
    
    def __init__(self,
                 m=None,
                 c=None,
                 k=None,
                 z=None,
                 z_s=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k = k
        if c is not None: self.c = c
        if z is not None: self.z = z
        if z_s is not None: self.z_s = z_s

        self.ivar = ivar
        self.qs = [self.z,self.z_s]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.damper = Damper(self.c, pos1=self.z_s, qs=self.qs)
        self.spring = Spring(self.k, pos1=self.z, pos2=self.z_s, qs=self.qs)

        components['material_point'] = self.material_point
        components['damper'] = self.damper
        components['spring'] = self.spring

        return components
    
class ROVandStructureImpact(ROVImpactDamped):

    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    m_n=Symbol('m_n', positive=True)
    c=Symbol('c',positive=True)
    k=Symbol('k', positive=True)
    ivar=Symbol('t')
    
    z, z_s=dynamicsymbols('z z_s')
    
    def __init__(self,
                 m=None,
                 c=None,
                 k=None,
                 z=None,
                 z_s=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if k is not None: self.k = k
        if c is not None: self.c = c
        if z is not None: self.z = z
        if z_s is not None: self.z_s = z_s

        self.ivar = ivar
        self.qs = [self.z,self.z_s]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring = Spring(self.k, pos1=self.z, pos2=self.z_s, qs=self.qs)
        self.damper = Damper(self.c, pos1=self.z_s, qs=self.qs)
        self.material_point_aux = MaterialPoint(self.m_n, self.z_s, qs=self.qs)

        components['material_point'] = self.material_point
        components['spring'] = self.spring
        components['damper'] = self.damper
        components['material_point1'] = self.material_point_aux

        return components

class ROVandStructureImpactMassModel(ROVandStructureImpact):

    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    m=Symbol('m', positive=True)
    m_n=Symbol('m_n', positive=True)
    k_n=Symbol('k_n', positive=True)
    c=Symbol('c', positive=True)
    k=Symbol('k', positive=True)
    ivar=Symbol('t')
    
    z, z_s=dynamicsymbols('z z_s')
    
    def __init__(self,
                 m=None,
                 m_n=None,
                 c=None,
                 k=None,
                 k_n=None,
                 z=None,
                 z_s=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if m_n is not None: self.m_n = m_n
        if k is not None: self.k = k
        if k_n is not None: self.k_n = k_n
        if c is not None: self.c = c
        if z is not None: self.z = z
        if z_s is not None: self.z_s = z_s

        self.ivar = ivar
        self.qs = [self.z,self.z_s]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.material_point = MaterialPoint(self.m, self.z, qs=self.qs)
        self.spring_rov = Spring(self.k, pos1=self.z, pos2=self.z_s, qs=self.qs)
        self.damper = Damper(self.c, pos1=self.z_s, qs=self.qs)
        self.material_point_aux = MaterialPoint(self.m_n, self.z_s, qs=self.qs)
        self.spring_str = Spring(self.k_n, pos1=self.z_s, pos2=0, qs=self.qs)

        components['material_point'] = self.material_point
        components['spring_rov'] = self.spring_rov
        components['damper'] = self.damper
        components['material_point1'] = self.material_point_aux
        components['spring_str'] = self.spring_str

        return components
    