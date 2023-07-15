from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, EngineMount
from ..continuous import ContinuousSystem, PlaneStressProblem
from dynpy.models.mechanics.tmd import TunedMassDamper

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin

#####

class SDOFWinchSystem(ComposedSystem):


    scheme_name = 'winch.png'
    real_name = 'winch_cross_sec.png'
    I_s=Symbol('I_s', positive=True)
    I_k=Symbol('I_k', positive=True)
    I_1=Symbol('I_1', positive=True)
    i_1=Symbol('i_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    i_2=Symbol('i_2', positive=True)
    I_3=Symbol('I_3', positive=True)
    I_4=Symbol('I_4', positive=True)
    I_b=Symbol('I_b', positive=True)
    A=Symbol('A', positive=True)
    B=Symbol('B', positive=True)
    ivar=Symbol('t')
    phi=dynamicsymbols('\\varphi')
    dphi=dynamicsymbols('\\varphi',1)
    D=Symbol('D', positive=True)
    mu=Symbol('\\mu', positive=True)
    M_s=Symbol('M_s', positive=True)
    M_T=Symbol('M_T', positive=True)
    G=Symbol('G', positive=True)
    g=Symbol('g',positive=True)
    alpha=Symbol('\\alpha')
    pos1=phi
    qs=[phi]
    phi_1=phi
    phi_2=phi_1*i_1
    phi_3=phi_2*i_2
    def __init__(self,
                 I_s=None,
                 I_k=None,
                 I_1=None,
                 i_1=None,
                 I_2=None,
                 i_2=None,
                 I_3=None,
                 I_4=None,
                 I_b=None,
                 A=None,
                 B=None,
                 ivar=Symbol('t'),
                 phi=None,
                 dphi=None,
                 D=None,
                 mu=None,
                 M_s=None,
                 M_T=None,
                 G=None,
                 g=None,
                 alpha=None,
                 phi_1=None,
                 phi_2=None,
                 phi_3=None,
                 qs=None,
                 **kwargs):
        

        
        if I_s is not None: self.I_s = I_s
        if I_k is not None: self.I_k = I_k
        if I_1 is not None: self.I_1 = I_1
        if i_1 is not None: self.i_1 = i_1
        if I_2 is not None: self.I_2 = I_2
        if I_3 is not None: self.I_3 = I_3
        if i_2 is not None: self.i_2 = i_2
        if I_4 is not None: self.I_4 = I_4
        if I_b is not None: self.I_b = I_b
        if D is not None: self.D = D
        if G is not None: self.G = G
        if g is not None: self.g = g
        if mu is not None: self.mu = mu
        if M_T is not None: self.M_T = M_T
        if M_s is not None: self.M_s = M_s
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if phi_3 is not None: self.phi_3 = phi_3
        if alpha is not None: self.alpha = alpha
        if A is not None: self.A=A
        if B is not None: self.B=B
        if dphi is not None: self.dphi = dphi
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if phi_3 is not None: self.phi_3 = phi_3
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._engine = Disk(self.I_s, self.phi_1, qs=self.qs)
#         self.disk_k = Disk(I_k, phi_1, qs=qs)
        self._disk_1 = Disk(self.I_1, self.phi_1, qs=self.qs)
        self._disk_2 = Disk(self.I_2, self.phi_2, qs=self.qs)
        self._disk_3 = Disk(self.I_3, self.phi_2, qs=self.qs)
        self._disk_4 = Disk(self.I_4, self.phi_3, qs=self.qs)
        self._disk_B = Disk(self.I_b, self.phi_3, qs=self.qs)
        self._load = MaterialPoint(self.G/self.g, pos1=self.phi_3* self.D/2, qs=self.qs)
        self._friction_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * cos(self.alpha)*self.mu , qs=self.qs)
        self._force = Force(-self.M_T,pos1=self.phi_3,qs=self.qs)
        self._drive = Force(self.M_s,pos1=self.phi_1,qs=self.qs)
        self._gravity_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * sin(self.alpha) , qs=self.qs)


        components['engine'] = self._engine
        components['disk_1'] = self._disk_1
        components['disk_2'] = self._disk_2
        components['disk_3'] = self._disk_3
        components['disk_4'] = self._disk_4
        components['disk_B'] = self._disk_B
        components['load'] = self._load
        components['friction_comp'] = self._friction_comp
        components['force'] = self._force
        components['drive'] = self._drive
        components['gravity_comp'] = self._gravity_comp
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
            self.I_s: [20,30,40],
            self.I_k: [20,30,40],
            self.I_1: [20,30,40],
            self.I_2: [20,30,40],
            self.I_3: [20,30,40],
            self.I_4: [20,30,40],
            self.I_b: [20,30,40],
            self.i_1: [20,30,40],
            self.i_2: [20,30,40],
            self.D: [10,20,30],
            self.G: [10,20,30],
            self.g: [10,20,30],
            self.mu : [10,20,30],
            self.M_T : [10,20,30],
            self.M_s : [10,20,30],
            self.phi_1: [10,20,30],
            self.phi_2: [10,20,30],
            self.phi_3  : [10,20,30],
            self.alpha  : [10,20,30],
        }

        return default_data_dict
        
    def steady_angular_velocity(self):
        obj=self
        eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)
        return omega_steady_state
    
    def reduced_torque(self):
        obj=self
        
        red_tor = -(obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())

        return delta_1

    def reduced_inertia(self):
        obj = self
        
        ans = obj.inertia_matrix()[0]

        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.02)-1)*I_r).n(5))
        return ik
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_3
        return st_val
    

    
class SDOFWinchSystemTest(ComposedSystem):


    scheme_name = 'winch_scheme_test.png'
    real_name = 'winch_cross_sec.png'
    I_s=Symbol('I_s', positive=True)
    I_k=Symbol('I_k', positive=True)
    I_1=Symbol('I_1', positive=True)
    i_1=Symbol('i_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    I_b=Symbol('I_b', positive=True)
    A=Symbol('A', positive=True)
    B=Symbol('B', positive=True)
    ivar=Symbol('t')
    phi=dynamicsymbols('\\varphi')
    dphi=dynamicsymbols('\\varphi',1)
    D=Symbol('D', positive=True)
    mu=Symbol('\\mu', positive=True)
    M_s=Symbol('M_s', positive=True)
    M_T=Symbol('M_T', positive=True)
    G=Symbol('G', positive=True)
    g=Symbol('g',positive=True)
    alpha=Symbol('\\alpha')
    pos1=phi
    qs=[phi]
    phi_1=phi
    phi_2=phi_1*i_1

    def __init__(self,
                 I_s=None,
                 I_k=None,
                 I_1=None,
                 i_1=None,
                 I_2=None,
                 I_b=None,
                 A=None,
                 B=None,
                 ivar=Symbol('t'),
                 phi=None,
                 dphi=None,
                 D=None,
                 mu=None,
                 M_s=None,
                 M_T=None,
                 G=None,
                 g=None,
                 alpha=None,
                 phi_1=None,
                 phi_2=None,
                 qs=None,
                 **kwargs):
        

        
        if I_s is not None: self.I_s = I_s
        if I_k is not None: self.I_k = I_k
        if I_1 is not None: self.I_1 = I_1
        if i_1 is not None: self.i_1 = i_1
        if I_2 is not None: self.I_2 = I_2
        if I_b is not None: self.I_b = I_b
        if D is not None: self.D = D
        if G is not None: self.G = G
        if g is not None: self.g = g
        if mu is not None: self.mu = mu
        if M_T is not None: self.M_T = M_T
        if M_s is not None: self.M_s = M_s
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if alpha is not None: self.alpha = alpha
        if A is not None: self.A=A
        if B is not None: self.B=B
        if dphi is not None: self.dphi = dphi
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._engine = Disk(self.I_s, self.phi_1, qs=self.qs)
#         self.disk_k = Disk(I_k, phi_1, qs=qs)
        self._disk_1 = Disk(self.I_1, self.phi_1, qs=self.qs)
        self._disk_2 = Disk(self.I_2, self.phi_2, qs=self.qs)
#         self._disk_3 = Disk(self.I_3, self.phi_2, qs=self.qs)
#         self._disk_4 = Disk(self.I_4, self.phi_3, qs=self.qs)
        self._disk_B = Disk(self.I_b, self.phi_2, qs=self.qs)
        self._load = MaterialPoint(self.G/self.g, pos1=self.phi_2* self.D/2, qs=self.qs)
#         self._friction_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * cos(self.alpha)*self.mu , qs=self.qs)
        self._force = Force(-self.M_T,pos1=self.phi_2,qs=self.qs)
        self._drive = Force(self.M_s,pos1=self.phi_1,qs=self.qs)
        self._gravity_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_2* self.D/2 * sin(self.alpha) , qs=self.qs)


        components['engine'] = self._engine
        components['disk_1'] = self._disk_1
        components['disk_2'] = self._disk_2
#         components['disk_3'] = self._disk_3
#         components['disk_4'] = self._disk_4
        components['disk_B'] = self._disk_B
        components['load'] = self._load
#         components['friction_comp'] = self._friction_comp
        components['force'] = self._force
        components['drive'] = self._drive
        components['gravity_comp'] = self._gravity_comp
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
            self.I_s: [20,30,40],
            self.I_k: [20,30,40],
            self.I_1: [20,30,40],
            self.I_2: [20,30,40],
            self.I_b: [20,30,40],
            self.i_1: [20,30,40],
            self.D: [10,20,30],
            self.G: [10,20,30],
            self.g: [10,20,30],
            self.mu : [10,20,30],
            self.M_T : [10,20,30],
            self.M_s : [10,20,30],
            self.phi_1: [10,20,30],
            self.phi_2: [10,20,30],
            self.alpha  : [10,20,30],
        }

        return default_data_dict
        
    def steady_angular_velocity(self):
        obj=self
        eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)
        return omega_steady_state
    
    def reduced_torque(self):
        obj=self
        
        red_tor = -(obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())

        return delta_1

    def reduced_inertia(self):
        obj = self
        
        ans = obj.inertia_matrix()[0]

        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.02)-1)*I_r).n(5))
        return ik
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_2
        return st_val