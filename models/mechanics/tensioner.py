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


class BlowerToothedBelt(ComposedSystem):

    scheme_name = 'blower_toothed_belt.png'
    real_name = 'blown_440_big_block.jpg'
    detail_scheme_name = 'blower_roller_bolt.png'
    detail_real_name = 'tensioner_pulley.jpg'

    m=Symbol('m', positive=True)
    k_belt=Symbol('k_b', positive=True)
    k_tensioner=Symbol('k_t', positive=True)
    ivar=Symbol('t')
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    z0 = Symbol('z_0', positive=True)
    z=dynamicsymbols('z')


    def __init__(self,
                 m=None,
                 k_belt=None,
                 k_tensioner=None,
                 ivar=Symbol('t'),
                 Omega=None,
                 F=None,
                 z0=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_belt is not None: self.k_belt = k_belt
        if k_tensioner is not None: self.k_tensioner = k_tensioner
        if F is not None: self.F = F
        if Omega is not None: self.Omega = Omega
        if z is not None: self.z = z
        if ivar is not None: self.ivar = ivar
        if z0 is not None: self.z0 = z0

        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z, qs=[self.z])(label = 'Material point')
        self._upper_belt = Spring(self.k_belt, self.z, pos2=0, qs=[self.z])(label = 'Upper belt stiffness')
        self._lower_belt = Spring(self.k_belt, self.z, pos2=0, qs=[self.z])(label = 'Lower belt stiffness')
        self._tensioner = Spring(self.k_tensioner, self.z, pos2=self.z0, qs=[self.z])(label = 'Tensioner stiffness')
        self._force = Force(self.F * sin(self.Omega * self.ivar), pos1=self.z, qs=[self.z])(label = 'Force')



        components['_mass'] = self._mass
        components['_upper_belt'] = self._upper_belt
        components['_lower_belt'] = self._lower_belt
        components['_tensioner'] = self._tensioner
        components['_force'] = self._force
        return components

    def get_default_data(self):


        k0, F0 = symbols('k0 F0', positive=True)

        default_data_dict = {
            self.F: [F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0],
            self.z0: [F0 / k0 * S.One / 4 * no for no in range(1, 2)],
        }

        return default_data_dict
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m: [2 * no for no in range(1, 8)],
            self.k_belt: [10e3 * no for no in range(1, 8)],
            self.k_tensioner: [10e5 * no for no in range(4, 8)],
            self.F: [500 * no for no in range(1, 8)],
            self.Omega: [0.02*no for no in range(1,4)],
            self.Q: [200* no for no in range(1,8)],
            self.z0: [0.2 * no for no in range(1, 2)],
        }

        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of the system on the spring',
            self.k_belt: r'Belt stiffness factor',
            self.k_tensioner: r'Tensioner spring stiffnessfactor',
            self.Omega: r'Angle of exciting function',
            self.F: r'Exciting force',
            self.Q: r'Load of the system',
            self.z: r'Displacement of the system'
        }

        return self.sym_desc_dict

    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST,
        mech_comp.TensionerForce,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticPinDiameter,
        mech_comp.DynamicPinDiameter
        ]

        return comp_list
    
    def tensioner_belt_force(self):
        return self.k_tensioner * self.steady_solution()[0]

    def left_belt_force(self):
        return self.k_belt * self.steady_solution()[0]

    def right_belt_force(self):
        return self.k_belt * self.steady_solution()[0]

    def max_static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
#        abs(self.static_load().doit()[0])
        return abs(free_coeff)

    
    def dynamic_force(self):

        amps = self._fodes_system.steady_solution.as_dict()
        force_in_tensioner = self.components['_tensioner'].force().doit().expand()#.doit()
        data=self._given_data

        ##### new force implementation
        ans = force_in_tensioner.subs(amps).expand().doit()

        return ans.subs(data)
    
    
    def max_dynamic_force(self):

        ans=self.dynamic_force()
        data=self._given_data
        sin_coeff=ans.coeff(sin(self.Omega*self.ivar)).subs(data)
        #display(sincoeff)
        cos_coeff=ans.coeff(cos(self.Omega*self.ivar)).subs(data)
        #display(coscoeff)
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        #display(freecoeff)
        return sqrt(sin_coeff**2 + cos_coeff**2) + abs(free_coeff)

    def static_force_pin_diameter(self):
        kt = Symbol('k_shear', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_shear', positive=True)
        Re = Symbol('R_e', positive=True)
        data=self._given_data
        load = abs(self.max_dynamic_force())
        return (((4 * load) / (pi * kt * Re))**(1 / 2)).subs(data)




class DampedBlowerToothedBelt(BlowerToothedBelt):

    scheme_name = 'damped_blower_toothed_belt.png'
    real_name = 'blown_440_big_block.jpg'
    detail_scheme_name = 'blower_roller_bolt.png'
    detail_real_name = 'tensioner_pulley.jpg'


    c_belt = Symbol('c_b', positive=True)
    c_tensioner = Symbol('c_t', positive=True)
    lam = Symbol('lambda', positive=True)
    z0 = Symbol('z0', positive=True)

    def __init__(
            self,
            m=None,
            k_belt=None,
            k_tensioner=None,
            ivar=Symbol('t'),
            Omega=None,
            F=None,
            z0=None,
            z=None,
            c_belt=None,
            c_tensioner=None,
            lam=None,
            qs=None,
            **kwargs):
        if lam is not None: self.lam = lam
        if c_belt is not None: self.c_belt = c_belt
        if c_tensioner is not None: self.c_tensioner = c_tensioner
        if z0 is not None: self.z0 = z0
        if qs is not None: self.qs = qs
        else: self.qs = [self.z]

        super().__init__(m=m, k_belt=k_belt, k_tensioner=k_tensioner, ivar=ivar, Omega=Omega, F=F, z0=z0, z=z)

    @property
    def components(self):

        components = {}

        self._blower_toothed_belt = BlowerToothedBelt(m=self.m, k_belt=self.k_belt, k_tensioner=self.k_tensioner, ivar=self.ivar, Omega=self.Omega, F=self.F, z0=self.z0, z=self.z)(label = 'Blower Toothed Belt')
        self._upper_belt_damping = Damper(self.c_belt, pos1=self.z, pos2=0, qs=[self.z])(label = 'Upper belt damping')
        self._lower_belt_damping = Damper(self.c_belt, pos1=self.z, pos2=0, qs=[self.z])(label = 'Lower belt damping')
        self._tensioner_damping = Damper(self.c_tensioner, pos1=self.z, pos2=self.z0, qs=[self.z])(label = 'Tensioner damping')

        components['_blower_toothed_belt'] = self._blower_toothed_belt
        components['_upper_belt_damping'] = self._upper_belt_damping
        components['_lower_belt_damping'] = self._lower_belt_damping
        components['_tensioner_damping'] = self._tensioner_damping

        return components

    def get_default_data(self):

        default_data_dict = {

            self.c_belt: [self.lam * (self.k_belt)],
            self.c_tensioner: [self.lam * (self.k_tensioner)],
        }

        return default_data_dict
    
    def reference_data_dict(self):

        lam, m0, k0, l0, F0 = symbols('lambda m_0 k_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.c_belt: [self.lam*(self.k_belt)],
            self.c_tensioner: [self.lam*(self.k_tensioner)],
        }

        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.c_belt: [(self.k_belt)*self.lam],
            self.c_tensioner: [(self.k_tensioner)*self.lam],

            self.m: [2 * no for no in range(1, 8)],
            self.k_belt: [1000 * no for no in range(1, 8)],
            self.k_tensioner: [10000 * no for no in range(4, 8)],
            self.F: [500 * no for no in range(1, 8)],
            self.Omega: [2*3.14*no/4 for no in range(1,4)],
            self.lam: [no/50 for no in range(7,13)],
            self.z0: [0.2 * no for no in range(1, 2)],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of the system on the spring',
            self.k_belt: r'Belt stiffness factor',
            self.k_tensioner: r'Tensioner spring stiffness factor',
            self.c_belt: r'Belt damping factor',
            self.c_tensioner: r'Tensioner damping factor',
            self.Omega: r'Angle of exciting function',
            self.F: r'Exciting force',
            self.z: r'Displacement of the system'
        }

        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list=[
        *REPORT_COMPONENTS_LIST,
        mech_comp.TensionerForce,
        mech_comp.TensionerDamperForce,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticPinDiameter,
        mech_comp.DynamicPinDiameter,
        mech_comp.DampedVibrationFrequency
        ]

        return comp_list


    def tensioner_damper_force(self):
        return self.c_tensioner * self.steady_solution()[0].diff(self.ivar)

    def left_belt_damper_force(self):
        return self.c_belt * self.steady_solution()[0].diff(self.ivar)

    def right_belt_damper_force(self):
        return self.c_belt * self.steady_solution()[0].diff(self.ivar)
