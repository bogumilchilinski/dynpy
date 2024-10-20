from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy import Heaviside
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, CombustionEngine
from .elements import Resistor, Inductor, Capacitor, VoltageSource

import base64
import random
import IPython as IP
import numpy as np
import inspect

from dynpy.models.mechanics.trolley import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin
from dynpy.models.electric.engine import ElectromotiveForce, RotorTorque

class ElectricMotorcycle2DoF(ComposedSystem):
    scheme_name = 'DC_motor.png'
    
    U_z = Symbol('U_z',positive=True)
    R_w = Symbol('R_w', positive=True)
    L_w = Symbol('L_w', positive=True)
    k_e = Symbol('k_e', positive=True)
    B = Symbol('B', positive=True)
    J = Symbol('J', positive=True)
    k_m = Symbol('k_m', positive=True)
    
    r = Symbol('r',positive=True)
    
    Cd=Symbol('C_d',positive=True) #dragcoefficient
    Ad=Symbol('A_d',positive=True) #air density
    Af=Symbol('A_f',positive=True) #front area of the vehicle
    i=Symbol('i',positive=True) #przełożenie
    
    Rpwm=Symbol('R_pwm',positive=True)
    Rpwm2=Symbol('R_pwm2',positive=True)
    
    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    ft=Symbol('f_t',positive=True)
    eta_p=Symbol('\eta_p',positive=True)
    v=Symbol('v',positive=True)
    
    charge=dynamicsymbols('q')
    x=dynamicsymbols('x')

    
    def __init__(self,
                 U_z=None,
                 R_w=None,
                 L_w=None,  
                 k_e=None,
                 B=None, 
                 J=None, 
                 k_m=None,
                 r=None,
                 Cd=None,
                 Ad=None,
                 Af=None,
                 Rpwm=None,
                 Rpwm2=None,
                 i=None,
                 m=None,
                 g=None,
                 ft=None,
                 eta_p=None,
                 v=None,
                 ivar=Symbol('t'),
                 **kwargs):
        
        if U_z is not None: self.U_z = U_z
        if R_w is not None: self.R_w = R_w
        if L_w is not None: self.L_w = L_w
        if k_e is not None: self.k_e = k_e
        if B is not None: self.B = B
        if J is not None: self.J = J
        if k_m is not None: self.k_m = k_m
        if r is not None: self.r = r
        
        if Cd is not None: self.Cd = Cd
        if Ad is not None: self.Ad = Ad
        if Af is not None: self.Af = Af
        if i is not None: self.i=i
        
        if Rpwm is not None: self.Rpwm=Rpwm
        if Rpwm2 is not None: self.Rpwm2=Rpwm2
        
        if m is not None: self.m=m
        if g is not None: self.g=g
        if ft is not None: self.ft=ft
        if eta_p is not None: self.eta_p=eta_p
        if v is not None: self.v=v
        
        self.ivar = ivar
        
        self.qs = [self.charge, self.x]
        self._init_from_components(**kwargs)
        
    @property
    def components(self):
        components = {}
        
        self.resistance=Resistor(-self.R_w,self.charge,qs=self.qs,ivar=self.ivar)
        self.inductor = Inductor(-self.L_w, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(-self.U_z, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.inertia=Inductor(-self.J, self.x/self.r, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.viscous_friction=Resistor(-self.B,self.x/self.r,qs=self.qs,ivar=self.ivar)
        self.electromotive_force = ElectromotiveForce(-self.k_m*self.charge.diff(self.ivar)*self.i, self.x/self.r, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(self.k_e*self.x.diff(self.ivar)/self.r, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.air_drag=Force(S.One/2*self.Ad*self.Cd*self.Af*(self.x.diff(self.ivar))**2,self.x,qs=self.qs,ivar=self.ivar)
        self.pwm=Force(self.Rpwm*Heaviside(-(self.ivar)+5)*Heaviside(self.charge.diff(self.ivar)-700)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.pwm2=Force(self.Rpwm2*Heaviside((self.ivar)-5)*Heaviside(self.charge.diff(self.ivar)-300)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)

        components['resistor'] = self.resistance
        components['inductor'] = self.inductor
        components['voltage_source'] = self.voltage_source
        components['inertia'] = self.inertia
        components['viscous_friction'] = self.viscous_friction
        components['electromotive_force'] = self.electromotive_force
        components['rotor_torque'] = self.rotor_torque
        components['air_drag'] = self.air_drag
        components['pwm']=self.pwm
        components['pwm2']=self.pwm2

        return components
    
    def required_power(self):
        preq=Symbol('P_req',positive=True)
        preq_eq=Eq(preq,(self.m*self.g*self.ft+(self.Cd*self.Ad*(self.v**2)/2*self.Af))/self.eta_p*self.v)
        return preq_eq
    
    def power_dict(self):
        return {self.m:350,self.g:9.81,self.ft:0.015,self.Cd:0.65,self.Ad:1.2047,self.Af:0.6,self.eta_p:0.98,self.v:41.6}
        
    def symbols_description(self):
        i_roc=self.i_w.diff(self.ivar)
        self.sym_desc_dict = {
            self.U_z: r'voltage supplying the rotor',
            self.R_w: r'equivalent resistance of the rotor windings',
            self.L_w: r'equivalent inductance of the rotor windings',
            self.E: r'electromotive force of induction',
            self.U_Rw: r'voltage across the rotor winding resistance',
            self.U_Lw: r'voltage related to the rotor inductance',
            self.k_e: r'back electromotive force constant',
            self.M_s: r'rotor torque',
            self.B: r'coefficient of viscous friction reduced to the rotor shaft',
            self.J: r'moment of inertia reduced to the rotor shaft',
            self.M_obc: r'engine load torque',
            self.k_m: r'torque constant',
            self.i_w: r'rotor winding current',
            i_roc:r'rotor winding current rate of change',
            self.omega_s:r'angular velocity of the rotor',
            self.M_a:r'rotor angular acceleration torque',
            self.M_r:r'rotor motion resistance torque',
            self.omega_s.diff(self.ivar):r'angular acceleration of the rotor',
        }
        return self.sym_desc_dict
    
    def get_default_data(self):
        default_data_dict = {
            self.R_w: [2],
            self.L_w: [0.1],
            self.k_e: [0.1],
            self.k_m: [0.1],
            self.J: [0.1],
            self.B: [0.5],
            self.M_obc: [0.2],
            self.r:[0.5],
            self.c:[1000],
            self.m:[250],
            self.Cd:[0.65],
            self.Ad:[1.2047],
            self.Af:[0.6]
        }
        return default_data_dict
    

    def _dimensionless_ode(self):
        
        from ...solvers.linear import ODESystem,FirstOrderLinearODESystem
        
        t = self.ivar
        tau = Symbol('tau')
        
        syms_dict = {
                    # dcmotor.U_z: 0,
                    self.R_w: 1,
                    self.L_w: 1,
                    # dcmotor.k_e: Symbol('kappa',positive=True),
                    # dcmotor.k_m: Symbol('kappa',positive=True),
                    
                    self.k_e: Symbol('kappa',negative=True),
                    self.k_m: Symbol('kappa',negative=True),
                    
                    self.J: 1,
                    self.B: 1,
                    t:tau,
                    # dcmotor.M_obc: 10,
                    }
        dvars = Matrix([[self.i_w], [self.omega_s]]).subs({t:tau})
        
        
        ode = ODESystem(odes=(self._eoms.subs(syms_dict))*-1, dvars = dvars, ode_order=1,ivar=tau)
        
        #return FirstOrderLinearODESystem.from_ode_system(ode)
        return ode