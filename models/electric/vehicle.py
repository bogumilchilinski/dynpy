from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy import Heaviside, DiracDelta
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
from dynpy.models.electric.motor import ElectromotiveForce, RotorTorque
from dynpy.solvers.linear import ODESystem

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
    i=Symbol('i_r',positive=True) #przełożenie
    
    Rpwm=Symbol('R_{pwm1}',positive=True)
    Rpwm2=Symbol('R_{pwm2}',positive=True)
    
    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)
    ft=Symbol('f_t',positive=True)
    eta_p=Symbol('\eta_p',positive=True)
    v=Symbol('v',positive=True)
    
    charge=dynamicsymbols('q')
    x=dynamicsymbols('x')

    current=dynamicsymbols('i')
    velocity=dynamicsymbols('v')
    
    ic1=Symbol('i_{c1}',positive=True)
    ic2=Symbol('i_{c2}',positive=True)
    
    time_cut=Symbol('t_c',positive=True)
    
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
                 ic1=None,
                 ic2=None,
                 time_cut=None,
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
        
        if ic1 is not None: self.ic1=ic1
        if ic2 is not None: self.ic2=ic2
        if time_cut is not None: self.time_cut=time_cut
        
        self.preq=Symbol('P_{req}',positive=True)
        
        self.ivar = ivar
        
        self.qs = [self.charge, self.x]
        self._init_from_components(**kwargs)
        
    @property
    def components(self):
        self.cut1=self.ic1
        self.cut2=self.ic2
        components = {}
        
        self.resistance=Resistor(-self.R_w,self.charge,qs=self.qs,ivar=self.ivar)
        self.inductor = Inductor(-self.L_w, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(-self.U_z, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.inertia=Inductor(-self.J, self.x/self.r, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.viscous_friction=Resistor(-self.B,self.x/self.r,qs=self.qs,ivar=self.ivar)
        self.electromotive_force = ElectromotiveForce(-self.k_m*self.charge.diff(self.ivar)*self.i, self.x/self.r, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(self.k_e*self.x.diff(self.ivar)/self.r, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.air_drag=Force(S.One/2*self.Ad*self.Cd*self.Af*(self.x.diff(self.ivar))**2,self.x,qs=self.qs,ivar=self.ivar)
        self.pwm=Force(self.Rpwm*Heaviside(-(self.ivar)+self.time_cut)*Heaviside(self.charge.diff(self.ivar)-self.cut1)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.pwm2=Force(self.Rpwm2*Heaviside((self.ivar)-self.time_cut)*Heaviside(self.charge.diff(self.ivar)-self.cut2)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)

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

        preq_eq=Eq(self.preq,(self.m*self.g*self.ft+(self.Cd*self.Ad*(self.v**2)/2*self.Af))/self.eta_p*self.v)
        return preq_eq
    
    def power_dict(self):
        return {self.m:350,self.g:9.81,self.ft:0.015,self.Cd:0.65,self.Ad:1.2047,self.Af:0.6,self.eta_p:0.98,self.v:41.6}
        
    def symbols_description(self):
        i_roc=self.charge.diff(self.ivar)
        self.sym_desc_dict = {
            self.U_z: r'voltage supplying the rotor',
            self.R_w: r'equivalent resistance of the rotor windings',
            self.L_w: r'equivalent inductance of the rotor windings',
            self.eta_p:r'efficiency',
            self.Cd:r'aerodynamic drag coefficient',
            self.Ad:r'air density',
            self.k_e: r'back electromotive force constant',
            self.v:r'maximum velocity',
            self.B: r'coefficient of viscous friction reduced to the rotor shaft',
            self.J: r'moment of inertia reduced to the rotor shaft',
            self.Af:r'frontal area',
            self.k_m: r'torque constant',
            self.charge:r'electric charge',
            self.charge.diff(self.ivar): r'current',
            self.charge.diff(self.ivar,self.ivar):r'current rate of change',
#             i_roc:r'rotor winding current rate of change',
            self.m:r'mass',
            self.g:r'gravitational acceleration',
            self.ft:r'friction coefficient',
#             self.omega_s.diff(self.ivar):r'angular acceleration of the rotor',
            self.preq:r'required power',
            self.time_cut:r'current swich time',
            self.Rpwm2:r'lower current control resistance',
            self.Rpwm:r'higher current control resistance',
            self.ic1:r'higher current cut value',
            self.ic2:r'lower current cut value',
            self.i:r'final drive ratio',
            self.r:r'wheel radius',
            self.x:r'linear displacement',
            self.x.diff(self.ivar):r'velocity',
            self.x.diff(self.ivar,self.ivar):r'acceleration',
            self.velocity.diff(self.ivar):r'acceleration',
            self.current:r'current',
            self.current.diff(self.ivar):r'current rate of change',
            self.velocity:r'velocity'
        }
        return self.sym_desc_dict
    
    def get_default_data(self):
        default_data_dict = {
            self.R_w: 0.0105,
            self.L_w: 0.000182,
            self.k_e: 0.224,
            self.k_m: 0.224,
            self.J: 20,
            self.B: 0.004794,
            self.r:0.3175,
            self.m:250,
            self.Cd:0.65,
            self.Ad:1.2047,
            self.Af:0.6,
#             self.eps:0.18,
#             self.delta:0,
            self.U_z:96,
            self.i:1.925,
            self.Rpwm:0.126,
            self.Rpwm2:0.28,
            self.time_cut:5,
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
    
    @property
    def as_first_order_em(self):
        inert_mat=self.mass_matrix
        damp_mat=self.damping_matrix().subs(self.charge,self.current).subs(DiracDelta(self.current.diff(self.ivar)-self.cut1),0).subs(Heaviside(self.current.diff(self.ivar)-self.cut1),Heaviside(self.current-self.cut1)).subs(self.Ad*self.Af*self.Cd*self.x.diff(self.ivar),self.Ad*self.Af*self.Cd*self.velocity*S.Half).subs(DiracDelta(self.current.diff(self.ivar)-self.cut2),0).subs(Heaviside(self.current.diff(self.ivar)-self.cut2),Heaviside(self.current-self.cut2))
        sec_to_first_comp=inert_mat*Matrix([[self.current.diff(self.ivar)],[self.velocity.diff(self.ivar)]])
        first_to_zero_comp=(damp_mat*Matrix([[self.current],[self.velocity]]))
        eoms_mat=sec_to_first_comp+first_to_zero_comp-self.external_forces()
        return ODESystem(-eoms_mat,Matrix([[self.current],[self.velocity]]),ivar=self.ivar,ode_order=None)
    
class ElectricMotorcycle2DoFMSM(ElectricMotorcycle2DoF):
    eps=Symbol('varepsilon',positive=True)
    delta=Symbol('delta',positive=True)
    
    @property
    def components(self):
        components = {}
        self.cut=self.ic
        r_pwm_1=self.Rpwm*(Heaviside(self.charge.diff(self.ivar)-self.cut)*self.eps+1)
        
        
        self.resistance=Resistor(-self.R_w,self.charge,qs=self.qs,ivar=self.ivar)
        self.inductor = Inductor(-self.L_w, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(-self.U_z, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.inertia=Inductor(-self.J, self.x/self.r, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.viscous_friction=Resistor(-self.B,self.x/self.r,qs=self.qs,ivar=self.ivar)
        self.electromotive_force = ElectromotiveForce(-self.k_m*self.charge.diff(self.ivar)*self.i, self.x/self.r, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(self.k_e*self.x.diff(self.ivar)/self.r, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.air_drag=Force(S.One/2*self.Ad*self.Cd*self.Af*(self.x.diff(self.ivar))**2,self.x,qs=self.qs,ivar=self.ivar)
        self.pwm=Force(r_pwm_1*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
#         self.pwm2=Force(self.Rpwm2*Heaviside((self.ivar)-5)*Heaviside(self.charge.diff(self.ivar)-300)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)

        components['resistor'] = self.resistance
        components['inductor'] = self.inductor
        components['voltage_source'] = self.voltage_source
        components['inertia'] = self.inertia
        components['viscous_friction'] = self.viscous_friction
        components['electromotive_force'] = self.electromotive_force
        components['rotor_torque'] = self.rotor_torque
        components['air_drag'] = self.air_drag
        components['pwm']=self.pwm
#         components['pwm2']=self.pwm2

        return components

    @property
    def as_first_order_em(self):
        inert_mat=self.mass_matrix
        damp_mat=self.damping_matrix().subs(self.charge,self.current).subs(DiracDelta(self.current.diff(self.ivar)-self.cut),0).subs(Heaviside(self.current.diff(self.ivar)-self.cut),Heaviside(self.current-self.cut)).subs(self.Ad*self.Af*self.Cd*self.x.diff(self.ivar),self.Ad*self.Af*self.Cd*self.velocity*S.Half)
        sec_to_first_comp=inert_mat*Matrix([[self.current.diff(self.ivar)],[self.velocity.diff(self.ivar)]])
        first_to_zero_comp=(damp_mat*Matrix([[self.current],[self.velocity]]))
        eoms_mat=sec_to_first_comp+first_to_zero_comp-self.external_forces()
        return ODESystem(eoms_mat,Matrix([[self.current],[self.velocity]]),ivar=self.ivar,ode_order=None)
    def get_default_data(self):
        default_data_dict = {
            self.R_w: 0.0105,
            self.L_w: 0.000182,
            self.k_e: 0.224,
            self.k_m: 0.224,
            self.J: 20,
            self.B: 0.004794,
            self.r:0.3175,
            self.m:250,
            self.Cd:0.65,
            self.Ad:1.2047,
            self.Af:0.6,
            self.eps:0.18,
            self.delta:0,
            self.U_z:96,
            self.i:1.925,
            self.Rpwm:0.18
        }
        return default_data_dict
    
class RegenerativeBrakingElectricMotorcycle2DoFMSM(ElectricMotorcycle2DoFMSM):
    eps=Symbol('varepsilon',positive=True)
    delta=Symbol('delta',positive=True)
    
    @property
    def components(self):
        components = {}
        r_pwm_1=self.Rpwm*(1+self.eps*Heaviside(self.charge.diff(self.ivar)+20))
        
        
        self.resistance=Resistor(-self.R_w,self.charge,qs=self.qs,ivar=self.ivar)
        self.inductor = Inductor(-self.L_w, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(-self.U_z, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.inertia=Inductor(-self.J, self.x/self.r, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.viscous_friction=Resistor(-self.B,self.x/self.r,qs=self.qs,ivar=self.ivar)
        self.electromotive_force = ElectromotiveForce(-self.k_m*self.charge.diff(self.ivar)*self.i, self.x/self.r, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(self.k_e*self.x.diff(self.ivar)/self.r, self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.air_drag=Force(S.One/2*self.Ad*self.Cd*self.Af*(self.x.diff(self.ivar))**2,self.x,qs=self.qs,ivar=self.ivar)
        self.pwm=Force(r_pwm_1*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)
#         self.pwm2=Force(self.Rpwm2*Heaviside((self.ivar)-5)*Heaviside(self.charge.diff(self.ivar)-300)*self.charge.diff(self.ivar), self.charge, qs=self.qs, ivar=self.ivar , frame=base_frame)

        components['resistor'] = self.resistance
        components['inductor'] = self.inductor
        components['voltage_source'] = self.voltage_source
        components['inertia'] = self.inertia
        components['viscous_friction'] = self.viscous_friction
        components['electromotive_force'] = self.electromotive_force
        components['rotor_torque'] = self.rotor_torque
        components['air_drag'] = self.air_drag
        components['pwm']=self.pwm
#         components['pwm2']=self.pwm2

        return components
    
    
    
    
class adtest(ComposedSystem):
    scheme_name = 'DC_motor.png'
    

    
    Cd=Symbol('C_d',positive=True) #dragcoefficient
    Ad=Symbol('A_d',positive=True) #air density
    Af=Symbol('A_f',positive=True) #front area of the vehicle

    m=Symbol('m',positive=True)
    g=Symbol('g',positive=True)

    x=dynamicsymbols('x')

    
    def __init__(self,
                 Cd=None,
                 Ad=None,
                 Af=None,
                 m=None,
                 ivar=Symbol('t'),
                 **kwargs):
        

        
        if Cd is not None: self.Cd = Cd
        if Ad is not None: self.Ad = Ad
        if Af is not None: self.Af = Af
        
        
        if m is not None: self.m=m


        
        self.ivar = ivar
        
        self.qs = [self.x]
        self._init_from_components(**kwargs)
        
    @property
    def components(self):

        components = {}
        
        self.mass=MaterialPoint(-self.m,self.x,qs=self.qs,ivar=self.ivar)
        self.air_drag=Force(S.One/2*self.Ad*self.Cd*self.Af*(self.x.diff(self.ivar))**2,self.x,qs=self.qs,ivar=self.ivar)
        self.friction=Force(0.015*self.m*9.81,self.x,qs=self.qs,ivar=self.ivar)

        components['resistor'] = self.mass

        components['air_drag'] = self.air_drag
        components['fric'] = self.friction

        return components
    
    def get_default_data(self):
        default_data_dict = {
            self.m:350,
            self.Cd:0.65,
            self.Ad:1.2047,
            self.Af:0.6,
        }
        return default_data_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    