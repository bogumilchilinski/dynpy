from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

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

class DCMotor(ComposedSystem):
    scheme_name = 'DC_motor.png'
    
    U_z = Symbol('U_z', positive=True)
    R_w = Symbol('R_w', positive=True)
    L_w = Symbol('L_w', positive=True)
    E = Symbol('E', positive=True)
    U_Rw = Symbol('U_Rw', positive=True)
    U_Lw = Symbol('U_Lw', positive=True)
    k_e = Symbol('k_e', positive=True)
    M_s = Symbol('M_s', positive=True)
    B = Symbol('B', positive=True)
    J = Symbol('J', positive=True)
    M_obc = Symbol('M_l', positive=True)
    k_m = Symbol('k_m', positive=True)
    M_a = Symbol('M_a', positive=True)
    M_r = Symbol('M_r', positive=True)
    i_w=dynamicsymbols('i_w')
    omega_s=dynamicsymbols('omega_s')
    

    
    def __init__(self,
                 U_z=None, 
                 R_w=None,
                 L_w=None, 
                 E=None, 
                 U_Rw=None, 
                 U_Lw=None, 
                 k_e=None,
                 M_s=None, 
                 B=None, 
                 J=None, 
                 M_obc=None, 
                 k_m=None,
                 M_a=None,
                 M_r=None,
                 ivar=Symbol('t'),
                 **kwargs):
        
        if U_z is not None: self.U_z = U_z
        if R_w is not None: self.R_w = R_w
        if L_w is not None: self.L_w = L_w
        if E is not None: self.E = E
        if U_Rw is not None: self.U_Rw = U_Rw
        if U_Lw is not None: self.U_Lw = U_Lw
        if k_e is not None: self.k_e = k_e
        if M_s is not None: self.M_s = M_s
        if B is not None: self.B = B
        if J is not None: self.J = J
        if M_obc is not None: self.M_obc = M_obc
        if k_m is not None: self.k_m = k_m
        if M_a is not None: self.M_a = M_a
        if M_r is not None: self.M_r = M_r
        
        self.ivar = ivar
        self.qs = [self.i_w, self.omega_s]
        self._init_from_components(**kwargs)
        
    @property
    def components(self):
        components = {}
        
        self.resistor = CurrentDependentResistor(-self.R_w, self.i_w,  qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.inductor = Resistor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
#         self.inductor = Inductor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(-self.U_z, self.i_w, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.electromotive_force = ElectromotiveForce(self.k_e*self.omega_s, self.i_w, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.engine_load_torque = EngineLoadTorque(self.M_obc, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(-self.k_m*self.i_w, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_angular_acceleration_torque = RotorAngularAccelerationTorque(-self.J, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_resistance_movement_torque = ResistanceMovementTorque(-self.B, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        
        components['resistor'] = self.resistor
        components['inductor'] = self.inductor
        components['voltage_source'] = self.voltage_source
        components['electromotive_force'] = self.electromotive_force
        components['engine_load_torque'] = self.engine_load_torque
        components['rotor_torque'] = self.rotor_torque
        components['rotor_angular_acceleration_torque'] = self.rotor_angular_acceleration_torque
        components['rotor_resistance_movement_torque'] = self.rotor_resistance_movement_torque
        
        return components
        
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
            self.U_z: [10],
            self.R_w: [2],
            self.L_w: [0.1],
            self.k_e: [0.1],
            self.k_m: [0.1],
            self.J: [0.1],
            self.B: [0.5],
            self.M_obc: [0.2]
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
    

class CurrentDependentResistor(Spring):
    """
    Model of Resistor
    
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    def __init__(self, resistance, q0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(stiffness=resistance, pos1=q0, pos2=0, qs=qs, ivar=Symbol('t'), frame=base_frame)
        
        
class ElectromotiveForce(Force):
    """
    Model of a Electromotive Force
    """
    def __init__(self, electric_constant, q0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(force=electric_constant, pos1=q0,  qs=qs, ivar=ivar, frame=frame)
        
        
class EngineLoadTorque(Force):
    """
    Model of a engine load torque
    """
    def __init__(self, torque, q0,  qs=None, ivar=Symbol('t'), frame = base_frame):
        super().__init__(force=torque, pos1=q0, qs=qs, ivar=ivar, frame=frame)
        
        
class RotorTorque(Force):
    """
    Model of a Electromotive Force
    """
    def __init__(self, mechanical_constant, q0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(force=mechanical_constant, pos1=q0, qs=qs, ivar=ivar, frame=frame)
        
        
class RotorAngularAccelerationTorque(Damper):
    """
    Model of Resistor
    
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    def __init__(self, moment_of_interia, q0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(c=moment_of_interia, pos1=q0, pos2=0, qs=qs, ivar=Symbol('t'), frame=base_frame)
        
class ResistanceMovementTorque(Spring):
    """
    Model of a Electromotive Force
    """
    def __init__(self, friction_coefficient, q0, q1=0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(stiffness=friction_coefficient, pos1=q0, pos2=q1,l_0=0 ,  qs=qs, ivar=ivar, frame=frame)