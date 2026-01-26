import base64
import inspect
import random
from dataclasses import dataclass, field

import IPython as IP
import numpy as np
from sympy import (
    Abs,
    Eq,
    Function,
    Heaviside,
    Matrix,
    N,
    Number,
    S,
    Subs,
    Symbol,
    Tuple,
    asin,
    cos,
    diag,
    diff,
    dsolve,
    factorial,
    flatten,
    fraction,
    hessian,
    im,
    latex,
    oo,
    pi,
    sin,
    solve,
    solveset,
    sqrt,
    symbols,
    integrate,
)
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint
from .battery import BatteryCell

from dynpy.models.mechanics.trolley import (
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)

from ...dynamics import HarmonicOscillator, LagrangesDynamicSystem, mech_comp
from ..elements import (
    PID,
    CombustionEngine,
    Damper,
    Disk,
    Excitation,
    Force,
    GravitationalForce,
    MaterialPoint,
    RigidBody2D,
    Spring,
    base_frame,
    base_origin,
)
from .elements import Capacitor, Inductor, Resistor, VoltageSource


@dataclass
class DCMotor(ComposedSystem):
    scheme_name = "DC_motor.png"
    R_a=Symbol('R_a',real=True)
    L_a=Symbol('L_a',real=True)
    L_w=Symbol('L_w',real=True)
    R_w=Symbol('R_w',real=True)
    J=Symbol('J',real=True)
    B=Symbol('B',real=True)
    K_e=Symbol('K_e',real=True)
    K_t=Symbol('K_t',real=True)
    U_z=Symbol('U_z',real=True)
    k_e=Symbol('k_e',real=True)
    M_obc=Symbol('M_obc',real=True)
    k_m=Symbol('k_m',real=True)
    M_a=Symbol('M_a',real=True)
    M_r=Symbol('M_r',real=True)
    i_w: Symbol = dynamicsymbols("i_w")
    omega_s: Symbol = dynamicsymbols("omega_s")
    ivar: Symbol = Symbol("t")  # wartość domyślna

    def __init__(
        self,
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
        ivar=Symbol("t"),
        **kwargs
    ):

        if U_z is not None:
            self.U_z = U_z
        if R_w is not None:
            self.R_w = R_w
        if L_w is not None:
            self.L_w = L_w
        if E is not None:
            self.E = E
        if U_Rw is not None:
            self.U_Rw = U_Rw
        if U_Lw is not None:
            self.U_Lw = U_Lw
        if k_e is not None:
            self.k_e = k_e
        if M_s is not None:
            self.M_s = M_s
        if B is not None:
            self.B = B
        if J is not None:
            self.J = J
        if M_obc is not None:
            self.M_obc = M_obc
        if k_m is not None:
            self.k_m = k_m
        if M_a is not None:
            self.M_a = M_a
        if M_r is not None:
            self.M_r = M_r

        self.ivar = ivar
        self.qs = [self.i_w, self.omega_s]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.resistor = CurrentDependentResistor(
            self.R_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.inductor = Resistor(
            self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        #         self.inductor = Inductor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(
            self.U_z, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.electromotive_force = ElectromotiveForce(
            -self.k_e * self.omega_s,
            self.i_w,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        self.engine_load_torque = EngineLoadTorque(
            -self.M_obc, self.omega_s, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.rotor_torque = RotorTorque(
            self.k_m * self.i_w,
            self.omega_s,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        self.rotor_angular_acceleration_torque = RotorAngularAccelerationTorque(
            self.J, self.omega_s, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.rotor_resistance_movement_torque = ResistanceMovementTorque(
            self.B, self.omega_s, qs=self.qs, ivar=self.ivar, frame=base_frame
        )

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor
        components["voltage_source"] = self.voltage_source
        components["electromotive_force"] = self.electromotive_force
        components["engine_load_torque"] = self.engine_load_torque
        components["rotor_torque"] = self.rotor_torque
        components["rotor_angular_acceleration_torque"] = (
            self.rotor_angular_acceleration_torque
        )
        components["rotor_resistance_movement_torque"] = (
            self.rotor_resistance_movement_torque
        )

        return components

    def symbols_description(self):
        i_roc = self.i_w.diff(self.ivar)
        self.sym_desc_dict = {
            self.U_z: r"voltage supplying the rotor",
            self.R_w: r"equivalent resistance of the rotor windings",
            self.L_w: r"equivalent inductance of the rotor windings",
            self.E: r"electromotive force of induction",
            self.U_Rw: r"voltage across the rotor winding resistance",
            self.U_Lw: r"voltage related to the rotor inductance",
            self.k_e: r"back electromotive force constant",
            self.M_s: r"rotor torque",
            self.B: r"coefficient of viscous friction reduced to the rotor shaft",
            self.J: r"moment of inertia reduced to the rotor shaft",
            self.M_obc: r"engine load torque",
            self.k_m: r"torque constant",
            self.i_w: r"rotor winding current",
            i_roc: r"rotor winding current rate of change",
            self.omega_s: r"angular velocity of the rotor",
            self.M_a: r"rotor angular acceleration torque",
            self.M_r: r"rotor motion resistance torque",
            self.omega_s.diff(self.ivar): r"angular acceleration of the rotor",
        }
        return self.sym_desc_dict

    def get_default_data(self):
        default_data_dict = {
            self.U_z: 10,
            self.R_w: 2,
            self.L_w: 0.1,
            self.k_e: 0.1,
            self.k_m: 0.1,
            self.J: 0.1,
            self.B: 0.5,
            self.M_obc: 0.2,
        }
        return default_data_dict

    def _dimensionless_ode(self):

        from ...solvers.linear import FirstOrderLinearODESystem, ODESystem

        t = self.ivar
        tau = Symbol("tau")

        syms_dict = {
            # dcmotor.U_z: 0,
            self.R_w: 1,
            self.L_w: 1,
            # dcmotor.k_e: Symbol('kappa',positive=True),
            # dcmotor.k_m: Symbol('kappa',positive=True),
            self.k_e: Symbol("kappa", negative=True),
            self.k_m: Symbol("kappa", negative=True),
            self.J: 1,
            self.B: 1,
            t: tau,
            self.M_obc: Symbol("mu", positive=True),
            self.U_z: Symbol("upsilon", negative=True),
        }
        dvars = Matrix([[self.i_w], [self.omega_s]]).subs({t: tau})

        ode = ODESystem(
            odes=(self._eoms.subs(syms_dict)) * -1, dvars=dvars, ode_order=1, ivar=tau
        )

        # return FirstOrderLinearODESystem.from_ode_system(ode)
        return ode


class DCMotorIIOrder(ComposedSystem):
    scheme_name = "DC_motor.png"

    U_z = Symbol("U_z", positive=True)
    R_w = Symbol("R_w", positive=True)
    L_w = Symbol("L_w", positive=True)
    E = Symbol("E", positive=True)
    U_Rw = Symbol("U_Rw", positive=True)
    U_Lw = Symbol("U_Lw", positive=True)
    k_e = Symbol("k_e", positive=True)
    M_s = Symbol("M_s", positive=True)
    B = Symbol("B", positive=True)
    J = Symbol("J", positive=True)
    M_obc = Symbol("M_l", positive=True)
    k_m = Symbol("k_m", positive=True)
    M_a = Symbol("M_a", positive=True)
    M_r = Symbol("M_r", positive=True)
    i_w = dynamicsymbols("i_w")
    omega_s = dynamicsymbols("omega_s")
    charge = dynamicsymbols("q_1")
    phi = dynamicsymbols("varphi")

    def __init__(
        self,
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
        ivar=Symbol("t"),
        **kwargs
    ):

        if U_z is not None:
            self.U_z = U_z
        if R_w is not None:
            self.R_w = R_w
        if L_w is not None:
            self.L_w = L_w
        if E is not None:
            self.E = E
        if U_Rw is not None:
            self.U_Rw = U_Rw
        if U_Lw is not None:
            self.U_Lw = U_Lw
        if k_e is not None:
            self.k_e = k_e
        if M_s is not None:
            self.M_s = M_s
        if B is not None:
            self.B = B
        if J is not None:
            self.J = J
        if M_obc is not None:
            self.M_obc = M_obc
        if k_m is not None:
            self.k_m = k_m
        if M_a is not None:
            self.M_a = M_a
        if M_r is not None:
            self.M_r = M_r

        self.ivar = ivar
        self.qs = [self.charge, self.phi]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.resistance = Resistor(-self.R_w, self.charge, qs=self.qs, ivar=self.ivar)
        self.inductor = Inductor(
            -self.L_w, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.inertia = Inductor(
            -self.J, self.phi, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.viscous_friction = Resistor(-self.B, self.phi, qs=self.qs, ivar=self.ivar)
        #         self.resistor = CurrentDependentResistor(-self.R_w, self.i_w,  qs=self.qs, ivar=self.ivar , frame=base_frame)
        # #         self.inductor = Resistor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
        #         self.inductor = Inductor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(
            -self.U_z, self.charge, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.electromotive_force = ElectromotiveForce(
            +self.k_e * self.phi.diff(self.ivar),
            self.charge,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        self.engine_load_torque = EngineLoadTorque(
            self.M_obc, self.phi, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.rotor_torque = RotorTorque(
            -self.k_m * self.charge.diff(self.ivar),
            self.phi,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        #         self.rotor_angular_acceleration_torque = RotorAngularAccelerationTorque(-self.J, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        #         self.rotor_resistance_movement_torque = ResistanceMovementTorque(-self.B, self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)

        components["resistor"] = self.resistance
        components["inductor"] = self.inductor
        components["voltage_source"] = self.voltage_source
        components["inertia"] = self.inertia
        components["viscous_friction"] = self.viscous_friction
        components["electromotive_force"] = self.electromotive_force
        components["engine_load_torque"] = self.engine_load_torque
        components["rotor_torque"] = self.rotor_torque
        #         components['rotor_angular_acceleration_torque'] = self.rotor_angular_acceleration_torque
        #         components['rotor_resistance_movement_torque'] = self.rotor_resistance_movement_torque

        return components

    def symbols_description(self):
        i_roc = self.i_w.diff(self.ivar)
        self.sym_desc_dict = {
            self.U_z: r"voltage supplying the rotor",
            self.R_w: r"equivalent resistance of the rotor windings",
            self.L_w: r"equivalent inductance of the rotor windings",
            self.E: r"electromotive force of induction",
            self.U_Rw: r"voltage across the rotor winding resistance",
            self.U_Lw: r"voltage related to the rotor inductance",
            self.k_e: r"back electromotive force constant",
            self.M_s: r"rotor torque",
            self.B: r"coefficient of viscous friction reduced to the rotor shaft",
            self.J: r"moment of inertia reduced to the rotor shaft",
            self.M_obc: r"engine load torque",
            self.k_m: r"torque constant",
            self.i_w: r"rotor winding current",
            i_roc: r"rotor winding current rate of change",
            self.omega_s: r"angular velocity of the rotor",
            self.M_a: r"rotor angular acceleration torque",
            self.M_r: r"rotor motion resistance torque",
            self.omega_s.diff(self.ivar): r"angular acceleration of the rotor",
        }
        return self.sym_desc_dict

    def get_default_data(self):
        default_data_dict = {
            self.U_z: 10,
            self.R_w: 2,
            self.L_w: 0.1,
            self.k_e: 0.1,
            self.k_m: 0.1,
            self.J: 0.1,
            self.B: 0.5,
            self.M_obc: 0.2,
        }
        return default_data_dict

    def _dimensionless_ode(self):

        from ...solvers.linear import FirstOrderLinearODESystem, ODESystem

        t = self.ivar
        tau = Symbol("tau")

        syms_dict = {
            # dcmotor.U_z: 0,
            self.R_w: 1,
            self.L_w: 1,
            # dcmotor.k_e: Symbol('kappa',positive=True),
            # dcmotor.k_m: Symbol('kappa',positive=True),
            self.k_e: Symbol("kappa", negative=True),
            self.k_m: Symbol("kappa", negative=True),
            self.J: 1,
            self.B: 1,
            t: tau,
            # dcmotor.M_obc: 10,
        }
        dvars = Matrix([[self.i_w], [self.omega_s]]).subs({t: tau})

        ode = ODESystem(
            odes=(self._eoms.subs(syms_dict)) * -1, dvars=dvars, ode_order=1, ivar=tau
        )

        # return FirstOrderLinearODESystem.from_ode_system(ode)
        return ode



class DCMotorHeaviside(DCMotor):
    scheme_name = "DC_motor.png"

    U_z = Symbol("U_z", positive=True)
    R_w = Symbol("R_w", positive=True)
    L_w = Symbol("L_w", positive=True)
    E = Symbol("E", positive=True)
    U_Rw = Symbol("U_Rw", positive=True)
    U_Lw = Symbol("U_Lw", positive=True)
    k_e = Symbol("k_e", positive=True)
    M_s = Symbol("M_s", positive=True)
    B = Symbol("B", positive=True)
    J = Symbol("J", positive=True)
    M_obc = Symbol("M_l", positive=True)
    k_m = Symbol("k_m", positive=True)
    M_a = Symbol("M_a", positive=True)
    M_r = Symbol("M_r", positive=True)
    i_w = dynamicsymbols("i_w")

    omega_s = dynamicsymbols("omega_s")
    T = Symbol("T", positive=True)

    def __init__(
        self,
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
        T=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if U_z is not None:
            self.U_z = U_z
        if R_w is not None:
            self.R_w = R_w
        if L_w is not None:
            self.L_w = L_w
        if E is not None:
            self.E = E
        if U_Rw is not None:
            self.U_Rw = U_Rw
        if U_Lw is not None:
            self.U_Lw = U_Lw
        if k_e is not None:
            self.k_e = k_e
        if M_s is not None:
            self.M_s = M_s
        if B is not None:
            self.B = B
        if J is not None:
            self.J = J
        if M_obc is not None:
            self.M_obc = M_obc
        if k_m is not None:
            self.k_m = k_m
        if M_a is not None:
            self.M_a = M_a
        if M_r is not None:
            self.M_r = M_r
        if T is not None:
            self.T = T

        self.ivar = ivar
        self.qs = [self.i_w, self.omega_s]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self.resistor = CurrentDependentResistor(
            self.R_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.inductor = Inductor(
            self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        #         self.inductor = Inductor(-self.L_w, self.i_w, qs=self.qs, ivar=self.ivar, frame=base_frame)
        self.voltage_source = VoltageSource(
            -self.U_z * Heaviside(sin(2 * pi * self.ivar / self.T)),
            self.i_w,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        #         self.voltage_source = VoltageSource(-self.U_z*(Heaviside(self.ivar)-Heaviside(self.ivar-2)+Heaviside(self.ivar-4)-Heaviside(self.ivar-6)+Heaviside(self.ivar-8)), self.i_w, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.electromotive_force = ElectromotiveForce(
            -self.k_e * self.omega_s,
            self.i_w,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        self.engine_load_torque = EngineLoadTorque(
            self.M_obc * Heaviside(sin(2 * pi * self.ivar / self.T - pi / 2)),
            self.omega_s,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        #         self.engine_load_torque = EngineLoadTorque(self.M_obc*(Heaviside(self.ivar-1)-Heaviside(self.ivar-3)+Heaviside(self.ivar-5)-Heaviside(self.ivar-7)+Heaviside(self.ivar-9)), self.omega_s, qs=self.qs, ivar=self.ivar , frame=base_frame)
        self.rotor_torque = RotorTorque(
            self.k_m * self.i_w,
            self.omega_s,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )
        self.rotor_angular_acceleration_torque = RotorAngularAccelerationTorque(
            self.J, self.omega_s, qs=self.qs, ivar=self.ivar, frame=base_frame
        )
        self.rotor_resistance_movement_torque = ResistanceMovementTorque(
            self.B, self.omega_s, qs=self.qs, ivar=self.ivar, frame=base_frame
        )

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor
        components["voltage_source"] = self.voltage_source
        components["electromotive_force"] = self.electromotive_force
        components["engine_load_torque"] = self.engine_load_torque
        components["rotor_torque"] = self.rotor_torque
        components["rotor_angular_acceleration_torque"] = (
            self.rotor_angular_acceleration_torque
        )
        components["rotor_resistance_movement_torque"] = (
            self.rotor_resistance_movement_torque
        )

        return components

    def symbols_description(self):
        i_roc = self.i_w.diff(self.ivar)
        self.sym_desc_dict = {
            self.U_z: r"voltage supplying the rotor",
            self.R_w: r"equivalent resistance of the rotor windings",
            self.L_w: r"equivalent inductance of the rotor windings",
            self.E: r"electromotive force of induction",
            self.U_Rw: r"voltage across the rotor winding resistance",
            self.U_Lw: r"voltage related to the rotor inductance",
            self.k_e: r"electric constant",
            self.M_s: r"rotor torque",
            self.B: r"coefficient of viscous friction reduced to the rotor shaft",
            self.J: r"moment of inertia reduced to the rotor shaft",
            self.M_obc: r"motor load torque",
            self.k_m: r"mechanical constant",
            self.i_w: r"rotor winding current",
            i_roc: r"rotor winding current rate of change",
            self.omega_s: r"angular velocity of the rotor",
            self.M_a: r"rotor angular acceleration torque",
            self.M_r: r"rotor motion resistance torque",
            self.omega_s.diff(self.ivar): r"angular acceleration of the rotor",
            self.T: "Heaviside cycle length",
            Symbol("\Theta"): "Heaviside function",
        }
        return self.sym_desc_dict


class CurrentDependentResistor(Spring):
    """
    Model of Resistor

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """

    def __init__(self, resistance, q0, qs=None, ivar=Symbol("t"), frame=base_frame):
        super().__init__(
            stiffness=resistance,
            pos1=q0,
            pos2=0,
            qs=qs,
            ivar=Symbol("t"),
            frame=base_frame,
        )


class ElectromotiveForce(Force):
    """
    Model of a Electromotive Force
    """

    def __init__(
        self, electric_constant, q0, qs=None, ivar=Symbol("t"), frame=base_frame
    ):
        super().__init__(
            force=electric_constant, pos1=q0, qs=qs, ivar=ivar, frame=frame
        )


class EngineLoadTorque(Force):
    """
    Model of a engine load torque
    """

    def __init__(self, torque, q0, qs=None, ivar=Symbol("t"), frame=base_frame):
        super().__init__(force=torque, pos1=q0, qs=qs, ivar=ivar, frame=frame)


class RotorTorque(Force):
    """
    Model of a Electromotive Force
    """

    def __init__(
        self, mechanical_constant, q0, qs=None, ivar=Symbol("t"), frame=base_frame
    ):
        super().__init__(
            force=mechanical_constant, pos1=q0, qs=qs, ivar=ivar, frame=frame
        )


class RotorAngularAccelerationTorque(Damper):
    """
    Model of Resistor

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """

    def __init__(
        self, moment_of_interia, q0, qs=None, ivar=Symbol("t"), frame=base_frame
    ):
        super().__init__(
            c=moment_of_interia,
            pos1=q0,
            pos2=0,
            qs=qs,
            ivar=Symbol("t"),
            frame=base_frame,
        )


class ResistanceMovementTorque(Spring):
    """
    Model of a Electromotive Force
    """

    def __init__(
        self,
        friction_coefficient,
        q0,
        q1=0,
        qs=None,
        ivar=Symbol("t"),
        frame=base_frame,
    ):
        super().__init__(
            stiffness=friction_coefficient,
            pos1=q0,
            pos2=q1,
            l_0=0,
            qs=qs,
            ivar=ivar,
            frame=frame,
        )
class DCMotorWithBattery(ComposedSystem):

    scheme_name = "thevenincircuit.png"
    real_name = "liioncell.PNG"

    R_1 = Symbol("R_1", positive=True)
    R_2 = Symbol("R_2", positive=True)

    L_1 = Symbol("L_1", positive=True)
    L_2 = Symbol("L_2", positive=True)

    C = Symbol("C", positive=True)
    U = Symbol("U", positive=True)
    q_1 = dynamicsymbols("q_1")
    q_2 = dynamicsymbols("q_2")
    t = Symbol("t")
    U_li = Function("U_li")(t)
    U_oc = Symbol("U_oc")
    R_0 = Symbol("R_0")
    I_li = Function("I_zad")(t)
    R_th = Symbol("R_th")
    C_th = Symbol("C_th")
    SOC = Function("SOC")(t)
    SOC_init = Symbol("SOC_init")
    C_rated = Symbol("C_rated")
    t_0 = Symbol("t_0")
    q_0 = Symbol("q_0")
    U_th = Function("U_th")(t)
    funcI = Function("funcI")(t)
    # U_th = Symbol('U_th')

    def __init__(
        self,
        R_1=None,
        R_2=None,
        C=None,
        U=None,
        q_1=None,
        q_2=None,
        t=None,
        U_th=None,
        R_th=None,
        C_th=None,
        L_1=None,
        L_2=None,
        I_li=None,
        q_0=None,
        ivar=Symbol("t"),
        **kwargs,
    ):

        if t is not None:
            self.t = t
        if R_1 is not None:
            self.R_1 = R_1
        if R_2 is not None:
            self.R_2 = R_2
        if C is not None:
            self.C = C
        if U is not None:
            self.U = U
        if q_1 is not None:
            self.q_1 = q_1
        if q_2 is not None:
            self.q_2 = q_2
        if U_th is not None:
            self.U_th = U_th
        if R_th is not None:
            self.R_th = R_th
        if C_th is not None:
            self.C_th = C_th
        if L_1 is not None:
            self.L_1 = L_1
        if L_2 is not None:
            self.L_2 = L_2
        if I_li is not None:
            self.I_li = I_li
        if q_0 is not None:
            self.q_0 = q_0
        self.qs = [self.q_1, self.q_2]
        self._init_from_components(**kwargs)

    @property
    def components(self):
        components = {}

        self._battery = BatteryCell(R_1=self.R_1, R_2=self.R_2, C=self.C, U=self.U, q_1=self.q_1, q_2=self.q_2, t=self.t, U_th=self.U_th, R_th=self.R_th, C_th=self.C_th, L_1=self.L_1, L_2=self.L_2, I_li=self.I_li)
        self._motor = DCMotorIIOrder(qs=[self.q_1])


        components["battery"] = self._battery
        components["motor"] = self._motor


        return components
    def soc(self):
        soc=(self.q_0-integrate(self.qs[0]))/self.q_0
        return soc
    def get_default_data(self):
        # Konfiguracja: Pakiet Li-Ion 3S (12.6V max) + Mocny silnik DC
        default_data_dict = {
            # --- BATERIA (Model Thevenina 2RC) ---
            self.R_1: 0.015,       # [Ohm] Rezystancja dynamiczna (szybka) - typowa dla pakietu
            self.R_2: 0.030,       # [Ohm] Rezystancja dynamiczna (wolna)
            self.C: 500.0,         # [F] Pojemność w gałęzi RC (modeluje relaksację, nie pojemność główną!)
            
            # Indukcyjności pasożytnicze (małe, ale niezerowe dla numeryki)
            self.L_1: 1e-6,        # [H] 1 uH
            self.L_2: 1e-6,        # [H] 1 uH
            
            # Parametry główne baterii
            self.U_oc: 12.6,       # [V] Napięcie jałowe (Open Circuit) dla pełnego naładowania (3x4.2V)
            self.U_th: 12.6,       # [V] Początkowe napięcie na kondensatorach (stan równowagi)
            self.U: 12.0,          # [V] Napięcie nominalne pod obciążeniem
            self.R_0: 0.15,        # [Ohm] Rezystancja wewnętrzna całego pakietu (szeregowa)
            
            # Termika baterii
            self.R_th: 2.0,        # [K/W] Rezystancja termiczna (chłodzenie pasywne jest słabsze niż 0.1)
            self.C_th: 1000.0,     # [J/K] Pojemność cieplna (masa ok. 1-2 kg)
            
            # Stan baterii
            self.SOC_init: 0.95,   # [0-1] Prawie pełna
            self.C_rated: 18000,   # [As] 5 Ah * 3600 = 18000 As (Realna pojemność małego pakietu)

            # --- SILNIK DC (Obciążenie) ---
            # Parametry dobrane tak, by silnik mógł ruszyć pod obciążeniem
            DCMotorIIOrder.U_z: 12.0,    # [V] Napięcie zasilania (nominalne)
            
            DCMotorIIOrder.R_w: 0.25,    # [Ohm] Rezystancja uzwojeń (Stall Current ~48A)
            DCMotorIIOrder.L_w: 0.002,   # [H] Indukcyjność (typowa dla silników tej mocy)
            
            # Stałe silnika (k_e * w = U_emf)
            # Dla k=0.1: przy 12V max prędkość to ~120 rad/s (1150 RPM)
            DCMotorIIOrder.k_e: 0.1,     # [V*s/rad]
            DCMotorIIOrder.k_m: 0.1,     # [N*m/A] (W układzie SI k_m = k_e)
            
            # Mechanika
            DCMotorIIOrder.J: 0.02,      # [kg*m^2] Bezwładność wirnika + przekładni
            DCMotorIIOrder.B: 0.05,      # [N*m*s/rad] Współczynnik tarcia lepkiego
            
            # Obciążenie
            # Moment startowy silnika: M = k_m * (U/R) = 0.1 * (12/0.25) = 4.8 Nm
            # Obciążenie musi być mniejsze niż 4.8 Nm, żeby ruszył!
            DCMotorIIOrder.M_obc: 2.5    # [N*m] Obciążenie nominalne (np. jazda pod górkę)
        }
        return default_data_dict

