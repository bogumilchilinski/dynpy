from numbers import Number

from sympy import (
    Derivative,
    Eq,
    Expr,
    Heaviside,
    Matrix,
    S,
    Symbol,
    cos,
    diag,
    diff,
    sin,
    sqrt,
    symbols,
)
from sympy.physics.mechanics import Point, ReferenceFrame, dynamicsymbols
from sympy.physics.vector import vlatex, vpprint

from ..elements import Damper, Force, MaterialPoint, Spring, base_frame, base_origin


class Capacitor(Spring):
    """
    Model of a Capacitor

    Creates a singular model, after inputing correct values of capacity - c and general coordinate(s), which analytically display the dynamics of displacing spring after cummulating PE.
    """

    scheme_name = "spring.png"
    real_name = "spring.png"

    def __init__(self, capacity, q0, q1=0, qs=None, ivar=Symbol("t"), frame=base_frame):
        super().__init__(
            stiffness=1 / capacity,
            pos1=q0,
            pos2=q1,
            l_0=0,
            qs=qs,
            ivar=ivar,
            frame=frame,
        )


class Inductor(MaterialPoint):
    """
    Model of a Material point with changing point of mass:
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """

    scheme_name = "material_point.png"
    real_name = "material_point.png"

    def __init__(self, inductance, q0, qs=None, frame=base_frame, ivar=Symbol("t")):
        super().__init__(m=inductance, pos1=q0, qs=qs, frame=frame, ivar=ivar)


class Resistor(Damper):
    """
    Model of Resistor

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """

    scheme_name = "material_point.png"
    real_name = "material_point.png"

    def __init__(self, resistance, q0, qs=None, ivar=None, frame=base_frame):
        """
        Linear electric resistor – analogia tłumika lepkiego.
        Jeżeli użytkownik nie poda `ivar`, przyjmujemy Symbol('t').
        """
        ivar = ivar or Symbol("t")
        super().__init__(c=resistance, pos1=q0, pos2=0, qs=qs, ivar=ivar, frame=frame)


class VoltageSource(Force):
    """
    Model of a voltage source

    Creates a singular model, after inputing correct values of capacity - c and general coordinate(s), which analytically display the dynamics of displacing spring after cummulating PE.
    """

    scheme_name = "voltgesrc.png"
    real_name = "spring.png"

    def __init__(self, voltage, q0, qs=None, ivar=Symbol("t"), frame=base_frame):
        super().__init__(force=voltage, pos1=q0, qs=qs, ivar=ivar, frame=frame)


class EMForce(Force):
    """
    Creates enforcement that shows the effect of EM coupling.
    """

    scheme_name = ""
    real_name = ""

    def __init__(self, voltage, q0, qs=None, ivar=Symbol("t"), frame=base_frame):
        super().__init__(force=voltage, pos1=q0, qs=qs, ivar=ivar, frame=frame)
