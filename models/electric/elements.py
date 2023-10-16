from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Derivative, Expr)
from numbers import Number
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

from ..elements import MaterialPoint, Damper, Spring , base_frame, base_origin




class Capacitor(Spring):

    """
    Model of a Capacitor

    Creates a singular model, after inputing correct values of capacity - c and general coordinate(s), which analytically display the dynamics of displacing spring after cummulating PE.
    """
    scheme_name = 'spring.png'
    real_name = 'spring.png'

    def __init__(self, capacity, q0,  qs=None,ivar=Symbol('t'), frame = base_frame):
        super().__init__(stiffness=1/capacity, pos1=q0, pos2=0,l_0=0 ,  qs=qs,ivar=ivar, frame = frame)


class Inductor(MaterialPoint):
    """
    Model of a Material point with changing point of mass:
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'

    def __init__(self, inductance, q0 , qs=None, frame=base_frame, ivar=Symbol('t')):
        super().__init__(m=inductance, pos1=q0 , qs=qs, frame=frame, ivar=ivar)

class Resistor(Damper):
    """
    Model of Resistor
    
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'
    def __init__(self, resistance, q0,  qs=None, ivar=Symbol('t'), frame=base_frame):
        super().__init__(c=resistance, pos1=q0, pos2=0, qs=qs, ivar=Symbol('t'), frame=base_frame)