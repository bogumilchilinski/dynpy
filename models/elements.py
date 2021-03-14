from sympy import S, Symbol, diff,cos,sin
from ..dynpy import LagrangesDynamicSystem
from sympy.physics.mechanics import ReferenceFrame,Point


class MaterialPoint(LagrangesDynamicSystem):
    """
    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion
    """
    def __init__(self, mass, pos, ivar=Symbol('t')):

        qs = [pos]
        Lagrangian = S.Half * mass * pos.diff(ivar)**2

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Spring(LagrangesDynamicSystem):
    """
    Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after       cummulating PE.
    """
    def __init__(self, stiffness, pos1, pos2=0, qs=None, ivar=Symbol('t')):
        self.stiffness = stiffness
        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        Lagrangian = -(S.One / 2 * (stiffness * (pos1 - pos2)**2))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Disk(LagrangesDynamicSystem):
    """
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a             rotating wheel.
    """
    def __init__(self,
                 I,
                 pos1,
                 qs=None,
                 ivar=Symbol('t'),
                 evaluate=True):

        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        Lagrangian = (S.One / 2 * (I * (pos1.diff(ivar))**2))
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Pendulum(LagrangesDynamicSystem):
    """
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a             rotating wheel.
    """
    def __init__(self,
                 m,
                 l,
                 pos1,
                 g=Symbol('g',positive=True),
                 qs=None,
                 ivar=Symbol('t'),
                 evaluate=True):

        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        Lagrangian = (S.One / 2 * (m*l**2 * (pos1.diff(ivar))**2)) + m*g*l*cos(pos1)
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

class RigidBody2d(LagrangesDynamicSystem):
    """
    Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which       analytically display the dynamics of rotating and propeling object: beam.
    """
    def __init__(self,
                 m,
                 I,
                 pos_lin,
                 pos_rot,
                 ivar=Symbol('t'),
                 evaluate=True):

        qs = [pos_lin, pos_rot]
        Lagrangian = S.One / 2 * m * pos_lin.diff(
            ivar)**2 + S.One / 2 * I * pos_rot.diff(ivar)**2
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Damper(LagrangesDynamicSystem):
    """
    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional       to the velocity, in time, of an inertial virbating element.
    """
    def __init__(self,
                 c,
                 pos1,
                 pos2=0,
                 qs=None,
                 ivar=Symbol('t'),
                 frame=ReferenceFrame('N')):

        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1, pos2]

        dpos1 = diff(pos1, ivar)
        dpos2 = diff(pos2, ivar)

        P = Point('P')
        P.set_vel(frame, (dpos1 - dpos2) * frame.x)
        Lagrangian = 0
        forcelist = [
            (P, -((S.One / 2) * c * (dpos1 - dpos2)**2).diff(dpos1) * frame.x)
        ]

        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
