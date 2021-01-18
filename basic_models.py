from sympy import S, Symbol
from .dynpy import LagrangesDynamicSystem


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
    def __init__(self, stiffness, pos1, pos2=0, qs=None, ivar=sym.Symbol('t')):
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
                 m,
                 r,
                 I,
                 pos1,
                 qs=None,
                 ivar=sym.Symbol('t'),
                 evaluate=True):

        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        I = S.One / 2 * m * r**2
        Lagrangian = (S.One / 2 * (I * (pos1.diff(ivar))**2))
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
                 ivar=sym.Symbol('t'),
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
                 ivar=sym.Symbol('t'),
                 frame=N):

        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1, pos2]

        dpos1 = diff(pos1, t)
        dpos2 = diff(pos2, t)
        N = ReferenceFrame('N')
        P = Point('P')
        P.set_vel(N, (dpos1 - dpos2) * N.x)
        Lagrangian = 0
        forcelist = [
            (P, -((S.One / 2) * c * (dpos1 - dpos2)**2).diff(dpos1) * N.x)
        ]

        super().__init__(0, qs=qs, forcelist=forcelist, frame=N, ivar=ivar)
