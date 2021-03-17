from sympy import *  # S, Symbol, diff
<<<<<<< HEAD
from dynamics import LagrangesDynamicSystem
=======
from ..dynamics import LagrangesDynamicSystem
>>>>>>> 23b3203dd953fa3867db922b5cb77ec6203508e9
from sympy.physics.mechanics import *
from sympy.physics.vector import *

class MaterialPoint(LagrangesDynamicSystem):
    """
    Model of a Material point with changing point of mass:
    """
    def __init__(self, m, pos1=0, pos2=0, pos_c=0, qs=None,  ivar=Symbol('t')):
        
        if pos1 == 0 and pos2 == 0:
            
            self.qs = [pos_c]

            Lagrangian = S.Half * m * (diff(pos_c,ivar))**2
        
        elif pos_c == 0:
            
            if qs == qs:
                qs = [pos1]
            else:
                qs = [pos1, pos2]
            
            Lagrangian = S.Half * m * (diff(pos1,ivar) + diff(pos2,ivar))**2
            
            """Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Spring(LagrangesDynamicSystem):
    """
    Model of a Spring:
    """
    def __init__(self, k, pos1=0, pos2=0, pos_c=0, qs=None, ivar=Symbol('t')):
            
        if pos1 == 0 and pos2 == 0:
            
            self.qs = [pos_c]

            Lagrangian = -S.Half * k * (pos_c)**2
        
        elif pos_c == 0:
            
            if qs == qs:
                qs = [pos1]
            else:
                qs = [pos1, pos2]
            
            Lagrangian = -S.Half * k * (pos1 - pos2)**2
            
            """
            Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after            cummulating PE.
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Centroid(LagrangesDynamicSystem):
    """
    Model of a changing centroid for potential energy:
    """
    def __init__(self, m, g, pos1=0, pos_c=0, qs=None, ivar=Symbol('t')):
        
        if pos1 == 0:
            
            self.qs = [pos_c]

            Lagrangian = -(m * g * (pos_c))
        
        elif pos_c == 0:
            
            qs = [pos1]
            
            Lagrangian = -(m * g * (pos1))
            
            """
            Creates a singular model, after inputing correct values of gravity field - g, mass of - m as well as additionaly the general coordiante
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Disk(LagrangesDynamicSystem):
    """
    Model of a Disk:
    """
    def __init__(self, I, pos1=0, pos_c=0, qs=None, ivar=Symbol('t')):
        
        if pos1 == 0:
            
            self.qs = [pos_c]

            Lagrangian = S.Half * I * diff(pos_c,ivar)**2
        
        elif pos_c == 0:
            
            qs = [pos1]
            
            Lagrangian = S.Half * I * diff(pos1,ivar)**2
            
            """
            Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a                 rotating wheel.
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))
        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class RigidBody2D(LagrangesDynamicSystem):
    """
    Model of a 2DoF Rigid body:
    """
    def __init__(self, m, I, pos_lin=0, pos_rot=0, pos_lin_c=0, pos_rot_c=0, qs=None, ivar=Symbol('t')):
        
        if pos_lin == 0 and pos_rot == 0:
            
            self.qs = [pos_lin_c, pos_rot_c]

            Lagrangian = S.Half * m * diff(pos_lin_c,ivar)**2 + S.Half * I * diff(pos_rot_c,ivar)**2
        
        elif pos_lin_c == 0 and pos_rot_c == 0:
            
            qs = [pos_lin, pos_rot]
            
            Lagrangian = S.Half * m * diff(pos_lin,ivar)**2 + S.Half * I * diff(pos_rot,ivar)**2
            
            """
                Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which             analytically display the dynamics of rotating and propeling object: beam.
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))
        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)
  
        
class Pendulum(LagrangesDynamicSystem):
    """
    Model of a sDoF mathematical Pendulum:
    """
    def __init__(self, m,  g, l, pos1=0, trig=0, pos_c=0, qs=None, ivar=Symbol('t')):

        
        if trig == 0:
            
            self.qs = [pos_c]
            
            if pos1 == 0:
                P = pos_c
            else:
                P = pos1

            Lagrangian = S.Half * m * l**2 * diff(P,ivar)**2 - m * g * l * (pos_c)
        
        elif pos_c == 0:
            
            qs = [pos1]
            
            if trig == cos(pos1):
                trig = cos(pos1)
            elif trig == sin(pos1):
                trig = sin(pos1)
            else:
                Lagrangian = 0
            
            Lagrangian = S.Half * m * l**2 * diff(pos1,ivar)**2 - m * g * l * (1-trig)
            
            """
            Creates a singular model, after inputing correct values of mass - m , gravitational field - g, length of a strong - l and general coordinate which estabilshes an analytical display of a mathematical model of a sDoF pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.
            """
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Damper(LagrangesDynamicSystem):
    """
    Model of a Damper:

    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional to the velocity, in time, of an inertial virbating element.
    """
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=ReferenceFrame('N')):
        

        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1, pos2]
        
        dpos1 = diff(pos1, ivar)
        dpos2 = diff(pos2, ivar)
        
        Pd = Point('Pd')
        Pd.set_vel(frame, (dpos1 - dpos2) * frame.x)
        
        D = (((S.Half) * c * (dpos1 - dpos2)**2).diff(dpos1))
        
        forcelist = [(Pd, -D*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class PID(LagrangesDynamicSystem):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    def __init__(self, kp, ki, kd, pos1, qs=None, ivar=Symbol('t'), frame=ReferenceFrame('N')):

        if qs == None:
            qs = [pos1]

        dpos1 = diff(pos1, ivar)
        d2pos1 = diff(dpos1, ivar)
        
        Pp = Point('Pp')
        Pp.set_vel(frame, dpos1 * frame.x)
        
        In = pos1  # In - error
        u = pos1*ki + dpos1*kp + d2pos1*kd
        G = In - u
        
        forcelist = [(Pp, G*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class Excitation(LagrangesDynamicSystem):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    def __init__(self, f, pos_rot, ivar=Symbol('t'), frame=ReferenceFrame('N')):
        
        qs = [pos_rot]
        
        dpos_rot = diff(pos_rot,ivar)
        
        Pe = Point('Pe')
        Pe.set_vel(frame, dpos_rot * frame.x)
        Pe.vel(frame)

        F = f * cos(dpos_rot*ivar) + f * sin(dpos_rot*ivar)
        
        forcelist = [(Pe, F*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

class Force(LagrangesDynamicSystem):
    """
    Creates enforcement.
    """
    def __init__(self,
                 F,
                 pos1 = None,
                 qs = None,
                 velocity= None,
                 ivar=Symbol('t'),
                 frame=None):

        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        if isinstance(pos1, Point):
            P = pos1
            P.set_vel(frame, pos1.vel(frame))
            
        else:
            frame=ReferenceFrame("N")
            P = Point('P')
            pos1 = pos1        
            dpos1 = diff(pos1, ivar)
            P.set_vel(frame, dpos1 * frame.x)

        forcelist = [(P, F* frame.x), (P, F * frame.y), (P, F*frame.z)]

        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

######################################################################################################################################

class SpringFrame(LagrangesDynamicSystem):
    """
    Model of a Spring based on a Particle & Frame instances
    """

    def __init__(self, k, pos1, pos2, qs=None, ivar=Symbol('t'), frame=ReferenceFrame('N')):
        
        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1, pos2]
        
        P = Point('P')
        P.set_vel(frame, (1) * frame.x)
        
        Pa = Particle('Pa', P, 1)
        
        Pa.potential_energy = S.Half * k * pos1.pos_from(pos2).magnitude()**2
        L = Lagrangian(frame, Pa)

        super().__init__(Lagrangian=L, qs=qs, frame=frame, ivar=ivar)

