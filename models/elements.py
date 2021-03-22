from sympy import *  # S, Symbol, diff

from ..dynamics import LagrangesDynamicSystem

from sympy.physics.mechanics import *
from sympy.physics.vector import *

base_frame=ReferenceFrame('N')

class MaterialPoint(LagrangesDynamicSystem):
    """
    Model of a Material point with changing point of mass:
    """
    """Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
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
            
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Spring(LagrangesDynamicSystem):
    """
    Model of a Spring:
    """
    """
    Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after            cummulating PE.
    """
    def __init__(self, stiffness, pos1, pos2=0, qs=None, ivar=Symbol('t')):
        self.stiffness = stiffness
        if qs == None:
            qs = [pos1]
        else:
            qs = qs

        Lagrangian = -(S.One / 2 * (stiffness * (pos1 - pos2)**2))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)



class Centroid(LagrangesDynamicSystem):
    """
    Model of a changing centroid for potential energy:
    """

    """
    Creates a singular model, after inputing correct values of gravity field - g, mass of - m as well as additionaly the general coordiante
    """
    def __init__(self, m, g, pos1=0, pos_c=0, qs=None, ivar=Symbol('t')):
        
        if pos1 == 0:
            
            self.qs = [pos_c]

            Lagrangian = -(m * g * (pos_c))
        
        elif pos_c == 0:
            
            qs = [pos1]
            
            Lagrangian = -(m * g * (pos1))

        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Disk(LagrangesDynamicSystem):
    """
    Model of a Disk:
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a rotating wheel.
    """

    def __init__(self, I, pos1=0, pos_c=0, qs=None, ivar=Symbol('t')):
        
        if pos1 == 0:
            
            self.qs = [pos_c]

            Lagrangian = S.Half * I * diff(pos_c,ivar)**2
        
        elif pos_c == 0:
            
            qs = [pos1]
            
            Lagrangian = S.Half * I * diff(pos1,ivar)**2
            

        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))
        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class RigidBody2D(LagrangesDynamicSystem):
    """
    Model of a 2DoF Rigid body:
    """
    """
    Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which             analytically display the dynamics of rotating and propeling object: beam.
    """
    def __init__(self, m, I, pos_lin=0, pos_rot=0, pos_lin_c=0, pos_rot_c=0, qs=None, ivar=Symbol('t')):
        
        if pos_lin == 0 and pos_rot == 0:
            
            self.qs = [pos_lin_c, pos_rot_c]

            Lagrangian = S.Half * m * diff(pos_lin_c,ivar)**2 + S.Half * I * diff(pos_rot_c,ivar)**2
        
        elif pos_lin_c == 0 and pos_rot_c == 0:
            
            qs = [pos_lin, pos_rot]
            
            Lagrangian = S.Half * m * diff(pos_lin,ivar)**2 + S.Half * I * diff(pos_rot,ivar)**2
            
        else:
            
            Lagrangian = 0
            print(str('YOU CANNOT APPLY BOTH METHODS AT ONCE'))
        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)
  
        
class Pendulum(LagrangesDynamicSystem):
    """
    Model of a sDoF mathematical Pendulum:
    """
    """
            Creates a singular model, after inputing correct values of mass - m , gravitational field - g, length of a strong - l and general coordinate which estabilshes an analytical display of a mathematical model of a sDoF pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.
    """
    def __init__(self, m,  g, l, angle=0, qs=None, ivar=Symbol('t')):


        if qs == None:
            qs = [angle]
        else:
            qs = qs


        Lagrangian = S.Half * m * l**2 * diff(angle,ivar)**2 - m * g * l * (1-cos(angle))



        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Damper(LagrangesDynamicSystem):
    """
    Model of a Damper:

    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional to the velocity, in time, of an inertial virbating element.
    """
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame):
        

        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1, pos2]
        
        dpos1 = diff(pos1, ivar)
        dpos2 = diff(pos2, ivar)
        
        P = Point('P')
        P.set_vel(frame, (dpos1 - dpos2) * frame.x)
        
        D = (((S.Half) * c * (dpos1 - dpos2)**2).diff(dpos1))
        
        forcelist = [(P, -D*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class PID(LagrangesDynamicSystem):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    def __init__(self, kp, ki, kd, pos1, qs=None, ivar=Symbol('t'), frame=base_frame):

        if qs == None:
            qs = [pos1]

        dpos1 = diff(pos1, ivar)
        d2pos1 = diff(dpos1, ivar)
        
        P = Point('P')
        P.set_vel(frame, dpos1 * frame.x)
        
        In = pos1  # In - error
        u = pos1*ki + dpos1*kp + d2pos1*kd
        G = In - u
        
        forcelist = [(P, G*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class Excitation(LagrangesDynamicSystem):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    def __init__(self, f, pos_rot, ivar=Symbol('t'), frame=base_frame):
        
        qs = [pos_rot]
        
        dpos_rot = diff(pos_rot,ivar)
        
        P = Point('P')
        P.set_vel(frame, dpos_rot * frame.x)
        P.vel(frame)

        F = f * cos(dpos_rot*ivar) + f * sin(dpos_rot*ivar)
        
        forcelist = [(P, F*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

class Force(LagrangesDynamicSystem):
    """
    Creates enforcement.
    """
    def __init__(self,
                 force,
                 pos1 = None,
                 qs = None,
                 ivar=Symbol('t'),
                 frame = base_frame):

        if not qs == None:
            qs = [pos1]
        else:
            qs = qs
    
        if isinstance(pos1, Point):
            P = pos1
            if not qs:
                diffs=P.vel(frame).magnitude().atoms(Derivative)
                
                qs= [deriv.args[0] for deriv in diffs]
                print(qs)
        else:
            P = Point('P')
            P.set_vel(frame,pos1.diff(ivar)*frame.x)
            force=force*frame.x


        forcelist = [(P, force)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

######################################################################################################################################