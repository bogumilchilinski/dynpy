from sympy import *
from ..dynamics import LagrangesDynamicSystem
from sympy.physics.mechanics import *
from sympy.physics.vector import *
import base64
import IPython as IP

base_frame=ReferenceFrame('N')

class Elements(LagrangesDynamicSystem):
    """Base class for all elements
    """
    @classmethod
    def preview(cls,real=False):
        if real:
            path = __file__.replace('elements.py', 'images/') + cls.real_name
            with open(f"{path}", "rb") as image_file:
                encoded_string = base64.b64encode(image_file.read())
            image_file.close()
            
        else:
            path = __file__.replace('elements.py', 'images/') + cls.scheme_name
            with open(f"{path}", "rb") as image_file:
                encoded_string = base64.b64encode(image_file.read())
            image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))

class MaterialPoint(Elements):
    """
    Model of a Material point with changing point of mass:
    """
    """Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'
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


class Spring(Elements):
    """
    Model of a Spring:
    """
    """
    Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after            cummulating PE.
    """
    scheme_name = 'spring.png'
    real_name = 'spring.png'
    def __init__(self, stiffness, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame = base_frame):
        if not qs:
            qs = [pos1]
        else:
            qs = qs

        if isinstance(pos1, Point): 
            P1 = pos1
            P2 = pos2
            if not qs:
                diffs = P1.vel(frame).magnitude().atoms(Derivative)
                diffs2 = P2.vel(frame).magnitude().atoms(Derivative)
                
                qs= [deriv.args[0] for deriv in diffs]
                qs = np.append(qs, [deriv.args[0] for deriv in diffs2])
                print(qs)

            Pa = Particle('Pa', P1, 1)
            Pa.potential_energy = S.Half * k * P1.pos_from(P2).magnitude()**2

            L = Lagrangian(frame, Pa)

        else:
            L = -S.Half * stiffness * (pos1 - pos2)**2


        super().__init__(Lagrangian=L, qs=qs, ivar=ivar)



class GravitationalForce(Elements):
    """
    Model of a changing centroid for potential energy:
    """
    """
    Creates a singular model, after inputing correct values of gravity field - g, mass of - m as well as additionaly the general coordiante
    """
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'
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


class Disk(Elements):
    """
    Model of a Disk:
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a rotating wheel.
    """
    scheme_name = 'disk.png'
    real_name = 'disk.png'
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

        
class RigidBody2D(Elements):
    """
    Model of a 2DoF Rigid body:
    """
    """
    Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which             analytically display the dynamics of rotating and propeling object: beam.
    """
    scheme_name = 'rigid_body2D.png'
    real_name = 'rigid_body2D.png'
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
  
        
class Damper(Elements):
    """
    Model of a Damper:

    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional to the velocity, in time, of an inertial virbating element.
    """
    scheme_name = 'damper.png'
    real_name = 'damper.png'
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame):
        

        if qs == None:
            qs = [pos1]
        else:
            qs = qs
        
        dpos1 = diff(pos1, ivar)
        dpos2 = diff(pos2, ivar)
        
        P = Point('P')
        P.set_vel(frame, (dpos1 - dpos2) * frame.x)
        
        D = (((S.Half) * c * (dpos1 - dpos2)**2).diff(dpos1))
        
        forcelist = [(P, -D*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class PID(Elements):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'
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

        
class Excitation(Elements):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'
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
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'
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