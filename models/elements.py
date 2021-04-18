from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Point, Derivative, Number, Expr)
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator

import base64
import IPython as IP

base_frame=ReferenceFrame('N')
base_origin=Point('O')


class GeometryOfPoint:
    def __init__(self, *args, frame=base_frame , ivar=Symbol('t')):

        if isinstance(args[0],Point):
            self._point=args[0]

        if isinstance(args[0],Number) or isinstance(args[0],Expr):
            P = Point('P')
            P.set_pos(base_origin, frame.x*args[0])
            P.set_vel(frame, frame.x*diff(args[0], ivar))
            self._point=P

        else:
            print('Unsupported data type')
            P = Point('P')
            P.set_pos(base_origin, frame.x*0)
            P.set_vel(frame, frame.x*diff(0, ivar))
            self._point=P
            


    def get_point(self):

        return self._point



class Element(LagrangesDynamicSystem):
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

class MaterialPoint(Element):
    """
    Model of a Material point with changing point of mass:
    """
    """Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'
    def __init__(self, m, pos1, qs=None,  ivar=Symbol('t')):
        
        if not qs:
            
            self.qs = [pos1]

        Lagrangian = S.Half * m * (pos1.diff(ivar))**2


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class Spring(Element):
    """
    Model of a Spring:
    """
    """
    Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after            cummulating PE.
    """
    scheme_name = 'spring.png'
    real_name = 'spring.png'
    def __init__(self, stiffness, pos1, pos2=0,  qs=None, l0 = 0, ivar=Symbol('t'), frame = base_frame):
        if not qs:
            qs = [pos1]

        pos1=GeometryOfPoint(pos1).get_point()
        pos2=GeometryOfPoint(pos2).get_point()

        if isinstance(pos1,Point):
            u = pos1.pos_from(pos2).magnitude()-l0
            
            L = -S.Half * stiffness * (u**2).expand().simplify()
            
        else:
            L = -S.Half * stiffness * (pos1 - pos2)**2

        super().__init__(Lagrangian=L, qs=qs, ivar=ivar)

        
class NonlinSpring__RefFrme_Pt(Element):
    """
    Model of a Nonlinear Spring with whole ReferenceFrame, Point packed in the class  -- Please Advise! ... qs = [ ] must be supplied in the calling sequence:
    """
    """
    Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s), which analytically display the dynamics of displacing spring after            cummulating PE. The class work in two modes. In coordinate mode and Point mode where the Point mode may be supplemented with parameter.
    """

    def __init__(self, k, l_0, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame):

        if pos1 == qs or pos2 == qs:

            if qs == None:
                self.qs = [pos1]
            elif not pos2 == 0:
                qs = [pos1,pos2]

        else:
            pos1 is type(Point) or pos2 is type(Point)

        P1 = Point('P1')
        P1.set_pos(P1 , frame.x*pos1)

        P2 = Point('P2')
        P2.set_pos(P1, frame.x*pos1 + frame.y*pos2)

        L = - S.Half * k *  (P1.pos_from(P2).magnitude()-l_0)**2

        super().__init__(Lagrangian=L, qs=qs, ivar=ivar, frame=frame)

        
class GravitationalForce(Element):
    """
    Model of a changing centroid for potential energy:
    """
    """
    Creates a singular model, after inputing correct values of gravity field - g, mass of - m as well as additionaly the general coordiante
    """
    scheme_name = ''
    real_name = ''
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


class Disk(Element):
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

        
class RigidBody2D(Element):
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
  
        
class Damper(Element):
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
        P.set_vel(frame, 1 * frame.x)
        
        D = ((S.Half) * c * (dpos1 - dpos2)**2)
        
        
             
        forcelist = [(P, -D.diff(ivar)*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)

        
class PID(Element):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    scheme_name = ''
    real_name = ''
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

        
class Excitation(Element):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    scheme_name = ''
    real_name = ''
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
    scheme_name = ''
    real_name = ''
    def __init__(self,
                 force,
                 pos1 = None,
                 qs = None,
                 ivar=Symbol('t'),
                 frame = base_frame):

        if qs == None:
            qs = [pos1]
            
        else:
            qs = qs
    
        if isinstance(pos1, Point):
            P = pos1
            if qs==None:
                diffs=P.vel(frame).magnitude().atoms(Derivative)
                
                qs= [deriv.args[0] for deriv in diffs]
                print(qs)
        else:
            P = Point('P')
            P.set_vel(frame,pos1.diff(ivar)*frame.x)
            force=-force*frame.x

        forcelist=[(P,force)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

###################################################################################################################################