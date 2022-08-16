from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Derivative,Integral, Expr,Function,latex, Heaviside, atan, pi)
from numbers import Number
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

import numpy as np

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, GeometryScene
from ..utilities.templates import tikz

#from dgeometry import GeometryScene,np

from pylatex import TikZ,TikZNode
import matplotlib.pyplot as plt

import base64
import random
import IPython as IP

base_frame=ReferenceFrame('N')
base_origin=Point('O')






class GeometryOfPoint:

    def __init__(self, *args , frame=base_frame, base_origin=base_origin , ivar=Symbol('t')):
        
        #print(type(args),args)
        #name = str(args[0])

        if isinstance(args[0], Point):
            self.Point = args[0]

        elif isinstance(args[0], Expr) or isinstance(args[0], Number):
            P1 = Point('P1')
            P1.set_pos(base_origin, frame.x*args[0] )
            P1.set_vel(frame, frame.x*diff(args[0], ivar) )

            self.Point = P1
            #print(self.Point)

            #print('Expression based on the generalized coordinate')

        elif isinstance(args, tuple):

            if len(args) == 2:
                P1 = Point('P1')
                P1.set_pos(base_origin, frame.x*args[0] + frame.y*args[1])
                P1.set_vel(frame, frame.x*diff(args[0], ivar) + frame.y*diff(args[1], ivar))
                self.Point = P1
                #print(self.Point)

            elif len(args) == 3:
                P2 = Point('P2')
                P2.set_pos(base_origin, frame.x*args[0] + frame.y*args[1] + frame.z*args[2])
                P2.set_vel(frame, frame.x*diff(args[0], ivar) + frame.y*diff(args[1], ivar) + frame.z*diff(args[2], ivar))
                self.Point = P2
                #print(self.Point)

            else:
                print('Tuple out of range')
                P_na = Point('P_na')
                P_na.set_pos(base_origin, frame.x*0)
                P_na.set_vel(frame, frame.x*diff(0, ivar))
                self.Point = P_na

    def get_point(self):
        return self.Point




class Element(LagrangesDynamicSystem):
    """Base class for all elements
    """

    scheme_name = 'damper.PNG'
    real_name = 'damper.PNG'
    _default_folder_path = "./dynpy/models/images/"
    

    
    @classmethod
    def _scheme(cls):

        path = cls._default_folder_path + cls.scheme_name
        
        return path

    @classmethod
    def _real_example(cls):
        path = cls._default_folder_path + cls.real_name

        return path
    
    @classmethod
    def _detail_real(cls):
        path = cls._default_folder_path + cls.detail_real_name

        return path
    
    @classmethod
    def _detail_scheme(cls):
        path = cls._default_folder_path + cls.detail_scheme_name

        return path

    
    def get_default_data(self):
        return {}

    def get_numerical_data(self):
        return {}
    
    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    def get_numerical_parameters(self):

        default_data_dict = self.get_numerical_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict
    
    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))


    def __str__(self):
        

        return f'{(self._label)}'    

        
#Pati
class MaterialPoint(Element):
    """
    Model of a Material point with changing point of mass:

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'
    def __init__(self, m, pos1 , qs=None, frame=base_frame, ivar=Symbol('t')):
        
        if not qs:
            self.qs = [pos1]

        pos1=GeometryOfPoint(pos1).get_point()
            
        if isinstance(pos1, Point):
            
            Lagrangian = S.Half * m * ((base_origin.pos_from(pos1).diff(ivar,frame).magnitude())**2)
            #Lagrangian = S.Half * m * ((base_origin.vel(frame).magnitude())**2)
            
        else:
            Lagrangian = S.Half * m * diff(pos1,ivar)**2


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar,frame=frame)
        
        self._kinetic_energy = Lagrangian
        
        

#Mateusz
class Spring(Element):
    """
    Model of a Spring: Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE.
    """
    scheme_name = 'spring.PNG'
    real_name = 'spring_real.png'
    
    stiffness=Symbol('k',positive=True)
    l_s=Symbol('l_s', positive=True)
    k0 =Symbol('k_0',positive=True)
    
    def __init__(self, stiffness, pos1, pos2=0, l_s=None ,  qs=None, ivar=Symbol('t'), frame = base_frame):

        if stiffness is not None: self.stiffness=stiffness
        if pos1 is not None: self._pos1=pos1
        self._pos2=pos2
        if l_s is not None: self.l_s=l_s
        else:
            l_s=0
            self.l_s=l_s

        if not qs:
            qs = [pos1]

        Lagrangian = (-S.Half * self.stiffness * ((self.pos2.pos_from(self.pos1)).magnitude() - self.l_s  )**2 )


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame)

        self._potential_energy = - Lagrangian

    @property
    def pos1(self):
        return GeometryOfPoint(self._pos1).get_point()

    @property
    def pos2(self):
        return GeometryOfPoint(self._pos2).get_point()

    def spring_force(self):
        if self.pos2==0:
            force = (self.stiffness * self.pos1)
        else:
            force = self.stiffness * (self.pos2.pos_from(self.pos1)).magnitude() - self.l_s
        return force

    def _plot_2d(self, language='en'):

        class_name = self.__class__.__name__

        coords=[(0,0),(0,3),(3,6),(-3,12),(0,15),(0,18)]
        x_coords = [ x+5 for x,y in coords]
        y_coords = [ y-10 for x,y in coords]

        res = GeometryScene.ax_2d.plot(x_coords,y_coords, label=class_name)
        print(res)
        
        
    def get_default_data(self):

        k0= self.k0
        self.default_data_dict={
            self.stiffness: [S.One / 100 *  k0 * no for no in range(80, 135)],
        }
        return {**super().get_default_data(),**self.default_data_dict}
        
    def get_numerical_data(self):

        self.default_numerical_dict={
            self.stiffness: [S.One / 100 *  50e3 * no for no in range(80, 135)],
        }
        
        return {**super().get_numerical_data(),**self.default_numerical_dict}
    
    @classmethod 
    def from_default(cls,
                 stiffness=None,
                 pos1=None,
                 frame=base_frame,
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if stiffness is not None: cls.stiffness=stiffness
        if pos1 is None: pos1= dynamicsymbols('z') 
        if not qs: 
            qs = [pos1] 
        ivar=cls.ivar

        return cls(stiffness=stiffness, pos1=pos1, qs=qs,ivar=ivar, frame = frame)

class EngineMount(Spring):

    """
    Model of a Engine Mount: Creates a singular model, after inputing correct values of stiffeness - k and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE.
    """
    real_path='buick_regal_3800.jpg'

    def pin_diameter(self):
        #średnica sworznia ze wzoru na ścinanie
        kt=Symbol('k_t', positive=True)
        Re=Symbol('R_e', positive=True)

        return sqrt((4*S.One*self.spring_force())/(pi*kt*Re)).doit()

    def get_numerical_data(self):
        self.default_data_dict={

        }
        return default_data_dict
    

class GravitationalForce(Element):
    """
    Model of a changing centroid for potential energy. Creates a singular model, after inputing correct values of gravity field - g,
     mass of - m as well as additionaly the general coordiante
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, m, g, pos1, qs=None, ivar=Symbol('t')):

        if not qs:
            qs=[pos1]

        Lagrangian = -(m * g * (pos1))


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        self._potential_energy = - Lagrangian


class Disk(Element):
    """
    Model of a Disk:
    Creates a singular model, after inputing correct values of moment of inertia - I and rotational general coordinate, which analytically displays the dynamics of a rotating wheel.
    """
    scheme_name = 'disk.png'
    real_name = 'disk.png'
    def __init__(self, I, pos1 , qs=None, frame=base_frame, ivar=Symbol('t')):
        
        if not qs:
            self.qs = [pos1]

        pos1=GeometryOfPoint(pos1).get_point()

        if isinstance(pos1, Point):
            Lagrangian = S.Half * I * ((base_origin.pos_from(pos1).diff(ivar,frame).magnitude())**2)
                
        else:
            Lagrangian = S.Half * I * diff(pos1,ivar)**2


        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame)
        self._kinetic_energy = Lagrangian

        

class RigidBody2D(Element):
    """
    Model of a 2DoF Rigid body:
    """
    """
    Creates a singular model, after inputing correct values of moment of inertia - I for rotational element and mass - m for linear element and general coordinates, which             analytically display the dynamics of rotating and propeling object: beam.
    """
    scheme_name = 'rigid_body2D.png'
    real_name = 'rigid_body2D.png'
    def __init__(self, m, I, pos_lin=0, pos_rot=0, qs=None, frame=base_frame, ivar=Symbol('t')):
        
        if pos_lin == 0 and pos_rot == 0:
            self.qs = [pos_lin, pos_rot]

        if isinstance(pos_lin, Point) and isinstance(pos_rot, Point):
            Lagrangian = S.Half * m * ((base_origin.pos_from(pos_lin).diff(ivar,frame).magnitude())**2) + S.Half * I *((base_origin.pos_from(pos_rot).diff(ivar,frame).magnitude())**2)
        
        else:
            Lagrangian = S.Half * m * diff(pos_lin,ivar)**2 + S.Half * I * diff(pos_rot,ivar)**2

        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar, frame=frame)
        
        self._kinetic_energy = Lagrangian
        
class Damper(Element):
    """
    Model of a Damper:

    Creates a singular model, after inputing correct values of the damping coefficient - c and general coordinates, which establishes a damping force directly proportional to the velocity, in time, of an inertial virbating element.
    """
    scheme_name = 'damper.PNG'
    real_name = 'damper.PNG'
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame):
        
        if not qs:
            self.qs = [pos1]

        if isinstance(pos1, Point):
            pos1=GeometryOfPoint(pos1,frame=frame).get_point()
            pos2=GeometryOfPoint(pos2,frame=frame).get_point()
            
            
            
            D = S.Half * c * ((pos1.vel(frame)-pos2.vel(frame)).magnitude()**2).doit()  # pos1.vel(frame).magntidue() behaves as diff(pos1,ivar)
        else:
            D = S.Half * c * (diff(pos1-pos2,ivar))**2

        points_dict={}
        for coord in qs:
            
            coord_vel=diff(coord,ivar)
            
            P_tmp = Point(f'P_{str(coord)}')
            P_tmp.set_vel(frame, coord_vel * frame.x)
            points_dict[coord_vel]=P_tmp

    
        forcelist = [ (point_dmp,-diff(D,coord_vel)*frame.x)  for coord_vel,point_dmp in points_dict.items() ]
        #print(forcelist)
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        self._dissipative_potential = D
        


        
class Excitation(Element):
    """
    Model of a harmonic extorsion applied onto the elemnt:
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, f, pos_rot, ivar=Symbol('t'), frame=base_frame):
        
        qs = [pos_rot]
        
        pos_rot = GeometryOfPoint(pos_rot).get_point()
        
        dpos_rot = diff(pos_rot,ivar)
        
        P = Point('P')
        P.set_vel(frame, dpos_rot * frame.x)
        P.vel(frame)

        F = f * cos(dpos_rot*ivar) + f * sin(dpos_rot*ivar)
        
        forcelist = [(P, F*frame.x)]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

class Force(Element):
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

        if qs is None:
            qs = [pos1]
            
        else:
            qs = qs
    
        if isinstance(pos1, Point):
            P = pos1
            if qs is None:
                diffs=P.vel(frame).magnitude().atoms(Derivative)
                
                qs= [deriv.args[0] for deriv in diffs]
                print(qs)
        else:
            
            P = Point('P')
            P.set_vel(frame,pos1.diff(ivar)*frame.x)
            force=+force*frame.x

        forcelist=[(P,force)]
        
        #display(forcelist)
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)
        

class CombustionEngine(Force):

    eta= Symbol('eta', positive = True)
    n_engine = Symbol('n')
    omega=Symbol('omega')

    def __init__(self, omega, characteristic=None, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

        if characteristic is None:
            self._characteristic= self.characteristic
        else:
            self._characteristic =  characteristic
        if eta is not None: self.eta = eta
        if omega is not None: self.omega = omega
            

        super().__init__(self.eta*self._characteristic, omega.integrate(ivar), qs=qs, ivar=ivar)

    @property
    def characteristic(self):
        n_engine = self.omega*60/(2*np.pi)
        return -4e-5*self.n_engine**2 + 0.2157*self.n_engine + 113.07

    @classmethod
    def from_data(cls, data, n, degree=5, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):

        omega = 2*3.14/60*n
        
        x = data.index.to_numpy()
        y = data.iloc[:,0].to_numpy()

        coeffs=reversed((np.polyfit(x, y, degree)))

        
        omega_min = 2*3.14/60*x[0]
        omega_max = 2*3.14/60*x[-1]
        
        char_slope = 10
        
        window   =   ((S.One/2+atan(char_slope*(omega-omega_min))/np.pi))    -  ((S.One/2+atan(char_slope*(omega-omega_max))/np.pi))
        #window = Heaviside(omega-omega_min) - Heaviside(omega-omega_max)
        
        omega = 2*3.14/60*n
        
        char_poly  =window *sum([ coeff * n ** no for no, coeff in enumerate(coeffs)])
        

        
        obj = cls(omega = omega, characteristic = char_poly, eta=eta, qs = qs, ivar = ivar, frame = frame)
        obj._char_poly  =  char_poly
        
        return obj

    @classmethod
    def from_rpm(cls, n, eta = None, qs=None, ivar=Symbol('t'),frame = base_frame):
        omega = 2*3.14/60*n
        return cls(omega = omega, eta = eta, qs=qs, ivar=ivar, frame = frame)
    
class IntegralElement(Force):

    k_I = Symbol('k_I', positive=True)
    error = dynamicsymbols('e')
    target = dynamicsymbols('target')
    reference = S.Zero
    
    def __init__(self, k_I, error, target=None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):
        
        if reference is not None: self.reference = reference
        
        if target is None and isinstance(error,Function):
               target=error

        super().__init__(k_I*(self.reference - error).integrate(ivar), target, qs=qs, ivar=ivar)


class DerivativeElement(Force):
    
    k_D = Symbol('k_D', positive=True)
    error = dynamicsymbols('e')
    target = dynamicsymbols('target')
    reference = S.Zero
    
    def __init__(self, k_D, error, target=None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):
        
        if reference is not None: self.reference = reference
        
        if target is None and isinstance(error,Function):
               target=error

        super().__init__(k_D*diff((self.reference - error),ivar), target, qs=qs, ivar=ivar)


class ProportionalElement(Force):

    k_P = Symbol('k_P', positive=True)
    error = Symbol('e', positive=True)
    target = Symbol('target', positive=True)
    reference = S.Zero
    
    def __init__(self, k_P, error, target=None, reference = None,  qs=None, ivar=Symbol('t'), frame = base_frame):
        
        if reference is not None: self.reference = reference
        
        if target is None and isinstance(error,Function):
               target=error

        super().__init__(k_P*(self.reference - error), target, qs=qs, ivar=ivar)

        
        
class PID(Element):
    """
    Model of a PID controller:

Creates a model of a PID controller (proportional , integral , derivative) which initates a smoother control over delivered system oscillatory system in form of Lagrange method. Taking into account that a PID regulator is ment for a SISO system the model is built for sDoF system, hence in case of building a mDoF the process will require using two PIDs.
    """
    scheme_name = ''
    real_name = ''
    def __init__(self, kp, ki, kd, pos1, qs=None, ivar=Symbol('t'), frame=base_frame):

        if qs:
            qs = [pos1]

        dpos1 = diff(pos1, ivar)

        In = pos1  # In - error
        u = pos1*ki + dpos1*kp + d2pos1*kd
        G = In - u
        
        forcelist = [(P, G*frame.x)]
        
        pos1=GeometryOfPoint(pos1,frame=frame).get_point()
        
        if isinstance(pos1, Point):
            D = S.Half * c * ((pos1.vel(frame)-pos2.vel(frame)).magnitude()**2).doit()  #pos1*ki + pos1.vel(frame)*kp + d2pos1*kd
        else:
            D = S.Half * c * (diff(pos1-pos2,ivar))**2

        points_dict={}
        for coord in qs:
            
            coord_vel=diff(coord,ivar)
            
            P_tmp = Point(f'P_{str(coord)}')
            P_tmp.set_vel(frame, coord_vel * frame.x)
            points_dict[coord_vel]=P_tmp

        forcelist = [ (point_dmp,-diff(G,coord_vel)*frame.x)  for coord_vel,point_dmp in points_dict.items() ]
        
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)