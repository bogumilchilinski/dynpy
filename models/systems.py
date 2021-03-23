from sympy import *
from sympy.physics.mechanics import *

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator
from .elements import MaterialPoint, Spring, RigidBody2D, Force

from sympy.physics.mechanics import *
from sympy.physics.vector import *

import base64
import IPython as IP


class ComposedSystem(HarmonicOscillator):
    """Base class for all systems
    
    """
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'
    
    
    @classmethod
    def _scheme(cls):
        
        path = __file__.replace('systems.py', 'images/') + cls.scheme_name

        
        return path
    
    @classmethod
    def _real_example(cls):
        
        path = __file__.replace('systems.py', 'images/') + cls.real_name

        
        return path
    
    @classmethod
    def preview(cls,example=False):
        if example:
            path=cls._real_example()
             
        else:
            path=cls._scheme()
            
        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()


        return IP.display.Image(base64.b64decode(encoded_string))


class SDoFHarmonicOscillator(ComposedSystem):
    """Ready to use sample Single Degree of Freedom System with mass on spring
        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            k = Spring coefficient
                -Spring carrying the system

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates
                
        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant k
        
        >>> t = symbols('t')
        >>> m, k = symbols('m k')
        >>> qs = dynamicsymbols('z') # Generalized Coordinates 
        >>> T = S.Half*m*z.diff(t)**2 # Kinetic Energy 
        >>> V = S.Half*k*z**2 # Potential Energy 
        >>> L = T - V # Lagrangian Calculation
        >>> N = ReferenceFrame('N') # Defining of reference frame for coordinate system
        >>> P = Point('P') # Defining point in space
        >>> P.set_vel(N, z.diff(t)*N.y) # Set velocity of point P in reference system N on axis z
        >>> Forcelist = [(P,f*sin(omega*t)*N.y)] # external forces on the system 
        >>> mass = dyn.HarmonicOscillator(dyn.LagrangesDynamicSystem(L, qs=[z], frame=N)) # Initialization of LagrangesDynamicSystem instance
        
        -We define the symbols and dynamicsymbols
        -Kinetic energy T and potential energy v are evaluated to calculate the lagrangian L
        -Reference frame was created with point P defining the position and the velocity determined on the z axis
        -external forces assigned 
        -finally we determine the instance of the system using class LagrangeDynamicSystem
    """
    scheme_name = 'engine.png'
    real_name = '???.png'
        
    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('z')):

        self.m = m
        self.k = k
        self.qs = qs

        self.mass = MaterialPoint(m, pos1=qs)
        self.spring = Spring(k, pos1=qs)
        system = self.mass + self.spring
        
        super().__init__(system)


class DDoFVehicleSuspension(ComposedSystem):
    scheme_name = 'car.png'
    real_name = 'car_real.jpg'
    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_rod', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_l=Symbol('k_2', positive=True),
                 k_r=Symbol('k_1', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('z, varphi')):

        z, phi = qs

        self.m = m  #mass of a rod
        self.l_l = l_l  #offset of left spring
        self.l_r = l_r  #offset of right spring
        self.l_rod = l_rod  #length of a rod
        self.k_l = k_l  #left spring
        self.k_r = l_r  #right spring
        self.I = I  #moment of inertia of a rod
        self.qs = qs

        self.body = RigidBody2D(m, I, pos_lin=z, pos_rot=phi,
                                qs=qs)  #rod ---->
        self.spring_1 = Spring(k_l, pos1=z + phi * l_l, qs=qs)  #left spring
        self.spring_2 = Spring(k_r, pos1=z - phi * l_r, qs=qs)  # right spring
        system = self.body + self.spring_1 + self.spring_2
        
        super().__init__(system)


class Pendulum(ComposedSystem):
    """
    Model of a sDoF mathematical Pendulum:

            Creates a singular model, after inputing correct values of mass - m , gravitational field - g, length of a strong - l and general coordinate which estabilshes an analytical display of a mathematical model of a sDoF pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.
    """
    scheme_name = 'pendulum.png'
    real_name = '???'
    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t')):

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        Lagrangian = S.Half * m * l**2 * diff(angle, ivar)**2 - m * g * l * (1 - cos(angle))
        
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)


class SDoFPendulum(ComposedSystem):
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'
    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 F=Symbol('F', positive=True),
                 ivar=Symbol('t'),
                 qs=[dynamicsymbols('varphi')]):
        phi = qs[0]

        self.pendulum = Pendulum(m, g, l, qs=qs)
        self.force = Force(F, pos1=phi,qs=[phi])
        system = self.pendulum + self.force
        
        super().__init__(system)


class DDoFDoublePendulum(ComposedSystem):
    scheme_name = '???.png'
    real_name = '???.png'
    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 k=Symbol('k', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('varphi, varphi2')):

        phi, phi2 = qs

        self.spring = Spring(k, pos1=phi * l, pos2=phi2 * l, qs=qs)
        self.pendulum_1 = Pendulum(m, g, l, angle=phi, qs=qs)
        self.pendulum_2 = Pendulum(m, g, l, angle=phi2, qs=qs)
        system = self.spring + self.pendulum_1 + self.pendulum_2
        
        super().__init__(system)

