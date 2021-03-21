from sympy import *
from sympy.physics.mechanics import *

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator
from .elements import MaterialPoint, Spring

from sympy.physics.mechanics import *
from sympy.physics.vector import *

import base64
import IPython as IP


class SDOFHarmonicOscillator(HarmonicOscillator):
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
    def __init__(self, m=Symbol('m',positive=True), k=Symbol('k',positive=True), ivar=Symbol('t'), qs=dynamicsymbols('z') ):
        
        self.m = m
        self.k = k
        self.qs = qs
        
        self.mass = MaterialPoint(m,pos1=qs)
        self.spring = Spring(k,pos1=qs)
        system = self.mass + self.spring
        
        super().__init__(system)
        
    @classmethod
    def preview(cls, img):
        with open(f'{img}', "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))
        
class SDOFcar(HarmonicOscillator):
    
    def __init__(self,
                 m=Symbol('m',positive=True),
                 l_l=Symbol('l_l',positive=True),
                 l_r=Symbol('l_r',positive=True),
                 l_rod=Symbol('l_rod',positive=True),
                 k_l=Symbol('k_l',positive=True),
                 k_r=Symbol('k_r',positive=True),
                 I=Symbol('I',positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('z') ):
        
        self.m = m #mass of a rod
        self.l_l = l_l #length of left spring
        self.l_r = l_r #length of right spring
        self.l_rod = l_rod #length of a rod
        self.k_l = k_l #left spring 
        self.k_r = l_r #right spring
        self.I = I #moment of inertia of a rod
        self.qs = qs
        
#         self.RigidBody2d = RigidBody2d(m, I, pos_lin = qs, pos_rot = ... ) #rod ----> nie wiem jak zorbic jego pozycje...
        self.spring_1 = Spring_2(k_l,pos1=qs) #left spring
        self.spring_2 = Spring_2(k_r,pos2=qs) # right spring
        system = self.spring_1 + self.spring_2 + self.RigidBody2d
        
        super().__init__(system)