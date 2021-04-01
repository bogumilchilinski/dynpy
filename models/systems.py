from sympy import *
from sympy.physics.mechanics import *

from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator
from .elements import *

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
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

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
        >>> m, k = symbols('m, k')
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
                 l_rod=Symbol('2l', positive=True),
                 l_l=Symbol('l', positive=True),
                 l_r=Symbol('l', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('z, varphi')):

        z, phi = qs

        self.m = m  # mass of a rod
        self.l_l = l_l  # offset of left spring
        self.l_r = l_r  # offset of right spring
        self.l_rod = l_rod  # length of a rod
        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.F_engine = F_engine
        self.qs = qs

        self.body = RigidBody2D(m, I, pos_lin=z, pos_rot=phi,
                                qs=qs)  # rod ---->
        self.spring_1 = Spring(k_1, pos1=z + phi * l_l, qs=qs)  # left spring
        self.spring_2 = Spring(k_2, pos1=z - phi * l_r, qs=qs)  # right spring
        self.force = Force(F_engine, pos1=z - l_r * phi, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force

        super().__init__(system)


class DDoFShaft(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents the Kinematicly excited shaft with two disks.
        =========

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
        >>> m, k = symbols('m, k')
        >>> qs = dynamicsymbols('z') # Generalized Coordinates 
        >>> DDoFShaft()

        -defines the symbols and dynamicsymbols
        -finally determines the instance of the system using class DDoFShaft
    """

    scheme_name = 'mdof_shaft_vibration.png'
    real_name = 'tandem_lift.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('2l', positive=True),
                 l_l=Symbol('l', positive=True),
                 l_r=Symbol('l', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 input_disp=dynamicsymbols('theta', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('z, varphi')):

        z, phi = qs

        self.m = m  # mass of a rod
        self.l_l = l_l  # offset of left spring
        self.l_r = l_r  # offset of right spring
        self.l_rod = l_rod  # length of a rod
        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.input_display = theta
        self.qs = qs

        self.disc_1 = Disk(I, pos1=phi, qs=qs)
        self.spring_1 = Spring(k_2, phi1, phi2, qs=qs)  # left spring
        self.disc_2 = Disk(I, pos1=phi2, qs=qs)
        self.spring_2 = Spring(
            k_1, pos1=phi2, pos2=theta, qs=qs)  # right spring
        system = self.disc_1 + self.disc_2 + self.spring_1 + self.spring_2 + self.force

        super().__init__(system)


class Pendulum(ComposedSystem):
    """
    Model of a sDoF mathematical Pendulum:

            Creates a singular model, after inputing correct values of mass - m , gravitational field - g, length of a strong - l and general coordinate which estabilshes an analytical display of a mathematical model of a sDoF pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.
    """
    scheme_name = 'pendulum.png'
    real_name = 'pendulum.png'

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

        Lagrangian = S.Half * m * l**2 * diff(
            angle, ivar)**2 - m * g * l * (1 - cos(angle))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)
# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie
# jeśli jest układ sdof, to nie indeksujemy parametrów


class SDoFPendulum(ComposedSystem):
    scheme_name = 'horizontal_forced_pendulum.png'
    real_name = 'pendulum2_real.jpg'

    def __init__(
            self,
            m1=Symbol('m_1', positive=True),
            g=Symbol('g', positive=True),
            l1=Symbol('l_1', positive=True),
            F=Symbol('F', positive=True),
            angle=dynamicsymbols('varphi_1'),
            qs=None,
            ivar=Symbol('t'),
    ):
        phi = angle

        self.pendulum = Pendulum(m1, g, l1, angle=angle)
        self.force = Force(-F * l1 * cos(phi), pos1=phi, qs=[phi])
        system = self.pendulum + self.force

        super().__init__(system)
# usunąć indeksy parametrów, zarówno w klasie jak i na rysunkach
# ujednolicić zmienną (angle albo phi)


class SDoFDampedPendulum(ComposedSystem):
    scheme_name = 'damped_pendulum.png'
    real_name = 'pendulum2_real.jpg'

    def __init__(
            self,
            m1=Symbol('m_1', positive=True),
            g=Symbol('g', positive=True),
            l1=Symbol('l_1', positive=True),
            c=Symbol('c', positive=True),
            angle=dynamicsymbols('varphi_1'),
            qs=None,
            ivar=Symbol('t'),
    ):
        phi = angle

        self.pendulum = Pendulum(m1, g, l1, angle=angle)
        self.force = Force(-c * diff(angle, ivar), pos1=phi, qs=[phi])
        system = self.pendulum + self.force

        super().__init__(system)
# to samo co w poprzedniej klasie


class DDoFCouplePendulum(ComposedSystem):
    scheme_name = 'mdof_dpendulum.png'
    real_name = 'tandem_lift.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 k=Symbol('k', positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('varphi1, varphi2')):

        phi, phi2 = qs

        self.spring = Spring(k, pos1=phi * l, pos2=phi2 * l, qs=qs)
        self.pendulum_1 = Pendulum(m, g, l, angle=phi, qs=qs)
        self.pendulum_2 = Pendulum(m, g, l, angle=phi2, qs=qs)
        system = self.spring + self.pendulum_1 + self.pendulum_2

        super().__init__(system)
# konsekwentnie używać indeksów dla mdofów - poprawić phi na phi1 (zgodnie z obrazkiem)

# class SDoFEngine(ComposedSystem):
#     scheme_name = 'engine.png'
#     real_name = 'engine_real.PNG'

#     def __init__(self, system=None,ivar=Symbol('t')):

#         t=ivar

#         m, m_0, k, M, k_m, g, F_1, F_2, Omega, F, R, e, m_e, J, k_m, beta, k_m = symbols(
#             'm,m_0,k,M,k_v,g,F_1,F_2,Omega, F_0, R, e, m_e, J, k_m, beta, k_m',
#             positive=True)

#         x_1, x_2, z, phi, theta = dynamicsymbols('x_1,x_2,z, varphi,theta')
#         dx_1, dx_2, dz = dynamicsymbols('x_1,x_2,z', 1)

#         T = S.Half * M * dz**2 + S.Half * m_e * (z + e * cos(phi)).diff(t)**2

#         V = S.Half * 2 * k_m * z**2

#         L_engine = (T - V)

#         engine_base = HarmonicOscillator(L_engine, qs=[z])

#         super().__init__(engine_base)


class SDoFEngine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 dz=dynamicsymbols('dz'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.MaterialPoint_1 = MaterialPoint(M, pos_c=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e, pos1=z+e*cos(phi), qs=[z])
        self.Spring = Spring(2*k_m, pos1=z, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2

        super().__init__(system)


class DDoFEngine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 k_TMD=Symbol('k_TMD', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 m_TMD=Symbol('m_TMD', positive=True),
                 e=Symbol('e', positive=True),
                 dz=dynamicsymbols('dz'),
                 z=dynamicsymbols('z'),
                 z_TMD=dynamicsymbols('z_TMD'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.MaterialPoint_1 = MaterialPoint(M, pos_c=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e, pos1=z+e*cos(phi), qs=[z])
        self.MaterialPoint_3 = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.Spring_1 = Spring(2*k_m, pos1=z, qs=[z])
        self.Spring_2 = Spring(k_TMD, pos2=z_TMD, qs=[z_TMD])

        system = self.Spring_1 + self.Spring_2 + self.MaterialPoint_1 + \
            self.MaterialPoint_2 + self.MaterialPoint_3

        super().__init__(system)


class SDoFNonlinearEngine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 beta=Symbol('beta', positive=True),
                 l0=Symbol('l_0', positive=True),
                 dz=dynamicsymbols('dz'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        N = ReferenceFrame('N')
        O = Point('O')

        P1 = Point('P1')
        P1.set_pos(O, 0*N.x + 0*N.y)

        P2 = Point('P2')
        P2.set_pos(O, l0*sin(beta)*N.x + (z + l0*cos(beta))*N.y)

        self.MaterialPoint_1 = MaterialPoint(M, z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e, z+e*cos(phi), qs=[z])
        self.Spring = Spring(2*k_m, pos1=P1, pos2=P2, l0=l0, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2

        super().__init__(system)
        
        
class MDoFTMD(ComposedSystem):
    scheme_name = 'mdof_tmd.png'
    real_name = 'mdof_tmd_real.png' # trzeba dodać zdjecie rela tmd bo nie mamy

    
    def __init__(self,
                 m=Symbol('m', positive=True),
                 me=Symbol('me', positive=True),
                 k=Symbol('k', positive=True),
                 ke=Symbol('ke', positive=True),
                 l0=Symbol('l_0', positive=True),
                 F=Symbol('F', positive=True),
                 xe=dynamicsymbols('xe'),
                 xb=dynamicsymbols('xb'),
                 angle=dynamicsymbols('Omega'),
                 ivar=Symbol('t', positive=True),
                 system=None):
    
    
        self.MaterialPoint_1 = MaterialPoint(m, pos1 = xe, qs=[xe])
        self.MaterialPoint_2 = MaterialPoint(me, pos1 = xb, qs=[xb])
        self.Spring_1 = Spring(k, pos1=xe,  qs=[xe])
        self.Spring_2 = Spring(ke, pos1=xe, pos2=xb, qs=[xe, xb])
        self.force = Force(F * sin(angle * ivar), pos1=xe, qs=[xe])
        
        
        system = self.Spring_1 + self.Spring_2 + self.MaterialPoint_1 + self.MaterialPoint_2 + self.Force

        super().__init__(system)
        

class MDoFWinch(ComposedSystem):
    #scheme_name = 'mdof_winch.png'
    #real_name = 'mdof_winch_real.png' # tez nei mam yzdjecia real wincha

    
    def __init__(self,
                 M=Symbol('M', positive=True),
                 M_engine=Symbol('M_{engine}', positive=True),
                 I = Symbol('I', positive=True),
                 k = Symbol('k', positive=True),
                 r = Symbol('r', positive=True),
                 l = Symbol('l', positive=True),
                 m = Symbol('m', positive=True),
                 g = Symbol('g', positive=True),
                 theta=dynamicsymbols('theta'),
                 phi=dynamicsymbols('phi'),
                 system=None):
        
        x = r*cos(theta) + l*sin(phi)
        y = r*sin(theta) - l*cos(phi)
        F = m*g*r
        
        
        self.disc_1 = Disk(I, pos_c = theta, qs=[phi, theta])
        self.Spring = Spring(k, theta, qs = [phi, theta])
        self.MaterialPoint_1 = MaterialPoint(m, x, qs=[phi, theta])
        self.MaterialPoint_2 = MaterialPoint(m, y, qs=[phi, theta])
        self.M_engine = Force(F, theta, qs = [phi, theta])
        
        system = self.MaterialPoint_1 + self.MaterialPoint_2 + self.disc_1 + self.Spring  + self.M_engine

        super().__init__(system.linearized())
        
class SDoFTrolleyWithNonlinearSpring(ComposedSystem):
    #scheme_name = ''
    #real_nanme = ''

    def __init__(self, m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t', positive=True),
                 F=Symbol('F_0', positive=True),
                 x=dynamicsymbols('x'),
                 Omega=Symbol('Omega', positive=True),
                 system=None):

        non_linear_spring_trolley = (MaterialPoint(m, x) +
                                     Spring(k, pos1=(sqrt(x**2 + l**2) - l_0), qs=[x]) + 
                                     Force(-F * cos(Omega*ivar), pos1=x, qs=[x]))

        super().__init__(non_linear_spring_trolley)


class MDoFTMD(ComposedSystem):
    #scheme_name = '...'
    #real_name = '...'

    def __init__(self, system=None, ivar=Symbol('t')):

        t = ivar

        xb, xe, z = dynamicsymbols('xb,xe,z')
        m_0, m, me, k_0, k, ke, F = symbols(
            'm_0, m, m_e, k_0, k, k_e, F', positive=True)

        T = S.Half * m * xb.diff(t)**2 + S.Half * m/10 * xe.diff(t)**2
        V = S.Half * k * xb**2 + S.Half * ke * (xe - xb)**2

        L_TMD = (T - V)

        tmd_base = HarmonicOscillator(
            L_TMD, qs=[xb, xe], forcelist=[])

        super().__init__(tmd_base)


class MDoFShaft(ComposedSystem):
    scheme_name = '...'
    real_name = '...'

    def __init__(self, system=None, ivar=Symbol('t')):

        t = ivar

        m, m_0, k, M, k_m, g, F_1, F_2, Omega, F, R, e, m_e, J, k_m, beta, k_m = symbols(
            'm,m_0,k,M,k_v,g,F_1,F_2,Omega, F_0, R, e, m_e, J, k_m, beta, k_m',
            positive=True)

        T = 0
        V = 0

        L_Shaft = (T - V)

        shaft_base = HarmonicOscillator(
            L_shaft, qs=[xb, xe], forcelist=[], frame=N)

        super().__init__(shaft_base)
