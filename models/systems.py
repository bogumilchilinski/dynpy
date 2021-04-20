from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, pi, latex, dsolve,
                   solve, fraction, factorial)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator
from .elements import MaterialPoint, Spring, NonlinSpring__RefFrme_Pt, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force

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

    def get_default_data(self):
        return None


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
        >>> mass = SDoFHarmonicOscillator(m,k, qs=[z],) # Initialization of LagrangesDynamicSystem instance

        -We define the symbols and dynamicsymbols
        -Kinetic energy T and potential energy v are evaluated to calculate the lagrangian L
        -Reference frame was created with point P defining the position and the velocity determined on the z axis
        -external forces assigned
        -Next we determine the instance of the system using class LagrangeDynamicSystem
        -We call out the instance of the class
        -If necessary assign values for the default arguments


    """
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z')):

        self.m = m
        self.k = k

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k, z, qs=[z])
        system = self.mass + self.spring

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict


class SDoFBeamBridge(ComposedSystem):
    """Ready to use model of bridge represented by the mass supported by elastic beam.
        Arguments:
        =========
            m = Symbol object
                -Mass embedded on beam.

            k = Symbol object
                -Bending stiffness of the beam

            g = Symbol object
                -Gravitational field acceleration

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass hanged on the elastic beam with the stiffness k in the gravitational field

        >>> t = symbols('t')
        >>> m, k = symbols('m, k')
        >>> qs = dynamicsymbols('z') # Generalized Coordinates
        >>> mass = SDoFHarmonicOscillator(m,k, qs=[z],) # Initialization of LagrangesDynamicSystem instance

        -We define the symbols and dynamicsymbols
        -Kinetic energy T and potential energy v are evaluated to calculate the lagrangian L
        -Reference frame was created with point P defining the position and the velocity determined on the z axis
        -external forces assigned
        -Next we determine the instance of the system using class LagrangeDynamicSystem
        -We call out the instance of the class
        -If necessary assign values for the default arguments


    """
    scheme_name = 'beam_bridge.PNG'
    real_name = 'beam_bridge_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 z=dynamicsymbols('z')):

        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)
        system = self.mass + self.spring + self.gravity_force + self.force

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict


class BeamBridgeTMD(ComposedSystem):

    scheme_name = 'beam_bridge.PNG'
    real_name = 'beam_bridge_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 m_TMD=Symbol('m_TMD', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 k_TMD=Symbol('k_TMD', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 z=dynamicsymbols('z'),
                 z_TMD=dynamicsymbols('z_TMD')):

        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0
        self.m_TMD = m_TMD
        self.k_TMD = k_TMD

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.gravity_TMD = GravitationalForce(self.m_TMD, self.g, z_TMD)
        self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)
        self.TMD = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.spring_TMD = Spring(k_TMD, z, z_TMD, qs=[z, z_TMD])
        system = (self.mass + self.spring + self.gravity_force + self.force +
                  self.TMD + self.spring_TMD + self.gravity_TMD)

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        E, I, l, m0, k0 = symbols('E I l_beam m_0 k_0', positive=True)

        default_data_dict = {
            self.m: [20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0],
            self.k_beam: [
                2 * 48 * E * I / l**3, 3 * 48 * E * I / l**3,
                4 * 48 * E * I / l**3, 5 * 48 * E * I / l**3,
                6 * 48 * E * I / l**3
            ],
            self.m_TMD: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.k_TMD: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0]
        }

        return default_data_dict


class SDoFDampedHarmonicOscillator(ComposedSystem):

    scheme_name = '???'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 c=Symbol('c', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z')):

        self.m = m
        self.k = k
        self.c = c
        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k, z, qs=[z])
        self.damper = Damper(c, z)
        system = self.mass + self.spring + self.damper

        super().__init__(system)


class DDoFSimplifyVehicleSuspension(ComposedSystem):
    """
    Ready to use sample Double Degree of Freedom System represents symmetrical kinematically excited beam with two springs.
        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            I = Moment of inertia
                -Moment of inertia of a rod

            l_rod = lenght of a rod
                -Dimension of a rod imitating vehicle's floor plate

            l_r = offset of right spring
                -Dimension of a right spring's offset

            l_l = offset of left spring
                -Dimension of a left spring's offset

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            F_engine = thrust of an engine
                -Force made by engine, which allows vehicle to move forward

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, k_1, k_2, l_rod, l_l, l_r, F_engine = symbols('m, k_1, k_2, l_rod, l_l, l_r, F_{engine}')
        >>> qs = dynamicsymbols('z, varphi') 
        >>> DDoFSymplifyVehicleSuspension()

        -We define the symbols and dynamicsymbols
        -and then we determine the instance of the system using class DDoFSymplifyVehicleSuspension()
    """

    scheme_name = 'car.png'
    real_name = 'car_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_rod', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t', positive=True),
                 qs=dynamicsymbols('z, phi')):

        z, phi = qs

        self.m = m  # mass of a rod
        self.l_l = l_l  # offset of left spring
        self.l_r = l_r  # offset of right spring
        self.l_rod = l_rod  # length of a rod
        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.F_engine = F_engine

        self.body = RigidBody2D(m, (m / 12) * (2 * l_rod)**2,
                                pos_lin=z,
                                pos_rot=phi,
                                qs=qs)
        self.spring_1 = Spring(k_1, pos1=z + phi * l_l, qs=qs)  # left spring
        self.spring_2 = Spring(k_2, pos1=z - phi * l_r, qs=qs)  # right spring
        self.force = Force(F_engine, pos1=z - l_r * phi, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.I: r'Moment of Inertia',
            self.l_rod: r'Length of the rod',
            self.l_l: r'offset of left spring',
            self.l_r: r'offset of right spring',
            self.k_1: r'Right spring stiffness coefficient',
            self.k_2: r'Left spring stiffness coefficient',
            self.F_engine: r'Force',
        }
        return self.sym_desc_dict


class DDoFVehicleSuspension(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents kinematically excited beam with two springs.
        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            I = Moment of inertia
                -Moment of inertia of a rod

            l_rod = lenght of a rod
                -Dimension of a rod imitating vehicle's floor plate

            l_r = offset of right spring
                -Dimension of a right spring's offset

            l_l = offset of left spring
                -Dimension of a left spring's offset

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            F_engine = thrust of an engine
                -Force made by engine, which allows vehicle to move forward

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, k_1, k_2, l_rod, l_l, l_r, F_engine = symbols('m, k_1, k_2, l_rod, l_l, l_r, F_{engine}')
        >>> qs = dynamicsymbols('z, varphi') # Generalized Coordinates
        >>> DDoFVehicleSuspension()

        -We define the symbols and dynamicsymbols
        -and then we determine the instance of the system using class DDoFVehicleSuspension()
    """

    scheme_name = 'car.png'
    real_name = 'car_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_rod', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t', positive=True),
                 qs=dynamicsymbols('z, \\varphi')):

        z, phi = qs

        self.m = m  # mass of a rod
        self.l_l = l_l  # offset of left spring
        self.l_r = l_r  # offset of right spring
        self.l_rod = l_rod  # length of a rod
        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.F_engine = F_engine

        self.body = RigidBody2D(m, I, pos_lin=z, pos_rot=phi, qs=qs)  # rod
        self.spring_1 = Spring(k_1, pos1=z + phi * l_l, qs=qs)  # left spring
        self.spring_2 = Spring(k_2, pos1=z - phi * l_r, qs=qs)  # right spring
        self.force = Force(F_engine, pos1=z - l_r * phi, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.I: r'Moment of Inertia',
            self.l_rod: r'Length of the rod',
            self.l_l: r'offset of left spring',
            self.l_r: r'offset of right spring',
            self.k_1: r'Right spring stiffness coefficient',
            self.k_2: r'Left spring stiffness coefficient',
            self.F_engine: r'Force',
        }
        return self.sym_desc_dict


class DDoFDampedVehicleSuspension(ComposedSystem):

    scheme_name = 'damped_car.png'
    real_name = 'car_real.jpg'

    def __init__(self,
                 non_damped_system,
                 c_l=Symbol('c_l', positive=True),
                 c_r=Symbol('c_r', positive=True),
                 l_cl=Symbol('l_{cl}', positive=True),
                 l_cr=Symbol('l_{cr}', positive=True),
                 k_1=DDoFVehicleSuspension().k_1,
                 k_2=DDoFVehicleSuspension().k_2,
                 l_l=DDoFVehicleSuspension().l_l,
                 l_r=DDoFVehicleSuspension().l_r,
                 qs=dynamicsymbols('z, \\varphi')):
        z, phi = qs
        self.k_1 = k_1
        self.k_2 = k_2
        self.c_l = c_l
        self.c_r = c_r
        self.l_cl = l_cl
        self.l_cr = l_cr
        self.l_l = l_l
        self.l_r = l_r
        self.nds = non_damped_system
        self.damper_l = Damper(c=c_l, pos1=z + phi * l_cl,
                               qs=qs)  # left damper
        self.damper_r = Damper(c=c_r, pos1=z - phi * l_cr,
                               qs=qs)  # right damper
        system = self.nds + self.damper_l + self.damper_r

        super().__init__(system)

    def get_default_data(self):

        c0, k_0, l_l0 = symbols('c_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.c_l: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.k_1: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            self.l_l: [2 * l_l0, 3 * l_l0, 4 * l_l0, 5 * l_l0, 6 * l_l0]
        }

        return default_data_dict


#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m: r'mass of system on the spring',
#             self.I: r'Moment of Inertia',
#             self.l_rod: r'Length of the rod',
#             self.l_l: r'offset of left spring',
#             self.l_r: r'offset of right spring',
#             self.k_1: r'Right spring stiffness coefficient',
#             self.k_2: r'Left spring stiffness coefficient',
#             self.F_engine: r'Force',
#         }
#         return self.sym_desc_dict


class DDoFShaft(ComposedSystem):
    """Ready to use sample Double Degree of Freedom System represents the Kinematicly excited shaft with two disks.
    =========
            I = Moment of Inertia
                -Moment of Inertia in case of both disc

            k_1 =Right spring coefficient
                -Right spring carrying the system

            k_2 =Left spring coefficient
                -Left spring carrying the system

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

    >>> t = symbols('t')
    >>> I, k1, k2 = symbols('I, k_1, k_2')
    >>> qs = dynamicsymbols('phi_1, phi_2') # Generalized Coordinates
    >>> DDoFShaft()

    -defines the symbols and dynamicsymbols
    -finally determines the instance of the system using class DDoFShaft
    """

    scheme_name = 'ddof_shaft.png'
    real_name = 'ddof_shaft_real.png'

    def __init__(self,
                 I=Symbol('I', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 input_displacement=dynamicsymbols('theta'),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('varphi_1, varphi_2')):

        phi1, phi2 = qs
        theta = input_displacement

        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.input_displacement = input_displacement
        self.qs = qs

        self.disc_1 = Disk(I, pos1=phi1, qs=qs)
        self.spring_1 = Spring(k_2, phi1, phi2, qs=qs)  # left spring
        self.disc_2 = Disk(I, pos1=phi2, qs=qs)
        self.spring_2 = Spring(k_1, pos1=phi2, pos2=theta,
                               qs=qs)  # right spring
        system = self.disc_1 + self.disc_2 + self.spring_1 + self.spring_2

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict


class Pendulum(ComposedSystem):
    """
    Model of a sDoF mathematical Pendulum. The "trig" arg follows up on defining the angle of rotation over a specific axis hence choosing apporperietly either sin or cos.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional's field acceleration

            l = lenght
                -Dimension of pendulum's strong

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """
    scheme_name = 'pendulum.png'
    real_name = 'pendulum_real.jpg'

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

        self.m = m
        self.g = g
        self.l = l

        Lagrangian = S.Half * m * l**2 * diff(
            angle, ivar)**2 - m * g * l * (1 - cos(angle))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict


# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie


class SDoFFreePendulum(ComposedSystem):
    """
    Model of a sDoF free pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> SDOFFreePendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SFoDFreePendulum()
    """
    scheme_name = 'free_sdof_pendulum.png'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t')):

        self.m = m
        self.g = g
        self.l = l

        self.pendulum = Pendulum(m, g, l, angle=angle)
        system = self.pendulum

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict


class SDoFExcitedPendulum(ComposedSystem):
    """
    Model of a sDoF Excited Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            F = Force
                -Pendulum's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l, F = symbols('m, g, l, F')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> SDoFExcitedPendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SDoFExcitedPendulum()
    """
    scheme_name = 'horizontal_forced_pendulum.png'
    real_name = 'pendulum2_real.jpg'

    def __init__(
            self,
            m=Symbol('m', positive=True),
            g=Symbol('g', positive=True),
            l=Symbol('l', positive=True),
            F=Symbol('F', positive=True),
            angle=dynamicsymbols('varphi'),
            qs=None,
            ivar=Symbol('t'),
    ):
        phi = angle

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.F = F

        self.pendulum = Pendulum(m, g, l, angle=phi)
        self.force = Force(-F * l * cos(phi), pos1=phi, qs=qs)
        system = self.pendulum + self.force

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.F: r'Force',
        }
        return self.sym_desc_dict


class SDoFDampedPendulum(ComposedSystem):
    """
    Model of a sDoF damped Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            c = damper coefficient
                -value of damper coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l, c = symbols('m, g, l, c')
        >>> qs = dynamicsymbols('varphi') # Generalized Coordinates
        >>> SDoFDampedPendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFDampedPendulum()
    """
    scheme_name = 'damped_pendulum.png'
    real_name = 'pendulum2_real.jpg'

    def __init__(
            self,
            m=Symbol('m', positive=True),
            g=Symbol('g', positive=True),
            l=Symbol('l', positive=True),
            c=Symbol('c', positive=True),
            angle=dynamicsymbols('varphi'),
            qs=None,
            ivar=Symbol('t'),
    ):
        phi = angle

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.c = c

        self.Pendulum = Pendulum(m, g, l, angle=phi)
        self.Force = Force(-c * diff(phi, ivar), pos1=phi, qs=qs)
        system = self.Pendulum + self.Force

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict


class DDoFCoupledPendulum(ComposedSystem):
    """
    Model of a DDoF Coupled Pendulum.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -Dimension of pendulum strong

            k = spring coefficient
                -value of spring coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, g, l, k = symbols('m, g, l, k')
        >>> qs = dynamicsymbols('varphi_1, varphi_2') # Generalized Coordinates
        >>> DDoFCouplePendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """
    scheme_name = 'mdof_dpendulum.png'
    real_name = 'lifting_tandem.png'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 k=Symbol('k', positive=True),
                 qs=dynamicsymbols('phi_1, phi_2')):

        phi1, phi2 = qs

        self.m = m
        self.g = g
        self.l = l
        self.k = k

        self.spring = Spring(k, pos1=(phi1 * (l)), pos2=(phi2 * (l)), qs=[qs])
        self.pendulum_1 = Pendulum(m, g, l, angle=phi1, qs=[qs])
        self.pendulum_2 = Pendulum(m, g, l, angle=phi2, qs=[qs])

        system = self.pendulum_1 + self.pendulum_2 + self.spring
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.k: r'Stifness coefficient',
        }
        return self.sym_desc_dict


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
    """
    Model of a SDoF engine.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            me = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            km = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A SDoF engine oscillating up and down while being held up by a two springs with the given stiffness

        >>> t = symbols('t')
        >>> M, me, e, km = symbols('M, m_e, e, k_m')
        >>  phi = dynamicsymbosl('phi')
        >>> qs = dynamicsymbols('z') #Generalized coordinate
        >>> SDoFEngine()

    """
    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.Spring = Spring(2 * k_m, pos1=z, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class SDoFDampedEngine(ComposedSystem):
    scheme_name = 'engine_with_damper.png'
    real_name = 'engine_real.PNG'
    """
    Model of a SDoF engine.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            me = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            k_m = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A SDoF engine oscillating up and down while being held up by a two springs with the given stiffness

        >>> t = symbols('t')
        >>> M, me, e, km = symbols('M, m_e, e, k_m')
        >>  phi = dynamicsymbosl('phi')
        >>> qs = dynamicsymbols('z') #Generalized coordinate
        >>> SDoFEngine()

    """
    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 c_m=Symbol('c_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.M = M
        self.k_m = k_m
        self.c_m = c_m
        self.m_e = m_e
        self.e = e
        self.phi = phi

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.Spring = Spring(2 * k_m, pos1=z, qs=[z])
        self.damper = Damper(2 * c_m, pos1=z)

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2 + self.damper
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'unbalanced rotating mass',
            self.e: r'radius of rotation',
        }
        return self.sym_desc_dict


class EngineWithTMD(ComposedSystem):
    """
    Model of a DDoF Coupled Pendulum.

        Arguments:
        =========
            M = Mass
                -Mass of system on spring

            m_e = mass
                -value of particle mass

            m_TMD = mass
                -value of TMD mass

            e = distance
                -motion radius of a particle

            k_m = spring coefficient
                -value of spring coefficient

            k_TMD = spring coefficient
                -value of spring coefficient

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> M, m_e, m_TMD, e, k_m, k_TMD = symbols('M, m_e, m_TMD, e, k_m, k_TMD')
        >>> qs = dynamicsymbols('z, z_TMD') # Generalized Coordinates
        >>> DDoFCouplePendulum()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = 'engine_TMD.PNG'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 k_TMD=Symbol('k_{TMD}', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 m_TMD=Symbol('m_{TMD}', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 z_TMD=dynamicsymbols('z_{TMD}'),
                 phi=dynamicsymbols('varphi'),
                 ivar=Symbol('t'),
                 system=None):

        self.M = M
        self.k_m = k_m
        self.k_TMD = k_TMD
        self.m_e = m_e
        self.m_TMD = m_TMD
        self.e = e
        self.z = z
        self.z_TMD = z_TMD
        self.phi = phi

        self.MaterialPoint_1 = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.MaterialPoint_3 = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.Spring_1 = Spring(2 * k_m, pos1=z, qs=[z])
        self.Spring_2 = Spring(k_TMD, pos1=z, pos2=z_TMD, qs=[z_TMD])

        system = self.Spring_1 + self.Spring_2 + self.MaterialPoint_1 + \
            self.MaterialPoint_2 + self.MaterialPoint_3
        super().__init__(system)

    def equilibrium_equation(self, static_disp_dict=None):
        static_disp_dict = {
            self.z: Symbol('z_0', positive=True),
            self.z_TMD: Symbol('z_{TMD0}', positive=True)
        }

        return super().equilibrium_equation(static_disp_dict=static_disp_dict)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class SDoFNonlinearEngine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'
    """
    Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            me = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            km = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            beta = angle
                -angle of springs that hold the system

            l0 = length
                -Initial length of non-linear springs

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        >>> t = symbols('t')
        >>> M, me, e, km, beta, l0 = symbols('M, m_e, e, k_m, beta, l_0')
        >>> qs = dynamicsymbols('z') 
        >>> SDoFNonlinearEngine()
    """
    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 beta=Symbol('beta', positive=True),
                 l0=Symbol('l_0', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.beta = beta
        self.e = e
        self.l0 = l0
        self.z = z
        self.phi = phi

        N = ReferenceFrame('N')
        O = Point('O')

        P1 = Point('P1')
        P1.set_pos(O, 0 * N.x + 0 * N.y)

        P2 = Point('P2')
        P2.set_pos(O, l0 * sin(beta) * N.x + (z + l0 * cos(beta)) * N.y)

        self.MaterialPoint_1 = MaterialPoint(M, z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e, z + e * cos(phi), qs=[z])
        self.Spring = Spring(2 * k_m, pos1=P1, pos2=P2, l0=l0, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
            self.beta: r'',
        }
        return self.sym_desc_dict


class MDoFTMD(ComposedSystem):
    scheme_name = 'mdof_tmd.png'
    real_name = 'mdof_tmd_real.png'
    """
    Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            me = Mass
                -Mass of tuned mass damper

            k = spring coefficient
                -value of spring coefficient

            ke = spring coefficient
                -value of spring coefficient that tuned mass damper is mounted

            F = Force
                -Trolley's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        >>> t = symbols('t')
        >>> m, me, k, ke, F = symbols('m, m_e, k, k_e, F')

        
        
    """
    def __init__(self,
                 m=Symbol('m', positive=True),
                 me=Symbol('m_e', positive=True),
                 k=Symbol('k', positive=True),
                 ke=Symbol('k_e', positive=True),
                 F=Symbol('F', positive=True),
                 xe=dynamicsymbols('x_e'),
                 xb=dynamicsymbols('x_b'),
                 angle=dynamicsymbols('Omega'),
                 ivar=Symbol('t', positive=True),
                 system=None):

        self.m = m
        self.me = me
        self.k = k
        self.ke = ke
        self.F = F

        self.MaterialPoint_1 = MaterialPoint(m, pos1=xe, qs=[xe])
        self.MaterialPoint_2 = MaterialPoint(me, pos1=xb, qs=[xb])
        self.Spring_1 = Spring(k, pos1=xe, qs=[xe])
        self.Spring_2 = Spring(ke, pos1=xe, pos2=xb, qs=[xe, xb])
        self.Force = Force(F * sin(angle * ivar), pos1=xe, qs=[xe])

        system = self.Spring_1 + self.Spring_2 + \
            self.MaterialPoint_1 + self.MaterialPoint_2 + self.Force
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of main object',
            self.me: r'Mass of the TMD',
            self.k: r'Stiffness coefficient',
            self.ke: r'Stiffness coefficient',
            self.F: r'Force',
        }
        return self.sym_desc_dict


class MDoFWinch(ComposedSystem):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = 'mdof_winch.png'
    real_name = 'mdof_winch_real.png'

    def __init__(self,
                 I=Symbol('I', positive=True),
                 k=Symbol('k', positive=True),
                 r=Symbol('r', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 theta=dynamicsymbols('theta'),
                 phi=dynamicsymbols('phi'),
                 system=None):

        self.I = I
        self.k = k
        self.r = r
        self.l = l
        self.m = m
        self.g = g
        self.theta = dynamicsymbols('theta')
        self.phi = dynamicsymbols('phi')

        x = r * cos(phi) + (l + r * (phi - theta)) * sin(phi)
        y = r * sin(phi) + (l + r * (phi - theta)) * cos(phi)
        F = m * g * r

        self.disc_1 = Disk(I, pos_c=theta, qs=[phi, theta])
        self.spring = Spring(k, theta, qs=[phi, theta])
        self.material_point_1 = MaterialPoint(m, x, qs=[phi, theta])
        self.material_point_2 = MaterialPoint(m, y, qs=[phi, theta])
        self.M_engine = Force(F, theta, qs=[phi, theta])
        self.gravity = GravitationalForce(m, g, pos1=-y, qs=[phi])

        system = self.material_point_1 + self.material_point_2 + \
            self.disc_1 + self.spring + self.M_engine + self.gravity

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia of the winch',
            self.k: r'Spring stiffness',
            self.r: r'Winch radius',
            self.l: r'Winch length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict


class MDoFElasticPendulum(ComposedSystem):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = 'mdof_winch.png'
    real_name = 'mdof_winch_real.png'

    def __init__(self,
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('varphi'),
                 system=None):

        self.k = k

        self.l = l
        self.m = m
        self.g = g
        self.phi = phi
        self.z = z

        x = (l + z) * sin(phi)
        y = (l + z) * cos(phi)

        self.frame = ReferenceFrame('N')

        self.payload = Point('payload')
        #self.payload.set_vel(
        #     frame, (sqrt((diff(x, ivar)**2 + diff(y, ivar)**2).simplify()))*frame.x)
        self.payload.set_vel(
            self.frame,
            sqrt(diff(z, ivar)**2 + (diff(phi, ivar) * (l + z))**2) * self.frame.x)

        #print(payload, 'try', type(payload))

        self.spring = Spring(k, z, qs=[phi, z])
        self.material_point_1 = MaterialPoint(m,
                                              self.payload,
                                              qs=[phi, z],
                                              frame=self.frame)
        # self.material_point_2 = MaterialPoint(m, y, qs=[phi, z])
        #self.M_engine = Force(F, theta, qs=[phi, z])
        self.gravity = GravitationalForce(m, g, pos1=-y, qs=[phi, z])
        system = (self.spring + self.gravity + self.material_point_1
                  )  # + self.material_point_2

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r'Spring stiffness',
            self.l: r'Winch length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict


class MDoFDampedElasticPendulum(ComposedSystem):
    """
    Model of a Double Degree of Freedom Involute Pendulum (Winch)

        Arguments:
        =========
            m = Mass
                -Mass of the payload

            I = Moment of Inertia
                -disc moment of inertia

            g = gravitional field
                -value of gravitional field acceleration

            l = lenght
                -initial length of the cable

            r = lenght
                -radius of the cylinder

            k = torsional stiffness coefficient
                -value of torsional spring coefficient

            ivar = symbol object
                -Independant time variable

            phi = dynamicsymbol object
                -pendulation angle of the mass m

            theta = dynamicsymbol object
                -oscillation angle of the cylinder

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass m pendulating on cable l_0 which is wounded on the cylinder with the radius R.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFCouplePendulum()
    """

    scheme_name = 'mdof_winch.png'
    real_name = 'mdof_winch_real.png'

    def __init__(self,
                 undamped_system,
                 c=Symbol('c', positive=True),
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('varphi')):

        self.c = c
        self.k = k
        self.l = l
        self.m = m
        self.g = g
        self.phi = phi
        self.z = z
        self.undamped = undamped_system
        self.damper = Damper(c=c,
                             pos1=self.undamped.payload,
                             qs=[phi, z],
                             frame=self.undamped.frame)

        display(self.damper._eoms)
        
        system = self.undamped + self.damper

        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r'Spring stiffness',
            self.l: r'Winch length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict


class Inverted_Pendulum(HarmonicOscillator):
    def __init__(self,
                 M=symbols('M', positive=True),
                 m=symbols('m', positive=True),
                 I=symbols('I', positive=True),
                 g=symbols('g', positive=True),
                 b=symbols('b', positive=True),
                 l=symbols('l', positive=True),
                 F=symbols('F', positive=True),
                 var=dynamicsymbols('x, phi'),
                 ivar=Symbol('t'),
                 system=None):

        x, phi = var

        self.rod = (
            RigidBody2D(m,
                        I,
                        pos_lin=0,
                        pos_rot=0,
                        pos_lin_c=(x + l * sin(phi)),
                        pos_rot_c=phi,
                        qs=[x, phi]) +
            MaterialPoint(m, l * cos(phi), qs=[phi]) +
            GravitationalForce(m, g, pos1=0, pos_c=l * cos(phi), qs=[phi]))

        self.cart = MaterialPoint(M, x, qs=[x])

        self.force = Force(F, x)

        self.friction = Damper(b, x)

        system = self.rod + self.cart + self.friction + self.force

        super().__init__(system)


class SDoFTrolleyWithNonlinearSpring(ComposedSystem):
    scheme_name = 'troleywithnonlinspring.PNG'
    real_name = 'trolleywithnonlinearspring_real.png'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t', positive=True),
                 F=Symbol('F_0', positive=True),
                 x=dynamicsymbols('x'),
                 Omega=Symbol('Omega', positive=True),
                 system=None):
        """
        Model of Single Degree of Freedom Trolley with nonlinear spring (type of inverted pendulum)

        Arguments:
        =========
            m = Mass
                -Mass of system on spring

            l = length
                -actual length of the non-linear spring

            l_0 = lenght
                -Initial length of non-linear spring

            k = spring coefficient
                -value of spring coefficient

            F = Force
                -Trolley's exciting force

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates

        Example
        =======
        A mass oscillating up and down while being held up by a spring with a spring constant kinematicly

        >>> t = symbols('t')
        >>> m, l, l_0, k, F = symbols('m, l, l_0, k, F')
        >>> qs = dynamicsymbols('x') # Generalized Coordinate
        >>> Omega = symbols('Omega')
        >>> SDoFTrolleyWithNonlinearSpring()

        -We define the symbols and dynamicsymbols
        -determine the instance of the pendulum by using class SDoFTrolleyWithNonlinearSpring()
    """

        self.m = m
        self.k = k
        self.l = l
        self.l_0 = l_0
        self.F = F

        self.MaterialPoint = MaterialPoint(m, x, qs=[x])
        self.Spring = Spring(k, pos1=(sqrt(x**2 + l**2) - l_0), qs=[x])
        self.Force = Force(-F * cos(Omega * ivar), pos1=x, qs=[x])

        system = self.MaterialPoint + self.Spring + self.Force
        super().__init__(system)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.k: 'Spring Stiffness',
            self.l: r'length',
            self.l_0: r'length',
            self.F: r'Force',
        }
        return self.sym_desc_dict


# class MDoFShaft(ComposedSystem):
#     scheme_name = '...'
#     real_name = '...'

#     def __init__(self, system=None, ivar=Symbol('t')):

#         t = ivar

#         m, m_0, k, M, k_m, g, F_1, F_2, Omega, F, R, e, m_e, J, k_m, beta, k_m = symbols(
#             'm,m_0,k,M,k_v,g,F_1,F_2,Omega, F_0, R, e, m_e, J, k_m, beta, k_m',
#             positive=True)

#         T = 0
#         V = 0

#         L_Shaft = (T - V)

#         shaft_base = HarmonicOscillator(
#             L_shaft, qs=[xb, xe], forcelist=[], frame=N)

#         super().__init__(shaft_base)
