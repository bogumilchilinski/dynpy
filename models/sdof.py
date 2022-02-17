from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, pi, latex, dsolve,
                   solve, fraction, factorial,Subs, Number, oo, Abs, N)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from ..continuous import ContinuousSystem, PlaneStressProblem

import base64
import random
import IPython as IP
import numpy as np
import inspect

class ComposedSystem(HarmonicOscillator):
    """Base class for all systems

    """
    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'

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

    def calculations_steps(self,preview=True,system=None,code=False):

#         latex_store=AutoBreak.latex_backend
#         AutoBreak.latex_backend = latex

        print('zlo')
        print(inspect.getsource(self.__class__))

        doc_model=super().calculations_steps(preview=True,code=code)


#         AutoBreak.latex_backend = latex_store
        return doc_model


    def get_default_data(self):
        return None

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()


        
        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict=None

        return parameters_dict


class SpringMassSystem(ComposedSystem):

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
                 z=dynamicsymbols('z'),
                 **kwargs
                 ):

        self.m = m
        self.k = k

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k, z, qs=[z])
        composed_system = self.mass + self.spring

        super().__init__(composed_system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict


class BeamBridge(ComposedSystem):
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
                 z=dynamicsymbols('z'),
                 **kwargs):

        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)
        composed_system = self.mass + self.spring + self.gravity_force + self.force

        super().__init__(composed_system,**kwargs)


    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict


class BeamBridgeTMD(ComposedSystem):

    scheme_name = 'bridge_tmd.png'
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
                 z_TMD=dynamicsymbols('z_TMD'),
                 **kwargs):

        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0
        self.m_TMD = m_TMD
        self.k_TMD = k_TMD
        self.z_TMD = z_TMD
        self.z = z

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.gravity_TMD = GravitationalForce(self.m_TMD, self.g, z_TMD)
        self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)
        self.TMD = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.spring_TMD = Spring(k_TMD, z, z_TMD, qs=[z, z_TMD])
        composed_system = (self.mass + self.spring + self.gravity_force + self.force +
                  self.TMD + self.spring_TMD + self.gravity_TMD)

        super().__init__(composed_system,**kwargs)



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
            self.k_TMD: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0],
        }

        return default_data_dict

<<<<<<< HEAD
class DampedHarmonicOscillator(ComposedSystem):

=======


    
class DampedSpringMassSystem(ComposedSystem):
>>>>>>> f61d6eec9ffde7f92f0d5ab570ea34046b6959de
    scheme_name = '???'
    real_name = 'engine_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 c=Symbol('c', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 **kwargs
                 ):

        self.m = m
        self.k = k
        self.c = c
        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k, z, qs=[z])
        self.damper = Damper(c,pos1=z, qs=[z])
        system = self.mass + self.spring + self.damper

        super().__init__(system,**kwargs)

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
    scheme_name = 'undamped_pendulum.png'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l

        Lagrangian = S.Half * m * l**2 * diff(
            angle, ivar)**2 - m * g * l * (1 - cos(angle))

        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar,**kwargs)

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.l: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0,7*l0, 8*l0, 9*l0,10*l0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie

class FreePendulum(ComposedSystem):
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
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.g = g
        self.l = l

        self.pendulum = Pendulum(m, g, l, angle=angle)
        system = self.pendulum

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict
    
    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, 1*m0, S.Half * m0, S.Half**2 * m0, 3*S.Half * m0],
            self.l: [2 * l0, 1*l0, S.Half * l0, S.Half**2 * l0, 3*S.Half * l0],
        }
        return default_data_dict

class ExcitedPendulum(ComposedSystem):
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
            **kwargs
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

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.F: r'Force',
        }
        return self.sym_desc_dict


class DampedPendulum(ComposedSystem):
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
            **kwargs
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
        self.Damper = Damper(c, l * phi, qs=qs)
        system = self.Pendulum + self.Damper

        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, l0, c0 = symbols('m_0 l_0 c_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.l: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.c: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0]
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict


class ExcitedDampedPendulum(ComposedSystem):

    scheme_name = 'damped_excited_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    def __init__(
            self,
            m=Symbol('m', positive=True),
            g=Symbol('g', positive=True),
            l=Symbol('l', positive=True),
            c=Symbol('c', positive=True),
            F=Symbol('F', positive=True),
            Omega=Symbol('Omega', positive=True),
            angle=dynamicsymbols('varphi'),
            qs=None,
            ivar=Symbol('t'),
            **kwargs
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
        self.F = F
        self.Omega = Omega

        self.Pendulum = Pendulum(m, g, l, angle=phi)
        self.Damper = Damper(c, l * phi, qs=qs)
        self.Force = Force(-F * sin(Omega * ivar), pos1=phi, qs=[phi])
        system = self.Pendulum + self.Damper + self.Force

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict


class PendulumKinematicExct(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'pendulum_real.jpg'

    def __init__(self,
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 phi=dynamicsymbols('\\varphi'),
                 x_e=dynamicsymbols('x_e'),
                 **kwargs):

        self.l = l
        self.m = m
        self.g = g
        self.phi = phi
        self.x_e = x_e

        x = l * sin(phi) + x_e
        y = l * cos(phi)

        self.material_point_1 = MaterialPoint(m, x, qs=[phi])
        self.material_point_2 = MaterialPoint(m, y, qs=[phi])
        self.gravity = GravitationalForce(m, g, pos1=-y, qs=[phi])

        system = self.material_point_1 + self.material_point_2 + self.gravity

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.x_e: r'Kinematic lateral excitation',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, S.Half * m0, 4 * m0, m0, S.Half**2 * m0, 8*m0, 16*m0],
            self.l: [2 * l0, S.Half * l0, 4 * l0, S.Half**2 * l0, 3 * l0, 3 * S.Half * l0, 9 * l0, 3*S.Half**2 * l0],
        }
        return default_data_dict


class Winch(ComposedSystem):

    scheme_name = 'sdof_winch.PNG'
    real_name = 'winch_mechanism_real.PNG'

    def __init__(self,
                 r=Symbol('r', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 phi=dynamicsymbols('\\varphi'),
                 **kwargs):

        self.r = r
        self.l = l
        self.m = m
        self.g = g
        self.phi = phi

        x = r * cos(phi) + (l + r * phi) * sin(phi)
        y = -r * sin(phi) + (l + r * phi) * cos(phi)

        self.material_point_1 = MaterialPoint(m, x, qs=[phi])
        self.material_point_2 = MaterialPoint(m, y, qs=[phi])
        self.gravity = GravitationalForce(m, g, pos1=-y, qs=[phi])

        system = self.material_point_1 + self.material_point_2 + self.gravity

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.r: r'Winch radius',
            self.l: r'Winch length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.r: [l0, 2 * l0, 4 * l0, l0, 8*l0],
            self.m: [2 * m0, S.Half * m0, 4 * m0, m0, S.Half**2 * m0, 8*m0, 16*m0],
            self.l: [2 * l0, S.Half * l0, 4 * l0, S.Half**2 * l0, 3 * l0, 3 * S.Half * l0, 9 * l0, 3*S.Half**2 * l0],
        }
        return default_data_dict

class Engine(ComposedSystem):
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
        >>> Engine()

    """

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

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
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


class DampedEngine(ComposedSystem):
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
                 **kwargs):

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
        self.damper = Damper(2 * c_m, pos1=z, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2 + self.damper
        super().__init__(system,**kwargs)

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
                 **kwargs):

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
        super().__init__(system,**kwargs)

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


class NonlinearEngine(ComposedSystem):
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'
    """
    Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

        Arguments:
        =========
            M = Mass
                -Mass of system (engine block) on spring

            m_e = Mass
                -Mass of particle

            e = distance
                -motion radius of a particle

            d = distance
                -distance between mounting point and engine

            k_m = spring coefficient
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
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.beta = beta
        self.d=d
        self.e = e
        self.l_0 = l_0
        self.z = z
        self.phi = phi

        N = ReferenceFrame('N')
        O = Point('O')

        P1 = Point('P1')
        P1.set_pos(O, 0 * N.x + 0 * N.y)

        P2 = Point('P2')
        P2.set_pos(O, d  * N.x + (z ) * N.y)

        self.MaterialPoint_1 = MaterialPoint(M, z, qs=[z])
        self.MaterialPoint_2 = MaterialPoint(m_e, z + e * cos(phi), qs=[z])
        self.Spring = Spring(2 * k_m, pos1=P1, pos2=P2, l_0=l_0, qs=[z])

        system = self.Spring + self.MaterialPoint_1 + self.MaterialPoint_2
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
            self.l_0: r'',
            self.beta: r'',
        }
        
        return self.sym_desc_dict
    def get_default_data(self):

        m0, k0, e0 = symbols('m_0 k_0 e_0', positive=True)

        default_data_dict = {
            self.M: [100*m0,300*m0,500*m0,700*m0,900*m0,200 * m0, 400 * m0,600*m0,800*m0],
            self.m_e: [m0,3*m0,5*m0,7*m0,9*m0,2 * m0, 4 * m0,6*m0,8*m0],
            self.k_m: [k0,2*k0,4*k0,6*k0,8*k0, 3 * k0,5*k0,7*k0,9*k0],
            self.e: [2 * e0, S.Half * e0, 4 * e0, S.Half**2 * e0,3 * e0,3* S.Half * e0, 9 * e0, 3*S.Half**2 * e0],
        }
        return default_data_dict
class StraightNonlinearEngine(NonlinearEngine):
    """
    Model of an exemplary Engine with nonlinear suspension aligned horizontally.

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
    scheme_name='non_linear_engine.png'

    def get_default_data(self):

        m0, k0, e0, l, omega = symbols('m_0 k_0 e_0 l Omega', positive=True)

        default_data_dict = {
            self.M: [100*m0,300*m0,500*m0,700*m0,900*m0,200 * m0, 400 * m0,600*m0,800*m0],
            self.m_e: [m0,3*m0,5*m0,7*m0,9*m0,2 * m0, 4 * m0,6*m0,8*m0],
            self.k_m: [k0,2*k0,4*k0,6*k0,8*k0, 3 * k0,5*k0,7*k0,9*k0],
            self.e: [2 * e0, S.Half * e0, 4 * e0, S.Half**2 * e0,3 * e0,3* S.Half * e0, 9 * e0, 3*S.Half**2 * e0],
            self.l_0:[S.Half*l,l,2*l],
            self.d:[4*l,8*l],
            self.beta:[S.Half*pi],
            self.phi:[omega*self.ivar]
        }
        return default_data_dict



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
                 **kwargs):

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

        super().__init__(system,**kwargs)


class TrolleyWithNonlinearSpring(ComposedSystem):
    scheme_name = 'sdof_nonlin_trolley.PNG'
    real_name = 'trolleywithnonlinearspring_real.png'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t', positive=True),
                 F=Symbol('F_0', positive=True),
                 x=dynamicsymbols('x'),
                 Omega=Symbol('Omega', positive=True),
                 **kwargs):
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
        self.d = d
        self.l_0 = l_0
        self.F = F

        self.MaterialPoint = MaterialPoint(m, x, qs=[x])
        self.Spring = Spring(k, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])
        self.Force = Force(F * cos(Omega * ivar), pos1=x, qs=[x])

        system = self.MaterialPoint + self.Spring + self.Force
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.k: 'Spring Stiffness',
            self.l: r'length',
            self.l_0: r'length',
            self.F: r'Force',
        }
        return self.sym_desc_dict
    
    def get_default_data(self):

        m0,l0,k0 = symbols('m_0 l_0 k_0', positive=True)

        default_data_dict = {
            self.m :[S.Half * m0, 1 * m0, 2 * m0, S.Half**2 * m0, 3*S.Half * m0],
            self.d :[S.Half * l0, 1 * l0, 2 * l0, S.Half**2 * l0, 3*S.Half * l0],
            self.k: [S.Half * k0, 1 * k0, 2 * k0, S.Half**2 * k0, 3*S.Half * k0],
        }

        return default_data_dict


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


class NonLinearTrolley(ComposedSystem):

    scheme_name = 'nonlin_trolley.PNG'
    real_name = 'nonlin_trolley_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 k=Symbol('k', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x'),
                 **kwargs):

        self.m = m
        self.k = k
        self.d = d
        self.l_0 = l_0
        self.x = x

        self.trolley = MaterialPoint(m, x, qs=[x]) + Spring(
            k, pos1=(sqrt(x**2 + d**2) - l_0), qs=[x])

        
        super().__init__(self.trolley,**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.d: [1 * l0, 2 * l0, S.Half * l0, 3 * S.Half * l0, 1 * l0],
            self.k:
            [S.Half * k0, S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }


        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Trolley Mass',
            self.k: 'Spring Stiffness',
            self.d: r'length',
            self.l_0: r'length',
        }
        return self.sym_desc_dict

class TriplePendulum(ComposedSystem):
    scheme_name = 'SDOFTriplePendulum.PNG'
    real_name = 'TriplePendulum_real.jpg'

    def __init__(self,
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 l_1=Symbol('l_1', positive=True),
                 l_2=Symbol('l_2', positive=True),
                 l_3=Symbol('l_3', positive=True),
                 g=Symbol('g', positive=True),
                 phi=dynamicsymbols('varphi'),
                 qs=dynamicsymbols('varphi'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.l_1 = l_1
        self.l_2 = l_2
        self.l_3 = l_3
        self.phi = phi
        self.g = g

        self.pendulum1 = Pendulum(m1, g, l_1, angle=phi, qs=[phi])
        self.pendulum2 = Pendulum(m2, g, (l_1+l_2), angle=phi, qs=[phi])
        self.pendulum3 = Pendulum(m3, g, (l_1+l_2+l_3), angle=phi, qs=[phi])

        system = self.pendulum1 + self.pendulum2 + self.pendulum3
        super().__init__(system(qs),**kwargs)


    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_3: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
        }

        return default_data_dict
    
    

    
    

class LagrangeIOnMathFunction(ComposedSystem):

    scheme_name = 'mat_point_parabola.PNG'
    real_name = 'tautochrone_curve_small.gif'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 x=dynamicsymbols('x'),
                 y=dynamicsymbols('y'),
                 a=symbols('a',positive=True),
                 R=symbols('R',positive=True),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('x,y'),
                 **kwargs):

        self.m = m
        self.x = x
        self.y = y
        self.a = a
        self.R = R
        self.g = g

        system = HarmonicOscillator(S.Half*m*x.diff(ivar)**2+S.Half*m*y.diff(ivar)**2-m*g*y,qs=[x,y])

        super().__init__(system(qs),**kwargs)

    def get_default_data(self):


        m0 = symbols('m_0', positive=True)
        x  = self.x
        a, Omega = symbols('a, Omega', positive=True)

        default_data_dict = {
            self.m :[S.Half * m0, 1 * m0, 2 * m0, 2**2 * m0, S.Half**2 * m0,8*m0,S.Half**3],
            self.y:[ a*x**2, a*(1-cos(x)),a*sin(x)**2,a*sin(x)**4,a*x**4]

        }

        return default_data_dict
        

class CrankSystem(ComposedSystem):

    scheme_name = 'crank_mechanism.png'
    real_name = 'crank_slider_real.jpg'

    def __init__(self,
                 I=Symbol('I', positive=True),
                 r=Symbol('r', positive=True),
                 h=Symbol('h', positive=True),
                 a=Symbol('a', positive=True),
                 b=Symbol('b', positive=True),
                 phi=dynamicsymbols('\\varphi'),
                 beta=dynamicsymbols('beta'),
#                  alpha=dynamicsymbols('alpha'),
                 **kwargs):

        self.I = I
        self.h=h
        self.r=r
        self.phi = phi
        self.beta = beta
#         self.alpha = alp
        self.a = a
        self.b = b

        self.crank = MaterialPoint(I, phi, qs=[phi])
        composed_system = (self.crank)

        super().__init__(composed_system,**kwargs)

    @property
    def _dbeta(self):
        beta=atan(self.r/self.l*self.phi) #it's probably wrong - has to be checked
        return beta.diff(self.ivar)
    
    @property
    def _displacement_d(self):

        return -(-self.b*sqrt((self.h**2 - self.h**2*self.a**2/self.b**2 + 2*self.h*self.r*cos(self.phi) - 2*self.h*self.r*self.a**2*cos(self.phi)/self.b**2 + self.r**2 - self.r**2*self.a**2*cos(self.phi)**2/self.b**2)/(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2)) - self.r*self.a*sin(self.phi)/sqrt(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2))
    
    @property
    def _velocity_b2(self):
        return (self.phi.diff(self.ivar)*self.r)
    @property
    def _velocity_b2b3(self):
        return (sqrt(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2)).diff(self.ivar)
    @property
    def _velocity_b3(self):
        gamma=asin(self._velocity_b2b3/self._velocity_b2)
        return self._velocity_b2*cos(gamma)
    @property
    def _velocity_d(self):
        return self._displacement_d.diff(self.ivar)
    @property
    def _velocity_c(self):
        omega3=self._velocity_b3/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi))
        return omega3*self.a
    @property
    def _acceleration_b2(self):
        return self.phi.diff(self.ivar)**2*self.r
    @property
    def _acceleration_b3n(self):
        return (self._velocity_b3**2/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))
    @property
    def _acceleration_cn(self):
        return (self._velocity_b3**2*(self.a/(sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))**2))
    @property
    def _acceleration_d(self):
        return self._velocity_d.diff(self.ivar)

    
    @property
    def _omega_3(self):
        return (self._velocity_b3/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))
    @property
    def linkage_ang_velocity(self):
        return self._dbeta

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I, l, m0, k0 = symbols('E I0 l_beam m_0 k_0', positive=True)
        
        default_data_dict = {
            self.I: [20 * I, 30 * I,60 * I,50 * I,40 * I,],
            self.h:[1,2,3],
            self.r:[0.1,0.2,0.3],
            self.a:[10],
            self.b:[20],
            self.phi : [10,20,30],
        }

        return default_data_dict

class WinchSystem(ComposedSystem):


    scheme_name = 'Winch_System.png'
    real_name = 'winch_mechanism_real.PNG'

    def __init__(self,
                 I_s=Symbol('I_s', positive=True),
                 I_k=Symbol('I_k', positive=True),
                 I_1=Symbol('I_1', positive=True),
                 i_1=Symbol('i_1', positive=True),
                 I_2=Symbol('I_2', positive=True),
                 i_2=Symbol('i_2', positive=True),
                 I_3=Symbol('I_3', positive=True),
                 I_4=Symbol('I_4', positive=True),
                 I_b=Symbol('I_b', positive=True),
                 A=Symbol('A', positive=True),
                 B=Symbol('B', positive=True),
                 ivar=Symbol('t'),
                 phi=dynamicsymbols('\\varphi'),
                 dphi=dynamicsymbols('\\varphi',1),
                 D=Symbol('D', positive=True),
                 mu=Symbol('\\mu', positive=True),
                 M_s=Symbol('M_s', positive=True),
                 M_T=Symbol('M_T', positive=True),
                 G=Symbol('G', positive=True),
                 g=Symbol('g',positive=True),
                 alpha=Symbol('\\alpha'),
                 **kwargs):
        
        pos1=phi
        qs=[phi]
        phi_1=phi
        phi_2=phi_1*i_1
        phi_3=phi_2*i_2

        self.I_s = I_s
        self.I_k = I_k
        self.I_1 = I_1
        self.i_1 = i_1
        self.I_2 = I_2
        self.I_3 = I_3
        self.i_2 = i_2
        self.I_4 = I_4
        self.I_b = I_b
        self.D = D
        self.G = G
        self.g = g
        self.mu = mu
        self.M_T = M_T
        self.M_s = M_s
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.alpha = alpha
        self.A=A
        self.B=B
        self.dphi = dphi

        
        self.engine = Disk(I_s, phi_1, qs=qs)
#         self.disk_k = Disk(I_k, phi_1, qs=qs)
        self.disk_1 = Disk(I_1, phi_1, qs=qs)
        self.disk_2 = Disk(I_2, phi_2, qs=qs)
        self.disk_3 = Disk(I_3, phi_2, qs=qs)
        self.disk_4 = Disk(I_4, phi_3, qs=qs)
        self.disk_B = Disk(I_b, phi_3, qs=qs)
        self.load = MaterialPoint(G/g, pos1=phi_3* D/2, qs=qs)
        self.friction_comp = GravitationalForce(G/g, g, pos1=phi_3* D/2 * cos(alpha)*mu , qs=qs)
        self.force = Force(-M_T,pos1=phi_3,qs=qs)
        self.drive = Force(M_s,pos1=phi_1,qs=qs)
        self.gravity_comp = GravitationalForce(G/g, g, pos1=phi_3* D/2 * sin(alpha) , qs=qs)

        system = self.engine + self.disk_1 + self.disk_2 + self.disk_3 + self.disk_4 + self.disk_B + self.load + self.friction_comp + self.gravity_comp + self.force + self.drive

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
            self.I_s: [20,30,40],
            self.I_k: [20,30,40],
            self.I_1: [20,30,40],
            self.I_2: [20,30,40],
            self.I_3: [20,30,40],
            self.I_4: [20,30,40],
            self.I_b: [20,30,40],
            self.i_1: [20,30,40],
            self.i_2: [20,30,40],
            self.D: [10,20,30],
            self.G: [10,20,30],
            self.g: [10,20,30],
            self.mu : [10,20,30],
            self.M_T : [10,20,30],
            self.M_s : [10,20,30],
            self.phi_1: [10,20,30],
            self.phi_2: [10,20,30],
            self.phi_3  : [10,20,30],
            self.alpha  : [10,20,30],
        }

        return default_data_dict
        
    def steady_angular_velocity(self):
        obj=self
        eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)
        return omega_steady_state
    
    def reduced_torque(self):
        obj=self
        
        red_tor = -(obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())

        return delta_1

    def reduced_inertia(self):
        obj = self
        
        ans = obj.inertia_matrix()[0]

        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.02)-1)*I_r).n(5))
        return ik
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_3
        return st_val
        
class DrivetrainVehicleSystem(ComposedSystem):


    #scheme_name = 'vehicle_drivetrain.png'
    #real_name = 'vehicle_drivetrain_real.PNG'

    def __init__(self,
                 R=Symbol('R'),
                 M_m=Symbol('M_m', positive=True),
                 omega_m=Symbol('\\omega_m', positive=True),
                 I_m=Symbol('I_m', positive=True),
                 I=Symbol('I', positive=True),
                 I_k=Symbol('I_k', positive=True),
                 omega_1=Symbol('\\omega_1', positive=True),
                 I_a=Symbol('I_a', positive=True),
                 z_1=Symbol('z_1', positive=True),
                 I_1=Symbol('I_1', positive=True),
                 z_2=Symbol('z_2', positive=True),
                 I_2=Symbol('I_2', positive=True),
                 I_b=Symbol('I_b', positive=True),
                 z_3=Symbol('z_3', positive=True),
                 I_3=Symbol('I_3', positive=True),
                 z_4=Symbol('z_4', positive=True),
                 I_4=Symbol('I_4', positive=True),
                 omega_3=Symbol('omega_3', positive=True),
                 I_c=Symbol('I_c',positive=True),
                 I_d=Symbol('I_d', positive=True),
                 I_w=Symbol('I_w', positive=True),
                 r=Symbol('r', positive=True),
                 m=Symbol('m', positive=True),
                 i_1=Symbol('i_1',positive=True),
                 i_2=Symbol('i_2',positive=True),
                 ivar=Symbol('t'),
                 A=Symbol('A', positive=True),
                 B=Symbol('B', positive=True),
                 phi=dynamicsymbols('\\varphi'),
                 dphi=dynamicsymbols('\\varphi',1),
                 alpha=Symbol('alpha'),
                 c=Symbol('c'),
                 g=Symbol('g'),
                 **kwargs):
        
        pos1=phi
        qs=[phi]
        phi_1=phi
#         i_1=z_2/z_1
#         i_2=z_4/z_3
        phi_2=phi_1/i_1
        phi_3=phi_2/i_2

        self.M_s = M_m
        self.omega_m = omega_m
        self.I_1 = I_1
        self.i_1 = i_1
        self.I_2 = I_2
        self.I_3 = I_3
        self.i_2 = i_2
        self.I_4 = I_4
        self.I_a = I_a
        self.I_b = I_b
        self.I_c = I_c
        self.I_d = I_d
        self.I_k = I_k
        self.I_m = I_m
        self.I_w = I_w
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.alpha = alpha
        self.A = A
        self.B = B
        self.dphi = dphi
        self.r = r
        self.D = 2*self.r
        self.m = m
        self.c = c
        self.g = g
        
        self.engine_force = Force(M_m, pos1=phi_1,qs=qs)
        self.engine_inertia = Disk(I_m, phi_1, qs=qs)
#         self.flywheel = Disk(I, phi_1, qs=qs)
        self.clutch = Disk(I_c, phi_1, qs=qs)
        self.driveshaft_1 = Disk(I_a, phi_1, qs=qs)
        self.transmission_i = Disk(I_1, phi_1, qs=qs)
        self.transmission_o = Disk(I_2, phi_2, qs=qs)
        self.driveshaft_2 = Disk(I_b, phi_2, qs=qs)
        self.diff_i = Disk(I_3, phi_2, qs=qs)
        self.diff_o = Disk(I_4, phi_3, qs=qs)
        self.axleshaft = Disk(I_d, phi_3, qs=qs)
        self.wheels = Disk(4*I_w, phi_3, qs=qs)
        self.mass = MaterialPoint(m, pos1=phi_3*r , qs=qs)
        self.gravity = Force(-m*g*sin(alpha), pos1=phi_3* r,qs=qs)#GravitationalForce(m, g, pos1=phi_3* r * sin(alpha) , qs=qs)
        self.air_force = Damper(c,pos1=phi_3*r,qs=qs) #(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame)

        system = self.engine_inertia + self.engine_force + self.clutch + self.driveshaft_1 + self.transmission_i + self.transmission_o + self.driveshaft_2 + self.diff_i + self.diff_o + self.axleshaft + self.wheels + self.mass + self.gravity + self.air_force

        super().__init__(system,**kwargs)
        
        self.__cached_solution = None

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I_c: r'moment of inertia',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
#             self.I_s: [20,30,40],
#             self.I_k: [20,30,40],
#             self.I_1: [20,30,40],
#             self.I_2: [20,30,40],
#             self.I_3: [20,30,40],
#             self.I_4: [20,30,40],
#             self.I_b: [20,30,40],
#             self.i_1: [20,30,40],
#             self.i_2: [20,30,40],
#             self.D: [10,20,30],
#             self.G: [10,20,30],
#             self.g: [10,20,30],
#             self.mu : [10,20,30],
#             self.M_T : [10,20,30],
#             self.M_s : [10,20,30],
#             self.phi_1: [10,20,30],
#             self.phi_2: [10,20,30],
#             self.phi_3  : [10,20,30],
        }

        return default_data_dict
    
    def linear_vehicle_velocity(self):
        obj=self
        lin_veh_vel=omega_steady_state*vel_ratio
        return lin_veh_vel
    
    def velocity_ratio(self):
        obj=self
        vel_ratio=self.phi_3*self.r/self.phi
        return vel_ratio
    

    def steady_angular_velocity(self):
        obj=self

        
        eom_sol = obj._eom_solution()
        #omega_steady_state =  (eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo))

        #display(obj._eoms[0])
        #eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        #display(eoms_num)
        #eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)

        return omega_steady_state
    
    def steady_angular_velocity_rev_per_min(self):
        obj=self
        ang_vel = obj.steady_angular_velocity()
        return ang_vel* 60/360
    
    def reduced_torque(self):
        obj=self
        red_tor = (obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())
#         0.4*(num_data.loc[case_no,'A']-Abs(M_Z))/(omega_ust**2*I_r)
        return delta_1

    def reduced_inertia(self):
        obj = self
        ans = obj.inertia_matrix()[0]
        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.004)-1)*I_r).n(5))
        return ik
    
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_3
        return st_val
    
    def reducted_mass(self):
        obj=self
        
        return obj.flywheel_moment_of_inertia()/(obj.startup_mass_velocity()*obj.phi1)**2
    
    def reducted_force(self):
        obj=self
        
        return obj.reduced_torque()/obj.startup_mass_velocity()


    
    def _eom_solution(self):
        obj=self
        
        if obj.__cached_solution is None:
        
            eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
            EOM_dsolve = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
            obj.__cached_solution = EOM_dsolve
        else:
            EOM_dsolve = obj.__cached_solution
            
       
        
        return EOM_dsolve