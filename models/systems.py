from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, pi, latex, dsolve,
                   solve, fraction, factorial,Subs, Number)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from ..continuous import ContinuousSystem, PlaneStressProblem

import base64
import random
import IPython as IP
import numpy as np

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


class ContinuousSystem(ContinuousSystem):
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
    
class PlaneStressProblem(PlaneStressProblem):
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


class BeamBridgeDampedTMD(ComposedSystem):

    scheme_name = 'bridge_tmd_dmp.png'
    real_name = 'beam_bridge_real.PNG'

    def __init__(self,
                 non_damped_system,
                 c_TMD=Symbol('c_TMD', positive=True),
                 c=Symbol('c', positive=True),
                 z_TMD=BeamBridgeTMD().z_TMD,
                 z=BeamBridgeTMD().z,
                 **kwargs):
        qs = z, z_TMD
        self.nds = non_damped_system
        self.c = c
        self.c_TMD = c_TMD

        self.damper_TMD = Damper(c=c_TMD, pos1=z - z_TMD, qs=qs)  # left damper
        self.damper = Damper(c=c, pos1=z, qs=qs)
        system = (self.nds + self.damper + self.damper_TMD)

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self,real=False):

        E, I, l, m0, k0, lamb = symbols('E I l_beam m_0 k_0 lambda',
                                        positive=True)

        if real is False:
            
            default_data_dict = {
                self.nds.m: [20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0],
                self.c: [lamb * self.nds.k_beam],
                self.nds.k_beam: [
                    20 * 48 * E * I / l**3, 30 * 48 * E * I / l**3,
                    40 * 48 * E * I / l**3, 50 * 48 * E * I / l**3,
                    60 * 48 * E * I / l**3
                ],
                self.c_TMD: [lamb * self.nds.k_TMD],
                self.nds.m_TMD: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
                self.nds.k_TMD: [
                    1 * 48 * E * I / l**3, 3 * 48 * E * I / l**3,
                    3 * 48 * E * I / l**3, 5 * 48 * E * I / l**3,
                    5 * 48 * E * I / l**3
                ],
            }
        else:
            numerized_dict = {m0:100 , lamb:2 , E:200000, I:0.48, l:3}
            default_data_dict = {
                self.nds.m: [20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0],
                self.c: [lamb * self.nds.k_beam],
                self.nds.k_beam: [
                    20 * 48 * E * I / l**3, 30 * 48 * E * I / l**3,
                    40 * 48 * E * I / l**3, 50 * 48 * E * I / l**3,
                    60 * 48 * E * I / l**3],
                self.c_TMD: [lamb * self.nds.k_TMD],
                self.nds.m_TMD: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
                self.nds.k_TMD: [
                    1 * 48 * E * I / l**3, 3 * 48 * E * I / l**3,
                    3 * 48 * E * I / l**3, 5 * 48 * E * I / l**3,
                    5 * 48 * E * I / l**3]
            }

            values_of_dict = flatten([default_data_dict[self.nds.m], default_data_dict[self.c],
                                      default_data_dict[self.nds.k_beam], default_data_dict[self.c_TMD],
                                      default_data_dict[self.nds.m_TMD], default_data_dict[self.nds.k_TMD]
                                     ])

            for i in range(len(values_of_dict)):
                    default_data_dict = display(values_of_dict[i].subs(numerized_dict))

        return default_data_dict


    def get_real_data(self):
        
        return self.get_default_data(True)
    
class SDoFDampedHarmonicOscillator(ComposedSystem):
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

    scheme_name = 'vehicle_suspension.PNG'
    real_name = 'vehicle_suspension_real.PNG'

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
                 qs=dynamicsymbols('z, \\varphi'),
                 **kwargs):

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
        self.force = Force(F_engine, pos1=z - phi * l_r, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force

        super().__init__(system,**kwargs)

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
    
    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)
        e,m_e,nu,t = symbols('e m_e \\nu t',positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.k_1: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, 3*S.Half * k0],
            self.k_2: [S.Half * k0, S.Half**2 * k0, 1* k0, 2 * k0, 3 * k0],
            self.l_l: [1 * l0, 2 * l0, S.Half * l0, 3 * l0, 4 * l0],
            self.l_r: [3*S.Half * l0, 2 * l0, 1 * l0, 3 * l0, 5 * l0],
            self.l_rod: [S.Half * l0, 2 * l0, 1 * l0, 4 * l0, 6 * l0],
            self.F_engine: [e*m_e*nu**2*cos(nu*t)],
            m_e: [S.Half**2 * m0, S.Half**3 * m0, S.Half**4 * m0, S.Half**5 * m0, S.Half**5 * m0],
            e: [S.Half**2 * l0, S.Half**3 * l0, S.Half**4 * l0, S.Half**5 * l0, S.Half**5 * l0],
        }

        return default_data_dict


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

    scheme_name = 'car.PNG'
    real_name = 'car_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_{rod}', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_2=Symbol('k_r', positive=True),
                 k_1=Symbol('k_l', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t', positive=True),
                 qs=dynamicsymbols('z, \\varphi'),
                 **kwargs):

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

        super().__init__(system,**kwargs)

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

    scheme_name = 'damped_car_new.PNG'
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
                 qs=dynamicsymbols('z, \\varphi'),
                 **kwargs):

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

        super().__init__(system,**kwargs)

    def get_default_data(self):

        c0, k_0, l_l0, omega, F_0 = symbols('c_0 k_0 l_0 omega F_0',
                                            positive=True)

        default_data_dict = {
            self.c_r: [self.c_l],
            self.k_2: [self.k_1],
            self.l_r: [self.l_l],
            self.l_cr: [self.l_l],
            self.l_cl: [self.l_l],
            self.c_l: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.k_1: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            self.l_l: [2 * l_l0, 3 * l_l0, 4 * l_l0, 5 * l_l0, 6 * l_l0],
            self.nds.F_engine: [
                2 * F_0 * sin(omega * self.nds.ivar),
                3 * F_0 * sin(omega * self.nds.ivar),
                4 * F_0 * sin(omega * self.nds.ivar),
                5 * F_0 * sin(omega * self.nds.ivar),
                6 * F_0 * sin(omega * self.nds.ivar)
            ]
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
                 l=Symbol('l', positive=True),
                 I=Symbol('I', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 input_displacement=dynamicsymbols('theta'),
                 phi_1=dynamicsymbols('\\varphi_1'),                 
                 phi_2=dynamicsymbols('\\varphi_2'),                 
                 phi=dynamicsymbols('\\varphi'),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2'),
                 **kwargs):


        theta = input_displacement
        
        self.phi_1=phi_1
        self.phi_2=phi_2        
        

        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.I = I  # moment of inertia of a rod
        self.input_displacement = input_displacement
        self.qs = qs
        
        self.phi=phi

        self.disc_1 = Disk(I, pos1=phi_1, qs=qs)
        self.spring_1 = Spring(k_2, phi_1, phi_2, qs=qs)  # left spring
        self.disc_2 = Disk(I, pos1=phi_2, qs=qs)
        self.spring_2 = Spring(k_1, pos1=phi_2, pos2=theta,
                               qs=qs)  # right spring
        system = self.disc_1 + self.disc_2 + self.spring_1 + self.spring_2

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict

    
    def get_default_data(self):


        m0, l0 , G, l = symbols('m_0 l_0 G l', positive=True)
        theta0, Omega = symbols('theta_0, Omega', positive=True)

        default_data_dict = {
            self.I: [S.Half*m0*l**2,S.One*m0*l**2,(S.Half**2)*m0*l**2],
            
            

            self.k_1: [1 * G * l**4 /(10* l), 2 *  G * l**4 /(10* l), S.Half *  G * l**4 /(10* l), 4 *  G * l**4 /(10* l), (S.Half**2) *  G * l**4 /(10* l)],
            self.k_2: [1 * G * l**4 /(10* l), 2 *  G * l**4 /(10* l), S.Half *  G * l**4 /(10* l), 4 *  G * l**4 /(10* l), (S.Half**2) *  G * l**4 /(10* l)],

            l:[1 * l0, 2 * l0, S.Half * l0, 4 * l0, (S.Half**2) * l0],
            
            self.phi_1:[self.phi, 0],
            self.phi_2:[self.phi],
            self.input_displacement:[theta0* cos(Omega * self.ivar) ],
        }

        return default_data_dict
    

class DDoFDampedShaft(ComposedSystem):


    scheme_name = 'ddof_damped_shaft.png'
    real_name = 'ddof_shaft_real.png'

    def __init__(self,
                 I=Symbol('I', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 c_1=Symbol('c_1', positive=True),
                 c_2=Symbol('c_1', positive=True),
                 input_displacement=dynamicsymbols('theta'),
                 ivar=Symbol('t'),
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2'),
                 **kwargs):

        phi1, phi2 = qs
        theta = input_displacement

        self.k_2 = k_2  # left spring
        self.k_1 = k_1  # right spring
        self.c_1 = c_1  # right spring
        self.c_2 = c_2  # right spring
        self.I = I  # moment of inertia of a rod
        self.input_displacement = input_displacement
        self.qs = qs

        self.disc_1 = Disk(I, pos1=phi1, qs=qs)
        self.spring_1 = Spring(k_2, phi1, phi2, qs=qs)  # left spring
        self.disc_2 = Disk(I, pos1=phi2, qs=qs)
        self.spring_2 = Spring(k_1, pos1=phi2, pos2=theta,
                               qs=qs)  # right spring
        self.damper_1 = Damper(c_2, phi1, phi2, qs=qs)  # left spring
        self.damper_2 = Damper(c_1, pos1=phi2, pos2=theta,
                               qs=qs)  # right spring
        system = self.disc_1 + self.disc_2 + self.spring_1 + self.spring_2 + self.damper_1 + self.damper_2

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'Moment of Inertia',
            self.k_1: r'',
            self.k_2: r'',
        }
        return self.sym_desc_dict
    def get_default_data(self):

        I0, k0, lamb = symbols('I_0 k_0 lambda', positive=True)

        default_data_dict = {
            self.k_2: [2 * k0, 4 * k0,6*k0,8*k0,10*k0],
            self.k_1: [k0, 3 * k0,5*k0,7*k0,9*k0],
            self.I: [2 * I0, S.Half * I0, 4 * I0, S.Half**2 * I0,3 * I0,3* S.Half * I0, 9 * I0, 3*S.Half**2 * I0],
            self.c_1: [lamb * self.k_1],
            self.c_2: [lamb * self.k_2],
        }
        return default_data_dict

    def get_random_parameters(self):



        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
            }
          
        return parameters_dict
    
    
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


class SDoFExcitedDampedPendulum(ComposedSystem):

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


class SDoFPendulumKinematicExct(ComposedSystem):

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


class SDoFWinch(ComposedSystem):


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
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2'),
                 **kwargs):

        phi1, phi2 = qs

        self.m = m
        self.g = g
        self.l = l
        self.k = k

        self.spring = Spring(k, pos1=(phi1 * (l/2)), pos2=(phi2 * (l/2)), qs=[qs])
        self.pendulum_1 = Pendulum(m, g, l, angle=phi1, qs=[qs])
        self.pendulum_2 = Pendulum(m, g, l, angle=phi2, qs=[qs])

        system = self.pendulum_1 + self.pendulum_2 + self.spring
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.k: r'Stifness coefficient',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            
            self.l: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],

#             self.phi_1: [self.phi_l, 0],
#             self.phi_2: [self.phi_l, self.phi_r, 0],
#             self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict
    
class DDoFLinearizedCoupledPendulum(ComposedSystem):
    """
    Model of a DDoF Linearized Coupled Pendulum.

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
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2'),
                 **kwargs):

        phi1, phi2 = qs

        self.m = m
        self.g = g
        self.l = l
        self.k = k

        self.spring = Spring(k, pos1=(phi1 * (l/2)), pos2=(phi2 * (l/2)), qs=[qs])
        self.pendulum_1 = Pendulum(m, g, l, angle=phi1, qs=[qs]).linearized()
        self.pendulum_2 = Pendulum(m, g, l, angle=phi2, qs=[qs]).linearized()

        system = self.pendulum_1 + self.pendulum_2 + self.spring
        super().__init__(system,**kwargs)
        
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.k: r'Stifness coefficient',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
#             self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            
            self.l: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],

#             self.phi_1: [self.phi_l, 0],
#             self.phi_2: [self.phi_l, self.phi_r, 0],
#             self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict
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


class SDoFNonlinearEngine(ComposedSystem):
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
class SDoFStraightNonlinearEngine(SDoFNonlinearEngine):
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
                 **kwargs):

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
        super().__init__(system,**kwargs)

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
                 **kwargs):

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

        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, l0, I0, k0, r0 = symbols('m_0 l_0 I_0 k_0 r_0', positive=True)

        default_data_dict = {
            self.m: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
            self.l: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.I: [2 * I0, 3 * I0, 4 * I0, 5 * I0, 6 * I0],
            self.k: [2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0],
            self.r: [2 * r0, 3 * r0, 4 * r0, 5 * r0, 6 * r0],
        }
        return default_data_dict

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

    scheme_name = 'elastic_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    def __init__(self,
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('\\varphi'),
                 **kwargs):

        self.k = k
        self.l = l
        self.m = m
        self.g = g
        self.phi = phi
        self.z = z

        x = (l + z) * sin(phi)
        y = (l + z) * cos(phi)

        self.frame = base_frame

        self.payload = Point('payload')
#         self.payload.set_vel(
#             self.frame, (sqrt((diff(x, ivar)**2 + diff(y, ivar)**2).simplify()))*self.frame.x)
        self.payload.set_vel(
            self.frame,
            sqrt(diff(z, ivar)**2 + (diff(phi, ivar) * (l + z))**2) *
            self.frame.x)

#         print(self.payload, 'try', type(self.payload))

        self.spring = Spring(k, z, qs=[phi, z])
        self.material_point_1 = HarmonicOscillator(
                                              S.Half*m * (diff(z, ivar)**2 + (diff(phi, ivar) * (l + z))**2),
                                              qs=[phi, z],
                                             )
        # self.material_point_2 = MaterialPoint(m, y, qs=[phi, z])
        self.gravity = GravitationalForce(m, g, pos1=-y, qs=[phi, z])
        system = (self.spring + self.gravity + self.material_point_1
                  )  # + self.material_point_2

        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.k: r'Spring stiffness',
            self.l: r'Pendulum length',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, k0 = symbols('m_0 l_0 k_0', positive=True)

        default_data_dict = {
            self.m: [1 * m0, S.Half**2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.l: [1 * l0, S.Half * l0, 2 * l0, 1 * l0, 2 * l0],
            self.k: [S.Half * k0, 1 * k0, 2 * k0, S.Half * k0, 2 * k0]
        }
        return default_data_dict
    
    def linearized(self, x0=None, op_point=True, hint=None, label=None):

        if hint is None:
            hint=[self.phi]
        
        return super().linearized(x0=x0, op_point=True, hint=hint, label=label)

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

    scheme_name = 'damped_elastic_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    def __init__(self,
                 undamped_system,
                 c=Symbol('c', positive=True),
                 k=Symbol('k', positive=True),
                 l=Symbol('l', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 ivar=Symbol('t'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('\\varphi'),
                 **kwargs):

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

        #         display(self.damper._eoms)

        system = self.undamped + self.damper

        super().__init__(system,**kwargs)

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


class SDoFTrolleyWithNonlinearSpring(ComposedSystem):
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
        self.Force = Force(-F * cos(Omega * ivar), pos1=x, qs=[x])

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


class DDoFTwoDisksWithThreeSprings(ComposedSystem):
    scheme_name = 'ddof_disks_3_springs_scheme.png'
    real_name = 'nonlin_trolley_real.PNG'

    def __init__(self,
                 d=Symbol('d', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 kl=Symbol('k_l', positive=True),
                 kc=Symbol('k_c', positive=True),
                 kr=Symbol('k_r', positive=True),
                 R=Symbol('R', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 xl=dynamicsymbols('x_l'),
                 xr=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l, x_r'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.kl = kl
        self.kc = kc
        self.kr = kr
        self.R = R
        self.l_0 = l_0
        self.d = d
        self.xl = xl
        self.xr = xr
        self.x = x

        self.Disk1 = MaterialPoint(m1, xl, qs=[xl]) + MaterialPoint(m1/2*R**2, xl/R, qs=[xl]) + Spring(kl, xl, qs=[xl])
        self.Disk2 = MaterialPoint(m2, xr, qs=[xr]) + MaterialPoint(m2/2*R**2, xr/R, qs=[xr]) + Spring(kr,xr, qs=[xr])
        self.Spring = Spring(kc, xl, xr, qs=[xl, xr])

        system = self.Disk1 + self.Spring + self.Disk2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

#             self.d: [1 * l0, 2 * l0, S.Half * l0, 4 * S.Half * l0,  S.Half**2 * l0],

            self.kl: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
            self.kr: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
            self.kc: [S.Half * k0, S.Half**2 * k0, 1 * k0, 4 * S.Half * k0, 2 * k0],
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
            self.m1: r'Disk Mass',
            self.m2: r'Disk Mass',
            self.kl: 'Left Spring Stiffness',
            self.kr: 'Right Spring Stiffness',
            self.kc: 'Central Spring Stiffness',
            self.l: r'Length',
            self.l_0: r'initial Spring Length',
        }
        return self.sym_desc_dict


class DDoFTwoNonLinearTrolleys(ComposedSystem):
    scheme_name = 'ddof_nonlin_trolleys.PNG'
    real_name = 'dwa_wozki_XD.PNG'

    def __init__(self,
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k1=Symbol('k_1', positive=True),
                 k2=Symbol('k_2', positive=True),
                 k3=Symbol('k_3', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 x1=dynamicsymbols('x1'),
                 x2=dynamicsymbols('x2'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x1, x2'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.d = d
        self.l_0 = l_0
        self.x1 = x1
        self.x2 = x2
        self.x = x

        self.Trolley1 = MaterialPoint(m1, x1, qs=[x1]) + Spring(
            k1, pos1=(sqrt(x1**2 + d**2) - l_0), qs=[x1])
        self.Trolley2 = MaterialPoint(m2, x2, qs=[x2]) + Spring(
            k2, pos1=(sqrt(x2**2 + d**2) - l_0), qs=[x2])
        self.Spring = Spring(k3, x1, x2,qs=[x1,x2])

        system = self.Trolley1 + self.Spring + self.Trolley2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.d: [1 * l0, 2 * l0, S.Half * l0, 3 * S.Half * l0, 1 * l0],
            self.k1:
            [S.Half * k0, S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0],
            self.k2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k3:
            [S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0, 5 * S.Half * k0],
            self.x1: [self.x, 0],
            self.x2: [self.x, S.Zero],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x1] == S.Zero:
            parameters_dict[self.x2] = self.x

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Trolley Mass',
            self.m2: r'Trolley Mass',
            self.k1: 'Spring Stiffness',
            self.k2: 'Spring Stiffness',
            self.k3: 'Spring Stiffness',
            self.d: r'length',
            self.l_0: r'length',
        }
        return self.sym_desc_dict

class SDoFNonLinearTrolley(ComposedSystem):

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

class MDoFForcedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'mdof_three_trolleys.PNG'
    real_name = 'three_carriages.PNG'

    def __init__(self,

                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_12=Symbol('k_12', positive=True),
                 k_c12=Symbol('k_c12', positive=True),
                 k_23=Symbol('k_23', positive=True),
                 k_c23=Symbol('k_c23', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 F=Symbol('F', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 
                 x_l=dynamicsymbols('x_l'),
                 x_c=dynamicsymbols('x_c'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x_3=dynamicsymbols('x_3'),
                 
                 qs=dynamicsymbols('x_l x_c x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_12 = k_12
        self.k_c12 = k_c12
        self.k_23 = k_23
        self.k_c23 = k_c23
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_c = x_c
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x_3 = x_3
        self.Omega = Omega

        self.Trolley_1 = MaterialPoint(m1, x_l, qs=[x_l]) + Spring(
            k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[
                x_l
            ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(
                -F * cos(Omega * ivar), pos1=x_l, qs=[x_l])
        
        
        
        self.Trolley_2 = MaterialPoint(m2, x_c, qs=qs) + Spring(
            k_12, x_l, x_c,qs=qs) + Spring(k_12, x_l, x_c,qs=qs) + Spring(
                k_23, x_c, x_r,qs=qs) + Spring(k_23, x_c, x_r,qs=[x_c,x_l]) + Spring(
                    k_c12, x_l, x_c,qs=qs) + Spring(k_c23, x_c, x_r,qs=qs)
        self.Trolley_3 = MaterialPoint(m3, x_r, qs=qs) + Spring(
            k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[
                x_r
            ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(
                -2 * F * cos(Omega * ivar), pos1=x_r, qs=[x_r])

        system = self.Trolley_1 + self.Trolley_2 + self.Trolley_3
        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2, 0],
            self.x_r: [self.x_2, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict



class CSBeam(ContinuousSystem):
    """
    Model of a Continuous Beam

        Arguments:
        =========
            m = Mass
                -Mass of the payload



        Example
        =======
        A vibrating beam for given Young modulus and moment of inertia can be created in the following way.

        >>> t = symbols('t')
        >>> R,l_0 = symbols('R, l_0',positive=True)
        >>> MDoFWinch(r=R,l=l_0)

        -define the symbols and dynamicsymbols
        -determine the instance of the beam by using class CSBeam()
    """

    scheme_name = 'supported_beam_scheme.PNG'
    real_name = 'supported_beam_real.PNG'

    def __init__(self,
                E=Symbol('E',positive=True),
                A=Symbol('A',positive=True),
                I=Symbol('I',positive=True),
                rho=Symbol('\\rho',positive=True),
                w=Function('w'),
                bc_dict=None,
                time=Symbol('t'),
                loc=Symbol('x'),
                **kwargs
                ):

        self.E = E
        self.A = A
        self.I = I
        self.rho = rho
        self.w=w(time,loc)
        self.time=time
        self.loc=loc
        self.l=Symbol('L',positive=True)
        

        L_beam=S.One/2*(A*rho*(self.w.diff(self.time))**2-E*I*(self.w.diff(self.loc,2))**2)

        super().__init__(L_beam,q=self.w,bc_dict=bc_dict,t_var=self.time, spatial_var=self.loc,**kwargs)
        
        
        self._sep_expr=2*self.L.subs({self.w.diff(self.time):0,(self.w.diff(self.loc,2)):1}) .doit()  *Symbol('k',positive=True)**4 
        
        print(self._sep_expr)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.E: r'Young modulus',
            self.E: r'Young modulus',
            self.I: r'Cross-section moment of inertia',
            self.A: r'Area of cross-section',
            self.rho: r'density',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        E_0, A_0, I_0, L_0 = symbols('E_0, A_0, I_0, L_0', positive=True)
        
        
        l=self.l
        x=self.loc
        X=Function('X')(x)
        
        
        sup_ls=Subs(X.diff(x,0),x,0)
        sup_rs=Subs(X.diff(x,0),x,l)

        fix_ls=Subs(X.diff(x,1),x,0)
        fix_rs=Subs(X.diff(x,1),x,l)

        mb_ls=Subs(X.diff(x,2),x,0)
        mb_rs=Subs(X.diff(x,2),x,l)

        v_ls=Subs(X.diff(x,3),x,0)
        v_rs=Subs(X.diff(x,3),x,l)
        
        

        default_data_dict = {
            self.E: [2.5*E_0,1.25*E_0,0.75*E_0,1.35*E_0,1*E_0,2*E_0,3*E_0,4*E_0,5*E_0,6*E_0,7*E_0,8*E_0,9*E_0 ],
            self.A: [1 * A_0, 2 * A_0, S.Half * A_0, 1 * A_0, 2 * A_0,3*A_0,4*A_0,5*A_0,6*A_0,7*A_0,8*A_0,9*A_0,10*A_0,11*A_0],
            self.I: [1 * I_0, 2 * I_0, S.Half * I_0, 1 * I_0, 2 * I_0,3*I_0,4*I_0,5*I_0,6*I_0,7*I_0,8*I_0,9*I_0,10*I_0,11*I_0],
           self.BC: [ {sup_ls:0,mb_ls:0,sup_rs:0,mb_rs:0},{fix_ls:0,v_ls:0,sup_rs:0,mb_rs:0}, ],
           self.l:[1 * L_0, 2 * L_0, S.Half * L_0, 1 * L_0, 3 * L_0,4*L_0,5*L_0,6*L_0,7*L_0,8*L_0,9*L_0,10*L_0,11*L_0],

        }

        return default_data_dict

    def get_random_parameters(self):

        E_0, A_0, I_0, L_0 = symbols('E_0, A_0, I_0, L_0', positive=True)

        data_dict=super().get_random_parameters()
        data_dict[self.A]  = (L_0**2 * random.choice([0.125, 1.25, 0.0125, 1.23, 0.128 ])/ (data_dict[self.E]/E_0) / (data_dict[self.I]/I_0) ).n(3)

        return data_dict
    
    
class CSRod(ContinuousSystem):
    """
    Model of a Continuous Rod

        Arguments:
        =========
            E = Young modulus
                -Young modulus of a rod's material
            A = Area of cross-section

            rho_L = Linear density of a rod's material

            u = Space-Time function
                -Allow to describe a rod vibration in two-variable state

        Example
        =======
        A vibrating rod for given Young modulus and area of cross-section can be created in the following way.

        >>> t = symbols('t')
        >>> E,A,rho_L = symbols('E,A,rho_L',positive=True)
        >>> CSRod()

        -define the symbols and dynamicsymbols
        -determine the instance of the beam by using class CSRod()
    """

    scheme_name = 'rod_scheme.PNG'
    real_name = 'rod_real.PNG'

    def __init__(self,
                E=Symbol('E',positive=True),
                A=Symbol('A',positive=True),
                rho_L=Symbol('\\rho_L',positive=True),
                u=Function('u'),
                bc_dict=None,
                time=Symbol('t'),
                loc=Symbol('x'),
                **kwargs
                ):

        self.E = E
        self.A = A
        self.rho_L = rho_L
        self.u=u(time,loc)
        self.time=time
        self.loc=loc
        self.l=Symbol('L',positive=True)
        
        

        L_rod=S.One/2*(rho_L*(self.u.diff(self.time))**2-A*E*(self.u.diff(self.loc))**2)

        super().__init__(L_rod,q=self.u,bc_dict=bc_dict,t_var=self.time, spatial_var=self.loc,**kwargs)

        self._sep_expr=2*self.L.subs({self.u.diff(self.time):0,(self.u.diff(self.loc)):1}) .doit()  *Symbol('k',positive=True)**2 
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.E: r'Young modulus',
            self.A: r'Area of cross-section',
            self.rho_L: r'Linear density',
            self.u: 'Space-Time function',
        }
        return self.sym_desc_dict


    def get_default_data(self):

        E_0, A_0, L_0 = symbols('E_0, A_0, L_0', positive=True)
        
        x=self.loc
        l=self.l
        X=Function('X')(x)
        
        fix_ls=Subs(X.diff(x,0),x,0)
        fix_rs=Subs(X.diff(x,0),x,l)
        free_ls=Subs(X.diff(x,1),x,0)
        free_rs=Subs(X.diff(x,1),x,l)


        default_data_dict = {
            self.E: [2.5*E_0,1.25*E_0,0.75*E_0,1.35*E_0 ],
            self.A: [1 * A_0, 2 * A_0, S.Half * A_0, 1 * A_0, 2 * A_0],
            
           self.BC: [ {fix_ls:0,free_rs:0},{free_ls:0,free_rs:0},{fix_ls:0,fix_rs:0} ],
           self.l:[1 * L_0, 2 * L_0, S.Half * L_0, 1 * L_0, 2 * L_0],
            }

        return default_data_dict

    def get_random_parameters(self):

        E_0, A_0, L_0 = symbols('E_0, A_0, L_0', positive=True)

        data_dict=super().get_random_parameters()
        data_dict[self.A]  = (L_0**2 * random.choice([0.125, 0.0125, 0.00125, 0.123, 0.0128 ])/ (data_dict[self.E]/E_0) ).n(3)

        return data_dict

class CSString(ContinuousSystem):


    scheme_name = 'string_scheme.png'
    real_name = 'double_bass_string.jpg'

    def __init__(self,
                Ty=Symbol('T_y',positive=True),
                A=Symbol('A',positive=True),
                rho=Symbol('\\rho',positive=True),
                w=Function('w'),
                bc_dict=None,
                time=Symbol('t'),
                loc=Symbol('x'),
                **kwargs
                ):

        self.Ty = Ty
        self.A = A
        self.rho = rho
        self.w=w(time,loc)
        self.time=time
        self.loc=loc
        self.l=Symbol('L',positive=True)
        
        

        L_rod=S.One/2*(A*rho*(self.w.diff(self.time))**2-Ty*(self.w.diff(self.loc))**2)

        super().__init__(L_rod,q=self.w,bc_dict=bc_dict,t_var=self.time, spatial_var=self.loc,**kwargs)

        self._sep_expr=2*self.L.subs({self.w.diff(self.time):0,(self.w.diff(self.loc)):1}) .doit()  *Symbol('k',positive=True)**2 
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.Ty: r'Initial tension',
            self.A: r'Area of cross-section',
            self.rho: r'Material density',
            self.w: 'Space-Time function',
        }
        return self.sym_desc_dict


    def get_default_data(self):

        T0, A_0, L_0 = symbols('T_0, A_0, L_0', positive=True)
        
        x=self.loc
        l=self.l
        X=Function('X')(x)
        
        fix_ls=Subs(X.diff(x,0),x,0)
        fix_rs=Subs(X.diff(x,0),x,l)
        free_ls=Subs(X.diff(x,1),x,0)
        free_rs=Subs(X.diff(x,1),x,l)


        default_data_dict = {
            self.Ty: [2.5*T0,1.25*T0,0.75*T0,1.35*T0 ],
            self.A: [1 * A_0, 2 * A_0, S.Half * A_0, 1 * A_0, 2 * A_0],
            
           self.BC: [{fix_ls:0,fix_rs:0} ],
           self.l:[1 * L_0, 2 * L_0, S.Half * L_0, 3 * L_0, 2 * L_0,4*L_0,5*L_0,6*L_0,7*L_0,8*L_0,9*L_0],

        }

        return default_data_dict

    def get_random_parameters(self):
        
        T0, A_0, L_0 = symbols('T_0, A_0, L_0', positive=True)
        
        data_dict=super().get_random_parameters()
        data_dict[self.A]  = (L_0**2 * random.choice([0.03, 0.031, 0.0031, 0.033, 3.1 ])/ (data_dict[self.Ty]/T0) ).n(3)
        
        
        return data_dict

                                                  
class CSShaft(ContinuousSystem):


    scheme_name = 'supported_shaft_scheme.PNG'
    real_name = 'shaft_real.PNG'

    def __init__(self,
                G=Symbol('G',positive=True),

                M=Symbol('M',positive=True),
                I=Symbol('I',positive=True),

                rho=Symbol('\\rho',positive=True),
                phi=Function('\\phi'),
                bc_dict=None,
                time=Symbol('t'),
                loc=Symbol('x'),
                **kwargs
                ):


        self.M = M
        self.G = G
        self.I = I

        self.rho = rho
        self.phi=phi(time,loc)
        self.time=time
        self.loc=loc
        self.l=Symbol('L',positive=True)
        
        


        L_shaft=S.One/2*(rho*I*(self.phi.diff(self.time))**2-G*I*(self.phi.diff(self.loc))**2)


        super().__init__(L_shaft,q=self.phi,bc_dict=bc_dict,t_var=self.time, spatial_var=self.loc,**kwargs)

        self._sep_expr=2*self.L.subs({self.phi.diff(self.time):0,(self.phi.diff(self.loc)):1}) .doit()  *Symbol('k',positive=True)**2 
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.G: r'''Kirchoff's modulus ''',
            
            self.rho: r'Material density',
            self.phi: 'Space-Time function',
        }
        return self.sym_desc_dict


    def get_default_data(self):

        G_0, L_0 = symbols('G_0, L_0', positive=True)
        
        x=self.loc
        l=self.l
        X=Function('X')(x)
        
        fix_ls=Subs(X.diff(x,0),x,0)
        fix_rs=Subs(X.diff(x,0),x,l)
        free_ls=Subs(X.diff(x,1),x,0)
        free_rs=Subs(X.diff(x,1),x,l)


        default_data_dict = {
            self.G: [2.5*G_0,1.25*G_0,0.75*G_0,1.35*G_0 ],
            

           self.BC: [{fix_ls:0,free_rs:0},{fix_ls:0,fix_rs:0} ],
           self.l:[1 * L_0, 2 * L_0, S.Half * L_0, 3 * L_0, 2 * L_0],


        }

        return default_data_dict

    def get_random_parameters(self):
        
        G_0, L_0 = symbols('G_0, L_0', positive=True)
        
        data_dict=super().get_random_parameters()

        data_dict[self.I]  = (L_0**4 * random.choice([0.1, 0.01, 0.011, 0.11, 1.1,11 ])/ (data_dict[self.G]/G_0) ).n(3)

        
        
        return data_dict


class CSCylinder(PlaneStressProblem):
    
    scheme_name = 'cylinder_scheme.PNG'
    real_name = 'cylinder_real_cannon.PNG'
    
    def __init__(self,
                 
                 disp_func=[Function('\\mathit{u}')(Symbol('r')), 0],
                 stress_tensor=Matrix(2, 2, [
                     Function('\\sigma_r')(Symbol('r')), 0, 0,
                     Function('\\sigma_\\varphi')(Symbol('r'))
                 ]),

                 bc_dict=None,
                 coords=[Symbol('r'), Symbol('\\varphi')],
                 E=Symbol('E', positive=True),
                 nu=Symbol('\\nu', positive=True),
                 D=Symbol('D', positive=True),
                 volumetric_load=0,
                 **kwargs
                ):
        


        
        super().__init__(disp_func=disp_func,stress_tensor=stress_tensor,bc_dict=bc_dict,coords=coords,E=E,nu=nu,D=D,volumetric_load=volumetric_load,**kwargs)
        

        
class CSPlate(PlaneStressProblem):
    
    scheme_name = 'plate_point_load.PNG'
    real_name = 'reservoir.jpg'
    
    def __init__(self,
                 disp_func=[Function('\\psi')(Symbol('r')), 0],
                 stress_tensor=Matrix(2, 2, [
                     Function('\\mathit{m}_r')(Symbol('r')), 0, 0,
                     Function('\\mathit{m}_\\varphi')(Symbol('r'))
                 ]),
                 bc_dict=None,
                 coords=[Symbol('r'), Symbol('\\varphi')],
                 E=Symbol('E', positive=True),
                 nu=Symbol('\\nu', positive=True),      
                 D=Symbol('D_h', positive=True),
                 h=Symbol('h', positive=True),
                 volumetric_load=0,
                 **kwargs
                ):

#         print('__init')
#         print(type(self))
#         display(disp_func)
        
        super().__init__(disp_func=disp_func,stress_tensor=stress_tensor,bc_dict=bc_dict,coords=coords,E=E,nu=nu,D=D,volumetric_load=volumetric_load,**kwargs)
        
#         print('__init - after super')
#         print(type(self))
#         display(self.u)
        self._h=h
        
        self._subs_dict = {E:D*(1-nu**2)*12/(h**3)} 
        self.E_module=(h**3)/12*E
#         self.u = Matrix([Function('\\psi')(Symbol('r')), 0])
        
#         print('__init - after super')
#         print(type(self))
#         display(self.u)

# <<<<<<< HEAD
# OK
class MDoFDoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'DoubleTrolley_real.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m_w1=Symbol('m_w1', positive=True),
                 m_w2=Symbol('m_w2', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 lam=Symbol('lambda', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.x_l = x_l
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x=x
        self.m1 = m1
        self.m2 = m2
        self.m_w1 = m_w1
        self.m_w2 = m_w2
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.R = R
        self.lam = lam

        self.trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m_w1, x_l/2, qs = [x_l]) + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]) + MaterialPoint(m_w1, x_l/2, qs = [x_l])
                          + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]))

        self.trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r])
                          + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]))

        system = self.trolley_1 + self.trolley_2
        super().__init__(system(qs),**kwargs)
        
        

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            
            
            self.lam:[0],


        }

        return default_data_dict



    
class SDoFDoubleTrolleyDifferentWheels(ComposedSystem):
    scheme_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'
    real_name = 'MDOFDoubleTrolleyDifferentWheels.PNG'

    def __init__(self,
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 m_w1=Symbol('m_w1', positive=True),
                 m_w2=Symbol('m_w2', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 R=Symbol('R', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 lam=Symbol('lambda', positive=True),

                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m_w1 = m_w1
        self.m_w2 = m_w2
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.R = R
        self.lam = lam

        self.trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m_w1, x_l/2, qs = [x_l]) + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]) + MaterialPoint(m_w1, x_l/2, qs = [x_l])
                          + MaterialPoint(m_w1/2, x_l/2, qs = [x_l]))

        self.trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]) + MaterialPoint(m_w2, x_r/2, qs = [x_r])
                          + MaterialPoint(m_w2/2, x_r/2, qs = [x_r]))

        system = self.trolley_1 + self.trolley_2
        super().__init__(system(qs),**kwargs)

        
    def get_default_data(self):

        m0, k0, l0, lam = symbols('m_0 k_0 l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w1: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.m_w2: [1 * m0, 2 * m0, S.Half * m0, S.Half**2 *  m0, 2**2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict

    
    

class MDoFDampedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'MDOF_Damped_Trolleys_With_Springs.PNG'
    real_name = 'two_damped_trolleys.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 c_cl=Symbol('c_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 c_cc=Symbol('c_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 c_cr=Symbol('c_cr', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m = m
        self.k_l = k_l
        self.c_cl = c_cl
        self.k_c = k_c
        self.c_cc = c_cc
        self.k_r = k_r
        self.c_cr = c_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_l, pos1=x_l, qs=[x_l]) + Damper(c_cl, pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))

        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Damper(c_cc, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) 
                          + Damper(c_cr, pos1=x_r, qs=[x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) 
                          + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0, lam = symbols('m k l_0 lambda', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.c_cr: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cc: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],
            self.c_cl: [lam *  k0,lam * 2 * k0,lam * S.Half * k0,lam * 4 * k0,lam * S.Half**2 * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict


    
    
# OK
class MDoFThreePendulumsWithSprings(ComposedSystem):
    scheme_name = 'MDOF_Three_Pendulum_Spring.png'
    real_name = 'three_carriages.PNG'

    def __init__(self,
                 phi_1=dynamicsymbols('varphi_1'),
                 phi_2=dynamicsymbols('varphi_2'),
                 phi_3=dynamicsymbols('varphi_3'),
                 phi_l=dynamicsymbols('varphi_l'),
                 phi_r=dynamicsymbols('varphi_r'),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_3=Symbol('k_3', positive=True),
                 k_4=Symbol('k_4', positive=True),
                 l=Symbol('l', positive=True),
                 g=Symbol('g', positive=True),
                 qs=dynamicsymbols('varphi_1, varphi_2, varphi_3'),
                 ivar=Symbol('t'),
                 system = None,
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_1 = k_1
        self.k_2 = k_2
        self.k_3 = k_3
        self.k_4 = k_4
        self.l = l
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.phi_l = phi_l
        self.phi_r = phi_r
        self.g = g

        if system == None:
            self.Pendulum1 = Pendulum(m1, g, l, angle=phi_1, qs=[phi_1])
            self.Pendulum2 = Pendulum(m2, g, l, angle=phi_2, qs=[phi_2])
            self.Pendulum3 = Pendulum(m3, g, l, angle=phi_3, qs=[phi_3])
            self.Spring1 = Spring(k_1, pos1=(sqrt((l/2*(sin(phi_1)-sin(phi_2)))**2+(l/2*(cos(phi_1)-cos(phi_2)))**2)), qs=[phi_1, phi_2])
            self.Spring2 = Spring(k_2, pos1=(sqrt((l/2*(sin(phi_2)-sin(phi_3)))**2+(l/2*(cos(phi_2)-cos(phi_3)))**2)), qs=[phi_2, phi_3])
            self.Spring3 = Spring(k_3, pos1=(sqrt((l*(sin(phi_1)-sin(phi_2)))**2+(l*(cos(phi_1)-cos(phi_2)))**2)), qs=[phi_1, phi_2])
            self.Spring4 = Spring(k_4, pos1=(sqrt((l*(sin(phi_2)-sin(phi_3)))**2+(l*(cos(phi_2)-cos(phi_3)))**2)), qs=[phi_2, phi_3])

            system = self.Pendulum1 + self.Pendulum2 + self.Pendulum3 + self.Spring1 + self.Spring2 + self.Spring3 + self.Spring4
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)


        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.k_1: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.phi_1: [self.phi_l, 0],
            self.phi_2: [self.phi_l, self.phi_r, 0],
            self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi_1] != self.phi_l or parameters_dict[self.phi_2] != self.phi_l:

            parameters_dict[self.phi_1] = self.phi_l

        if parameters_dict[self.phi_2] != self.phi_r or parameters_dict[self.phi_3] != self.phi_r:

            parameters_dict[self.phi_2] = self.phi_r


        return parameters_dict

    
class MDoFLinearizedThreePendulumsWithSprings(ComposedSystem):
    scheme_name = 'three_pendulums_forced.PNG'
    real_name = 'lifting_tandem.png'

    def __init__(self,
                 phi_1=dynamicsymbols('\\varphi_1'),
                 phi_2=dynamicsymbols('\\varphi_2'),
                 phi_3=dynamicsymbols('\\varphi_3'),
                 phi_l=dynamicsymbols('\\varphi_l'),
                 phi_c=dynamicsymbols('\\varphi_c'),
                 phi_r=dynamicsymbols('\\varphi_r'),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_1=Symbol('k_1', positive=True),
                 k_2=Symbol('k_2', positive=True),
                 k_3=Symbol('k_3', positive=True),
                 k_4=Symbol('k_4', positive=True),
                 l=Symbol('l', positive=True),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('\Omega', positive=True),
                 F=Symbol('F',positive=True),
                 qs=dynamicsymbols('\\varphi_l, \\varphi_c, \\varphi_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.k_1 = k_1
        self.k_2 = k_2
        self.k_3 = k_3
        self.k_4 = k_4
        self.l = l
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.phi_l = phi_l
        self.phi_c = phi_c
        self.phi_r = phi_r
        self.g = g
        self.Omega = Omega
        self.F = F

        self.Pendulum1 = Pendulum(m1, g, l, angle=phi_l, qs=[phi_l]).linearized() + Force(-F * l * cos(Omega * ivar), pos1=phi_l, qs=[phi_l])
        self.Pendulum2 = Pendulum(m2, g, l, angle=phi_c, qs=[phi_c]).linearized()
        self.Pendulum3 = Pendulum(m3, g, l, angle=phi_r, qs=[phi_r]).linearized() + Force(-2* F * l* cos(Omega * ivar), pos1=phi_r, qs=[phi_r])
        self.Spring1 = Spring(k_1, pos1=(phi_l * (l/2)), pos2=(phi_c * (l/2)), qs=[phi_l, phi_c])
        self.Spring2 = Spring(k_2, pos1=(phi_c * (l/2)), pos2=(phi_r * (l/2)), qs=[phi_c, phi_r])
        self.Spring3 = Spring(k_3, pos1=(phi_l * l), pos2=(phi_c * l), qs=[phi_l, phi_c])
        self.Spring4 = Spring(k_4, pos1=(phi_c * l), pos2=(phi_r * l), qs=[phi_c, phi_r])

        system = self.Pendulum1 + self.Pendulum2 + self.Pendulum3 + self.Spring1 + self.Spring2 + self.Spring3 + self.Spring4
        super().__init__(system(qs),**kwargs)


    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k_1: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_3: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_4: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.phi_1: [self.phi_l, 0],
            self.phi_2: [self.phi_l, self.phi_r, 0],
            self.phi_3: [self.phi_r, 0],
        }

        return default_data_dict





    
class DDoFTwoNonLinearDisks(ComposedSystem):
    scheme_name = 'MDOF_Double_Disk.png'
    real_name = 'dwa_wozki_XD.PNG'

    def __init__(self,
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 kl=Symbol('k_l', positive=True),
                 kc=Symbol('k_c', positive=True),
                 kr=Symbol('k_r', positive=True),
                 R=Symbol('R', positive=True),
                 d=Symbol('d', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 ivar=Symbol('t'),
                 xl=dynamicsymbols('x_l'),
                 xr=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l, x_r'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.kl = kl
        self.kc = kc
        self.kr = kr
        self.R = R
        self.l_0 = l_0
        self.d = d
        self.xl = xl
        self.xr = xr
        self.x = x

        self.Disk1 = MaterialPoint(m1, xl, qs=[xl]) + MaterialPoint(m1/2*R**2, xl/R, qs=[xl]) + Spring(kl, pos1=(sqrt(xl**2 + d**2) - l_0), qs=[xl])
        self.Disk2 = MaterialPoint(m2, xr, qs=[xr]) + MaterialPoint(m2/2*R**2, xr/R, qs=[xr]) + Spring(kr, pos1=(sqrt(xr**2 + d**2) - l_0), qs=[xr])
        self.Spring = Spring(kc, xl, xr, qs=[xl, xr])

        system = self.Disk1 + self.Spring + self.Disk2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.d: [1 * l0, 2 * l0, S.Half * l0, 3 * S.Half * l0, 1 * l0],

            self.kl: [S.Half * k0, S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0],
            self.kr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.kc: [S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0, 5 * S.Half * k0],

            self.xl: [self.x, 0],
            self.xr: [self.x, S.Zero],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }


        if parameters_dict[self.xl] == S.Zero:
            parameters_dict[self.xr] = self.x

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Disk Mass',
            self.m2: r'Disk Mass',
            self.kl: 'Left Spring Stiffness',
            self.kr: 'Right Spring Stiffness',
            self.kc: 'Central Spring Stiffness',
            self.l: r'Length',
            self.l_0: r'initial Spring Length',
        }
        return self.sym_desc_dict



class MDoFDoubleTrolleyWithNonlinearSprings(ComposedSystem):
    scheme_name = 'MDOF_Double_Trolley_With_Springs.PNG'
    real_name = 'sdof_nonlinspring_trolleys_real.PNG'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 k_cc=Symbol('k_cc', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 mu=Symbol('mu', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_r=dynamicsymbols('x_r'),
                 x=dynamicsymbols('x'),
                 qs=dynamicsymbols('x_l x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_c = k_c
        self.k_cc = k_cc
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_r = x_r
        self.x = x
        self.R = R
        self.mu = mu

        self.Trolley_1 = (MaterialPoint(m1, x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) 
                          + Spring(k_cl*(1-mu*x_l**2), pos1=x_l, qs=[x_l]) 
                          + MaterialPoint(m, x_l/2, qs = [x_l]) + MaterialPoint(m/2, x_l/2, qs = [x_l]) + MaterialPoint(m, x_l/2, qs = [x_l]) 
                          + MaterialPoint(m/2, x_l/2, qs = [x_l]))
        
        
       
        self.Trolley_2 = (MaterialPoint(m2, x_r, qs=[x_r]) + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_c, pos1=x_l, pos2=x_r, qs=[x_l, x_r]) + Spring(k_cc*(1-mu*(x_r - x_l)**2), pos1=x_l, pos2=x_r, qs=[x_l, x_r]) 
                          + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_cr*(1-mu*x_r**2), pos1=x_r, qs=[x_r]) 
                          + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]) + MaterialPoint(m, x_r/2, qs = [x_r]) + MaterialPoint(m/2, x_r/2, qs = [x_r]))

        system = self.Trolley_1 + self.Trolley_2
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cc: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x, 0],
            self.x_r: [self.x, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] == 0 and parameters_dict[self.x_r]==0:

            parameters_dict[self.x_l] = self.x


        return parameters_dict




class MDoFForcedDisksWithParallelSprings(ComposedSystem):
    scheme_name = 'MDOF_Forced_Disks_With_Parallel_Springs.PNG'
    real_name = 'three_rollers_real.png'

    def __init__(self,
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 k_cl=Symbol('k_cl', positive=True),
                 k_12=Symbol('k_12', positive=True),
                 k_c12=Symbol('k_c12', positive=True),
                 k_23=Symbol('k_23', positive=True),
                 k_c23=Symbol('k_c23', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_cr=Symbol('k_cr', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_c=dynamicsymbols('x_c'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 qs=dynamicsymbols('x_l x_c x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.R = R
        self.k_l = k_l
        self.k_cl = k_cl
        self.k_12 = k_12
        self.k_c12 = k_c12
        self.k_23 = k_23
        self.k_c23 = k_c23
        self.k_r = k_r
        self.k_cr = k_cr
        self.x_l = x_l
        self.x_c = x_c
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.Omega = Omega

        self.Disk1 = MaterialPoint(m, x_l, qs=[x_l]) + MaterialPoint(
            m / 2, x_l, qs=[x_l]) + MaterialPoint(m1, x_l, qs=[x_l]) + Spring(
                k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[
                    x_l
                ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(
                    -F_0 * cos(Omega * ivar), pos1=x_l, qs=[x_l])
        self.Disk2 = MaterialPoint(m, x_c, qs=[x_c]) + MaterialPoint(
            m / 2, x_c, qs=[x_c]) + MaterialPoint(m2, x_c, qs=[
                x_c
            ]) + Spring(k_12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(
                k_c12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(
                    k_12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(
                        k_23, pos1=x_c, pos2=x_r, qs=[x_c, x_r]) + Spring(
                            k_c23, pos1=x_c, pos2=x_r, qs=[x_c, x_r]) + Spring(
                                k_23, pos1=x_c, pos2=x_r, qs=[x_c, x_r])
        self.Disk3 = MaterialPoint(m, x_r, qs=[x_r]) + MaterialPoint(
            m / 2, x_r, qs=[x_r]) + MaterialPoint(m3, x_r, qs=[x_r]) + Spring(
                k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[
                    x_r
                ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(
                    -2 * F_0 * cos(Omega * ivar), pos1=x_r, qs=[x_r])


        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m k l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],

            self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2],
            self.x_r: [self.x_2, 0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[
                self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[
                self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict



    
    
# OK
class MDoFForcedDisksWithSerialSprings(ComposedSystem):
    scheme_name = 'MDOF_Forced_Disks_With_Serial_Springs.PNG'
    real_name = 'three_carriages.PNG'

    def __init__(self,
                 r=Symbol('r', positive=True), #!!! Important - it's dummy variable which is to remove when the LagrangesDynamicSystem inits will be improved
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 k_ll=Symbol('k_ll', positive=True),
                 k_lr=Symbol('k_lr', positive=True),
                 k_12l=Symbol('k_12l', positive=True),
                 k_12r=Symbol('k_12r', positive=True),
                 k_23l=Symbol('k_23l', positive=True),
                 k_23r=Symbol('k_23r', positive=True),
                 k_rl=Symbol('k_rl', positive=True),
                 k_rr=Symbol('k_rr', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 x_l=dynamicsymbols('x_l'),
                 x_c=dynamicsymbols('x_c'),
                 x_r=dynamicsymbols('x_r'),
                 x_1=dynamicsymbols('x_1'),
                 x_2=dynamicsymbols('x_2'),
                 x_3=dynamicsymbols('x_3'),
                 qs=dynamicsymbols('x_l x_c x_r'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.R = R
        self.k_ll = k_ll
        self.k_lr = k_lr
        self.k_12l = k_12l
        self.k_12r = k_12r
        self.k_23l = k_23l
        self.k_23r = k_23r
        self.k_rl = k_rl
        self.k_rr = k_rr
        self.x_l = x_l
        self.x_c = x_c
        self.x_r = x_r
        self.x_1 = x_1
        self.x_2 = x_2
        self.x_3 = x_3
        self.Omega = Omega

        self.Disk1 =  MaterialPoint(m, x_l, qs = [x_l]) + MaterialPoint(m/2, x_l, qs = [x_l]) + MaterialPoint(m1, x_l, qs = [x_l]) + Spring((k_ll*k_lr)/(k_ll+k_lr), pos1 = x_l, qs = [x_l]) + Force(-2*F_0 * cos(Omega * ivar), pos1 = x_l, qs = [x_l])
        self.Disk2 =  MaterialPoint(m, x_c, qs = [x_c]) + MaterialPoint(m/2, x_c, qs = [x_c]) + MaterialPoint(m2, x_c, qs = [x_c]) + Spring((k_12l*k_12r)/(k_12l+k_12r), pos1 = x_l, pos2 = x_c, qs = [x_l, x_c]) + Spring((k_23l*k_23r)/(k_23l+k_23r), pos1 = x_c, pos2 = x_r, qs = [x_c, x_r])
        self.Disk3 =  MaterialPoint(m, x_r, qs = [x_r]) + MaterialPoint(m/2, x_r, qs = [x_r]) + MaterialPoint(m3, x_r, qs = [x_r]) + Spring((k_rl*k_rr)/(k_rl+k_rr), pos1 = x_r, qs = [x_r]) + Force(-F_0 * cos(Omega * ivar), pos1 = x_r, qs = [x_r])


        
        
        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            
            self.k_ll: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_lr: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_12l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_12r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_23l: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_23r: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_rl: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],
            self.k_rr: [1 * k0, 2 * k0, S.Half * k0, 4 * k0, S.Half**2 * k0],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2, 0],
            self.x_r: [self.x_2, 0],
            }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict
    
    

    
class MDoFTriplePendulum(ComposedSystem):
    scheme_name = 'MDOFTriplePendulum.PNG'
    real_name = 'TriplePendulum_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 l_1=Symbol('l_1', positive=True),
                 l_2=Symbol('l_2', positive=True),
                 l_3=Symbol('l_3', positive=True),
                 g=Symbol('g', positive=True),
                 phi_1=dynamicsymbols('varphi_1'),
                 phi_2=dynamicsymbols('varphi_2'),
                 phi_3=dynamicsymbols('varphi_3'),
                 phi_u=dynamicsymbols('varphi_u'),
                 phi_l=dynamicsymbols('varphi_l'),
                 qs=dynamicsymbols('varphi_1 varphi_2 varphi_3'),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.l_1 = l_1
        self.l_2 = l_2
        self.l_3 = l_3
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.phi_u = phi_u
        self.phi_l = phi_l
        self.g = g

        x_2 = sin(phi_1)*l_1 + sin(phi_2)*l_2
        y_2 = cos(phi_1)*l_1 + cos(phi_2)*l_2
        x_3 = x_2 + sin(phi_3)*l_3
        y_3 = y_2 + cos(phi_3)*l_3

        self.Pendulum1 = Pendulum(m1, g, l_1, angle=phi_1, qs=[phi_1])
        self.material_point_11 = MaterialPoint(m2, x_2, qs=[phi_1, phi_2])
        self.material_point_21 = MaterialPoint(m2, y_2, qs=[phi_1, phi_2])
        self.gravity_1 = GravitationalForce(m2, g, pos1=-y_2, qs=[phi_2])
        self.material_point_12 = MaterialPoint(m3, x_3, qs=[phi_1, phi_2, phi_3])
        self.material_point_22 = MaterialPoint(m3, y_3, qs=[phi_1, phi_2, phi_3])
        self.gravity_2 = GravitationalForce(m3, g, pos1=-y_3, qs=[phi_3])
        
        system = self.Pendulum1 + self.material_point_11 + self.material_point_12 + self.material_point_21 + self.material_point_22 + self.gravity_1 + self.gravity_2
        super().__init__(system(qs).linearized(),**kwargs)

    def get_default_data(self):


        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_3: [2 * l0, 4 * l0, S.Half * l0, 2 * l0, S.Half * l0],

            self.phi_1:[self.phi_u,0],
            self.phi_2:[self.phi_u,self.phi_l],
            self.phi_3:[self.phi_l],
        }

        return default_data_dict    

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi_2] == parameters_dict[self.phi_3]:

            parameters_dict[self.phi_2] = self.phi_u


        return parameters_dict

    

    
    
class SDoFTriplePendulum(ComposedSystem):
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
    
    
class MDoFTripleShaft(ComposedSystem):
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

    scheme_name = 'MDoFTripleShaft.PNG'
    real_name = 'MDoFTripleShaft_real.jpg'

    def __init__(self,
                 m = symbols('m', positive=True),
                 m_1=Symbol('m_1', positive=True),
                 m_2=Symbol('m_2', positive=True),
                 m_3=Symbol('m_3', positive=True),
                 l_1=Symbol('l_1', positive=True),
                 l_2=Symbol('l_2', positive=True),
                 l_3=Symbol('l_3', positive=True),
                 d=Symbol('d', positive=True),
                 d_1=Symbol('d_1', positive=True),
                 d_2=Symbol('d_2', positive=True),
                 d_3=Symbol('d_3', positive=True),
                 G=Symbol('G', positive=True),
                 M_1=Symbol('M_1', positive=True),
                 M_2=Symbol('M_2', positive=True),
                 M_3=Symbol('M_3', positive=True),
                 Omega=Symbol('\\Omega', positive=True),
                 delta=Symbol('\\delta', positive=True),

                 input_displacement=dynamicsymbols('theta'),
                 phi_1=dynamicsymbols('\\varphi_1'),
                 phi_2=dynamicsymbols('\\varphi_2'),
                 phi_3=dynamicsymbols('\\varphi_3'),

                 phi_l=dynamicsymbols('\\varphi_L'),
                 phi_r=dynamicsymbols('\\varphi_R'),

                 ivar=Symbol('t'),
                 qs=dynamicsymbols('\\varphi_1, \\varphi_2, \\varphi_3,'),
                 **kwargs):


        theta = input_displacement
        
        self.m_1=m_1
        self.m_2=m_2
        self.m_3=m_3
        self.G=G

        self.l_1=l_1
        self.l_2=l_2
        self.l_3=l_3

        self.d_1=d_1
        self.d_2=d_2
        self.d_3=d_3
        self.d=d

        self.M_1=M_1
        self.M_2=M_2
        self.M_3=M_3

        self.Omega=Omega
        self.delta=delta

        self.phi_1=phi_1
        self.phi_2=phi_2
        self.phi_3=phi_3
        
        self.phi_l=phi_l
        self.phi_r=phi_r

        self.input_displacement = input_displacement
        self.qs = qs

        I_0=S.Half*pi*(d/2)**4
        I_1=S.Half*m_1*(d_1/2)**2
        I_2=S.Half*m_2*(d_2/2)**2
        I_3=S.Half*m_3*(d_3/2)**2

        k_1=(G*I_0)/l_1
        k_2=(G*I_0)/l_1
        k_3=(G*I_0)/l_1

        self.disk_1 = Disk(I_1, pos1=phi_1, qs=qs) + Force(-M_1 * cos(Omega * ivar + delta), pos1 = phi_1, qs = [phi_1])
        self.spring_1 = Spring(k_1, pos1=theta, pos2=phi_1, qs=qs)  # left rod
        self.disk_2 = Disk(I_2, pos1=phi_2, qs=qs) + Force(-M_2 * cos(Omega * ivar + delta), pos1 = phi_2, qs = [phi_2])
        self.spring_2 = Spring(k_2, pos1=phi_1, pos2=phi_2, qs=qs)  # central rod
        self.disk_3 = Disk(I_3, pos1=phi_3, qs=qs) + Force(-M_3 * cos(Omega * ivar + delta), pos1 = phi_3, qs = [phi_3])
        self.spring_3 = Spring(k_3, pos1=phi_2, pos2=phi_3, qs=qs)  # right rod

        system = self.disk_1 + self.disk_2 + self.disk_3 + self.spring_1 + self.spring_2 + self.spring_3

        super().__init__(system(qs),**kwargs)

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.I: r'Moment of Inertia',
#             self.k_1: r'',
#             self.k_2: r'',
#         }
#         return self.sym_desc_dict

    
    def get_default_data(self):


        m0, M0, l0 , d0 ,M0= symbols('m_0 M_0 l_0 d_0 M_0', positive=True)
        theta0, Omega = symbols('theta_0, Omega', positive=True)

        default_data_dict = {
            
            Symbol('I_m'):[S.Half*Symbol('m')*Symbol('R')**2],
            Symbol('I_S'):[S.Half*pi*Symbol('R_S')**4],
            Symbol('k_S'):[Symbol('G')*Symbol('I_S')/Symbol('L_S')],
            
            
            self.m_1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m_2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m_3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_3: [2 * l0, 4 * l0, S.Half * l0, 2 * l0, S.Half * l0],

            self.d_1: [1 * d0, 2 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d_2: [1 * d0, 2 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d_3: [2 * d0, 4 * d0, S.Half * d0, 2 * d0, S.Half * d0],
            self.d: [2 * d0, 4 * d0, S.Half * d0, 2 * d0, S.Half * d0],

            self.phi_1:[self.phi_l,0],
            self.phi_2:[self.phi_l,self.phi_r],
            self.phi_3:[self.phi_r,0],
            
            self.delta:[0],

            self.input_displacement:[0],
            
            self.M_1:[M0],
            self.M_2:[M0],
            self.M_3:[M0],
            

        }

        return default_data_dict


    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi_1] != self.phi_l or parameters_dict[self.phi_2] != self.phi_l:

            parameters_dict[self.phi_1] = self.phi_l

        if parameters_dict[self.phi_2] != self.phi_r or parameters_dict[self.phi_3] != self.phi_r:

            parameters_dict[self.phi_3] = self.phi_r
            

        return parameters_dict
    
    

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

class SDOFWinchSystem(ComposedSystem):


    scheme_name = 'Winch_System.png'
    real_name = 'winch_mechanism_real.PNG'

    def __init__(self,
                 I_k=Symbol('I_k', positive=True),
                 I_1=Symbol('I_1', positive=True),
                 i_1=Symbol('i_1', positive=True),
                 I_2=Symbol('I_2', positive=True),
                 i_2=Symbol('i_2', positive=True),
                 I_3=Symbol('I_3', positive=True),
                 I_4=Symbol('I_4', positive=True),
                 I_b=Symbol('I_b', positive=True),
                 ivar=Symbol('t'),
                 phi=dynamicsymbols('\\varphi'),
                 D=Symbol('D', positive=True),
                 v=Symbol('v', positive=True),
                 M_s=Symbol('M_s', positive=True),
                 M_T=Symbol('M_T', positive=True),
                 G=Symbol('G', positive=True),
                 g=Symbol('g',positive=True),
                 alpha=Symbol('\\alpha'),
                 **kwargs):
        
        pos1=phi
        qs=[phi]
        phi_1=phi
        phi_2=phi_1/i_1
        phi_3=phi_2/i_2

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
        self.v = v
        self.M_T = M_T
        self.M_s = M_s
        self.phi = phi
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        

        

        self.disk_k = Disk(I_k, phi_1, qs=qs)
        self.disk_1 = Disk(I_1, phi_1, qs=qs)
        self.disk_2 = Disk(I_2, phi_2, qs=qs)
        self.disk_3 = Disk(I_3, phi_2, qs=qs)
        self.disk_4 = Disk(I_4, phi_3, qs=qs)
        self.disk_B = Disk(I_b, phi_3, qs=qs)
        self.load = MaterialPoint(G/g, pos1=phi_3* D/2, qs=qs)
        self.gravity = GravitationalForce(G/g, g, pos1=phi_3* D/2 * cos(alpha) , qs=qs)
        self.force = Force(-M_T,pos1=phi_3,qs=qs)
        self.engine = Force(M_s,pos1=phi_3,qs=qs)

        system = self.disk_k + self.disk_1 + self.disk_2 + self.disk_3 + self.disk_4 + self.disk_B + self.load + self.gravity + self.force + self.engine

        super().__init__(system,**kwargs)
        
        

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
        