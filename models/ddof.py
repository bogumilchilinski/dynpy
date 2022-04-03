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
    detail_scheme_name = 'damped_car_new.PNG'
    detail_real_name = 'car_real.jpg'
    @classmethod
    def _scheme(cls):
        if 'systems.py' in __file__: 
            path = __file__.replace('systems.py', 'images/') + cls.scheme_name
        if 'sdof.py' in __file__: 
            path = __file__.replace('sdof.py', 'images/') + cls.scheme_name
        if 'ddof.py' in __file__: 
            path = __file__.replace('ddof.py', 'images/') + cls.scheme_name
        if 'mdof.py' in __file__: 
            path = __file__.replace('mdof.py', 'images/') + cls.scheme_name
        return path

    @classmethod
    def _real_example(cls):
        if 'systems.py' in __file__: 
            path = __file__.replace('systems.py', 'images/') + cls.real_name
        if 'sdof.py' in __file__: 
            path = __file__.replace('sdof.py', 'images/') + cls.real_name
        if 'ddof.py' in __file__: 
            path = __file__.replace('ddof.py', 'images/') + cls.real_name
        if 'mdof.py' in __file__: 
            path = __file__.replace('mdof.py', 'images/') + cls.real_name

        return path
    
    @classmethod
    def _detail_real(cls):
        if 'systems.py' in __file__: 
            path = __file__.replace('systems.py', 'images/') + cls.detail_real_name
        if 'sdof.py' in __file__: 
            path = __file__.replace('sdof.py', 'images/') + cls.detail_real_name
        if 'ddof.py' in __file__: 
            path = __file__.replace('ddof.py', 'images/') + cls.detail_real_name
        if 'mdof.py' in __file__: 
            path = __file__.replace('mdof.py', 'images/') + cls.detail_real_name

        return path
    
    @classmethod
    def _detail_scheme(cls):
        if 'systems.py' in __file__: 
            path = __file__.replace('systems.py', 'images/') + cls.detail_scheme_name
        if 'sdof.py' in __file__: 
            path = __file__.replace('sdof.py', 'images/') + cls.detail_scheme_name
        if 'ddof.py' in __file__: 
            path = __file__.replace('ddof.py', 'images/') + cls.detail_scheme_name
        if 'mdof.py' in __file__: 
            path = __file__.replace('mdof.py', 'images/') + cls.detail_scheme_name

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


class BeamBridgeTMD(ComposedSystem):

    scheme_name = 'bridge_tmd.png'
    real_name = 'beam_bridge_real.PNG'
    detail_real_name = 'real_bidge_tmd.jpg'

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
        self.force = Force(F_0 * sin(Omega * ivar), pos1=z)
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
            self.m: [10*m0, 20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0, 70 * m0, 80 * m0, 90 * m0],
            self.k_beam: [1 * 48 * E * I / l**3,
                2 * 48 * E * I / l**3, 3 * 48 * E * I / l**3,
                4 * 48 * E * I / l**3, 5 * 48 * E * I / l**3,
                6 * 48 * E * I / l**3,
                7 * 48 * E * I / l**3,
                8 * 48 * E * I / l**3,
                9 * 48 * E * I / l**3
            ],
            self.m_TMD: [1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0],
            self.k_TMD: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
        }

        return default_data_dict


class BeamBridgeDampedTMD(ComposedSystem):

    scheme_name = 'bridge_tmd_dmp.png'
    real_name = 'beam_bridge_real.PNG'
    detail_real_name = 'real_bidge_tmd.jpg'
    
    def __init__(self,
                 m=Symbol('m', positive=True),
                 m_TMD=Symbol('m_TMD', positive=True),
                 k_beam=Symbol('k_beam', positive=True),
                 k_TMD=Symbol('k_TMD', positive=True),
                 ivar=Symbol('t'),
                 g=Symbol('g', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 F_0=Symbol('F_0', positive=True),
                 c_TMD=Symbol('c_TMD', positive=True),
                 c=Symbol('c', positive=True),
                 z_TMD=BeamBridgeTMD().z_TMD,
                 z=BeamBridgeTMD().z,
                 **kwargs):
        qs = z, z_TMD
        
        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F_0 = F_0
        self.m_TMD = m_TMD
        self.k_TMD = k_TMD
        self.z_TMD = z_TMD
        self.z = z
        self.c = c
        self.c_TMD = c_TMD

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.gravity_TMD = GravitationalForce(self.m_TMD, self.g, z_TMD)
        self.force = Force(F_0 * sin(Omega * ivar), pos1=z)
        self.TMD = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.spring_TMD = Spring(k_TMD, z, z_TMD, qs=[z, z_TMD])
        
        self.damper_TMD = Damper(c=c_TMD, pos1=z - z_TMD, qs=qs)  # left damper
        self.damper = Damper(c=c, pos1=z, qs=qs)
        composed_system = self.mass + self.spring + self.gravity_force + self.force + self.TMD + self.spring_TMD + self.gravity_TMD + self.damper + self.damper_TMD

        super().__init__(composed_system,**kwargs)

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
                self.m: [10 * m0, 20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0, 70*m0, 80*m0, 90*m0],
                self.c: [lamb * self.k_beam],
                self.k_beam: [10 * 48 * E * I / l**3,
                    20 * 48 * E * I / l**3, 30 * 48 * E * I / l**3,
                    40 * 48 * E * I / l**3, 50 * 48 * E * I / l**3,
                    60 * 48 * E * I / l**3, 70 * 48 * E * I / l**3, 80 * 48 * E * I / l**3, 90 * 48 * E * I / l**3
                ],
                self.c_TMD: [lamb * self.k_TMD],
                self.m_TMD: [m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7*m0, 8*m0, 9*m0],
                self.k_TMD: [
                    1 * 48 * E * I / l**3, 2 * 48 * E * I / l**3,
                    3 * 48 * E * I / l**3, 4 * 48 * E * I / l**3,
                    5 * 48 * E * I / l**3, 6 * 48 * E * I / l**3, 7 * 48 * E * I / l**3, 8 * 48 * E * I / l**3, 9 * 48 * E * I / l**3
                ],
            }
        else:
            numerized_dict = {m0:100 , lamb:2 , E:200000, I:0.48, l:3}
            default_data_dict = {
                self.m: [20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0],
                self.c: [lamb * self.k_beam],
                self.k_beam: [
                    20 * 48 * E * I / l**3, 30 * 48 * E * I / l**3,
                    40 * 48 * E * I / l**3, 50 * 48 * E * I / l**3,
                    60 * 48 * E * I / l**3],
                self.c_TMD: [lamb * self.k_TMD],
                self.m_TMD: [2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0],
                self.k_TMD: [
                    1 * 48 * E * I / l**3, 3 * 48 * E * I / l**3,
                    3 * 48 * E * I / l**3, 5 * 48 * E * I / l**3,
                    5 * 48 * E * I / l**3]
            }

            values_of_dict = flatten([default_data_dict[self.m], default_data_dict[self.c],
                                      default_data_dict[self.k_beam], default_data_dict[self.c_TMD],
                                      default_data_dict[self.m_TMD], default_data_dict[self.k_TMD]
                                     ])

            for i in range(len(values_of_dict)):
                    default_data_dict = display(values_of_dict[i].subs(numerized_dict))

        return default_data_dict


    def get_real_data(self):
        
        return self.get_default_data(True)

class SimplifySuspension(ComposedSystem):
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


class VehicleSuspension(ComposedSystem):
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


# class DampedVehicleSuspension(ComposedSystem):

#     scheme_name = 'damped_car_new.PNG'
#     real_name = 'car_real.jpg'

#     def __init__(self,
#                  non_damped_system,
#                  c_l=Symbol('c_l', positive=True),
#                  c_r=Symbol('c_r', positive=True),
#                  l_cl=Symbol('l_{cl}', positive=True),
#                  l_cr=Symbol('l_{cr}', positive=True),
#                  k_1=VehicleSuspension().k_1,
#                  k_2=VehicleSuspension().k_2,
#                  l_l=VehicleSuspension().l_l,
#                  l_r=VehicleSuspension().l_r,
#                  qs=dynamicsymbols('z, \\varphi'),
#                  **kwargs):

#         z, phi = qs

#         self.k_1 = k_1
#         self.k_2 = k_2
#         self.c_l = c_l
#         self.c_r = c_r
#         self.l_cl = l_cl
#         self.l_cr = l_cr
#         self.l_l = l_l
#         self.l_r = l_r
#         self.nds = non_damped_system
#         self.damper_l = Damper(c=c_l, pos1=z + phi * l_cl,
#                                qs=qs)  # left damper
#         self.damper_r = Damper(c=c_r, pos1=z - phi * l_cr,
#                                qs=qs)  # right damper
#         system = self.nds + self.damper_l + self.damper_r

#         super().__init__(system,**kwargs)

#     def get_default_data(self):

#         c0, k_0, l_l0, omega, F_0 = symbols('c_0 k_0 l_0 omega F_0',
#                                             positive=True)

#         default_data_dict = {
#             self.c_r: [self.c_l],
#             self.k_2: [self.k_1],
#             self.l_r: [self.l_l],
#             self.l_cr: [self.l_l],
#             self.l_cl: [self.l_l],
#             self.c_l: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
#             self.k_1: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
#             self.l_l: [2 * l_l0, 3 * l_l0, 4 * l_l0, 5 * l_l0, 6 * l_l0],
#             self.nds.F_engine: [
#                 2 * F_0 * sin(omega * self.nds.ivar),
#                 3 * F_0 * sin(omega * self.nds.ivar),
#                 4 * F_0 * sin(omega * self.nds.ivar),
#                 5 * F_0 * sin(omega * self.nds.ivar),
#                 6 * F_0 * sin(omega * self.nds.ivar)
#             ]
#         }

#         return default_data_dict
class DampedVehicleSuspension(ComposedSystem):

    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_{rod}', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t', positive=True),
                 qs=dynamicsymbols('z, \\varphi'),
                 c_l=Symbol('c_l', positive=True),
                 c_r=Symbol('c_r', positive=True),
                 l_cl=Symbol('l_{cl}', positive=True),
                 l_cr=Symbol('l_{cr}', positive=True),

                 **kwargs):

        z, phi = qs

        self.m = m  # mass of a rod

        self.l_rod = l_rod  # length of a rod
        self.k_r = k_r  # left spring
        self.k_l = k_l  # right spring
        self.I = I  # moment of inertia of a rod
        self.F_engine = F_engine


        self.c_l = c_l
        self.c_r = c_r
        self.l_cl = l_cl
        self.l_cr = l_cr
        self.l_l = l_l
        self.l_r = l_r
        
        self.damper_l = Damper(c=c_l, pos1=z + phi * l_l,
                               qs=qs)  # left damper
        self.damper_r = Damper(c=c_r, pos1=z - phi * l_r,
                               qs=qs)  # right damper
        
        self.body = RigidBody2D(m, I, pos_lin=z, pos_rot=phi, qs=qs)  # rod
        self.spring_1 = Spring(k_l, pos1=z + phi * l_l, qs=qs)  # left spring
        self.spring_2 = Spring(k_r, pos1=z - phi * l_r, qs=qs)  # right spring
        self.force = Force(F_engine, pos1=z - l_r * phi, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force + self.damper_l + self.damper_r

        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0',
                                            positive=True)

        default_data_dict = {
            self.I: [self.m*self.l_rod**2],
            self.l_rod:[l0,2*l0,3*l0,4*l0,5*l0,6*l0,7*l0,8*l0,9*l0],
            self.m: [m0,2*m0,3*m0,4*m0,5*m0,6*m0,7*m0,8*m0,9*m0],
            self.c_r: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.k_l: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            self.k_r: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            
            
            self.c_l: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.l_r: [0.2*self.l_rod,0.25 * self.l_rod, 0.3 * self.l_rod, 0.4 * self.l_rod, 0.5 * self.l_rod,0.35 * self.l_rod, 0.45 * self.l_rod],
            self.l_l: [0.2*self.l_rod,0.25 * self.l_rod, 0.3 * self.l_rod, 0.4 * self.l_rod, 0.5 * self.l_rod,0.35 * self.l_rod, 0.45 * self.l_rod],

            self.F_engine: [
                2 * F_0 * sin(omega * self.ivar),
                3 * F_0 * sin(omega * self.ivar),
                4 * F_0 * sin(omega * self.ivar),
                5 * F_0 * sin(omega * self.ivar),
                6 * F_0 * sin(omega * self.ivar)
            ]
        }

        return default_data_dict
    
class UndampedVehicleSuspension(ComposedSystem):
    scheme_name = 'car_undamped.png'
    real_name = 'car_real.jpg'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 I=Symbol('I', positive=True),
                 l_rod=Symbol('l_{rod}', positive=True),
                 l_l=Symbol('l_l', positive=True),
                 l_r=Symbol('l_r', positive=True),
                 k_r=Symbol('k_r', positive=True),
                 k_l=Symbol('k_l', positive=True),
                 F_engine=Symbol('F_{engine}', positive=True),
                 ivar=Symbol('t', positive=True),
                 qs=dynamicsymbols('z, \\varphi'),
                 c_l=Symbol('c_l', positive=True),
                 c_r=Symbol('c_r', positive=True),
                 l_cl=Symbol('l_{cl}', positive=True),
                 l_cr=Symbol('l_{cr}', positive=True),

                 **kwargs):

        z, phi = qs

        self.m = m  # mass of a rod

        self.l_rod = l_rod  # length of a rod
        self.k_r = k_r  # left spring
        self.k_l = k_l  # right spring
        self.I = I  # moment of inertia of a rod
        self.F_engine = F_engine


        self.c_l = c_l
        self.c_r = c_r
        self.l_cl = l_cl
        self.l_cr = l_cr
        self.l_l = l_l
        self.l_r = l_r
        
        self.damper_l = Damper(c=c_l, pos1=z + phi * l_l,
                               qs=qs)  # left damper
        self.damper_r = Damper(c=c_r, pos1=z - phi * l_r,
                               qs=qs)  # right damper
        
        self.body = RigidBody2D(m, I, pos_lin=z, pos_rot=phi, qs=qs)  # rod
        self.spring_1 = Spring(k_l, pos1=z + phi * l_l, qs=qs)  # left spring
        self.spring_2 = Spring(k_r, pos1=z - phi * l_r, qs=qs)  # right spring
        self.force = Force(F_engine, pos1=z - l_r * phi, qs=qs)
        system = self.body + self.spring_1 + self.spring_2 + self.force

        super().__init__(system,**kwargs)

    def get_default_data(self):

        m0, l0, c0, k_0, l_l0, omega, F_0 = symbols('m_0 l_0 c_0 k_0 l_l0 Omega F_0',
                                            positive=True)

        default_data_dict = {
            self.I: [self.m*self.l_rod**2],
            self.l_rod:[l0,2*l0,3*l0,4*l0,5*l0,6*l0,7*l0,8*l0,9*l0],
            self.m: [m0,2*m0,3*m0,4*m0,5*m0,6*m0,7*m0,8*m0,9*m0],
            self.c_r: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.k_l: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            self.k_r: [2 * k_0, 3 * k_0, 4 * k_0, 5 * k_0, 6 * k_0],
            
            
            self.c_l: [2 * c0, 3 * c0, 4 * c0, 5 * c0, 6 * c0],
            self.l_r: [0.2*self.l_rod,0.25 * self.l_rod, 0.3 * self.l_rod, 0.4 * self.l_rod, 0.5 * self.l_rod,0.35 * self.l_rod, 0.45 * self.l_rod],
            self.l_l: [0.2*self.l_rod,0.25 * self.l_rod, 0.3 * self.l_rod, 0.4 * self.l_rod, 0.5 * self.l_rod,0.35 * self.l_rod, 0.45 * self.l_rod],

            self.F_engine: [
                2 * F_0 * sin(omega * self.ivar),
                3 * F_0 * sin(omega * self.ivar),
                4 * F_0 * sin(omega * self.ivar),
                5 * F_0 * sin(omega * self.ivar),
                6 * F_0 * sin(omega * self.ivar)
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


class Shaft(ComposedSystem):
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
        self.spring_1 = Spring(k_1, phi_1, phi_2, qs=qs)  # left spring
        self.disc_2 = Disk(I, pos1=phi_2, qs=qs)
        self.spring_2 = Spring(k_2, pos1=phi_2, pos2=theta,
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
            self.I: [S.Half*m0*l**2,S.One*m0*l**2,(S.Half**2)*m0*l**2,2*m0*l**2,3*m0*l**2,4*m0*l**2,5*m0*l**2,6*m0*l**2,7*m0*l**2,8*m0*l**2,9*m0*l**2],
            
            

            self.k_1: [1 * G * l**4 /(10* l), 2 *  G * l**4 /(10* l), S.Half *  G * l**4 /(10* l), 4 *  G * l**4 /(10* l), (S.Half**2) *  G * l**4 /(10* l), 3 *  G * l**4 /(10* l), 5 *  G * l**4 /(10* l), 6 *  G * l**4 /(10* l), 7 *  G * l**4 /(10* l), 8 *  G * l**4 /(10* l), 9 *  G * l**4 /(10* l)],
            self.k_2: [1 * G * l**4 /(10* l), 2 *  G * l**4 /(10* l), S.Half *  G * l**4 /(10* l), 4 *  G * l**4 /(10* l), (S.Half**2) *  G * l**4 /(10* l), 3 *  G * l**4 /(10* l), 5 *  G * l**4 /(10* l), 6 *  G * l**4 /(10* l), 7 *  G * l**4 /(10* l), 8 *  G * l**4 /(10* l), 9 *  G * l**4 /(10* l)],

            l:[1 * l0, 2 * l0, S.Half * l0, 4 * l0, (S.Half**2) * l0, 3 * l0, (S.Half**2) * l0, 5 * l0, (S.Half**2) * l0, 6 * l0, (S.Half**2) * l0, 7 * l0, (S.Half**2) * l0, 8 * l0, (S.Half**2) * l0, 9 * l0, (S.Half**2) * l0],
            
#             self.phi_1:[self.phi, 0],
#             self.phi_2:[self.phi],
            self.input_displacement:[theta0* cos(Omega * self.ivar) ],
        }

        return default_data_dict
    

class DampedShaft(ComposedSystem):


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



# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie


class CoupledPendulum(ComposedSystem):
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
    phi_1,phi_2=dynamicsymbols('\\varphi_1, \\varphi_2')
    _hint = [phi_1,phi_2]

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 k=Symbol('k', positive=True),
                 qs=[phi_1,phi_2],
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
    
class LinearizedCoupledPendulum(ComposedSystem):
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

    scheme_name = 'tmd_engine_vertical_spring_gravity.png'
    real_name = 'tmd_engine_real.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 k_TMD=Symbol('k_{TMD}', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 m_TMD=Symbol('m_{TMD}', positive=True),
                 e=Symbol('e', positive=True),
                 z=dynamicsymbols('z'),
                 z_TMD=dynamicsymbols('z_{TMD}'),
                 Omega=Symbol('\Omega',positive=True),
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
        self.Omega=Omega
        self.t=ivar
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
    def get_default_data(self):

        m0, k0, e0 = symbols('m_0 k_0 e_0', positive=True)

        default_data_dict = {
            self.M: [10 * m0, 20 * m0, 30 * m0, 40 * m0, 50 * m0, 60 * m0, 70 * m0, 80 * m0, 90 * m0],
            self.k_m: [k0, 2*k0, 3*k0, 4*k0, 5*k0, 6*k0, 7*k0, 8*k0, 9*k0],
            self.m_TMD: [m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0],
            self.k_TMD: [k0, 2*k0, 3*k0, 4*k0, 5*k0, 6*k0, 7*k0, 8*k0, 9*k0],
            self.phi:[self.Omega*self.t],
            self.m_e: [0.1*m0, 0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0, 0.7 * m0, 0.8 * m0, 0.9 * m0],
            self.e:[e0, 2 * e0, 3 * e0, 4 * e0, 5 * e0, 6 * e0, 7*e0, 8*e0, 9* e0]
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


class TwoDisksWithThreeSprings(ComposedSystem):
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

class TwoNonLinearTrolleys(ComposedSystem):
    scheme_name = 'ddof_nonlin_trolleys.PNG'
    real_name = 'dwa_wozki_XD.PNG'

    def __init__(self,
                 g=Symbol('g', positive=True),
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
            self.k1: [S.Half * k0, S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0],
            self.k2: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
            self.k3: [S.Half * k0, 1 * k0, 3 * S.Half * k0, 2 * k0, 5 * S.Half * k0],
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


class TwoNonLinearDisks(ComposedSystem):
    scheme_name = 'sdof_nonlin_disc.png'
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

class VeeEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'vee_engine_perpendicular_springs.png'
    real_name = '440_magnum_v8.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l',positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1x = MaterialPoint(M, pos1=x, qs=[x])
        self.MaterialPoint_1z = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2x = MaterialPoint(m_e,
                                             pos1=x + e * sin(phi),
                                             qs=[x])
        self.MaterialPoint_2z = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.SpringVer = Spring(2 * k_m, pos1=z, qs=[z])
        self.SpringHor = Spring(2 * k_m, pos1=x, qs=[x])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1x + self.MaterialPoint_2x + self.MaterialPoint_1z + self.MaterialPoint_2z
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict



    
class NonLinearVeeEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'nonlin_vee_engine_perpendicular_springs.png'
    real_name = '440_magnum_v8.jpg'

    def __init__(self,
                 M=Symbol('M', positive=True),
                 k_m=Symbol('k_m', positive=True),
                 m_e=Symbol('m_e', positive=True),
                 e=Symbol('e', positive=True),
                 l=Symbol('l',positive=True),
                 x=dynamicsymbols('x'),
                 z=dynamicsymbols('z'),
                 phi=dynamicsymbols('phi'),
                 ivar=Symbol('t', positive=True),
                 **kwargs):

        self.M = M
        self.k_m = k_m
        self.m_e = m_e
        self.e = e

        self.MaterialPoint_1x = MaterialPoint(M, pos1=x, qs=[x])
        self.MaterialPoint_1z = MaterialPoint(M, pos1=z, qs=[z])
        self.MaterialPoint_2x = MaterialPoint(m_e,
                                             pos1=x + e * sin(phi),
                                             qs=[x])
        self.MaterialPoint_2z = MaterialPoint(m_e,
                                             pos1=z + e * cos(phi),
                                             qs=[z])
        self.SpringVer = Spring(2 * k_m, pos1=(z**2+x**2)**(1/2)-x, qs=[z])
        self.SpringHor = Spring(2 * k_m, pos1=(z**2+x**2)**(1/2)-z, qs=[x])
        system = self.SpringVer + self.SpringHor + self.MaterialPoint_1x + self.MaterialPoint_2x + self.MaterialPoint_1z + self.MaterialPoint_2z
        super().__init__(system,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }

        return self.sym_desc_dict
