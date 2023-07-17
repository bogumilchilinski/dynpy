from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, pi, latex, dsolve,
                   solve, fraction, factorial,Subs, Number, oo, Abs, N)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ..dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from .elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin
from ..continuous import ContinuousSystem, PlaneStressProblem



import base64
import random
import IPython as IP
import numpy as np

import inspect

from .mechanics.principles import ComposedSystem, NonlinearComposedSystem,  base_frame, base_origin,cached_property, lru_cache

class ComposedSystem(HarmonicOscillator):
    """Base class for all systems

    """
    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'
    _default_args = ()
    _default_folder_path = "./dynpy/models/images/"


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

    @lru_cache   
    def linearized(self, x0=None, op_point=False, hint=[], label=None):

        return type(self).from_system(super().linearized(x0=x0,op_point=op_point,hint=hint,label=label))



# class TMD(ComposedSystem):
#     scheme_name = 'mdof_tmd.png'
#     real_name = 'mdof_tmd_real.png'
#     """
#     Model of an exemplary Tuned Mass Damper (TMD) simulated as Double Degree of Freedom of coupled trolleys.

#         Arguments:
#         =========
#             m = Mass
#                 -Mass of system on spring

#             me = Mass
#                 -Mass of tuned mass damper

#             k = spring coefficient
#                 -value of spring coefficient

#             ke = spring coefficient
#                 -value of spring coefficient that tuned mass damper is mounted

#             F = Force
#                 -Trolley's exciting force

#             ivar = symbol object
#                 -Independant time variable

#             qs = dynamicsymbol object
#                 -Generalized coordinates

#         Example
#         =======
#         >>> t = symbols('t')
#         >>> m, me, k, ke, F = symbols('m, m_e, k, k_e, F')

        
        
#     """

#     def __init__(self,
#                  m=Symbol('m', positive=True),
#                  me=Symbol('m_e', positive=True),
#                  k=Symbol('k', positive=True),
#                  ke=Symbol('k_e', positive=True),
#                  F=Symbol('F', positive=True),
#                  xe=dynamicsymbols('x_e'),
#                  xb=dynamicsymbols('x_b'),
#                  angle=dynamicsymbols('Omega'),
#                  ivar=Symbol('t', positive=True),
#                  **kwargs):

#         self.m = m
#         self.me = me
#         self.k = k
#         self.ke = ke
#         self.F = F

#         self.MaterialPoint_1 = MaterialPoint(m, pos1=xe, qs=[xe])
#         self.MaterialPoint_2 = MaterialPoint(me, pos1=xb, qs=[xb])
#         self.Spring_1 = Spring(k, pos1=xe, qs=[xe])
#         self.Spring_2 = Spring(ke, pos1=xe, pos2=xb, qs=[xe, xb])
#         self.Force = Force(F * sin(angle * ivar), pos1=xe, qs=[xe])

#         system = self.Spring_1 + self.Spring_2 + \
#             self.MaterialPoint_1 + self.MaterialPoint_2 + self.Force
#         super().__init__(system,**kwargs)

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m: r'Mass of main object',
#             self.me: r'Mass of the TMD',
#             self.k: r'Stiffness coefficient',
#             self.ke: r'Stiffness coefficient',
#             self.F: r'Force',
#         }
#         return self.sym_desc_dict

#DO ZROBIENIA
class DDoFWinch(ComposedSystem):
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



#DO ZROBIENIA
class DampedElasticPendulum(ComposedSystem):
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

#DO SPRAWDZENIA
class ForcedTrolleysWithSprings(ComposedSystem):
    scheme_name = 'mdof_three_trolleys.PNG'
    real_name = 'three_carriages.PNG'

    def __init__(self,
                 dum=Symbol('dumy'),
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
                 F=Symbol('F_0', positive=True),
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
        super().__init__(system(x_l,x_c),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0,],
            self.k_l: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_cl: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_12: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_23: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_r: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],
            self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,],

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
        print(parameters_dict)
        if parameters_dict[self.x_l] != self.x_1 or parameters_dict[self.x_c] != self.x_1:

            parameters_dict[self.x_l] = self.x_1

        if parameters_dict[self.x_c] != self.x_2 or parameters_dict[self.x_r] != self.x_2:

            parameters_dict[self.x_r] = self.x_2

        return parameters_dict
    
#DO ZROBIENIA //wydane Maja
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
    
#DO ZROBIENIA //wydane Michał Szewczyk
class MDoFForcedDisksWithParallelSprings(ComposedSystem):

    _default_subs_method='direct'
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
                ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(2*F_0 * cos(Omega * ivar), pos1=x_l, qs=[x_l])

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
                ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(F_0 * cos(Omega * ivar), pos1=x_r, qs=[x_r])


        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

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
    
#   def get_default_data(self):

#         m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

#         default_data_dict = {
#             self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
#             self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

#             self.x_l: [self.x_1, 0],
#             self.x_c: [self.x_1, self.x_2, 0],
#             self.x_r: [self.x_2, 0],
#         }

#         return default_data_dict

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
    
#DO SPRAWDZENIA
class DoubleTrolleyDifferentWheels(ComposedSystem):
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

#DO SPRAWDZENIA
class DampedTrolleysWithSprings(ComposedSystem):
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

#DO ZROBIENIA
class ThreePendulumsWithSprings(ComposedSystem):
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

    
#DO SPRAWDZENIA
class DoubleTrolleyWithNonlinearSprings(ComposedSystem):
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
    
#DO ZROBIENIA
class MDoFForcedSimpleDisksWithParallelSprings(ComposedSystem):

    _default_subs_method='direct'
    scheme_name = 'three_simple_disks.png'
    real_name = 'three_rollers_real.png'

    def __init__(self,
                 dum=Symbol('dum'),
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
                 F=Symbol('F', positive=True),
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
        self.R=R
        self.Disk1 = MaterialPoint(m, x_l, qs=[x_l]) + MaterialPoint(m/2*R**2, x_l/R, qs=[x_l])  + Spring(k_l, pos1=x_l, qs=[x_l]) + Spring(k_l, pos1=x_l, qs=[x_l]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(2*F * cos(Omega * ivar), pos1=x_l, qs=[x_l])

        self.Disk2 = MaterialPoint(m, x_c, qs=[x_c]) + MaterialPoint(m/2*R**2, x_c/R, qs=[x_c]) + Spring(k_12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(k_c12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(k_12, pos1=x_l, pos2=x_c, qs=[x_l, x_c]) + Spring(k_23, pos1=x_c, pos2=x_r, qs=[x_c, x_r]) + Spring(k_c23, pos1=x_c, pos2=x_r, qs=[x_c, x_r]) + Spring(k_23, pos1=x_c, pos2=x_r, qs=[x_c, x_r])

        self.Disk3 = MaterialPoint(m, x_r, qs=[x_r]) + MaterialPoint(m/2*R**2, x_r/R, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_r, pos1=x_r, qs=[x_r]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(F * cos(Omega * ivar), pos1=x_r, qs=[x_r])


        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0, = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {

            self.m: [1 * m0, 2 * m0, 4 * m0, 3 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0, 9 * m0],

            self.k_l: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_cl: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_12: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_c12: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_23: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_c23: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_r: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],
            self.k_cr: [1 * k0, 2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0],

            self.x_l: [self.x_1, 0],
            self.x_c: [self.x_1, self.x_2],
            self.x_r: [self.x_2, 0],
        }

        return default_data_dict
    
#   def get_default_data(self):

#         m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

#         default_data_dict = {
#             self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
#             self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

#             self.x_l: [self.x_1, 0],
#             self.x_c: [self.x_1, self.x_2, 0],
#             self.x_r: [self.x_2, 0],
#         }

#         return default_data_dict

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

#DO ZROBIENIA
class MDoFForcedSimpleDisksWithSerialSprings(ComposedSystem):
    scheme_name = 'three_simple_disks_serial.png'
    real_name = 'three_carriages.PNG'

    def __init__(self,
                 r=Symbol('r', positive=True), #!!! Important - it's dummy variable which is to remove when the LagrangesDynamicSystem inits will be improved
                 R=Symbol('R', positive=True),
                 m=Symbol('m', positive=True),
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

        self.Disk1 =  MaterialPoint(m, x_l, qs=[x_l]) + MaterialPoint(m/2*R**2, x_l/R, qs=[x_l]) + Spring((k_ll*k_lr)/(k_ll+k_lr), pos1 = x_l, qs = [x_l]) + Force(-2*F_0 * cos(Omega * ivar), pos1 = x_l, qs = [x_l])
        self.Disk2 =  MaterialPoint(m, x_c, qs=[x_c]) + MaterialPoint(m/2*R**2, x_c/R, qs=[x_c]) + Spring((k_12l*k_12r)/(k_12l+k_12r), pos1 = x_l, pos2 = x_c, qs = [x_l, x_c]) + Spring((k_23l*k_23r)/(k_23l+k_23r), pos1 = x_c, pos2 = x_r, qs = [x_c, x_r])
        self.Disk3 =  MaterialPoint(m, x_r, qs=[x_r]) + MaterialPoint(m/2*R**2, x_r/R, qs=[x_r]) + Spring((k_rl*k_rr)/(k_rl+k_rr), pos1 = x_r, qs = [x_r]) + Force(-F_0 * cos(Omega * ivar), pos1 = x_r, qs = [x_r])


        
        
        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

        default_data_dict = {

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
    
#DO SPRAWDZENIA
class ForcedDisksWithParallelSprings(ComposedSystem):

    _default_subs_method='direct'
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
                ]) + Spring(k_cl, pos1=x_l, qs=[x_l]) + Force(2*F_0 * cos(Omega * ivar), pos1=x_l, qs=[x_l])

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
                ]) + Spring(k_cr, pos1=x_r, qs=[x_r]) + Force(F_0 * cos(Omega * ivar), pos1=x_r, qs=[x_r])


        system = self.Disk1 + self.Disk2 + self.Disk3
        super().__init__(system(qs),**kwargs)

    def get_default_data(self):

        m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

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

            self.x_l: [self.x_1],
            self.x_c: [0],
            self.x_r: [self.x_2],
        }

        return default_data_dict
    
#   def get_default_data(self):

#         m0, k0, l0 = symbols('m_0 k_0 l_0', positive=True)

#         default_data_dict = {
#             self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
#             self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
#             self.k_l: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cl: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c12: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_c23: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_r: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],
#             self.k_cr: [1 * k0, 2 * k0, S.Half * k0, 2 * k0, S.Half * k0],

#             self.x_l: [self.x_1, 0],
#             self.x_c: [self.x_1, self.x_2, 0],
#             self.x_r: [self.x_2, 0],
#         }

#         return default_data_dict

#     def get_random_parameters(self):

#         default_data_dict = self.get_default_data()

#         parameters_dict = {
#             key: random.choice(items_list)
#             for key, items_list in default_data_dict.items()
#         }

#         if parameters_dict[self.x_l] != self.x_1 or parameters_dict[
#                 self.x_c] != self.x_1:

#             parameters_dict[self.x_l] = self.x_1

#         if parameters_dict[self.x_c] != self.x_2 or parameters_dict[
#                 self.x_r] != self.x_2:

#             parameters_dict[self.x_r] = self.x_2

#         return parameters_dict



    
    
#DO ZROBIENIA?, SPRAWDZIC CZY TO NIE POWTORKA ZADAN O PODOBNEJ NAZWIE!
class ForcedDisksWithSerialSprings(ComposedSystem):
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


    def get_default_data(self):


        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],

            self.phi_1:[self.phi_u],
            self.phi_2:[self.phi_l],
        }

        return default_data_dict    

#     def get_random_parameters(self):

#         default_data_dict = self.get_default_data()

#         parameters_dict = {
#             key: random.choice(items_list)
#             for key, items_list in default_data_dict.items()
#         }

#         if parameters_dict[self.phi_2] == parameters_dict[self.phi_3]:

#             parameters_dict[self.phi_2] = self.phi_u


#         return parameters_dict

    


#DO ZROBIENIA
class LinearizedTriplePendulum(ComposedSystem):
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

