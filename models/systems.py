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
    _default_args = ()

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
    
    def linearized(self, x0=None, op_point=False, hint=[], label=None):

        return type(self).from_system(super().linearized(x0=x0,op_point=op_point,hint=hint,label=label))


#ZAPYTAC SIE BOGUSIA CO Z TYM
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


#DO SPRAWDZENIA
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



#DO SPRAWDZENIA
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



#DO SPRAWDZENIA
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

    
    
#DO SPRAWDZENIA
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


    
    

#DO SPRAWDZENIA
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

#DO ZROBIENIA
class SDOFWinchSystem(ComposedSystem):


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

    
#DO ZROBIENIA
class SDOFDrivetrainVehicleSystem(ComposedSystem):


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

#DO ZROBIENIA
class LagrangeIBlocksOnInclinedPlane(ComposedSystem):
    scheme_name = 'ddof_disks_3_springs_scheme.png'
    real_name = 'nonlin_trolley_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_1', positive=True),
                 m2=Symbol('m_2', positive=True),
                 m3=Symbol('m_3', positive=True),
                 m4=Symbol('m_4', positive=True),
                 R=Symbol('R', positive=True),
                 g=Symbol('g', positive=True),
                 alpha=Symbol('alpha',positive=True),
                 beta=Symbol('beta',positive=True),
                 ivar=Symbol('t'),
                 x1=dynamicsymbols('x_1'),
                 x2=dynamicsymbols('x_2'),
                 x3=dynamicsymbols('x_3'),
                 x4=dynamicsymbols('x_4'),
                 phi=dynamicsymbols('\\varphi'),
                 qs=dynamicsymbols('x_1, x_2, x_3, x_4, \\varphi'),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4

        self.Mass1 = MaterialPoint(m1, pos1=x1, qs=[x1]) + GravitationalForce(m1, g, pos1=-x1*sin(alpha), qs=[x1])
        self.Mass2 = MaterialPoint(m2, pos1=x2, qs=[x2]) + GravitationalForce(m2, g, pos1=-x2*sin(alpha), qs=[x2])
        self.Mass3 = MaterialPoint(m3, pos1=x3, qs=[x3]) + GravitationalForce(m3, g, pos1=-x3*sin(beta), qs=[x3])
        self.Mass4 = MaterialPoint(m4, pos1=x4, qs=[x4]) + GravitationalForce(m4, g, pos1=-x4*sin(beta), qs=[x4])
        self.Pulley = MaterialPoint(1/2*m*R**2, pos1=phi, qs=[phi])

        system = self.Mass1 + self.Mass2 + self.Mass3 + self.Mass4 + self.Pulley
        super().__init__(system.lagrangian(),qs=qs, hol_coneqs=[x1-x2,phi*R-x2,x2-x3,phi*R-x3],**kwargs)
    def get_default_data(self):

        m0 = symbols('m_0', positive=True)

        default_data_dict = {
            self.m: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m2: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m3: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
            self.m4: [S.Half * m0, 1 * m0, 2 * m0, 4 * m0, S.Half**2 * m0],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        return parameters_dict