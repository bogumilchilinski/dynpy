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