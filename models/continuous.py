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

        path = __file__.replace('continuous.py', 'images/') + cls.scheme_name

        return path

    @classmethod
    def _real_example(cls):

        path = __file__.replace('continuous.py', 'images/') + cls.real_name

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