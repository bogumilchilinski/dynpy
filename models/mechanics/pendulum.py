from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin


import base64
import random
import IPython as IP
import numpy as np
import inspect

from .trolley import ComposedSystem, NonlinearComposedSystem

#DONE
class Pendulum(NonlinearComposedSystem):
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
    #scheme_name = 'undamped_pendulum.png'
    #real_name = 'pendulum_real.jpg'

    
    m=Symbol('m', positive=True)
    g=Symbol('g', positive=True)
    l=Symbol('l', positive=True)
    angle=dynamicsymbols('\\varphi')
    qs=None
    
    def __init__(self,
                 m=None,
                 g=None,
                 l=None,
                 angle=None,
                 ivar=None,
                 **kwargs):

        #if qs == None:
        #    qs = [angle]
        #else:
        #    qs = qs

        if m is not None: self.m = m
        if g is not None: self.g = g
        if l is not None: self.l = l
        if angle is not None: self.angle = angle
        if ivar is not None: self.ivar = ivar

        
        self.qs = [self.angle]

        self._init_from_components(**kwargs)
        
        
        
    @property
    def components(self):

        components = {}
        
        self.gravitationalforce = GravitationalForce(self.m, self.g, self.l * (1 - cos(self.angle)), qs = self.qs)
        self.material_point = MaterialPoint(self.m * self.l**2, pos1=self.angle, qs=self.qs)
        
        components['gravitationalforce'] = self.gravitationalforce        
        components['material_point'] = self.material_point

        
        return components
        
    def get_default_data(self):

       
        m0, l0 = self.m0, self.l0
        
        
        default_data_dict = {

            self.m: [S.One * no * m0 * 10 for no in range(5, 8)],
            self.l: [S.One * no * l0 for no in range(10, 20)],
            
        }
        return default_data_dict

    def get_numerical_data(self):

       
        m0, l0 = self.m0, self.l0
        
        
        default_data_dict = {

            self.m: [S.One * no * 10 for no in range(5, 8)],
            self.l: [S.One * no for no in range(10, 20)],
            
        }
        return default_data_dict
    
    
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            #         mech_comp.LinearizationComponent,
            #         mech_comp.GoverningEquationComponent,
            #         mech_comp.FundamentalMatrixComponent,
            #         mech_comp.GeneralSolutionComponent,
            #         mech_comp.SteadySolutionComponent,
        ]

        return comp_list

#DONE
class PulledPendulum(ComposedSystem):
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
        >>> qs = dynamicsymbols('varphi') --> Generalized Coordinates
        >>> Pendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class Pendulum()
    """
    scheme_name = 'undamped_pendulum.png'
    real_name = 'pendulum_real.jpg'

    
    m=Symbol('m', positive=True),
    g=Symbol('g', positive=True),
    l=Symbol('l', positive=True),
    angle=dynamicsymbols('\\varphi'),
    qs=None,
    ivar=Symbol('t')

    
    def __init__(self,
                 m=None,
                 g=None,
                 l=None,
                 angle=None,
                 ivar=None,
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if l is not None: self.l = l
        if angle is not None: self.angle = angle
        if ivar is not None: self.ivar = ivar
        

        
        self.qs = [self.angle]

        self._init_from_components(**kwargs)
        
    @property
    def components(self):
        
        components = {}
        
        self._pendulum = Pendulum(self.m,self.phi,qs=[phi])
        self._force = Force(-2 * self.m * self.l * (self.g / self.l * cos(pi)),
                           angle,qs=[phi])
        
        
        components['_pendulum'] = self._pendulum
        components['force'] = self.force

        
        return components

    def get_default_data(self):

        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
        }
        return self.sym_desc_dict

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            #         mech_comp.LinearizationComponent,
            mech_comp.GoverningEquationComponent,
            mech_comp.FundamentalMatrixComponent,
            mech_comp.GeneralSolutionComponent,
            mech_comp.SteadySolutionComponent,
        ]
        return comp_list


# wymienić obrazek na taki, gdzie nie ma wymuszenia i symbole na obrazku będą zgodne z tymi w klasie

#DONE
class FreePendulum(Pendulum):
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
        >>> FreePendulum()

        -We define the symbols and dynamicsymbols
        -if dynamicsymbols is not defined that parameter would be set as "varphi" as a default
        -determine the instance of the pendulum by using class SFoDFreePendulum()
    """
    scheme_name = 'free_sdof_pendulum.png'
    real_name = 'pendulum_real.jpg'


#TO_DO
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
    scheme_name = 'damped_excited_pendulum.PNG'
    real_name = 'pendulum2_real.jpg'

    def __init__(self,
                 dummy=Symbol('dummy', positive=True),
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 F=Symbol('F', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        phi = angle
        self.phi = phi

        if qs == None:
            qs = [angle]
        else:
            qs = qs

        self.m = m
        self.g = g
        self.l = l
        self.F = F

        Omega = Symbol('Omega', positive=True)
        self.Omega = Omega
        self.pendulum = Pendulum(m, g, l, angle=phi)
        self.force = Force(-F * l * sin(Omega * ivar), pos1=phi, qs=qs)
        system = self.pendulum + self.force

        super().__init__(system, **kwargs)

    def get_default_data(self):

        m0, l0, g0, F0 = symbols('m_0 l_0 g_0 F_0', positive=True)

        default_data_dict = {
            self.m:
            [2 * m0, 1 * m0, S.Half * m0, S.Half**2 * m0, 3 * S.Half * m0],
            self.l:
            [2 * l0, 1 * l0, S.Half * l0, S.Half**2 * l0, 3 * S.Half * l0],
            self.g: [g0],
            self.F:
            [2 * F0, 1 * F0, S.Half * F0, S.Half**2 * F0, 3 * S.Half * F0],
        }
        return default_data_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m1: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.F: r'Force',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, F0 = symbols('m_0 l_0 F_0', positive=True)

        default_data_dict = {
            self.m1: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
            self.F: [
                1 * F0, 2 * F0, 3 * F0, 4 * F0, 5 * F0, 6 * F0, 7 * F0, 8 * F0,
                9 * F0, 10 * F0, 11 * F0, 12 * F0, 13 * F0, 14 * F0, 15 * F0,
                16 * F0, 17 * F0, 18 * F0, 19 * F0, 20 * F0, 21 * F0, 22 * F0,
                23 * F0, 24 * F0, 25 * F0, 26 * F0, 27 * F0, 28 * F0, 29 * F0,
                30 * F0
            ],
        }
        return default_data_dict

#TO_DO
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

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 c=Symbol('c', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
        phi = angle
        self.phi = phi

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

        super().__init__(system, **kwargs)

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

#TO_DO
class ExcitedDampedPendulum(ComposedSystem):

    scheme_name = 'damped_excited_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    def __init__(self,
                 m=Symbol('m', positive=True),
                 g=Symbol('g', positive=True),
                 l=Symbol('l', positive=True),
                 c=Symbol('c', positive=True),
                 F=Symbol('F', positive=True),
                 Omega=Symbol('Omega', positive=True),
                 angle=dynamicsymbols('varphi'),
                 qs=None,
                 ivar=Symbol('t'),
                 **kwargs):
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

        super().__init__(system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of pendulum',
            self.g: r'Gravity constant',
            self.l: r'Pendulum length',
            self.c: r'Damping coefficient',
        }
        return self.sym_desc_dict

#TO_DO
class PendulumKinematicExct(ComposedSystem):

    scheme_name = 'kin_exct_pendulum.PNG'
    real_name = 'elastic_pendulum_real.PNG'

    l = Symbol('l', positive=True)
    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    phi = dynamicsymbols('\\varphi')
    x_e = dynamicsymbols('x_e')

    def __init__(self,
                 l=None,
                 m=None,
                 g=None,
                 phi=None,
                 x_e=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if l is not None: self.l = l
        if m is not None: self.m = m
        if g is not None: self.g = g
        if phi is not None: self.phi = phi
        if x_e is not None: self.x_e = x_e
        self.ivar = ivar
        self.qs = [self.phi]

        self.x = self.l * sin(self.phi) + self.x_e
        self.y = self.l * cos(self.phi)

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._material_point_1 = MaterialPoint(self.m, self.x, qs=self.qs)
        self._material_point_2 = MaterialPoint(self.m, self.y, qs=self.qs)
        self._gravity = GravitationalForce(self.m, self.g, pos1=-self.y, qs=self.qs)

        components['_material_point_1'] = self._material_point_1
        components['_material_point_2'] = self._material_point_2
        components['_gravity'] = self._gravity

        return components


    def symbols_description(self):
        self.sym_desc_dict = {
            self.l: r'Pendulum length',
            self.x_e: r'Kinematic lateral excitation',
            self.m: r'Mass',
            self.g: 'Gravity constant',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, l0, x0, Omega = symbols('m_0 l_0 x_0 Omega', positive=True)

        default_data_dict = {
            self.m: [
                1 * m0, 2 * m0, 3 * m0, 4 * m0, 5 * m0, 6 * m0, 7 * m0, 8 * m0,
                9 * m0, 10 * m0, 11 * m0, 12 * m0, 13 * m0, 14 * m0, 15 * m0,
                16 * m0, 17 * m0, 18 * m0, 19 * m0, 20 * m0, 21 * m0, 22 * m0,
                23 * m0, 24 * m0, 25 * m0, 26 * m0, 27 * m0, 28 * m0, 29 * m0,
                30 * m0
            ],
            self.l: [
                1 * l0, 2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0, 7 * l0, 8 * l0,
                9 * l0, 10 * l0, 11 * l0, 12 * l0, 13 * l0, 14 * l0, 15 * l0,
                16 * l0, 17 * l0, 18 * l0, 19 * l0, 20 * l0, 21 * l0, 22 * l0,
                23 * l0, 24 * l0, 25 * l0, 26 * l0, 27 * l0, 28 * l0, 29 * l0,
                30 * l0
            ],
            self.x_e: [x0 * sin(self.Omega * self.ivar)]
        }
        return default_data_dict

    def max_static_cable_force(self):
        return (self.m * self.g).subs(self._given_data)

    def max_dynamic_cable_force(self):

        omg_amp = ComposedSystem(
            self.linearized()).frequency_response_function() * self.Omega

        return (self.m * self.l * (omg_amp)**2 + self.max_static_cable_force())

    def static_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_cable_force()) / (pi * kr * Re))**(1 / 2)

    def dynamic_cable_diameter(self):
        kr = Symbol('k_r', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_cable_force()) / (pi * kr * Re))**(1 / 2)

