from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, EngineMount
from ..continuous import ContinuousSystem, PlaneStressProblem
from .tmd import TunedMassDamper

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin
from functools import cached_property

#####

class EngineHousing(MaterialPoint):
    """
    Model of a EngineHousing composed of material point with changing point of mass:

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'

    def get_default_data(self):
        m0=Symbol('m_0',positive=True)
        default_data_dict={
            self.m:[S.One *10 * m0 * no  for no in range(75, 120)],
        }
        return default_data_dict

    def get_numerical_data(self):

        self.default_data_dict={
            self.m:[10* no/100 for no in range(80,120)],
        }
        return self.default_data_dict

    
#DONE #Mateusz
############################### Nowa podstawowa klasa dla silników o 1 stopniu swobody
class FreeEngine(ComposedSystem):
    """
    Base class for all engine systems. This class consists only of the engine block and crank system, involving the gravity factor.
    =========
    M = Symbol object
        - Mass of the engine housing
    m_e = Symbol object
        - Reduced mass of the crank system
    e = Symbol object
        - Length of the crank - radius of the circular motion
    g = Symbol object
        - Gravitational field acceleration
    z = Symbol object
        - Vertical coordinate z
    ivar = Symbol object
        - Independant time variable
    qs = Dynamicsymbol object
        - Generalized coordinates
    
    """
    scheme_name = 'engine_block.png'
    real_name = 'engine_real.PNG'
    
    M = Symbol('M', positive=True)  #Mass of the engine housing
    g = Symbol('g', positive=True)  #Gravity constant
    z = dynamicsymbols('z')  #Displacement coordinate
    m_e = Symbol('m_e', positive=True)  #Reduced mass of the crank system
    e = Symbol('e', positive=True)  #Length of the crank - radius of the circular motion
    phi = dynamicsymbols('\\varphi')
    ivar = Symbol('t')
    m0 = Symbol('m_0', positive=True)

    e0 = Symbol('e_0', positive=True)
    omega = Symbol('Omega', positive=True)    

    def __init__(self,
                 M=None,
                 m_e=None,
                 e=None,
                 g=None,
                 z=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if g is not None: self.g = g
        if z is not None: self.z = z
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        M, m_e, e, z, phi = self.M, self.m_e, self.e, self.z, self.phi

        self._engine_housing = EngineHousing(M, pos1=z, qs=[z])(label='Engine housing')
        self._crank = MaterialPoint(m_e, pos1=z + e * cos(phi), qs=[z])(label='Position of reduced mass of the crank system')
        self._housing_gravity = GravitationalForce(m=self.M, g=self.g, pos1=self.z, qs=[self.z])(label='Gravitational Force of housing')
        self._crank_gravity = GravitationalForce(m=self.m_e, g=self.g, pos1=self.z+self.e*cos(self.phi), qs=[self.z])(label='Gravitational Force of crank')

        components['_engine_housing'] = self._engine_housing
        components['_crank'] = self._crank
        components['_housing_gravity'] = self._housing_gravity
        components['_crank_gravity'] = self._crank_gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of the engine block',
            self.m_e: r'Reduced mass of the crank system',
            self.e: r'Length of the crank - radius of the circular motion',
            self.g: r'Gravity constant',
            self.z: r'Displacement coordinate'
        }
        return self.sym_desc_dict

    def get_default_data(self):

        t = self.ivar
        
        default_data_dict = {
#             self.M: [10 * self.m0 * no / 100 for no in range(10, 150)],
#             self.m_e: [self.m0 * no / 100 for no in range(80, 120)],
            self.e: [S.One *self.e0 * no  for no in range(8, 12)],
            self.phi: [self.omega * self.ivar]
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
#             self.M: [4 * no for no in range(10, 150)],
#             self.m_e: [0.5 * no for no in range(80, 120)],
            self.e: [2/100 * no for no in range(5, 15)],
            self.phi: [2 * 3.14 * 10 * self.ivar],
#             self.g: [9.81]
        }
        return default_data_dict
    
    
#DONE #Mateusz
############################### Nowa podstawowa klasa dla silników o 1 stopniu swobody
class FreeEngineWithConstrains(ComposedSystem):
    """
    Base class for all engine systems. This class consists only of the engine block and crank system, involving the gravity factor.
    =========
    M = Symbol object
        - Mass of the engine housing
    m_e = Symbol object
        - Reduced mass of the crank system
    e = Symbol object
        - Length of the crank - radius of the circular motion
    g = Symbol object
        - Gravitational field acceleration
    z = Symbol object
        - Vertical coordinate z
    ivar = Symbol object
        - Independant time variable
    qs = Dynamicsymbol object
        - Generalized coordinates
    
    """
    scheme_name = 'engine_block.png'
    real_name = 'engine_real.PNG'
    
    M = Symbol('M', positive=True)  #Mass of the engine housing
    g = Symbol('g', positive=True)  #Gravity constant
    z = dynamicsymbols('z')  #Displacement coordinate
    z_e = dynamicsymbols('z_e')  #Displacement coordinate
    m_e = Symbol('m_e', positive=True)  #Reduced mass of the crank system
    e = Symbol('e', positive=True)  #Length of the crank - radius of the circular motion
    phi = dynamicsymbols('\\varphi')
    ivar = Symbol('t')
    m0 = Symbol('m_0', positive=True)

    e0 = Symbol('e_0', positive=True)
    omega = Symbol('Omega', positive=True)    

    def __init__(self,
                 M=None,
                 m_e=None,
                 e=None,
                 g=None,
                 z=None,
                 z_e=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if g is not None: self.g = g
        if z is not None: self.z = z
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z]

        self._init_from_components({**kwargs,'hol_coneqs':[self.z_e - self.z -10]})

    @cached_property
    def components(self):
        components = {}
        M, m_e, e, z, phi, z_e = self.M, self.m_e, self.e, self.z, self.phi, self.z_e

        self._engine_housing = EngineHousing(M, pos1=z, qs=[z])(label='Engine housing')
        self._crank = MaterialPoint(m_e, pos1=z_e, qs=[z_e])(label='Position of reduced mass of the crank system')
        self._housing_gravity = GravitationalForce(m=self.M, g=self.g, pos1=self.z, qs=[self.z])(label='Gravitational Force of housing')
        self._crank_gravity = GravitationalForce(m=self.m_e, g=self.g, pos1=z_e, qs=[self.z_e])(label='Gravitational Force of crank')

        components['_engine_housing'] = self._engine_housing
        components['_crank'] = self._crank
        components['_housing_gravity'] = self._housing_gravity
        components['_crank_gravity'] = self._crank_gravity

        return components
    
    
#DONE #Mateusz
################### Podstawowa klasa dla silników o 2 stopniach swobody
class FreeEngineDDOF(FreeEngine):
    """Base class for all engine systems with 2 degrees of freedom. This class consists only of the engine block and crank system, involving the gravity factor.
    =========
            M = Symbol object
                - Mass of the engine housing
            m_e = Symbol object
                - Reduced mass of the crank system
            e = Symbol object
                - Length of the crank - radius of the circular motion
            g = Symbol object
                - Gravitational field acceleration
            x = Symbol object
                - Vertical coordinate x,,
            z = Symbol object
                - Vertical coordinate z
            ivar = Symbol object
                - Independant time variable
            qs = Dynamicsymbol object
                - Generalized coordinates
    
    """
    scheme_name = 'engine_block_DDOF.png'
    real_name = 'engine_real.PNG'

    M = Symbol('M', positive=True)  #Mass of the engine housing
    g = Symbol('g', positive=True)  #Gravity constant
    z = dynamicsymbols('z')  #Displacement coordinate z
    x = dynamicsymbols('x')  #Displacement coordinate x
    m_e = Symbol('m_e', positive=True)  #Reduced mass of the crank system
    e = Symbol('e', positive=True)  #Length of the crank - radius of the circular motion
    phi = dynamicsymbols('\\varphi')
    ivar = Symbol('t')
    m0 = Symbol('m0', positive=True)
    m_e0 = Symbol('m_e0', positive=True)
    e0 = Symbol('e_0', positive=True)

    def __init__(self,
                 M=None,
                 m_e=None,
                 e=None,
                 g=None,
                 z=None,
                 x=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if g is not None: self.g = g
        if z is not None: self.z = z
        if x is not None: self.x = x
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z, self.x]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        M, m_e, e, g, z, x, phi = self.M, self.m_e, self.e, self.g, self.z, self.x, self.phi

        self._engine_block_1DOF = FreeEngine(M, m_e, e, g, z, phi)(label='Engine block with 1 DOF')
        self._engine_housing_horizontal_component = EngineHousing(M, pos1=x, qs=[x])(label='Engine hounsing horizontal component')
        self._crank_horizontal_component = MaterialPoint(m_e, pos1=x + e * sin(phi), qs=[x])(label='Position of reduced mass of the crank system in direction x')

        components['_engine_block_1DOF'] = self._engine_block_1DOF
        components['_engine_housing_horizontal_component'] = self._engine_housing_horizontal_component
        components['_crank_horizontal_component'] = self._crank_horizontal_component

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of the engine block',
            self.m_e: r'Reduced mass of the crank system',
            self.e: r'Length of the crank - radius of the circular motion',
            self.g: r'Gravity constant',
            self.z: r'Displacement coordinate z',
            self.x: r'Displacement coordinate x'
        }
        return self.sym_desc_dict

#     def get_default_data(self):

#         t = self.ivar
        
#         default_data_dict = {
#             self.M: [10 * self.m0 * no / 100 for no in range(10, 150)],
#             self.m_e: [self.m0 * no / 100 for no in range(80, 120)],
#             self.e: [self.e0 * no / 100 for no in range(80, 120)],
#             self.phi: [2 * 3.14 * self.ivar]
#         }
#         return default_data_dict

#     def get_numerical_data(self):

#         default_data_dict = {
#             self.M: [4 * no for no in range(10, 150)],
#             self.m_e: [0.5 * no for no in range(80, 120)],
#             self.e: [2/100 * no for no in range(5, 15)],
#             self.phi: [2 * 3.14 * 10 * self.ivar],
#             self.g: [9.81]
#         }
#         return default_data_dict
#dobrać dane numeryczne
#####
#DONE #Sav 
class Engine(ComposedSystem):
    """
    Ready to use model of engine represented by the rotating mass of a crankshaft and mass of the engine.
    ==========
    
    Arguments:
    ==========

        M = Symbol object
            -Mass of an engine.

        m_e = Symbol object
            -Mass of a crankshaft.

        k_m = Symbol object
            -Stifness of the engine mounts
            
        phi = Dynamicsymbol object
            -rotation angle
            
        e =  Symbol object
            -offset of the rotating mass

        g = Symbol object
            -Gravitational field acceleration
            
        z = Symbol object
            -vertical Z coordinate
            
        q = Symbol object
            -Directional coefficent

        ivar = symbol object
            -Independant time variable

        qs = dynamicsymbol object
            -Generalized coordinates


    """
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'

    M = Symbol('M', positive=True)  #Mass of system (engine block) on spring
    k_m = Symbol('k_m', positive=True)  #stiffness
    m_e = Symbol('m_e', positive=True)  #Mass of particle
    g = Symbol('g', positive=True)
    e = Symbol('e', positive=True)  #distance -motion radius of a particle
    z = dynamicsymbols('z')  #generalized coordinate
    phi = dynamicsymbols('varphi')
    ivar = Symbol('t')
    m0 = Symbol('m_0', positive=True)
    k_m0 = Symbol('k_m0', positive=True)
    m_e0 = Symbol('m_e0', positive=True)
    e0 = Symbol('e_0', positive=True)
    omega = Symbol('Omega', positive=True)

    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 g=None,
                 e=None,
                 z=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if g is not None: self.g = g
        if e is not None: self.e = e
        if z is not None: self.z = z
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}
        M, k_m, m_e, e, z, phi = self.M, self.k_m, self.m_e, self.e, self.z, self.phi

        self._engine = FreeEngine(self.M, self.m_e,self.e ,self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine')
        
        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount',scheme_options={'at':(1,-8)})
        
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount',scheme_options={'at':(4,-8)})

        components['_engine'] = self._engine
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount

        return components

    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
#        mech_comp.FundamentalMatrixComponent,
#        mech_comp.GeneralSolutionComponent,
#        mech_comp.SteadySolutionComponent,
#        mech_comp.SpringForce,
            
        ]
        
        return comp_list
    
#DONE 
class EngineVerticalSpringGravity(Engine):

    """Ready to use model of engine represented by the rotating mass of a crankshaft and mass of the engine.
        Arguments:
        =========
            M = Symbol object
                -Mass of an engine.

            m_e = Symbol object
                -Mass of a crankshaft.

            k_m = Symbol object
                -Stifness of the engine mounts
                
                
            phi = Dynamicsymbol object
                -rotation angle
                
            e =  Symbol object
                -offset of the rotating mass

            g = Symbol object
                -Gravitational field acceleration
                
            z = Symbol object
                -vertical Z coordinate
                
            q = Symbol object
                -Directional coefficent

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates


    """

    scheme_name = 'engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'

    M = Symbol('M', positive=True)
    m_e = Symbol('m_e', positive=True)
    phi = dynamicsymbols('varphi')
    g = Symbol('g', positive=True)
    k_m = Symbol('k_m', positive=True)
    c_m = Symbol('c_m', positive=True)
    e = Symbol('e', positive=True)
    z = dynamicsymbols('z')

    Omega = Symbol('Omega', positive=True)

    #m0 = Symbol('m_0', positive=True)

    def __init__(self,
                 M=None,
                 m_e=None,
                 k_m=None,
                 e=None,
                 g=None,
                 phi=None,
                 z=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if phi is not None: self.phi = phi
        if g is not None: self.g = g
        if k_m is not None: self.k_m = k_m
        if e is not None: self.e = e
        if z is not None: self.z = z

        self.qs = [self.z]
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._engine_block = FreeEngine(self.M, self.m_e,self.e, self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine', scheme_options = {'at':(0,0)})

        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount', scheme_options = {'at':(0,-8)})
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount', scheme_options = {'at':(5,-8)})
        

        components['_engine_block'] = self._engine_block
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        return components

class EngineConstantVelocityVerticalSpringGravity(Engine):

    """Ready to use model of engine represented by the rotating mass of a crankshaft and mass of the engine.
        Arguments:
        =========
            M = Symbol object
                -Mass of an engine.

            m_e = Symbol object
                -Mass of a crankshaft.

            k_m = Symbol object
                -Stifness of the engine mounts
                
                
            phi = Dynamicsymbol object
                -rotation angle
                
            e =  Symbol object
                -offset of the rotating mass

            g = Symbol object
                -Gravitational field acceleration
                
            z = Symbol object
                -vertical Z coordinate
                
            q = Symbol object
                -Directional coefficent

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates


    """

    scheme_name = 'engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'

    M = Symbol('M', positive=True)
    m_e = Symbol('m_e', positive=True)
    g = Symbol('g', positive=True)
    k_m = Symbol('k_m', positive=True)
    c_m = Symbol('c_m', positive=True)
    e = Symbol('e', positive=True)
    z = dynamicsymbols('z')

    Omega = Symbol('Omega', positive=True)

    #m0 = Symbol('m_0', positive=True)

    def __init__(self,
                 M=None,
                 m_e=None,
                 k_m=None,
                 e=None,
                 g=None,
                 z=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if g is not None: self.g = g
        if k_m is not None: self.k_m = k_m
        if e is not None: self.e = e
        if z is not None: self.z = z
        if Omega is not None: self.Omega=Omega

        self.qs = [self.z]
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._engine_block = FreeEngine(self.M, self.m_e,self.e, self.g, z= self.z, phi = self.Omega*self.ivar, qs=[self.z])(label='Engine', scheme_options = {'at':(0,0)})

        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount', scheme_options = {'at':(0,-8)})
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount', scheme_options = {'at':(5,-8)})
        

        components['_engine_block'] = self._engine_block
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        return components
    
    # def get_default_data(self):

    #     m0, e0, Omega0 = symbols('m_0 e_0 Omega_0',positive=True)
        
    #     default_data_dict = {
    #         self.M: [10 * m0 * no / 100 for no in range(10, 150)],
    #         self.m_e: [m0 * no / 100 for no in range(80, 120)],
    #         self.e: [e0 * no / 100 for no in range(80, 120)],
    #     }
    #     return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.M: [4 * no for no in range(10, 100)],
            self.m_e: [0.5 * no for no in range(80, 120)],
            self.e: [2/100 * no for no in range(5, 15)],
            self.Omega: [2 * 3.14 * no for no in range(1,20)],
            self.g: [9.81]
        }
        return default_data_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        mech_comp.SteadySolutionComponent,
        mech_comp.SpringForce,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticPinDiameter,
        mech_comp.DynamicPinDiameter,

        ]

        return comp_list

    
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def dynamic_force(self):
        return self.spring_force().expand().doit().n(6)
    
    def spring_force(self):
        k_m=self.k_m
        z=self.z
        sol_dict=self._fodes_system.steady_solution.as_dict()
        F_km=(k_m*z).subs(sol_dict).subs(self._given_data)
        
        return F_km
    
    
    def left_engine_mount_force(self):
        return self._left_mount.stiffness * self.steady_solution()[0]

    def right_engine_mount_force(self):
        return self._right_mount.stiffness * self.steady_solution()[0]

    def max_static_force(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force(self):
        return self.frequency_response_function() * self._left_mount.stiffness + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)
    

class EngineConstantWithReactionForce(EngineConstantVelocityVerticalSpringGravity):

    """Ready to use model of engine represented by the rotating mass of a crankshaft and mass of the engine.
        Arguments:
        =========
            M = Symbol object
                -Mass of an engine.

            m_e = Symbol object
                -Mass of a crankshaft.

            k_m = Symbol object
                -Stifness of the engine mounts
                
                
            phi = Dynamicsymbol object
                -rotation angle
                
            e =  Symbol object
                -offset of the rotating mass

            g = Symbol object
                -Gravitational field acceleration
                
            z = Symbol object
                -vertical Z coordinate
                
            q = Symbol object
                -Directional coefficent

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates


    """

    scheme_name = 'engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'

    M = Symbol('M', positive=True)
    m_e = Symbol('m_e', positive=True)
    g = Symbol('g', positive=True)
    k_m = Symbol('k_m', positive=True)
    c_m = Symbol('c_m', positive=True)
    e = Symbol('e', positive=True)
    z = dynamicsymbols('z')

    Omega = Symbol('Omega', positive=True)

    #m0 = Symbol('m_0', positive=True)

    def __init__(self,
                 M=None,
                 m_e=None,
                 k_m=None,
                 e=None,
                 g=None,
                 z=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if m_e is not None: self.m_e = m_e
        if g is not None: self.g = g
        if k_m is not None: self.k_m = k_m
        if e is not None: self.e = e
        if z is not None: self.z = z
        if Omega is not None: self.Omega=Omega

        self.qs = [self.z]
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._engine_block = FreeEngine(self.M, self.m_e,self.e, self.g, z= self.z, phi = self.Omega*self.ivar, qs=[self.z])(label='Engine', scheme_options = {'at':(0,0)})

        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount', scheme_options = {'at':(0,-8)})
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount', scheme_options = {'at':(5,-8)})
        

        components['_engine_block'] = self._engine_block
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        return components
    
    # def get_default_data(self):

    #     m0, e0, Omega0 = symbols('m_0 e_0 Omega_0',positive=True)
        
    #     default_data_dict = {
    #         self.M: [10 * m0 * no / 100 for no in range(10, 150)],
    #         self.m_e: [m0 * no / 100 for no in range(80, 120)],
    #         self.e: [e0 * no / 100 for no in range(80, 120)],
    #     }
    #     return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.M: [4 * no for no in range(10, 100)],
            self.m_e: [0.5 * no for no in range(80, 120)],
            self.e: [2/100 * no for no in range(5, 15)],
            self.Omega: [2 * 3.14 * no for no in range(1,20)],
            self.g: [9.81]
        }
        return default_data_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        mech_comp.SteadySolutionComponent,
        mech_comp.SpringForce,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticPinDiameter,
        mech_comp.DynamicPinDiameter,

        ]

        return comp_list
    
    def left_engine_mount_force(self):
        return self._left_mount.stiffness * self.steady_solution()[0]

    def right_engine_mount_force(self):
        return self._right_mount.stiffness * self.steady_solution()[0]

    def max_static_force(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force(self):
        return self.frequency_response_function() * self._left_mount.stiffness + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)
    
    
    

#TODO #DDOF #DO POPRAWY !!!
#class BoxerEnginePerpendicularSprings(Engine):
#    scheme_name = 'boxer_engine_perpendicular_springs.png'
#    real_name = 'f6c_valkyrie.jpg'
#
#    M = Symbol('M', positive=True),
#    k_m = Symbol('k_m', positive=True),
#    m_e = Symbol('m_e', positive=True),
#    e = Symbol('e', positive=True),
#    l = Symbol('l', positive=True),
#    x = dynamicsymbols('x'),
#    z = dynamicsymbols('z'),
#    phi = dynamicsymbols('phi'),
#    ivar = Symbol('t'),
#
#    def __init__(self,
#                 M=None,
#                 k_m=None,
#                 m_e=None,
#                 e=None,
#                 l=None,
#                 x=None,
#                 z=None,
#                 phi=None,
#                 ivar=Symbol('t'),
#                 **kwargs):
#
#        if M is not None: self.M = M
#        if k_m is not None: self.k_m = k_m
#        if m_e is not None: self.m_e = m_e
#        if e is not None: self.e = e
#        if l is not None: self.l = l
#        if x is not None: self.x = x
#        if z is not None: self.z = z
#        if phi is not None: self.phi = phi
#
#        self._init_from_components(**kwargs)
#
#        @property
#        def components(self):
#
#            components = {}
#
#            self._engine_block = MaterialPoint(self.M,
#                                               pos1=self.x,
#                                               qs=[self.x
#                                                   ])(label='Engine block')
#            self._unbalanced_mass = MaterialPoint(
#                self.m_e, pos1=self.x + self.e * sin(self.phi),
#                qs=[self.x])(label='Unbalanced mass in cranksystem')
#            self._spring_vertical = Spring(
#                2 * self.k_m, pos1=self.x,
#                qs=[self.x])(label='Vertical engine mount')
#            self._spring_horizontal = Spring(
#                2 * self.k_m, pos1=self.z,
#                qs=[self.x])(label='Horizontal engine mount')
#
#            components['_engine_block'] = self._engine_block
#            components['_unbalanced_mass'] = self._unbalanced_mass
#            components['_spring_vertical'] = self._spring_vertical
#            components['_spring_horizontal'] = self._spring_horizontal
#
#            return components
#
#        def symbols_description(self):
#            self.sym_desc_dict = {
#                self.M: r'Mass of engine block',
#                self.k_m: r'Spring stiffness coefficient',
#                self.m_e: r'Mass of unbalanced mass',
#                self.e: r'radius of unbalanced mass',
#            }
#            return self.sym_desc_dict

#DONE
class DampedEngineVerticalSpringGravity(Engine):
    scheme_name = 'damped_engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'
    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    l = Symbol('l', positive=True)
    z = dynamicsymbols('z')
    Omega = Symbol('Omega', positive=True)
    phi = dynamicsymbols('varphi')
    ivar = Symbol('t')
    c_m = Symbol('c_m', positive=True)
    g = Symbol('g', positive=True)
    lam = Symbol('lambda', positive=True)

    m0 = Symbol('m_0', positive=True)

    def __init__(
            self,
            M=None,
            k_m=None,
            m_e=None,
            e=None,
            l=None,
            z=None,
            Omega=None,
            phi=None,
            c_m=None,
            g=None,
            lam=None,
            ivar=Symbol('t'),
            qs=None,
            **kwargs):

        if phi is not None: self.phi = phi
        if z is not None: self.z = z
        if Omega is not None: self.Omega = Omega
        self.t = ivar
        if M is not None: self.M = M
        if k_m is not None: self.k_m
        if c_m is not None: self.c_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if g is not None: self.g = g
        if lam is not None: self.lam = lam
        self.qs = [self.phi, self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        #print('print', self.M, self.m_e, self.k_m, self.e, self.l, self.z,self.Omega)
        
        self._engine = FreeEngine(self.M, self.m_e, self.e, self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine')
        
        self._left_mount = Spring(self.k_m, pos1=self.z, qs=[self.z])(label='Left engine mount')
        
        self._right_mount = Spring(self.k_m, pos1=self.z, qs=[self.z])(label='Right engine mount')

        self._left_damper = Damper((self.c_m), pos1=self.z, qs=[self.z])(label='Damping in left engine mount')
        
        self._right_damper = Damper((self.c_m), pos1=self.z, qs=[self.z])(label='Damping in right engine mount')
        
        components['_engine'] = self._engine
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        components['_left_damper'] = self._left_damper
        components['_right_damper'] = self._right_damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)
        m0, k0, Omega = self.m0, self.k0, self.Omega

        default_data_dict = {
            self.c_m: [lam * self.k_m],
            self.M: [S.One * m0 * no * 10 for no in range(5, 8)],
            self.m_e: [S.One / 2 * m0 * no for no in range(1, 8)],
            self.k_m: [S.One / 100 * k0 * no for no in range(80, 135)],
            self.e: [S.One / 10 * e0 * no for no in range(7, 14)],
            self.phi: [S.One * Omega * self.ivar],
        }

        return default_data_dict

    def get_numerical_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)
        m0, k0 = self.m0, self.k0

        default_data_dict = {
            self.c_m: [0.01 * self.k_m],
            self.M: [2 * no * 10 for no in range(5, 8)],
            self.m_e: [S.One / 2 * 2 * no for no in range(1, 8)],
            self.k_m: [S.One / 100 * 50e3 * no for no in range(80, 135)],
            self.e: [0.01 * no for no in range(7, 14)],
            self.phi: [S.One / 10  * self.ivar * no for no in range(7, 14)],
            self.g: [9.81],
        }

        return default_data_dict

class DampedEngineConstantVelocityVerticalSpringGravity(EngineConstantVelocityVerticalSpringGravity):
    scheme_name = 'damped_engine_vertical_spring_gravity.png'
    real_name = 'paccar.jpg'
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'
    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    l = Symbol('l', positive=True)
    z = dynamicsymbols('z')
    Omega = Symbol('Omega', positive=True)
    phi = dynamicsymbols('varphi')
    ivar = Symbol('t')
    c_m = Symbol('c_m', positive=True)
    g = Symbol('g', positive=True)
    lam = Symbol('lambda', positive=True)

    m0 = Symbol('m_0', positive=True)

    def __init__(
            self,
            M=None,
            k_m=None,
            m_e=None,
            e=None,
            l=None,
            z=None,
            Omega=None,
            phi=None,
            c_m=None,
            g=None,
            lam=None,
            ivar=Symbol('t'),
            qs=None,
            **kwargs):

        if phi is not None: self.phi = phi
        if z is not None: self.z = z
        if Omega is not None: self.Omega = Omega
        self.t = ivar
        if M is not None: self.M = M
        if k_m is not None: self.k_m
        if c_m is not None: self.c_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if g is not None: self.g = g
        if lam is not None: self.lam = lam
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        #print('print', self.M, self.m_e, self.k_m, self.e, self.l, self.z,self.Omega)
        
        self._engine = FreeEngine(self.M, self.m_e, self.e, self.g, z= self.z, phi = self.Omega*self.ivar, qs=[self.z])(label='Engine')
        
        self._left_mount = Spring(self.k_m, pos1=self.z, qs=[self.z])(label='Left engine mount')
        
        self._right_mount = Spring(self.k_m, pos1=self.z, qs=[self.z])(label='Right engine mount')

        self._left_damper = Damper((self.c_m), pos1=self.z, qs=[self.z])(label='Damping in left engine mount')
        
        self._right_damper = Damper((self.c_m), pos1=self.z, qs=[self.z])(label='Damping in right engine mount')
        
        components['_engine'] = self._engine
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        components['_left_damper'] = self._left_damper
        components['_right_damper'] = self._right_damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)
        m0, k0, Omega = self.m0, self.k0, self.Omega

        default_data_dict = {
            self.c_m: [lam * self.k_m],
            self.M: [S.One * m0 * no * 10 for no in range(5, 8)],
            self.m_e: [S.One / 2 * m0 * no for no in range(1, 8)],
            self.k_m: [S.One / 100 * k0 * no for no in range(80, 135)],
            self.e: [S.One / 10 * e0 * no for no in range(7, 14)],
            self.phi: [S.One * Omega * self.ivar],
        }

        return default_data_dict

    def get_numerical_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)
        m0, k0 = self.m0, self.k0

        default_data_dict = {
            self.c_m: [0.01 * self.k_m],
            self.M: [2 * no * 10 for no in range(5, 8)],
            self.m_e: [S.One / 2 * 2 * no for no in range(1, 8)],
            self.k_m: [S.One / 100 * 50e3 * no for no in range(80, 135)],
            self.e: [0.01 * no for no in range(7, 14)],
            self.phi: [S.One / 10 * self.ivar * no for no in range(7, 14)],
            self.g: [9.81],
        }

        return default_data_dict
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        mech_comp.SteadySolutionComponent,
        mech_comp.SpringForce,
        mech_comp.MaxStaticForce,
        mech_comp.MaxDynamicForce,
        mech_comp.StaticPinDiameter,
        mech_comp.DynamicPinDiameter,

        ]

        return comp_list


##DONE #DDOF #Mateusz
class InlineEnginePerpendicularSprings(Engine):
    scheme_name = 'inline_engine_perpendicular_springs_DDOF.png'
    real_name = 'paccar.jpg'

    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    l = Symbol('l', positive=True)
    x = dynamicsymbols('x')
    z = dynamicsymbols('z')
    phi = dynamicsymbols('varphi')
    ivar = Symbol('t')

    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 e=None,
                 l=None,
                 x=None,
                 z=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if l is not None: self.l = l
        if x is not None: self.x = x
        if z is not None: self.z = z
        if phi is not None: self.phi = phi

        self.qs = [self.x, self.z]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._engine_block_DDOF = FreeEngineDDOF(self.M, self.m_e, self.e, self.g, self.z, self.x, self.phi)
        
        self._spring_ver_left = Spring(self.k_m,pos1=(self.z**2 + self.x**2)**(1 / 2) - self.z,qs=[self.z])(label='Left vertical engine mount')
        self._spring_hor_left = Spring(self.k_m,pos1=(self.z**2 + self.x**2)**(1 / 2) - (self.x + self.l),qs=[self.x])(label='Left horizontal engine mount')
        self._spring_ver_right = Spring(self.k_m,pos1=(self.z**2 + (self.x + self.l)**2)**(1 / 2) - self.z,qs=[self.z])(label='Right vertical engine mount')
        self._spring_hor_right = Spring(self.k_m,pos1=(self.z**2 +(self.x + self.l)**2)**(1 / 2) - (self.x + self.l),qs=[self.x])(label='Right horizontal engine mount')

        components['_engine_block_DDOF'] = self._engine_block_DDOF
        components['_spring_ver_left'] = self._spring_ver_left
        components['_spring_hor_left'] = self._spring_hor_left
        components['_spring_ver_right'] = self._spring_ver_right
        components['_spring_hor_right'] = self._spring_hor_right

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of the engine block',
            self.m_e: r'Reduced mass of the crank system',
            self.k_m: r'Spring stiffness coefficient',
            self.e: r'Length of the crank - radius of the circular motion',
            self.g: r'Gravity constant',
            self.z: r'Displacement coordinate z',
            self.x: r'Displacement coordinate x'
        }
        return self.sym_desc_dict


#Pioter #DONE - PATI zrobione
class NonLinearInlineEnginePerpendicularSpringsGravity(Engine):

    scheme_name = 'nonlinear_inline_engine_perpendicular_springs.png'
    real_name = 'paccar.jpg'
    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    l = Symbol('l', positive=True)
    g = Symbol('g', positive=True)
    z = dynamicsymbols('z')
    phi = dynamicsymbols('phi')
    Omega = Symbol('Omega', positive=True)
    ivar = Symbol('t')
    d = Symbol('d', positive=True)

    e0 = Symbol('e_0', positive=True)

    def __init__(
            self,
            d1=None,
            #d2=None,
            M=None,
            k_m=None,
            m_e=None,
            e=None,
            l=None,
            g=None,
            z=None,
            phi=None,
            Omega=None,
            ivar=None,
            d=None,
            **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if phi is not None: self.phi = phi
        if Omega is not None: self.Omega = Omega
        if d is not None: self.d = d
        if l is not None: self.l = l
        if g is not None: self.g = g
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z]

        self._init_from_components(**kwargs)
    @property
    def components(self):
        components = {}
        M, k_m, m_e, e, z, phi = self.M, self.k_m, self.m_e, self.e, self.z, self.phi

        self._engine = FreeEngine(self.M, self.m_e, self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine')
        
        self._springVer  = Spring(self.k_m, self.z, qs=[self.z])(label='Left engine mount')
        
        self._springHor = Spring(self.k_m, self.z, qs=[self.z])(label='Right engine mount')

        components['_engine'] = self._engine
        components['_springVer'] = self._springVer
        components['_springHor'] = self._springHor

        return components


    def linearized(self):

        return type(self).from_system(super().linearized())

    def get_default_data(self):

        m0, k0, e0, l0 = symbols('m_0 k_0 e_0 l_0', positive=True)

        default_data_dict = {
            self.d: [l0 * (S.One + no / 20) for no in range(5, 10)],
            self.l: [l0],
        }

        return {**super().get_default_data(), **default_data_dict}

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if 4 * parameters_dict[self.k_m] - 2 * parameters_dict[
                self.k_m] * parameters_dict[self.l] / parameters_dict[
                    self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

#DONE-PATI ZROBI
class NonLinearBoxerEnginePerpendicularSprings(Engine):
    scheme_name = 'nonlin_boxer_engine_perpendicular_springs.png'
    real_name = 'f6c_valkyrie.jpg'
    M = Symbol('M', positive=True)
    k_m=Symbol('k_m', positive=True)
    m_e=Symbol('m_e', positive=True)
    e=Symbol('e', positive=True)
    l=Symbol('l', positive=True)
    x=dynamicsymbols('x')
    phi=Symbol('phi', positive=True)
    Omega=Symbol('Omega', positive=True)
    ivar=Symbol('t')
    d=Symbol('d', positive=True)
        
    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 e=None,
                 l=None,
                 x=None,
                 phi=None,
                 Omega=None,
                 ivar=None,
                 d=None,
                 **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if l is not None: self.l = l
        if x is not None: self.x = x            
        if phi is not None: self.phi = phi
        if Omega is not None: self.Omega = Omega
        if d is not None: self.d = d
        if ivar is not None: self.ivar = ivar
        self.qs = [self.x]
        
        self._init_from_components(**kwargs)
        
    @property
    def components(self):
        components = {}
        self.g=0
        M, k_m, m_e, e, x, phi = self.M, self.k_m, self.m_e, self.e, self.x, self.phi
        
        self._engine = FreeEngine(self.M, self.m_e,g=0, z= self.x, phi = self.phi, qs=[self.x])(label='Engine')
        
        self._springVer  = Spring(self.k_m, self.x, qs=[self.x])(label='Left engine mount')
        
        self._springHor = Spring(self.k_m, self.x, qs=[self.x])(label='Right engine mount')

        components['_engine'] = self._engine
        components['_springVer'] = self._springVer
        components['_springHor'] = self._springHor

        return components



    def get_default_data(self):

        m0, k0, e0, l0 = symbols('m_0 k_0 e_0 l_0', positive=True)

        default_data_dict = {
            self.M: [
                200 * m0, 350 * m0, 400 * m0, 550 * m0, 650 * m0, 700 * m0,
                800 * m0
            ],
            self.k_m: [
                2 * k0, 3 * k0, 4 * k0, 5 * k0, 6 * k0, 7 * k0, 8 * k0, 9 * k0,
                10 * k0
            ],
            self.m_e: [
                0.2 * m0, 0.3 * m0, 0.4 * m0, 0.5 * m0, 0.6 * m0, 0.7 * m0,
                0.8 * m0, 0.9 * m0
            ],
            self.e: [2 * e0, 3 * e0, 4 * e0, 5 * e0, 6 * e0],
            self.d: [2 * l0, 3 * l0, 4 * l0, 5 * l0, 6 * l0],
            self.l: [l0],
            #             self.g:[g],
            #             self.phi:[self.Omega*self.t],
            self.phi: [self.Omega * self.t]
            
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }
        if 4 * parameters_dict[self.k_m] - 2 * parameters_dict[
                self.k_m] * parameters_dict[self.l] / parameters_dict[
                    self.d] == 0:
            parameters_dict[self.d] = 2 * parameters_dict[self.d]

        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


#Sav #DONE
class DampedEngine(Engine):
    """Ready to use model of engine represented by the rotating mass of a crankshaft and mass of the engine.
        Arguments:
        =========
            M = Symbol object
                -Mass of an engine.

            m_e = Symbol object
                -Mass of a crankshaft.

            k_m = Symbol object
                -Stifness of the engine mounts
                
            c_m = Symbol object
                -Damping of the engine mounts
                
                
            phi = Dynamicsymbol object
                -rotation angle
                
            e =  Symbol object
                -offset of the rotating mass

            g = Symbol object
                -Gravitational field acceleration
                
            z = Symbol object
                -vertical Z coordinate
                
            q = Symbol object
                -Directional coefficent

            ivar = symbol object
                -Independant time variable

            qs = dynamicsymbol object
                -Generalized coordinates


    """

    scheme_name = 'engine_with_damper.png'
    real_name = 'engine_real.PNG'

    c_m = Symbol('c_m', positive=True)

    def __init__(self,
                 M=None,
                 k_m=None,
                 c_m=None,
                 m_e=None,
                 g=None,
                 e=None,
                 z=None,
                 phi=None,
                 ivar=None,
                 **kwargs):

        if c_m is not None: self.c_m = c_m
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if g is not None: self.g = g
        if e is not None: self.e = e
        if z is not None: self.z = z
        if phi is not None: self.phi = phi
        if ivar is not None: self.ivar = ivar
        self.qs = [self.z]

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components={}
        self._engine = FreeEngine(self.M, self.m_e, self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine')
        
        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount')
        
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount')


        self.damper_left = Damper(2 * self.c_m, pos1=self.z,
                                  qs=[self.z])(label='Damping in left mount')
        self.damper_right = Damper(2 * self.c_m, pos1=self.z,
                                   qs=[self.z])(label='Damping in right mount')

        
        components['_engine'] = self._engine
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount
        components['damper_left'] = self.damper_left
        components['damper_right'] = self.damper_right

        return components

#DONE #Amadi #Sav
class NonlinearEngine(Engine, NonlinearComposedSystem):
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
    scheme_name = 'non_linear_engine.png'

    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    m_e = Symbol('m_e', positive=True)
    e = Symbol('e', positive=True)
    d = Symbol('d', positive=True)
    l_0 = Symbol('l_0', positive=True)
    z = dynamicsymbols('z')
    phi = dynamicsymbols('\\varphi')
    Omega = Symbol('Omega', positive=True)
    g = Symbol('g', positive=True)

    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 e=None,
                 d=None,
                 l_0=None,
                 z=None,
                 phi=None,
                 ivar=Symbol('t'),
                 Omega=None,
                 g=None,
                 **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if d is not None: self.d = d
        if phi is not None: self.phi = phi
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if l_0 is not None: self.l_0 = l_0
        if z is not None: self.z = z
        if Omega is not None: self.Omega = Omega
        if g is not None: self.g = g
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._engine = FreeEngine(self.M, self.m_e, self.g, z= self.z, phi = self.phi, qs=[self.z])(label='Engine')
        
        self._left_mount  = EngineMount(self.k_m, self.z, qs=[self.z])(label='Left engine mount')
        
        self._right_mount = EngineMount(self.k_m, self.z, qs=[self.z])(label='Right engine mount')
        self._spring_left = EngineMount(2 * self.k_m,
                             pos1=(self.z**2 + self.d**2)**0.5 - self.l_0,
                             qs=[self.z])(label='Nonlinear engine mount left')
        self._spring_right = EngineMount(2 * self.k_m,
                             pos1=(self.z**2 + self.d**2)**0.5 - self.l_0,
                             qs=[self.z])(label='Nonlinear engine mount right')

            

        components['_engine'] = self._engine  
        components['_spring_left'] = self._spring_left
        components['_spring_right'] = self._spring_right
        components['_left_mount'] = self._left_mount
        components['_right_mount'] = self._right_mount


        return components


    
    
# DONE
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
    scheme_name = 'non_linear_engine.png'

    def get_default_data(self):

        m0, k0, e0, l, omega = symbols('m_0 k_0 e_0 l Omega', positive=True)

        default_data_dict = {
            self.M: [
                100 * m0, 300 * m0, 500 * m0, 700 * m0, 900 * m0, 200 * m0,
                400 * m0, 600 * m0, 800 * m0
            ],
            self.m_e: [
                m0, 3 * m0, 5 * m0, 7 * m0, 9 * m0, 2 * m0, 4 * m0, 6 * m0,
                8 * m0
            ],
            self.k_m: [
                k0, 2 * k0, 4 * k0, 6 * k0, 8 * k0, 3 * k0, 5 * k0, 7 * k0,
                9 * k0
            ],
            self.e: [
                2 * e0, S.Half * e0, 4 * e0, S.Half**2 * e0, 3 * e0,
                3 * S.Half * e0, 9 * e0, 3 * S.Half**2 * e0
            ],
            self.l_0: [S.Half * l, l, 2 * l],
            self.d: [4 * l, 8 * l],
            self.beta: [S.Half * pi],
            self.phi: [omega * self.ivar]
        }
        return default_data_dict


#TODO #DDOF # DO POPRAWY !!!
#class NonLinearVeeEnginePerpendicularSprings(Engine):
#    scheme_name = 'nonlin_vee_engine_perpendicular_springs.png'
#    real_name = '440_magnum_v8.jpg'
#
#    M = Symbol('M', positive=True)
#    k_m = Symbol('k_m', positive=True)
#    m_e = Symbol('m_e', positive=True)
#    e = Symbol('e', positive=True)
#    l = Symbol('l', positive=True)
#    x = dynamicsymbols('x')
#    z = dynamicsymbols('z')
#    phi = dynamicsymbols('phi')
#
#    def __init__(self,
#                 M=None,
#                 k_m=None,
#                 m_e=None,
#                 e=None,
#                 l=None,
#                 x=None,
#                 z=None,
#                 phi=None,
#                 ivar=Symbol('t'),
#                 **kwargs):
#
#        if M is not None: self.M = M
#        if k_m is not None: self.k_m = k_m
#        if m_e is not None: self.m_e = m_e
#        if e is not None: self.e = e
#        if l is not None: self.l = l
#        if x is not None: self.x = x
#        if z is not None: self.z = z
#        if phi is not None: self.phi = phi
#
#        self.qs = [self.x, self.z]  #phi nie jest wspolrzedna tego przypadku
#        self.ivar = ivar
#
#        self._init_from_components(**kwargs)
#
#    @property
#    def components(self):
#
#        components = {}
#
#        self._engine_block_horizontal_component = MaterialPoint(
#            self.M, pos1=self.x,
#            qs=[self.x])(label='Engine block (horizontal component)')
#        self._engine_block_vertical_component = MaterialPoint(
#            self.M, pos1=self.z,
#            qs=[self.z])(label='Engine block (vertical component)')
#        self._rotating_unbalanced_mass_horizontal_component = MaterialPoint(
#            self.m_e, pos1=self.x + self.e * sin(self.phi),
#            qs=[self.x
#                ])(label='Rotating unbalanced mass (horizontal component)')
#        self._rotating_unbalanced_mass_vertical_component = MaterialPoint(
#            self.m_e, pos1=self.z + self.e * cos(self.phi),
#            qs=[self.z])(label='Rotating unbalanced mass (vertical component)')
#        self._engine_mount_left_vertical_component = Spring(
#            self.k_m,
#            pos1=(self.z**2 + self.x**2)**(1 / 2) - self.x,
#            qs=[self.z])(label='Engine mount left (vertical component)')
#        self._engine_mount_left_horizontal_component = Spring(
#            self.k_m,
#            pos1=(self.z**2 + self.x**2)**(1 / 2) - self.z,
#            qs=[self.x])(label='Engine mount left (horizontal component)')
#        self._engine_mount_right_vertical_component = Spring(
#            self.k_m,
#            pos1=(self.z**2 + self.x**2)**(1 / 2) - self.x,
#            qs=[self.z])(label='Engine mount right (vertical component)')
#        self._engine_mount_right_horizontal_component = Spring(
#            self.k_m,
#            pos1=(self.z**2 + self.x**2)**(1 / 2) - self.z,
#            qs=[self.x])(label='Engine mount right (horizontal component)')
#
#        components[
#            '_engine_block_horizontal_component'] = self._engine_block_horizontal_component
#        components[
#            '_engine_block_vertical_component'] = self._engine_block_vertical_component
#        components[
#            '_rotating_unbalanced_mass_horizontal_component'] = self._rotating_unbalanced_mass_horizontal_component
#        components[
#            '_rotating_unbalanced_mass_vertical_component'] = self._rotating_unbalanced_mass_vertical_component
#        components[
#            '_engine_mount_left_vertical_component'] = self._engine_mount_left_vertical_component
#        components[
#            '_engine_mount_left_horizontal_component'] = self._engine_mount_left_horizontal_component
#        components[
#            '_engine_mount_right_vertical_component'] = self._engine_mount_right_vertical_component
#        components[
#            '_engine_mount_right_horizontal_component'] = self._engine_mount_right_horizontal_component
#
#        return components
#
#        super().__init__(system, **kwargs)
#
#    def symbols_description(self):
#        self.sym_desc_dict = {
#            self.M: r'Mass of engine block',
#            self.k_m: r'Spring stiffness coefficient',
#            self.m_e: r'Mass in eccentric displacement',
#            self.e: r'Eccentric displacement',
#            #            self.l: r'?????',
#            self.phi: r'Angle of deviation from the vertical',
#            self.x: r'Horizontal displacement',
#            self.z: r'Vertical displacement',
#        }
#
#        return self.sym_desc_dict

# DONE #Amadi Sprawdza
class EngineWithTMD(Engine):
    """
    Model of a DDoF Engine with Tuned Mass Damper attached

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
    detail_scheme_name = 'sruba_pasowana.png'

    M = Symbol('M', positive=True)
    k_m = Symbol('k_m', positive=True)
    k_E = Symbol('k_TMD', positive=True)
    m_e = Symbol('m_e', positive=True)
    m_E = Symbol('m_TMD', positive=True)
    e = Symbol('e', positive=True)
    z = dynamicsymbols('z')
    z_E = dynamicsymbols('z_TMD')
    phi = dynamicsymbols('varphi')
    Omega = Symbol('Omega', positive=True)
    ivar = Symbol('t')
    Re = Symbol('R_e', positive=True)
    g = Symbol('g', positive=True)
    m0 = Symbol('m_0', positive=True)
    k0 = Symbol('k_0', positive=True)
    e0 = Symbol('e_0', positive=True)

    def __init__(self,
                 M=None,
                 k_m=None,
                 k_E=None,
                 m_e=None,
                 m_E=None,
                 e=None,
                 z=None,
                 z_E=None,
                 phi=None,
                 g=None,
                 Omega=None,
                 ivar=Symbol('t'),
                 **kwargs):
        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if k_E is not None: self.k_TMD = k_TMD
        if phi is not None: self.phi = phi
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if m_E is not None: self.m_TMD = m_TMD
        if z_E is not None: self.z_TMD = z_TMD
        if z is not None: self.z = z
        if Omega is not None: self.Omega = Omega
        if g is not None: self.g = g
        self.ivar = ivar
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
        components = {}

        
        
        
        #self._engine = EngineVerticalSpringGravity(
        self._engine = EngineConstantVelocityVerticalSpringGravity(
            M=self.M,
            m_e=self.m_e,
            k_m=self.k_m,
            e=self.e,
            g=self.g,
            phi=self.phi,
            z=self.z)(label='Engine (damped object)')

        self._TMD = TunedMassDamper(self.m_E, self.k_E, self.z_E,self.z)(label='TMD (damper)')
        self._TMD_gravity = GravitationalForce(self.m_E,self.g,self.z_E)(label='Gravitational force for TMD')
        
        components['_engine'] = self._engine
        components['_TMD'] = self._TMD
        components['_TMD_gravity'] = self._TMD_gravity


        return components

    def equilibrium_equation(self, static_disp_dict=None, **kwargs):
        
        if static_disp_dict is None:
            static_disp_dict = {
                self.z: Symbol('z_0', positive=True),
                self.z_E: Symbol('z_{TMD0}', positive=True)
            }

        return super().equilibrium_equation(static_disp_dict=static_disp_dict,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, e0 = self.m0, self.k0, self.e0

        default_data_dict = {
            self.phi: [self.Omega * self.ivar],
            self.M: [S.One *m0 * no for no in range(75, 100)],
            self.k_m: [S.One * k0 * no for no in range(1, 20)],
            self.m_E: [S.One * m0 * 0.1 * no for no in range(75, 125)],
            self.k_E: [S.One * k0 * no/10 for no in range(1, 20)],
            self.m_e: [S.One *m0 * no for no in range(1, 25)],
            self.e: [S.One *e0 * no / 10 for no in range(1, 20)],
        }

        return default_data_dict

    def get_numerical_data(self):

        m0, k0, e0, g, lam = symbols('m_0 k_0 e_0 g lambda', positive=True)
        m0, k0 = self.m0, self.k0
        m0 = 10
        k0 = 10
        default_data_dict = {
            self.M: [no for no in range(75, 100)],
            self.m_e: [no for no in range(1, 10)],
            self.k_m: [1000 * no for no in range(1, 20)],
            self.m_E: [0.1 * no for no in range(75, 125)],
            self.k_E: [100 * no for no in range(1, 20)],
            self.e: [no / 10 for no in range(1, 20)],
            self.phi: [self.Omega * self.ivar],
            self.Omega:[S.One * no * 3.1415 for no in range(1, 10)],
        }

        return default_data_dict
    
    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        parameters_dict[self.m_E] = parameters_dict[self.M]+parameters_dict[self.m_e]/10
        parameters_dict[self.k_E] = parameters_dict[self.k_m]/10

        return parameters_dict


    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.MDoFGeneralSolutionComponent,
        mech_comp.MDoFSteadySolutionComponent,
#         mech_comp.SpringForce,
#         mech_comp.MaxStaticForce,
#         mech_comp.MaxDynamicForce,
#         mech_comp.StaticPinDiameter,
#         mech_comp.DynamicPinDiameter,

        ]

        return comp_list

#     def max_static_force_pin(self):

#         force = self._engine._left_mount.force().subs(self._op_points()[0])

#         return force
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def dynamic_force(self):
        return self.spring_force().expand().doit().n(6)
    
    def spring_force(self):
        
        '''This method evaluates dynamic force in single spring of the engine. It works by substituting coordinates in spring force with steady solutions.'''
        force_km=EngineWithTMD().components['_engine'].components['_left_mount'].force()
        sol_dict=self._fodes_system.steady_solution.as_dict()
        F_km=(force_km).subs(sol_dict).subs(self._given_data)
        
        return F_km
    
    def spring_force_tmd(self):
        
        '''This method evaluates dynamic force in spring between tuned mass damper and engine. It works by substituting coordinates in spring force with steady solutions.'''
        
        force_TMD=EngineWithTMD().components['_TMD'].components['_spring'].force()
        sol_dict=self._fodes_system.steady_solution.as_dict()
        F_ktmd=force_TMD.subs(sol_dict).subs(self._given_data)
        
        return F_ktmd
    
       
    def z_tmd_tune_coefficient(self):
        
        '''This method calculates the tuned mass damper tune coefficient. It assumes that nominator of steady solution is equal to 0 and then solves the equation with respect to k_TMD (stiffness of the tune mass damper spring).'''

        phi=self.phi
        Omega=self.Omega
        t=self.ivar
        k_E=self.k_E
        steady_z_tmd=type(self)()._fodes_system.steady_solution.subs(phi, Omega*t).doit()
        coef_z_tmd=list(fraction(steady_z_tmd[0].coeff(cos(Omega*t))))
        sol_z_tmd=solve(Eq(coef_z_tmd[0], 0), k_E)
        
        return sol_z_tmd[0]
    
class VeeEnginePerpendicularSprings(ComposedSystem):
    scheme_name = 'vee_engine_perpendicular_springs.png'
    real_name = '440_magnum_v8.jpg'
    
    M=Symbol('M', positive=True)
    k_m=Symbol('k_m', positive=True)
    m_e=Symbol('m_e', positive=True)
    e=Symbol('e', positive=True)
    l=Symbol('l',positive=True)
    x=dynamicsymbols('x')
    z=dynamicsymbols('z')
    phi=dynamicsymbols('phi')
    ivar=Symbol('t')

    def __init__(self,
                 M=None,
                 k_m=None,
                 m_e=None,
                 e=None,
                 l=None,
                 x=None,
                 z=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if M is not None: self.M = M
        if k_m is not None: self.k_m = k_m
        if m_e is not None: self.m_e = m_e
        if e is not None: self.e = e
        if l is not None: self.l = l
        if x is not None: self.x = x
        if z is not None: self.z = z
        if phi is not None: self.phi = phi
            
        self.ivar=ivar
        self.qs=[self.x, self.z]
        
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):
            
        components = {}
            
        self._MaterialPoint_1x = MaterialPoint(self.M, pos1=self.x, qs=[self.x])
        self._MaterialPoint_1z = MaterialPoint(self.M, pos1=self.z, qs=[self.z])
        self._MaterialPoint_2x = MaterialPoint(self.m_e,
                                             pos1=self.x + self.e * sin(self.phi),
                                             qs=[self.x])
        self._MaterialPoint_2z = MaterialPoint(self.m_e,
                                             pos1=self.z + self.e * cos(self.phi),
                                             qs=[self.z])
        self._SpringVer = Spring(2 * self.k_m, pos1=self.z, qs=[self.z])
        self._SpringHor = Spring(2 * self.k_m, pos1=self.x, qs=[self.x])
        
        components['material_point_1x'] = self._MaterialPoint_1x
        components['material_point_1z'] = self._MaterialPoint_1z
        components['material_point_2x'] = self._MaterialPoint_2x
        components['material_point_2z'] = self._MaterialPoint_2z
        components['vertical_spring'] = self._SpringVer
        components['horizontal_spring'] = self._SpringHor
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.M: r'Mass of engine block',
            self.k_m: r'Spring stiffness coefficient',
            self.m_e: r'',
            self.e: r'',
        }
        return self.sym_desc_dict


    