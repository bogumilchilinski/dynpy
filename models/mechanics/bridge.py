from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from  ..continuous import ContinuousSystem, PlaneStressProblem
from dynpy.models.mechanics.tmd import TunedMassDamper

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST

from functools import cached_property

# class ComposedSystem(HarmonicOscillator):
#     """Base class for all systems

#     """
#     scheme_name = 'damped_car_new.PNG'
#     real_name = 'car_real.jpg'
#     detail_scheme_name = 'sruba_pasowana.png'
#     detail_real_name = 'buick_regal_3800.jpg'
#     _default_args = ()
#     _default_folder_path = "./dynpy/models/images/"
#     z = dynamicsymbols('z')
#     m0 = Symbol('m_0', positive=True)
#     k0 = Symbol('k_0', positive=True)
#     F0 = Symbol('F_0', positive=True)
#     Omega0 = Symbol('Omega_0', positive=True)

#     @classmethod
#     def _scheme(cls):

#         path = cls._default_folder_path + cls.scheme_name

#         return path

#     @classmethod
#     def _real_example(cls):
#         path = cls._default_folder_path + cls.real_name

#         return path

#     @classmethod
#     def _detail_real(cls):
#         path = cls._default_folder_path + cls.detail_real_name

#         return path

#     @classmethod
#     def _detail_scheme(cls):
#         path = cls._default_folder_path + cls.detail_scheme_name

#         return path

#     def _init_from_components(self,*args,system=None,**kwargs):
        
#         if system is None:
#             composed_system = self._elements_sum
#         else:
#             composed_system = system
            
#         #print('CS',composed_system._components)
#         super().__init__(None,system = composed_system)
        
#         #print('self',self._components)
#         if self._components is None:
#             comps = {}
#         else:
#             comps=self._components

#         self._components = {**comps,**self.components}

#     def __init__(self,
#                  Lagrangian=None,
#                  m0 = None,
#                  qs=None,
#                  forcelist=None,
#                  bodies=None,
#                  frame=None,
#                  hol_coneqs=None,
#                  nonhol_coneqs=None,
#                  label=None,
#                  ivar=None,
#                  evaluate=True,
#                  system=None):


#         if ivar is not None: self.ivar = ivar
#         if m0 is not None: self.m0 = m0

#         if qs is not None:
#             self.qs = qs
#         else:
#             self.qs = [self.z]


#         self._init_from_components(system=system)


#     @property
#     def components(self):

#         components = {}

#         self._material_point = MaterialPoint(self.m0, self.qs[0],
#                                              self.qs)('Material Point')
#         components['_material_point'] = self._material_point

#         return components

#     @property
#     def elements(self):
        
        
#         return {**super().components,**self.components}
    
#     @classmethod
#     def preview(cls, example=False):
#         if example:
#             path = cls._real_example()

#         elif example == 'detail_scheme_name':
#             path = cls._detail_scheme()
#         elif example == 'detail_real_name':
#             path = cls._detail_real()
#         else:
#             path = cls._scheme()
#         print(path)
#         with open(f"{path}", "rb") as image_file:
#             encoded_string = base64.b64encode(image_file.read())
#         image_file.close()

#         return IP.display.Image(base64.b64decode(encoded_string))

#     def get_default_data(self):
#         return None

#     def get_numerical_data(self):
#         return None
    
#     def get_random_parameters(self):

#         default_data_dict = self.get_default_data()

#         if default_data_dict:
#             parameters_dict = {
#                 key: random.choice(items_list)
#                 for key, items_list in default_data_dict.items()
#             }
#         else:
#             parameters_dict = None

#         return parameters_dict

#     def get_numerical_parameters(self):

#         default_data_dict = self.get_numerical_data()

#         if default_data_dict:
#             parameters_dict = {
#                 key: random.choice(items_list)
#                 for key, items_list in default_data_dict.items()
#             }
#         else:
#             parameters_dict = None

#         return parameters_dict
    
#     @property
#     def _report_components(self):

#         comp_list = [
#             mech_comp.TitlePageComponent,
#             mech_comp.SchemeComponent,
#             mech_comp.ExemplaryPictureComponent,
#             mech_comp.KineticEnergyComponent,
#             mech_comp.PotentialEnergyComponent,
#             mech_comp.LagrangianComponent,
#             mech_comp.GoverningEquationComponent,
#             mech_comp.FundamentalMatrixComponent,
#             mech_comp.GeneralSolutionComponent,
#             mech_comp.SteadySolutionComponent,
#         ]

#         return comp_list

    
#     def linearized(self):

#         return type(self).from_system(super().linearized())

# class CompoundSystem(ComposedSystem):

#     z = dynamicsymbols('z')
#     _p = Symbol('p')


#     @property
#     def components(self):

#         components = {}

#         self._material_point = MaterialPoint(self._p, self.qs[0],
#                                              self.qs)('Material Point')
#         components['_material_point'] = self._material_point

#         return components


# class NonlinearComposedSystem(ComposedSystem):

#     def frequency_response_function(self,
#                                     frequency=Symbol('Omega', positive=True),
#                                     amplitude=Symbol('a', positive=True)):

#         omega = ComposedSystem(self.linearized()).natural_frequencies()[0]
#         eps = self.small_parameter()

#         exciting_force = self.external_forces()[0]

#         comps = exciting_force.atoms(sin, cos)
#         exciting_amp = sum([exciting_force.coeff(comp) for comp in comps])
#         inertia = self.inertia_matrix()[0]

#         return amplitude * (-frequency**2 + omega**2) * inertia + S(
#             3) / 4 * eps * amplitude**3 - exciting_amp

#     def amplitude_from_frf(self, amplitude=Symbol('a', positive=True)):

#         return solveset(self.frequency_response_function(), amplitude)

#     @property
#     def _report_components(self):

#         comp_list = [
#             mech_comp.TitlePageComponent,
#             mech_comp.SchemeComponent,
#             mech_comp.ExemplaryPictureComponent,
#             mech_comp.KineticEnergyComponent,
#             mech_comp.PotentialEnergyComponent,
#             mech_comp.LagrangianComponent,
#             mech_comp.LinearizationComponent,
#             mech_comp.GoverningEquationComponent,
#             mech_comp.FundamentalMatrixComponent,
#             mech_comp.GeneralSolutionComponent,
#             mech_comp.SteadySolutionComponent,
#         ]

#         return comp_list


class BeamStructure(MaterialPoint):
    """
    Model of a Beam Structure composed of material point with changing point of mass:

    Creates a material point of an inertial body, after inputing correct values of mass -m and general coordinates, which follows a linear motion.
    """
    scheme_name = 'material_point.png'
    real_name = 'material_point.png'

    def get_default_data(self):
        m0=Symbol('m_0',positive=True)
        default_data_dict={
            self.m:[S.One * m0 * no for no in range(100, 200)],
        }
        return default_data_dict

    def get_numerical_data(self):

        self.default_data_dict={
            self.m:[no * 10000 for no in range(100, 200)],
        }
        return self.default_data_dict

class BeamElasticity(Spring):
    """
    Model of a Beam Elasticity: Creates a singular model, after inputing correct values of stiffeness - k_beam and general coordinate(s),
    which analytically display the dynamics of displacing spring after cummulating PE.
    """
    scheme_name =  'spring.PNG'
    real_name = 'spring_real.png'

    def get_default_data(self):
        E=Symbol('E',positive=True)
        I=Symbol('I',positive=True)
        l=Symbol('l',positive=True)
        default_data_dict={
            self.stiffness:[S.One * 48 * E * I / l**3 * no for no in range(50, 80)],

        }
        return default_data_dict

    def get_numerical_data(self):

        self.default_data_dict={
            self.stiffness:[no * 1000 for no in range(500, 800)],
        }
        return self.default_data_dict

# class BridgeExcitation(Force):
#     """
#     Creates enforcement.
#     """
#     scheme_name = ''
#     real_name = ''

#     def get_default_data(self):
#         F0=Symbol('F0',positive=True)
#         Omega0=Symbol('Omega0',positive=True)
#         default_data_dict={
#             self.force: [S.One * no * F0 * 1000 for no in range(10, 25)],
#         }
#         return default_data_dict

#     def get_numerical_data(self):

#         self.default_data_dict={
#             self.force:[no * 1000 for no in range(10, 25)],
#             self.Omega: [2 * 3.14 * 100],
#         }
#         return self.default_data_dict

#Amadi
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
     #test

    """
    scheme_name = 'beam_bridge.PNG'
    real_name = 'beam_bridge_real.PNG'

    m = Symbol('m', positive=True)
    k_beam = Symbol('k_beam', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    z = dynamicsymbols('z')


    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if z is not None: self.z = z
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._mass = BeamStructure(self.m, self.z, qs=[self.z])(label='Mass of a beam')
        self._spring = BeamElasticity(self.k_beam, self.z, qs=[self.z])(label='Beam stiffness')
        self._gravity_force = GravitationalForce(self.m, self.g,self.z)(label='Gravity field')
        self._force = Force(-self.F * sin(self.Omega * self.ivar),pos1=self.z)(label='External force')

        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravity_force'] = self._gravity_force
        components['_force'] = self._force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'Gravitational field acceleration',
            self.Omega: r'Excitation frequency',
            self.F: r'Force acting on a bridge',
            self.l: r'Lenght of a beam',
            self.module: r'Youngs modulus',
            self.inertia: r'Interia of a beam'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        
        m0, k0, E, I, l, E0, I0= symbols('m_0 k_0 E I l E_0 I_0', positive=True)
        
        F0=Symbol('F_0',positive=True)
        Omega0=Symbol('Omega_0',positive=True)
        E0=Symbol('E_0',positive=True)
        I0=Symbol('I_0',positive=True)
        l0=Symbol('l_0',positive=True)
        module=Symbol('E',positive=True)
        inertia=Symbol('I',positive=True)
        l=Symbol('l',positive=True)
        
        default_data_dict={
            l: [S.One * no * l0 for no in range(1, 4)],
            #module: [S.One * E0 * no for no in range(10, 20)],
            #inertia: [S.One * no * I0 for no in range(10, 20)],
            self.F: [S.One * no * F0  for no in range(10, 25)],
            #self.m :[S.One * no * m0 for no in range(10, 20)]
#            self.Omega: [S.One * Omega0]
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict={
            self.F:[no * 1000 for no in range(10, 25)],
            self.Omega: [3.14 * no for no in range(1,10)],
            
        }
        return default_data_dict

    
    def dynamic_force(self):
        data=self._given_data
        amps=self._fodes_system.steady_solution.as_dict()
        dyn_force=(self.components['_spring'].force().subs(amps)).subs(data).expand().doit()
        
        return dyn_force
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    
#Amadi
class BeamBridgeDamped(BeamBridge):
    """Ready to use model of damped bridge represented by the mass supported by elastic beam.
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
    scheme_name = 'bridge_dmp.png'
    real_name = 'beam_bridge_real.PNG'

    m = Symbol('m', positive=True)
    k_beam = Symbol('k_beam', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    z = dynamicsymbols('z')
    c = Symbol('c', positive=True)


    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 c=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if z is not None: self.z = z
        if c is not None: self.c = c
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._beam_bridge = BeamBridge(self.m, self.k_beam, self.ivar, self.g, self.Omega, self.F)(label='Beam Bridge')
        self._damper = Damper(self.c, pos1=self.z, qs=[self.z])(label='Damping of a beam')

        components['_beambridge'] = self._beam_bridge
        components['_damper'] = self._damper

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'Gravitational field acceleration',
            self.Omega: r'Excitation frequency',
            self.F: r'Force acting on a bridge',
            self.l: r'Lenght of a beam',
            self.module: r'Youngs modulus',
            self.inertia: r'Interia of a beam',
            self.c: r'Damping coefficient'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        lam=Symbol('lambda', positive = True)
        E, I, l=symbols('E I l',positive=True)
        F0=Symbol('F_0',positive=True)
        m0=Symbol('m_0',positive=True)
        k0=Symbol('k_0',positive=True)
        
        default_data_dict = {
#             self.k_beam [S.One * 48 * E * I / l**3],
            self.m: [S.One * m0 * no for no in range(100, 150)],
            self.c: [S.One * self.k_beam * lam],
            self.k_beam: [S.One * k0 * no for no in range(50, 100)],
            self.F: [S.One * F0 * no for no in range(10, 25)],
#             lam: [S.One * no * lam0 /1000 for no in range(1, 5)],
        }
        return default_data_dict
    

    def get_numerical_data(self):

        default_data_dict = {
            self.c:[S.One * no for no in range(5,20)],
        }
        return default_data_dict
    
    def dynamic_force(self):
        data=self._given_data
        amps=self._fodes_system.steady_solution.as_dict()
        dyn_force=(self._beam_bridge.components['_spring'].force().subs(amps)).subs(data).expand().doit()

        return dyn_force
#Amadi
class BeamBridgeTMD(BeamBridge):

    scheme_name = 'bridge_tmd.png'
    real_name = 'beam_bridge_real.PNG'
    detail_real_name = 'real_bidge_tmd.jpg'

    m=Symbol('m', positive=True)
    m_TMD=Symbol('m_TMD', positive=True)
    k_beam=Symbol('k_beam', positive=True)
    k_TMD=Symbol('k_TMD', positive=True)
    g=Symbol('g', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    z=dynamicsymbols('z')
    z_TMD=dynamicsymbols('z_TMD')

    
    def __init__(self,
                 m=None,
                 m_TMD=None,
                 k_beam=None,
                 k_TMD=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 z=None,
                 z_TMD=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if m_TMD is not None: self.m_TMD = m_TMD
        if k_TMD is not None: self.k_TMD = k_TMD
        if z_TMD is not None: self.z_TMD = z_TMD
        if z is not None: self.z = z
        self.ivar=ivar
        
        self._init_from_components(**kwargs)

    @cached_property
    def components(self):

        components = {}

        self._beam_bridge = BeamBridge(self.m, self.k_beam, self.ivar, self.g, self.Omega, self.F)(label='Beam Bridge')
        self._TMD = TunedMassDamper(self.m_TMD, self.k_TMD, self.z_TMD,self.z)(label='TMD (damper)')
        self._TMD_gravity = GravitationalForce(self.m_TMD, self.g,self.z_TMD)(label='Gravity field for TMD')
        
        components['_beam_bridge'] = self._beam_bridge
        components['_TMD'] = self._TMD
        components['_TMD_gravity'] = self._TMD_gravity

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
        self.m: r'Mass of beam',
        self.k_beam: r'Beam stiffness',
        self.g: r'Gravitational field acceleration',
        self.Omega: r'Excitation frequency',
        self.F: r'Force acting on a bridge',
        self.l: r'Lenght of a beam',
        self.module: r'Youngs modulus',
        self.inertia: r'Interia of a beam',
        self.m_TMD: r'Mass of TMD',
        self.k_TMD: r'Spring stiffness of a TMD'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, k0, E, I, l, E0, I0= symbols('m_0 k_0 E I l E_0 I_0', positive=True)

        default_data_dict = {
            self.m: [S.One * m0  * no for no in range(100, 200)],
            self.m_TMD: [S.One * m0/10  * no for no in range(100, 200)],
            self.k_TMD: [S.One * 48/10 * E * I / l**3 * no for no in range(50, 80)],
            self.k_beam: [S.One * 48 * E * I / l**3 * no for no in range(50, 80)],
#             E: [S.One * no * E0 for no in range(10,50)],
#             I: [S.One * no * I0 for no in range(20,80)],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        parameters_dict[self.m_TMD] = parameters_dict[self.m]/10
        parameters_dict[self.k_TMD] = parameters_dict[self.k_beam]/10

        return parameters_dict
    
    
    def get_numerical_data(self):

        default_data_dict = {
            self.m_TMD: [no * 1000 for no in range(100, 200)],
            self.k_TMD: [no * 100 for no in range(500, 800)],
            self.g: [9.81],

        }
        return default_data_dict

    def dynamic_force(self):
        
        '''This method counts dynamic force in in the beam of the investigated bridge without TMD. It works by substituting coordinates in beam force with steady solutions.'''
        
        force=BeamBridgeTMD().components['_beam_bridge'].components['_spring'].force()
        sol_dict=self._fodes_system.steady_solution.as_dict()
        F_kbeam=force.subs(sol_dict).subs(self._given_data)
        
        return F_kbeam.subs(self._given_data)
    
    def static_force(self):
        data=self._given_data
        ans=self.dynamic_force()
        free_coeff=ans.subs({cos(self.Omega*self.ivar):0, sin(self.Omega*self.ivar):0}).subs(data)
        return (free_coeff)
    
    def dynamic_force_TMD(self):
        '''This method counts dynamic force in spring between tuned mass damper and beam. It works by substituting coordinates in spring force with steady solutions.'''
        force=BeamBridgeTMD().components['_TMD'].components['_spring'].force()
        sol_dict=self._fodes_system.steady_solution.as_dict()
        F_tmd=force.subs(sol_dict).subs(self._given_data)
        
        return F_tmd
    
    def tmd_tune_coefficient(self):
        '''This method calculates the TMD tune coefficient. It assumes that nominator of steady solution is equal to 0 and then solves the equation with respect to k_TMD.'''
        
        pure_system = type(self)()
        
        Omega=pure_system.Omega
        t=self.ivar
        k_TMD=pure_system.k_TMD
        steady_z_tmd=pure_system._fodes_system.steady_solution.doit()
        coef_z_tmd=list(fraction(steady_z_tmd[0].coeff(sin(Omega*t))))
        sol_z_tmd=solve(Eq(coef_z_tmd[0], 0), k_TMD)
        
        return sol_z_tmd[0].subs(self._given_data)
    
    @property
    def _report_components(self):
        
        comp_list=[
        *REPORT_COMPONENTS_LIST,
#         mech_comp.FundamentalMatrixComponent,
#         mech_comp.MDoFGeneralSolutionComponent,
#         mech_comp.MDoFSteadySolutionComponent,
        mech_comp.TMDForceComponent

        ]

        return comp_list
    
class BeamBridgeDampedTMD(BeamBridge):

    scheme_name = 'bridge_tmd_dmp.png'
    real_name = 'beam_bridge_real.PNG'
    detail_real_name = 'real_bidge_tmd.jpg'

    c_TMD=Symbol('c_TMD', positive=True)
    c=Symbol('c', positive=True)
    m=Symbol('m', positive=True)
    m_TMD=Symbol('m_TMD', positive=True)
    k_beam=Symbol('k_beam', positive=True)
    k_TMD=Symbol('k_TMD', positive=True)
    g=Symbol('g', positive=True)
    Omega=Symbol('Omega', positive=True)
    F=Symbol('F', positive=True)
    z=dynamicsymbols('z')
    z_TMD=dynamicsymbols('z_TMD')

    def __init__(self,
                 m=None,
                 m_TMD=None,
                 k_beam=None,
                 k_TMD=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 c_TMD=None,
                 c=None,
                 z_TMD=None,
                 z=None,
                 **kwargs):
        
        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if m_TMD is not None: self.m_TMD = m_TMD
        if k_TMD is not None: self.k_TMD = k_TMD
        if z_TMD is not None: self.z_TMD = z_TMD
        if z is not None: self.z = z
        if c is not None: self.c = c
        if c_TMD is not None: self.c_TMD = c_TMD
        self.ivar=ivar
        
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        
        self._beam_bridge_damped = BeamBridgeDamped(self.m, self.k_beam, self.ivar, self.g, self.Omega, self.F, self.c)(label='Beam Bridge Damped')
        self._TMD = TunedMassDamper(self.m_TMD, self.k_TMD, self.z_TMD,self.z)(label='TMD (damper)')
        self._TMD_gravity = GravitationalForce(self.m_TMD, self.g,self.z_TMD)(label='Gravity field for TMD')
        self._damper_TMD = Damper(self.c_TMD, pos1=self.z - self.z_TMD, qs=[self.z, self.z_TMD])(label='Damper of a TMD')
        
        components['_beam_bridge_damped'] = self._beam_bridge_damped
        components['_TMD'] = self._TMD
        components['_TMD_gravity'] = self._TMD_gravity
        components['_damper_TMD'] = self._damper_TMD

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
        self.m: r'Mass of beam',
        self.k_beam: r'Beam stiffness',
        self.g: r'Gravitational field acceleration',
        self.Omega: r'Excitation frequency',
        self.F: r'Force acting on a bridge',
        self.l: r'Lenght of a beam',
        self.module: r'Youngs modulus',
        self.inertia: r'Interia of a beam',
        self.m_TMD: r'Mass of TMD',
        self.k_TMD: r'Spring stiffness of a TMD',
        self.c: r'Damping coefficient of a bridge',
        self.c_TMD: r'Damping coefficient of a TMD'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        lam=Symbol('lambda', positive = True)
        m0=Symbol('m_0', positive = True)
        k0=Symbol('k_0', positive = True)
        E, I = symbols('E I', positive=True)
        default_data_dict = {
            self.c_TMD: [S.One * self.k_TMD * lam],
            self.m_TMD: [S.One * no * m0 * 1000 for no in range(20, 80)],
            self.k_TMD: [S.One * no * (48 * E * I /self.l**3)* 1000 for no in range(10, 20)]
        }

        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.c_TMD:[S.One * no for no in range(5,10)],
            self.m_TMD: [S.One * no * 1000 for no in range(20, 80)],
            self.k_TMD: [S.One * no * 1000 for no in range(15, 20)]
        }
        return default_data_dict