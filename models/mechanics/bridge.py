from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from  ..continuous import ContinuousSystem, PlaneStressProblem

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
    detail_scheme_name = 'sruba_pasowana.png'
    detail_real_name = 'buick_regal_3800.jpg'
    _default_args = ()
    _default_folder_path = "./dynpy/models/images/"

    
    z = dynamicsymbols('z')
    
    m0 = Symbol('m_0', positive=True)
    k0 = Symbol('k_0', positive=True)
    F0 = Symbol('F_0', positive=True)
    Omega0 = Symbol('Omega_0', positive=True)

    @classmethod
    def _scheme(cls):

        path = cls._default_folder_path + cls.scheme_name

        return path

    @classmethod
    def _real_example(cls):
        path = cls._default_folder_path + cls.real_name

        return path

    @classmethod
    def _detail_real(cls):
        path = cls._default_folder_path + cls.detail_real_name

        return path

    @classmethod
    def _detail_scheme(cls):
        path = cls._default_folder_path + cls.detail_scheme_name

        return path


    def _init_from_components(self,*args,system=None,**kwargs):
        
        if system is None:
            composed_system = self._elements_sum
        else:
            composed_system = system
            
        #print('CS',composed_system._components)
        super().__init__(None,system = composed_system)
        
        #print('self',self._components)
        if self._components is None:
            comps = {}
        else:
            comps=self._components

        self._components = {**comps,**self.components}

    def __init__(self,
                 Lagrangian=None,
                 m0 = None,
                 qs=None,
                 forcelist=None,
                 bodies=None,
                 frame=None,
                 hol_coneqs=None,
                 nonhol_coneqs=None,
                 label=None,
                 ivar=None,
                 evaluate=True,
                 system=None):


        if ivar is not None: self.ivar = ivar
        if m0 is not None: self.m0 = m0

        if qs is not None:
            self.qs = qs
        else:
            self.qs = [self.z]


        self._init_from_components(system=system)


    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self.m0, self.qs[0],
                                             self.qs)('Material Point')
        components['_material_point'] = self._material_point

        return components

    @property
    def elements(self):
        
        
        return {**super().components,**self.components}



    
    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        elif example == 'detail_scheme_name':
            path = cls._detail_scheme()
        elif example == 'detail_real_name':
            path = cls._detail_real()
        else:
            path = cls._scheme()
        print(path)
        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))

    def get_default_data(self):
        return None

    def get_numerical_data(self):
        return None
    
    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict

    def get_numerical_parameters(self):

        default_data_dict = self.get_numerical_data()

        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict = None

        return parameters_dict
    
    
    @property
    def _report_components(self):

        comp_list = [
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
        ]

        return comp_list


    
    
    def linearized(self):

        return type(self).from_system(super().linearized())

#     def tensioner_belt_force(self):
#         return self.k_tensioner * self.steady_solution()

#     def left_belt_force(self):
#         return self.k_belt * self.steady_solution()

#     def right_belt_force(self):
#         return self.k_belt * self.steady_solution()

#     def max_static_force_pin(self):
#         return abs(self.static_load().doit()[0])

#     def max_dynamic_force_pin(self):
#         return self.frequency_response_function() * self.stiffness_matrix(
#         )[0] + self.max_static_force_pin()

#     def static_force_pin_diameter(self):
#         kt = Symbol('k_t', positive=True)
#         Re = Symbol('R_e', positive=True)
#         return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

#     def dynamic_force_pin_diameter(self):
#         kt = Symbol('k_t', positive=True)
#         Re = Symbol('R_e', positive=True)
#         return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


class CompoundSystem(ComposedSystem):

    z = dynamicsymbols('z')
    _p = Symbol('p')


    @property
    def components(self):

        components = {}

        self._material_point = MaterialPoint(self._p, self.qs[0],
                                             self.qs)('Material Point')
        components['_material_point'] = self._material_point

        return components


    

class NonlinearComposedSystem(ComposedSystem):

    def frequency_response_function(self,
                                    frequency=Symbol('Omega', positive=True),
                                    amplitude=Symbol('a', positive=True)):

        omega = ComposedSystem(self.linearized()).natural_frequencies()[0]
        eps = self.small_parameter()

        exciting_force = self.external_forces()[0]

        comps = exciting_force.atoms(sin, cos)
        exciting_amp = sum([exciting_force.coeff(comp) for comp in comps])
        inertia = self.inertia_matrix()[0]

        return amplitude * (-frequency**2 + omega**2) * inertia + S(
            3) / 4 * eps * amplitude**3 - exciting_amp

    def amplitude_from_frf(self, amplitude=Symbol('a', positive=True)):

        return solveset(self.frequency_response_function(), amplitude)

    @property
    def _report_components(self):

        comp_list = [
            mech_comp.TitlePageComponent,
            mech_comp.SchemeComponent,
            mech_comp.ExemplaryPictureComponent,
            mech_comp.KineticEnergyComponent,
            mech_comp.PotentialEnergyComponent,
            mech_comp.LagrangianComponent,
            mech_comp.LinearizationComponent,
            mech_comp.GoverningEquationComponent,
            mech_comp.FundamentalMatrixComponent,
            mech_comp.GeneralSolutionComponent,
            mech_comp.SteadySolutionComponent,
        ]

        return comp_list

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


    """
    scheme_name = 'beam_bridge.PNG'
    real_name = 'beam_bridge_real.PNG'

    m = Symbol('m', positive=True)
    k_beam = Symbol('k_beam', positive=True)
    g = Symbol('g', positive=True)
    Omega = Symbol('Omega', positive=True)
    F = Symbol('F', positive=True)
    l = Symbol('l', positive=True)
    module = Symbol('E', positive=True)
    inertia = Symbol('I', positive=True)
    z = dynamicsymbols('z')
    E0 = Symbol('E_0', positive=True)
    I0 = Symbol('I_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    m0 = Symbol('m_0', positive=True)
    lam0 = Symbol('lam_0', positive=True)
    F0 = Symbol('F_0', positive=True)

    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 l=None,
                 module=None,
                 inertia=None,
                 z=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if l is not None: self.l = l
        if z is not None: self.z = z
        if module is not None: self.module = module
        if inertia is not None: self.inertia = inertia
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z,
                                   qs=[self.z])(label='Material point')
        self._spring = Spring(self.k_beam, self.z, qs=[self.z])(label='Beam')
        self._gravity_force = GravitationalForce(self.m, self.g,
                                                 self.z)(label='Gravity field')
        self._force = Force(-self.F * sin(self.Omega * self.ivar),
                            pos1=self.z)(label='External force')

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

        #E0, I0, l0, m0, k0,c0, lam0= symbols('E_0 I_0 l_0 m_0 k_0 c_0 lambda_0', positive=True)
        E0, I0, l0, m0, lam0, F0 = self.E0, self.I0, self.l0, self.m0, self.lam0, self.F0
        default_data_dict = {

            #self.lam:[10],
            #self.c:[self.k_beam*self.lam],
            self.m: [S.One * no * m0 * 10 for no in range(2, 8)],
            self.module: [S.One * E0 * no for no in range(10, 20)],
            self.inertia: [S.One * no * I0 for no in range(10, 20)],
            self.l: [S.One * no * l0 for no in range(10, 25)],
            #self.lam:[lam0,2*lam0,3*lam0,4*am0,5*lam0,6*lam0,7*lam0,8*lam0,9*lam0],
            self.F: [S.One * no * F0 for no in range(10, 25)],
            self.k_beam: [S.One * 48 * self.module * self.inertia / self.l**3]
        }
        return default_data_dict


#Amadi
class BeamBridgeDamped(ComposedSystem):
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
    l = Symbol('l', positive=True)
    module = Symbol('E', positive=True)
    inertia = Symbol('I', positive=True)
    z = dynamicsymbols('z')
    c = Symbol('c', positive=True)
    E0 = Symbol('E_0', positive=True)
    I0 = Symbol('I_0', positive=True)
    l0 = Symbol('l_0', positive=True)
    m0 = Symbol('m_0', positive=True)
    lam0 = Symbol('lambda_0', positive=True)
    F0 = Symbol('F_0', positive=True)
    lam = Symbol('lambda', positive=True)

    def __init__(self,
                 m=None,
                 k_beam=None,
                 ivar=Symbol('t'),
                 g=None,
                 Omega=None,
                 F=None,
                 l=None,
                 module=None,
                 inertia=None,
                 z=None,
                 c=None,
                 lam=None,
                 **kwargs):

        if m is not None: self.m = m
        if k_beam is not None: self.k_beam = k_beam
        if g is not None: self.g = g
        if Omega is not None: self.Omega = Omega
        if F is not None: self.F = F
        if l is not None: self.l = l
        if z is not None: self.z = z
        if module is not None: self.module = module
        if inertia is not None: self.inertia = inertia
        if c is not None: self.c = c
        if lam is not None: self.lam = lam
        self.ivar = ivar

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z,
                                   qs=[self.z])(label='Material point')
        self._spring = Spring(self.k_beam, self.z,
                              qs=[self.z])(label='Beam stiffness')
        self._gravity_force = GravitationalForce(self.m, self.g,
                                                 self.z)(label='Gravity field')
        self._force = Force(-self.F * sin(self.Omega * self.ivar),
                            pos1=self.z)(label='External force')
        self._damper = Damper(self.c, pos1=self.z,
                              qs=[self.z])(label='Damping of a beam')

        components['_mass'] = self._mass
        components['_spring'] = self._spring
        components['_gravity_force'] = self._gravity_force
        components['_force'] = self._force
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

        #E0, I0, l0, m0, k0,c0, lam0= symbols('E_0 I_0 l_0 m_0 k_0 c_0 lambda_0', positive=True)
        E0, I0, l0, m0, lam0, F0 = self.E0, self.I0, self.l0, self.m0, self.lam0, self.F0
        default_data_dict = {

            #self.lam:[10,
            self.m: [S.One * no * m0 * 10 for no in range(2, 8)],
            self.module: [S.One * E0 * no for no in range(10, 20)],
            self.inertia: [S.One * no * I0 for no in range(10, 20)],
            self.l: [S.One * no * l0 for no in range(10, 25)],
            self.lam: [S.One * no * lam0 for no in range(1, 5)],
            self.F: [S.One * no * F0 for no in range(10, 25)],
            self.k_beam: [S.One * 48 * self.module * self.inertia / self.l**3],
            self.c: [self.k_beam * self.lam]
        }
        return default_data_dict

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
                 F=Symbol('F', positive=True),
                 z=dynamicsymbols('z'),
                 z_TMD=dynamicsymbols('z_TMD'),
                 **kwargs):

        self.m = m
        self.k_beam = k_beam
        self.g = g
        self.Omega = Omega
        self.F = F
        self.m_TMD = m_TMD
        self.k_TMD = k_TMD
        self.z_TMD = z_TMD
        self.z = z

        self.mass = MaterialPoint(m, z, qs=[z])
        self.spring = Spring(k_beam, z, qs=[z])
        self.gravity_force = GravitationalForce(self.m, self.g, z)
        self.gravity_TMD = GravitationalForce(self.m_TMD, self.g, z_TMD)
        self.force = Force(F * sin(Omega * ivar), pos1=z)
        self.TMD = MaterialPoint(m_TMD, pos1=z_TMD, qs=[z_TMD])
        self.spring_TMD = Spring(k_TMD, z, z_TMD, qs=[z, z_TMD])
        composed_system = (self.mass + self.spring + self.gravity_force + self.force +
                  self.TMD + self.spring_TMD )

        super().__init__(composed_system,**kwargs)



    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k_beam: r'Beam stiffness',
            self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):

        E, I, l, m0, k0,F0 = symbols('E I l_beam m_0 k_0 F_0', positive=True)

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
            self.F:[F0*no for no in range(0,10)]
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
