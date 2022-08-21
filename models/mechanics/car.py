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

    def tensioner_belt_force(self):
        return self.k_tensioner * self.steady_solution()

    def left_belt_force(self):
        return self.k_belt * self.steady_solution()

    def right_belt_force(self):
        return self.k_belt * self.steady_solution()

    def max_static_force_pin(self):
        return abs(self.static_load().doit()[0])

    def max_dynamic_force_pin(self):
        return self.frequency_response_function() * self.stiffness_matrix(
        )[0] + self.max_static_force_pin()

    def static_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_static_force_pin()) / (pi * kt * Re))**(1 / 2)

    def dynamic_force_pin_diameter(self):
        kt = Symbol('k_t', positive=True)
        Re = Symbol('R_e', positive=True)
        return ((4 * self.max_dynamic_force_pin()) / (pi * kt * Re))**(1 / 2)


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
    

class CarMovementConstantThrottle(ComposedSystem):

    x = dynamicsymbols('x')
    e = dynamicsymbols('e')

    F_t = Symbol('F_t', positive=True)
    F_p= Symbol('F_p', positive=True)
    F_w = Symbol('F_w', positive=True)
    F_b = Symbol('F_b', positive=True)
    F_n = Symbol('F_n', positive=True)

    T = Symbol('T', positive=True)
    i= Symbol('i', positive=True)
    i_c = Symbol('i_c', positive=True)
    eta = Symbol('eta', positive=True)
    r_d = Symbol('r_d', positive=True)

    m= Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    f_0 = Symbol('f_0', positive=True)
    kappa = Symbol('kappa', positive=True)

    A = Symbol('A', positive=True)
    c_x = Symbol('c_x', positive=True)
    rho = Symbol('rho', positive=True)

    m_red = Symbol('m_red', positive=True)
    J_k= Symbol('J_k', positive=True)
    J_s = Symbol('J_s', positive=True)
    alpha = Symbol('alpha', positive=True)

    a_h = Symbol('a_h', positive=True)
    u = Symbol('u', positive=True)
    
    n = Symbol('n', positive=True)
    
    P = Symbol('P', positive=True)
    I = Symbol('I', positive=True)
    D = Symbol('D', positive=True)
    
    n_min = Symbol('n_min', positive=True)
    n_max = Symbol('n_max', positive=True)
    
    def __init__(self,
                x = None,
                T = None,
                i = None,
                i_c = None,
                eta = None,
                r_d = None,
                m = None,
                g = None,
                f_0 = None,
                kappa = None,
                A = None,
                c_x = None,
                rho = None,
                m_red = None,
                J_k= None,
                J_s = None,
                alpha = None,
                a_h = None,
                u = None,
                n = None,
                P = None,
                I = None,
                D = None,
                e = None,
                n_min = None,
                n_max = None,
                ivar=None,
                **kwargs):

        if x is not None: self.x = x
        if T is not None: self.T = T
        if i is not None: self.i = i
        if eta is not None: self.eta = eta
        if i_c is not None: self.i_c = i_c
        if r_d is not None: self.r_d = r_d
        if m is not None: self.m = m
        if g is not None: self.g = g
        if f_0 is not None: self.f_0 = f_0
        if kappa is not None: self.kappa = kappa
        if A is not None: self.A = A
        if c_x is not None: self.c_x = c_x
        if rho is not None: self.rho = rho
        if m_red is not None: self.m_red = m_red
        if J_k is not None: self.J_k = J_k
        if J_s is not None: self.J_s = J_s
        if alpha is not None: self.alpha = alpha
        if a_h is not None: self.a_h = a_h
        if u is not None: self.u = u
        if n is not None: self.n = n
        if P is not None: self.P = P
        if I is not None: self.I = I
        if D is not None: self.D = D
        if e is not None: self.e = e
        if n_min is not None: self.n_min = n_min
        if n_max is not None: self.n_max = n_max
        if ivar is not None: self.ivar = ivar
            
        self.qs = [self.x]
        
        #self.m_red = (self.m + (4*self.J_k+self.J_s*self.i_c**2)/self.r_d**2)
        
        self.c = S.One/2*self.A*self.c_x*self.rho*diff(self.x,self.ivar)*2/3
        self.mu = self.m_red*self.g*self.f_0*self.kappa
        self.f = self.m_red*self.g*self.f_0
        self.F = (self.T*self.i*self.i_c*self.eta)/self.r_d
        self.n = (diff(self.x,self.ivar)*self.i*self.i_c*60 )/(2*np.pi*self.r_d)
        
        self.error = -7.21914962310842e-24*self.ivar**10 + 4.09285356849447e-20*self.ivar**9 - 9.80641183997517e-17*self.ivar**8 + 1.29051452578393e-13*self.ivar**7 - 1.01472654714565e-10*self.ivar**6 + 4.85741460232082e-8*self.ivar**5 - 1.37892965900947e-5*self.ivar**4 +             0.00213941785802763*self.ivar**3 - 0.148427605950237*self.ivar**2 + 7.13838237925456*self.ivar - 67.8916383683618 - self.x + 67

        self._init_from_components(**kwargs)

        
    @property
    def components(self):

        components = {}

        v = self.x.diff(self.ivar)
        omega = v

        self._car_model = MaterialPoint(self.m_red,
                                     pos1=self.x,
                                     qs=self.qs)(label='car model')

        self._drag = Damper(self.c,
                         pos1=self.x,
                         qs=self.qs)(label='drag')

        self._friction1 = Damper(self.mu,
                         pos1=self.x,
                         qs=self.qs)(label='friction1')

        self._friction2 = Force(-self.f,
                         pos1=self.x,
                         qs=self.qs)(label='friction2')
        
        self._throttle = Force(0.5*self.F+0.5*self.F*sin(3.14*omega/(15)),
                         pos1=self.x,
                         qs=self.qs)(label='throttle')

        components['_car_model'] = self._car_model
        components['_drag'] = self._drag
        components['_friction1'] = self._friction1
        components['_friction2'] = self._friction2
        components['_throttle'] = self._throttle
        return components
    
    @property
    def current_gear(self):
        return 1.024

    def symbols_description(self):
        self.sym_desc_dict = {
            self.T: r'Moment obrotowy',
            self.i: r'Przełożenie obecnego biegu',
            self.eta: r'Sprawność przekładni',
            self.i_c: r'Przełożenie przekładni głównej',
            self.r_d: r'Promień dynamiczny',
            self.m: r'Masa pojazdu',
            self.g: r'Przyśpieszenie ziemskie',
            self.f_0: r'Współczynnik oporu toczenia',
            self.kappa: r'Współczynnik proporcjonalności współczynnika oporu toczenia zależnego od prędkości',
            self.A: r'Powierzchnia czołowa pojazdu',
            self.c_x: r'Współczynnik oporu powietrza',
            self.rho: r'Gęstość powietrza',
            self.m_red: r'Masa redukowana',
            self.J_k: r'Moment bezwładności koła',
            self.J_s: r'Moment bezwładności mas wirujących silnika',
            self.alpha: r'Kąt nachylenia drogi',
            self.a_h: r'Podstawowe opóźnienie siły hamowania',
            self.u: r'Siła hamowania',
            self.x: r'Przemieszczenie w osi x',
            self.ivar: r'Czas',
        }

        return self.sym_desc_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.T: [140],
            self.i: [self.current_gear],
            self.eta: [0.95],
            self.i_c: [4.105],
            self.r_d: [0.97*(0.235*0.45+8.5*0.0254)],
            self.m: [1760],
            self.g: [9.81],
            self.f_0: [0.01],
            self.kappa: [0.00001],
            self.A: [2.14],
            self.c_x: [0.28],
            self.rho: [1.2],
            self.m_red: [1831.6],
            self.J_k: [0.9],
            self.J_s: [0.2],
            self.alpha: [0],
            self.a_h: [4],
            self.P: [350], #2100
            self.I: [10],
            self.D: [20], #100
        }

        return default_data_dict
    
    
    def get_default_data(self):

        default_data_dict = {
            self.T: [140],
            self.i: [self.current_gear],
            self.eta: [0.95],
            self.i_c: [4.105],
            self.r_d: [0.97*(0.235*0.45+8.5*0.0254)],
            self.m: [1760],
            self.g: [9.81],
            self.f_0: [0.01],
            self.kappa: [0.00001],
            self.A: [2.14],
            self.c_x: [0.28],
            self.rho: [1.2],
            self.m_red: [1831.6],
            self.J_k: [0.9],
            self.J_s: [0.2],
            self.alpha: [0],
            self.a_h: [4],
            self.P: [210],
            self.I: [0],
            self.D: [100],
        }

        return default_data_dict
    
    
class SecondGearCarMovement(CarMovementConstantThrottle):
    @property
    def current_gear(self):
        return 2.08
    
class ThirdGearCarMovement(CarMovementConstantThrottle):
    @property
    def current_gear(self):
        return 1.361
class FourthGearCarMovement(CarMovementConstantThrottle):
    @property
    def current_gear(self):
        return 0.84

class FifthGearCarMovement(CarMovementConstantThrottle):
    @property
    def current_gear(self):
        return 0.686
