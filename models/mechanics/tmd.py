from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from  ..continuous import ContinuousSystem, PlaneStressProblem
from .trolley import SpringMassSystem

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem


    
class TMD(SpringMassSystem):
    """
    Ready to use sample Single Degree of Freedom System of Tuned Mass Damper
    

    Arguments:
    =========
        m = Mass
            -Mass of TMD

        k = Spring coefficient
            -Spring carrying the TMD

        ivar = symbol object
            -Independant time variable

        qs = dynamicsymbol object
            -Generalized coordinates

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

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
    scheme_name = 'tmd_engine_vertical_spring_nogravity.png'
    real_name = 'tmd_engine_real.jpg'

    m_E=Symbol('m_E', positive=True)
    k_E=Symbol('k_E', positive=True)
    z=dynamicsymbols('z')
    z_E=dynamicsymbols('z_E')
    
    def __init__(self,
                 m_E=None,
                 k_E=None,
                 z_E=None,
                 z=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m_E is not None: self.m_E = m_E
        if k_E is not None: self.k_E = k_E
        if z is not None: self.z = z
        if z_E is not None: self.z_E = z_E
        self.ivar = ivar
        
   
        self.qs = [self.z_E]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self._material_point = MaterialPoint(self.m_E, self.z_E, qs=self.qs)(label='Material point - mass of the TMD')
        self._spring = Spring(self.k_E, pos1=self.z_E, pos2=self.z, qs=self.qs)(label='Spring - stiffness of the spring')
        
        components['_material_point'] = self._material_point
        components['_spring'] = self._spring
        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of system on the spring',
            self.k: r'Spring coefficient ',
        }

        return self.sym_desc_dict

