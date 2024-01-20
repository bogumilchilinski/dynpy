import sympy
import pandas
from sympy import *
from scipy import integrate
import numpy as np
from sympy import lambdify
#sprawdzenie poprawności importów
init_printing()


from pylatex import Document,Section
from dynpy.utilities import *
from dynpy.solvers.linear import ODESystem
import numpy as np
from dynpy.dynamics import *
from sympy import *
from sympy.physics.mechanics import *
from dynpy.models.mechanics.trolley import *
from dynpy.models.mechanics.principles import *
from dynpy.models.electric import *
import pandas as pd
from dynpy.utilities.report import ReportText, Markdown, Picture, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, Markdown, SympyFormula, TikZPlot, CurrentContainer, AutoMarker
from dynpy.utilities.report import SymbolsDescription, DescriptionsRegistry
from dynpy.utilities.adaptable import *
from dynpy.models.electric.elements import *
#mechanics_printing(pretty_print=True)


from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)
import numpy as np
from numpy import *

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, CombustionEngine

from .elements import Resistor, Inductor, Capacitor, VoltageSource
from ..mechanics.principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST

from sympy import symbols, Eq, diff

class BatteryCell(ComposedSystem):
    scheme_name = 'thevenincircuit.PNG'
    real_name = 'liioncell.PNG'

    R_1=Symbol('R_1', positive=True)
    R_2=Symbol('R_2', positive=True)
    C=Symbol('C', positive=True)
    U=Symbol('U', positive=True)
    q_1=dynamicsymbols('q_1')
    q_2=dynamicsymbols('q_2')
    t = Symbol('t')
    U_li = Function('U_li')(t)
    U_oc = Symbol('U_oc')
    R_0 = Symbol('R_0')
    I_li = Function('I_zad')(t)
    R_th = Symbol('R_th')
    C_th = Symbol('C_th')
    SOC = Function('SOC')(t)
    SOC_init = Symbol('SOC_init')
    C_rated = Symbol('C_rated')
    t_0 = Symbol('t_0')
    U_th = Function('U_th')(t)
    funcI= Function('funcI')(t)
    U_th = Symbol('U_th')
    R_th = Symbol('R_th')
    C_th = Symbol('C_th')
    I_li = Symbol('I_li')

    def __init__(self,
                 R_1=None,
                 R_2=None,
                 C=None,
                 U=None,
                 q_1=None,
                 q_2=None,
                 t=None,
                 U_th=None,
                 R_th=None,
                 C_th=None,
                 I_li=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if t is not None: self.t = t
        if R_1 is not None: self.R_1 = R_1
        if R_2 is not None: self.R_2 = R_2
        if C is not None: self.C= C
        if U is not None: self.U = U
        if q_1 is not None: self.q_1 = q_1
        if q_2 is not None: self.q_2 = q_2
        if U_th is not None: self.U_th = U_th
        if R_th is not None: self.R_th = R_th
        if C_th is not None: self.C_th = C_th
        if I_li is not None: self.I_li = I_li
        self.qs = [self.q_1, self.q_2]
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._resistor_1 = Resistor(self.R_1, self.q_1, qs=[self.q_1])
        self._resistor_2 = Resistor(self.R_2, q0=self.q_2, qs=[self.q_2])
        self._voltagesource = VoltageSource(self.U, q0 = self.q_1, qs=[self.q_1])
        self._capacitor = Capacitor(self.C, q0 = self.q_1-self.q_2, qs=[self.q_1, self.q_2])



        components['resistor_1'] = self._resistor_1
        components['resistor_2'] = self._resistor_2
        components['voltagesource'] = self._voltagesource
        components['capacitor'] = self._capacitor

        return components
    
    
    def voltage_ode(self):

        t = self.ivar
        U_th = self.U_th
        R_th = self.R_th
        C_th = self.C_th
        I_li = self.I_li

        
#         t =Symbol('t')
#         U_th=Symbol('U_th')
#         R_th=Symbol('R_th')
#         C_th=Symbol('C_th')
#         I_li=Symbol('I_li')
        
        
        exprI=0.084*Heaviside(t-100)+0*Heaviside(t-400)-0*Heaviside(t-600)-0*Heaviside(t-800)
        #plot(exprI,(t,0,1001))
        funcI = lambdify(t, exprI)
        
        #Derivative(Function('U')(t),t, evaluate=False))
        #tworzenie odesystemu z drugiego wzoru z publikacji
        #uth_eq=Eq(Derivative(U_th,t, evaluate=False),((U_th)/(R_th*C_th))+I_li/C_th) to nie działa
        
        uth_eq=Eq(diff(U_th,t, evaluate=False),((U_th)/(R_th*C_th))+I_li/C_th)
        #display(uth_eq)
        
        return uth_eq