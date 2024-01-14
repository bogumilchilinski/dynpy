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

    def __init__(self,
                 R_1=None,
                 R_2=None,
                 C=None,
                 U=None,
                 q_1=None,
                 q_2=None,
                 t=None,
                 U_th=None,
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
        
        t=self.ivar
        
        exprI=0.084*Heaviside(t-100)+0*Heaviside(t-400)-0*Heaviside(t-600)-0*Heaviside(t-800)
        #plot(exprI,(t,0,1001))
        funcI = lambdify(t, exprI)
        
        
        #tworzenie odesystemu z drugiego wzoru z publikacji
        uth_eq=Eq(diff(U_th,t),((U_th)/(R_th*C_th))+I_li/C_th)
        display(uth_eq)
        ode_th = ODESystem(odes=uth_eq.rhs+uth_eq.lhs, dvars=U_th, ode_order=1)
        analytical_sol1=ode_th.solution
        #display(analytical_sol1)
        
        
        napiecie_eq=Eq(U_li,U_oc-R_0*I_li-U_th)
        display(napiecie_eq)
        
        
        calka=integrate(I_li, (t, t_0, t))
        SOC_eq_bezcalki = Eq(SOC,((1)/(C_rated)))
        SOC_eq = Eq(SOC,SOC_init+((1)/(C_rated))*calka)
        display(SOC_eq)
        
        
        #pojemnosc w amperogodzinach

        data_dict1 = {
            U_oc:4.25,
            R_0:0.02423,
            R_th:0.1,
            C_th:100,
            SOC_init:0.99,
            C_rated:100,
            t_0:0,
        }

        U_oc_dict = {
            '100':4.25,
            '90':4.25,
            '80':4.25,
            '70':4.25,
            '60':3.80,
            '50':3.70,
            '40':3.60,
            '30':3.50,
            '20':3.40,
            '10':3.3,
            '0':3.2
        }

        n_points=2000 #rozdzielczosc symulacji
        t_array=np.linspace(0,1001,n_points)#czas symulacji i prad symulacji


        prad_list=[]
        u0_list=[]
        for i in t_array:
            prad_list.append(funcI(i))
        display(prad_list)

        u0_list = [i * data_dict1[R_0] for i in prad_list]


        #obliczanie rownania rozniczkowego
        ode_th.subs(I_li,exprI).subs(data_dict1).numerized(backend='numpy').compute_solution(t_array,[0])


        #tworzenie dataframe z wynikami
        tabelka=ode_th.subs(I_li,exprI).subs(data_dict1).numerized(backend='numpy').compute_solution(t_array,[0])


        wynsym = list(tabelka.iloc[:,1])
        U_th_list=wynsym


        u_array=[]
        soc_array=[]

        for i in range(n_points):
            soc_array.append(float((SOC_eq_bezcalki.rhs.subs(data_dict1)*(-1)*(np.trapz(prad_list[0:i], x=t_array[0:i]))/3600+SOC_init.subs(data_dict1))))
        #tutaj jest dopasowanie ocv w zaleznosci od soc, zamiast ifow
            soc_calk = int(soc_array[i]*10)
            soc_calk = soc_calk*10
            data_dict1["U_oc"]=U_oc_dict[str(soc_calk)]


            u_array.append(data_dict1["U_oc"]-U_th_list[i]-u0_list[i])



        df_u = pd.DataFrame({t : t_array})
        df_u[U_li] = u_array


        df_I = pd.DataFrame({t : t_array})
        df_I[I_li] = prad_list

        df_SOC = pd.DataFrame({t : t_array})
        df_SOC[SOC] = soc_array


        tdf_u=TimeDataFrame(df_u).set_index(t)
        tdf_I=TimeDataFrame(df_I).set_index(t)
        tdf_SOC=TimeDataFrame(df_SOC).set_index(t)