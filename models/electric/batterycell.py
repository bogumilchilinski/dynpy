import sympy
import sympy as sym
import pandas
from sympy import *
#sprawdzenie poprawności importów
init_printing()
sympy.Symbol, sym.Symbol, Symbol
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
#from dynpy.models.electric import batterycell 
#mechanics_printing(pretty_print=True)




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
    #U_th = Symbol('U_th')
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
    
    
    def voltage_ode(self):

        t = self.ivar
        U_th = self.U_th
        R_th = self.R_th
        C_th = self.C_th
        I_li = self.I_li
        
        
        uth_eq=Eq(diff(U_th,t, evaluate=False),((U_th)/(R_th*C_th))+I_li/C_th)
        ode_th = ODESystem(odes=uth_eq.rhs+uth_eq.lhs, dvars=U_th, ode_order=1)
        
        return ode_th
    
    
    def _uth_simulation(self,data_dict1,U_oc_dict,R_0_dict,step_time,step_current):
        
        t = self.ivar
        U_th = self.U_th
        R_th = self.R_th
        C_th = self.C_th
        I_li = self.I_li
        U_li = self.U_li
        U_oc = self.U_oc
        R_0 = self.R_0
        t_0 = self.t_0
        SOC = self.SOC
        C_rated = self.C_rated
        SOC_init = self.SOC_init
        ode_th = self.voltage_ode()
        
        

     
        #wzor na napiecie
        napiecie_eq=Eq(U_li,U_oc-R_0*I_li-U_th)

        
        #wzor na soc
        calka=integrate(I_li, (t, t_0, t))
        SOC_eq_bezcalki = Eq(SOC,((1)/(C_rated)))
        SOC_eq = Eq(SOC,SOC_init+((1)/(C_rated))*calka)
        
        
        #dane symulacji
        n_points=step_time[-1]*4 #rozdzielczosc symulacji probki na sekunde
        t_array=np.linspace(0,step_time[-1]*1.1,n_points)#czas symulacji i prad symulacji

        #tworzenie pustych list
        prad_list=[]
        u0_list=[]
        u_array=[]
        soc_array=[]
        
        #tworzenie wymuszenia
        exprI=step_current[0]*Heaviside(t-step_time[0])+step_current[1]*Heaviside(t-step_time[1])+step_current[2]*Heaviside(t-step_time[2])+step_current[3]*Heaviside(t-step_time[3])
        funcI = lambdify(t, exprI)
        

        #tworzenie dataframe z wynikami rownania rozniczkowego+tworzenie z tego listy
        tabelka=ode_th.subs(I_li,exprI).subs(data_dict1).numerized(backend='numpy').compute_solution(t_array,[0])
        U_th_list = list(tabelka.iloc[:,1])




        for i in range(n_points):
            prad_list.append(funcI(t_array[i]))
            soc_array.append(float((SOC_eq_bezcalki.rhs.subs(data_dict1)*(-1)*(np.trapz(prad_list[0:i], x=t_array[0:i]))/3600+SOC_init.subs(data_dict1))))

            #tutaj jest dopasowanie ocv w zaleznosci od soc, zamiast ifow
            soc_calk = int(soc_array[i]*10)
            soc_calk = soc_calk*10
            data_dict1["U_oc"]=U_oc_dict[str(soc_calk)]
            data_dict1["R_0"]=R_0_dict[str(soc_calk)]

            u0_list.append( prad_list[-1] * data_dict1[R_0])
            u_array.append(data_dict1["U_oc"]-U_th_list[i]-u0_list[i])


        wartosci = [u_array,prad_list,soc_array,t_array]
#komentuje zeby sprawdzic czy dobrze licze
        return wartosci
    
        #return prad_list



    def VoltageResponse(self,data_dict1,U_oc_dict,R_0_dict,step_time,step_current):

        t=self.t
        U_li= self.U_li
        wyn = self._uth_simulation(data_dict1,U_oc_dict,R_0_dict,step_time,step_current)
        t_array = wyn[3]
        u_array = wyn[0]

        df_u = pd.DataFrame({t : t_array})
        df_u[U_li] = u_array
        df_u = df_u.set_index(t)

        return df_u
    
    def SocLevelResponse(self,data_dict1,U_oc_dict,R_0_dict,step_time,step_current):
        
        t=self.t
        SOC= self.SOC
        wyn = self._uth_simulation(data_dict1,U_oc_dict,R_0_dict,step_time,step_current)
        t_array = wyn[3]
        soc_array = wyn[2]

        df_SOC = pd.DataFrame({t : t_array})
        df_SOC[SOC] = soc_array
        df_SOC = df_SOC.set_index(t)

        return df_SOC
    
    
    def CurrentForcing(self,data_dict1,U_oc_dict,R_0_dict,step_time,step_current):
        
        t=self.t
        I_li= self.I_li
        wyn = self._uth_simulation(data_dict1,U_oc_dict,R_0_dict,step_time,step_current)
        t_array = wyn[3]
        prad_list = wyn[1]

        df_I = pd.DataFrame({t : t_array})
        df_I[I_li] = prad_list
        df_I = df_I.set_index(t)

        return df_I
        
        

#         tdf_u=TimeDataFrame(df_u).set_index(t)
#         tdf_I=TimeDataFrame(df_I).set_index(t)
#         tdf_SOC=TimeDataFrame(df_SOC).set_index(t)