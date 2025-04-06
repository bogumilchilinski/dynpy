from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset, integrate, lambdify, Heaviside, integrate, exp, atan, sign, sign)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, CombustionEngine


import pandas as pd

from .elements import *

import base64
import random
import IPython as IP
import numpy as np
import inspect


from ..mechanics.trolley import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin

t = Symbol('t')

#from dynpy.models.electric import batterycell 
#mechanics_printing(pretty_print=True)




class BatteryCell(ComposedSystem):
    
    
    scheme_name = 'thevenincircuit.png'
    real_name = 'liioncell.PNG'

    R_1=Symbol('R_1', positive=True)
    R_2=Symbol('R_2', positive=True)

    L_1=Symbol('L_1', positive=True)
    L_2=Symbol('L_2', positive=True)
    
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
        
        self._coil_1 = Inductor(self.L_1, self.q_1, qs=[self.q_1])
        self._coil_2 = Inductor(self.L_2, self.q_2, qs=[self.q_2])
        
        self._resistor_1 = Resistor(self.R_1, self.q_1, qs=[self.q_1])
        self._resistor_2 = Resistor(self.R_2, q0=self.q_2, qs=[self.q_2])
        self._voltagesource = VoltageSource(self.U, q0 = self.q_1, qs=[self.q_1])
        self._capacitor = Capacitor(self.C, q0 = self.q_1-self.q_2, qs=[self.q_1, self.q_2])

        
        components['coil_1'] = self._coil_1
        components['coil_2'] = self._coil_2
        
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
        
        
        
        
        uth_eq=Eq(diff(U_th,t, evaluate=False),((U_th)/(R_th*C_th))+I_li/C_th)
        ode_th = ODESystem(odes=uth_eq.rhs+uth_eq.lhs, dvars=U_th, ode_order=1)
        
        return ode_th
    
    
    def uth_simulation(self,data_dict1,celldata, step_time,step_current):
        
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
        t_array=np.linspace(0,step_time[-1]*1,n_points)#czas symulacji i prad symulacji

        #tworzenie pustych list
        prad_list=[]
        u0_list=[]
        u_array=[]
        soc_array=[]
        
        #tworzenie wymuszenia
        exprI=0*Heaviside(t-0)
        for i in range(len(step_time)):
            exprI=exprI + step_current[i]*Heaviside(t-step_time[i])
        funcI = lambdify(t, exprI)
        #exprI=step_current[0]*Heaviside(t-step_time[0])+step_current[1]*Heaviside(t-step_time[1])+step_current[2]*Heaviside(t-step_time[2])+step_current[3]*Heaviside(t-step_time[3])
        #funcI = lambdify(t, exprI)
        

        #tworzenie dataframe z wynikami rownania rozniczkowego+tworzenie z tego listy
        tabelka=ode_th.subs(I_li,exprI).subs(data_dict1).numerized(backend='numpy').compute_solution(t_array,[0])
        
        #display(ode_th)
        #display(tabelka)
        
        U_th_list = list(tabelka.iloc[:,1])




        for i in range(n_points):
            prad_list.append(funcI(t_array[i]))
            soc_array.append(float((SOC_eq_bezcalki.rhs.subs(data_dict1)*(-1)*(np.trapz(prad_list[0:i], x=t_array[0:i]))/3600+SOC_init.subs(data_dict1))))

            #tutaj jest dopasowanie ocv w zaleznosci od soc, zamiast ifow
            soc_calk = int(soc_array[i]*100)
            #soc_calk = soc_calk*100
            data_dict1[self.U_oc]= celldata.at[soc_calk, 'OCV']
            data_dict1[self.R_0]=celldata.at[soc_calk, 'R0']

            u0_list.append( prad_list[-1] * data_dict1[self.R_0])
            u_array.append(data_dict1[self.U_oc]-U_th_list[i]-u0_list[i])


        wartosci = [u_array,prad_list,soc_array,t_array]
        
        return wartosci
    

    def _voltage_response(self,data_dict1,celldata,step_time,step_current):

        t=self.t
        U_li= self.U_li
        wyn = self.uth_simulation(data_dict1,celldata,step_time,step_current)
        t_array = wyn[3]
        u_array = wyn[0]

        df_u = pd.DataFrame({t : t_array})
        df_u[U_li] = u_array
        df_u = df_u.set_index(t)

        return df_u
    
    def _soc_level_response(self,data_dict1,celldata,step_time,step_current):
        
        t=self.t
        SOC= self.SOC
        wyn = self.uth_simulation(data_dict1,celldata,step_time,step_current)
        t_array = wyn[3]
        soc_array = wyn[2]

        df_SOC = pd.DataFrame({t : t_array})
        df_SOC[SOC] = soc_array
        df_SOC = df_SOC.set_index(t)

        return df_SOC
    
    
    def _current_forcing(self,data_dict1,celldata,step_time,step_current):
        
        t=self.t
        I_li= self.I_li
        wyn = self.uth_simulation(data_dict1,celldata,step_time,step_current)
        t_array = wyn[3]
        prad_list = wyn[1]

        df_I = pd.DataFrame({t : t_array})
        df_I[I_li] = prad_list
        df_I = df_I.set_index(t)

        return df_I
        
        

#         tdf_u=TimeDataFrame(df_u).set_index(t)
#         tdf_I=TimeDataFrame(df_I).set_index(t)
#         tdf_SOC=TimeDataFrame(df_SOC).set_index(t)


class BatteryCharging(ComposedSystem):
    
    V_OCV = Symbol('V_OCV',positive=True)
    Vd = Symbol('Vd',positive=True)
    Vo = Symbol('Vo',positive=True)
    Ro = Symbol('Ro',positive=True)
    current = Symbol('I_ch', positive=True)
    t0 = Symbol('t_0',positive=True)
    resistance_0 = Symbol('R_0', positive=True)
    resistance_coefficient_alpha = Symbol('R_Alpha', positive=True) #symbol alfy
    resistance_coefficient_beta = Symbol('R_Beta', positive=True) #symbol bety
    resistance_coefficient_gamma = Symbol('R_Gamma', positive=True) #symbol gammy
    resistance_coefficient_blocking = Symbol('R_Omega', positive=True) #symbol omegi
    charger_blocking_resistance = Symbol('R_b', positive=True) 
    resistance  = Symbol('R' , positive= True)
    coil_coefficient = Symbol('L',positive=True)
    capacitor_coefficient= Symbol('C',positive=True)
    electric_charge = dynamicsymbols('q')
    k_p = Symbol ('k_p',positive=True)
    k_d = Symbol ('k_d',positive=True)
    PD = Symbol ('PD',positive=True)
    q_char = Symbol ('q_char',positive=True)

    def __init__(self,
                 V_OCV = V_OCV,
                 Vo = Vo,
                 Vd = Vd,
                 Ro = Ro,
                 current = current,
                 t0 = t0,
                 resistance_0 = resistance_0,
                 resistance_coefficient_alpha = resistance_coefficient_alpha,
                 resistance_coefficient_beta = resistance_coefficient_beta,
                 resistance_coefficient_gamma = resistance_coefficient_gamma,
                 resistance_coefficient_blocking = resistance_coefficient_blocking,
                 charger_blocking_resistance = charger_blocking_resistance,
                 resistance = resistance,
                 coil_coefficient = coil_coefficient,
                 capacitor_coefficient=capacitor_coefficient,
                 electric_charge = electric_charge,
                 k_p = k_p,
                 k_d = k_d,
                 PD = PD,
                 q_char = q_char,
                 qs=[electric_charge],
                 ivar=Symbol('t'),
                 **kwargs):
    
        t=ivar
        self.resistance_0 = resistance_0
        self.Vo = Vo
        self.Vd = Vd
        self.Ro = Ro
        self.resistance_coefficient_alpha =  resistance_coefficient_alpha
        self.resistance_coefficient_beta =  resistance_coefficient_beta
        self.resistance_coefficient_gamma =  resistance_coefficient_gamma
        self.resistance_coefficient_blocking =  resistance_coefficient_blocking
        self.coil_coefficient = coil_coefficient
        self.capacitor_coefficient = capacitor_coefficient
        self.electric_charge = electric_charge
        self.current = current
        self.qs = [self.electric_charge]
        
        ### krzywa U(SOC) w procesie ładowania
        x_char = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,     2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6650, 6800, 6850, 6900, 6950, 7000, 7050, 7100, 7150, 7200, 7250, 7300, 7350, 7400, 7450, 7500, 7550, 7600, 7650, 7700, 7750, 7800, 7850, 7900, 7950, 8000, 8050, 8100, 8150, 8200, 8250, 8300, 8350, 8400, 8450, 8500, 8550, 8600, 8650, 8800, 8900, 9400, 9800, 10200, 10600]

        y_char = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.55, 2.62, 2.74, 2.86, 2.98, 3.1, 3.22, 3.34, 3.46, 3.58, 3.7, 3.75, 3.801, 3.802, 3.803, 3.804, 3.805, 3.806, 3.807, 3.808, 3.809, 3.81, 3.811, 3.812, 3.813, 3.814, 3.815, 3.816, 3.817, 3.818, 3.819, 3.82, 3.821, 3.822, 3.823, 3.824, 3.825, 3.826, 3.827, 3.828, 3.829, 3.83, 3.831, 3.832, 3.833, 3.834, 3.835, 3.836, 3.837, 3.838, 3.839, 3.84, 3.841, 3.842, 3.843, 3.844, 3.845, 3.846, 3.847, 3.848, 3.849, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4, 4.01, 4.02, 4.05, 4.1, 4.15, 4.2, 4.22, 4.24, 4.26, 4.28, 4.3, 4.32, 4.34, 4.36, 4.38, 4.4, 4.42, 4.44, 4.46, 4.48, 4.5, 4.52, 4.6, 4.8, 5, 5.5, 6, 6.5, 7]
        char = np.polyfit(x_char,y_char,11)
        
        self.V_OCV  = char[0]*self.electric_charge**11 + char[1]*self.electric_charge**10 + char[2]*self.electric_charge**9 + char[3]*self.electric_charge**8 + char[4]*self.electric_charge**7 + char[5]*self.electric_charge**6 + char[6]*self.electric_charge**5 + char[7]*self.electric_charge**4 + char[8]*self.electric_charge**3 + char[9]*self.electric_charge**2 + char[10]*self.electric_charge + char[11]  
        
        
        self.charger_blocking_resistance = exp((self.electric_charge + (-1 * self.q_char - 250))/50) #self.resistance_coefficient_blocking
        
        self.resistance = self.resistance_0 + exp((self.resistance_coefficient_alpha - self.electric_charge)/self.resistance_coefficient_beta) +  exp((self.electric_charge + self.resistance_coefficient_gamma)/5550 ) 
        
        self.material_point_L = MaterialPoint(self.coil_coefficient,self.electric_charge, qs=[self.electric_charge])
        self.damper_R= Damper(self.resistance,self.electric_charge, qs=[self.electric_charge])
        self.damper_R_blocking= Damper(self.charger_blocking_resistance,self.electric_charge, qs=[self.electric_charge])
        self.force_OCV = Force(self.V_OCV,self.electric_charge,qs=[self.electric_charge])     
        self.spring_C = Spring(1/self.capacitor_coefficient, self.electric_charge, qs=[self.electric_charge])
        
        
        system = self.damper_R + self.damper_R_blocking + self.force_OCV + self.material_point_L + self.spring_C
        #super().__init__(system(qs),**kwargs)

        ###### regulator PD
        
        self.q_char = q_char
        self.t0 = ((self.q_char*0.0002777)/self.current)*3600 # ładunek żądany/stała(zamiana C na AH) dzielone przez prad ladowania  i zamiana na sekundy
        self.charger_slope = 0.001
        self.k_p = k_p 
        self.k_d = k_d   
        
        pd = Force((self.k_p* (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) *(self.current - self.electric_charge.diff(t)), self.electric_charge, qs = [self.electric_charge]  ) + Force((self.k_d * (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) * (self.current - self.electric_charge.diff(t)).diff(t),self.electric_charge ,qs = [self.electric_charge])

        system_pid = system+pd(qs)
        super().__init__(system_pid,**kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance_0: r'Rezystancja podstawowa',
            self.resistance: r'Rezystor (SOC)',
            self.resistance_coefficient_alpha: r'Współczynnik rezystancji /alpha',
            self.resistance_coefficient_beta: r'Współczynnik rezystancji /beta',
            self.resistance_coefficient_gamma: r'Współczynnik rezystancji /gamma',
            self.charger_blocking_resistance: r'Rezystancja blokująca',
            self.resistance_coefficient_blocking: r'Współczynnik rezystancji blokującej', 
            self.capacitor_coefficient: r'Kondensator',
            self.coil_coefficient: r'Cewka indukcyjna',
            self.V_OCV: r'Napięcie',
            self.Vo:r'Napięcie',
            self.electric_charge: r'Ładunek',
            self.electric_charge.diff(t): r'Prąd',
            self.current: r'Prąd ładowania',
            self.t0: r'Czas ładowania',
            self.pd: r'Regulator PD',
            self.k_p: r'Człon proporcjonalny regulatora',
            self.k_d: r'Człon różniczkujący regulatora',
            self.ivar: r'Czas',

        }
        return self.sym_desc_dict


class BatteryModeling(ComposedSystem):
    
    V_OCV = Symbol('V_OCV',positive=True)
    current = Symbol('I_ch', positive=True)
    t0 = Symbol('t_0',positive=True)
    resistance_0 = Symbol('R_0', positive=True)
    resistance_coefficient_alpha = Symbol('R_Alpha', positive=True) #symbol alfy
    resistance_coefficient_beta = Symbol('R_Beta', positive=True) #symbol bety
    resistance_coefficient_gamma = Symbol('R_Gamma', positive=True) #symbol gammy
    resistance_coefficient_blocking = Symbol('R_Omega', positive=True) #symbol omegi
    charger_blocking_resistance = Symbol('R_b', positive=True) 
    resistance  = Symbol('R' , positive= True)
    coil_coefficient = Symbol('L',positive=True)
    capacitor_coefficient= Symbol('C',positive=True)
    electric_charge = dynamicsymbols('q')
    pid = Symbol ('PD',positive=True)
    q_char = Symbol ('q_char',positive=True)
    k_p = Symbol ('k_p',positive=True)
    k_d = Symbol ('k_d',positive=True)
    PD = Symbol ('PD',positive=True)
    n_0 = Symbol('n_0',positive=True)
    n_1 = Symbol('n_0',positive=True)
    n_2 = Symbol('n_0',positive=True)
    n_3 = Symbol('n_0',positive=True)
    n_4 = Symbol('n_0',positive=True)
    n_5 = Symbol('n_0',positive=True)
    n_6 = Symbol('n_0',positive=True)
    n_7 = Symbol('n_0',positive=True)
    n_8 = Symbol('n_0',positive=True)
    n_9= Symbol('n_0',positive=True)
    n_10 = Symbol('n_0',positive=True)
    n_11 = Symbol('n_0',positive=True)
    charger_slope = Symbol('chs',positive=True)

    def __init__(self,
                 V_OCV = V_OCV,
                 current = current,
                 t0 = t0,
                 resistance_0 = resistance_0,
                 resistance_coefficient_alpha = resistance_coefficient_alpha,
                 resistance_coefficient_beta = resistance_coefficient_beta,
                 resistance_coefficient_gamma = resistance_coefficient_gamma,
                 resistance_coefficient_blocking = resistance_coefficient_blocking,
                 charger_blocking_resistance = charger_blocking_resistance,
                 resistance = resistance,
                 coil_coefficient = coil_coefficient,
                 capacitor_coefficient=capacitor_coefficient,
                 electric_charge = electric_charge,
                 pid = pid,
                 q_char = q_char,
                 k_p = k_p,
                 k_d = k_d,
                 PD = PD,
                 n_0 = n_0,
                 n_1 = n_1,
                 n_2 = n_2,
                 n_3 = n_3,
                 n_4 = n_4,
                 n_5 = n_5,
                 n_6 = n_6,
                 n_7 = n_7,
                 n_8 = n_8,
                 n_9 = n_9,
                 n_10 = n_10,
                 n_11 = n_11,
                 charger_slope = charger_slope,
                 qs=[electric_charge],
                 ivar=Symbol('t'),
                 **kwargs):
    
        t=ivar
        self.resistance_0 = resistance_0
        self.resistance_coefficient_alpha =  resistance_coefficient_alpha
        self.resistance_coefficient_beta =  resistance_coefficient_beta
        self.resistance_coefficient_gamma =  resistance_coefficient_gamma
        self.resistance_coefficient_blocking =  resistance_coefficient_blocking
        self.coil_coefficient = coil_coefficient
        self.capacitor_coefficient = capacitor_coefficient
        self.electric_charge = electric_charge
        self.current = current
        self.qs = [self.electric_charge]
        
        ### krzywa U(SOC) w procesie ładowania
#         x_char = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,     2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6650, 6800, 6850, 6900, 6950, 7000, 7050, 7100, 7150, 7200, 7250, 7300, 7350, 7400, 7450, 7500, 7550, 7600, 7650, 7700, 7750, 7800, 7850, 7900, 7950, 8000, 8050, 8100, 8150, 8200, 8250, 8300, 8350, 8400, 8450, 8500, 8550, 8600, 8650, 8800, 8900, 9400, 9800, 10200, 10600]

#         y_char = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.55, 2.62, 2.74, 2.86, 2.98, 3.1, 3.22, 3.34, 3.46, 3.58, 3.7, 3.75, 3.801, 3.802, 3.803, 3.804, 3.805, 3.806, 3.807, 3.808, 3.809, 3.81, 3.811, 3.812, 3.813, 3.814, 3.815, 3.816, 3.817, 3.818, 3.819, 3.82, 3.821, 3.822, 3.823, 3.824, 3.825, 3.826, 3.827, 3.828, 3.829, 3.83, 3.831, 3.832, 3.833, 3.834, 3.835, 3.836, 3.837, 3.838, 3.839, 3.84, 3.841, 3.842, 3.843, 3.844, 3.845, 3.846, 3.847, 3.848, 3.849, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4, 4.01, 4.02, 4.05, 4.1, 4.15, 4.2, 4.22, 4.24, 4.26, 4.28, 4.3, 4.32, 4.34, 4.36, 4.38, 4.4, 4.42, 4.44, 4.46, 4.48, 4.5, 4.52, 4.6, 4.8, 5, 5.5, 6, 6.5, 7]
#         char = np.polyfit(x_char,y_char,11)
        
        self.V_OCV  = n_0*self.electric_charge**11 + n_1*self.electric_charge**10 + n_2*self.electric_charge**9 + n_3*self.electric_charge**8 + n_4*self.electric_charge**7 + n_5*self.electric_charge**6 + n_6*self.electric_charge**5 + n_7*self.electric_charge**4 + n_8*self.electric_charge**3 + n_9*self.electric_charge**2 + n_10*self.electric_charge + n_11  
        
        
        self.charger_blocking_resistance = exp((self.electric_charge + (-1 * self.q_char - 250))/50) #self.resistance_coefficient_blocking
        self.resistance = self.resistance_0 + exp((self.resistance_coefficient_alpha - self.electric_charge)/self.resistance_coefficient_beta) +  exp((self.electric_charge + self.resistance_coefficient_gamma)/self.resistance_coefficient_beta) 
        self.material_point_L = MaterialPoint(self.coil_coefficient,self.electric_charge, qs=[self.electric_charge])
        self.damper_R= Damper(self.resistance,self.electric_charge, qs=[self.electric_charge])
        self.force_OCV = Force(self.V_OCV,self.electric_charge,qs=[self.electric_charge])     
        self.spring_C = Spring(1/self.capacitor_coefficient, self.electric_charge, qs=[self.electric_charge])
        
        
        system =self.damper_R  +  self.force_OCV + self.material_point_L + self.spring_C 
        #super().__init__(system(qs),**kwargs)

        ###### regulator PD
        self.damper_Rb= Damper(self.charger_blocking_resistance,self.electric_charge, qs=[self.electric_charge])
        self.q_char = q_char
        self.t0 = ((self.q_char)/self.current) # ładunek żądany/stała(zamiana C na AH) dzielone przez prad ladowania  i zamiana na sekundy
        self.charger_slope = charger_slope
        self.k_p = k_p 
        self.k_d = k_d
        
        self.pid = ((self.k_p* (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) *(self.current - self.electric_charge.diff(t))) + ((self.k_d * (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) * (self.current - self.electric_charge.diff(t)).diff(t))
        
        pd = Force((self.k_p* (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) *(self.current - self.electric_charge.diff(t)), self.electric_charge, qs = [self.electric_charge]  ) + Force((self.k_d * (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) * (self.current - self.electric_charge.diff(t)).diff(t),self.electric_charge ,qs = [self.electric_charge]) + self.damper_Rb
        
        system_pid = system+pd(qs)
        super().__init__(system_pid,**kwargs)
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance_0: r'Rezystancja podstawowa',
            self.resistance: r'Rezystor (SOC)',
            self.resistance_coefficient_alpha: r'Współczynnik rezystancji \\alpha',
            self.resistance_coefficient_beta: r'Współczynnik rezystancji \\beta',
            self.resistance_coefficient_gamma: r'Współczynnik rezystancji \\gamma',
            self.resistance_coefficient_blocking: r'Współczynnik rezystancji blokującej', 
            self.capacitor_coefficient: r'Kondensator',
            self.coil_coefficient: r'Cewka indukcyjna',
            self.V_OCV: r'Napięcie',
            self.electric_charge: r'Ładunek',
            self.electric_charge.diff(t): r'Prąd',
            self.current: r'Prąd ładowania',
            self.ivar: r'czas',
            
        }
        return self.sym_desc_dict
    
    
    
class CircutRl(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RL system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    inductance=Symbol('L', positive=True)
    q0=dynamicsymbols('q_c')
    qs=dynamicsymbols('qs')
    frame=Symbol('frame', positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 inductance=None,
                 ivar=None,
                 q0=None,
                 qs=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if inductance is not None: self.inductance = inductance
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if q0 is not None: self.q0 = q0
        if qs is not None: self.qs = qs
        if frame is not None: self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self.resistor = Resistor(self.resistance, self.q0,  qs=self.qs,ivar=self.ivar , frame=base_frame)('resistor')
        self.inductor = Inductor(self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('inductor')
        
        components['resistor'] = self.resistor
        components['inductor'] = self.inductor

        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r'resistance of resistor',
            self.inductance: r'inductance of inductor',
        }

        return self.sym_desc_dict
    
    
    
class CircutRC(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RC system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    capacity=Symbol('C', positive=True)
    q0=dynamicsymbols('q_c')
    qs=dynamicsymbols('qs')
    frame=Symbol('frame', positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 capacity=None,
                 ivar=None,
                 q0=None,
                 qs=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if capacity is not None: self.capacity = capacity
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if q0 is not None: self.q0 = q0
        if qs is not None: self.qs = qs
        if frame is not None: self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self.resistor = Resistor(self.resistance, self.q0,  qs=self.qs,ivar=self.ivar , frame=base_frame)('resistor')
        self.capacitor = Capacitor(self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('capacitor')
        
        components['resistor'] = self.resistor
        components['capacitor'] = self.capacitor

        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r'resistance of resistor',
            self.capacity: r'capacity of capacitor',
        }

        return self.sym_desc_dict
    
    
    
class CircuitRLC(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    inductance=Symbol('L', positive=True)
    capacity=Symbol('C', positive=True)
    q0=dynamicsymbols('q')
    frame=Symbol('frame', positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 inductance=None,
                 capacity=None,
                 ivar=None,
                 q0=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if inductance is not None: self.inductance = inductance
        if capacity is not None: self.capacity = capacity
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if q0 is not None: self.q0 = q0
        if frame is not None: self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        
        self.resistor = Resistor(self.resistance, self.q0,  qs=self.qs,ivar=self.ivar , frame=base_frame)('resistor')
        self.inductor = Inductor(self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('inductor')
        self.capacitor = Capacitor(self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('capacitor')
        
        components['resistor'] = self.resistor
        components['inductor'] = self.inductor
        components['capacitor'] = self.capacitor

        
        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r'resistance of resistor',
            self.inductance: r'inductance of inductor',
            self.capacity: r'capacity of capacitor',
        }

        return self.sym_desc_dict
    
class CircuitRLCWithPWM(CircuitRLC):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    inductance=Symbol('L', positive=True)
    capacity=Symbol('C', positive=True)
    U=Symbol('U',positive=True)
    eps=Symbol('varepsilon',postive=True)
    rho=Symbol('rho',positive=True)
    delta=Symbol('delta',positive=True)
    T=Symbol('T',positive=True)
    q0=dynamicsymbols('q')
    omega=Symbol('omega',positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 inductance=None,
                 capacity=None,
                 U=None,
                 eps=None,
                 rho=None,
                 delta=None,
                 T=None,
                 ivar=None,
                 q0=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if inductance is not None: self.inductance = inductance
        if capacity is not None: self.capacity = capacity
        if U is not None: self.U=U
        if eps is not None: self.eps=eps
        if rho is not None: self.rho=rho
        if delta is not None: self.delta=delta
        if T is not None: self.T=T
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if q0 is not None: self.q0 = q0
        if frame is not None: self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        self.r_var = self.resistance*(1+self.eps*self.rho  ) 
        
        
        self.resistor = Resistor(self.r_var, self.q0,  qs=self.qs,ivar=self.ivar , frame=base_frame)('resistor')
        self.inductor = Inductor(self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('inductor')
        self.capacitor = Capacitor(self.capacity/self.delta/self.eps, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('capacitor')
        self.voltage= VoltageSource(self.U,self.q0,ivar=self.ivar,qs=self.qs)
        
        components['resistor'] = self.resistor
        components['inductor'] = self.inductor
        components['capacitor'] = self.capacitor
        components['voltage'] = self.voltage

        
        return components
        
    def symbols_description(self):
        
        self.sym_desc_dict = {
            self.resistance: r'resistance',
            self.inductance: r'inductance',
            self.capacity: r'capacitance',
            self.U:r'source voltage',
            self.eps:r'small parameter',
            self.rho:r'formula representing resistance variation',
            self.q0:r'electric charge',
            self.q0.diff(self.ivar):r'current',
            self.q0.diff(self.ivar,self.ivar):r'current rate of change'

        }

        return self.sym_desc_dict
    
    def para(self, angle=2*pi):
        para_series_heave =  10*Heaviside(sin(2*pi/(self.T)*self.ivar),0.5)
        new_eq=self.subs(self.rho,para_series_heave)
        return new_eq
    
    def wave(self):
        
        wave1=(1*(sin(2*pi*self.ivar/self.T))+2.0)*(1/2 + 1/2*sign(sin(2*pi*self.ivar/self.T)))
        wave2=(1*(-sin(2*pi*self.ivar/self.T))-3.0)*(1/2+ 1/2*sign(-sin(2*pi*self.ivar/self.T)))
        waves=wave1+wave2
#         wave1=(100*(sin(self.ivar/self.T))+1000)*Heaviside(sin(self.ivar/self.T))
#         wave2=(100*(-sin(self.ivar/self.T)))*Heaviside(-sin(self.ivar/self.T))
        waves=wave1+wave2
        new_eq = self.subs(self.rho,waves)

        return new_eq
        
    def rect(self,no=6,numerical=False):
        t=self.ivar
        omg = self.omega
        rectangular=5*((1/2+1/2*sign(sin(2*pi*self.ivar/self.self.T)))-S.Half)
#         rectangular_func = lambdify((t, T), rectangular)
#         rectangular_values = rectangular_func(t_values, T_value)
        
#         trig=sum([Heaviside(omg*t-1) + 0* Heaviside(omg*t-2)  for ind in range(no)])
#         new_eq=self.subs(self.rho,trig)
        new_eq=self.subs(self.rho,rectangular)
        
        return new_eq
    
    def approx_rect(self,no=6,numerical=False):
        
        if numerical is True:
            amps_list = [2.27348466531425,0, 0.757408805199249,0, 0.453942816897038,0, 0.323708002807428,0, 0.25121830779797,0, 0.204977919963796,0, 0.172873394602606,0, 0.149252079729775,0, 0.131121653619234,0, 0.116749954968057,0]
        else:
            amps_list = symbols(f'a_0:{no}')
        rectangular_approx = sum([amp*2/sqrt(2)*sin((2*(no)+1)*2*pi*self.ivar/self.T) for no,amp in enumerate(amps_list[0:])])
        new_eq=self.subs({self.rho:rectangular_approx})
        
        return new_eq
    
#     def ode_with_delta(self):
#         delta=Symbol('delta', positive=True)
#         eps =Symbol('varepsilon', positive=True)
        
#         with_delta = self._eoms[0]+self.resistance*eps*delta*self.q0.diff(self.ivar)
#         delta_sys = type(self._ode_system)(with_delta, Matrix([self.q0]), ivar=self.ivar, ode_order=2)

#         return delta_sys

class CircuitRLCWithHeavisidePWM(CircuitRLCWithPWM):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    inductance=Symbol('L', positive=True)
    capacity=Symbol('C', positive=True)
    R_min=Symbol('R_min',positive=True)
    U=Symbol('U',positive=True)
    eps=Symbol('varepsilon',postive=True)
    rho=Symbol('rho',positive=True)
    delta=Symbol('delta',positive=True)
    T=Symbol('T',positive=True)
    q0=dynamicsymbols('q')
    omega=Symbol('omega',positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 R_min=None,
                 inductance=None,
                 capacity=None,
                 U=None,
                 eps=None,
                 rho=None,
                 delta=None,
                 T=None,
                 ivar=None,
                 q0=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if R_min is not None: self.R_min=R_min
        if inductance is not None: self.inductance = inductance
        if capacity is not None: self.capacity = capacity
        if U is not None: self.U=U
        if eps is not None: self.eps=eps
        if rho is not None: self.rho=rho
        if delta is not None: self.delta=delta
        if T is not None: self.T=T
        if ivar is not None: self.ivar = ivar
        if z is not None: self.z = z
        if q0 is not None: self.q0 = q0
        if frame is not None: self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        self.r_var = self.resistance*(1+self.eps*self.rho  ) 
        
        
        self.resistor = Resistor(self.R_min+self.resistance*Heaviside(sin(2*pi*self.ivar/self.T)), self.q0,  qs=self.qs,ivar=self.ivar , frame=base_frame)('resistor')
        self.inductor = Inductor(self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('inductor')
        self.capacitor = Capacitor(self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame)('capacitor')
        self.voltage= VoltageSource(self.U,self.q0,ivar=self.ivar,qs=self.qs)
        
        components['resistor'] = self.resistor
        components['inductor'] = self.inductor
        components['capacitor'] = self.capacitor
        components['voltage'] = self.voltage

        
        return components
    
    def symbols_description(self):
        
        self.sym_desc_dict = {
            self.resistance: r'amplitude of the variable part of the resistance',
            self.R_min:r'constant value of the resistance',
            self.inductance: r'inductance',
            self.capacity: r'capacitance',
            self.U:r'source voltage',
            self.q0:r'electric charge',
            self.q0.diff(self.ivar):r'current',
            self.q0.diff(self.ivar,self.ivar):r'current rate of change'

        }

        return self.sym_desc_dict
        

class BatteryCharging(ComposedSystem):
    
    V_OCV = Symbol('V_OCV',positive=True)
    current = Symbol('I_ch', positive=True)
    t0 = Symbol('t_0',positive=True)
    Tk = Symbol('T_k',positive=True)
    A = Symbol('A',positive=True)
    Phi = Symbol('\Phi', positive = True)
    resistance_0 = Symbol('R_0', positive=True)
    resistance_coefficient_alpha = Symbol('R_Alpha', positive=True) #symbol alfy
    resistance_coefficient_beta = Symbol('R_Beta', positive=True) #symbol bety
    resistance_coefficient_gamma = Symbol('R_Gamma', positive=True) #symbol gammy
    resistance_coefficient_blocking = Symbol('R_Omega', positive=True) #symbol omegi
    charger_blocking_resistance = Symbol('R_b', positive=True) 
    resistance  = Symbol('R' , positive= True)
    coil_coefficient = Symbol('L',positive=True)
    capacitor_coefficient= Symbol('C',positive=True)
    electric_charge = dynamicsymbols('q')
    k_p = Symbol ('k_p',positive=True)
    k_d = Symbol ('k_d',positive=True)
    PD = Symbol ('PD',positive=True)
    q_char = Symbol ('q_char',positive=True)

    def __init__(self,
                 V_OCV = V_OCV,
                 current = current,
                 t0 = t0,
                 Tk = Tk,
                 A = A,
                 Phi = Phi,
                 resistance_0 = resistance_0,
                 resistance_coefficient_alpha = resistance_coefficient_alpha,
                 resistance_coefficient_beta = resistance_coefficient_beta,
                 resistance_coefficient_gamma = resistance_coefficient_gamma,
                 resistance_coefficient_blocking = resistance_coefficient_blocking,
                 charger_blocking_resistance = charger_blocking_resistance,
                 resistance = resistance,
                 coil_coefficient = coil_coefficient,
                 capacitor_coefficient=capacitor_coefficient,
                 electric_charge = electric_charge,
                 k_p = k_p,
                 k_d = k_d,
                 PD = PD,
                 q_char = q_char,
                 qs=[electric_charge],
                 ivar=Symbol('t'),
                 **kwargs):
    
        t=ivar
        self.resistance_0 = resistance_0
        self.Tk= Tk
        self.A = A
        self.Phi = Phi
        self.resistance_coefficient_alpha =  resistance_coefficient_alpha
        self.resistance_coefficient_beta =  resistance_coefficient_beta
        self.resistance_coefficient_gamma =  resistance_coefficient_gamma
        self.resistance_coefficient_blocking =  resistance_coefficient_blocking
        self.coil_coefficient = coil_coefficient
        self.capacitor_coefficient = capacitor_coefficient
        self.electric_charge = electric_charge
        self.current = current
        self.qs = [self.electric_charge]
        
        ### krzywa U(SOC) w procesie ładowania
        x_char = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,     2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6650, 6800, 6850, 6900, 6950, 7000, 7050, 7100, 7150, 7200, 7250, 7300, 7350, 7400, 7450, 7500, 7550, 7600, 7650, 7700, 7750, 7800, 7850, 7900, 7950, 8000, 8050, 8100, 8150, 8200, 8250, 8300, 8350, 8400, 8450, 8500, 8550, 8600, 8650, 8800, 8900, 9400, 9800, 10200, 10600]

        y_char = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.55, 2.62, 2.74, 2.86, 2.98, 3.1, 3.22, 3.34, 3.46, 3.58, 3.7, 3.75, 3.801, 3.802, 3.803, 3.804, 3.805, 3.806, 3.807, 3.808, 3.809, 3.81, 3.811, 3.812, 3.813, 3.814, 3.815, 3.816, 3.817, 3.818, 3.819, 3.82, 3.821, 3.822, 3.823, 3.824, 3.825, 3.826, 3.827, 3.828, 3.829, 3.83, 3.831, 3.832, 3.833, 3.834, 3.835, 3.836, 3.837, 3.838, 3.839, 3.84, 3.841, 3.842, 3.843, 3.844, 3.845, 3.846, 3.847, 3.848, 3.849, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4, 4.01, 4.02, 4.05, 4.1, 4.15, 4.2, 4.22, 4.24, 4.26, 4.28, 4.3, 4.32, 4.34, 4.36, 4.38, 4.4, 4.42, 4.44, 4.46, 4.48, 4.5, 4.52, 4.6, 4.8, 5, 5.5, 6, 6.5, 7]
        char = np.polyfit(x_char,y_char,11)
        
        self.V_OCV  = char[0]*self.electric_charge**11 + char[1]*self.electric_charge**10 + char[2]*self.electric_charge**9 + char[3]*self.electric_charge**8 + char[4]*self.electric_charge**7 + char[5]*self.electric_charge**6 + char[6]*self.electric_charge**5 + char[7]*self.electric_charge**4 + char[8]*self.electric_charge**3 + char[9]*self.electric_charge**2 + char[10]*self.electric_charge + char[11]  
        
        
        self.charger_blocking_resistance = exp((self.electric_charge + (-1 * self.q_char - 250))/50) #self.resistance_coefficient_blocking
        
        self.resistance = self.resistance_0 + exp((self.resistance_coefficient_alpha - self.electric_charge)/self.resistance_coefficient_beta) +  exp((self.electric_charge + self.resistance_coefficient_gamma)/5550 ) 
        
        self.material_point_L = MaterialPoint(self.coil_coefficient,self.electric_charge, qs=[self.electric_charge])
        self.damper_R= Damper(self.resistance,self.electric_charge, qs=[self.electric_charge])
        self.damper_R_blocking= Damper(self.charger_blocking_resistance,self.electric_charge, qs=[self.electric_charge])
        self.force_OCV = Force(self.V_OCV,self.electric_charge,qs=[self.electric_charge])     
        self.spring_C = Spring(1/self.capacitor_coefficient, self.electric_charge, qs=[self.electric_charge])
        
        
        system = self.damper_R + self.damper_R_blocking + self.force_OCV + self.material_point_L + self.spring_C
        #super().__init__(system(qs),**kwargs)

        ###### regulator PD
        
        self.q_char = q_char
        self.t0 = ((self.q_char*0.0002777)/self.current)*3600 # ładunek żądany/stała(zamiana C na AH) dzielone przez prad ladowania  i zamiana na sekundy
        self.charger_slope = 0.001
        self.k_p = k_p 
        self.k_d = k_d   
        
        pd = Force((self.k_p* (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) *(self.current - self.electric_charge.diff(t)), self.electric_charge, qs = [self.electric_charge]  ) + Force((self.k_d * (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) * (self.current - self.electric_charge.diff(t)).diff(t),self.electric_charge ,qs = [self.electric_charge])
        
        system_pid = system+pd(qs)
        super().__init__(system_pid,**kwargs)
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance_0: r'Rezystancja podstawowa',
            self.resistance: r'Rezystor (SOC)',
            self.resistance_coefficient_alpha: r'Współczynnik rezystancji /alpha',
            self.resistance_coefficient_beta: r'Współczynnik rezystancji /beta',
            self.resistance_coefficient_gamma: r'Współczynnik rezystancji /gamma',
            self.charger_blocking_resistance: r'Rezystancja blokująca',
            self.resistance_coefficient_blocking: r'Współczynnik rezystancji blokującej', 
            self.capacitor_coefficient: r'Kondensator',
            self.coil_coefficient: r'Cewka indukcyjna',
            self.V_OCV: r'Napięcie',
            self.electric_charge: r'Ładunek',
            self.electric_charge.diff(t): r'Prąd',
            self.current: r'Prąd ładowania',
            self.t0: r'Czas ładowania',
            self.PD: r'Regulator PD',
            self.k_p: r'Człon proporcjonalny regulatora',
            self.k_d: r'Człon różniczkujący regulatora',
            self.ivar: r'Czas',
            
        }
        return self.sym_desc_dict


