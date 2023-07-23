from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, CombustionEngine

from .elements import Resistor, Inductor, Capacitor

import base64
import random
import IPython as IP
import numpy as np
import inspect


from dynpy.models.mechanics.trolley import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin

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
    
    
    
    class CircutRLC(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """
    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'
    

    resistance=Symbol('R', positive=True)
    inductance=Symbol('L', positive=True)
    capacity=Symbol('C', positive=True)
    q0=dynamicsymbols('q_c')
    qs=dynamicsymbols('qs')
    frame=Symbol('frame', positive=True)
    ivar=Symbol('t')
    
    
    
    def __init__(self,
                 resistance=None,
                 inductance=None,
                 capacity=None,
                 ivar=None,
                 q0=None,
                 qs=None,
                 frame=None,
                 z=None,
                 **kwargs):

        
        
        if resistance is not None: self.resistance = resistance
        if inductance is not None: self.inductance = inductance
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
# class BatteryCharging(ComposedSystem):
    
#     V_OCV = Symbol('V_OCV',positive=True)
#     current = Symbol('I_ch', positive=True)
#     t0 = Symbol('t_0',positive=True)
#     Tk = Symbol('T_k',positive=True)
#     A = Symbol('A',positive=True)
#     Phi = Symbol('\Phi', positive = True)
#     resistance_0 = Symbol('R_0', positive=True)
#     resistance_coefficient_alpha = Symbol('R_Alpha', positive=True) #symbol alfy
#     resistance_coefficient_beta = Symbol('R_Beta', positive=True) #symbol bety
#     resistance_coefficient_gamma = Symbol('R_Gamma', positive=True) #symbol gammy
#     resistance_coefficient_blocking = Symbol('R_Omega', positive=True) #symbol omegi
#     charger_blocking_resistance = Symbol('R_b', positive=True) 
#     resistance  = Symbol('R' , positive= True)
#     coil_coefficient = Symbol('L',positive=True)
#     capacitor_coefficient= Symbol('C',positive=True)
#     electric_charge = dynamicsymbols('q')
#     k_p = Symbol ('k_p',positive=True)
#     k_d = Symbol ('k_d',positive=True)
#     PD = Symbol ('PD',positive=True)
#     q_char = Symbol ('q_char',positive=True)

#     def __init__(self,
#                  V_OCV = V_OCV,
#                  current = current,
#                  t0 = t0,
#                  Tk = Tk,
#                  A = A,
#                  Phi = Phi,
#                  resistance_0 = resistance_0,
#                  resistance_coefficient_alpha = resistance_coefficient_alpha,
#                  resistance_coefficient_beta = resistance_coefficient_beta,
#                  resistance_coefficient_gamma = resistance_coefficient_gamma,
#                  resistance_coefficient_blocking = resistance_coefficient_blocking,
#                  charger_blocking_resistance = charger_blocking_resistance,
#                  resistance = resistance,
#                  coil_coefficient = coil_coefficient,
#                  capacitor_coefficient=capacitor_coefficient,
#                  electric_charge = electric_charge,
#                  k_p = k_p,
#                  k_d = k_d,
#                  PD = PD,
#                  q_char = q_char,
#                  qs=[electric_charge],
#                  ivar=Symbol('t'),
#                  **kwargs):
    
#         t=ivar
#         self.resistance_0 = resistance_0
#         self.Tk= Tk
#         self.A = A
#         self.Phi = Phi
#         self.resistance_coefficient_alpha =  resistance_coefficient_alpha
#         self.resistance_coefficient_beta =  resistance_coefficient_beta
#         self.resistance_coefficient_gamma =  resistance_coefficient_gamma
#         self.resistance_coefficient_blocking =  resistance_coefficient_blocking
#         self.coil_coefficient = coil_coefficient
#         self.capacitor_coefficient = capacitor_coefficient
#         self.electric_charge = electric_charge
#         self.current = current
#         self.qs = [self.electric_charge]
        
#         ### krzywa U(SOC) w procesie ładowania
#         x_char = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,     2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6650, 6800, 6850, 6900, 6950, 7000, 7050, 7100, 7150, 7200, 7250, 7300, 7350, 7400, 7450, 7500, 7550, 7600, 7650, 7700, 7750, 7800, 7850, 7900, 7950, 8000, 8050, 8100, 8150, 8200, 8250, 8300, 8350, 8400, 8450, 8500, 8550, 8600, 8650, 8800, 8900, 9400, 9800, 10200, 10600]

#         y_char = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.55, 2.62, 2.74, 2.86, 2.98, 3.1, 3.22, 3.34, 3.46, 3.58, 3.7, 3.75, 3.801, 3.802, 3.803, 3.804, 3.805, 3.806, 3.807, 3.808, 3.809, 3.81, 3.811, 3.812, 3.813, 3.814, 3.815, 3.816, 3.817, 3.818, 3.819, 3.82, 3.821, 3.822, 3.823, 3.824, 3.825, 3.826, 3.827, 3.828, 3.829, 3.83, 3.831, 3.832, 3.833, 3.834, 3.835, 3.836, 3.837, 3.838, 3.839, 3.84, 3.841, 3.842, 3.843, 3.844, 3.845, 3.846, 3.847, 3.848, 3.849, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4, 4.01, 4.02, 4.05, 4.1, 4.15, 4.2, 4.22, 4.24, 4.26, 4.28, 4.3, 4.32, 4.34, 4.36, 4.38, 4.4, 4.42, 4.44, 4.46, 4.48, 4.5, 4.52, 4.6, 4.8, 5, 5.5, 6, 6.5, 7]
#         char = np.polyfit(x_char,y_char,11)
        
#         self.V_OCV  = char[0]*self.electric_charge**11 + char[1]*self.electric_charge**10 + char[2]*self.electric_charge**9 + char[3]*self.electric_charge**8 + char[4]*self.electric_charge**7 + char[5]*self.electric_charge**6 + char[6]*self.electric_charge**5 + char[7]*self.electric_charge**4 + char[8]*self.electric_charge**3 + char[9]*self.electric_charge**2 + char[10]*self.electric_charge + char[11]  
        
        
#         self.charger_blocking_resistance = exp((self.electric_charge + (-1 * self.q_char - 250))/50) #self.resistance_coefficient_blocking
        
#         self.resistance = self.resistance_0 + exp((self.resistance_coefficient_alpha - self.electric_charge)/self.resistance_coefficient_beta) +  exp((self.electric_charge + self.resistance_coefficient_gamma)/5550 ) 
        
#         self.material_point_L = MaterialPoint(self.coil_coefficient,self.electric_charge, qs=[self.electric_charge])
#         self.damper_R= Damper(self.resistance,self.electric_charge, qs=[self.electric_charge])
#         self.damper_R_blocking= Damper(self.charger_blocking_resistance,self.electric_charge, qs=[self.electric_charge])
#         self.force_OCV = Force(self.V_OCV,self.electric_charge,qs=[self.electric_charge])     
#         self.spring_C = Spring(1/self.capacitor_coefficient, self.electric_charge, qs=[self.electric_charge])
        
        
#         system = self.damper_R + self.damper_R_blocking + self.force_OCV + self.material_point_L + self.spring_C
#         #super().__init__(system(qs),**kwargs)

#         ###### regulator PD
        
#         self.q_char = q_char
#         self.t0 = ((self.q_char*0.0002777)/self.current)*3600 # ładunek żądany/stała(zamiana C na AH) dzielone przez prad ladowania  i zamiana na sekundy
#         self.charger_slope = 0.001
#         self.k_p = k_p 
#         self.k_d = k_d   
        
#         pd = Force((self.k_p* (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) *(self.current - self.electric_charge.diff(t)), self.electric_charge, qs = [self.electric_charge]  ) + Force((self.k_d * (S.One/2-atan(self.charger_slope*(self.ivar-self.t0))/pi)) * (self.current - self.electric_charge.diff(t)).diff(t),self.electric_charge ,qs = [self.electric_charge])
        
#         system_pid = system+pd(qs)
#         super().__init__(system_pid,**kwargs)
        
#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.resistance_0: r'Rezystancja podstawowa',
#             self.resistance: r'Rezystor (SOC)',
#             self.resistance_coefficient_alpha: r'Współczynnik rezystancji /alpha',
#             self.resistance_coefficient_beta: r'Współczynnik rezystancji /beta',
#             self.resistance_coefficient_gamma: r'Współczynnik rezystancji /gamma',
#             self.charger_blocking_resistance: r'Rezystancja blokująca',
#             self.resistance_coefficient_blocking: r'Współczynnik rezystancji blokującej', 
#             self.capacitor_coefficient: r'Kondensator',
#             self.coil_coefficient: r'Cewka indukcyjna',
#             self.V_OCV: r'Napięcie',
#             self.electric_charge: r'Ładunek',
#             self.electric_charge.diff(t): r'Prąd',
#             self.current: r'Prąd ładowania',
#             self.t0: r'Czas ładowania',
#             self.PD: r'Regulator PD',
#             self.k_p: r'Człon proporcjonalny regulatora',
#             self.k_d: r'Człon różniczkujący regulatora',
#             self.ivar: r'Czas',
            
#         }
#         return self.sym_desc_dict

