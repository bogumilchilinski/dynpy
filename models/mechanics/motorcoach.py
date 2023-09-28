from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from dynpy import * # enables mechanical models for mathematical modelling


from sympy import * # provides mathematical interface for symbolic calculations


from sympy.physics.mechanics import *

import sympy as sym
import numpy as np
import numpy.fft as fft

import matplotlib.pyplot as plt
import math

from dynpy.utilities.report import SystemDynamicsAnalyzer

import pandas as pd

from dynpy.utilities.adaptable import *


from pint import UnitRegistry
ureg = UnitRegistry()

ix = pd.IndexSlice



mechanics_printing(pretty_print=True)

t,f= symbols('t, f')


from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST


#TODO 151
class FourDOFTrolleySuspension(ComposedSystem):

    scheme_name = 'car.PNG'
    real_name = 'car_real.jpg'
    
    z,phi,z_pw,z_lw=dynamicsymbols('z, \\varphi z_pw z_lw')
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_k=Symbol('m_k', positive=True)
    m_p=Symbol('m_p', positive=True)
    m_l=Symbol('m_l', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_k=Symbol('l_k', positive=True)
    l_p=Symbol('l_p', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_l=Symbol('k_l1', positive=True)
    k_lw=Symbol('k_lw', positive=True)
    k_p=Symbol('k_p1', positive=True)
    k_pw=Symbol('k_pw', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_p=Symbol('c_p', positive=True)
    c_lw=Symbol('c_lw', positive=True)
    c_pw=Symbol('c_pw', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    ivar=Symbol('t')
    f=Symbol('f', positive=True)
    
    
    def __init__(self,
                 m=None,
                 m_k=None,
                 m_p=None,
                 m_l=None,
                 m_b=None,
                 I=None,
                 l_rod=None,
                 l_l=None,
                 l_p=None,
                 l_b=None,
                 l_k=None,
                 k_l=None,
                 k_lw=None,
                 k_p=None,
                 k_pw=None,
                 c_l=None,
                 c_p=None,
                 c_lw=None,
                 c_pw=None,
                 f=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 s=None,
                 phi=None,
                 z_pw=None,
                 z_lw=None,
                 A=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
            
                
            
                 **kwargs):
        
        if z is not None: self.z=z
        if z_lw is not None: self.z_lw=z_lw
        if z_pw is not None: self.z_pw=z_pw
        if z_l is not None: self.z_lw=z_lw
        if z_p is not None: self.z_pw=z_pw
        if z_b is not None: self.z_b=z_b
        if phi is not None: self.phi=phi

        if m is not None: self.m = m  # mass of a rod
        if m_p is not None: self.m_p = m_p
        if m_l is not None: self.m_l = m_l
        if m_b is not None: self.m_b = m_b
        if m_k is not None: self.m_k = m_k
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_p is not None: self.l_p = l_p #offset of right spring
        if l_b is not None: self.l_b = l_b  
        if l_k is not None: self.l_k = l_k  
        if l_rod is not None: self.l_rod = l_rod
        if k_l is not None: self.k_l1 = k_l1
        if k_lw is not None: self.k_lw = k_lw
        if k_p is not None: self.k_p1 = k_p1
        if k_pw is not None: self.k_pw = k_pw
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_lw is not None: self.c_lw=c_lw
        if c_pw is not None: self.c_pw=c_pw
        if A is not None: self.A=A
        if f is not None: self.f=f
        if omega is not None: self.omega=omega
        if I is not None: self.I = I  # moment of inertia of a rod
        if F_engine is not None: self.F_engine = F_engine

        
        self.s=self.A*sin(ivar*self.omega)
        
       
        
        self.right_wheel = MaterialPoint(self.m_p, pos1=self.z_pw, qs=[self.z_pw])
        self.left_wheel = MaterialPoint(self.m_l, pos1=self.z_lw, qs=[self.z_lw])
        self.spring_pw = Spring(self.k_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z_pw])
        self.spring_lw = Spring(self.k_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z_lw])
        self.damper_pw = Damper(self.c_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z_pw])
        self.damper_lw = Damper(self.c_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z_lw])
        
   
        
        self.body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=[self.z, self.phi])  # rod
        self.battery= MaterialPoint(self.m_b, pos1=self.z-self.phi*self.l_b, qs=[self.z, self.phi])
#         self.pasazer= MaterialPoint(self.m_k, pos1=self.z+self.phi*self.l_k, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_l, pos1=self.z+self.phi*self.l_l , pos2 = self.z_lw , qs=[self.z, self.phi, self.z_lw])  # left spring
        self.spring_2 = Spring(self.k_p, pos1=self.z-self.phi*self.l_p , pos2 = self.z_pw , qs=[self.z, self.phi, self.z_pw])
        self.damper_1= Damper(self.c_l, pos1=self.z+self.phi*self.l_l , pos2 = self.z_lw , qs=[self.z, self.phi, self.z_lw])
        self.damper_2 = Damper(self.c_p, pos1=self.z-self.phi*self.l_p , pos2 = self.z_pw , qs=[self.z, self.phi, self.z_pw])
        #self.force_1 = Force(self.A*sin(ivar*self.omega), pos1=self.z_lw , qs=[self.z_lw])
        #self.force_2 = Force(self.A*sin(ivar*self.omega), pos1=self.z_pw, qs=[self.z_pw])
        
        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body +self.battery+self.spring_1 + self.spring_2 +self.damper_1+self.damper_2 + self.right_wheel + self.left_wheel + self.spring_pw+ self.spring_lw+self.damper_lw+self.damper_pw

        super().__init__(**{'system':system,**kwargs})
        


    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Masa resorowana',
            self.m_p: r'Masa przedniej osi',
            self.m_l: r'Masa tylnej osi',
            self.m_b: r'Masa baterii trakcyjnej',
            self.m_b: r'Masa pasazera',
            self.I: r'Moment bezwładności',
            self.l_rod: r'Długość osi',
            self.l_p: r'Odległość środka ciężkości masy nadwozia autobusu do punktu styku przednich kół w kierunku poziomym',
            self.l_l: r'Odległość środka ciężkości masy nadwozia autobusu do punktu styku tylnych kół w kierunku poziomym',
            self.l_b: r'Odległość środka ciężkości masy nadwozia autobusu do środka ciężkości baterii',
            self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
            self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
            self.k_pw: r'Współczynnik sztywności dla przedniej opony',
            self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
            self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
            self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
            self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
            self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
            self.phi:  r'Kąt obrotu masy resorowanej',
            self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
            self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
            self.z_lw: r'Przemieszczenie pionowe tylnego koła',
            self.z_b: r'Przemieszczenie pionowe baterii',
            self.A: r'Amplituda siły wymuszającej',
            self.omega: r'Częstość siły wymuszającej',
            self.F_engine: r'Siła',
            self.ivar: r'Czas',
            self.f: r'Częstotliwość wymuszenia',
            
            
        }
        return self.sym_desc_dict
    
sys4=FourDOFTrolleySuspension()
y=[sys4.q,sys4.z, sys4.z.diff(t,t)]
ic_list = [0.1,0.1,0,0,0.1,0.1,0,0]


units_dict = {sys4.m:ureg.kilogram,
              sys4.m_p:ureg.kilogram,
              sys4.m_l:ureg.kilogram,
              sys4.m_b:ureg.kilogram,
              sys4.m_k:ureg.kilogram,
              sys4.A:ureg.meter,
              sys4.k_l:(ureg.newton / ureg.meter),
              sys4.k_lw:(ureg.newton / ureg.meter),
              sys4.k_p:(ureg.newton / ureg.meter),
              sys4.k_pw:(ureg.newton / ureg.meter),
              sys4.c_l:(ureg.newton*ureg.second / ureg.meter),
              sys4.c_lw:(ureg.newton*ureg.second / ureg.meter),
              sys4.c_p:(ureg.newton*ureg.second / ureg.meter),
              sys4.c_pw:(ureg.newton*ureg.second / ureg.meter),
              sys4.omega:ureg.radian,
              sys4.l_l:ureg.meter,
              sys4.l_p:ureg.meter,
              sys4.l_b:ureg.meter,
              sys4.l_rod:ureg.meter,
              sys4.I:(ureg.kilogram*ureg.meter*ureg.meter),
              sys4.phi:ureg.radian,
              sys4.z:ureg.meter,
              sys4.z.diff(t,t):ureg.meter/ureg.second**2,
              sys4.z_lw:ureg.meter,
              sys4.z_pw:ureg.meter,
              sys4.z_b:ureg.meter,
              sys4.F_engine:ureg.newton,
              t:ureg.second,
              f:ureg.hertz,
              
             }

unit=units_dict
#TODO        
class DDOFTrolleySuspension(ComposedSystem):

   
    
    z,phi,z_pw,z_lw=dynamicsymbols('z, \\varphi z_pw z_lw')
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_p=Symbol('m_p', positive=True)
    m_l=Symbol('m_l', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_p=Symbol('l_p', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_l=Symbol('k_l1', positive=True)
    k_lw=Symbol('k_lw', positive=True)
    k_p=Symbol('k_p1', positive=True)
    k_pw=Symbol('k_pw', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_p=Symbol('c_p', positive=True)
    c_lw=Symbol('c_lw', positive=True)
    c_pw=Symbol('c_pw', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    ivar=Symbol('t')
    k_mr=Symbol('k_mr',positive=True)
    m_nr=Symbol('m_nr')
    
    def __init__(self,
                  m=None,
                 m_p=None,
                 m_l=None,
                 m_b=None,
                 I=None,
                 l_rod=None,
                 l_l=None,
                 l_p=None,
                 l_b=None,
                 k_l=None,
                 k_lw=None,
                 k_p=None,
                 k_pw=None,
                 c_l=None,
                 c_p=None,
                 c_lw=None,
                 c_pw=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 s=None,
                 phi=None,
                 z_pw=None,
                 z_lw=None,
                 A=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
                 k_mr=None,
                 m_nr=None,
                 **kwargs):
        
        if z is not None: self.z=z
        if phi is not None: self.phi=phi
        

        if m is not None: self.m = m  # mass of a rod
        if m_p is not None: self.m_w = m_w
        if m_l is not None: self.m_w = m_w
        if m_b is not None: self.m_b = m_b
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_p is not None: self.l_p = l_p #offset of right spring
        if l_b is not None: self.l_p = l_b  
        if l_rod is not None: self.l_rod = l_rod
        if k_l is not None: self.k_l1 = k_l1
        if k_lw is not None: self.k_lw = k_lw
        if k_p is not None: self.k_p1 = k_p1
        if k_pw is not None: self.k_pw = k_pw
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_lw is not None: self.c_lw=c_lw
        if c_pw is not None: self.c_pw=c_pw
        if A is not None: self.A=A
        if omega is not None: self.omega=omega
        if I is not None: self.I = I  # moment of inertia of a rod
        if F_engine is not None: self.F_engine = F_engine     
        if k_mr is not None: self.k_mr=k_mr
        if m_nr is not None: self.m_nr=m_nr
            
            
        self.s=self.A*sin(ivar*self.omega)
        
        self.z_l=self.z+self.phi*self.l_l
        self.z_p=self.z-self.phi*self.l_p
        self.z_lw=(self.k_l*self.z_l)/(self.k_lw+self.k_l)
        self.z_pw=(self.k_p*self.z_p)/(self.k_pw+self.k_p)
        
        
        self.right_wheel = MaterialPoint(self.m_p, pos1=self.z_pw, qs=[self.z, self.phi])
        self.left_wheel = MaterialPoint(self.m_l, pos1=self.z_lw, qs=[self.z, self.phi])
        self.spring_pw = Spring(self.k_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        self.spring_lw = Spring(self.k_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        self.damper_pw = Damper(self.c_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        self.damper_lw = Damper(self.c_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring=Spring(k,pos1,po2,qs)
        
        self.body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=[self.z, self.phi])  # rod
        self.battery= MaterialPoint(self.m_b, pos1=self.z+self.phi*self.l_b, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_l, pos1=self.z_l , pos2 = self.z_lw , qs=[self.z, self.phi])  # left spring
        self.spring_2 = Spring(self.k_p, pos1=self.z_p, pos2=self.z_pw, qs=[self.z, self.phi])
        self.damper_1=Damper(self.c_l, pos1=self.z_p , pos2 = self.z_lw , qs=[self.z, self.phi])
        self.damper_2 = Damper(self.c_p, pos1=self.z_l , pos2 = self.z_pw , qs=[self.z, self.phi])

        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body +self.battery+ self.spring_1 + self.spring_2 + self.damper_1+ self.damper_2 + self.right_wheel + self.left_wheel + self.spring_pw + self.spring_lw+self.damper_lw+self.damper_pw
        
#         display(type(system))
        
#         system_new = system.subs({self.z_lw:self.z_l2,self.z_pw:self.z_p2}) 
        
#         display(type(system_new))
#         display(type(self.body))
        
        #super().__init__(system_new,**kwargs)

        super().__init__(**{'system':system,**kwargs})
        


        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Masa resorowana',
            self.m_p: r'Masa przedniej osi',
            self.m_l: r'Masa tylnej osi',
            self.m_b: r'Masa baterii trakcyjnej',
            self.I: r'Moment bezwładności',
            self.l_rod: r'Długość osi',
            self.l_p: r'Odległość od przedniej osi do środka ciężkości autobusu',
            self.l_l: r'Odległość od tylnej osi do środka ciężkości autobusu',
            self.l_b: r'Odległość od środka cieżkości baterii do środka ciężkości autobusu',
            self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
            self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
            self.k_pw: r'Współczynnik sztywności dla przedniej opony',
            self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
            self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
            self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
            self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
            self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
            self.phi:  r'Kąt obrotu masy resorowanej',
            self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
            self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
            self.z_lw: r'Przemieszczenie pionowe tylnego koła',
            self.z_b: r'Przemieszczenie pionowe baterii',
            self.A: r'Amplituda siły wymuszającej',
            self.omega: r'Częstość siły wymuszającej',
            self.ivar: r'Czas',
            self.k_mr: r'Zastępcza sztywność mas resorowanych',
            self.m_nr: r'Zastępcza masa nierosowana',
            
        }
        return self.sym_desc_dict
    
sys2=DDOFTrolleySuspension()



units_dict1 = {sys2.m:ureg.kilogram,
              sys2.m_p:ureg.kilogram,
              sys2.m_l:ureg.kilogram,
              sys2.m_b:ureg.kilogram,
              sys2.A:ureg.meter,
              sys2.k_l:(ureg.newton / ureg.meter),
              sys2.k_lw:(ureg.newton / ureg.meter),
              sys2.k_p:(ureg.newton / ureg.meter),
              sys2.k_pw:(ureg.newton / ureg.meter),
              sys2.c_l:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_lw:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_p:(ureg.newton*ureg.second / ureg.meter),
              sys2.c_pw:(ureg.newton*ureg.second / ureg.meter),
              sys2.omega:ureg.radian,
              sys2.l_l:ureg.meter,
              sys2.l_p:ureg.meter,
              sys2.l_b:ureg.meter,
              sys2.l_rod:ureg.meter,
              sys2.I:(ureg.kilogram*ureg.meter*ureg.meter),
              sys2.phi:ureg.radian,
              sys2.z:ureg.meter,
              sys2.z.diff(t,t):ureg.meter/ureg.second**2,
              sys2.z_lw:ureg.meter,
              sys2.z_pw:ureg.meter,
              sys2.z_b:ureg.meter,
              sys2.F_engine:ureg.newton,
              t:ureg.second,
              f:ureg.hertz,
              sys2.k_mr:(ureg.newton / ureg.meter),
              sys2.m_nr:ureg.kilogram,
             }

unit1=units_dict1
        #TODO 153
class DDOFTrolleySuspension2(ComposedSystem):

   
    
    z,phi,z_pw,z_lw=dynamicsymbols('z, \\varphi z_pw z_lw')
    m=Symbol('m', positive=True)
    m_b=Symbol('m_b', positive=True)
    m_p=Symbol('m_p', positive=True)
    m_l=Symbol('m_l', positive=True)
    I=Symbol('I', positive=True)
    l_rod=Symbol('l_{rod}', positive=True)
    l_l=Symbol('l_l', positive=True)
    l_p=Symbol('l_p', positive=True)
    l_b=Symbol('l_b', positive=True)
    k_l=Symbol('k_l1', positive=True)
    k_lw=Symbol('k_lw', positive=True)
    k_p=Symbol('k_p1', positive=True)
    k_pw=Symbol('k_pw', positive=True)
    c_l=Symbol('c_l', positive=True)
    c_p=Symbol('c_p', positive=True)
    c_lw=Symbol('c_lw', positive=True)
    c_pw=Symbol('c_pw', positive=True)
    F_engine=Symbol('F_{engine}', positive=True)
    k_p=Symbol('k_p', positive=True)
    k_l=Symbol('k_l', positive=True)
    A=Symbol('A', positive=True)
    omega=Symbol('omega', positive=True)
    z_l=Symbol('z_l', positive=True)
    z_p=Symbol('z_p',positive=True)
    z_b=Symbol('z_b',positive=True)
    s=Symbol('s', positive=True)
    k_l_zas=Symbol('k_l_zas', positive=True)
    k_p_zas=Symbol('k_p_zas', positive=True)
    c_l_zas=Symbol('c_l_zas', positive=True)
    c_p_zas=Symbol('c_lp_zas', positive=True)
    ivar=Symbol('t')
    
    def __init__(self,
                 m=None,
                 m_p=None,
                 m_l=None,
                 m_b=None,
                 I=None,
                 l_rod=None,
                 l_l=None,
                 l_p=None,
                 l_b=None,
                 k_l=None,
                 k_lw=None,
                 k_p=None,
                 k_pw=None,
                 c_l=None,
                 c_p=None,
                 c_lw=None,
                 c_pw=None,
                 F_engine=None,
                 ivar=Symbol('t'),
                 z=None,
                 phi=None,
                 z_pw=None,
                 z_lw=None,
                 A=None,
                 omega=None,
                 z_l=None,
                 z_p=None,
                 z_b=None,
                 s=None,
                 k_l_zas=None,
                 k_p_zas=None,
                 c_l_zas=None,
                 c_p_zas=None,
                 
                 **kwargs):
        if z is not None: self.z=z
        if z_lw is not None: self.z_lw=z_lw
        if z_pw is not None: self.z_pw=z_pw
        if z_l is not None: self.z_lw=z_lw
        if z_p is not None: self.z_pw=z_pw
        if z_b is not None: self.z_b=z_b
        if phi is not None: self.phi=phi

        if m is not None: self.m = m  # mass of a rod
        if m_p is not None: self.m_w = m_w
        if m_l is not None: self.m_w = m_w
        if m_b is not None: self.m_b = m_b
        if l_l is not None: self.l_l = l_l  # offset of left spring
        if l_p is not None: self.l_p = l_p #offset of right spring
        if l_b is not None: self.l_p = l_b  
        if l_rod is not None: self.l_rod = l_rod
        if k_l is not None: self.k_l1 = k_l1
        if k_lw is not None: self.k_lw = k_lw
        if k_p is not None: self.k_p1 = k_p1
        if k_pw is not None: self.k_pw = k_pw
        if c_l is not None: self.c_l=c_l
        if c_p is not None: self.c_p=c_p
        if c_lw is not None: self.c_lw=c_lw
        if c_pw is not None: self.c_pw=c_pw
        if A is not None: self.A=A
        if omega is not None: self.omega=omega
        if I is not None: self.I = I  # moment of inertia of a rod
        if F_engine is not None: self.F_engine = F_engine
        
                 
        self.s=self.A*sin(ivar*self.omega)
        
        self.k_l_zas=(self.k_l*self.k_lw)/(self.k_l+self.k_lw)
        self.k_p_zas=(self.k_p*self.k_pw)/(self.k_p+self.k_pw)
        self.c_l_zas=(self.c_l*self.c_lw)/(self.c_l+self.c_lw)
        self.c_p_zas=(self.c_p*self.c_pw)/(self.c_p+self.c_pw)
        
        self.z_l=self.z+self.phi*self.l_l
        self.z_p=self.z-self.phi*self.l_p
        self.z_lw=(self.k_l*self.z_l)/(self.k_lw+self.k_l)
        self.z_pw=(self.k_p*self.z_p)/(self.k_pw+self.k_p)
        
        #self.right_wheel = MaterialPoint(self.m_w, pos1=self.z_pw, qs=[self.z, self.phi])
        #self.left_wheel = MaterialPoint(self.m_w, pos1=self.z_lw, qs=[self.z, self.phi])
        #self.spring_pw = Spring(self.k_pw, pos1=self.z_pw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring_lw = Spring(self.k_lw, pos1=self.z_lw, pos2=self.s, qs=[self.z, self.phi])
        #self.spring=Spring(k,pos1,po2,qs)
        
        self.body = RigidBody2D(self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=[self.z, self.phi])  # rod
        self.battery= MaterialPoint(self.m_b, pos1=self.z+self.phi*self.l_b, qs=[self.z, self.phi])
        self.spring_1 = Spring(self.k_l_zas, pos1=self.z_l , pos2 = self.s , qs=[self.z, self.phi])  # left spring
        self.spring_2 = Spring(self.k_p_zas, pos1=self.z_p, pos2=self.s, qs=[self.z, self.phi])
        self.damper_1=Damper(self.c_l_zas, pos1=self.z_p , pos2 = self.s, qs=[self.z, self.phi])
        self.damper_2 = Damper(self.c_p_zas, pos1=self.z_l , pos2 = self.s , qs=[self.z, self.phi])

        #self.force = Force(-self.F_engine, pos1=self.z - self.l_p * self.phi, qs=[self.z, self.phi])
        system = self.body +self.battery+ self.spring_1 + self.spring_2 + self.damper_1+ self.damper_2 
        
#         display(type(system))
        
#         system_new = system.subs({self.z_lw:self.z_l2,self.z_pw:self.z_p2}) 
        
#         display(type(system_new))
#         display(type(self.body))
        
        #super().__init__(system_new,**kwargs)

        super().__init__(**{'system':system,**kwargs})        

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Masa resorowana',
            self.m_p: r'Masa przedniej osi',
            self.m_l: r'Masa tylnej osi',
            self.m_b: r'Masa baterii trakcyjnej',
            self.I: r'Moment bezwładności',
            self.l_rod: r'Długość osi',
            self.l_p: r'Odległość od przedniej osi do środka ciężkości autobusu',
            self.l_l: r'Odległość od tylnej osi do środka ciężkości autobusu',
            self.l_b: r'Odległość od środka cieżkości baterii do środka ciężkości autobusu',
            self.k_p: r'Współczynnik sztywności dla sprężyny przedniej nadwozia',
            self.k_l: r'Współczynnik sztywności dla sprężyny tylnej nadwozia',
            self.k_pw: r'Współczynnik sztywności dla przedniej opony',
            self.k_lw:  r'Współczynnik sztywności dla tylnej opony',
            self.c_l:  r'Współczynnik tłumienia dla amortzatora przedniego n',
            self.c_l: r'Współczynnik tłumienia dla amortzatora tylnego',
            self.c_lw:  r'Współczynnik tłumienia dla tylnej opony',
            self.c_pw:  r'Współczynnik tłumienia dla przedniej opony',
            self.phi:  r'Kąt obrotu masy resorowanej',
            self.z: r'Przemieszczenie pionowe środka cieżkości masy resorowanej',
            self.z_pw:  r'Przemieszczenie pionowe przedniego koła',
            self.z_lw: r'Przemieszczenie pionowe tylnego koła',
            self.z_b: r'Przemieszczenie pionowe baterii',
            self.A: r'Amplituda siły wymuszającej',
            self.omega: r'Częstość siły wymuszającej',
            self.ivar: r'Czas',
            self.k_l_zas: r'Zastępcza wartość współczynnika sztywności zawieszenia tylnej osi',
            self.k_p_zas: r'Zastępcza wartość współczynnika sztywności zawieszenia przedniej osi',
            self.c_l_zas: r'Zastępcza wartość współczynnika tłumienia zawieszenia tylnej osi',
            self.c_p_zas: r'Zastępcza wartość współczynnika tłumienia zawieszenia przedniej osi',
        }
        return self.sym_desc_dict
                 
sys22=DDOFTrolleySuspension2()      
                 
units_dict2 = {sys22.m:ureg.kilogram,
              sys22.m_p:ureg.kilogram,
              sys22.m_l:ureg.kilogram,
              sys22.m_b:ureg.kilogram,
              sys22.A:ureg.meter,
              sys22.k_l:(ureg.newton / ureg.meter),
              sys22.k_lw:(ureg.newton / ureg.meter),
              sys22.k_p:(ureg.newton / ureg.meter),
              sys22.k_pw:(ureg.newton / ureg.meter),
              sys22.c_l:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_lw:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_p:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_pw:(ureg.newton*ureg.second / ureg.meter),
              sys22.omega:ureg.radian,
              sys22.l_l:ureg.meter,
              sys22.l_p:ureg.meter,
              sys22.l_b:ureg.meter,
              sys22.l_rod:ureg.meter,
              sys22.I:(ureg.kilogram*ureg.meter*ureg.meter),
              sys22.phi:ureg.radian,
              sys22.z:ureg.meter,
              sys22.z.diff(t,t):ureg.meter/ureg.second**2,
              sys22.z_lw:ureg.meter,
              sys22.z_pw:ureg.meter,
              sys22.z_b:ureg.meter,
              t:ureg.second,
              f:ureg.hertz,
              sys22.F_engine:ureg.newton,
              sys22.k_p_zas:(ureg.newton / ureg.meter),
              sys22.k_l_zas:(ureg.newton / ureg.meter),
              sys22.c_l_zas:(ureg.newton*ureg.second / ureg.meter),
              sys22.c_p_zas:(ureg.newton*ureg.second / ureg.meter),
             }

unit2=units_dict2
   

