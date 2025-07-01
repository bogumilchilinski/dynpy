import importlib

from . import dynamics as dyn

importlib.reload(dyn)
import numpy as np
import pint
import pylatex
import sympy as sym
from pylatex import Command, NoEscape
from sympy import *
from sympy.physics.mechanics import *
from sympy.physics.vector import dynamicsymbols

from dynpy.dynamics import HarmonicOscillator, LagrangesDynamicSystem, mech_comp
from dynpy.models.elements import (
    PID,
    Damper,
    Disk,
    Excitation,
    Force,
    GravitationalForce,
    MaterialPoint,
    RigidBody2D,
    Spring,
    base_frame,
    base_origin,
)
from dynpy.models.mechanics.trolley import SpringDamperMassSystem, SpringMassSystem

from . import dynamics as dyn

# from dynpy.models.sdof import ComposedSystem, SpringMassSystem, DampedSpringMassSystem, BeamBridgeDamped
ureg = pint.UnitRegistry()
from dynpy.models.mechanics.principles import (
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)

mechanics_printing()


t = Symbol("t")  # independent variable - time


class LEWEL16_simple(SpringMassSystem):
    scheme_name = "train_sdof.png"
    real_name = "lewel16.jpg"


class LEWEL16_simple_damped(SpringDamperMassSystem):
    scheme_name = "train_sdof_damped.png"
    real_name = "lewel16.jpg"


class LeafSpring(ComposedSystem):

    m = Symbol("m", positive=True)
    k_beam = Symbol("k_beam", positive=True)
    ivar = Symbol("t")
    g = Symbol("g", positive=True)
    l = Symbol("l", positive=True)
    n = Symbol("n", positive=True)
    t = Symbol("s", positive=True)
    w = Symbol("w", positive=True)
    module = Symbol("E", positive=True)
    z = dynamicsymbols("z")

    scheme_name = "bridge_dmp.png"
    real_name = "leaf_spring.jpg"

    def __init__(
        self,
        m=None,
        k_beam=None,
        ivar=None,
        g=None,
        l=None,
        n=None,
        t=None,
        w=None,
        module=None,
        z=None,
        **kwargs
    ):

        if m is not None:
            self.m = m
        if k_beam is not None:
            self.k_beam = k_beam
        if g is not None:
            self.g = g
        if l is not None:
            self.l = l
        if n is not None:
            self.n = n
        if t is not None:
            self.t = t
        if w is not None:
            self.w = w
        if z is not None:
            self.z = z
        if module is not None:
            self.module = module

        self.qs = [self.z]
        self._init_from_components(**kwargs)

    #         self.force = Force(-F_0 * sin(Omega * ivar), pos1=z)

    @property
    def components(self):

        components = {}

        self._mass = MaterialPoint(self.m, self.z, qs=self.qs)
        self._spring = Spring(self.k_beam, self.z, qs=self.qs)
        self._gravitational_force = GravitationalForce(
            self.m, self.g, self.z, qs=self.qs
        )
        components["_mass"] = self._mass
        components["_spring"] = self._spring
        components["_gravitational_force"] = self._gravitational_force

        return components

    def set_equivalent_data(self):
        # E - Young's modulus, n - number of leaves, w - width of leaves, t - thickness of leaves, l - span; https://www.piping-designer.com/index.php/disciplines/mechanical/stationary-equipment/fastener/2486-leaf-spring-stiffness

        equivalent_data_dict = {
            self.k_beam: (8 * self.module * self.n * self.w * self.t**3)
            / (3 * self.l**3),
        }

        return equivalent_data_dict

    def get_default_data(self):
        default_data_dict = {
            self.m: 3000,
            self.module: 2.1 * 10**11,
            self.l: 0.9,
            self.g: 9.81,
            self.n: 8,
            self.t: 0.0135,
            self.w: 0.08,
        }

        return default_data_dict

    def get_real_data(self):

        real_data_dict = {
            self.m: 3000,
            self.module: 2.1 * 10**11,
            self.l: 0.9,
            self.g: 9.81,
            self.n: 8,
            self.t: 0.0135,
            self.w: 0.08,
        }

        return real_data_dict


class DDoFTrain(ComposedSystem):

    m = Symbol("m", positive=True)
    I = Symbol("I", positive=True)
    k = Symbol("k", positive=True)
    c = Symbol("c", positive=True)
    l_l = Symbol("l_l", positive=True)
    l_r = Symbol("l_r", positive=True)
    w = Symbol("w", positive=True)
    k_t = Symbol("k_t", positive=True)
    #     Omega1=Symbol('Omega1',postive=True)
    Omega2 = Symbol("Omega2", postive=True)
    S = Symbol("S", postive=True)
    #     A=Symbol('A',postive=True)
    lam = Symbol("lambda", postive=True)
    ivar = Symbol("t")
    s = Symbol("s", positive=True)
    v = Symbol("v", positive=True)
    d = Symbol("d", positive=True)
    h = Symbol("h", positive=True)
    l_0 = Symbol("l_0", positive=True)
    f = Symbol("f", positive=True)
    # Dynamic Symbols
    z = dynamicsymbols("z")
    phi = dynamicsymbols("varphi")

    def __init__(
        self,
        m=None,
        I=None,
        k=None,
        c=None,
        l_l=None,
        l_r=None,
        w=None,
        #                      Omega1=None,
        Omega2=None,
        S=None,
        #                      A=None,
        lam=None,
        z=None,
        phi=None,
        s=None,
        v=None,
        d=None,
        h=None,
        k_t=None,
        l_0=None,
        f=None,
        ivar=Symbol("t"),
        **kwargs
    ):

        if m is not None:
            self.m = m
        if I is not None:
            self.I = I
        if k is not None:
            self.k = k
        if c is not None:
            self.c = c
        if l_l is not None:
            self.l_l = l_l
        if l_r is not None:
            self.l_r = l_r
        if w is not None:
            self.w = w
        #             if Omega1 is not None: self.Omega1 = Omega1
        if Omega2 is not None:
            self.Omega2 = Omega2
        if S is not None:
            self.S = S
        #             if A is not None: self.A = A
        if lam is not None:
            self.lam = lam
        if z is not None:
            self.z = z
        if phi is not None:
            self.phi = phi
        if s is not None:
            self.s = s
        if v is not None:
            self.v = v
        if d is not None:
            self.d = d
        if h is not None:
            self.h = h
        if k_t is not None:
            self.k_t = k_t
        if l_0 is not None:
            self.l_0 = l_0
        if f is not None:
            self.f = f
        self.ivar = ivar

        self.qs = [self.z, self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._train_body = RigidBody2D(
            self.m, self.I, pos_lin=self.z, pos_rot=self.phi, qs=self.qs
        )(label="body")
        #             self._harmonic_force = Force(self.A*sin(self.Omega1*self.ivar),
        #                                           pos1=self.phi,
        #                                           qs=self.qs)(label='harmonic force')

        #             self._spring_left = Spring(self.k,
        #                                  pos1=self.z + self.phi*self.l_l,
        #                                        pos2=Heaviside(self.ivar-50)*self.S1*sin(self.Omega2*(self.ivar-self.t0)),
        #                                   qs=self.qs)(label='left spring')

        #             self._spring_right = Spring(self.k,
        #                                  pos1=self.z - self.phi*self.l_r,
        #                                         pos2=Heaviside(self.ivar-50)*self.S2*sin(self.Omega2*self.ivar),
        #                                   qs=self.qs)(label='right spring')

        #             self._damper_left = Damper(self.c,
        #                                        pos1= self.z + self.phi*self.l_l,
        #                                        pos2=Heaviside(self.ivar-50)*self.S1*sin(self.Omega2*(self.ivar-self.t0)),
        #                                        qs=self.qs) (label = 'left damper')

        #             self._damper_right = Damper(self.c,
        #                                        pos1= self.z - self.phi*self.l_r,
        #                                        pos2=Heaviside(self.ivar-50)*self.S2*sin(self.Omega2*self.ivar),
        #                                        qs=self.qs) (label = 'right damper')
        ####################################################### Wymuszneie zwykłe ##################################################################
        #             self._spring_left = Spring(self.k,
        #                                  pos1=self.z + self.phi*self.l_l ,
        #                                        pos2=self.S1*sin(self.Omega2*(self.ivar-self.s/self.v)), #i tak samo niżej W TŁUMIKACH TEŻ ;0
        #                                  qs=self.qs)(label='left spring')

        #             self._spring_right = Spring(self.k,
        #                                  pos1=self.z - self.phi*self.l_r,
        #                                         pos2=self.S1*sin(self.Omega2*(self.ivar)),
        #                                   qs=self.qs)(label='right spring')

        #             self._damper_left = Damper(self.c,
        #                                        pos1= self.z + self.phi*self.l_l,
        #                                        pos2=self.S1*sin(self.Omega2*(self.ivar-self.s/self.v)),
        #                                        qs=self.qs) (label = 'left damper')

        #             self._damper_right = Damper(self.c,
        #                                        pos1= self.z - self.phi*self.l_r,
        #                                        pos2=self.S1*sin(self.Omega2*(self.ivar)),
        #                                        qs=self.qs) (label = 'right damper')

        ####################################################### Wymuszneie zależne od prędkości ########################################################

        self._spring_left = Spring(
            self.k,
            pos1=self.z + self.phi * self.l_l,
            pos2=self.S * sin(2 * pi * self.v / self.d * (self.ivar - self.s / self.v)),
            qs=self.qs,
        )(label="left spring")

        self._spring_right = Spring(
            self.k,
            pos1=self.z - self.phi * self.l_r,
            pos2=self.S * sin(2 * pi * self.v / self.d * (self.ivar)),
            qs=self.qs,
        )(label="right spring")

        self._damper_left = Damper(
            self.c,
            pos1=self.z + self.phi * self.l_l,
            pos2=self.S * sin(2 * pi * self.v / self.d * (self.ivar - self.s / self.v)),
            qs=self.qs,
        )(label="left damper")

        self._damper_right = Damper(
            self.c,
            pos1=self.z - self.phi * self.l_r,
            pos2=self.S * sin(2 * pi * self.v / self.d * (self.ivar)),
            qs=self.qs,
        )(label="right damper")

        #######################################################TORSIONAL I HORIZONTAL########################################################

        #             self._torsional_spring = Spring(self.k_t,
        #                                  pos1=self.z/self.l_r,
        # #                                         pos2=self.S1*sin(2*pi*self.v/self.d*(self.ivar)),
        #                                   qs=self.qs)(label='torsional spring')
        #            self._horizontal_spring = Spring(self.k_t,
        #                                 pos1=((self.z +self.phi*self.l_l)**2+self.l_0**2)-self.l_0 ,
        #                                         pos2=self.S1*sin(2*pi*self.v/self.d*(self.ivar)),
        #                                  qs=self.qs)(label='torsional spring')

        components["_train_body"] = self._train_body
        components["_spring_left"] = self._spring_left
        components["_spring_right"] = self._spring_right
        components["_damper_left"] = self._damper_left
        components["_damper_right"] = self._damper_right
        #            components['_horizontal_spring']=self._horizontal_spring
        #             components['_harmonic_force'] = self._harmonic_force

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.z: r"przemieszczenie pionowe nadwozia",
            self.phi: r"przemieszczenie kątowe nadwozia",
            self.phi.diff(self.ivar): r"prędkość kątowa nadwozia",
            self.z.diff(self.ivar): r"prędkość pionowa nadwozia",
            self.phi.diff(self.ivar, self.ivar): r"przyspieszenie kątowe nadwozia",
            self.z.diff(self.ivar, self.ivar): r"przyspieszenie pionowe nadwozia",
            self.m: r"masa modelu połówkowego",
            self.I: r"moment bezwładności modelu połówkowego",
            self.k: r"sztywność zastępcza usprężynowania pierwszego i drugiego stopnia",
            self.c: r"tłumienie zastępcze zestawów tłumiących pierwszego i drugiego stopnia",
            self.l_r: r"odległość układu sprężysto-tłumiącego po stronie prawej od środka obrotu",
            self.l_l: r"odległość układu sprężysto-tłumiącego po stronie lewej od środka obrotu",
            self.v: r"pozioma prędkość pojazdu",
            self.d: r"długość nierówności powstającej na szynie",
            self.S: r"amplituda wymuszenia kinematycznego",
            self.s: r"odległość pomiędzy wózkami jezdnymi",
            self.ivar: r"czas",
            self.w: r"długość modelu połówkowego",
            self.h: r"wysokość modelu połówkowego",
        }

        return self.sym_desc_dict

    def get_default_data(self):

        default_data_dict = {
            self.lam: S.One * self.lam,
            self.k: S.One * self.k,
            self.c: S.One * self.k * self.lam,
            #             self.c: S.One * self.c,
            self.m: S.One * self.m,
            self.I: S.One * self.I,
            self.l_l: S.One * self.l_l,
            self.l_r: S.One * self.l_r,
            self.w: S.One * self.w,
            #            self.Omega1: S.One * self.Omega1,
            self.Omega2: S.One * self.Omega2,
            self.S: S.One * self.S,
            #            self.A: S.One * self.A,
            self.s: S.One * self.s,
            self.v: S.One * self.v,
            self.d: S.One * self.d,
            self.h: S.One * self.h,
            self.f: S.One * self.f,
        }
        return default_data_dict

    def get_numerical_data(self):

        default_data_dict = {
            self.c: 15000,
            #            self.lam: 0.01,
            self.m: 11025,
            self.k: 500000,
            self.l_l: 3.063,
            self.l_r: 2.663,
            self.w: 6.125,
            #             self.Omega1 : 3.14*10,
            self.Omega2: 3.14 * 10,
            self.S: 0.01,
            #             self.A : 0.1,
            self.s: 5.725,
            self.v: 16.6,  # 60 km/h prędkość przejazdów
            self.d: 30,  # Zakres z artykułu 40-200m
            self.h: 1.8,  # Środek cięzkości w środku geometrycznym
            self.I: self.m * (self.w**2 + self.h**2) / 12,
        }
        return default_data_dict

    def default_data_table(self):
        default_data_dict = {
            self.c: 15000,
            #            self.lam: 0.01,
            self.m: 11025,
            self.k: 500000,
            self.l_l: 3.063,
            self.l_r: 2.663,
            self.w: 6.125,
            #             self.Omega1 : 3.14*10,
            self.Omega2: 3.14 * 10,
            self.S: 0.01,
            #             self.A : 0.1,
            self.s: 5.725,
            self.v: 16.6,  # 60 km/h prędkość przejazdów
            self.d: 30,  # Zakres z artykułu 40-200m
            self.h: 1.8,  # Środek cięzkości w środku geometrycznym
            self.I: self.m * (self.w**2 + self.h**2) / 12,
        }
        default_data_dict[self.I] = (
            (self.m * (self.w**2 + self.h**2) / 12)
            .subs(DDoFTrain().get_numerical_data())
            .n(3)
        )
        return default_data_dict

    def unit_registry(self):

        unit_dict = {
            self.phi: ureg.rad,
            self.z: ureg.meter,
            self.z.diff(DDoFTrain.ivar, 2): ureg.meter / ureg.sec**2,
            self.phi.diff(DDoFTrain.ivar, 2): ureg.rad / ureg.sec**2,
            self.c: ureg.newton * ureg.sec / ureg.meter,
            self.m: ureg.kilogram,
            self.I: ureg.kilogram * ureg.meter**2,
            self.k: ureg.newton / ureg.meter,
            self.l_l: ureg.meter,
            self.l_r: ureg.meter,
            #            self.Omega1 : ureg.rad/ureg.sec,
            self.Omega2: ureg.rad / ureg.sec,
            self.S: ureg.meter,
            #            self.A : ureg.meter,
            self.lam: ureg.sec,
            self.s: ureg.meter,
            self.v: ureg.meter / ureg.sec,
            self.d: ureg.meter,
            self.h: ureg.meter,
            self.f: ureg.hertz,
            self.w: ureg.meter,
        }
        return unit_dict

    def resonance_velocities(self):

        self.v_z1 = self.natural_frequencies()[0] * self.d / (2 * pi)
        self.v_z2 = self.natural_frequencies()[3] * self.d / (2 * pi)
        display(
            self.v_z1.subs(self.get_numerical_data()),
            self.v_z2.subs(self.get_numerical_data()),
        )
        return (
            self.v_z1.subs(self.get_numerical_data()),
            self.v_z2.subs(self.get_numerical_data()),
        )
