"""
circuits.py – wspólne mix-iny R-L-C dla wszystkich modeli elektrycznych.

Zamiast powielać klasy „CircutRl/CircutRC/…” w wielu modułach,
definiujemy jedno źródło prawdy:

* BaseCircuit – abstrakcyjny kontener na elementy,
* CircuitRL   – rezystor + cewka,
* CircuitRC   – rezystor + kondensator,
* CircuitRLC  – rezystor + cewka + kondensator.

Każda klasa zwraca słownik komponentów zgodny z interfejsem dynpy –
można go bezpośrednio włączyć do ComposedSystem.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List

from sympy import *
from sympy.physics.mechanics import dynamicsymbols

from .elements import Capacitor, Inductor, Resistor


from ..mechanics.trolley import (
    ComposedSystem,
    NonlinearComposedSystem,
    base_frame,
    base_origin,
)

# ──────────────────────────────────
# 1.  Klasa bazowa
# ──────────────────────────────────
@dataclass
class BaseCircuit:
    _counter = itertools.count()  #  <── poprawiona nazwa

    def _next_q(self):
        """Zwraca kolejną współrzędną q_i(t)."""
        return dynamicsymbols(f"q_{next(self._counter)}")

    def components(self) -> Dict[str, object]:
        """Zwraca dict {nazwa: element}.  Implementują go klasy potomne."""
        raise NotImplementedError


# ──────────────────────────────────
# 2.  Gotowe mix-iny
# ──────────────────────────────────
@dataclass
class CircuitRL(BaseCircuit):
    R: Symbol = Symbol("R", positive=True)
    L: Symbol = Symbol("L", positive=True)

    def components(self) -> Dict[str, object]:
        q = self._next_q() if self.qs is None else self.qs[0]
        return {
            "res": Resistor(R=self.R, q=q, ivar=self.ivar),
            "ind": Inductor(L=self.L, q=q, ivar=self.ivar),
        }


@dataclass
class CircuitRC(BaseCircuit):
    R: Symbol = Symbol("R", positive=True)
    C: Symbol = Symbol("C", positive=True)

    def components(self) -> Dict[str, object]:
        q = dynamicsymbols(f"q_{id(self)}") if self.qs is None else self.qs[0]
        return {
            "res": Resistor(R=self.R, q=q, ivar=self.ivar),
            "cap": Capacitor(C=self.C, q=q, ivar=self.ivar),
        }


@dataclass
class CircuitRLC(BaseCircuit):
    R: Symbol = Symbol("R", positive=True)
    L: Symbol = Symbol("L", positive=True)
    C: Symbol = Symbol("C", positive=True)

    def components(self) -> Dict[str, object]:
        q = dynamicsymbols(f"q_{id(self)}") if self.qs is None else self.qs[0]
        return {
            "res": Resistor(R=self.R, q=q, ivar=self.ivar),
            "ind": Inductor(L=self.L, q=q, ivar=self.ivar),
            "cap": Capacitor(C=self.C, q=q, ivar=self.ivar),
        }


from dynpy.utilities.report import TikZPicture


class ElecricScheme(TikZPicture):
    def _scheme_desc(self):
        return r"""
\begin{scope}[american, every node/.style={font=\small}, line width=0.6pt]
  % --- bloki logiczne -------------------------------------------------
  \tikzset{block/.style={
      draw, rounded corners, minimum width=3cm, minimum height=1.2cm,
      font=\sffamily\small
  }}
  \node[block] (rpi) at (0,0)               {RPi\,4};
  \node[block, right=3.5cm of rpi] (stm)    {RPI PICO};

  % --- linia UART -----------------------------------------------------
  \draw (rpi.east) -- node[above]{UART} (stm.west);

  % --- linia PWM z R₂ przed silnikiem --------------------------------
  \coordinate (pwm) at ($(stm.east)+(0.9,0)$);
  \draw (stm.east) -- node[above]{GP26} (pwm)
        to[R,l=$R_2$,*-*] ++(2,0) coordinate (motin);

  % --- blok silnika DC -----------------------------------------------
  \node[block, right=0.8cm of motin] (drv) {DC Motor};
  \draw (motin) -- (drv.west);

  % --- rezystor R₁ do masy -------------------------------------------
  \draw (pwm) to[R,l=$R_1$] ++(0,-2) node[ground]{};

  % --- wyłącznik oraz bateria ----------------------------------------
  \draw (drv.east) -- ++(1.2,0)
        to[spst,l=$S$,*-*] ++(1.5,0)
        coordinate (bat+) |- ++(0,-3)
        to[battery1,l=18650] ++(0,-2) node[ground]{};
\end{scope}
"""



class CircutRL(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RL system
    """

    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'

    resistance = Symbol("R", positive=True)
    inductance = Symbol("L", positive=True)
    q0 = dynamicsymbols("q_c")
    qs = dynamicsymbols("qs")
    frame = Symbol("frame", positive=True)
    ivar = Symbol("t")

    def __init__(
        self,
        resistance=None,
        inductance=None,
        ivar=None,
        q0=None,
        qs=None,
        frame=None,
        z=None,
        **kwargs,
    ):

        if resistance is not None:
            self.resistance = resistance
        if inductance is not None:
            self.inductance = inductance
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if q0 is not None:
            self.q0 = q0
        if qs is not None:
            self.qs = qs
        if frame is not None:
            self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.resistor = Resistor(
            self.resistance, self.q0, qs=self.qs, ivar=self.ivar, frame=base_frame
        )("resistor")
        self.inductor = Inductor(
            self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("inductor")

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r"resistance of resistor",
            self.inductance: r"inductance of inductor",
        }

        return self.sym_desc_dict


class CircutRC(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RC system
    """

    scheme_name = "thevenincircuit.png"
    real_name = "liioncell.PNG"

    resistance = Symbol("R", positive=True)
    capacity = Symbol("C", positive=True)
    q0 = dynamicsymbols("q_c")
    qs = dynamicsymbols("qs")
    frame = Symbol("frame", positive=True)
    ivar = Symbol("t")

    def __init__(
        self,
        resistance=None,
        capacity=None,
        ivar=None,
        q0=None,
        qs=None,
        frame=None,
        z=None,
        **kwargs,
    ):

        if resistance is not None:
            self.resistance = resistance
        if capacity is not None:
            self.capacity = capacity
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if q0 is not None:
            self.q0 = q0
        if qs is not None:
            self.qs = qs
        if frame is not None:
            self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.resistor = Resistor(
            self.resistance, self.q0, qs=self.qs, ivar=self.ivar, frame=base_frame
        )("resistor")
        self.capacitor = Capacitor(
            self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("capacitor")

        components["resistor"] = self.resistor
        components["capacitor"] = self.capacitor

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r"resistance of resistor",
            self.capacity: r"capacity of capacitor",
        }

        return self.sym_desc_dict


class CircuitRLC(ComposedSystem):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """

    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'

    resistance = Symbol("R", positive=True)
    inductance = Symbol("L", positive=True)
    capacity = Symbol("C", positive=True)
    q0 = dynamicsymbols("q")
    frame = Symbol("frame", positive=True)
    ivar = Symbol("t")

    def __init__(
        self,
        resistance=None,
        inductance=None,
        capacity=None,
        ivar=None,
        q0=None,
        frame=None,
        z=None,
        **kwargs,
    ):

        if resistance is not None:
            self.resistance = resistance
        if inductance is not None:
            self.inductance = inductance
        if capacity is not None:
            self.capacity = capacity
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if q0 is not None:
            self.q0 = q0
        if frame is not None:
            self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self.resistor = Resistor(
            self.resistance, self.q0, qs=self.qs, ivar=self.ivar, frame=base_frame
        )("resistor")
        self.inductor = Inductor(
            self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("inductor")
        self.capacitor = Capacitor(
            self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("capacitor")

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor
        components["capacitor"] = self.capacitor

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.resistance: r"resistance of resistor",
            self.inductance: r"inductance of inductor",
            self.capacity: r"capacity of capacitor",
        }

        return self.sym_desc_dict


class CircuitRLCWithPWM(CircuitRLC):
    """
    A class that determines the equation of an electrical circuit in an RLC system
    """

    # scheme_name = 'engine.png'
    # real_name = 'engine_real.PNG'

    resistance = Symbol("R", positive=True)
    inductance = Symbol("L", positive=True)
    capacity = Symbol("C", positive=True)
    U = Symbol("U", positive=True)
    eps = Symbol("varepsilon", postive=True)
    rho = Symbol("rho", positive=True)
    delta = Symbol("delta", positive=True)
    T = Symbol("T", positive=True)
    q0 = dynamicsymbols("q")
    omega = Symbol("omega", positive=True)
    ivar = Symbol("t")

    def __init__(
        self,
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
        **kwargs,
    ):

        if resistance is not None:
            self.resistance = resistance
        if inductance is not None:
            self.inductance = inductance
        if capacity is not None:
            self.capacity = capacity
        if U is not None:
            self.U = U
        if eps is not None:
            self.eps = eps
        if rho is not None:
            self.rho = rho
        if delta is not None:
            self.delta = delta
        if T is not None:
            self.T = T
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if q0 is not None:
            self.q0 = q0
        if frame is not None:
            self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        self.r_var = self.resistance * (1 + self.eps * self.rho)

        self.resistor = Resistor(
            self.r_var, self.q0, qs=self.qs, ivar=self.ivar, frame=base_frame
        )("resistor")
        self.inductor = Inductor(
            self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("inductor")
        self.capacitor = Capacitor(
            self.capacity / self.delta / self.eps,
            self.q0,
            ivar=self.ivar,
            qs=self.qs,
            frame=base_frame,
        )("capacitor")
        self.voltage = VoltageSource(self.U, self.q0, ivar=self.ivar, qs=self.qs)

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor
        components["capacitor"] = self.capacitor
        components["voltage"] = self.voltage

        return components

    def symbols_description(self):

        self.sym_desc_dict = {
            self.resistance: r"resistance",
            self.inductance: r"inductance",
            self.capacity: r"capacitance",
            self.U: r"source voltage",
            self.eps: r"small parameter",
            self.rho: r"formula representing resistance variation",
            self.q0: r"electric charge",
            self.q0.diff(self.ivar): r"current",
            self.q0.diff(self.ivar, self.ivar): r"current rate of change",
        }

        return self.sym_desc_dict

    def para(self, angle=2 * pi):
        para_series_heave = 10 * Heaviside(sin(2 * pi / (self.T) * self.ivar), 0.5)
        new_eq = self.subs(self.rho, para_series_heave)
        return new_eq

    def wave(self):

        wave1 = (1 * (sin(2 * pi * self.ivar / self.T)) + 2.0) * (
            1 / 2 + 1 / 2 * sign(sin(2 * pi * self.ivar / self.T))
        )
        wave2 = (1 * (-sin(2 * pi * self.ivar / self.T)) - 3.0) * (
            1 / 2 + 1 / 2 * sign(-sin(2 * pi * self.ivar / self.T))
        )
        waves = wave1 + wave2
        #         wave1=(100*(sin(self.ivar/self.T))+1000)*Heaviside(sin(self.ivar/self.T))
        #         wave2=(100*(-sin(self.ivar/self.T)))*Heaviside(-sin(self.ivar/self.T))
        waves = wave1 + wave2
        new_eq = self.subs(self.rho, waves)

        return new_eq

    def rect(self, no=6, numerical=False):
        t = self.ivar
        omg = self.omega
        rectangular = 5 * (
            (1 / 2 + 1 / 2 * sign(sin(2 * pi * self.ivar / self.T))) - S.Half
        )
        #         rectangular_func = lambdify((t, T), rectangular)
        #         rectangular_values = rectangular_func(t_values, T_value)

        #         trig=sum([Heaviside(omg*t-1) + 0* Heaviside(omg*t-2)  for ind in range(no)])
        #         new_eq=self.subs(self.rho,trig)
        new_eq = self.subs(self.rho, rectangular)

        return new_eq

    def approx_rect(self, no=6, numerical=False):

        if numerical is True:
            amps_list = [
                2.27348466531425,
                0,
                0.757408805199249,
                0,
                0.453942816897038,
                0,
                0.323708002807428,
                0,
                0.25121830779797,
                0,
                0.204977919963796,
                0,
                0.172873394602606,
                0,
                0.149252079729775,
                0,
                0.131121653619234,
                0,
                0.116749954968057,
                0,
            ]
        else:
            amps_list = symbols(f"a_0:{no}")
        rectangular_approx = sum(
            [
                amp * 2 / sqrt(2) * sin((2 * (no) + 1) * 2 * pi * self.ivar / self.T)
                for no, amp in enumerate(amps_list[0:])
            ]
        )
        new_eq = self.subs({self.rho: rectangular_approx})

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

    resistance = Symbol("R", positive=True)
    inductance = Symbol("L", positive=True)
    capacity = Symbol("C", positive=True)
    R_min = Symbol("R_min", positive=True)
    U = Symbol("U", positive=True)
    eps = Symbol("varepsilon", postive=True)
    rho = Symbol("rho", positive=True)
    delta = Symbol("delta", positive=True)
    T = Symbol("T", positive=True)
    q0 = dynamicsymbols("q")
    omega = Symbol("omega", positive=True)
    ivar = Symbol("t")

    def __init__(
        self,
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
        **kwargs,
    ):

        if resistance is not None:
            self.resistance = resistance
        if R_min is not None:
            self.R_min = R_min
        if inductance is not None:
            self.inductance = inductance
        if capacity is not None:
            self.capacity = capacity
        if U is not None:
            self.U = U
        if eps is not None:
            self.eps = eps
        if rho is not None:
            self.rho = rho
        if delta is not None:
            self.delta = delta
        if T is not None:
            self.T = T
        if ivar is not None:
            self.ivar = ivar
        if z is not None:
            self.z = z
        if q0 is not None:
            self.q0 = q0
        if frame is not None:
            self.frame = frame

        self.qs = [self.q0]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}
        self.r_var = self.resistance * (1 + self.eps * self.rho)

        self.resistor = Resistor(
            self.R_min + self.resistance * Heaviside(sin(2 * pi * self.ivar / self.T)),
            self.q0,
            qs=self.qs,
            ivar=self.ivar,
            frame=base_frame,
        )("resistor")
        self.inductor = Inductor(
            self.inductance, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("inductor")
        self.capacitor = Capacitor(
            self.capacity, self.q0, ivar=self.ivar, qs=self.qs, frame=base_frame
        )("capacitor")
        self.voltage = VoltageSource(self.U, self.q0, ivar=self.ivar, qs=self.qs)

        components["resistor"] = self.resistor
        components["inductor"] = self.inductor
        components["capacitor"] = self.capacitor
        components["voltage"] = self.voltage

        return components

    def symbols_description(self):

        self.sym_desc_dict = {
            self.resistance: r"amplitude of the variable part of the resistance",
            self.R_min: r"constant value of the resistance",
            self.inductance: r"inductance",
            self.capacity: r"capacitance",
            self.U: r"source voltage",
            self.q0: r"electric charge",
            self.q0.diff(self.ivar): r"current",
            self.q0.diff(self.ivar, self.ivar): r"current rate of change",
        }

        return self.sym_desc_dict
