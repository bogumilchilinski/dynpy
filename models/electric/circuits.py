from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from sympy import Symbol
from sympy.physics.mechanics import dynamicsymbols

# elementy bazowe z Twojego modułu elements.py
from .elements import Capacitor, Inductor, Resistor, VoltageSource


# ---------- 1. klasa bazowa -------------------------------------------------
@dataclass
class BaseCircuit:
    """
    Abstrakcyjny kontener na elementy R-L-C.
    Konkretny obwód nadpisuje tylko metodę `components()`.
    """

    ivar: Symbol = Symbol("t")  # zmienna niezależna
    qs: List[Symbol] | None = None  # lista współrzędnych uog.

    def components(self) -> Dict[str, object]:  # <- musi zwrócić dict elementów
        raise NotImplementedError


# ---------- 2. gotowe mix-iny ----------------------------------------------
@dataclass
class CircuitRL(BaseCircuit):
    R: Symbol = Symbol("R", positive=True)
    L: Symbol = Symbol("L", positive=True)

    def components(self):
        q = dynamicsymbols(f"q_{id(self)}") if self.qs is None else self.qs[0]
        return {
            "res": Resistor(R=self.R, q=q, ivar=self.ivar),
            "ind": Inductor(L=self.L, q=q, ivar=self.ivar),
        }


@dataclass
class CircuitRC(BaseCircuit):
    R: Symbol = Symbol("R", positive=True)
    C: Symbol = Symbol("C", positive=True)

    def components(self):
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

    def components(self):
        q = dynamicsymbols(f"q_{id(self)}") if self.qs is None else self.qs[0]
        return {
            "res": Resistor(self.R, q, ivar=self.ivar),
            "ind": Inductor(self.L, q, ivar=self.ivar),
            "cap": Capacitor(self.C, q, ivar=self.ivar),
        }
