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

from sympy import Symbol
from sympy.physics.mechanics import dynamicsymbols

from .elements import Capacitor, Inductor, Resistor


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
