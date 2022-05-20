import importlib
from . import dynamics as dyn
importlib.reload(dyn)
from . import dynamics as dyn

from sympy.physics.vector import dynamicsymbols
from sympy import *
from sympy.physics.mechanics import *
import sympy as sym
import numpy as np
import pylatex
from pylatex import Command, NoEscape
import pint
from dynpy.models.sdof import SpringMassSystem, DampedSpringMassSystem
ureg = pint.UnitRegistry()

mechanics_printing()



t=Symbol('t') #independent variable - time

class LEWEL16_simple(SpringMassSystem):
    scheme_name = 'train_sdof.png'
    real_name = 'lewel16.jpg'

class LEWEL16_simple_damped(DampedSpringMassSystem):
    scheme_name = 'train_sdof_damped.png'
    real_name = 'lewel16.jpg'
