from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Derivative, Expr)
from numbers import Number
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex

from ..elements import MaterialPoint, Damper, Spring




base_frame=ReferenceFrame('N')
base_origin=Point('O')

