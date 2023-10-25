from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, re, pi, latex,
                   dsolve, solve, fraction, factorial, Add, Mul, exp, zeros, shape,
                   numbered_symbols, integrate, ImmutableMatrix,Expr,Dict,Subs,Derivative,Dummy,
                   lambdify, Pow, Integral, init_printing)

from sympy.matrices.matrices import MatrixBase
from numbers import Number

###  exemplary comment
from sympy.physics.mechanics import dynamicsymbols, init_vprinting
from sympy.physics.vector.printing import vpprint, vlatex
from ..utilities.components.ode import en as ode

import pandas as pd
from ..utilities.adaptable import *

ODE_COMPONENTS_LIST = [
            ode.ODESystemComponent,
            ode.ODEInitComponent,
            ode.ODESystemCodeComponent,
            ode.VariablesComponent,
            ode.HomEquationComponent,
            ode.HomEquationCodeComponent,
            ode.GoverningEquationComponent,
            ode.FundamentalMatrixComponent,
#             ode.ODECharecteristicPolynomialComponent,
            ode.ODECharecteristicPolynomialCodeComponent,
            ode.GeneralSolutionComponent
]

class CommonFactorDetector:
    
    _default_symbol = Symbol('\\tau^2',positive=True)
    
    def __init__(self,fundamental_matrix,default_symbol=None):
        self._fundametnal_matrix  = fundamental_matrix
        
        if default_symbol is not None: self._default_symbol = default_symbol
            
        
    def _const_spotter(self):
        
        
        nums_list = [list(fraction(elem)[0].atoms(Symbol)) for elem in self._fundametnal_matrix]
        denoms_list = [list(fraction(elem)[1].atoms(Symbol)) for elem in self._fundametnal_matrix]
       
        num_common=set(sum(nums_list,[]))
        for elem in nums_list:
            #print(elem,elem != set())
            if set(elem) != set():
                
                num_common &= set(elem)

                
        denom_common=set(sum(denoms_list,[]))
        for elem in denoms_list:
            #print(elem,elem != set())
            if set(elem) != set():
                
                denom_common &= set(elem)
        
        
        
        
        return (list(num_common)+[None])[0],(list(denom_common)+[None])[0]


    def subs_dict(self):
        
        num, denom  = self._const_spotter()
        
        
        
        if num is not None and denom is not None:
            return {num: denom / self._default_symbol}
        else:
            return {}
    
    def callback_dict(self):
        
        num, denom  = self._const_spotter()
        
        if num is not None and denom is not None:
            return {self._default_symbol: denom / num}
        else:
            return {}
        
        
class TwoSystemsSimulation:

    def __init__(self,odes, t_span=None):
        self._odes = odes
        self._t_span = t_span
        
        
    @property
    def odes(self):

        odes = self._odes

        if not isinstance(odes,tuple):
            odes = (odes,odes)

        return odes

    def numerized(self):
        return self
    
    def compute_solution(self,t_span,ic_list=[0.0,0.0,0.0,0.0]):

        from .linear import FirstOrderLinearODESystemWithHarmonics
        
        ode1 , ode2 = self.odes

        t_span0 = t_span[0]
        t_span1 = t_span[1]

        fode_1=ode1.as_first_ode_linear_system()
        fodes_1=FirstOrderLinearODESystemWithHarmonics(fode_1,fode_1.dvars)
#         ode2.as_first_ode_linear_system()
#         fodes_2=FirstOrderLinearODESystemWithHarmonics(ode2,ode2.dvars)

        fode_2=ode2.as_first_ode_linear_system()
        fodes_2=FirstOrderLinearODESystemWithHarmonics(fode_2,fode_2.dvars)

        ana_sol_sym_1=fodes_1.solution.with_ics(ic_list)
        ana_sym_wyn_1=pd.DataFrame(ana_sol_sym_1.numerized().compute_solution(t_span=t_span0))

        ana_sol_sym_2=fodes_2.solution.with_ics(sol0=ana_sol_sym_1, ivar0=t_span1[0])
        ana_sym_wyn_2=pd.DataFrame(ana_sol_sym_2.numerized().compute_solution(t_span=t_span1))

        return TimeDataFrame(pd.concat([ana_sym_wyn_1,ana_sym_wyn_2]))
        
        