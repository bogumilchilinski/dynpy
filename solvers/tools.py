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
        
        
        
        
        return (list(num_common)+[S.One])[0],(list(denom_common)+[S.One])[0]


    def subs_dict(self):
        
        num, denom  = self._const_spotter()
        
        return {num: denom / self._default_symbol}