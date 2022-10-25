from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, re, pi, latex,
                   dsolve, solve, fraction, factorial, Add, Mul, exp, zeros, shape,
                   numbered_symbols, integrate, ImmutableMatrix,Expr,Dict,Subs,Derivative,Dummy,
                  lambdify)

from sympy.matrices.matrices import MatrixBase

from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from ..utilities.timeseries import TimeSeries, TimeDataFrame

from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7,TR6, TR5, TR3
from collections.abc import Iterable
from sympy.solvers.ode.systems import linodesolve
from functools import cached_property

from .numerical import OdeComputationalCase



class MultivariableTaylorSeries(Expr):
    """_summary_

    Args:
        Expr (_type_): _description_
    """
    def __new__(cls,expr, variables,*args, n=2, x0=None):
        """_summary_

        Parameters
        ----------
        expr : _type_
            _description_
        variables : _type_
            _description_
        n : int, optional
            _description_, by default 2
        x0 : _type_, optional
            _description_, by default None

        Returns
        -------
        _type_
            _description_
        """
        
        obj=super().__new__(cls,expr,variables,*args)
        obj._vars = variables
        obj._order = n
        obj._op_point = x0
        
        obj._expr_symbol = None
        
        return obj

    @property
    def order(self):
        return self._order
    
    def _set_default_op_point(self):
        '''
        It sets 0 as op_point if x0 (look __new__) is not provided. For lack of op_point, the result is the same as MacLaurina series 
        '''
        
        if self._op_point is None:
            self._op_point = {coord:0 for coord in self._vars}
        
        return self._op_point
    
    def _args_shifted(self,*args):
        
        self._set_default_op_point()
        
        args_shifted = {
            arg: arg - arg_shift
            for arg, arg_shift in self._op_point.items()
        }
        
        return args_shifted

    def _diff_orders_dict(self):
            
        order_max = self._order
        args = self._vars
        args_shifted = self._args_shifted()
        
            
        diff_orders_list = sum([
            list(itools.combinations_with_replacement(args, order))
            for order in range(1, order_max + 1, 1)
        ], [])
        
        
        diff_orders_dict = {
            comp: (sym.Mul(*comp).subs(args_shifted) / sym.Mul(*[
                sym.factorial(elem) for elem in sym.Poly(sym.Mul(*comp), *args).terms()[0][0]
            ])).doit()
            for comp in diff_orders_list
        }

        return diff_orders_dict

    
    def _diff_symbols_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        expr=self.args[0]
        op_point=self._op_point
        
        
        
        return {S.Zero:Subs(expr,list(op_point.keys()),list(op_point.values())),**{args_tmp:Subs(Derivative(expr,*args_tmp,evaluate=False),list(op_point.keys()),list(op_point.values())) 
            for args_tmp, poly in diff_orders_dict.items()}}

    def _diff_expr_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        expr=self.args[0]
        op_point=self._op_point
        
        return {S.Zero:expr.subs(op_point),**{args_tmp:expr.diff(*args_tmp).subs(op_point)
            for args_tmp, poly in diff_orders_dict.items()}}
    
    def _components_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        derivatives_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return {
                Subs(expr,list(op_point.keys()),list(op_point.values())):expr.subs(op_point),
                **{derivatives_dict[args_tmp] : expr.diff(*args_tmp).subs(op_point) for args_tmp, poly in diff_orders_dict.items()}
        }
        
    def _series(self):
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return expr.subs(op_point).doit()+Add(*[expr.diff(*args_tmp).subs(op_point).doit() * poly for args_tmp, poly in diff_orders_dict.items()],evaluate=False)
    
    def _symbolic_sum(self):
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return Subs(expr,list(op_point.keys()),list(op_point.values()))+Add(*[Mul(diff_dict[args_tmp]  ,poly,evaluate=True) for args_tmp, poly in diff_orders_dict.items()],evaluate=False)
    
    def _latex(self,*args):
        
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return '+'.join([latex(Mul(diff_dict[args_tmp]  ,poly,evaluate=True)) for args_tmp, poly in diff_orders_dict.items()])

    
    def calculation_steps(self,expr_symbol=None,form=None):

        obj = self
        
        if expr_symbol is None:
            obj._expr_symbol = self.args[0]
            
        obj_sym = self.__class__(expr_symbol,self._vars, n=self._order, x0=self._op_point)
        
        expr_dict=(self._diff_expr_dict())
        diffs_dict=(obj_sym._diff_symbols_dict())
        

        
        
        return [Eq(diffs_dict[key],expr_dict[key].doit())   for  key  in diffs_dict.keys()]
    
    def __str__(self,*args):
        return (self.args[0]).__str__()
    
    def __repr__(self,*args):
        return (self.args[0]).__repr__()

    
    
    
    def doit(self,**hints):
        return self._series().doit()




class AnalyticalSolution(Matrix):
    def __new__(cls, data,rhs=None ,evaluate=True, **options):

        
        
        if isinstance(data,(dict,Dict)):
            return cls.from_dict(data,**options)
        if isinstance(data,Eq):
            return cls.from_eq(data,**options)
        if isinstance(data,Matrix) and isinstance(data[0],Eq):
            return cls.from_eqns_matrix(data,**options)
        else:
            return cls._constructor(data,rhs ,evaluate=evaluate, **options)

        return obj



    
    @classmethod
    def _constructor(cls, lhs,rhs=None ,evaluate=True, **options):
        
        if not isinstance(lhs,Iterable): 
            lhs = Matrix([lhs])
            
        # it should be reimplemented with @property if None as odes_rhs is spotted
        if rhs is None:
            rhs = 0*lhs
        elif not isinstance(rhs,Iterable):
            rhs = Matrix([rhs])

        obj = super().__new__(cls,rhs,evaluate=evaluate, **options)

        obj._lhs=lhs

        return obj
    
    @classmethod
    def from_dict(cls,dictionary,**options): # działa
        
        return cls._constructor( Matrix(list(dictionary.keys())),Matrix(list(dictionary.values())) ,**options   )

    
    @classmethod
    def from_eq(cls,matrix_eq,**options): #
        if isinstance(matrix_eq,Eq):
            
            eq = matrix_eq
        else:
            print ('TypeError: matrix_eq is not Eq')
            return None
        
        obj = cls._constructor(eq.lhs,eq.rhs ,evaluate=True ,**options   )
        
        return obj

    @classmethod
    def from_eqns_matrix(cls,eqns_matrix,**options):

        dvars = []
        values = []
        if isinstance(eqns_matrix,Matrix):
            matrix = eqns_matrix
        else:
            print ('TypeError: eqns_matrix is not matrix')
            return None
        
        for eqs in shape(matrix):
            if isinstance(matrix[eqs-1],Eq):
                eq = matrix[eqs-1]
            else:
                print ('TypeError: eqprestion is not Eq')
                return None

            dvar = eq.lhs
            dvars.append(dvar)
            
            value = eq.rhs
            values.append(value)
            
        obj = cls._constructor(Matrix(dvars),Matrix(values),evaluate=True ,**options)

        return obj
    
    
    def subs(self,*args,**kwargs):
        

        
        obj = super().subs(*args,**kwargs)
        obj._lhs=self._lhs
        
        return obj
    
    
    def __add__(self,other):
        
        if isinstance(other,AnalyticalSolution):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
        obj = super().__add__(other)
        obj._lhs=self._lhs
        
        return obj

    
    def __rsub__(self,other):

        if isinstance(other,AnalyticalSolution):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])

        obj = super().__rsub__(other)
        obj._lhs=self._lhs

        return obj


    def __sub__(self,other):

        if isinstance(other,AnalyticalSolution):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])

        obj = super().__sub__(other)
        obj._lhs=self._lhs

        return obj


    def __mul__(self,other):

        obj = super().__mul__(other)
        obj._lhs=self._lhs

        return obj

    def doit(self,**hints):

        obj = super().doit(**hints)
        obj._lhs=self._lhs

        return obj
    
    def _eval_applyfunc(self, f):
        
        obj = super()._eval_applyfunc(f)
        obj._lhs=self._lhs
        
        return obj
    
    def __call__(self,t,params={}):
        
        solution = (self.nth_order_solution(order).rhs.subs(self.extra_params).subs(self.params_values))
        
        solution = solution.applyfunc(lambda x: x.subs(self.extra_params).subs(self.params_values))

#         print('solution after extra params')
#         display(-(self.nth_order_solution(order).rhs.subs(self.extra_params)))
            
        if isinstance(ivar, Symbol):

            ics_dict=self._ic_from_sol(order=order,formula=True)
            display(ics_dict)

            return solution.subs(ics_dict).subs(self.extra_params).subs(self.ivar, ivar)
        else:

            ics_dict=self._ic_from_sol(order=order,formula=False)
#             display(Dict(ics_dict).subs(self.params_values))
            
            ics_symbols_dict={sym:val for sym, val in zip(self.ics_symbols,self.ics)  }
            

            nth_order_solution_fun = lambdify(self.ivar, (solution.subs(ics_dict).subs(self.extra_params).subs(ics_symbols_dict)).subs(self.params_values).n(), 'numpy')

            #return nth_order_solution_fun(ivar)
            solution  = TimeDataFrame(data={
                dvar: data[0]
                for dvar, data in zip(self.dvars, nth_order_solution_fun(ivar))
            },
                                 index=ivar)
            
            
            for dvar in self.dvars:
                solution[dvar.diff(self.ivar,1)]=solution[dvar].gradient()
            
            for dvar in self.dvars:
                solution[dvar.diff(self.ivar,2)]=solution[dvar.diff(self.ivar,1)].gradient()
            
            return solution
    
    
    @property
    def lhs(self):
        return self._lhs

    
    @property
    def rhs(self):
        return Matrix(list(self))
    

    def as_matrix(self):
        return Matrix( self.lhs-self.rhs )
    
    def as_eq(self):
        return Eq(self.lhs,self.rhs)
    
    def as_iterable(self):

        return [(lhs,comp) for lhs,comp  in zip(self._lhs,self)]
        
    
    def as_dict(self):

        return Dict({lhs:rhs  for  lhs,rhs in self.as_iterable()})
    
    def as_eq_list(self):

        return [ Eq(lhs,comp,evaluate=False) for lhs,comp  in zip(self._lhs,self)]
    
    @property
    def _lhs_repr(self):
        return self.lhs
    
    def __repr__(self):

        return f'{self._lhs_repr} = {self.rhs}'
    
    
    def _latex(self,*args):

        return latex(Eq(self._lhs_repr,self.rhs,evaluate=False))

    def numerized(self,t_span=[],parameters={}):
        
        ivar = Symbol('t')
        
        solution = self.rhs.subs(parameters)
        
        sol_func = lambdify(ivar, solution, 'numpy')
        
        numerized_data = (sol_func(t_span))
        
        dvars = self.lhs
        
        numerized_sol  = TimeDataFrame(data={dvar:data[0] for dvar,data in zip(dvars,numerized_data)},index=t_span)
        numerized_sol.index.name = ivar
        
        return numerized_sol
    
    
class ODESolution(AnalyticalSolution):
    
    def constant(self,initial_conditions_matrix, msm_params):
        
        self.initial_conditions_matrix = initial_conditions_matrix
        self.msm_params = msm_params

        constant_dict=solve((self.subs(t,0).subs(eps,0).subs(msm_params).doit().applyfunc(lambda row: row.doit()).rhs-initial_conditions_matrix).subs(msm_params).n())
        
        return constant_dict
    
    def generate_ana_sol(self, msm_params, t_span, l_w, l_wnum=None):
        self.msm_params = msm_params
        self.t_span = t_span
        self.l_w = l_w
        self.l_wnum = l_wnum
        
        if isinstance(l_wnum,None):
            l_wnum = l_w
        
        ana_sol = self.subs(msm_params).subs(l_w,l_wnum).subs(self.constant()).n().numerized(t_span)
        
        return ana_sol

    def generate_num_sys(self, msm_params, t_span, l_w, l_wnum=None):
        self.msm_params = msm_params
        self.t_span = t_span
        self.l_w = l_w
        self.l_wnum = l_wnum
        
        if isinstance(l_wnum,None):
            l_wnum = l_w
        
        num_sys=self.subs(msm_params).subs(l_w,l_wnum).numerized().compute_solution(t_span ,self.generate_ana_sol().iloc[0,:].to_numpy(),params_values=msm_params)
        
        return num_sys

class ODESystem(AnalyticalSolution):
    
    _ivar = Symbol('t')
    _parameters = None
    _ode_order = 1
    
    @classmethod
    def from_ode_system(cls,ode_system):

        sys = ode_system
        odes=sys.lhs
        odes_rhs = sys.rhs
        dvars = sys.dvars
        ivar = sys.ivar
        parameters = sys._parameters
        ode_order = sys.ode_order

        return cls._constructor(odes , dvars, odes_rhs , ivar,ode_order=ode_order ,parameters=parameters)

    

    
    @classmethod
    def set_default_parameters(cls,parameters):
        cls._parameters = parameters
        return cls


    
    
    def __new__(cls, odes,dvars,odes_rhs=None , ivar=None ,ode_order=None,evaluate=True, parameters = None, **options):
        


        if ivar is None:
            ivar = cls._ivar
            
        if ode_order is None:
            ode_order = cls._ode_order

        
        if isinstance(odes,dict):
            return cls.from_dict(odes,**options)
        if isinstance(odes,Eq):
            return cls.from_eq(odes,**options)
        if isinstance(odes,MatrixBase) and odes.jacobian(dvars.diff(ivar,ode_order)).det()==0:
            return cls.from_rhs(odes_rhs=odes,dvars=dvars,ivar=ivar,ode_order=ode_order,parameters = parameters)
        else:
            return cls._constructor(odes,dvars,odes_rhs=odes_rhs , ivar=ivar,ode_order=ode_order ,evaluate=evaluate, parameters = parameters, **options)




    @classmethod
    def from_rhs(cls,odes_rhs,dvars,ivar=None,ode_order=1,parameters=None):
        
        if ivar is None:
            ivar = cls._ivar
        
        vels = dvars.diff(ivar)
        

        return cls._constructor( odes=vels , dvars=dvars  , odes_rhs = odes_rhs   ,ivar=ivar, ode_order=ode_order,parameters=parameters)
        
        
    @classmethod
    def _constructor(cls,odes,dvars,odes_rhs=None , ivar=None,ode_order=None ,evaluate=True, parameters = None, **options):
        

        if not isinstance(dvars,Iterable):
            dvars = Matrix([dvars])
        
        obj = super()._constructor(odes,odes_rhs,evaluate=evaluate, **options)
        
        obj._parameters = parameters
        obj._dvars = dvars
        obj._const_list = []
        
        
        if ivar is not None:
            obj._ivar = ivar

        if ode_order is not None:
            obj._ode_order = ode_order

            
        return obj
    
    @cached_property
    def ivar(self):
        return self._ivar
    
    @cached_property
    def dvars(self):
        return self._dvars
    
    @cached_property
    def _fode_dvars(self):
        return Matrix([self._dvars.diff(self.ivar,order) for order  in range(self.ode_order)])
    
    @cached_property
    def _highest_diff(self):
        return self.dvars.diff(self.ivar,self.ode_order)
    
    @cached_property
    def odes_rhs(self):
        
        diffs = self._highest_diff
        
        if diffs == self.lhs:
            return self.rhs
        else:
            diff_coeffs_mat=self.lhs.jacobian(diffs)
            
            return diff_coeffs_mat.inv()*(self.rhs  - self.lhs.subs({elem:0 for elem in diffs }))

                                          
    @cached_property
    def parameters(self):
        
        return self._parameters

    @cached_property    
    def ode_order(self):
        
        return self._ode_order
    
    @cached_property
    def odes(self):
        return Matrix([Add(lhs,-rhs,evaluate=False) if lhs-rhs == 0 else lhs-rhs   for lhs,rhs   in zip(self.lhs,self.rhs) ])


    def as_eq(self):
        return Eq(self.odes,self.dvars*0,evaluate=False)
    
    
    @cached_property
    def _lhs_repr(self):
        #return (self.lhs).applyfunc(lambda obj: Derivative(obj,self.ivar,evaluate=False)  )
        return (self.lhs).applyfunc(lambda obj: obj.diff(self.ivar)  )

    
    
    def _latex(self,*args):

        
        return f'{latex(self.as_eq())}~~for~~{latex(self.dvars)}' 

    def as_matrix(self):
        return Matrix(self._lhs_repr - self.rhs) 
    

    def as_first_ode_linear_system(self):
        
        
        return FirstOrderLinearODESystem.from_ode_system(self)

    def approximated(self, n=3, x0=None, op_point=False, hint=[], label=None):
        """
        Returns approximated N-th order function calculated with Taylor series method as an instance of the class
        """

        # print('x0',x0)
        if not x0:
            x0 = {coord: 0 for coord in self.dvars}

        # #display(self._op_points(hint=hint, subs=True))
        # if op_point:
        #     x0.update(self._op_points(hint=hint, subs=True)[0])
        #     #print('current op')
        #     #display(self._op_points(hint=hint, subs=True))

        
        ivar = self.ivar
        diff_matrix = self.lhs.jacobian(self.dvars.diff(ivar))
        
        
        
        rhs_eqns=diff_matrix.inv()*self.rhs
        
        lin_eqns = Matrix([MultivariableTaylorSeries(ode, self.dvars, n=n, x0=x0).doit() for ode in rhs_eqns   ])
    
    
        if n==1:
            return FirstOrderLinearODESystem(lin_eqns,dvars=self.dvars,ivar=self.ivar)
        else:
            return FirstOrderODESystem(lin_eqns,dvars=self.dvars,ivar=self.ivar)

    def numerized(self,parameters={},ic_list=[]):
        '''
        Takes values of parameters, substitutes it into the list of parameters and changes it into a Tuple. Returns instance of class OdeComputationalCase.
        '''
        return OdeComputationalCase(odes_system=self.odes_rhs,dvars=self.dvars,ivar=self.ivar)
        

    def linearized(self, x0=None, op_point=False, hint=[], label=None):
        """
        Returns approximated first order function calculated with Taylor series method as an instance of the class. It enables to obtain linearized output.
        Arguments:
        =========
            System = Created system based on symbolical represent of mechanical parts of it
            
            op_point - boolean, which points out if the operating point will be evaluated
            
            x0 - setting operating point
            
            hint - (optional) Adds additional equation to equilibrium condition and calculate op_point as equilibrium system.
            
            label=None (optional): string
                Label of the class instance. Default label: '{Class name} with {length of qs} DOF'

        Example:
        =======
        Creating the examplary system. A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') 
        >>> Pendulum()

        Creating linerized system in symbolic pattern
        >>> System_linearized = Sytem.linearized()

        """

        return self.approximated(n=1, x0=x0, op_point=op_point, hint=hint, label=label)
    
    
    
    
    
    
    def subs(self,*args,**kwargs):
        

        
        obj = super().subs(*args,**kwargs)

        obj._lhs = self.lhs.subs(*args,**kwargs)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj
    
    
    def __add__(self,other):
        
        if isinstance(other,self.__class__):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
        obj = super().__add__(other)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj

    
    def __rsub__(self,other):
        
        if isinstance(other,self.__class__):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
        obj = super().__rsub__(other)
        obj._lhs=self._lhs
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj

    
    def __sub__(self,other):
        
        if isinstance(other,self.__class__):
            other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
        obj = super().__sub__(other)
        
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj

    
    
    def __mul__(self,other):
        
        obj = super().__mul__(other)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj
    
    def doit(self,**hints):
        
        obj = super().doit(**hints)
        obj._lhs = self._lhs.doit(**hints)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj
    
    def copy(self):
        
        
        
        return copy.copy(self)
    
    
    def _eval_applyfunc(self, f):
        
        obj = super()._eval_applyfunc(f)
        obj._lhs = self.lhs._eval_applyfunc(f)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        
        return obj
    

    def is_linear(self):
        
        if isinstance(self,self.linearized()):
            return print('True')
        else:
            return print('False')


    def numerized(self,parameters={},ic_list=[]):
        '''
        Takes values of parameters, substitutes it into the list of parameters and changes it into a Tuple. Returns instance of class OdeComputationalCase.
        '''

        
        ode=self.as_first_ode_linear_system()
        
        return OdeComputationalCase(odes_system=ode.rhs,dvars=ode.dvars,ivar=ode.ivar)
    
    

class FirstOrderODESystem(ODESystem):
    

    
    def solution(self,ics=None):
        return self.linearized().solution()


        

    def linearized(self, x0=None, op_point=False, hint=[], label=None):
        """
        Returns approximated first order function calculated with Taylor series method as an instance of the class. It enables to obtain linearized output.
        Arguments:
        =========
            System = Created system based on symbolical represent of mechanical parts of it
            
            op_point - boolean, which points out if the operating point will be evaluated
            
            x0 - setting operating point
            
            hint - (optional) Adds additional equation to equilibrium condition and calculate op_point as equilibrium system.
            
            label=None (optional): string
                Label of the class instance. Default label: '{Class name} with {length of qs} DOF'

        Example:
        =======
        Creating the examplary system. A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') 
        >>> Pendulum()

        Creating linerized system in symbolic pattern
        >>> System_linearized = Sytem.linearized()

        """

        return self.approximated(n=1, x0=x0, op_point=op_point, hint=hint, label=label)


    def subs(self,*args,**kwargs):
        

        
        obj = super().subs(*args,**kwargs)
        obj._lhs=self._lhs.subs(*args,**kwargs)
        
        return obj
    


    
    
class FirstOrderLinearODESystem(FirstOrderODESystem):
    
    @classmethod
    def from_odes(cls,odes_system,dvars,ivar=None,ode_order=1,parameters=None):
        
        
        diffs = dvars.diff(ivar,ode_order)

        diffs_coeffs_mat  = odes_system.jacobian(diffs) #it should be reimplemented with .jacobian in Analytical Solution class
        
        if diffs_coeffs_mat == Matrix.zeros(len(diffs)):
        
            odes = odes_system
        else:
            odes = diffs_coeffs_mat.inv() * (-odes_system.subs({ d_dvar:0     for d_dvar in diffs } ))
            
        f_ord_dvars=Matrix(sum([list(dvars.diff(ivar,no))   for no in  range(ode_order)],[]))
        
        display(odes)
        display(f_ord_dvars)
            
        return cls._constructor( odes=f_ord_dvars.diff(ivar) , dvars=f_ord_dvars  , odes_rhs = Matrix([f_ord_dvars[len(dvars):],odes] )  ,ivar=ivar,ode_order=ode_order,parameters=parameters)

    @classmethod
    def from_ode_system(cls,ode_system):

        
        
        sys = ode_system
        odes=sys.lhs
        odes_rhs = sys.rhs
        dvars = sys.dvars
        ivar = sys.ivar
        ode_order = 1
        parameters = sys._parameters
        


        
        
        f_ord_dvars=sys._fode_dvars
        ode_rhs = Matrix(list(f_ord_dvars[len(dvars):]) + list(sys.odes_rhs))

        
        return cls._constructor(f_ord_dvars.diff(ivar) , f_ord_dvars,ode_rhs, ivar = ivar,ode_order=ode_order ,parameters=parameters)
    
    
    @cached_property
    def odes_rhs(self):
        return self.rhs
    
    @cached_property
    def _fundamental_matrix(self):
        
        return self.odes_rhs.jacobian(self.dvars)
    
    
    
    @cached_property
    def _free_terms(self):
        
        return self.odes_rhs.subs({coord:0 for coord in self.dvars})


    @cached_property
    def _auxiliary_fundamental_matrix(self):
        
        odes_list = list(self.odes_rhs)
        dvars_list = list(self.dvars)
        
        if len(self.dvars) >= 2:
            no = int(len(self.dvars)/2)
        
            odes=  odes_list[no:] + odes_list[:no] 
            dvars = dvars_list[no:] + dvars_list[:no]
        else:
            odes=odes_list
            dvars=dvars_list
        
        odes = Matrix(odes)
        dvars = Matrix(dvars)
        
        return odes.jacobian(dvars)
    
    @cached_property
    def _auxiliary_free_terms(self):
        
        odes_list = list(self.odes_rhs)
        dvars = self.dvars
        
        if len(self.dvars) >= 2:
            no = int(len(self.dvars)/2)
        
            odes=  odes_list[no:] + odes_list[:no] 

        else:
            odes=odes_list

        odes = Matrix(odes)
        

        
        return odes.subs({coord:0 for coord in dvars})   
    
    
    
    def _const_mapper(self,const_set,parameters=None):
        
        if parameters is None:
            parameters = self.parameters
        
        const_base_dict={dummy_symbol : Symbol(f'C_{no+1}')   for  no,dummy_symbol  in enumerate(const_set)}
        
        
        
        
        if parameters is None:
            const_dict=const_base_dict
            
            #const_dict=solve([key-val  for  key,val in const_base_dict.items()] , dummies_set )
            
        else:
            const_fun_dict = {dummy_symbol : Function(str(c_symbol))   for  dummy_symbol,c_symbol  in const_base_dict.items()}
            const_dict = {dummy_symbol : c_fun(*parameters)   for  dummy_symbol,c_fun  in const_fun_dict.items()}
            
            #display([key-val  for  key,val in const_dict.items()],list(const_dict.values()))
            
            #const_dict=solve([key-val  for  key,val in const_dict.items()] , dummies_set  )
            
            

        #display(const_dict)
            
        self._const_list = list(const_dict.values())
        

        
        return const_dict

    @cached_property
    def solution_map(self):
        return lambda obj: TR8(TR10(TR8(obj.doit(conds='none')).expand()).expand())
        
    
    @cached_property
    def _general_solution(self):

        
        
        
        A = self._auxiliary_fundamental_matrix
        
        sol = AnalyticalSolution(self.dvars,linodesolve(A,t=self.ivar,b=0*self.dvars))#.applyfunc(self.solution_map)
        
        dummies_set= sol.atoms(Dummy)
        const_dict = self._const_mapper(dummies_set)

        if len(self.dvars) >= 2:
            no = int(len(self.dvars)/2)
        
            sol=  sol[no:] + sol[:no]

        return AnalyticalSolution(self.dvars,sol).subs( const_dict )
                                 

                                 
    @cached_property
    def _steady_solution(self):
                                 
        
        
        
        A = self._auxiliary_fundamental_matrix
        b = self._auxiliary_free_terms
                                 

        sol = AnalyticalSolution(self.dvars,linodesolve(A,t=self.ivar,b=b)).applyfunc(self.solution_map)
        
        dummies_set= sol.atoms(Dummy)


        if len(self.dvars) >= 2:
            no = int(len(self.dvars)/2)
        
            sol=  sol[no:] + sol[:no] 

        return AnalyticalSolution(self.dvars,sol).subs({dum_sym:0 for dum_sym in dummies_set})

    @cached_property
    def const_set(self):

        return self._const_list
    
    
    @cached_property
    def general_solution(self):

        return self._general_solution


    @cached_property
    def steady_solution(self):

        return self._steady_solution
    
    @property
    def solution(self):

        return self.general_solution + self.steady_solution
    

    def _latex(self,*args):

        
        return latex(Eq(self.lhs,self.rhs ,evaluate=False   ))    
    
    
    
class FirstOrderODE:
    """
    This class represents the system of first order defferential equations. Allows to get the solution of given system.
    """

    _const_list = set()

    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 equation_type=None):
        '''
        Supply the following arguments for the initialization of FirstOrderODE
        
        Args:
        '''
        #         if isinstance(odes_system, LinearDynamicSystem):
        #             self.dvars = odes_system.q
        #             odes_system = odes_system._eoms

        self.odes_system = odes_system

        #         print('solver init')

        #         display(self.odes_system)

        self.governing_equations = self.odes_system
        self.ivar = ivar

        if len(dvars) > 0:
            self.dvars = dvars

        self.ic_point = ic_point

        if isinstance(params, dict):
            self.params_values = params
            self.params = list(self.params_values.keys())

        elif (all(isinstance(elem, Symbol) for elem in params)):
            self.params = params
            self.params_values = params_values

        elif (all(
                isinstance(elem, tuple) and len(elem) == 2
                for elem in params)):
            self.params_values = {var: value for var, value in params}
            self.params = list(self.params_values.keys())

        else:
            self.params = self.odes_system.free_symbols
            self.params.remove(self.ivar)
            self.params_values = params_values

        self.t_span = t_span
        self.__numerical_odes = None

        self.eq_type = equation_type

    def stiffness_matrix(self):
        '''
        Returns the system stiffness matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized coordinates vector.
        '''

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        return odes.jacobian(reg_vars).subs({coord: 0
                                             for coord in reg_vars}).doit()

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        dvars_ddot = list(sym.Matrix(reg_vars).diff(self.ivar, 2))

        result = odes.jacobian(dvars_ddot).subs(
            {coord: 0
             for coord in reg_vars}).doit()

        #         print('------------------- linear mat ------------------')
        #         display(self.odes_system,self.dvars,dvars_ddot,result)

        #         print('------------------- linear mat ------------------')

        return result

    def damping_matrix(self):
        '''
        Returns the system damping matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        dvars_dot = list(sym.Matrix(reg_vars).diff(self.ivar, 1))

        result = odes.jacobian(dvars_dot).subs(
            {coord: 0
             for coord in reg_vars}).doit()

        return result

    def external_forces(self):
        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        return odes.subs({gen_coord: 0 for gen_coord in self.dvars}).doit()

    def diagonalize(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        ode_sys = self.odes_system

        #display(ode_sys)

        main_matrix = ode_sys.jacobian(self.dvars).subs(
            {dvar: 0
             for dvar in self.dvars})

        inertia_matrix = ode_sys.jacobian(self.dvars).subs(
            {dvar.diff(self.ivar): 0
             for dvar in self.dvars}) #.diff(self.ivar)

        #         display(inertia_matrix,main_matrix)

        linear_odes = inertia_matrix.inv() * (
            main_matrix * sym.Matrix(self.dvars) +
            ode_sys.subs({dvar: 0
                          for dvar in self.dvars}).doit())

        linear_odes_dict = {
            coord: eq
            for coord, eq in zip(
                self.dvars,
                linear_odes,
            )
        }

        const_odes_list = {}
        regular_odes_list = {}

        for coord, ode in (linear_odes_dict.items()):
            if ode == 0:
                const_odes_list[coord] = ode
            else:
                regular_odes_list[coord] = ode


#         print('regular odes')
#         display(*list(regular_odes_list.values()))

        if isinstance(self.ivar, Function):
            regular_vars = Matrix(list(
                regular_odes_list.values())).atoms(Function) - {self.ivar}
        else:
            regular_vars = Matrix(list(
                regular_odes_list.values())).atoms(Function)

        regular_vars = regular_odes_list.keys()

        #         display(regular_vars)

        #display( row  if row  in Matrix(regular_odes_list).rows  )
        #         display( list(regular_odes_list.values()) )
        regular_main_matrix = Matrix(list(
            regular_odes_list.values())).jacobian(list(regular_vars))
        #singular_odes=[no  for no in  const_odes_list]
        #         print('diagonalize')
        #         print([Matrix(list(regular_odes_list.values())) ,list(regular_vars) ,const_odes_list])
        #         print('diagonalize')

        return [
            Matrix(list(regular_odes_list.values())),
            list(regular_vars), const_odes_list
        ]

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        odes, vars = self.diagonalize()[0:2]

        return odes.jacobian(vars).diagonalize()[1]

    def eigenmodes(self):
        '''
        Determines the system eigenmodes matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        odes, vars = self.diagonalize()[0:2]
        #         display(odes,vars)

        return odes.jacobian(vars).diagonalize()[0]

    def damped_natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        natural_freqs = list({(im((eigen))).doit()
                              for eigen in self.eigenvalues()
                              if not eigen == 0})

        return diag(*natural_freqs)

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''

        #         print('------------------- linear gen sol ----------------')
        #         display(self.odes_system,self.dvars)

        #         display(self.inertia_matrix(),self.stiffness_matrix())
        #         print('------------------- linear gen sol ----------------')

        C = numbered_symbols('C', start=1)

        C_list = []

        for i in range(len(self.dvars) * 2):

            C_list += [next(C)]

        args_list = self.dvars[0].args

        if len(self.dvars[0].args) > 1:

            params = {*args_list} - {self.ivar}

            C_list = [Function(str(C_tmp)) for C_tmp in C_list]
            C_list = [(C_tmp)(*params) for C_tmp in C_list]

        #         print('o tu')
        self.__class__._const_list |= set(C_list)

        modes, eigs = self.eigenmodes(), self.eigenvalues()
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        #         display( (const_part))

        #Y_mat = Matrix(self.dvars)

        #diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        # t_sol = self.ivar

#         display(modes,eigs)
        solution = [
            C_list[i] * modes[:, i] * exp(eigv * self.ivar)
            for i, eigv in enumerate([eigv for eigv in eigs if not eigv == 0])
        ]
        #         print('sol')
        #         display(solution)
        #         print('sol')

        const_dict = ({
            coord: C_list[no + len(solution)]
            for no, coord in enumerate(const_part.keys())
        })

        sol_dict = {
            coord: eq
            for coord, eq in zip(
                reg_vars, list(sum(solution, Matrix([0] * len(solution[0])))))
        }

        return {**sol_dict, **const_dict}

    def steady_solution(self, initial_conditions=None):

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        ext_forces = self.external_forces().expand().applyfunc(
            lambda row: (TR8(row).expand()))

        #         sin_components=ext_forces.atoms(sin)
        #         cos_components=ext_forces.atoms(cos)

        #        display('sin cos reco',)
        components = ext_forces.atoms(sin, cos)

#         display('ext_forces', ext_forces)
#         display(odes)
#         display(reg_vars)

        steady_sol = Matrix([0 for gen_coord in self.dvars])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
#             print('omg')
#             display(omg)
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

#             display(self.inertia_matrix())
#             display(self.damping_matrix())

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (
                sym.re(fund_mat.inv()).doit() * amp_vector) * comp + (
                    sym.im(fund_mat.inv()).doit() * amp_vector) * TR3(
                        comp.subs(omg * self.ivar, omg * self.ivar - pi / 2))

            #steady_sol += ((fund_mat.inv()) * amp_vector) * comp

            ext_forces -= (amp_vector * comp).expand()
        #print(ext_forces.doit().expand())

        # const_elems = lambda expr: [
        #     expr for expr in comp.expand().args if not expr.has(self.ivar)
        #     if isinstance(comp, Add)
        # ]

        const_mat = Matrix([
            sum((expr for expr in comp.expand().args if not expr.has(self.ivar)
                 if isinstance(comp, Add)), 0) for comp in ext_forces.doit()
        ])

        #        display('const',const_mat)

        steady_sol += (self.stiffness_matrix().inv() * const_mat).doit()

        ext_forces -= const_mat

        #        display('res',ext_forces)

        if ext_forces.doit().expand() != sym.Matrix(
            [0 for gen_coord in self.dvars]):
            #             print('o tu')
            #             display((self.governing_equations - self.external_forces() +
            #                      ext_forces).expand().doit())
            eqns_res = (self.governing_equations - self.external_forces() +
                        ext_forces).expand().doit()

            if len(self.dvars) == 1:
                steady_sol += Matrix([
                    sym.dsolve(eqns_res[0], self.dvars[0]).rhs.subs({
                        Symbol('C1'):
                        0,
                        Symbol('C2'):
                        0
                    })
                ])
            else:
                steady_sol += sym.dsolve(eqns_res, self.dvars)

        return {coord: -eq for coord, eq in zip(reg_vars, steady_sol)}

    def solution(self, initial_conditions=None):
        return {
            key:
            self.general_solution(initial_conditions=initial_conditions)[key] +
            self.steady_solution(initial_conditions=initial_conditions)[key]
            for key in self.general_solution(
                initial_conditions=initial_conditions).keys()
        }


class LinearODESolution:

    _const_list = set()
    _cache={}
    

    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 equation_type=None):
        '''
        Supply the following arguments for the initialization of OdeComputationalCase:
        
        Args:
        '''
        #         if isinstance(odes_system, LinearDynamicSystem):
        #             self.dvars = odes_system.q
        #             odes_system = odes_system._eoms

        self.odes_system = odes_system

        #         print('solver init')

        #         display(self.odes_system)

        self.governing_equations = self.odes_system
        self.ivar = ivar

        if len(dvars) > 0:
            self.dvars = dvars

        self.ic_point = ic_point

        if isinstance(params, dict):
            self.params_values = params
            self.params = list(self.params_values.keys())

        elif (all(isinstance(elem, Symbol) for elem in params)):
            self.params = params
            self.params_values = params_values

        elif (all(
                isinstance(elem, tuple) and len(elem) == 2
                for elem in params)):
            self.params_values = {var: value for var, value in params}
            self.params = list(self.params_values.keys())

        else:
            self.params = self.odes_system.free_symbols
            self.params.remove(self.ivar)
            self.params_values = params_values

        self.t_span = t_span
        self.__numerical_odes = None

        self.eq_type = equation_type

    def stiffness_matrix(self):
        '''
        Returns the system stiffness matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized coordinates vector.
        '''
        return self.governing_equations.jacobian(self.dvars).subs(
            {coord: 0
             for coord in self.dvars}).doit()

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_ddot = list(sym.Matrix(self.dvars).diff(self.ivar, 2))

        result = self.governing_equations.jacobian(dvars_ddot).subs(
            {coord: 0
             for coord in self.dvars}).doit()

        #         print('------------------- linear mat ------------------')
        #         display(self.odes_system,self.dvars,dvars_ddot,result)

        #         print('------------------- linear mat ------------------')

        return result
    
    def fundamental_matrix(self, freq=Symbol('omega', positive=True)):
        '''
        Method returns a fundamental matrix of the system built from inertia and stiffness matrices. Takes one optional argument.

        Args:
            freq (optional, obj:Symbol): This argument lets user to choose a symbol for frequency representation. Set to 'omega' by default, 'positive=True' ensures its affiliation to real numbers domain.
        '''
        return -freq**2 * self.inertia_matrix() + self.stiffness_matrix()

    def damping_matrix(self):
        '''
        Returns the system damping matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_dot = list(sym.Matrix(self.dvars).diff(self.ivar, 1))

        return self.governing_equations.jacobian(dvars_dot).subs(
            {coord: 0
             for coord in self.dvars}).doit()

    def external_forces(self):
        return self.odes_system.subs(
            {gen_coord: 0
             for gen_coord in self.dvars}).doit()

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        #         display(self.governing_equations)

        q_dot = (Matrix(self.dvars).diff(self.ivar))

        ode_sys = Matrix([
            q_dot,
            self.inertia_matrix().inv() *
            (-self.stiffness_matrix() * Matrix(self.dvars) -
             self.damping_matrix() * q_dot)
        ])

        #         display(ode_sys)

        main_matrix = ode_sys.jacobian(list(self.dvars) + list(q_dot))

        #         display(main_matrix)

        return (main_matrix).diagonalize()[1]

    def damped_natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        natural_freqs = list({(im((eigen))).doit()
                              for eigen in self.eigenvalues()
                              if not eigen == 0})

        return diag(*natural_freqs)
    
    def natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        
        natural_freqs=list({(sqrt(abs(eigen**2))) for eigen  in self.eigenvalues() if not eigen==0 })
        
        return diag(*natural_freqs)

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''

        #         print('------------------- linear gen sol ----------------')
        #         display(self.odes_system,self.dvars)

        #         display(self.inertia_matrix(),self.stiffness_matrix())
        #         print('------------------- linear gen sol ----------------')

        C = numbered_symbols('C', start=1)

        C_list = []

        for i in range(len(self.dvars) * 2):

            C_list += [next(C)]

        args_list = self.dvars[0].args

        if len(self.dvars[0].args) > 1:

            params = {*args_list} - {self.ivar}

            C_list = [Function(str(C_tmp)) for C_tmp in C_list]
            C_list = [(C_tmp)(*params) for C_tmp in C_list]

        #         print('o tu')
        #         display(self.odes_system)

        self.__class__._const_list |= set(C_list)

        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        eigs = self.eigenvalues()

        Y_mat = Matrix(self.dvars)

        #diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        t_sol = self.ivar

        solution = [(C_list[2 * i] * modes[:, i] *
                     sin(im(eigs[2 * i + 1, 2 * i + 1]).doit() * t_sol) +
                     C_list[2 * i + 1] * modes[:, i] *
                     cos(im(eigs[2 * i + 1, 2 * i + 1]).doit() * t_sol)) *
                    exp(re(eigs[i, i]).doit() * t_sol)
                    for i, coord in enumerate(self.dvars)]

        return sum(solution, Matrix([0] * len(Y_mat)))

    def steady_solution(self, initial_conditions=None):
        if (tuple(self.dvars),ImmutableMatrix(self.governing_equations)) in self.__class__._cache:


            #print('cache for Linear ODE')
            steady_sol=self.__class__._cache[(tuple(self.dvars),ImmutableMatrix(self.governing_equations))]
            
        else:
            #print('new solution')
            ext_forces = self.external_forces().expand().applyfunc(
                lambda row: (TR8(row).expand()))

            #         sin_components=ext_forces.atoms(sin)
            #         cos_components=ext_forces.atoms(cos)

            #        display('sin cos reco',)
            components = [
                comp for comp in ext_forces.atoms(sin, cos) if comp.has(self.ivar)
            ]

            #        display('ext_forces',ext_forces)

            steady_sol = Matrix([0 for gen_coord in self.dvars])

            for comp in components:
                #display(comp)

                omg = (comp.args[0].diff(self.ivar)).doit()
                #            display(omg)
                amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

                #display(amp_vector)

                fund_mat = -self.inertia_matrix(
                ) * omg**2 + sym.I * omg * self.damping_matrix(
                ) + self.stiffness_matrix()

                #print('++linea.py++'*10)
                #display(fund_mat)
                #print('++linea.py++'*10)

                steady_sol += (
                    sym.re(fund_mat.inv()).doit() * amp_vector) * comp + (
                        sym.im(fund_mat.inv()).doit() * amp_vector) * TR3(
                            comp.subs(omg * self.ivar, omg * self.ivar - pi / 2))

                #steady_sol += ((fund_mat.inv()) * amp_vector) * comp

                ext_forces -= (amp_vector * comp).expand()
            #print(ext_forces.doit().expand())

            # const_elems = lambda expr: [
            #     expr for expr in comp.expand().args if not expr.has(self.ivar)
            #     if isinstance(comp, Add)
            # ]

            const_mat = Matrix([
                sum((expr for expr in comp.expand().args if not expr.has(self.ivar)
                     if isinstance(comp, Add)), 0) for comp in ext_forces.doit()
            ])

            #        display('const',const_mat)

            steady_sol += (self.stiffness_matrix().inv() * const_mat).doit()

            ext_forces -= const_mat

            #        display('res',ext_forces)

            if ext_forces.doit().expand() != sym.Matrix(
                [0 for gen_coord in self.dvars]):
                #             print('o tu')
                #             display((self.governing_equations - self.external_forces() +
                #                      ext_forces).expand().doit())
                eqns_res = (self.governing_equations - self.external_forces() +
                            ext_forces).expand().doit()

                if len(self.dvars) == 1:
                    steady_sol += Matrix([
                        sym.dsolve(eqns_res[0].expand().doit(), self.dvars[0]).rhs.subs({
                            Symbol('C1'):
                            0,
                            Symbol('C2'):
                            0
                        })
                    ])
                else:
                    steady_sol += sym.dsolve(eqns_res, self.dvars)
                    
            self.__class__._cache[(tuple(self.dvars),ImmutableMatrix(self.governing_equations))] = steady_sol

                
        return -steady_sol

    def solution(self, initial_conditions=None):
        return self.general_solution(
            initial_conditions=initial_conditions) + self.steady_solution(
                initial_conditions=initial_conditions)
