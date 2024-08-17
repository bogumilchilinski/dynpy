from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, re, pi, latex,
                   dsolve, solve, fraction, factorial, Add, Mul, exp, zeros, shape,
                   numbered_symbols, integrate, ImmutableMatrix,Expr,Dict,Subs,Derivative,Dummy,
                   lambdify, Pow, Integral, init_printing, I, N,eye, zeros, det, Integer,separatevars,Heaviside,simplify)

from sympy.matrices.matrices import MatrixBase
from sympy.solvers.ode.systems import matrix_exp, matrix_exp_jordan_form
from sympy.solvers.deutils import ode_order
from sympy import classify_ode

from numbers import Number

###  exemplary comment
from sympy.physics.mechanics import dynamicsymbols, init_vprinting
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from ..utilities.adaptable import TimeSeries, TimeDataFrame, NumericalAnalysisDataFrame

from collections import ChainMap

from IPython.display import display, Markdown

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7,TR6, TR5, TR3
from collections.abc import Iterable
from sympy.solvers.ode.systems import linodesolve
from functools import cached_property, lru_cache
cached_property = property

from .numerical import OdeComputationalCase

import time
import pandas as pd

from ..utilities.report import (SystemDynamicsAnalyzer,DMath,ReportText,SympyFormula, AutoBreak, PyVerbatim)
from ..utilities.templates.document import *
from ..utilities.templates import tikz
from ..utilities.components.ode import en as ode
from ..utilities.components.ode import pl as ode_comp_pl



import copy

from .tools import CommonFactorDetector, ODE_COMPONENTS_LIST, CodeFlowLogger

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




class AnalyticalSolution(ImmutableMatrix):

    _default_doctype = ExampleTemplate

    def __new__(cls, elements, vars = None, rhs=None, evaluate=True, **options):

        
        
        if isinstance(elements,(dict,Dict)):
            return cls.from_dict(elements,vars=vars,**options)
        if isinstance(elements,Eq):
            return cls.from_eq(elements,vars=vars,**options)
        if isinstance(elements,Iterable) and all([ isinstance(elem,Eq) for elem in  elements ]):
            return cls.from_eqns_list(elements,vars=vars,**options)
        else:
            return cls._constructor(elements,vars, rhs, evaluate=evaluate, **options)



    
    @classmethod
    def _constructor(cls, elements, vars= None ,rhs=None, evaluate=True, ics=None, **options):
        
        if not isinstance(elements,Iterable): 
            elems_lhs = Matrix([elements])
        elif isinstance(elements,Iterable):
            elems_lhs = Matrix(elements)
            
            
        # it should be reimplemented with @property if None as odes_rhs is spotted            
        if isinstance(rhs,Iterable):
            # display(rhs)
            # display(type(rhs))
            elems_rhs = Matrix(rhs)
            
        elif rhs is None:
            elems_rhs = None
        elif not isinstance(rhs,Iterable):
            elems_rhs = Matrix([rhs])            
            

        if elems_rhs is not None:
            elems = elems_lhs - elems_rhs
        else:
            elems = elems_lhs

        obj = super().__new__(cls, elems ,evaluate=evaluate, **options)

        obj._rhs= rhs    
        obj._vars = vars
        obj._ics = ics
        
        if obj.vars is not None:
            if len(obj.vars)!=len(obj):
            
                warn_str = '!!!!!!!!!!!!!!!!! WRONG SIZE of VARS and ELEMENTS'
                print(warn_str)
                return warn_str

        return obj
    
    @classmethod
    def from_dict(cls,dictionary,vars=None,**options): # działa
        '''
        This constructor creates object from map in the following form:
        
        vars -> right hand sides
        '''
        vars=Matrix(list(dictionary.keys()))
        elems=Matrix(list(dictionary.values()))
        
        return cls._constructor( vars ,vars = vars, rhs = elems ,**options   )

    @classmethod
    def from_vars_and_rhs(cls,vars,rhs, ics=None,**options): # działa
        '''
        This constructor creates object from map in the following form:
        
        vars -> right hand sides
        '''
            
        return cls._constructor( vars ,vars = vars, rhs = rhs, ics=ics ,**options   )


    
    @classmethod
    def from_eq(cls,matrix_eq,vars=None,**options): #
        if isinstance(matrix_eq,Eq):
            
            eq = matrix_eq
            
        else:
            print ('TypeError: matrix_eq is not Eq')
            return None
        
        obj = cls._constructor(eq.lhs,vars=vars,rhs=eq.rhs ,evaluate=True ,**options   )
        
        return obj

    @classmethod
    def from_eqns_list(cls,eqns_list,vars=None,**options):

        elems=[eqn.lhs for eqn in eqns_list]
        values=[eqn.rhs for eqn in eqns_list]
            
        obj = cls._constructor(Matrix(elems),vars=vars,rhs = Matrix(values),evaluate=True ,**options)
    
        return obj


    @classmethod
    def from_eqns_matrix(cls,eqns_matrix,vars=None,**options):

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
            
        obj = cls._constructor(Matrix(dvars),vars=vars,rhs=Matrix(values),evaluate=True ,**options)
    
        return obj
    
    def copy(self):
        
        obj=self._constructor(list(self),self._vars, self._rhs)
        obj = self._assign_properties(obj)
        
        
        return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    
    
    def _assign_properties(self,obj):
        
        #obj._rhs=self._rhs     
        obj._vars = self._vars
        #obj._const_list = self._const_list
        
        return obj
    
    
    def subs(self,*args,**kwargs):
        

        #obj = super().subs(*args,**kwargs)
        #obj = self._assign_properties(obj)
        result = self.rhs.subs(*args,**kwargs)
        
        
        #obj = self._assign_properties(obj)
        
        return type(self).from_vars_and_rhs(self.lhs,result)

    @cached_property 
    def _is_rhs_form(self):
        
        return self.lhs == self.vars


    def _fix_other(self,other):
        if isinstance(other,AnalyticalSolution):
            dict_elems = other.as_explicit_dict()
            
            display(dict_elems)
            display(self.vars)
            other = Matrix([other.as_explicit_dict()[coord]  for  coord  in self.vars ])
            #other = other.rhs        
        
        
    def _add_rhs_form(self,other):
        if isinstance(other,AnalyticalSolution):
            dict_elems = other.as_explicit_dict()
        else:
            dict_elems = AnalyticalSolution.from_vars_and_rhs(self.vars,other).as_explicit_dict()
            
        result = Matrix([self.as_explicit_dict()[coord]+ dict_elems[coord] for  coord  in self.lhs ])
        return type(self).from_vars_and_rhs(self.lhs,result)


    def _op_TYPE_rhs_form(self,other,op_func):
        if isinstance(other,AnalyticalSolution):
            dict_elems = other.as_explicit_dict()
        else:
            dict_elems = AnalyticalSolution.from_vars_and_rhs(self.vars,other).as_explicit_dict()
            
        
        elems_list = [op_func(self.as_explicit_dict()[coord],dict_elems[coord]) for  coord  in self.lhs ]
        
        result = Matrix(elems_list)
        obj=type(self).from_vars_and_rhs(self.lhs,result)
        self._assign_properties(obj)
        
        return self._assign_properties(obj)
    
    
    def _op_TYPE_sides_form(self,other,op_func):

        CodeFlowLogger(other,'other',self)

        if isinstance(other,AnalyticalSolution):
            other_rhs = other.rhs
            other_lhs = other.lhs
        else:
            other_rhs = other
            other_lhs = other*0
        
        
        #rhs = self.rhs.__add__(other_rhs)
        rhs = op_func(self.rhs,other_rhs)
        #lhs = self.lhs.__add__(other_lhs)
        lhs = op_func(self.lhs,other_lhs)
        
        CodeFlowLogger(rhs,'op rhs',self)
        CodeFlowLogger(lhs,'op lhs',self)
        
        obj = type(self)._constructor(lhs,self.vars,rhs)

        CodeFlowLogger(obj,'result',self)

        return self._assign_properties(obj)
    
    
    
    def __add__(self,other):
        
        op_func = lambda obj,other: obj.__add__(other)
        
        CodeFlowLogger(self._is_rhs_form,'__add__',self)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)


    def __radd__(self,other):

        op_func = lambda obj,other: obj.__radd__(other)
        
        CodeFlowLogger(self._is_rhs_form,'__radd__',self)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)


    
    def __rsub__(self,other):

        op_func = lambda obj,other: obj.__rsub__(other)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)


    def __sub__(self,other):

        op_func = lambda obj,other: obj.__sub__(other)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)


    def __mul__(self,other):

        op_func = lambda obj,other: obj.__mul__(other)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)



    def __rmul__(self,other):

        op_func = lambda obj,other: obj.__rmul__(other)
        
        if self._is_rhs_form:
            
            return self._op_TYPE_rhs_form(other,op_func)

        else:
            
            return self._op_TYPE_sides_form(other,op_func)


    
    
    def doit(self,**hints):

        obj = super().doit(**hints)
        obj = self._assign_properties(obj)
        obj._rhs = self.rhs.doit(**hints)

        return obj
    
    def _eval_applyfunc(self, f):
        
        obj = super()._eval_applyfunc(f)
        obj = self._assign_properties(obj)
        obj._rhs = self.rhs._eval_applyfunc(f)
        
        
        return obj

    def expand(self,deep=True, modulus=None, power_base=True, power_exp=True,mul=True, log=True, multinomial=True, basic=True, **hints):
        

        obj = super().expand(deep=deep, modulus=modulus, power_base=power_base, power_exp=power_exp,mul=mul, log=log, multinomial=multinomial, basic=basic,**hints)
        obj = self._assign_properties(obj)        
        obj._rhs = obj.rhs.expand(deep=deep, modulus=modulus, power_base=power_base, power_exp=power_exp,mul=mul, log=log, multinomial=multinomial, basic=basic,**hints)
        #obj._dvars=self._dvars
        #obj._ivar = self.ivar
        
        return obj
    
    
    def __call__(self,t,params={}):
        
        solution = (self.nth_order_solution(
            ).rhs.subs(self.extra_params).subs(self.params_values))
        
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
    def vars(self):
        if self._vars is not None:
            if isinstance(self._vars,Iterable):
                return self._vars
            else:
                return Matrix([self._vars])
                            
        else:
            return [Symbol(f'Eq_{no}')  for no in range( len(self.as_list())  )]
    
    
    @property
    def lhs(self):
        if self._rhs is None:

            return Matrix(list(self))
        else:
            return Matrix(list(self)) + self.rhs

    @property
    def rhs(self):
        if self._rhs is None:
            return Matrix(list(self))*0
        
        elif isinstance(self._rhs,Iterable):
            return Matrix(list(self._rhs)) 
        
        else:
            return Matrix([self._rhs])

    def as_list(self):
        
        return list(self)

    def as_matrix(self):
        '''
        Returns the right-hand side of the equation as a matrix.
        It returns Sympy Matrix object
        '''
        return Matrix(self.as_list())
    
    def as_eq(self):
        ''' Creates an equation using the Eq method from the Sympy library.
        It returns Sympy Equation class.
        '''
        return Eq(self.lhs,self.rhs)
    
    def as_iterable(self):

        return [(lhs,comp) for lhs,comp  in zip(self.vars,self.as_list())]
        
    
    def as_dict(self):

        return Dict({lhs:rhs  for  lhs,rhs in self.as_iterable()})
    
    def as_explicit_dict(self):

        return Dict({lhs:rhs  for  lhs,rhs in zip(self.lhs,self.rhs)})
    
    def as_eq_list(self):
        ''' Creates a zip object consisting of the left side of the self instance and self itself.
        It returns zip object.'''
        return [ Eq(lhs,comp,evaluate=False) for lhs,comp  in zip(self.lhs,self.rhs)]
    def as_eq_diff_list(self,ivar,diff_order):
        ''' Creates a zip object consisting of the left side of the self instance and self itself.
        It returns zip object.'''
        return [ Eq(lhs.diff(ivar,diff_order),comp.diff(ivar,diff_order),evaluate=False) for lhs,comp  in zip(self.lhs,self.rhs)]
    
    def system_parameters(self, parameter_values=None):
        '''
        Recognises system parameters as symbols which are not independent variable or its relations and returns the Tuple (Sympy object) containing these elements. It does not take into account undefined functions (instances of Function class) of independent variable.
        '''
        expr = (self.lhs-self.rhs).free_symbols
        ivar = self.ivar
        return SystemParameter(expr, ivar, parameter_values).system_parameters
#         params = (self.lhs-self.rhs).free_symbols
#         params.remove(self.ivar)
#         if parameter_values == None:
#             return list(params)
#         else:
#             return {
#                 param: parameter_values[no]
#                 for no, param in enumerate(params)
#             }
    
    @property
    def _lhs_repr(self):
        return self.lhs
    
    def __repr__(self):

        return f'{self._lhs_repr} = {self.rhs}'

    # @property
    # def _dvars(self):
    #     return self.lhs

    # @_dvars.setter
    # def _dvars(self,dvars):
    #     self._lhs = dvars

    @cached_property
    def dvars(self):
        return self._vars
    
    @property
    def _fode_dvars(self):
        return self.dvars
    
    def _latex(self,*args):
        # print('abc')
        # display(self.lhs)
        # display(self.rhs)
        return latex(Eq(self.lhs,self.rhs,evaluate=False)) + f'~~for~~{latex(self.vars)}'


    def numerized(self,parameters={},t_span=[],**kwargs):
        
        #ivar = Symbol('t')
        
        solution = self.subs(parameters).doit()
        

        
        return solution
    
    def compute_solution(self,
                         t_span=None,
                         ic_list=None,
                         t_eval=None,
                         params_values=None,
                         method='RK45',
                         derivatives=False):
        '''
        Returns the result of the computations of solve_ivp integrator from scipy.integrate module.
        '''


#         if not self._evaluated:
#             self.form_numerical_rhs()
#             self._evaluated=True
        
        t_0=time.time()

        solution = self.rhs

        ivar = list(self.dvars[0].args)[0]
        solution_fixed=solution+Matrix([exp(-ivar)*exp(-(1e6+t_span[0])) if isinstance(expr,Number) else 0  for expr in solution])



        #             print('num'*3)
        #             display(solution)

        sol_func = lambdify(ivar, solution_fixed, 'numpy')

        numerized_data = (sol_func(t_span))

        dvars = self.lhs

        numerized_sol  = TimeDataFrame(data={dvar:data[0] for dvar,data in zip(dvars,numerized_data)},index=t_span)
        numerized_sol.index.name = ivar 

        solution_tdf = numerized_sol



        velocities = self.dvars[int(len(self.dvars)/2) :]
        velocities=self.dvars
        for vel in velocities:
            solution_tdf[vel].to_numpy()
            gradient = np.gradient(solution_tdf[vel].to_numpy(),t_span)
            solution_tdf[vel.diff(ivar)] = gradient
            
        t_e=time.time()
        
        t_d = t_e - t_0
        print('_'*100,t_d)
        
        comp_time=t_d
        #if derivatives:

        solution_tdf._set_comp_time(comp_time)
        solution_tdf.index.name = ivar

        return solution_tdf


    def _as_na_df(self,parameter=None,param_span=None,dependencies_dict=None,coordinates=None, t_span=None):

        parameters = self.system_parameters()
        
        
        if len(parameters)<1:
            parameters = [Symbol('a')]

        if parameter is None:
            params_included = 1
            params_list = parameters[0:params_included]
        elif isinstance(parameter,Number):
            params_included = parameter     
            params_list = parameters[0:params_included]  
        else:
            params_included = 1
            params_list = parameters[0:params_included]

        param_list_el = params_list[0]
            
        if param_span is None:
            #param_span = [0.8,1,1.2]
            param_span = [0.9,1.0]

        if dependencies_dict is None:
            dependencies_dict = {}
            
        if t_span is None:
            t_span = [0.0]

        reference_data = {ref_val: 1 for ref_val in  self.system_parameters()[params_included:] }
#         reference_data = {}
#         display(reference_data)

        system = self

        if coordinates is None:
            Y = list(system._fode_dvars) + list(dependencies_dict.keys())
        elif coordinates=='acc' or coordinates==2:
            Y = list(system._fode_dvars) + list(system._fode_dvars.diff(self.ivar)) + list(dependencies_dict.keys())    
        else:
            Y = list(coordinates) + list(dependencies_dict.keys())
            
        Y_fixed = list({coord:coord   for coord  in  Y}.keys()    )    
        
        index = pd.Index(t_span,name=self.ivar)
        
        df_num = NumericalAnalysisDataFrame.from_model(system,
                                                        #parameter=params_list,
                                                        parameter=param_list_el,
                                                        span=param_span,
                                                        reference_data=reference_data,
                                                        coordinates=Y_fixed,
                                                        index=index)
        


        results_num = df_num#.perform_simulations(model_level_name=0,dependencies=dependencies_dict)
        #results = TimeDataFrame(results_num).droplevel(0,axis=1)
        results= results_num

        return results
    
    def numerical_analysis(self,parameter=None,param_span=None,dependencies_dict=None,coordinates=None, t_span=None):
                       
        return self._as_na_df(parameter=parameter, param_span=param_span, dependencies_dict=dependencies_dict,coordinates=coordinates, t_span=t_span)

    @property
    def _report_components(self):

        comp_list=[
        ode.ODESystemComponent,

        ]

        return comp_list

    @property
    def default_doctype(self):

        doctype=self._default_doctype

        return doctype

    @default_doctype.setter
    def default_doctype(self,value):

        self._default_doctype=value


    @property
    def report(self):

        sys=self
        doc = self.default_doctype()


        for comp in self._report_components:
            doc.append(comp(sys))

        return doc

    
    def __str__(self):
        return "Analytical Solution of ODE"
    def __repr__(self):
        return self.__str__()
    
    def _format_str(self, printer=None):

        return self.__str__()

    

class ODESolution(AnalyticalSolution):

    _ics = None    
    _default_ics = None
    _integration_consts = None #container for integration constants
    _ivar=Symbol('t')
    _dvars_str = None
    _ivar0 = 0 
    _sol0 = 0


    def _assign_properties(self,obj):
        
        #obj._lhs=self._lhs     
        obj._ivar = self.ivar
        #obj._rhs = self._rhs
        #obj._lhs = self._lhs
        obj._vars = self._vars
                
        return obj

    def expand(self,deep=True, modulus=None, power_base=True, power_exp=True,mul=True, log=True, multinomial=True, basic=True, **hints):
        

        obj = super().expand(deep=deep, modulus=modulus, power_base=power_base, power_exp=power_exp,mul=mul, log=log, multinomial=multinomial, basic=basic,**hints)
        

        obj._ivar = self.ivar
        
        return obj
    
    @property
    def ivar(self):
        return self._ivar
    
    @ivar.setter
    def ivar(self,ivar):
        if ivar is not None:
            self._ivar = ivar
        
        return self._ivar

    def append_integration_consts(self,const_list):
        
        if self._integration_consts is None:
            self._integration_consts = const_list
        else:
            self._integration_consts += const_list
            
        return self

    def _spot_constant(self):
        """
        Auxiliary method to spot or get integration constants from the system

        Returns
        -------
        list
            Containter with integration constants
        """
        
        const_list = self._integration_consts
        if const_list is None:
            const_list =  [expr for expr  in self.rhs.atoms(Symbol) if 'C' in str(expr) ]
        
        return const_list


    @property
    def ics_dvars(self):



        return self.dvars

    def default_ics(self,critical_point=False):
        
        if isinstance(self._default_ics,dict):
            ics_instance={coord:self._default_ics[coord] for coord in self.dvars if coord in self._default_ics}
            
            return {**{coord:0 for coord in self.dvars},**ics_instance}
        else:
            return {coord:0 for coord in self.dvars}
    




    def set_ics(self,ics):
        
        obj = copy.copy(self)
        
        if isinstance(ics,dict):
            obj._ics = ics
            
        return obj


    def _ics_dynamic_symbols(self):
        symbols_list = ['v', 'a']
        ics_dynamic_symbols = [Symbol(f'{coor}_{self._dvars_str}0') for coor in symbols_list]
        ics_dynamic_symbols = Matrix([Symbol(f'{self._dvars_str}0')] + ics_dynamic_symbols)
        return ics_dynamic_symbols



    @property
    def _ics_dict(self):
        """
        Property which manages the values of solutions intial conditions

        Returns
        -------
        dict
            descrption of returned data - its structure
        """
        
        
        #self._ics = {coord:1 for coord in self.dvar}
        
        if self._ics is None:
            return self._ics
        else:
            return self.default_ics

    @property
    def ivar_0(self):
        return self._ivar0

    @ivar_0.setter
    def ivar_0(self,ivar0):

        if isinstance(ivar0,(Symbol,Number)):
        
            self._ivar0 = ivar0
        elif not ivar0.has(self.ivar):
            self._ivar0 = ivar0
        
        else:
            print(f'something went wrong - ivar0 = {ivar0} which is not proper type ')

    def _get_constant_eqns(self,ics=None, sol0=None):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        if isinstance(sol0, ODESolution):
            ics_list = sol0.subs(self.ivar, self.ivar_0)  
        elif ics is None:
            temp = self._ics_dynamic_symbols()
            ics_list = [temp[index] for index in range(len(self.dvars))]
        elif isinstance(ics,dict):
            ics_list = [ics[coord]  for coord  in self.dvars]
        elif isinstance(ics,(list,tuple)):
            ics_list = ics
        #ics_list = [ics[coord]  for coord  in self.dvars]

        # display(self.rhs.subs(self.ivar,self.ivar_0))
        # display(Matrix(ics_list))
        
        
        const_eqns=Matrix(ics_list) - self.rhs.subs(self.ivar,self.ivar_0)

        CodeFlowLogger(const_eqns,'const eqns',self)
        
        return const_eqns
        
    def _calculate_constant(self,ics=None, sol0=None):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """

        const_eqns=self._get_constant_eqns(ics, sol0)

        CodeFlowLogger(const_eqns,'const eqns',self)
        
        const_list = self._spot_constant()

        CodeFlowLogger(const_list,'const list',self)
        
        return solve(const_eqns,const_list)
    
    def _get_dvars_str(self):
        return str(type(self.dvars[0]))
    
    def get_ics(self):
        return self._ics

    def with_ics(self, ics=None, ivar0=0, sol0=None):
        self._ics = ics
        self._dvars_str = self._get_dvars_str()
        self.ivar_0 = ivar0
        const_dict=self._calculate_constant(ics, sol0)
        CodeFlowLogger(const_dict,'const dict',self)
        
        return ODESolution.from_vars_and_rhs(self.vars,self.rhs.subs(const_dict), ics=self._ics) #temporary workaround that replaces subs
    
    def __call__(self,ivar=None,ics=None,params={}):
        
        
        #self.ivar = ivar
        
        
        if ivar is None or isinstance(ivar, Symbol):

            
            return self.with_ics(ics=ics)
            
        else:
            
            const_dict=self._calculate_constant(ics)
            
            return self.subs(const_dict)

    
class ODESystem(AnalyticalSolution):

    """
    This class creates an object as ordinary differential equation (ODE) in general form and provides methods for basic operations.

    Example
    =======
    
        >>>from dynpy.solvers.linear import *
        >>>import sympy

        >>>t = Symbol('t')
        >>>Omega, omega = symbols('Omega omega',positive)
        >>>x = Function('x')(t)

        >>>ode = ODESystem(odes=Matrix([x.diff(t, t) + omega**2 * x - sin(Omega*t)]), dvars=Matrix([x]), ode_order=2)
        >>>display(ode)

        >>>fode = ode.as_first_ode_linear_system()
        >>>display(fode)

        >>>sym_sol = ode.solution
        >>>display(sym_sol)

        >>>num_sol = sym_sol.subs({Symbol('C_1'): 2, Symbol('C_2'): 2, Symbol('omega'): 2, Symbol('fi'): 2}).numerized().compute_solution(np.linspace(0, 5, 100))
        >>>display(num_sol)
        >>>display(num_sol.plot())
    """
    
    _ivar = Symbol('t')
    _simp_dict = None
    _callback_dict = None

    _simp_dict = {}
    _callback_dict = {}
    
    _default_ics = None
    _const_list = []
    
    _default_detector = CommonFactorDetector #None
    #_default_detector = None
    
    _parameters = None
    _ode_order = 2
    
    _profcheck = False

    @classmethod
    def from_ode_system(cls,ode_system):

        sys = ode_system
        odes=sys.lhs
        odes_rhs = sys.rhs
        dvars = sys.dvars
        ivar = sys.ivar
        parameters = sys._parameters
        ode_order = sys.ode_order
        
        new_sys = cls._constructor(odes , dvars, odes_rhs , ivar,ode_order=ode_order ,parameters=parameters)
        new_sys._const_list = sys._const_list
        
        #display('subs dict',sys._simp_dict,cls)
        
        return new_sys.set_simp_deps(sys._simp_dict,sys._callback_dict,inplace=True)

    @classmethod
    def from_dynamic_system(cls,dyn_system, ode_order = None, parameters = None):
        """Class method that creates an ODESystem object from a dynamic system object. It also checks parameters and takes default class variables if it is needed.

        Arguments:
            dyn_system - dynamic system object;
            ode_order - order of a differential equation. Default value is None;
            parameters - parameters that could be implemented into newly created ODESystem object. Default value is None;
        """

        sys = dyn_system
        ds_lhs = sys._eoms

        dvars = sys.q
        ivar = sys.ivar
        
        if parameters == None:
            parameters = cls._parameters

        if ode_order == None:
            ode_order = cls._ode_order

        return cls._constructor(ds_lhs, dvars, ivar=ivar, ode_order = ode_order, parameters = parameters)

    
    @classmethod
    def set_default_parameters(cls,parameters):
        cls._parameters = parameters
        return cls


    
    @classmethod
    def from_random_data(cls):

        system=cls.from_default_data()
        subs_dict=system.get_random_parameters()
        
        return system.subs(subs_dict)
    
    def __new__(cls, odes, dvars, odes_rhs=None , ivar=None, ode_order=None, evaluate=True, parameters = None, **options):
        """
        Arguments
        =========
        
            odes : Any type
                An expression representing left-hand side of ODE. Prefered type is Martix.
                
            dvars : The same type as odes
                Represents the dependent variable.

        Default Arguments
        =================
        
            odes_rhs : The same type as odes
                An expression representing right-hand side of ODE, if odes_rhs is provided ODESystem converts equation into general form.
                
            ivar: Symbol
                Represents the independent variable and it's sympy.symbol('t') by default.
             
            ode_order : int 
                Value that represents order of ODE and it's 1 by default. If order of the ODE is NOT 1, ode_order MUST be changed.
            
            evaluate=True
            parameters=None
            **options
        """


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
    def _constructor(cls,odes,dvars=None,odes_rhs=None , ivar=None,ode_order=None ,evaluate=True, parameters = None, **options):
        

        if not isinstance(dvars,Iterable):
            dvars = Matrix([dvars])
        
        obj = super()._constructor(elements=odes,vars=dvars,rhs=odes_rhs,evaluate=evaluate, **options)
        
        obj._parameters = parameters
        obj._dvars = dvars
        obj._const_list = []
        
        
        if ivar is not None:
            obj._ivar = ivar

        if ode_order is not None:
            obj._ode_order = ode_order

            
        return obj


    def _as_msm(self):
        ode=self
        from dynpy.solvers.nonlinear import MultiTimeScaleSolution
        
        msm_eq=MultiTimeScaleSolution.from_ode_system(ode)
        return msm_eq
    
    def default_ics(self,critical_point=False):
        
        
        
        ics_init_dict = {coord:0.0 for coord in self._fode_dvars}
        
        if isinstance(self._default_ics,dict):
            ics_instance={coord:self._default_ics[coord] for coord in self.dvars if coord in self._default_ics}
            
            return {**ics_init_dict,**ics_instance}
        else:
            return ics_init_dict



    

    
    
    def set_simp_deps(self,dependencies,callback=None,inplace=False):
        

        
        if inplace is False:
            obj = self.copy()
        else:
            obj = self
        
        obj._simp_dict  = dependencies
        
        if callback is not None:
            obj._callback_dict = callback
        
        #print(obj._simp_dict)
        return obj
    
    @property
    def _simp_deps(self):
        
        #display('_simp_deps',self._simp_dict)
        
        if self._simp_dict is None:
            
            Detector = self._default_detector
            
            if Detector is None:
                self._simp_dict = {}
            else:
                self._simp_dict = Detector(self._fundamental_matrix).subs_dict()


            return self._simp_dict
        else:
            return self._simp_dict

    @cached_property
    def _callback_deps(self):
        
        if self._callback_dict is None:
            
            Detector = self._default_detector
            
            if Detector is None:
                self._callback_dict = {}
            else:
                self._callback_dict = Detector(self._fundamental_matrix).callback_dict()
            
            return self._callback_dict
        else:
            return self._callback_dict
    
    @cached_property
    def ivar(self):
        return self._ivar
    
    @cached_property
    def dvars(self):
        return self._dvars
    
    @cached_property
    def _fode_dvars(self):

        orders_dict = self.spot_order()
        sys_order=max(orders_dict.values())

        # a quick workaround - error to be solved
        #print(sys_order)
        if sys_order == 1:
            dvars_span = list(self.dvars)
        else:
            dvars_span = ([dvar.diff(self.ivar,order) for order  in range(sys_order) for dvar in self.dvars if order < orders_dict[dvar]])
            dvars_span = list(({dvar:1 for dvar in dvars_span}).keys())
            #dvars_span = list(np.unique(dvars_span))

        return Matrix(dvars_span)
    
    @cached_property
    def _highest_diff(self):
        
        orders = self.spot_order()
        dvars = self.dvars
        
        return Matrix([dvar.diff(self.ivar,orders[dvar])  for dvar in dvars])
    
    @cached_property
    def _reduction_eqns(self):
        
        return [eqn for eqn  in list(self._fode_dvars) if eqn not in self.dvars]
    
    
    
    @cached_property
    def odes_rhs(self):
        
        diffs = self._highest_diff
        
        if diffs == self.lhs:
            return self.rhs
        else:
            diff_coeffs_mat=self.lhs.jacobian(diffs)

            return diff_coeffs_mat.inv()*(self.rhs  - self.lhs.subs({elem:0 for elem in diffs }))

    #@lru_cache
    def _to_rhs_ode(self):
        

        diffs = self._fode_dvars.diff(self.ivar)
    
        
        if diffs == self.lhs:
            return self.copy()
        else:
            
            free_terms = self._free_component()
            
            diffs = self._highest_diff
            
            if diffs == self.lhs:
                odes_rhs =  self.rhs
            else: #nonlinear cases have to be implemented
                diff_coeffs_mat=self.lhs.jacobian(diffs)

                odes_rhs = diff_coeffs_mat.inv()*free_terms
            
            rhs_odes_mat = odes_rhs
        
            display(self.lhs - self.rhs,diffs)
            # rhs_odes_dict=solve(self.lhs - self.rhs,diffs)
            # rhs_odes_mat = Matrix([rhs_odes_dict[coord]  for coord in diffs])
            display(rhs_odes_mat)
            

            
            
            new_sys = type(self)._constructor( diffs,
                                            dvars=self.dvars,
                                            odes_rhs=rhs_odes_mat,
                                            ivar=self.ivar,
                                            ode_order=self.ode_order,
                                            #evaluate=self._evaluate,
                                            parameters=self.parameters
                                            )
                    
            new_sys._const_list = self._const_list
            
            return new_sys


    @cached_property
    def parameters(self):
        
        return self._parameters

    @cached_property    
    def ode_order(self):
        
        #if self._ode_order is None:
        
        return self._ode_order
    
    @cached_property
    def odes(self):
        return Matrix([Add(lhs,-rhs,evaluate=False) if lhs-rhs == 0 else lhs-rhs   for lhs,rhs   in zip(self.lhs,self.rhs) ])


    def as_eq(self):
        return Eq(self.lhs,self.rhs,evaluate=False)
    
    
    @cached_property
    def _lhs_repr(self):
        #return (self.lhs).applyfunc(lambda obj: Derivative(obj,self.ivar,evaluate=False)  )
        return (self.lhs).applyfunc(lambda obj: obj.diff(self.ivar)  )

    @property
    def _default_solver(self):
        
        return FirstOrderLinearODESystemWithHarmonics
        #return FirstOrderLinearODESystem

    def _latex(self,*args):

        if self.default_ics is None:
            latex_str =  f'{latex(self.as_eq())}~~for~~{latex(self.dvars)}'
        else:
            ics = self.default_ics()
            latex_str = f'{latex(self.as_eq())}~~for~~{latex(self.dvars)}~~with~~{latex(ics)}'
    #
        if len(self._callback_dict) == 0:
            return f'{latex_str}'
        else:
            cannonical_coeffs_list=[Eq(thing,self._callback_dict[thing],evaluate=False) for thing in self._callback_dict]
            
            return f'{latex_str}~~where~~{latex(cannonical_coeffs_list[0])}, {latex(cannonical_coeffs_list[1])}'

    def as_matrix(self):
        return self.odes
    
    
    def _as_fode(self):
        """Creates an object of FirstOrderLinearODESystem class."""
        
        #display('subs dict',self._simp_dict,type(self))
        
        fode = FirstOrderLinearODESystem.from_ode_system(self)
        
        fund_mat=fode._fundamental_matrix
        # # # print(type(self))
        # display(aux_mat)
        # display(aux_mat.det())
        
        # # print('det != 0 ')
        # display(det(aux_mat).atoms(Function))
        
        #if Heaviside in det(aux_mat).atoms():
        mat_elems = (list(fund_mat))
        dummy_list = symbols(f'A0:{len(mat_elems)}')

        
        
        aux_mat = Matrix(*fund_mat.shape,[dummy_list[no] if elem !=0 else 0  for no,elem in enumerate(mat_elems)])
        
        if det(aux_mat).has(self.ivar):
            return fode
        elif det(aux_mat)!=0:
            _solver = self._default_solver
            # print(_solver.__name__, ' used as solver')
            return _solver.from_ode_system(self)
        else:
            
            # print('FirstOrderLinearODESystem used as solver')
            return fode
            
        
        
    def as_first_ode_linear_system(self):
        
        return self._as_fode()

    
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
        

        
        #obj = super().subs(*args,**kwargs)

        lhs = self.lhs.subs(*args,**kwargs)
        rhs = self.rhs.subs(*args,**kwargs)
        
        obj = type(self)._constructor(lhs,self.vars,rhs)
        
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    
    
    # def __add__(self,other):
        
    #     if isinstance(other,self.__class__):
    #         other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
    #     obj = super().__add__(other)
    #     obj._dvars=self._dvars
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)

    
    # def __rsub__(self,other):
        
    #     if isinstance(other,self.__class__):
    #         other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
    #     obj = super().__rsub__(other)
    #     obj._lhs=self._lhs
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)

    
    # def __sub__(self,other):
        
    #     if isinstance(other,self.__class__):
    #         other = Matrix([other.as_dict()[coord]  for  coord  in self._lhs ])
        
    #     obj = super().__sub__(other)
        
    #     obj._dvars=self._dvars
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)

    
    
    # def __mul__(self,other):
        
    #     obj = super().__mul__(other)
    #     obj._dvars=self._dvars
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    
    # def doit(self,**hints):
        
    #     obj = super().doit(**hints)
    #     obj._rhs = self.rhs.doit(**hints)
    #     obj._dvars=self._dvars
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)

    # def expand(self,deep=True, modulus=None, power_base=True, power_exp=True,mul=True, log=True, multinomial=True, basic=True, **hints):
        

    #     obj = super().expand(deep=deep, modulus=modulus, power_base=power_base, power_exp=power_exp,mul=mul, log=log, multinomial=multinomial, basic=basic,**hints)
    #     obj=self.from_ode_system(self)
    #     obj._rhs = self.rhs.expand(deep=deep, modulus=modulus, power_base=power_base, power_exp=power_exp,mul=mul, log=log, multinomial=multinomial, basic=basic,**hints)
    #     obj._dvars=self._dvars
    #     obj._ivar = self.ivar
    #     obj._ode_order = self.ode_order
        
    #     return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    
    
    def copy(self):
        
        obj=self.from_ode_system(self)
        #obj = super().copy()
        obj._lhs = copy.copy(self.lhs)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        obj._const_list = self._const_list
        
        return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    
    
    def _eval_applyfunc(self, f):
        
        obj = super()._eval_applyfunc(f)
        obj._lhs = self.lhs._eval_applyfunc(f)
        obj._dvars=self._dvars
        obj._ivar = self.ivar
        obj._ode_order = self.ode_order
        
        
        return obj.set_simp_deps(self._simp_dict,self._callback_dict,inplace=True)
    

    def is_linear(self):
        
        if isinstance(self,self.linearized()):
            return print('True')
        else:
            return print('False')

# #     @lru_cache
#     def numerized(self,parameters=None,ic_list=[],backend='fortran',**kwrags):
#         '''
#         Takes values of parameters, substitutes it into the list of parameters and changes it into a Tuple. Returns instance of class OdeComputationalCase.
#         '''

        
#         ode=self.as_first_ode_linear_system()
#         if parameters is None:
#             parameters={}
        
#         return OdeComputationalCase(odes_system=ode.rhs,ivar=ode.ivar,dvars=ode.dvars,params= parameters,backend=backend)
    
    # def numerized(self,parameters={},ic_list=[]):
    #     '''
    #     Takes values of parameters, substitutes it into the list of parameters and changes it into a Tuple. Returns instance of class OdeComputationalCase.
    #     '''
    #     return OdeComputationalCase(odes_system=self.odes_rhs,dvars=self.dvars,ivar=self.ivar)

    def numerized(self,parameter_values=None,ic_list=None,backend='fortran',**kwrags):
        '''
        Takes values of parameters, substitute it into the sympy Dict. Redirects the numerizing to exectution method _numerized which has lru cache.
        '''

        if parameter_values is None: parameters_sympy_dict = Dict({})
        elif isinstance(parameter_values, dict): parameters_sympy_dict = Dict(parameter_values)
        elif isinstance(parameter_values, Dict): parameters_sympy_dict = parameter_values
        else: print("Podano zły typ danych - podaj słownik na parameter_values")

        if ic_list is None: ic_tuple=()
        elif isinstance(ic_list, list): ic_tuple = tuple(ic_list)
        elif isinstance(ic_list, tuple): ic_tuple = ic_list
        else: print("Podano zły typ danych - podaj list lub tuple na ic_list")

        return self._numerized(parameters=parameters_sympy_dict, ic_tuple=ic_tuple, backend=backend)

    #@lru_cache
    def _numerized(self,parameters=Dict(),ic_tuple=(),backend='fortran',**kwrags):
        '''
        Execution method that assumes all Args are hashable due to application of lru cache. Takes values of parameters and returns instance of class OdeComputationalCase.
        '''

        ode=self.as_first_ode_linear_system()._to_rhs_ode()

        return OdeComputationalCase(odes_system=ode.rhs,ivar=ode.ivar,dvars=ode.dvars,params= parameters,backend=backend)
    
    
    @property
    def _report_components(self):
        
        comp_list=[
        *ODE_COMPONENTS_LIST
        ]
        
        return comp_list

    @cached_property
    def _general_solution(self):

        fode = self._as_fode()
        #solver changing due to the computational simplicity reasons

        fode_sys= fode #FirstOrderLinearODESystem(fode,fode.dvars)        

        return fode_sys.general_solution

    @cached_property
    def _steady_solution(self):

        fode = self._as_fode()
        #solver changing due to the computational simplicity reasons

        fode_sys= fode #FirstOrderLinearODESystem(fode,fode.dvars)        

        return fode_sys.steady_solution



    
    
    @cached_property
    def general_solution(self):
        '''
        Solves the problem in the symbolic way and returns matrix of solution (in the form of equations (objects of ODESolution class)).
        '''
        return self._general_solution


    @cached_property
    def steady_solution(self):
        '''
        Provides the particular solution of a differential equation and returns matrix of solution (in the form of equations (objects of ODESolution class)).
        '''
        return self._steady_solution
    
    @property
    def solution(self):
        '''
        Provides the final solution of a differential equation and returns matrix of solution (in the form of equations (objects of ODESolution class)).
        '''
        return self.general_solution + self.steady_solution    
    

    # Derivative reset Method
    def static_equation(self):
        obj = self
        derivatives = list(obj.lhs.atoms(Derivative))
        derivatives_values = {item: 0 for item in derivatives}
        static_eqs = obj.subs(derivatives_values)
        return static_eqs


    def critical_point(self, params):
        """
        Calculating a critical point
        """
        obj = self
        critical_point = solve(Eq(obj.lhs, obj.rhs), self.dvars[0])
        critical_point = critical_point[self.dvars[0]]
        return Eq(self.dvars[0], critical_point.subs(params))


    # Designate Order Method
    def spot_order(self):
        
        
        odes = self.odes
        dvars = self.dvars
        
        
        orders={coord:max(odes.applyfunc(lambda elem: ode_order(elem,coord)))  for coord in dvars }

        return orders
    
    
    @property
    def details(self):
        obj = self
        fode = self._as_fode()

        display(obj)
        display(Markdown(f'Ivar: {self._ivar}'))
        display(Markdown(f'ODE order: {self._ode_order}'))
        display(Markdown(f'System reduced to first order:'), fode)
        
        
    def _free_component(self):
        #type(self.odes)
        return self.odes.subs({var:0 for var in self.dvars}).doit()
        #return (self.lhs - self.rhs).subs({var:0 for var in self.dvars}).doit()

    def _free_comps_by_sides(self):
        #type(self.odes)
        return self.subs({var:0 for var in self.dvars}).doit()
    
    def _hom_equation(self):
        
        # free_terms=self._free_component()
        ft_sides = self._free_comps_by_sides()
        # print('_hom_eqation start \n free terms')
        
        # display(free_terms)
        # display(ft_sides)
        
        # print('_hom_eqation : lhs')
        # display(self.lhs)
        # print('_hom_eqation : rhs')
        # display(self.rhs )
        # print('_hom_eqation - free : rhs')
        # display(self.rhs + free_terms)
        
        # print('ft_sides: lhs + rhs')
        # display(ft_sides.lhs,ft_sides.rhs)
        
        
        #hom_ode = ODESystem( self.lhs , dvars = self.dvars,odes_rhs=self.rhs + free_terms, ode_order = self._ode_order, ivar = self._ivar)
        hom_ode_sides = ODESystem( self.lhs-ft_sides.lhs , dvars = self.dvars,odes_rhs=self.rhs - ft_sides.rhs, ode_order = self._ode_order, ivar = self._ivar)
        
        # print('home composed')
        # display(hom_ode)
        # display(hom_ode_sides)
        # print('home composed end')
        
        return hom_ode_sides #dzieki Madi<3
    
    def char_polynomial(self):
        r=Symbol('r')
        stiffness_mat=self.lhs.jacobian(self.dvars)
        damping_mat=self.lhs.jacobian(diff(self.dvars))
        inertia_mat=self.lhs.jacobian(diff(self.dvars, self.ivar, 2))
        matrix_poly = inertia_mat*r**2 + damping_mat*r + stiffness_mat
        return matrix_poly.det()

    def _stiffness_matrix(self):
        return self.lhs.jacobian(self.dvars)
    
    def _inertia_matrix(self):
        return self.lhs.jacobian(diff(self.dvars,self.ivar,2))
    
    def _eigenvals_simplified(self):
        
        base_matrix = self._inertia_matrix().inv()*self._stiffness_matrix()
        base = simplify(base_matrix)
        base_matrix = base/base[0]

        nom = fraction(((base_matrix.trace()/2)**2 - base_matrix.det()).simplify().radsimp())[0]
        denom = fraction(((base_matrix.trace()/2)**2 - base_matrix.det()).simplify().radsimp())[1]
        
        delta_omg = sqrt(nom)/sqrt(denom)

        omg_mean = (base_matrix.trace()/2).simplify()
        omg_mean_nom = fraction(omg_mean)[0]
        omg_mean_denom = fraction(omg_mean)[1]

        omg12_omg22 = base_matrix[1]*-base_matrix[3]

        omg_1 = Symbol('\omega_{n_{1}}')
        omg_2 = Symbol('\omega_{n_{2}}')
        omg_n = Symbol('\omega_{n}')


        delta_omg_simp = omg_mean**4 - simplify(base_matrix)[1]*simplify(base_matrix)[2]


#         eigenval1 = sqrt(omg_mean-delta_omg)
#         eigenval2 = sqrt(omg_mean+delta_omg)

#         eigenval1 = omg_n*sqrt(1 - sqrt(1-(omg12_omg22)/omg_n**2))
#         eigenval2 = omg_n*sqrt(1 + sqrt(1-(omg12_omg22)/omg_n**2))
#        display(eigenval1)

#         display(omg_mean_nom)
#         display(omg_mean_denom)
        eigenval1 = Mul(sqrt(1 - sqrt(1-(omg12_omg22)/omg_mean**2)),Mul(omg_mean_nom,1/omg_mean_denom,evaluate=False),evaluate=False)
        eigenval2 = Mul(sqrt(1 + sqrt(1-(omg12_omg22)/omg_mean**2)),Mul(omg_mean_nom,1/omg_mean_denom,evaluate=False),evaluate=False)

        return Matrix([[eigenval1,0],[0,eigenval2]])



        
    def fundamental_matrix(self):
        r=Symbol('r')
        stiffness_mat=self.lhs.jacobian(self.dvars)
        damping_mat=self.lhs.jacobian(diff(self.dvars))
        inertia_mat=self.lhs.jacobian(diff(self.dvars, self.ivar, 2))
        fundamental_mat = stiffness_mat + r**2 * inertia_mat + damping_mat * r
        return fundamental_mat

    def _swept_analysis(self, subs_method=True, dvar=None, freq_symbol=None, amp_symbol=None, amplitude=None, ramp=None):

        '''
        Performs swept analysis on numerical ODESystem.

        subs_method - Boolean, switches methods for swept analysis, if True it substitutes values, if False it uses _hom_equation method.
        
        dvar - Symbol or Integer, coordinate which will be taken into consideration for performing swept analysis. This coordinate indicates which equation will be under excitation.

        freq_symbol - Symbol for substituting the excitation frequency value in equation.

        amp_symbol - Symbol for substituting the excitation amplitude value in equation.

        amplitude - Integer or Float, amplitude of excitation.

        ramp - Integer or Float, sets the steep coefficient of the ramp.
        '''
        if ramp is None:
            frequency = 0.005*self._ivar
        elif isinstance(ramp, (float, int)):
            frequency = ramp*self._ivar
#             display(frequency)
        else:
            return "Error, ramp should be a float or int type"

        if amplitude is None:
            amplitude=1
        elif isinstance(amplitude, (float, int)):
            pass
        else:
            return "Error, Amplitude should be a float or int type"

        if subs_method is True:

            if freq_symbol is None:
                freq_symbol = Symbol('Omega', positive=True)
            if amp_symbol is None:
                amp_symbol = Symbol('F', positive=True)
#             display(freq_symbol,amp_symbol)
            return type(self)(odes=self.odes.subs({freq_symbol: frequency, amp_symbol:amplitude}), dvars=self.dvars, ode_order=self._ode_order, ivar=self._ivar)#.numerized(backend='numpy')

        else:

            hom_ode = self._hom_equation()

            if dvar is None:
                dvar = self.dvars[0]

            if isinstance(dvar, Symbol):
                for no, sym in enumerate(self.dvars):
                    if sym == Symbol:
                        number = no

            elif isinstance(dvar, int):
                number = dvar

            else:
                return "Error"


            exct_mat = zeros(len(self.dvars),1)
            exct_mat[number,0] = amplitude*sin(Omega*self._ivar)


            return type(self)(odes=hom_ode.odes, odes_rhs=exct_mat, dvars=self.dvars, ode_order=self._ode_order, ivar=self._ivar)#.numerized(backend='numpy')


#         type(self)(odes=hom_eq, odes_rhs=zeros(len(self.dvars),1), dvars=self.dvars, ode_order=self._ode_order, ivar=self._ivar)

    def swept_analysis(self, subs_method=True, dvar=None, freq_symbol=None, amp_symbol=None, amplitude=None, ramp=None):

        return self._swept_analysis(subs_method=subs_method, dvar=dvar, freq_symbol=freq_symbol, amp_symbol=amp_symbol, amplitude=amplitude, ramp=ramp)

    def to_canonical_form(self):
        if self._ode_order == 2 and len(self) == 1:
            
            coord = self.dvars[0]
            ivar = self.ivar
            
            inertia = self.odes[0].coeff(diff(coord,self.ivar,2))
            
            new_eq = (self.odes[0]/inertia).expand()
            rhs_comp=self._free_component()/inertia
            
            omega_coeff = new_eq.coeff(self.dvars[0])
            h_coeff = new_eq.coeff(diff(self.dvars[0], self.ivar))
            
            omega_sym = Symbol('omega',positive=True)
            h_sym  = Symbol('h',positive=True)
            
            subs_dict = {omega_sym**2:omega_coeff, h_sym:S.Half * h_coeff}
            
            #subs_eq = new_eq.subs(subs_dict)
            canonical_ode = Matrix([coord.diff(ivar,2) + 2*h_sym * coord.diff(ivar) + omega_sym**2 * coord  ] )  +rhs_comp
            
            new_odesys = ODESystem(odes = canonical_ode, dvars = Matrix([coord]), ode_order = 2)
            
            new_odesys._callback_dict = subs_dict
            
            return new_odesys
        else:
            return self.copy()

    def _classify_ode(self):
        odes=Eq(self.lhs[0]-self.rhs[0],0)
        return classify_ode(odes,self.dvars[0])
### METODY - Karolina, do przejrzenia i wgl do zastanowienia czy chcemy to czy nie, na pewno dvars symbols sa must have no i te dwie ostatnie one ciagna z fode wsm##
    def _eqn_coeffs(self):

        system = self
        ivar = self.ivar
        dvar = system.dvars[0]

        coeffs_dict = {dvar.diff(ivar,2): system.lhs[0].coeff(dvar.diff(ivar,2)), #2nd derivative
                       dvar.diff(ivar): system.lhs[0].coeff(dvar.diff(ivar)), #1st derivative
                       dvar : system.lhs[0].coeff(dvar), #y
                       'free': system.lhs[0].subs({dvar:0})
                        }

        return coeffs_dict

    def _ode_solution(self):
        dvars=self.dvars
        dsol=dsolve(self.odes[0],self.dvars[0]).rhs
        return ODESolution(dvars,dvars,dsol)

    def _get_dvars_symbols(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar

        dvars_str=f"{str(dvars)}".replace(f'({str(ivar)})','')

        return {dvars:Symbol(dvars_str),
                ivar:ivar,
                dvars.diff(ivar):Symbol(f"{dvars_str}'"),
                dvars.diff(ivar,2):Symbol(f"{dvars_str}''")}

    def matrix_fundamental(self):
        return self._as_fode()._fundamental_matrix

    def _eigenvalues(self):
        return self._as_fode()._auxiliary_fundamental_matrix.diagonalize()[1]

    def _roots_list(self):
        return list(self._eigenvalues())
    
## koniec metod karoliny pozdraiwam ##

    @staticmethod
    def _get_code():
        from dynpy.utilities.report import ObjectCode
        
        code = ObjectCode('''
            from dynpy.solvers.linear import ODESystem
            from sympy import *

            F = Symbol('F', positive=True)
            Omega = Symbol('Omega', positive=True)
            t = Symbol('t')
            c = Symbol('c', positive=True)
            k = Symbol('k', positive=True)
            m = Symbol('m', positive=True)
            x = Function('x')
            eq = -F*sin(Omega*t) + c*Derivative(x(t), t) + k*x(t) + m*Derivative(x(t), (t, 2))

            dvars = Matrix([[x(t)]])
            odes = Matrix([[-F*sin(Omega*t) + c*Derivative(x(t), t) + k*x(t) + m*Derivative(x(t), (t, 2))]])
            ode_order = 2

            odesys = ODESystem(odes = odes , dvars = dvars, ode_order = ode_order)

            display(odesys)
                         
                         ''')
        return code
    
    def as_type(self, ode_type):        
        return ode_type.from_ode_system(self) 
    
class FirstOrderODESystem(ODESystem):
    
    
    @classmethod
    def from_ode_system(cls,ode_system):

        
        
        sys = ode_system
        odes=sys.lhs
        odes_rhs = sys.rhs
        dvars = sys.dvars
        ivar = sys.ivar
        ode_order = 1
        parameters = sys._parameters
        

        if sys.ode_order == 1:
            new_sys = cls._constructor(odes , dvars,odes_rhs, ivar = ivar,ode_order=1 ,parameters=parameters)
        else:
            f_ord_dvars=sys._fode_dvars
            aux_odes = sys._reduction_eqns
            ode_rhs = Matrix(aux_odes + list(sys.odes_rhs))
            new_sys = cls._constructor(f_ord_dvars.diff(ivar) , f_ord_dvars,ode_rhs, ivar = ivar,ode_order=1 ,parameters=parameters)
        # print(f'f_ord_dvars from from_ode_system by {cls}')
        # print(f'base system {type(sys)}')
        # display(f_ord_dvars)
        
        # print(f'sys rhs')
        # display(ode_rhs)
    
        
        
        #display('subs dict',sys._simp_dict,cls)
        new_sys._const_list = sys._const_list
        
        return new_sys.set_simp_deps(sys._simp_dict,sys._callback_dict,inplace=True)
    
    def _as_fode(self):
        """Creates an object of FirstOrderLinearODESystem class."""
        
        return self.copy()
    
    @cached_property
    def solution(self):
        return self.linearized().solution


        

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


#     def subs(self,*args,**kwargs):
        

        
#         obj = super().subs(*args,**kwargs)
#         obj._lhs=self._lhs.subs(*args,**kwargs)
#         obj.
        
#         return obj
    

    
class FirstOrderLinearODESystem(FirstOrderODESystem):
    #_const_list = []
    
    @classmethod
    def from_odes(cls,odes_system,dvars,ivar=None,ode_order=1,parameters=None):
        
        
        diffs = dvars.diff(ivar,ode_order)

        diffs_coeffs_mat  = odes_system.jacobian(diffs) #it should be reimplemented with .jacobian in Analytical Solution class
        
        if diffs_coeffs_mat == Matrix.zeros(len(diffs)):
        
            odes = odes_system
        else:
            odes = diffs_coeffs_mat.inv() * (-odes_system.subs({ d_dvar:0     for d_dvar in diffs } ))
            
        f_ord_dvars=Matrix(sum([list(dvars.diff(ivar,no))   for no in  range(ode_order)],[]))
        
        #display(odes)
        #display(f_ord_dvars)
            
        return cls._constructor( odes=f_ord_dvars.diff(ivar) , dvars=f_ord_dvars  , odes_rhs = Matrix([f_ord_dvars[len(dvars):],odes] )  ,ivar=ivar,ode_order=ode_order,parameters=parameters)


    
    
    @cached_property
    def odes_rhs(self): #hmmmm DAmian/Madi możemy zaryzykować tu w sumie
        return self.rhs
        #   return self._to_rhs_ode().rhs
    
    @cached_property
    def _fundamental_matrix(self):
        
        CodeFlowLogger(self._to_rhs_ode().rhs,"_to_odes_rhs",self)  
        
        CodeFlowLogger(self.rhs,"odes matrix lhs",self)       
        
        A = -self._to_rhs_ode().rhs.jacobian(self.dvars)

        CodeFlowLogger(A,"fundamental_matrix",self)
        
        return A
        #return self.odes_rhs.jacobian(self.dvars)
    
    
    @cached_property
    def _free_terms(self):
        
#         return self._to_rhs_ode().odes_rhs.subs({coord:0 for coord in self.dvars})
        return self.odes_rhs.subs({coord:0 for coord in self.dvars})



    @cached_property
    def _auxiliary_dvars(self):
        
        dvars = list(reversed(list(self.dvars)))
    
        
        return Matrix(dvars)  

    @cached_property
    def _auxiliary_fundamental_matrix(self):
        
        dvars = list(reversed(list(self.dvars)))
#         odes = list(reversed(list(self._to_rhs_ode().odes_rhs)))
        odes = list(reversed(list(self.odes_rhs)))
    
        
        return Matrix(odes).jacobian(dvars)  
    

    
    @cached_property
    def _auxiliary_free_terms(self):


        dvars = self.dvars

#         odes_list = list(reversed(list(self._to_rhs_ode().odes_rhs)))
        odes_list = list(reversed(list(self.odes_rhs)))

        odes = Matrix(odes_list)

        
        return odes.subs({coord:0 for coord in dvars})   
    
    def _const_mapper(self,const_set,parameters=None):
        
        if parameters is None:
            parameters = self.parameters
        
        const_base_dict={dummy_symbol : Symbol(f'C{no+1}')   for  no,dummy_symbol  in enumerate(const_set)}
        
        
        
        
        if parameters is None:
            const_dict=const_base_dict
            
            #const_dict=solve([key-val  for  key,val in const_base_dict.items()] , dummies_set )
            
        else:
            const_fun_dict = {dummy_symbol : Function(str(c_symbol))   for  dummy_symbol,c_symbol  in const_base_dict.items()}
            const_dict = {dummy_symbol : c_fun(*parameters)   for  dummy_symbol,c_fun  in const_fun_dict.items()}
            
            #display([key-val  for  key,val in const_dict.items()],list(const_dict.values()))
            
            #const_dict=solve([key-val  for  key,val in const_dict.items()] , dummies_set  )
            
            
        # print('const mapper is worong')
        # display(const_dict)
            
        self._const_list = list(const_dict.values())
        
        # print('const list status')
        # print(self._const_list)
        
        

        
        return const_dict

    @cached_property
    def solution_map(self):
        return lambda obj: TR8(TR10(TR8(obj.doit(conds='none')).expand()).expand())
        
    
    @cached_property
    def _general_solution(self):

        
        
        
        A = self._auxiliary_fundamental_matrix#.applyfunc(lambda elem: elem.subs(self._simp_deps,simultaneous=True))
        
        CodeFlowLogger(A,"_auxiliary_fundamental_matrix",self)
        
        dvars = self._auxiliary_dvars
        CodeFlowLogger(dvars,"aux dvars",self)
        
        
        #sol = AnalyticalSolution.from_vars_and_rhs(dvars,linodesolve(A,t=self.ivar,b=0*self.dvars)  )#.applyfunc(self.solution_map)
        sol = Matrix(linodesolve(A,t=self.ivar,b=0*self.dvars))
        
        CodeFlowLogger(sol,"sol AS",self)
        
        dummies_set= sol.atoms(Dummy)
        const_dict = self._const_mapper(dummies_set)


        sol = list(reversed(list(sol)))

        #print('self._callback_deps')
        #display(self._callback_deps)
        
        mapper = lambda elem: elem.subs(self._simp_deps,simultaneous=True).subs(self._callback_deps,simultaneous=True) 
        
        CodeFlowLogger(sol,"sol list",self)
        
        ode_sol = ODESolution.from_vars_and_rhs(self.dvars,sol).applyfunc(mapper).subs( const_dict )
        
        ode_sol.ivar=self.ivar
        ode_sol.append_integration_consts(list(const_dict.values()))
        
        CodeFlowLogger(ode_sol,"ODE sol with Sympy",self)
        
        return ode_sol
    
    
    @cached_property
    def _steady_solution(self):
        
        
        
        A = self._auxiliary_fundamental_matrix
        b = self._auxiliary_free_terms
        dvars = self._auxiliary_dvars                                 

        #print('steady sol')
        #display(A)
        #display(b)
        
        sol = Matrix(linodesolve(A,t=self.ivar,b=b)).applyfunc(self.solution_map)

        #display('ss',sol)
        
        dummies_set= sol.atoms(Dummy)


#         if len(self.dvars) >= 2:
#             no = int(len(self.dvars)/2)
        
#             sol=  sol[no:] + sol[:no] 
        sol = list(reversed(list(sol)))
    


        return ODESolution.from_vars_and_rhs(self.dvars,sol).subs({dum_sym:0 for dum_sym in dummies_set})

    @cached_property
    def const_set(self):

        return self._const_list
    
    
    @cached_property
    def general_solution(self):

        rhs_ode=self._to_rhs_ode() # it brings some errors - some _const_list was mapped
                                    # it should be reimplemented within scope of internal methods
        
        general_solution = rhs_ode._general_solution
        
        self._const_list = rhs_ode._const_list

        return general_solution



    @cached_property
    def steady_solution(self):

        rhs_ode=self._to_rhs_ode() # it brings some errors - some _const_list was mapped
                                    # it should be reimplemented within scope of internal methods
                                    
        steady_solution = rhs_ode._steady_solution
        
        self._const_list = rhs_ode._const_list

        return steady_solution
    
    @property
    def solution(self):
        
        gen_sol = self.general_solution
        steady_sol = self.steady_solution
        
        sol = gen_sol+steady_sol
        sol.append_integration_consts(gen_sol._spot_constant())
        
        return sol
    

    def _latex(self,*args):

#         f'{latex(self.as_eq())}~~for~~{latex(self.dvars)}'
        return latex(Eq(self.lhs,self.rhs ,evaluate=False   ))    
    
class FirstOrderLinearODESystemWithHarmonics(FirstOrderLinearODESystem):

    @cached_property
    def _auxiliary_fundamental_matrix(self):
        
        dvars = list(reversed(list(self.dvars)))
#         odes = list(reversed(list(self._to_rhs_ode().odes_rhs)))
        odes = list(reversed(list(self.odes_rhs)))
    
        
        return Matrix(odes).jacobian(dvars)  
    
    @cached_property
    def _reduced_fundamental_matrix(self):

        fund_mat=Matrix(self._auxiliary_fundamental_matrix)


        damping = self._is_proportional_damping
        if damping:
            #display('the damping is',damping)
            sl = int(fund_mat.shape[0]/2)
            fund_mat[:sl,:sl] = 0*fund_mat[:sl,:sl]
            
            return fund_mat#.subs(damping,0)
        else:
            #display(f'there is {damping} damping here xD')
            return fund_mat
        
    @cached_property
    def _trig_fundamental_matrix(self):

        fund_mat=Matrix(self._reduced_fundamental_matrix)
        sl = int(fund_mat.shape[0]/2)
        
        #display(fund_mat)
        
        fund_mat[sl:,:] = (-1)*fund_mat[sl:,:]
        
        #display(fund_mat)
        

        return fund_mat
        

    @cached_property
    def _reduced_eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''


        return self._reduced_fundamental_matrix.diagonalize()[1]
  
    @cached_property
    def _trig_eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        return self._trig_fundamental_matrix.diagonalize()[1]


    @cached_property
    def _reduced_modes(self):
        '''
        Returns reversed modes matrix (computed by Sympy's .diagonalize() method) by changing order of all rows.
        '''

        modes=self._reduced_fundamental_matrix.diagonalize()[0]
        rows_no = modes.shape[0]
        
        rows_list = [modes[row,:] for row in reversed(range(rows_no))]

        return Matrix(rows_list)

    
    @cached_property
    def _combined_modes(self):
        '''
        Returns reversed modes matrix (computed by Sympy's .diagonalize() method) by changing order of all rows.
        '''

        r_modes=Matrix(self._reduced_modes)
        t_modes=Matrix(self._trig_modes)
        
        eig_fun=self.eigenfunctions()
        
        for no,efun in enumerate(eig_fun):
            if efun.atoms(sin,cos) != set():
                
                r_modes[:,no] = t_modes[:,no]
                
        return r_modes
    
    
    @cached_property
    def _trig_modes(self):
        '''
        Returns reversed modes matrix (computed by Sympy's .diagonalize() method) by changing order of all rows.
        '''

        modes=self._trig_fundamental_matrix.diagonalize()[0]
        rows_no = modes.shape[0]
        
        rows_list = [modes[row,:] for row in reversed(range(rows_no))]

        return Matrix(rows_list)
    
    @cached_property
    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''


        return self._auxiliary_fundamental_matrix.diagonalize()[1]
  
    @cached_property
    def modes(self):
        '''
        Returns reversed modes matrix (computed by Sympy's .diagonalize() method) by changing order of all rows.
        '''

        modes=self._auxiliary_fundamental_matrix.diagonalize()[0]
        rows_no = modes.shape[0]
        
        rows_list = [modes[row,:] for row in reversed(range(rows_no))]

        return Matrix(rows_list)

    
    @cached_property
    def _is_proportional_damping(self):
        
        A_mat=self._fundamental_matrix
        
        
        aux_size=N((len(self.dvars))/2)
        
        #display(aux_size)
        sl = int(aux_size)
        
        if aux_size==round(aux_size):

            
            block_mats = A_mat[:sl,sl:],A_mat[:sl,:sl]
            
            is_reduced=block_mats[0]== eye(sl) and block_mats[1]== zeros(sl)
        else:
            is_reduced = False

        if is_reduced:
            K_mat=A_mat[sl:,:sl]
            C_mat=A_mat[sl:,sl:]
            
            lam=(C_mat[0,0]/K_mat[0,0])
            
            #display(K_mat,C_mat)
            
            if C_mat == K_mat*lam and lam != 0:
                return lam
            else:
                return False

        else:
            return False
        

    
    
    @cached_property
    def _general_solution(self):

        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        
        
        
        C = numbered_symbols('C', start=1)
        C_list = []

        for i in range(len(self.dvars) ):
            C_list += [next(C)]

        args_list = self.dvars[0].args

        if len(self.dvars[0].args) > 1:

            params = {*args_list} - {self.ivar}

            C_list = [Function(str(C_tmp)) for C_tmp in C_list]
            C_list = [(C_tmp)(*params) for C_tmp in C_list]
                
                
        
        fund_det=self._fundamental_matrix.det()

        
        if fund_det!=0:


            #         print('o tu')
            #         display(self.odes_system)

            #self.__class__._const_list |= set(C_list)
            damping = self._is_proportional_damping
            if damping is not False:
                #print('rayleight damping')
                modes = self._combined_modes

            else:
                #print('rayleight damping zero')
                modes = self.modes




            Y_mat = Matrix(self.dvars)

            
            
            
            #display(self.eigenfunctions())
            solution = self._combined_modes*Matrix([C_list[no]*self.eigenfunctions()[no] for   no   in range(len(self.dvars))]  )#.applyfunc(lambda elem: elem.rewrite(sin))
            
            if damping is not False:
                sl=int(N((len(self.dvars))/2))
                solution[sl:,:] = solution[:sl,:].diff(self.ivar)
            
            
            
        else:
            
            A = self._fundamental_matrix
            solution = (matrix_exp(A, self.ivar)*Matrix( C_list  ))

        ode_sol = ODESolution.from_vars_and_rhs(self.dvars,solution)
        ode_sol.ivar=self.ivar
        self._const_list = C_list
        ode_sol.append_integration_consts(C_list)    
        
        
        return ode_sol
    
    def eigenfunctions(self):
        
        
        damping = self._is_proportional_damping
        if damping is not False:
            #print('rayleight damping - new code')
            modes = self._reduced_modes
            modes_trig = self._trig_modes
            eigs = self._reduced_eigenvalues
        else:
            #print('rayleight damping zero')
            modes = self.modes
            eigs = self.eigenvalues

        #display(eigs)

        Y_mat = Matrix(self.dvars)

        eig_fun = []
        for   no   in range(len(self.dvars)):

            
            if eigs[no,no].n().is_imaginary:
                
                eigs_trig=self._trig_eigenvalues
                
                if damping:
                    h = damping/2*eigs_trig[no,no]**2
                else:
                    h = 0
                

                
                omg  = sqrt((eigs_trig[no,no])**2 - h**2)
                
                #display(h)
                #display(omg)
                
                if sin(no*pi/2) == 0:
                    gen_fun = cos
                else:
                    gen_fun = sin
                
                fun = (exp(-h*self.ivar)*(gen_fun(omg*self.ivar) ))#.expand()
            else:
                fun = exp(eigs[no,no]*self.ivar)
        
            eig_fun += [fun]
        
        return eig_fun

    ############################### NEW code
    
    def _sin_comp(self,omega,amp): 
        '''
        It applies generic form solution for the following differential equation
        \dot Y  =  A Y  + F \sin(\Omega t)
        
        The generic form is:
        
        D = -(A^{-1} \Omega^2 + A)^{-1}  F 
        C =   \Omega A^{-1} * D
        '''

        
        damping = self._is_proportional_damping
        #print('damping = ',damping)
        
        if damping is not False:
            
            
            sl=int(N((len(self.dvars))/2))
            #print('rayleight damping - new code stedy')
            modes = self._trig_modes[:sl,::2]
            modes_trig = self._trig_modes
            eigs = [self._trig_eigenvalues[no,no] for no  in range(len(self.dvars))[::2]]
            
            #display(eigs)
            
            
            
            b=Matrix(amp)[sl:,:]
            
            base_sin_list=[ (eigen**2 - omega**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            base_cos_list=[ -(damping*omega*eigen**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]

            
            
            modes_inv=modes.inv()
            modes_x_b=modes*b
            
            cos_comp_base = modes_inv*diag(*base_cos_list)*modes_x_b
            sin_comp_base = modes_inv*diag(*base_sin_list)*modes_x_b            
            
            cos_comp=Matrix(cos_comp_base).row_insert(sl,sin_comp_base*omega)
            sin_comp=Matrix(sin_comp_base).row_insert(sl,-cos_comp_base*omega)
            
        else:
            #print('rayleight damping zero - new code stedy')
            modes = self.modes
            eigs = self.eigenvalues
        


            A = self._fundamental_matrix
            b = amp


            sin_comp = -(A.inv() * omega**2 + A).inv()*b
            cos_comp = +omega*A.inv() * sin_comp
            
        return (cos_comp*cos(omega*self.ivar) +  sin_comp*sin(omega*self.ivar) ) 
    
    
    def _cos_comp(self,omega,amp):
        
        '''
        It applies generic form solution for the following differential equation
        \dot Y  =  A Y  + F \cos(\Omega t)
        
        The generic form is:
        
        C = -(A^{-1} \Omega^2 + A)^{-1}  F 
        D =  -\Omega A^{-1} * C
        '''

        
        damping = self._is_proportional_damping
        #print('damping = ',damping)
        
        if damping is not False:

            sl=int(N((len(self.dvars))/2))
            #print('rayleight damping - new code stedy')
            modes = self._reduced_modes[:sl,::2]
            modes_trig = self._trig_modes
            eigs = [self._trig_eigenvalues[no,no] for no  in range(len(self.dvars))[::2]]
            
            #display(eigs)
            

            b=Matrix(amp)[sl:,:]

            
            base_sin_list=[ (damping*omega*eigen**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            base_cos_list=[ (eigen**2 - omega**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            
            #display(self._reduced_modes)
            #display(base_cos_list)
            
            modes_inv=modes.inv()
            modes_x_b=modes*b
            
            cos_comp_base = modes_inv*diag(*base_cos_list)*modes_x_b
            sin_comp_base = modes_inv*diag(*base_sin_list)*modes_x_b            
            
            cos_comp=Matrix(cos_comp_base).row_insert(sl,sin_comp_base*omega)
            sin_comp=Matrix(sin_comp_base).row_insert(sl,-cos_comp_base*omega)
            
        else:
            #print('rayleight damping zero - new code stedy')
            modes = self.modes
            eigs = self.eigenvalues
        
            A = self._fundamental_matrix
            b = amp

            cos_comp = -(A.inv() * omega**2 + A).inv()*b
            sin_comp = -omega*A.inv() * cos_comp
          
        return cos_comp*cos(omega*self.ivar) +  sin_comp*sin(omega*self.ivar)  

    
    def _heaviside_sin_comp(self, exct, n=None):
    
        if type(exct) == Heaviside:
    
            for arg in exct.args:
                if type(arg) == sin:
                    sin_coeff = arg.args[0].diff(self.ivar)
    
        T = 2*pi/sin_coeff
    
        new_expr = 0
    
        if n is None:
            n=6
    
        for no in range(n):
    
            new_expr = new_expr + (-1)**no * Heaviside(self.ivar - T*no/2)
    
        return new_expr
    
    def _exp_sin_comp(self,omega,amp,a):

        damping = self._is_proportional_damping
        if damping is not False:
    
            sl=int(N((len(self.dvars))/2))
            #print('rayleight damping - new code stedy')
            modes = self._reduced_modes[:sl,::2]
            modes_trig = self._trig_modes
            eigs = [self._trig_eigenvalues[no,no] for no  in range(len(self.dvars))[::2]]
            
            #display(eigs)
    
            b=Matrix(amp)[sl:,:]
    
            
            base_exp_sin_list=[ (damping*omega*eigen**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            base_exp_cos_list=[ (eigen**2 - omega**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            
            
            modes_inv=modes.inv()
            modes_x_b=modes*b
            
            exp_cos_comp_base = modes_inv*diag(*base_exp_cos_list)*modes_x_b
            exp_sin_comp_base = modes_inv*diag(*base_exp_sin_list)*modes_x_b
            
            exp_cos_comp=Matrix(exp_cos_comp_base).row_insert(sl,exp_sin_comp_base*omega)
            exp_sin_comp=Matrix(exp_sin_comp_base).row_insert(sl,-exp_cos_comp_base*omega)
            
        else:
            #print('rayleight damping zero - new code stedy')
            modes = self.modes
            eigs = self.eigenvalues
        
            A = self._fundamental_matrix
            b = amp
    
            exp_cos_comp = -(A.inv() * omega**2 + A).inv()*b
            exp_sin_comp = -omega*A.inv() * exp_cos_comp
          
        return exp_cos_comp*exp(a*self.ivar)*cos(omega*self.ivar) +  exp_sin_comp*exp(a*self.ivar)*sin(omega*self.ivar)
    
    
    def _exp_cos_comp(self,omega,amp,a):

        damping = self._is_proportional_damping
        if damping is not False:
    
            sl=int(N((len(self.dvars))/2))
            #print('rayleight damping - new code stedy')
            modes = self._reduced_modes[:sl,::2]
            modes_trig = self._trig_modes
            eigs = [self._trig_eigenvalues[no,no] for no  in range(len(self.dvars))[::2]]
            
            #display(eigs)
    
            b=Matrix(amp)[sl:,:]
    
            
            base_exp_sin_list=[ (damping*omega*eigen**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            base_exp_cos_list=[ (eigen**2 - omega**2)/((eigen**2 - omega**2)**2 + (damping*omega*eigen**2)**2)  for eigen in eigs]
            
            
            modes_inv=modes.inv()
            modes_x_b=modes*b
            
            exp_cos_comp_base = modes_inv*diag(*base_exp_cos_list)*modes_x_b
            exp_sin_comp_base = modes_inv*diag(*base_exp_sin_list)*modes_x_b
            
            exp_cos_comp=Matrix(exp_cos_comp_base).row_insert(sl,exp_sin_comp_base*omega)
            exp_sin_comp=Matrix(exp_sin_comp_base).row_insert(sl,-exp_cos_comp_base*omega)
            
        else:
            #print('rayleight damping zero - new code stedy')
            modes = self.modes
            eigs = self.eigenvalues
        
            A = self._fundamental_matrix
            b = amp
    
            exp_cos_comp = -(A.inv() * omega**2 + A).inv()*b
            exp_sin_comp = -omega*A.inv() * exp_cos_comp
          
        return exp_cos_comp*exp(a*self.ivar)*cos(omega*self.ivar) +  exp_sin_comp*exp(a*self.ivar)*sin(omega*self.ivar)
    

    
    @cached_property
    def _get_excitation_comps(self):
        '''
        It expands the free terms vector to sequnece of harmonic (sin and cos) components.
        '''
        base_terms = self._free_terms.expand().applyfunc(lambda row: (TR8(row).expand()))
        terms = base_terms

        #base = sin,cos,exp,
        base = Function,

        # components = [
        #     (comp,terms.applyfunc(lambda row: row.coeff(comp))) for comp in terms.atoms(*base) if comp.has(self.ivar)
        # ]
        components = []
        #display(terms)
        
        base_comps_H = terms.atoms(Heaviside)
        
        for comp in base_comps_H:
            # print('_get_excitation_comps -> Heaviside only')
            CodeFlowLogger(comp,'_get_excitation_comps -> Heaviside only',self)
            # print('its called by:')
            # display(self)

            # if type(comp.args[0]) in (cos,sin) and comp.args[0].has(self.ivar): #and not type(comp)==Mul:
            #     #display(terms)
            #     coeff_of_comp = terms.applyfunc(lambda row: row.coeff(comp))
            #     components += [(comp,coeff_of_comp)]
            #     terms  = terms -  coeff_of_comp*comp
            
            if comp.has(self.ivar): #and not type(comp)==Mul:
                #display(terms)
                coeff_of_comp = terms.applyfunc(lambda row: row.coeff(comp))
                components += [(comp,coeff_of_comp)]
                terms  = terms -  coeff_of_comp*comp

        base_comps_sc = terms.atoms(sin,cos)

        #display(terms)
        for comp in base_comps_sc:
            # print('_get_excitation_comps - sin/cos')
            # display(comp,type(comp))
            # print('its called by:')
            # display(self)
            
            if comp.has(self.ivar): #and not type(comp)==Mul:
                #display(terms)
                coeff_of_comp = terms.applyfunc(lambda row: row.coeff(comp))
                components += [(comp,coeff_of_comp)]
                terms  = terms -  coeff_of_comp*comp
  
        base_comps_fun = terms.atoms(Function)

        #display(terms)
        for comp in base_comps_fun:
            # print('_get_excitation_comps - sin/cos')
            # display(comp,type(comp))
            # print('its called by:')
            # display(self)
            
            if comp.has(self.ivar): #and not type(comp)==Mul:
                #display(terms)
                coeff_of_comp = terms.applyfunc(lambda row: row.coeff(comp))
                components += [(comp,coeff_of_comp)]
                terms  = terms -  coeff_of_comp*comp
        
        
        rest = (base_terms - sum([coeff*comp  for comp,coeff in components],Matrix([0]*len(self.dvars)))).doit()
        
        #self._profcheck = True
        if self._profcheck:
            print('Terms collected')
            display(rest)
        

        # #        display('ext_forces',ext_forces)

        # steady_sol = Matrix([0 for gen_coord in self.dvars])

        # for comp in components:
        #     #display(comp)

        #     omg = (comp.args[0].diff(self.ivar)).doit()
        #     #            display(omg)
        #     amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

        #     #display(amp_vector)
        #display(rest)
        return components+[(S.One,rest)]


                                 
    @cached_property
    def _steady_solution(self):
        '''
        It applies generic form solution for the following differential equation
        \dot Y + A Y = F \cos(\Omega t)
        
        The generic form is:
        
        C = (A^{-1} \Omega^2 + A)^{-1}  F 
        D =  \Omega A^{-1} * C
        '''
        

        A = self._fundamental_matrix
        b = self._free_terms
        
        sol = 0*b
        
        #print([elem for elem,coeff in self._get_excitation_comps]) 
        
        for elem,coeff in self._get_excitation_comps:
            if type(elem) == cos:
                omg = (elem.args[0].diff(self.ivar)).doit()
                sol += self._cos_comp(omg,coeff) + self._sin_comp(omg,0*b)
            elif type(elem) == sin:
                omg = (elem.args[0].diff(self.ivar)).doit()
                sol += self._cos_comp(omg,0*b) + self._sin_comp(omg,coeff)

            elif type(elem) == exp:
                a = (elem.args[0].diff(self.ivar)).doit()
                if coeff[1,0].has(sin) == True:
                    for arg in coeff[1,0].args:
                        if type(arg) == sin:
                            omg = arg.args[0].diff(self.ivar).doit()
                    sol += self._exp_cos_comp(omg,0*b,a) + self._exp_sin_comp(omg,coeff,a)
                elif coeff[1,0].has(cos) == True:
                    for arg in coeff[1,0].args:
                        if type(arg) == cos:
                            omg = arg.args[0].diff(self.ivar).doit()
                    sol += self._exp_cos_comp(omg,coeff,a) + self._exp_sin_comp(omg,0*b,a)

            elif type(elem) == Heaviside and type(elem.args[0]) not in (sin,cos):
                
                # print('heaviside check')
                # display(coeff)
                # display(elem)
                ivar0 = -elem.args[0].subs(self.ivar,0)
                #   display(ivar0)
                
                ics_list = list(0*b)
                # print('homo eq')
                odes_hg =self._hom_equation()._as_fode()
                odes = ODESystem(odes_hg.lhs,self.dvars,odes_hg.rhs+coeff,ivar=self.ivar,ode_order=1)
                
                
                # print(f'homo eq of type {type(self._hom_equation())}')
                # display(self._hom_equation())
                
                # print(f'homo eq + coeff as type {type(odes)}')
                # display(odes)
                #sol_H = ODESystem.from_ode_system(odes)._as_fode().solution.with_ics(ics_list,ivar0)
                sol_H = FirstOrderLinearODESystemWithHarmonics.from_ode_system(odes).solution.with_ics(ics_list,ivar0)
                
                
                CodeFlowLogger(sol_H,'Heaviside result')
                CodeFlowLogger(sol_H.rhs*elem)
                
                sol += sol_H.rhs*elem

            elif type(elem) == Heaviside and type(elem.args[0]) in (sin,cos):
                
                # print('heaviside check')
                # display(coeff)
                # display(elem)
                arg = elem.args[0]
                
                if type(arg) == sin:
                    ivar0 =  0
                    gap =(2*pi/ arg.args[0].diff(self.ivar))#.n(6)
    
                elif type(arg) == cos:
                    ivar0 =  -(S.Half*pi/ arg.args[0].diff(self.ivar))#.n(6)
                    gap =(2*pi/ arg.args[0].diff(self.ivar))#.n(6)
                
                ics_list = list(0*b)
                # print('homo eq')
                odes_hg =self._hom_equation()._as_fode()
                odes = ODESystem(odes_hg.lhs,self.dvars,odes_hg.rhs+coeff,ivar=self.ivar,ode_order=1)
                odes_out = ODESystem(odes_hg.lhs,self.dvars,odes_hg.rhs-coeff,ivar=self.ivar,ode_order=1)
                
                # print(f'homo eq of type {type(self._hom_equation())}')
                # display(self._hom_equation())
                
                # print(f'homo eq + coeff as type {type(odes)}')
                # display(odes)
                #sol_H = ODESystem.from_ode_system(odes)._as_fode().solution.with_ics(ics_list,ivar0)
                
                n = Symbol('n')
                CodeFlowLogger(ivar0+gap*n,'ivar value',self)
                
                ivar_H_start = ivar0+gap*n
                ivar_H_end = ivar0+gap/2+gap*n
                
                #ivar0=Symbol('t0')
                sol_H = FirstOrderLinearODESystemWithHarmonics.from_ode_system(odes).solution.with_ics(ics_list,ivar_H_start)
                sol_H_out = FirstOrderLinearODESystemWithHarmonics.from_ode_system(odes_out).solution.with_ics(ics_list,ivar_H_end)
                
                from sympy import Sum,oo
                n_num = Symbol('N')
                
                sol_series_start =  sol_H.rhs.applyfunc(lambda row:  Sum(row*Heaviside(self.ivar - ivar_H_start ),(n,-n_num,n_num,)   ) )
                sol_series_end =  sol_H_out.rhs.applyfunc(lambda row:  Sum(row*Heaviside(self.ivar -  ivar_H_end    ),(n,-n_num,n_num)   ) )
                
                sol_series = sol_series_start + sol_series_end
                
                CodeFlowLogger(sol_H,'Heaviside result')
                CodeFlowLogger(sol_H.rhs*elem)
                
                CodeFlowLogger(sol_series,'SolSeries')
                #CodeFlowLogger(sol_H.rhs*elem)
                
                sol += sol_series#*elem


            elif elem == S.One:
                sol += self._cos_comp(0,coeff) + self._sin_comp(0,coeff)
            
            else:
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                print('element skiped - improve code')
                print(f'called by {type(self)}')
                print(f'elem - base of type {type(elem)}')
                display(elem)
                print('coeff')
            
        ode_sol = ODESolution.from_vars_and_rhs(self.dvars,sol)
        ode_sol.ivar=self.ivar
        
        return ode_sol
    
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

class BernoullisODE(FirstOrderODESystem):
    

    
    # def solution(self,ics=None):
    #     return self.linearized().solution()


        

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
    
    

    def solution(self):
        dvars = self.dvars
        for n in dsolve(self.odes[0]):    
            return AnalyticalSolution(dvars, n.rhs)
        
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
        Solves the problem in the symbolic way and returns matrix of solution (in the form of equations (objects of Eq class)).
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
                
#                 display(eqns_res)

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


# class NumericalODE(ODESystem):

# GOTOWE ale trzeba najpierw ogarnac komponenty raportujace

class SeparableODE(ODESystem):

    def _hom_equation(self):
        ode = super()._hom_equation()
        
        return HomogeneousSeparableODE.from_ode_system(ode)


#     @property
#     def _report_components(self):
        
#         comp_list=[
#         ode_comp_pl.SeparableODEIntroComponent,
#         ode_comp_pl.EquationDefinitionComponent,
#         ode_comp_pl.LinearTransformation,
#         ode_comp_pl.VariablesSeparationComponent,
#         ode_comp_pl.SeparatedVariablesIntegrationComponent,
#         ode_comp_pl.SolutionComponent,
#         ode_comp_pl.VariationOfConstant,
#         ]
        
#         return comp_list

    @property
    def _report_components(self):
        
        comp_list=[
        ode.SeparableODEIntroComponent,
        ode.EquationDefinitionComponent,
        ode.LinearTransformation,
        ode.VariablesSeparationComponent,
        ode.SeparatedVariablesIntegrationComponent,
        ode.SolutionComponent,
        ode.VariationOfConstant,
        ]
        
        return comp_list
    
    def _get_dvars_symbols(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar
        
        dvars_str=f"{str(dvars)}".replace(f'({str(ivar)})','')
#         dvars_str=f"{latex(dvars)[:-17]}"
    
        return {dvars:Symbol(dvars_str),
                ivar:ivar,
                dvars.diff(ivar):Symbol(f"{dvars_str}'")}

    def _variable_name(self, form='expr'):
        system = self
        
        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = f"{latex(fun)[:-17]}'"

        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{latex(fun)[:-17]}')

        y_sym=Symbol('y')
        
        return d_ivar, d_var, y_sym, fun_str
    
    def _function_separation(self,form='expr'):
        
        system = self
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        dvars = system.dvars[0]
        ivar = system.ivar
        dvar_sym= dvars.subs(system._get_dvars_symbols())
        
        
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]

        
        sep_vars_dict = (separatevars(ode_rhs.subs(dvars,dvar_sym),(ivar,dvar_sym),dict=True))

        
        g_fun_ivar = sep_vars_dict[ivar]*sep_vars_dict['coeff']
        h_fun_dvar = sep_vars_dict[dvar_sym]
        
        return g_fun_ivar,h_fun_dvar

    def _function_integration(self,form='expr'):
        system = self
        #y=Function('y')(x)
        
        ode_sys = system
        
        fun = system.dvars[0]
        ivar = system.ivar
        fun_str=f"{latex(fun)[:-17]}'"
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{latex(fun)[:-17]}')
        y_sym=Symbol('y')
        sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))
        g_fun_ivar = sep_vars_dict[ivar]*sep_vars_dict['coeff']
        h_fun_dvar = sep_vars_dict[y_sym]

        first_eq=Eq(Symbol(fun_str),g_fun_ivar*h_fun_dvar * y_sym)

        symbol_1=Symbol(latex(d_var/y_sym)+ '=' + latex(g_fun_ivar) + latex(d_ivar) )
        
        symbol_2=Symbol('\\int ' + latex(d_var/y_sym)+ '=  \\int ' + latex(g_fun_ivar) + latex(d_ivar) )

        symbol_3=Symbol(latex( Integral(1/h_fun_dvar,y_sym) )+ '= ' + latex(Integral(g_fun_ivar,ivar))  )


        const = Symbol('C')
        
        symbol_4=Symbol(latex( Integral(1/h_fun_dvar,y_sym).doit() )+ '= ' + latex(Integral(g_fun_ivar,ivar).doit()  )  + ' + \\ln '  + latex(const)   )
        return{
            'first_eq':first_eq,
            'symbol_1':symbol_1,
            'symbol_2':symbol_2,
            'symbol_3':symbol_1,
            'symbol_4':symbol_4
            
        }

    def _rename_variables(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar
        
        dvars_str=Symbol(f"{str(dvars)}".replace(f'({str(ivar)})',''))
        dvars_g=Function(f"{dvars_str}_o")(ivar)
    
        return {dvars:dvars_g,
                dvars.diff(ivar):Function(f"{dvars_g}'")(ivar)}
    
    def _ode_solution(self):
        dvars=self.dvars
        return ODESolution.from_vars_and_rhs(dvars,dsolve(self.odes[0],self.dvars[0]).rhs)
    def _get_dvars_symbols(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar
        
        dvars_str=f"{str(dvars)}".replace(f'({str(ivar)})','')
#         dvars_str=f"{latex(dvars)[:-17]}"
    
        return {dvars:Symbol(dvars_str),
                ivar:ivar,
                dvars.diff(ivar):Symbol(f"{dvars_str}'")}
    
    @cached_property
    def _general_solution(self):

        return self._ode_solution()
    
class HomogeneousSeparableODE(SeparableODE):

#     @property
#     def _report_components(self):
        
#         comp_list=[
#         ode_comp_pl.SeparableODEIntroComponent,
#         ode_comp_pl.EquationDefinitionComponent,
#         ode_comp_pl.VariablesSeparationComponent,
#         ode_comp_pl.SeparatedVariablesIntegrationComponent,
#         ode_comp_pl.SolutionComponent,
#         ]
        
#         return comp_list

    @property
    def _report_components(self):
        
        comp_list=[
        ode.SeparableODEIntroComponent,
        ode.EquationDefinitionComponent,
        ode.VariablesSeparationComponent,
        ode.SeparatedVariablesIntegrationComponent,
        ode.SolutionComponent,
        ]
        
        return comp_list
    
class LinearWithConstCoeffODE(ODESystem):
#     system = self
#     fun = system.dvars[0]
#     ivar = system.ivar

#     @property
#     def _report_components(self):

#         comp_list=[
#         ode_comp_pl.LinearODEIntroComponent,
#         ode_comp_pl.EquationDefinitionComponent,
#         ode_comp_pl.LinearTransformation,
#         ode_comp_pl.LinearToSeparable,
#         ode_comp_pl.VariationOfConstant
#         ]
        
#         return comp_list
    @property
    def _report_components(self):

        comp_list=[
        ode.LinearODEIntroComponent,
        ode.EquationDefinitionComponent,
        ode.LinearTransformation,
        ode.LinearToSeparable,
        ode.VariationOfConstant
        ]
        
        return comp_list

    def _hom_equation(self):
        ode = super()._hom_equation()
        
        return HomogeneousSeparableODE.from_ode_system(ode)
    
    def _eqn_coeffs(self):
        
        system = self
        ivar = self.ivar
        dvar = system.dvars[0]
        
        coeffs_dict = {dvar.diff(ivar): system.lhs[0].coeff(dvar.diff(ivar)), # coefficient of highest derivative
                       dvar : system.lhs[0].coeff(dvar),
                       'free': system.lhs[0].subs({dvar:0})
                        }
        
        return coeffs_dict
    
    def _variable_name(self, form='expr'):
        system = self
        
        fun = system.dvars[0]
        ivar = system.ivar
        
        fun_str = str(fun).replace(f'({str(ivar)})','' )
        y_prim = f"{fun_str}'"

        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{fun_str}')

        y_sym=Symbol(str(y))
        
        return {'ivar':d_ivar, 'dvar':d_var, 'y_fun':y_sym, 'y':fun_str} ##slownik zrobic z tego
    
    def _function_separation(self,form='expr'):
        
        system = self
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        dvars = system.dvars[0]
        ivar = system.ivar
        dvar_sym= dvars.subs(system._get_dvars_symbols())
        
        
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        
        sep_vars_dict = (separatevars(ode_rhs.subs(dvars,dvar_sym),(ivar,dvar_sym),dict=True))
        
        g_fun_ivar = sep_vars_dict[ivar]*sep_vars_dict['coeff']
        h_fun_dvar = sep_vars_dict[dvar_sym]
        
        return g_fun_ivar,h_fun_dvar

    def _function_integration(self,form='expr'):
        system = self
        
        ode_sys = system
        
        fun = system.dvars[0]
        ivar = system.ivar
        fun_str=self._variable_name()[3]
        
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        
        display(ode_rhs)
        
        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{latex(fun)[:-17]}')
        y_sym=Symbol(str(y))
        sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))
        g_fun_ivar = sep_vars_dict[ivar]*sep_vars_dict['coeff']
        h_fun_dvar = sep_vars_dict[y_sym]

        first_eq=Eq(Symbol(fun_str),g_fun_ivar*h_fun_dvar * y_sym)

        symbol_1=Symbol(latex(d_var/y_sym)+ '=' + latex(g_fun_ivar) + latex(d_ivar) )
        
        symbol_2=Symbol('\\int ' + latex(d_var/y_sym)+ '=  \\int ' + latex(g_fun_ivar) + latex(d_ivar) )

        symbol_3=Symbol(latex( Integral(1/h_fun_dvar,y_sym) )+ '= ' + latex(Integral(g_fun_ivar,ivar))  )

#         const = (system.solution._spot_constant())[0]
        const = Symbol('D')
        
        symbol_4=Symbol(latex( Integral(1/h_fun_dvar,y_sym).doit() )+ '= ' + latex(Integral(g_fun_ivar,ivar).doit()  )  + '+' + latex(const)   )
        return{
            'first_eq':first_eq,
            'symbol_1':symbol_1,
            'symbol_2':symbol_2,
            'symbol_3':symbol_1,
            'symbol_4':symbol_4
            
        }
    def _variation_of_constant(self):
        
        system = self
        ode_sys =system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        ode_sys2 = ode_sys
        x=ivar
        szur=ode_sys2.general_solution
        ode_syso=Eq(ode_sys2.lhs[0],ode_sys2.rhs[0])
        Cx = Function(szur._spot_constant()[0])(x)
        c_prim=Function("C'")(x)
        szur2=szur.subs(szur._spot_constant()[0], Cx)

        return szur2
    
    def _ode_solution(self):
        dvars=self.dvars
        return ODESolution.from_vars_and_rhs(dvars,dsolve(self.odes[0],self.dvars[0]).rhs)
    
    @cached_property
    def _general_solution(self):

        return self._ode_solution()

    def _variable_name(self, form='expr'):
        system = self
        
        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = f"{latex(fun)[:-17]}'"

        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{latex(fun)[:-17]}')

        y_sym=Symbol('y')
        
        return d_ivar, d_var, y_sym, fun_str

    def _get_dvars_symbols(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar
        
        dvars_str=f"{str(dvars)}".replace(f'({str(ivar)})','')
#         dvars_str=f"{latex(dvars)[:-17]}"
    
        return {dvars:Symbol(dvars_str),
                ivar:ivar,
                dvars.diff(ivar):Symbol(f"{dvars_str}'")}

class BernoulliODE(ODESystem):    
#     system = self
#     fun = system.dvars[0]
#     ivar = system.ivar

#     @property
#     def _report_components(self):

#         comp_list=[
#         ode_comp_pl.BernoulliODEIntroComponent,
#         ode_comp_pl.EquationDefinitionComponent,
#         ode_comp_pl.BernoulliTransformation,
#         ode_comp_pl.BernoulliLinearTransformation,
#         ]
        
#         return comp_list

    @property
    def _report_components(self):

        comp_list=[
        ode.BernoulliODEIntroComponent,
        ode.EquationDefinitionComponent,
        ode.BernoulliTransformation,
        ode.BernoulliLinearTransformation,
        ]
        
        return comp_list


    def get_nonlinear_power(self):
        
            
        
            
        power_fun = list(self.odes[0].atoms(Pow))[0]
        power_value = power_fun.atoms(Integer)
        
        
        return power_fun,list(power_value)[0]
    
    @property
    def power_Bernoulli(self):
        return self.get_nonlinear_power()[1]
        
        
    def _bernooullicos(self):
        alpha=Symbol('alpha',positive=True)
        system = self
        ivar = self.ivar
        dvar = system.dvars[0]
        
        
        coeffs_dict = {dvar.diff(ivar): system.lhs[0].coeff(dvar.diff(ivar)), # coefficient of highest derivative
                       dvar : system.lhs[0].coeff(dvar),
                       'free': system.lhs[0].subs({dvar:0})
                        }
        
        return coeffs_dict

    def _bernoulli_substitution(self,subs_fun=None):
        dvars=self.dvars
        ivar=self.ivar
        n_sub=self.power_Bernoulli

        y_pow=self.get_nonlinear_power()[0]
        n=Symbol('n')
        if subs_fun is None:
            z_fun=Function('z')(ivar)
        else:
            z_fun = subs_fun

        z_eq=Eq(z_fun,(dvars**(1-n).subs(n,n_sub))[0])
        z_eq_pow=Eq(z_eq.lhs**n_sub,z_eq.rhs**n_sub)
        z=z_eq.rhs
        y_sol_z=solve(z_eq,dvars[0])[0]

        z_prim=(z.diff(ivar))
        y_prim=y_sol_z.diff(ivar)
        y_pow=y_sol_z**n
        #z_sol=solve(z_eq_pow,y_pow)[

        
        return {
                'z_fun':Matrix([z_fun]),
                'z_diff':((dvars[0])**(S.One-n)).diff(ivar),
                dvars[0] : z_fun**(S.One/1-n)*z_fun**(n*0),
                dvars[0].diff(ivar):y_pow*z_fun.diff(ivar)*z_fun**(n*0),
                n:n_sub,
               }

    def _function_separation(self,form='expr'):
        
        system = self
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        dvars = system.dvars[0]
        ivar = system.ivar
        dvar_sym= dvars.subs(system._get_dvars_symbols())
        
        
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        
        sep_vars_dict = (separatevars(ode_rhs.subs(dvars,dvar_sym),(ivar,dvar_sym),dict=True))
        
        g_fun_ivar = sep_vars_dict[ivar]*sep_vars_dict['coeff']
        h_fun_dvar = sep_vars_dict[dvar_sym]
        
        return g_fun_ivar,h_fun_dvar        
        
    def _variable_name(self, form='expr'):
        system = self

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = f"{latex(fun)[:-17]}'"

        d_ivar = Symbol(f'd{latex(system.ivar)}')
        d_var = Symbol(f'd{latex(fun)[:-17]}')

#        y_sym=Symbol(str(y))

        return d_ivar, d_var, fun_str #y_sym

    def _ode_solution(self):
        dvars=self.dvars
        #return ODESolution.from_vars_and_rhs(dvars,dsolve(self.odes[0],self.dvars[0])[0].rhs)
        return ODESolution(dvars,dvars,dsolve(self.odes[0],self.dvars[0])[0].rhs)

    @cached_property
    def _general_solution(self):

        return self._ode_solution()

    def _get_dvars_symbols(self):
        system=self
        dvars = system.dvars[0]
        ivar = system.ivar
        
        dvars_str=f"{str(dvars)}".replace(f'({str(ivar)})','')
#         dvars_str=f"{latex(dvars)[:-17]}"
    
        return {dvars:Symbol(dvars_str),
                ivar:ivar,
                dvars.diff(ivar):Symbol(f"{dvars_str}'")}
    
    def _bernoulli_subs(self):
        
        #ode = self.subs(self._bernoulli_substitution())
        
        return self._as_reduced_ode()
    
    def _as_reduced_ode(self):
        reduction_dict=self._bernoulli_substitution()

        red_dvar=reduction_dict['z_fun']
        n_sub=self.power_Bernoulli
        
        red_expr=red_dvar[0]**n_sub
        
        ode_rhs=(((self.subs(reduction_dict).lhs*red_expr-self.subs(reduction_dict).rhs*red_expr  ))).expand()

        sep_ode = LinearWithConstCoeffODE(ode_rhs,dvars=red_dvar,ivar=self.ivar,ode_order=1)

        return sep_ode
    
class SystemParameter:
    
    def __init__(self, expr, ivar, parameter_values=None):
        
        
        params = expr
        
        if ivar in params:
            params.remove(ivar)
            
        if parameter_values == None:
            self.system_parameters = list(params)
        else:
            self.system_parameters =  {
                param: parameter_values[no]
                for no, param in enumerate(params)
            }
class HomoODE2ndOrderPL(ODESystem):

    @property
    def _report_components(self):

        comp_list=[
        ode_comp_pl.HomoPredictionIntroComponent,
        ode_comp_pl.EquationDefinitionComponent,
        ode_comp_pl.MainPredictionComponent,
        ode_comp_pl.RootsAnalysisComponent,
        ]
        return comp_list
class HomoODE2ndOrderEN(ODESystem):

    @property
    def _report_components(self):

        comp_list=[
        ode.HomoPredictionIntroComponent,
        ode.EquationDefinitionComponent,
        ode.MainPredictionComponent,
        ode.RootsAnalysisComponent,
        ]
        return comp_list