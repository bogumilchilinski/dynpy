from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   Function, lambdify, factorial, solve, Dict, Number, N, Add,
                   Mul, expand,zoo,exp,Dummy, det)
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import copy
import scipy.integrate as solver
from ..utilities.timeseries import DataMethods, SpectralMethods, TimeDomainMethods, SpectrumSeries, SpectrumFrame, TimeSeries, TimeDataFrame

from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3, TR0

from .linear import LinearODESolution, FirstOrderODE, AnalyticalSolution, FirstOrderLinearODESystem, FirstOrderODESystem, ODESystem, ODESolution,FirstOrderLinearODESystemWithHarmonics
from ..utilities.components.ode import en as ode_comp

#from timer import timer
from time import time
from collections.abc import Iterable
from functools import cached_property,lru_cache

from .tools import CommonFactorDetector, ODE_COMPONENTS_LIST, CodeFlowLogger

def const_no_gen():
    num = 0
    while True:
        yield num
        num += 1

class SimplifiedExpr:
    
    """

    This class provides methods which allow to obtain a differential equation in simplified form, or extract desired simplification parameters. 

    """

    _const_no=const_no_gen()
    
    _subs_container={}
    def __init__(self,expr,ivar,parameters = None ):

        self._expr = expr
        self._ivar = ivar
        self._parameters = parameters

    @property
    def _fun_list(self):
        
        """
        
        Method extracts functions from expression and returns them in a list.
        
        """
        
        funlist=[expr  for  expr in  (self._expr).atoms(sin,cos) if expr.has(self._ivar)]
        
        return funlist
    
    @property
    def simplified_expr(self):
        
        """
        
        This method returns a simplified expression, prepared in sum_expr method.
        
        """
        
        return self.sum_expr
    @property
    def sum_expr(self):
        
        """
        
        This method provides simplified expression that consist of main part of an equation and functions with additional coefficents.
        
        """

        data = self._expr_collector
        expr_dict={key:values[0]  for  key,values  in  data.items() }
        self.__class__._subs_container={**self.__class__._subs_container, **expr_dict}
        
        simplified_sum=sum([key*values[1]  for  key,values  in  data.items()  ],S.Zero)
        expanded_sum=sum([values[0]*values[1]  for  key,values  in  data.items()  ],S.Zero)
        
        expr_rest = (self._expr - expanded_sum).expand()
        # display(expr_rest)
        
        # print(len(self._expr.args))
        # print(len(simplified_sum.args))
        return simplified_sum + expr_rest

    
    @property
    def _expr_dict(self):
        
        """

        Method returns a dictionary that contains additional coefficients with their values.

        """
        
        expr_dict={key:values[0]  for  key,values  in  self._expr_collector.items() }
        return expr_dict
    @property
    def _expr_collector(self):
        
        """

        This method combines functions and coefficients and then returns a dictionary with segregated values.

        """
        
        parameters = self._parameters
        
        if parameters:
            aux_symbol = lambda no: Function(f'DummyFun_{next(self._const_no)}')(*parameters)
        else:
            aux_symbol = lambda no: Dummy(f'D_{no}')
        
        bucket = {aux_symbol(no):(self._expr.coeff(fun),fun) for no,fun in enumerate(self._fun_list)}
        return bucket
    
    @property
    def full_expr(self):
        
        """

        Main method that returns a full expression with functions and proper coefficients.

        """
        
        #display(self.__class__._subs_container)
    
        return self._expr.doit().expand().subs(self.__class__._subs_container)
    
    def __repr__(self,*args):
        return 'instance of SimplifiedExpr'

class PerturbationODESolution(ODESolution):
    _eps = Symbol('varepsilon')
    
    def set_small_parameter(self,eps,inplace = False):
        obj = copy.copy(self)
        obj._eps = eps
        
        return obj
    
    @property
    def small_parameter(self):
        return self._eps

    @small_parameter.setter
    def small_parameter(self,eps):
        
        self._eps = eps
        
        return self._eps
    
    @property
    def eps(self):
        return self.small_parameter
    
    def _calculate_constant(self,ics=None,*args,**kwargs):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """

        const_eqns=self._get_constant_eqns(ics).subs(self.eps,0).doit()
        const_list = self._spot_constant()
        
        #display('ivar',self.ivar)
        
        return solve(const_eqns,const_list)
    

class ODESystemApproximation(ODESystem):

    _simp_dict = {}
    _callback_dict = {}

    @property
    def _default_solver(self):
        
        return FirstOrderLinearODESystemWithHarmonics
        #return FirstOrderLinearODESystem
    
class NthOrderODEsApproximation(ODESystem):
    _ode_order=2

    _simp_dict = {}
    _callback_dict = {}
    
    # _const_list = []
    
    @property
    def _default_solver(self):
        
        return FirstOrderLinearODESystemWithHarmonics
        #return FirstOrderLinearODESystem
    
    @property
    def _secular_funcs(self):

        # weird error occurs, when .remove_secular_terms is called. if secular terms are missing, method returns None,
        # if statement should be considered
        
        #const_to_add=[eqn for eqn in self._const_list ]
        
        if self._parameters is None:
            params = []
        else:
            params = self._parameters
        
        secular_funcs = set() | (self.general_solution.rhs.atoms(Function) - set(
            self._const_list) - set(params) - {self.ivar})
        #display(secular_funcs)
        
        return {sec_fun for sec_fun in secular_funcs if sec_fun.has(self.ivar)}

    @property
    def secular_terms(self):

        #print('secsec')
        sec_funcs = self._secular_funcs
        #self._const_list


        CodeFlowLogger(sec_funcs,'sec_funcs',self)
        
        sec_conditions = Matrix(list(set(sum([list(self.odes.applyfunc(lambda entry: entry.coeff(func))) for func in sec_funcs],[]))-{0}))


        ivar = self._parameters[0]

        
        #sec_conditions =FirstOrderODESystem.from_ode_system(ODESystem([sec_conditions[1][1],sec_conditions[3][1]],dvars = self._const_list ,ivar=ivar  )).linearized().solution
        
        # print('some check fof sec eq')
        

        CodeFlowLogger(sec_conditions,'sec_conditions',self)
        
        const_list = [fun_cons for fun_cons in list(sec_conditions.atoms(Function)) if  'C' in str(fun_cons) ]
        

        CodeFlowLogger(const_list,'const_list',self)
        
        sec_odes=ODESystem((sec_conditions),
                            dvars=Matrix(const_list),
                            ivar=ivar,  
                            ode_order=None)
        

        CodeFlowLogger(sec_odes.as_type(FirstOrderLinearODESystem).solution,'sec_odes.as_type(FirstOrderLinearODESystem).solution',self)
        #fode_harm_rhs = FirstOrderLinearODESystemWithHarmonics.from_ode_system(sec_odes)
        
        # print('transformed')
        # display(fode_harm_rhs)
        # display(type(fode_harm_rhs))
        
        return sec_odes

    def remove_secular_terms(self):


        obj = ODESystem.from_ode_system(self)

        # subs_dict = {comp: 0 for comp in self._secular_funcs}
        
        # return obj.subs(subs_dict).doit().subs(zoo,0).doit()
        
        def eqns_map(obj):

            oper_expr = lambda obj: (TR10(
                TR8(TR10(TR0(obj.expand())).expand())).expand().doit().expand())

            #oper_expr = lambda obj: TR8(obj).expand()

            #oper_expr = lambda obj: obj.doit().expand()

            if isinstance(obj, Add):
                elems = obj.args
                return sum(oper_expr(elem) for elem in elems)
            else:
                return oper_expr(obj)
            

        obj = obj.applyfunc(eqns_map)
        
        secular_expr_list = [obj.applyfunc(lambda row: row.coeff(comp)*comp).rhs  for comp in self._secular_funcs]
        secular_expr = sum(secular_expr_list,Matrix( [0]*len(obj) ))

        return ODESystem((obj.as_matrix() - secular_expr).expand().subs({comp:0 for comp in self._secular_funcs}).doit(),obj.dvars,ivar=obj.ivar)



        
#     def find_approximation(self):

#         zeroth_ord = self.odes.applyfunc(lambda ode: ode.general_solution.self.remove_secular_terms())

#         for order,approx in enumerate(self.odes):

#             approx._parameters = self.t_list[order+1:]

#             nth_ord = approx.applyfunc(lambda ode: ode.steady_solution.self.remove_secular_terms())

#         return {order:}




class MultiTimeScaleSolution(ODESystem):
    _ode_order=2
    _stored_solution = None
    _saved_solution = {}

    _saved_numerized = {}
    _order = 2
    _scales_no = None

    def __new__(cls,
                odes_system,
                dvars,
                ivar=None,
                ics=None,
                eps=Symbol('varepsilon'),
                omega=None,
                order=2,
                scales_no=None,
                t_span=[],
                params=[],
                params_values={},
                extra_params={},
                ic_point={},
                equation_type=None,
                label=None,
                evaluate=True,
                **options):

        if not isinstance(odes_system, Iterable):
            odes_system = Matrix([odes_system])
        if not isinstance(dvars, Iterable):
            dvars = Matrix([dvars])

        obj = super().__new__(cls,
                              odes_system,
                              dvars=dvars,
                              ivar=ivar,
                              evaluate=evaluate,
                              **options)

        obj._const_list = []

        if ivar is not None:
            obj._ivar = ivar

        if label == None:
            label = obj._label = obj.__class__.__name__ + ' with ' + str(
                len(obj.dvars)) + ' equations'
            
        obj._label = 123
                
        obj.eps = eps
        obj._order = order
        
        if scales_no is None:
            obj._scales_no = obj.order+1

        obj._stored_solution = None
        obj._saved_solution = {}
        obj.extra_params = extra_params

        obj.secular_eq = {}
        obj.int_const = set()

        obj.omega = 1
        if omega:
            obj.omega = omega
        
        obj.ics = ics

        obj._const_sol = {}

        
        obj._eoms = odes_system
        obj.odes_system = odes_system
        #obj._odes = odes_system

        return obj
    
    @property
    def eps(self):
        return self._eps
    
    @eps.setter
    def eps(self,eps):
        self._eps=eps
        
    @property
    def order(self):
        return self._order
    
    @property
    def scales_no(self):
        if self._scales_no is None:
            return self.order+1
        else:
            return self._scales_no
    
    @scales_no.setter
    def scales_no(self,no):
        
        self._scales_no = no
    
    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, omega):
        self._omega = omega
        
    @property
    def ics(self):
        return self._ics

    @ics.setter
    def ics(self, ics):
        self._ics = ics
    @property
    def extra_params(self):
        return self._extra_params

    @ics.setter
    def extra_params(self, extra_params):
        self._ics = extra_params
        
    def __str__(self):
        return '123'

    def __repr__(self):
        return '123'

    def set_solution_order(self, order=1):

        self._order = order

    @property
    def _report_components(self):
        
        comp_list=[
        ode_comp.ODESystemComponent,
        ode_comp.VariablesComponent,

        ode_comp.ODESystemCodeComponent,
        ode_comp.MSMCalculationsOrderComponent,        
        ode_comp.PredictedSolutionComponent,
        ode_comp.DetailsOfPredictedSolutionComponent,

        ode_comp.GoverningEquationComponent,

        #ode_comp.ApproximatedNonlinearGoverningEquationComponent,
#         ode_comp.FundamentalMatrixComponent,
        ode_comp.ZerothOrderApproximatedEqComponent,
        ode_comp.FirstOrderApproximatedEqComponent,

        # ode_comp.GeneralSolutionComponent,
        ode_comp.ZerothOrderSolutionComponent,
        ode_comp.SecularTermsEquationsComponent,
            
        ]
        
        return comp_list


    
    @property
    def order(self):
        return self._order

    @order.setter
    def order(self,order):
        self._order = order

    def set_order(self,order):
        
        
        new_msm = MultiTimeScaleSolution(self.as_matrix(),self.vars,ivar=self.ivar,omega=S.One,order=order,scales_no=self.scales_no,eps=self.eps )
        new_msm._scales_no = self._scales_no
        
        return new_msm
    
    def _assign_properties(self,obj):
        
        obj._eps=self.eps     
        obj._vars = self._vars
        #obj._const_list = self._const_list
        
        return obj
    
    
    @property
    def ics_dvars(self):

        q = list(self.dvars)
        dq = list(Matrix(self.dvars).diff(self.ivar))

        return Matrix(q + dq)

#     @order.setter
#     def order(self, order):
#         self._order = order

    @cached_property
    def t_list(self):
        r"""
        Returns the list of time domanin slow varying function
        """

        self._t_list = [
            t_i(self.ivar)
            for t_i in symbols(f't_0:{self.scales_no}', cls=Function, real=True)
        ]

        return self._t_list
    
    @cached_property
    def _scales_formula(self):
        r"""
        Returns the list of time domanin slow varying function
        """

        t_list = self.t_list
        expr_list = [self.ivar*self.eps**no for no in range(len(t_list))]

        return AnalyticalSolution.from_vars_and_rhs(t_list, expr_list)
    
    
    def approx_function_symbol(self,order=0):
      
        dvars_no = len(self.dvars)
        
        if dvars_no==1:
        
            func_str = str(self.vars[0]).replace(f"({self.ivar})","")
        else:
            func_str = 'Y'
        
        return func_str + f"_{{{str(order)}}}"

    def approximation_function(self, order, max_order=1, ivar_list = None ):

        dvars_no = len(self.vars)

        if ivar_list is None:
            t_list = self.t_list
        else:
            t_list = ivar_list
        
        fun_str = self.approx_function_symbol(order)
        
        if dvars_no==1:

            aprrox_fun = [
                sym.Function(fun_str)(*t_list)
                for dvar_no, dvar in enumerate(self.vars)
            ]
        else:


            aprrox_fun = [
                sym.Function(fun_str + '_{' + str(order) + '}')(*t_list)
                for dvar_no, dvar in enumerate(self.vars)
            ]

        return Matrix(aprrox_fun)

    def approximation_function_without_scales(self, order, max_order=1):

        return self.approximation_function(order=order,max_order=max_order,ivar_list=[self.ivar])
    
    def part_derivative(self, order=1):
        return Eq(self.approximation_function_without_scales(order).diff(self.ivar)[0],self.approximation_function(order).diff(self.ivar)[0])
    
    def part_derivative_subs(self, order=1):
        deriv_dict = {}

        for index, value in enumerate(self._scales_formula):
            scales_lhs = self._scales_formula.lhs[index]
            scales_rhs = self._scales_formula.rhs[index]

            scales_lhs_deriv = self.first_order_subs().lhs[index]
            scales_rhs_deriv = self.first_order_subs().rhs[index]

            scales_lhs_second_deriv = self.second_order_subs().lhs[index]
            scales_rhs_second_deriv = self.second_order_subs().lhs[index]

            deriv_dict[scales_lhs] = scales_rhs
            deriv_dict[scales_lhs_deriv] = scales_rhs_deriv
            deriv_dict[scales_lhs_second_deriv] = scales_rhs_second_deriv

        return Eq(self.approximation_function_without_scales(order).diff(self.ivar)[0],self.part_derivative(order).rhs.subs(deriv_dict).simplify().expand())

    
    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.vars)
        order = self.order

        solution = sum(
            (self.approximation_function(comp_ord, order) * self.eps**comp_ord
            for comp_ord in range(order + 1)), sym.zeros(dvars_no, 1))

        return ODESolution.from_vars_and_rhs(self.vars, solution)

    def predicted_solution_without_scales(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)

        solution = sum(
             (self.approximation_function_without_scales(comp_ord, order) * self.eps**comp_ord
             for comp_ord in range(order + 1)), sym.zeros(dvars_no, 1))

        return ODESolution.from_vars_and_rhs(self.dvars, solution)
    
    def derivative_without_scale(self, degree, order=1, dict=False, equation=False):
        
        derivative_rhs = self.predicted_solution_without_scales(self.order).as_eq_list()[0].n(3).rhs.diff(self.ivar,degree)
        derivative_lhs = self.dvars[0].diff(self.ivar,degree)
        
        return AnalyticalSolution(derivative_lhs, rhs = derivative_rhs, vars = self.dvars[0])

    def first_order_subs(self):

        first_ord_subs = {
            t_i.diff(self.ivar): self.eps**t_ord
            for t_ord, t_i in enumerate(self.t_list)
        }

        return AnalyticalSolution.from_dict(first_ord_subs)
    
    def second_order_subs(self):

        sec_ord_subs = {
            t_i.diff(self.ivar, 2): 0
            for t_ord, t_i in enumerate(self.t_list)
        }

        return AnalyticalSolution.from_dict(sec_ord_subs)
    
    #@lru_cache
    def eoms_approximation(self, order=None, odes_system=None):

        
        order = self.order
        
        if not odes_system:
            odes_system = self.as_matrix()

        first_ord_subs = self.first_order_subs().as_explicit_dict()
        sec_ord_subs = self.second_order_subs().as_explicit_dict()

        #display(self.predicted_solution(order).as_dict().subs(sec_ord_subs).doit().subs(first_ord_subs).doit())

        
        odes_rhs = self.as_matrix().subs(self.predicted_solution(
            order).as_explicit_dict()).subs(sec_ord_subs).doit().subs(first_ord_subs).doit()

        #display()
        #eoms_approximated = ODESystem(odes_system,dvars=self.dvars,parameters = self.t_list).subs(
        #    self.predicted_solution(order).as_dict())
        eoms_approximated = - odes_rhs

        eoms_approximated = eoms_approximated.subs(sec_ord_subs).doit().subs(
            first_ord_subs).doit()

        t_fun_subs_dict = ({
            expr_tmp: expr_tmp.subs(self.ivar, self.t_list[0])
            for expr_tmp in (
                eoms_approximated.atoms(Function) - {*self.t_list} - {
                    *sum([
                        list(self.approximation_function(ord_tmp, order))
                        for ord_tmp in range(order + 1)
                    ], [])
                })
        })

        CodeFlowLogger(eoms_approximated,'eoms for approximation before subs',self)
        CodeFlowLogger(t_fun_subs_dict,'dict with time subs',self)

        return eoms_approximated.subs(t_fun_subs_dict).doit()

    #@lru_cache
    def eoms_approximation_list(self,
                                max_order=3,
                                min_order=0,
                                odes_system=None):

        eoms_approximated = self.eoms_approximation(
            order=max_order, odes_system=odes_system).expand()

        order_get = lambda obj, order: (obj.diff(self.eps, order) / factorial(
            order)).subs(self.eps, 0).doit()

        #NthOrderODEsApproximation
        
        
        approx_list=[
            NthOrderODEsApproximation(
                eoms_approximated.applyfunc(lambda obj: order_get(obj, order)),
                dvars=self.approximation_function(order),
                ivar=self.t_list[0],
                ode_order=2,
                parameters=self.t_list)
            for order in range(min_order, max_order + 1)
        ]
        
        
        return approx_list
        #return [NthOrderODEsApproximation.from_ode_system(approx_ode)   for approx_ode in approx_list]

    #@lru_cache
    def eoms_approx_with_const(self, order=1):

        self.secular_eq = {}

        approx_eoms_list = self.eoms_approximation_list(self.order)

        approx_eoms_list[0]._parameters = self._t_list[1:]
        
        # display(approx_eoms_list)
        # display(approx_eoms_list[0])
        
        sol = approx_eoms_list[0].solution
        # print('spot_const')
        # display(sol._spot_constant())
        
        #### !!!! 
        sol_list = [sol]
        
        
        sol_subs_list = [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit() for sol in sol_list]
        sol_subs_dict  = {  dvar:eqn     for  sol   in sol_subs_list  for dvar, eqn in sol.as_explicit_dict().items()}

        
        # print('_gen_sol')
        # display(sol_list)
        # display(sol_subs_dict)
        
        # display(approx_eoms_list[0]._const_list)
        approx_with_const = [approx_eoms_list[0]]

        for order, approx in enumerate(approx_eoms_list[1:]):

            approx._parameters = self._t_list[1:]
            #approx = approx.as_first_ode_linear_system()

            eqns_map = lambda obj: (TR10(TR8(TR10(obj.expand()).expand())).
                                    expand().doit().expand())

            def eqns_map(obj):

                oper_expr = lambda obj: (TR10(
                    TR8(TR10(TR0(obj.expand())).expand())).expand().doit().expand())

                #oper_expr = lambda obj: TR8(obj).expand()

                #oper_expr = lambda obj: obj.doit().expand()

                if isinstance(obj, Add):
                    elems = obj.args
                    return sum(oper_expr(elem) for elem in elems)
                else:
                    return oper_expr(obj)
            
            # print('approx eqn')
            # display(approx)
            
            sol_subs_list = [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit() for sol in sol_list]
            sol_subs_dict  = {  dvar:eqn     for  sol   in sol_subs_list  for dvar, eqn in sol.as_explicit_dict().items()}
            
            sol_subs_dict  = sol.as_explicit_dict()
        
            approx_subs = approx.applyfunc(eqns_map).subs(
                sol_subs_dict).applyfunc(eqns_map)
            
            # display(approx_subs)

            approx_subs = NthOrderODEsApproximation.from_ode_system(approx_subs)
            approx_subs._parameters = self.t_list[1:]
            
            approx_with_const += [approx_subs]
            
        return approx_with_const


    def _general_sol(self, order=1):

        self.secular_eq = {}

        
        approx_with_const = self.eoms_approx_with_const(order)

        sol = approx_with_const[0].solution
        sol_list = [sol]

        for order, approx in enumerate(approx_with_const[1:]):

            approx_subs = approx
            self.secular_eq[self.eps**(order +1)] = (approx_subs.secular_terms.as_type(FirstOrderLinearODESystem))
            
            #nonlin_ode.eoms_approx_with_const()[0].remove_secular_terms().steady_solution.rhs
            approx_subs = approx_subs#.steady_solution.rhs
            
            #aux_mat=approx_subs.jacobian(approx_subs.dvars)
            # if det(aux_mat)==0:

            #     SolverClass = FirstOrderLinearODESystem
            # else:
            #     SolverClass = FirstOrderLinearODESystemWithHarmonics

            # SolverClass = FirstOrderLinearODESystemWithHarmonics
            # #SolverClass = FirstOrderLinearODESystem
            
            if len(sol) > 0:
                ode_2_solve = approx_subs.subs(sol_list[-1].as_explicit_dict()).remove_secular_terms()
            else:
                ode_2_solve = approx_subs.remove_secular_terms()
            
            # ode_2_solve._ivar = self.t_list[0]
            
            # ode_2_solve = SolverClass.from_ode_system(ode_2_solve)

            

            sol = ode_2_solve.steady_solution.applyfunc(
                lambda obj: obj.expand())#.applyfunc(eqns_map)#.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).sum_expr)


            sol_list += [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit()]
            #sol_list += [sol]

        CodeFlowLogger(sol_list,'sol list',self)

        return sol_list

    @property
    def general_solution(self):

        
        order = self.order
        
        #print('my order',order)
        sol_list = []
        gen_sol = self._general_sol(order)
        
        # print('gen_sol')
        # display(*gen_sol)
        
        if self.eps in self.secular_eq:

            ode2check = self.secular_eq[self.eps]#.as_first_ode_linear_system()

            
            CodeFlowLogger(ode2check,'sec_eq',self)
            
            # display(ode2check._as_fode())
            # display(ode2check.solution)
            #ode2check._as_fode()
            #ode2check.solution
            C_const_sol = ode2check.solution.as_explicit_dict()
            

            CodeFlowLogger(C_const_sol,'tut C_const_sol',self)
            
        else:
            C_const_sol={}
        # print('tut')
        for order, approx_sol in enumerate(gen_sol):
            
            # print('stale: ')
            # display(C_const_sol)
            # print('aproox bez subsa')
            
            
            # display(approx_sol.rhs.applyfunc(lambda obj: obj.expand().doit().subs(C_const_sol).subs(C_const_sol)
            #                                ))
            solution = approx_sol.applyfunc(lambda obj: obj.expand().doit().subs(C_const_sol).subs(C_const_sol)).subs(
                {self.t_list[1]: self.eps  * self.ivar})

            # print('aproox ze subsa')
            # display(solution)

            sol_list += [(self.eps**order) *
                            solution.rhs.subs({
                                self.t_list[1]: self.eps  * self.ivar,
                                self.t_list[0]: self.ivar
                            })   ]
            
        #display(*list(SimplifiedExpr._subs_container.values()))
        result = (sum(sol_list, Matrix(2*len(self.dvars)*[0])  )).applyfunc(lambda obj: obj.expand().doit())
        
        new_res = PerturbationODESolution.from_vars_and_rhs(Matrix(list(self.dvars) +  list(self.dvars.diff(self.ivar)) ) , result)#.set_small_parameter(self.eps)
        new_res.small_parameter = self.eps
        new_res.ivar = self.ivar

        return new_res

    @cached_property
    def steady_solution(self):

        result = Matrix(2*len(self.dvars)*[0])

        new_res = PerturbationODESolution(Matrix(list(self.dvars) +  list(self.dvars.diff(self.ivar)) ) , result)#.set_small_parameter(self.eps)
        new_res.small_parameter = self.eps
        new_res.ivar = self.ivar
        return new_res

    @cached_property
    def solution(self):

        return self.general_solution


    def compute_solution(   
                        self,
                        t_span=None,
                        ic_list=None,
                        t_eval=None,
                        params_values=None,
                        method='RK45'
                        ):

        #         print('compute_solution')
        #         display(params_values)
        #         display(ic_list)

        if ic_list:
            print('class ics has been taken')
            self.ics = ic_list

        return self.__call__(ivar=t_span,
                             order=self._order,
                             params_values=params_values)


