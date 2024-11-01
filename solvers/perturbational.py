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
from functools import cached_property

from .tools import CommonFactorDetector, ODE_COMPONENTS_LIST, CodeFlowLogger

def const_no_gen():
    num = 0
    while True:
        yield num
        num += 1

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

    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.vars)

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
    

    def eoms_approximation(self, order=None, odes_system=None):

        if order is None:
            order = self.order
        
        if not odes_system:
            odes_system = self.as_matrix()

        first_ord_subs = self.first_order_subs()
#         {
#             t_i.diff(self.ivar): self.eps**t_ord
#             for t_ord, t_i in enumerate(self.t_list)
#         }
        sec_ord_subs = self.second_order_subs()
#         {
#             t_i.diff(self.ivar, 2): 0
#             for t_ord, t_i in enumerate(self.t_list)
#         }

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


    def eoms_approx_with_const(self, order=1):

        self.secular_eq = {}

        approx_eoms_list = self.eoms_approximation_list(order)

        approx_eoms_list[0]._parameters = self._t_list[1:]
        
        # display(approx_eoms_list)
        # display(approx_eoms_list[0])
        
        sol = approx_eoms_list[0].as_type(FirstOrderLinearODESystem).solution
        # print('spot_const')
        # display(sol._spot_constant())
        
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

                oper_expr = lambda obj: obj.doit().expand()

                if isinstance(obj, Add):
                    elems = obj.args
                    return sum(oper_expr(elem) for elem in elems)
                else:
                    return oper_expr(obj)
            
            # print('approx eqn')
            # display(approx)
            
            sol_subs_list = [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit() for sol in sol_list]
            sol_subs_dict  = {  dvar:eqn     for  sol   in sol_subs_list  for dvar, eqn in sol.as_explicit_dict().items()}
        
            approx_subs = approx.applyfunc(eqns_map).subs(
                sol_subs_dict).applyfunc(eqns_map)
            
            # display(approx_subs)

            approx_subs = NthOrderODEsApproximation.from_ode_system(approx_subs)
            approx_subs._parameters = self._t_list[1:]
            
            approx_with_const += [approx_subs]
            
        return approx_with_const


    def _general_sol(self, order=1):

        self.secular_eq = {}

        approx_eoms_list = self.eoms_approximation_list(order)

        approx_eoms_list[0]._parameters = self._t_list[1:]
        
        # display(approx_eoms_list)
        # display(approx_eoms_list[0])
        
        sol = approx_eoms_list[0].as_type(FirstOrderLinearODESystem).solution
        # print('spot_const')
        # display(sol._spot_constant())
        
        sol_list = [sol]
        sol_subs_list = [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit() for sol in sol_list]
        sol_subs_dict  = {  dvar:eqn     for  sol   in sol_subs_list  for dvar, eqn in sol.as_explicit_dict().items()}

        
        # print('_gen_sol')
        # display(sol_list)
        # display(sol_subs_dict)
        
        # display(approx_eoms_list[0]._const_list)

        for order, approx in enumerate(approx_eoms_list[1:]):

            approx._parameters = self._t_list[1:]
            #approx = approx.as_first_ode_linear_system()

            eqns_map = lambda obj: (TR10(TR8(TR10(obj.expand()).expand())).
                                    expand().doit().expand())

            def eqns_map(obj):

                oper_expr = lambda obj: (TR10(
                    TR8(TR10(TR0(obj.expand())).expand())).expand().doit().expand())

                oper_expr = lambda obj: obj.doit().expand()

                if isinstance(obj, Add):
                    elems = obj.args
                    return sum(oper_expr(elem) for elem in elems)
                else:
                    return oper_expr(obj)
            
            # print('approx eqn')
            # display(approx)
            
            sol_subs_list = [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit() for sol in sol_list]
            sol_subs_dict  = {  dvar:eqn     for  sol   in sol_subs_list  for dvar, eqn in sol.as_explicit_dict().items()}
        
            approx_subs = approx.applyfunc(eqns_map).subs(
                sol_subs_dict).applyfunc(eqns_map)
            
            # display(approx_subs)

            approx_subs = NthOrderODEsApproximation.from_ode_system(approx_subs)
            approx_subs._parameters = self._t_list[1:]
            
            # display(approx_subs)
            # display(approx_subs.ivar)
            # display(approx_subs.secular_terms)
            # display(type(approx_subs.secular_terms))

            self.secular_eq[self.eps**(order +
                                       1)] = (approx_subs.secular_terms.as_type(FirstOrderLinearODESystem))
            
            
            #print('with secular terms')
            #display(approx_subs)
            approx_subs = approx_subs.remove_secular_terms()
            
            #print('without secular terms')
            #display(approx_subs)

            #display(approx_subs.lhs,approx_subs.rhs)
            #display(approx_subs)
            #display(FirstOrderLinearODESystemWithHarmonics.from_ode_system(approx_subs).steady_solution)
            #display(type(approx_subs))
            
            #aux_mat=self.odes_system.jacobian(self.dvars)
            aux_mat=approx_subs.jacobian(approx_subs.dvars)
            if det(aux_mat)==0:

                SolverClass = FirstOrderLinearODESystem
            else:
                SolverClass = FirstOrderLinearODESystemWithHarmonics

            #SolverClass = FirstOrderLinearODESystemWithHarmonics
            #SolverClass = FirstOrderLinearODESystem
            
            
            ode_2_solve = ODESystem.from_ode_system(approx_subs)
            
            ode_2_solve._ivar = self._t_list[0]
            
            ode_2_solve = SolverClass.from_ode_system(ode_2_solve)
            # print('gen sol part')
            #display(*FirstOrderLinearODESystemWithHarmonics.from_ode_system(ode_2_solve)._get_excitation_comps)
            
            #fm_mat = ode_2_solve._as_fode()._fundamental_matrix
            #display(fm_mat)
            #display(fm_mat.diagonalize())
            

            
            sol = ode_2_solve.steady_solution.applyfunc(
                lambda obj: obj.expand()).applyfunc(eqns_map)#.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).sum_expr)

            # print('and its sol')
            # display(sol)
            #sol_subs_dict = {**sol_subs_dict, **sol.as_dict()}
            #display(*list(sol.as_matrix()))
            sol_list += [sol.applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).full_expr).doit()]
            #sol_list += [sol]
        print('sol list')
        display(*sol_list)

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
            print('sec_eq')
            ode2check = self.secular_eq[self.eps]#.as_first_ode_linear_system()
            print(type(ode2check))
            display(ode2check)
            # display(ode2check._as_fode())
            # display(ode2check.solution)
            #ode2check._as_fode()
            #ode2check.solution
            C_const_sol = ode2check.solution.as_explicit_dict()
            
            print('tut')
            display(C_const_sol)
            
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
                            solution.subs({
                                self.t_list[1]: self.eps  * self.ivar,
                                self.t_list[0]: self.ivar
                            }).rhs   ]
            
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




class WeakNonlinearProblemSolution(LinearODESolution):

    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 eps=Symbol('varepsilon'),
                 omega=None,
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 equation_type=None):

        super().__init__(odes_system=odes_system,
                         ivar=ivar,
                         dvars=dvars,
                         t_span=t_span,
                         params=params,
                         params_values=params_values,
                         ic_point=ic_point,
                         equation_type=equation_type)
        self.eps = eps

        #        display(omega)
        self.omega = 1
        if omega:
            self.omega = omega

    def __call__(self, time, order=1, params_values=None):

        if params_values:
            self.params_values = params_values

        if isinstance(time, Symbol):
            return self.nth_order_solution(order).subs(self.ivar, time).subs(
                self.params_values)
        else:
            nth_order_solution_fun = lambdify(
                self.ivar,
                self.nth_order_solution(order).subs(self.params_values),
                'numpy')

            #            return nth_order_solution_fun(time)
            result = TimeDataFrame(index=time,
                                   data={
                                       'solution':
                                       (nth_order_solution_fun(time))
                                   }).copy().apply(np.real)

            display(result)
            display(result.apply(np.real))

            return result

    def _format_solution(self, dvars, solution, dict=False, equation=False):

        if equation:
            solution = Matrix(
                [Eq(lhs, rhs) for lhs, rhs in zip(dvars, list(solution))])

        if dict:
            solution = {lhs: rhs for lhs, rhs in zip(dvars, list(solution))}

        return solution

    def approximation_function(self, order, max_order=1):

        dvars_no = len(self.dvars)

        t_list = self.t_list

        aprrox_fun = [
            sym.Function('Y_{' + str(dvar_no) + str(order) + '}')(*t_list)
            for dvar_no, dvar in enumerate(self.dvars)
        ]

        return Matrix(aprrox_fun)

    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)

        solution = sum(
            (self.approximation_function(comp_ord) * self.eps**comp_ord
             for comp_ord in range(order + 1)), sym.zeros(dvars_no, 1))

        return self._format_solution(dvars=self.dvars,
                                     solution=solution,
                                     dict=dict,
                                     equation=equation)

    def eigenvalues(self):
        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        return [eig for eig in eigs if eig != 0]

    def eigenvalue_approximation(self,
                                 order=1,
                                 eigenvalue=None,
                                 eigenvalue_no=None,
                                 eigenvalue_symbol=Symbol('omega'),
                                 series_coefficients=None):

        if not eigenvalue:
            eigenvalue = self.omega
        self.omega = eigenvalue

        if not series_coefficients:
            series_coefficients = [
                symbols('A_{}{}_1:{}'.format(eom_no, no, order + 1))
                for eom_no, eom in enumerate(self.dvars)
                for no, coord in enumerate(self.dvars)
            ]

        #approximation={eigenvalue:eigenvalue_symbol**2 - sum(coeff**(eps_ord+1)  for eps_ord,coeff  in enumerate(series_coefficients) )}
        approximations = [
            eigenvalue**2 - sum([
                coeff * self.eps**(eps_ord + 1)
                for eps_ord, coeff in enumerate(coeffs)
            ]) for coeffs in series_coefficients
        ]
        return approximations

    def eoms_approximation(self, order=3, odes_system=None):

        if not odes_system:
            odes_system = self.odes_system

        stiffness_mat = (self.stiffness_matrix())

        eps_stiffness_mat = Matrix(*stiffness_mat.shape, [
            stiff_comp * approx for stiff_comp, approx in zip(
                stiffness_mat, self.eigenvalue_approximation(order=order))
        ])

        odes_system = odes_system + (eps_stiffness_mat -
                                     stiffness_mat) * Matrix(self.dvars)

        eoms_approximated = Matrix(odes_system).subs(
            self.predicted_solution(order, dict=True))

        return eoms_approximated

    def eoms_approximation_list(self,
                                max_order=3,
                                min_order=0,
                                odes_system=None):

        eoms_approximated = self.eoms_approximation(
            order=max_order, odes_system=odes_system).expand()

        return [
            eoms_approximated.applyfunc(
                lambda obj: obj.diff(self.eps, order) / factorial(order)).subs(
                    self.eps, 0).doit()
            for order in range(min_order, max_order + 1)
        ]

    def nth_eoms_approximation(self, order=3):

        eoms_approximated = self.eoms_approximation_list(max_order=order,
                                                         min_order=order)

        return eoms_approximated[0]

    def _determine_secular_terms(self, zeroth_approx):

        sol_zeroth = self._find_nth_solution(zeroth_approx,
                                             order=0,
                                             secular_comps={})

        #eig_min=sqrt(self.eigenvalues()[0])

        #secular_comps = {sin(eig_min*self.ivar),cos(eig_min*self.ivar)}
        secular_comps = sum(sol_zeroth).atoms(sin, cos)

        return secular_comps

    def nth_order_solution(self, order=3):

        eoms_list = self.eoms_approximation_list(max_order=order)

        stiffness_mat = (self.stiffness_matrix())

        eps_stiffness_mat = Matrix(*stiffness_mat.shape, [
            stiff_comp * approx for stiff_comp, approx in zip(
                stiffness_mat, self.eigenvalue_approximation(order=order))
        ])

        secular_components = {
            comp: 0
            for comp in self._determine_secular_terms(
                zeroth_approx=eoms_list[0])
        }
        #        display('secular comps',secular_components)
        approx_dict = {}

        solution = []

        for comp_ord, eoms_nth in enumerate(eoms_list):

            #            display('trutututu',eoms_nth)

            nth_solution = self._find_nth_solution(
                eoms_nth.subs(approx_dict).doit().expand(),
                order=comp_ord,
                secular_comps=secular_components,
                dict=True)
            #            display(nth_solution)
            approx_dict.update(nth_solution)

        return Matrix(list(approx_dict.values()))

    def zeroth_approximation(self, dict=False, equation=False):

        #eoms = Matrix(self.odes_system).subs(self.eps, 0)

        eoms = (self.nth_eoms_approximation(0))

        #        display(eoms, self.approximation_function(order=0))

        solution = LinearODESolution(
            eoms, ivar=self.ivar,
            dvars=self.approximation_function(order=0)).solution()

        #         if equation:
        #             solution = Matrix([
        #                 Eq(lhs, rhs)
        #                 for lhs, rhs in zip(self.approximation_function(
        #                     order=0), list(solution))
        #             ])

        #         if dict:
        #             solution = {
        #                 lhs: rhs
        #                 for lhs, rhs in zip(self.approximation_function(
        #                     order=0), list(solution))
        #             }

        return self._format_solution(
            dvars=self.approximation_function(order=0),
            solution=solution,
            dict=dict,
            equation=equation)

    def _find_nth_solution(self,
                           eoms_nth,
                           order,
                           secular_comps,
                           dict=False,
                           equation=False):

        eoms = eoms_nth
        eoms = (eoms.expand()).applyfunc(
            lambda eqn: TR10(TR8(TR10(eqn).expand()).expand()).expand())

        #         print('='*100)
        #         display(*[row.coeff(comp)  for row in eoms  for comp in secular_comps])
        #         print('='*100)

        eoms = eoms.subs({comp: 0 for comp in secular_comps})

        #        display(eoms)

        solution = LinearODESolution(
            eoms,
            ivar=self.ivar,
            dvars=self.approximation_function(order=order)).solution()

        #         print('Const const')
        #         print(LinearODESolution._const_list)

        return self._format_solution(
            dvars=self.approximation_function(order=order),
            solution=solution,
            dict=dict,
            equation=equation)

    def nth_approximation(self, order=1, dict=False, equation=False):

        #         eoms = Matrix(self.odes_system).subs(
        #             self.predicted_solution(2, dict=True)).diff(self.eps,order).subs(
        #                 self.eps, 0).expand()

        eoms = self.nth_eoms_approximation(order)

        zeroth_approximation = self.zeroth_approximation(dict=True)

        secular_comps = sum(zeroth_approximation.values()).atoms(cos, sin)

        approx = [zeroth_approximation] + [
            self.nth_approximation(order=i, dict=True)
            for i in range(1, order)
        ]

        approx_dict = type({})(ChainMap(*approx))

        #        display('abc',approx_dict)

        return self._find_nth_solution(eoms.subs(approx_dict),
                                       order,
                                       secular_comps,
                                       dict=dict,
                                       equation=equation)

    def first_approximation(self, dict=False, equation=False):
        return self.nth_approximation(order=1, dict=dict, equation=equation)

