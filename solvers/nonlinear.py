from sympy import Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Function, lambdify, factorial,solve
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from ..utilities.timeseries import DataMethods, SpectralMethods, TimeDomainMethods, SpectrumSeries, SpectrumFrame, TimeSeries, TimeDataFrame

from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3

from .linear import LinearODESolution, FirstOrderODE

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
            return TimeDataFrame(
                index=time, data={'solution': nth_order_solution_fun(time)})

    def _format_solution(self, dvars, solution, dict=False, equation=False):

        if equation:
            solution = Matrix(
                [Eq(lhs, rhs) for lhs, rhs in zip(dvars, list(solution))])

        if dict:
            solution = {lhs: rhs for lhs, rhs in zip(dvars, list(solution))}

        return solution

    def approximation_function(self, order):

        dvars_no = len(self.dvars)

        return Matrix([
            me.dynamicsymbols('Y_{' + str(dvar_no) + str(order) +'}')
            for dvar_no, dvar in enumerate(self.dvars)
        ])

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

    
class MultiTimeScaleMethod(LinearODESolution):
    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 ics=None,
                 eps=Symbol('varepsilon'),
                 omega=None,
                 order=1,
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
        self._order=order
        
        self._stored_solution=None
        
        self.secular_eq=set()
        self.int_const=set()

        #        display(omega)
        self.omega = 1
        if omega:
            self.omega = omega
            
        self.ics=ics
        
        self._const_sol={}

                
    def set_solution_order(self,order=1):

        self._order=order

    @property
    def order(self):

        return self._order

    @order.setter
    def order(self,order):
        self._order=order



    @property
    def t_list(self):
        r"""
        Returns the list of time domanin slow varying function
        """

        self._t_list=[t_i(self.ivar) for t_i in symbols(f't_0:{self._order+1}',cls=Function, real=True)]
        return self._t_list


    def __call__(self, ivar, order=1, params_values=None):

        if params_values:
            self.params_values = params_values

        
        if not self._stored_solution:
            self._stored_solution = self.nth_order_solution(order)

        solution=self._stored_solution.subs(self.params_values).subs(self.params_values).n()
            
        if isinstance(ivar, Symbol):
            return solution.subs(self.ivar, ivar)
        else:
            #display(solution)
            #display(solution[1])
            nth_order_solution_fun = lambdify(self.ivar, (solution), 'numpy')

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

    def _format_solution(self, dvars, solution, dict=False, equation=False):

        if equation:
            solution = Matrix(
                [Eq(lhs, rhs) for lhs, rhs in zip(dvars, list(solution))])

        if dict:
            solution = {lhs: rhs for lhs, rhs in zip(dvars, list(solution))}

        return solution



    def approximation_function(self, order,max_order=1):

        dvars_no = len(self.dvars)

        
        t_list=self.t_list
        

# Function('X')(t0(t),t1(t),t2(t)).diff(t,2).subs({t0(t).diff(t,1):1,t1(t).diff(t,1):eps,t2(t).diff(t,1):eps**2}).doit().expand().subs(Function('X')(t0(t),t1(t),t2(t)) ,Function('X0')(t0(t),t1(t),t2(t)) + Function('X1')(t0(t),t1(t),t2(t))*eps + Function('X2')(t0(t),t1(t),t2(t))  *eps**2 ).doit().expand().coeff(eps)
        
        #print(t_list)
    
        return Matrix([
            sym.Function('Y_{' + str(dvar_no) + str(order) +'}')(*t_list) 
            for dvar_no, dvar in enumerate(self.dvars)
        ])

    
    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)


        
        
        solution = sum(
            (self.approximation_function(comp_ord,order) * self.eps**comp_ord
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
        
        
        first_ord_subs={t_i.diff(self.ivar):self.eps**t_ord  for t_ord,t_i  in enumerate(self.t_list)}
        sec_ord_subs={t_i.diff(self.ivar,2):0 for t_ord,t_i  in enumerate(self.t_list)}
        
        
        odes_system = odes_system
        


        

        eoms_approximated = Matrix(odes_system).subs(
            self.predicted_solution(order, dict=True)).subs(sec_ord_subs).doit().subs(first_ord_subs)

        t_fun_subs_dict=({
                expr_tmp:expr_tmp.subs(self.ivar,self.t_list[0]) for expr_tmp in
                (eoms_approximated.atoms(Function) -  {*self.t_list} - { *sum([list(self.approximation_function(ord_tmp,order)) for ord_tmp in range(order+1)],[])  })
            
            })
        
        return eoms_approximated.subs(t_fun_subs_dict)

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

    def _determine_secular_terms(self, zeroth_approx,order,ivar=None):
        
#         print('zeroth approx')
#         display(*zeroth_approx)
        
        
        if not ivar:
            ivar=self.t_list[0]

            
#         print('===')
        sol_zeroth = self._find_nth_solution(zeroth_approx,
                                             order=0,
                                             secular_comps={},ivar=ivar)

#         print('+++ sol of zeroth')
#         display('+'*100,sol_zeroth,'+'*100)
#         print('+++ sol of zeroth')
#         #eig_min=sqrt(self.eigenvalues()[0])

#         #secular_comps = {sin(eig_min*self.ivar),cos(eig_min*self.ivar)}
        secular_comps = sum(sol_zeroth).atoms(sin, cos)
     
        self.int_const|={ eq.coeff(comp)   for eq in sol_zeroth for comp in secular_comps}
    
        

#         display('+'*100,secular_comps,'+'*100)
        
        return secular_comps

    def _nth_order_solution_with_C(self, order=1):

        t_list=self.t_list


        eoms_list = self.eoms_approximation_list(max_order=order)

        stiffness_mat = (self.stiffness_matrix())

        eps_stiffness_mat = Matrix(*stiffness_mat.shape, [
            stiff_comp * approx for stiff_comp, approx in zip(
                stiffness_mat, self.eigenvalue_approximation(order=order))
        ])

        secular_components = {
            comp: 0
            for comp in self._determine_secular_terms(
                zeroth_approx=eoms_list[0],order=order)
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
                dict=True,ivar=self.t_list[0])
            #            display(nth_solution)
            approx_dict.update(nth_solution)
            
        sec_eq=Matrix(list(self.secular_eq - {S.Zero}))
#         display(sec_eq)
#         display(self.int_const)
        
#         const_sol=FirstOrderODE((sec_eq),
#                         self.t_list[1],
#                         dvars=Matrix(list(self.int_const))
#                         ).general_solution()
#         display(const_sol)

#         print('Const const')
#         print(FirstOrderODE._const_list)    
        
        t_back_subs={t_i:self.ivar*self.eps**t_ord for t_ord,t_i  in enumerate(self.t_list)}
        
#         display(t_back_subs)
        
        
        general_form_dict={var: eqn.subs(approx_dict)  for var, eqn in self.predicted_solution(order=order,dict=True).items()}
        
        
#         print('ics')
#         display(*list({var: eqn.subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
#         display(*list({var: eqn.diff(self.ivar).subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
        

                  
#         display(const_vals)
                  
        return general_form_dict

    
    def nth_order_solution(self, order=1):

        general_form_dict=self._nth_order_solution_with_C(order=order)
        
        
        
        
        sec_eq=Matrix(list(self.secular_eq - {S.Zero}))
#         display(sec_eq)
#         display(self.int_const)
        
        const_sol=FirstOrderODE((sec_eq),
                        self.t_list[1],
                        dvars=Matrix(list(self.int_const))
                        ).solution()
        
#         print('trial for SS')
        steady_sol=(FirstOrderODE((sec_eq),
                        self.t_list[1],
                        dvars=Matrix(list(self.int_const))
                        ).steady_solution()) 
        
        
        self._const_sol=const_sol
#         display(const_sol)

#         print('Const const')
#         print(FirstOrderODE._const_list)    
        
        t_back_subs={t_i:self.ivar*self.eps**t_ord for t_ord,t_i  in enumerate(self.t_list)}
        
#         display(t_back_subs)
        
        general_form_dict={key:eq.subs(const_sol).subs(t_back_subs)   for key,eq in general_form_dict.items()}

#         print('ics')
#         display(*list({var: eqn.subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
#         display(*list({var: eqn.diff(self.ivar).subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
        
        ics_eqns=(list({var: eqn.subs({self.ivar:0})  for var, eqn in general_form_dict.items()}.values())+
                    list({var: eqn.diff(self.ivar).subs({self.ivar:0})  for var, eqn in general_form_dict.items()}.values())
                      )
        
        
        
        
        eqns=Matrix([eq.expand() for eq in ics_eqns])
        
        
        
        const_from_eqn=[var for  var in list(eqns.atoms(Symbol,Function)) if var in FirstOrderODE._const_list]
        
        #display(Matrix(ics_eqns).jacobian(const_from_eqn))
        
        #display(const_from_eqn)
        
        #display(ics_eqns)
        const_vals_lin=Matrix([eq.expand() for eq in ics_eqns]).jacobian(const_from_eqn).subs({const:0 for const in const_from_eqn}).subs(self.eps,0)#*Matrix(self.ics)
        const_vals_zth=Matrix([eq.expand() for eq in ics_eqns]).subs({const:0 for const in const_from_eqn})
        
        
#         print('+++++++++  eqns for ics +++++++++++ ')
#         display(const_vals_lin)
#         display(const_vals_zth)
        
        self._ics_formula={const:val for const, val in zip(const_from_eqn,const_vals_lin.inv()*(Matrix(self.ics)-const_vals_zth))}
        
#         display(self._ics_formula)
#         print('+++++++++  eqns for ics +++++++++++ ')      
            
            
        return Matrix(
            list(general_form_dict.values())).subs(const_sol).subs(t_back_subs).subs(self._ics_formula)
        
    
    def zeroth_approximation(self, dict=False, equation=False):

        #eoms = Matrix(self.odes_system).subs(self.eps, 0)

        eoms = (self.nth_eoms_approximation(0))

        #        display(eoms, self.approximation_function(order=0))

        solution = LinearODESolution(
            eoms,
            ivar=self.ivar,
            dvars=self.approximation_function(order=len(self.t_list))).solution()

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
                           equation=False,
                           ivar=None):
        
        if not ivar:
            ivar=self.ivar

        eoms = eoms_nth
        
        # print('='*100,'eoms_nth for solve')
        # display(eoms )
        # print('='*100,'eoms_nth for solve')
        
        eoms = (eoms.expand()).applyfunc(
            lambda eqn: ((TR10(TR8(TR10(eqn).expand()).expand()).expand())))

        # print('='*100,'eoms_nth for solve')
        # display(eoms )
        # print('='*100,'eoms_nth for solve')
        
        self.secular_eq|={row.coeff(comp)  for row in eoms  for comp in secular_comps} #test

#         print('='*100)
#         display(*self.secular_eq)
#         print('='*100)

       
        

        eoms = eoms.subs({comp: 0 for comp in secular_comps})

#         print('='*100)
#         display(eoms)
#         display(self.approximation_function(order=order))
#         display(ivar)
#         print('='*100)

    
        
        
        solution = LinearODESolution(
            eoms,
            ivar=ivar,
            dvars=self.approximation_function(order=order)).solution()

#         print('=eoms and its linear sol'*100)
#         display(eoms)
#         display(solution)

#         print('=eoms and its sol'*100)
        
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

        zeroth_approximation = self.zeroth_approximation(dict=True,order=order)

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
    
    def numerized(self,
                 params_values={},
                 **kwargs):

        
        return self.__class__(
                odes_system=self.governing_equations,
                 ivar=self.ivar,
                 dvars= self.dvars,
                 ics=self.ics,
                 eps=self.eps,
                 omega=self.omega,
                 order=self._order,
                 t_span=[],
                 params=[],
                 params_values=params_values,
                 ic_point=self.ics,
                 equation_type=None)
    
    
    
    def compute_solution(self,
                         t_span=None,
                         ic_list=None,
                         t_eval=None,
                         params_values=None,
                         method='RK45'):
        
        return self.__call__( ivar=t_span, order=self._order, params_values=params_values)
    