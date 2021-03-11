from dynpy import dynamic 
from dynpy import solvers

from sympy import *
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
from sympy.simplify.fu import TR8, TR10, TR7, TR3

import numpy as np
import itertools as itools
import scipy.integrate as solver
from .timeseries import *
from collections import ChainMap
from IPython.display import display


def multivariable_taylor_series(expr, args, n=2, x0=None):
    '''
    Computes the multivariable Taylor series of expresion expr for the given order n in the neighbourhood of x0 point.
    
    Args:
        expr (:obj:Expr): Expresion describing the function under consideration.
        args (array_like): Set of the arguments.
        to be continued ...

    
    Examples:
        Taylor series expansion of the second order for the function $y+x**2$ is as follows:
        
        >>> multivariable_taylor_series(y+x**2,args=[x,y],n=2,x0=[0,0])
        y+x**2
    '''

    order_max = n
    op_point = x0

    args_shifted = {
        arg: arg - arg_shift
        for arg, arg_shift in op_point.items()
    }

    diff_orders = lambda var_list, order_arg: [
        args for tmp in range(order_arg)
    ]

    #     term_tmp_list=sum([list(itools.combinations_with_replacement(var_list,ord_tmp)) for ord_tmp in range(1,4)],[])

    #     term_tmp_dict= {comp:Mul(*comp)/Mul(*[factorial(elem)  for elem  in Poly(Mul(*comp),*var_list).terms()[0][0]]).doit()  for comp in term_tmp_list}
    #     term_tmp_dict

    diff_orders_list = sum([
        list(itools.combinations_with_replacement(args, order))
        for order in range(1, order_max + 1, 1)
    ], [])
    diff_orders_dict = {
        comp: (Mul(*comp).subs(args_shifted) / Mul(*[
            factorial(elem) for elem in Poly(Mul(*comp), *args).terms()[0][0]
        ])).doit()
        for comp in diff_orders_list
    }

    #print(diff_orders_list)

    return (sum([
        expr.diff(*args_tmp).subs(op_point).doit() * poly
        for args_tmp, poly in diff_orders_dict.items()
    ]) + expr.subs(op_point)).doit()

def scalar_fun_quadratic_form(expr, coordinates, op_point):
    '''Vector of deformation '''
    u = (Matrix(coordinates) -
         Matrix(coordinates).subs(op_point, simultaneous=True)).doit()

    constant_term = (expr.subs(op_point, simultaneous=True))
    linear_term = sum((Matrix([expr]).jacobian(coordinates).subs(
        op_point, simultaneous=True)) * u)
    quad_tem = sum(
        S.One / 2 * u.T *
        ((hessian(expr, coordinates).subs(op_point, simultaneous=True))) * u)

    return constant_term + linear_term + quad_tem

class LinearODESolution:
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
        if isinstance(odes_system, LinearDynamicSystem):
            self.dvars = odes_system.q
            odes_system = odes_system._eoms

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

        return self.governing_equations.jacobian(dvars_ddot).subs(
            {coord: 0
             for coord in self.dvars}).doit()

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

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        C = numbered_symbols('C', start=1)

        #         print('o tu')
        #         display(self.odes_system)

        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        Y_mat = Matrix(self.dvars)

        diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        t_sol = self.ivar

        solution = [
            next(C) * modes[:, i] * sin(sym.sqrt(eigs[i, i]) * t_sol) +
            next(C) * modes[:, i] * cos(sym.sqrt(eigs[i, i]) * t_sol)
            for i, coord in enumerate(self.dvars)
        ]

        return sum(solution, Matrix([0] * len(Y_mat)))

    def steady_solution(self, initial_conditions=None):

        ext_forces = self.external_forces().expand().applyfunc(
            lambda row: (TR8(row).expand()))

        #         sin_components=ext_forces.atoms(sin)
        #         cos_components=ext_forces.atoms(cos)

        #        display('sin cos reco',)
        components = ext_forces.atoms(sin, cos)

        #        display('ext_forces',ext_forces)

        steady_sol = Matrix([0 for gen_coord in self.dvars])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
            #            display(omg)
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

            #display(amp_vector)

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (sym.re(fund_mat.inv()).doit() * amp_vector) * comp + (sym.im(fund_mat.inv()).doit() * amp_vector) * TR3(comp.subs(omg*self.ivar,omg*self.ivar-pi/2))
            
            #steady_sol += ((fund_mat.inv()) * amp_vector) * comp 

            ext_forces -= (amp_vector * comp).expand()
        #print(ext_forces.doit().expand())

        const_elems = lambda expr: [
            expr for expr in comp.expand().args if not expr.has(self.ivar)
            if isinstance(comp, Add)
        ]

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

        return steady_sol

    def solution(self, initial_conditions=None):
        return self.general_solution(
            initial_conditions=initial_conditions) + self.steady_solution(
                initial_conditions=initial_conditions)

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


        display(self.omega)

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
            me.dynamicsymbols('Y_' + str(dvar_no) + str(order))
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

