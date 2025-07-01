import copy
import itertools as itools
from collections import ChainMap
from collections.abc import Iterable
from functools import cached_property

# from timer import timer
from time import time

import numpy as np
import scipy.integrate as solver
import sympy as sym
import sympy.physics.mechanics as me
from IPython.display import display
from sympy import (
    Add,
    Dict,
    Dummy,
    Eq,
    Function,
    Matrix,
    Mul,
    N,
    Number,
    S,
    Symbol,
    cos,
    det,
    diag,
    diff,
    exp,
    expand,
    factorial,
    lambdify,
    sin,
    solve,
    sqrt,
    symbols,
    zoo,
)
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vlatex, vpprint
from sympy.simplify.fu import TR0, TR3, TR7, TR8, TR10
from sympy.utilities.autowrap import autowrap, ufuncify

from ..utilities.components.ode import en as ode_comp
from ..utilities.timeseries import (
    DataMethods,
    SpectralMethods,
    SpectrumFrame,
    SpectrumSeries,
    TimeDataFrame,
    TimeDomainMethods,
    TimeSeries,
)
from .linear import (
    AnalyticalSolution,
    FirstOrderLinearODESystem,
    FirstOrderLinearODESystemWithHarmonics,
    FirstOrderODE,
    FirstOrderODESystem,
    LinearODESolution,
    ODESolution,
    ODESystem,
)
from .tools import ODE_COMPONENTS_LIST, CodeFlowLogger, CommonFactorDetector


def const_no_gen():
    num = 0
    while True:
        yield num
        num += 1


class SimplifiedExpr:
    """

    This class provides methods which allow to obtain a differential equation in simplified form, or extract desired simplification parameters.

    """

    _const_no = const_no_gen()

    _subs_container = {}

    def __init__(self, expr, ivar, parameters=None):

        self._expr = expr
        self._ivar = ivar
        self._parameters = parameters

    @property
    def _fun_list(self):
        """

        Method extracts functions from expression and returns them in a list.

        """

        funlist = [
            expr for expr in (self._expr).atoms(sin, cos) if expr.has(self._ivar)
        ]

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
        expr_dict = {key: values[0] for key, values in data.items()}
        self.__class__._subs_container = {**self.__class__._subs_container, **expr_dict}

        simplified_sum = sum([key * values[1] for key, values in data.items()], S.Zero)
        expanded_sum = sum(
            [values[0] * values[1] for key, values in data.items()], S.Zero
        )

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

        expr_dict = {key: values[0] for key, values in self._expr_collector.items()}
        return expr_dict

    @property
    def _expr_collector(self):
        """

        This method combines functions and coefficients and then returns a dictionary with segregated values.

        """

        parameters = self._parameters

        if parameters:
            aux_symbol = lambda no: Function(f"DummyFun_{next(self._const_no)}")(
                *parameters
            )
        else:
            aux_symbol = lambda no: Dummy(f"D_{no}")

        bucket = {
            aux_symbol(no): (self._expr.coeff(fun), fun)
            for no, fun in enumerate(self._fun_list)
        }
        return bucket

    @property
    def full_expr(self):
        """

        Main method that returns a full expression with functions and proper coefficients.

        """

        # display(self.__class__._subs_container)

        return self._expr.doit().expand().subs(self.__class__._subs_container)

    def __repr__(self, *args):
        return "instance of SimplifiedExpr"


class PerturbationODESolution(ODESolution):
    _eps = Symbol("varepsilon")

    def set_small_parameter(self, eps, inplace=False):
        obj = copy.copy(self)
        obj._eps = eps

        return obj

    @property
    def small_parameter(self):
        return self._eps

    @small_parameter.setter
    def small_parameter(self, eps):

        self._eps = eps

        return self._eps

    @property
    def eps(self):
        return self.small_parameter

    def _calculate_constant(self, ics=None, *args, **kwargs):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """

        const_eqns = self._get_constant_eqns(ics).subs(self.eps, 0).doit()
        const_list = self._spot_constant()

        # display('ivar',self.ivar)

        return solve(const_eqns, const_list)


class ODESystemApproximation(ODESystem):

    _simp_dict = {}
    _callback_dict = {}

    @property
    def _default_solver(self):

        return FirstOrderLinearODESystemWithHarmonics
        # return FirstOrderLinearODESystem


class ODESystemDsolve(ODESystem):
    _ode_order = 2
    _simp_dict = {}
    _callback_dict = {}

    @property
    def _default_solver(self):

        return FirstOrderLinearODESystemWithHarmonics
        # return FirstOrderLinearODESystem


class NthOrderODEsApproximation(ODESystem):
    _ode_order = 2

    _simp_dict = {}
    _callback_dict = {}

    # _const_list = []

    @property
    def _default_solver(self):

        return FirstOrderLinearODESystemWithHarmonics
        # return FirstOrderLinearODESystem

    @property
    def _secular_funcs(self):

        # weird error occurs, when .remove_secular_terms is called. if secular terms are missing, method returns None,
        # if statement should be considered

        # const_to_add=[eqn for eqn in self._const_list ]

        secular_funcs = set() | (
            self.general_solution.rhs.atoms(Function)
            - set(self._const_list)
            - set(self._parameters)
            - {self.ivar}
        )
        # display(secular_funcs)

        return {sec_fun for sec_fun in secular_funcs if sec_fun.has(self.ivar)}

    @property
    def secular_terms(self):

        # print('secsec')
        sec_funcs = self._secular_funcs
        # self._const_list

        display(sec_funcs)

        sec_conditions = Matrix(
            list(
                set(
                    sum(
                        [
                            list(self.odes.applyfunc(lambda entry: entry.coeff(func)))
                            for func in sec_funcs
                        ],
                        [],
                    )
                )
                - {0}
            )
        )

        ivar = self._parameters[0]

        # sec_conditions =FirstOrderODESystem.from_ode_system(ODESystem([sec_conditions[1][1],sec_conditions[3][1]],dvars = self._const_list ,ivar=ivar  )).linearized().solution

        # print('some check fof sec eq')

        display(sec_conditions)

        const_list = [
            fun_cons
            for fun_cons in list(sec_conditions.atoms(Function))
            if "C" in str(fun_cons)
        ]

        display(const_list)

        sec_odes = ODESystem(
            (sec_conditions), dvars=Matrix(const_list), ivar=ivar, ode_order=None
        )

        display(sec_odes.as_type(FirstOrderLinearODESystem).solution)
        # fode_harm_rhs = FirstOrderLinearODESystemWithHarmonics.from_ode_system(sec_odes)

        # print('transformed')
        # display(fode_harm_rhs)
        # display(type(fode_harm_rhs))

        return sec_odes

    def remove_secular_terms(self):

        obj = ODESystem.from_ode_system(self)

        # subs_dict = {comp: 0 for comp in self._secular_funcs}

        # return obj.subs(subs_dict).doit().subs(zoo,0).doit()

        secular_expr_list = [
            obj.applyfunc(lambda row: row.coeff(comp) * comp).rhs
            for comp in self._secular_funcs
        ]
        secular_expr = sum(secular_expr_list, Matrix([0] * len(obj)))

        return ODESystem(
            (obj.as_matrix() - secular_expr)
            .expand()
            .subs({comp: 0 for comp in self._secular_funcs})
            .doit(),
            obj.dvars,
            ivar=obj.ivar,
        )


#     def find_approximation(self):

#         zeroth_ord = self.odes.applyfunc(lambda ode: ode.general_solution.self.remove_secular_terms())

#         for order,approx in enumerate(self.odes):

#             approx._parameters = self.t_list[order+1:]

#             nth_ord = approx.applyfunc(lambda ode: ode.steady_solution.self.remove_secular_terms())

#         return {order:}


class MultiTimeScaleSolution(ODESystem):
    _ode_order = 2
    _stored_solution = None
    _saved_solution = {}

    _saved_numerized = {}
    _order = 2
    _scales_no = None

    def __new__(
        cls,
        odes_system,
        dvars,
        ivar=None,
        ics=None,
        eps=Symbol("varepsilon"),
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
        **options,
    ):

        if not isinstance(odes_system, Iterable):
            odes_system = Matrix([odes_system])
        if not isinstance(dvars, Iterable):
            dvars = Matrix([dvars])

        obj = super().__new__(
            cls, odes_system, dvars=dvars, ivar=ivar, evaluate=evaluate, **options
        )

        obj._const_list = []

        if ivar is not None:
            obj._ivar = ivar

        if label == None:
            label = obj._label = (
                obj.__class__.__name__ + " with " + str(len(obj.dvars)) + " equations"
            )

        obj._label = 123

        obj.eps = eps
        obj._order = order

        if scales_no is None:
            obj._scales_no = obj.order + 1

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
        # obj._odes = odes_system

        return obj

    @property
    def eps(self):
        return self._eps

    @eps.setter
    def eps(self, eps):
        self._eps = eps

    @property
    def order(self):
        return self._order

    @property
    def scales_no(self):
        if self._scales_no is None:
            return self.order + 1
        else:
            return self._scales_no

    @scales_no.setter
    def scales_no(self, no):

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
        return "123"

    def __repr__(self):
        return "123"

    def set_solution_order(self, order=1):

        self._order = order

    @property
    def _report_components(self):

        comp_list = [
            ode_comp.ODESystemComponent,
            ode_comp.VariablesComponent,
            ode_comp.ODESystemCodeComponent,
            ode_comp.MSMCalculationsOrderComponent,
            ode_comp.PredictedSolutionComponent,
            ode_comp.DetailsOfPredictedSolutionComponent,
            ode_comp.GoverningEquationComponent,
            # ode_comp.ApproximatedNonlinearGoverningEquationComponent,
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
    def order(self, order):
        self._order = order

    def set_order(self, order):

        new_msm = MultiTimeScaleSolution(
            self.as_matrix(),
            self.vars,
            ivar=self.ivar,
            omega=S.One,
            order=order,
            scales_no=self.scales_no,
            eps=self.eps,
        )
        new_msm._scales_no = self._scales_no

        return new_msm

    def _assign_properties(self, obj):

        obj._eps = self.eps
        obj._vars = self._vars
        # obj._const_list = self._const_list

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
            for t_i in symbols(f"t_0:{self.scales_no}", cls=Function, real=True)
        ]

        return self._t_list

    @cached_property
    def _scales_formula(self):
        r"""
        Returns the list of time domanin slow varying function
        """

        t_list = self.t_list
        expr_list = [self.ivar * self.eps**no for no in range(len(t_list))]

        return AnalyticalSolution.from_vars_and_rhs(t_list, expr_list)

    def approx_function_symbol(self, order=0):

        dvars_no = len(self.dvars)

        if dvars_no == 1:

            func_str = str(self.vars[0]).replace(f"({self.ivar})", "")
        else:
            func_str = "Y"

        return func_str + f"_{{{str(order)}}}"

    def approximation_function(self, order, max_order=1, ivar_list=None):

        dvars_no = len(self.vars)

        if ivar_list is None:
            t_list = self.t_list
        else:
            t_list = ivar_list

        fun_str = self.approx_function_symbol(order)

        if dvars_no == 1:

            aprrox_fun = [
                sym.Function(fun_str)(*t_list) for dvar_no, dvar in enumerate(self.vars)
            ]
        else:

            aprrox_fun = [
                sym.Function(fun_str + "_{" + str(order) + "}")(*t_list)
                for dvar_no, dvar in enumerate(self.vars)
            ]

        return Matrix(aprrox_fun)

    def approximation_function_without_scales(self, order, max_order=1):

        return self.approximation_function(
            order=order, max_order=max_order, ivar_list=[self.ivar]
        )

    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.vars)

        solution = sum(
            (
                self.approximation_function(comp_ord, order) * self.eps**comp_ord
                for comp_ord in range(order + 1)
            ),
            sym.zeros(dvars_no, 1),
        )

        return ODESolution.from_vars_and_rhs(self.vars, solution)

    def predicted_solution_without_scales(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)

        solution = sum(
            (
                self.approximation_function_without_scales(comp_ord, order)
                * self.eps**comp_ord
                for comp_ord in range(order + 1)
            ),
            sym.zeros(dvars_no, 1),
        )

        return ODESolution.from_vars_and_rhs(self.dvars, solution)

    def derivative_without_scale(self, degree, order=1, dict=False, equation=False):

        derivative_rhs = (
            self.predicted_solution_without_scales(self.order)
            .as_eq_list()[0]
            .n(3)
            .rhs.diff(self.ivar, degree)
        )
        derivative_lhs = self.dvars[0].diff(self.ivar, degree)

        return AnalyticalSolution(
            derivative_lhs, rhs=derivative_rhs, vars=self.dvars[0]
        )

    def first_order_subs(self):

        first_ord_subs = {
            t_i.diff(self.ivar): self.eps**t_ord
            for t_ord, t_i in enumerate(self.t_list)
        }

        return AnalyticalSolution.from_dict(first_ord_subs)

    def second_order_subs(self):

        sec_ord_subs = {
            t_i.diff(self.ivar, 2): 0 for t_ord, t_i in enumerate(self.t_list)
        }

        return AnalyticalSolution.from_dict(sec_ord_subs)

    def eoms_approximation(self, order=None, odes_system=None):

        if order is None:
            order = self.order

        if not odes_system:
            odes_system = self.as_matrix()

        first_ord_subs = self.first_order_subs().as_explicit_dict()
        sec_ord_subs = self.second_order_subs().as_explicit_dict()

        # display(self.predicted_solution(order).as_dict().subs(sec_ord_subs).doit().subs(first_ord_subs).doit())

        odes_rhs = (
            self.as_matrix()
            .subs(self.predicted_solution(order).as_explicit_dict())
            .subs(sec_ord_subs)
            .doit()
            .subs(first_ord_subs)
            .doit()
        )

        # display()
        # eoms_approximated = ODESystem(odes_system,dvars=self.dvars,parameters = self.t_list).subs(
        #    self.predicted_solution(order).as_dict())
        eoms_approximated = -odes_rhs

        eoms_approximated = (
            eoms_approximated.subs(sec_ord_subs).doit().subs(first_ord_subs).doit()
        )

        t_fun_subs_dict = {
            expr_tmp: expr_tmp.subs(self.ivar, self.t_list[0])
            for expr_tmp in (
                eoms_approximated.atoms(Function)
                - {*self.t_list}
                - {
                    *sum(
                        [
                            list(self.approximation_function(ord_tmp, order))
                            for ord_tmp in range(order + 1)
                        ],
                        [],
                    )
                }
            )
        }

        CodeFlowLogger(eoms_approximated, "eoms for approximation before subs", self)
        CodeFlowLogger(t_fun_subs_dict, "dict with time subs", self)

        return eoms_approximated.subs(t_fun_subs_dict).doit()

    def eoms_approximation_list(self, max_order=3, min_order=0, odes_system=None):

        eoms_approximated = self.eoms_approximation(
            order=max_order, odes_system=odes_system
        ).expand()

        order_get = (
            lambda obj, order: (obj.diff(self.eps, order) / factorial(order))
            .subs(self.eps, 0)
            .doit()
        )

        # NthOrderODEsApproximation

        approx_list = [
            NthOrderODEsApproximation(
                eoms_approximated.applyfunc(lambda obj: order_get(obj, order)),
                dvars=self.approximation_function(order),
                ivar=self.t_list[0],
                ode_order=2,
                parameters=self.t_list,
            )
            for order in range(min_order, max_order + 1)
        ]

        return approx_list
        # return [NthOrderODEsApproximation.from_ode_system(approx_ode)   for approx_ode in approx_list]

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
        sol_subs_list = [
            sol.applyfunc(
                lambda row: SimplifiedExpr(
                    row, ivar=self._t_list[0], parameters=self._t_list[1:]
                ).full_expr
            ).doit()
            for sol in sol_list
        ]
        sol_subs_dict = {
            dvar: eqn
            for sol in sol_subs_list
            for dvar, eqn in sol.as_explicit_dict().items()
        }

        # print('_gen_sol')
        # display(sol_list)
        # display(sol_subs_dict)

        # display(approx_eoms_list[0]._const_list)
        approx_with_const = [approx_eoms_list[0]]

        for order, approx in enumerate(approx_eoms_list[1:]):

            approx._parameters = self._t_list[1:]
            # approx = approx.as_first_ode_linear_system()

            eqns_map = lambda obj: (
                TR10(TR8(TR10(obj.expand()).expand())).expand().doit().expand()
            )

            def eqns_map(obj):

                oper_expr = lambda obj: (
                    TR10(TR8(TR10(TR0(obj.expand())).expand())).expand().doit().expand()
                )

                oper_expr = lambda obj: obj.doit().expand()

                if isinstance(obj, Add):
                    elems = obj.args
                    return sum(oper_expr(elem) for elem in elems)
                else:
                    return oper_expr(obj)

            # print('approx eqn')
            # display(approx)

            sol_subs_list = [
                sol.applyfunc(
                    lambda row: SimplifiedExpr(
                        row, ivar=self._t_list[0], parameters=self._t_list[1:]
                    ).full_expr
                ).doit()
                for sol in sol_list
            ]
            sol_subs_dict = {
                dvar: eqn
                for sol in sol_subs_list
                for dvar, eqn in sol.as_explicit_dict().items()
            }

            approx_subs = (
                approx.applyfunc(eqns_map).subs(sol_subs_dict).applyfunc(eqns_map)
            )

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
        sol_subs_list = [
            sol.applyfunc(
                lambda row: SimplifiedExpr(
                    row, ivar=self._t_list[0], parameters=self._t_list[1:]
                ).full_expr
            ).doit()
            for sol in sol_list
        ]
        sol_subs_dict = {
            dvar: eqn
            for sol in sol_subs_list
            for dvar, eqn in sol.as_explicit_dict().items()
        }

        # print('_gen_sol')
        # display(sol_list)
        # display(sol_subs_dict)

        # display(approx_eoms_list[0]._const_list)

        for order, approx in enumerate(approx_eoms_list[1:]):

            approx._parameters = self._t_list[1:]
            # approx = approx.as_first_ode_linear_system()

            eqns_map = lambda obj: (
                TR10(TR8(TR10(obj.expand()).expand())).expand().doit().expand()
            )

            def eqns_map(obj):

                oper_expr = lambda obj: (
                    TR10(TR8(TR10(TR0(obj.expand())).expand())).expand().doit().expand()
                )

                oper_expr = lambda obj: obj.doit().expand()

                if isinstance(obj, Add):
                    elems = obj.args
                    return sum(oper_expr(elem) for elem in elems)
                else:
                    return oper_expr(obj)

            # print('approx eqn')
            # display(approx)

            sol_subs_list = [
                sol.applyfunc(
                    lambda row: SimplifiedExpr(
                        row, ivar=self._t_list[0], parameters=self._t_list[1:]
                    ).full_expr
                ).doit()
                for sol in sol_list
            ]
            sol_subs_dict = {
                dvar: eqn
                for sol in sol_subs_list
                for dvar, eqn in sol.as_explicit_dict().items()
            }

            approx_subs = (
                approx.applyfunc(eqns_map).subs(sol_subs_dict).applyfunc(eqns_map)
            )

            # display(approx_subs)

            approx_subs = NthOrderODEsApproximation.from_ode_system(approx_subs)
            approx_subs._parameters = self._t_list[1:]

            # display(approx_subs)
            # display(approx_subs.ivar)
            # display(approx_subs.secular_terms)
            # display(type(approx_subs.secular_terms))
            print("HERE")
            display(approx_subs)
            self.secular_eq[self.eps ** (order + 1)] = (
                approx_subs.secular_terms.as_type(FirstOrderLinearODESystem)
            )

            # print('with secular terms')
            # display(approx_subs)
            approx_subs = approx_subs.remove_secular_terms()

            # print('without secular terms')
            # display(approx_subs)

            # display(approx_subs.lhs,approx_subs.rhs)
            # display(approx_subs)
            # display(FirstOrderLinearODESystemWithHarmonics.from_ode_system(approx_subs).steady_solution)
            # display(type(approx_subs))

            # aux_mat=self.odes_system.jacobian(self.dvars)
            aux_mat = approx_subs.jacobian(approx_subs.dvars)
            if det(aux_mat) == 0:

                SolverClass = FirstOrderLinearODESystem
            else:
                SolverClass = FirstOrderLinearODESystemWithHarmonics

            # SolverClass = FirstOrderLinearODESystemWithHarmonics
            # SolverClass = FirstOrderLinearODESystem

            ode_2_solve = ODESystem.from_ode_system(approx_subs)

            ode_2_solve._ivar = self._t_list[0]

            ode_2_solve = SolverClass.from_ode_system(ode_2_solve)
            # print('gen sol part')
            # display(*FirstOrderLinearODESystemWithHarmonics.from_ode_system(ode_2_solve)._get_excitation_comps)

            # fm_mat = ode_2_solve._as_fode()._fundamental_matrix
            # display(fm_mat)
            # display(fm_mat.diagonalize())

            sol = ode_2_solve.steady_solution.applyfunc(
                lambda obj: obj.expand()
            ).applyfunc(
                eqns_map
            )  # .applyfunc(lambda row: SimplifiedExpr(row,ivar=self._t_list[0],parameters=self._t_list[1:]).sum_expr)

            # print('and its sol')
            # display(sol)
            # sol_subs_dict = {**sol_subs_dict, **sol.as_dict()}
            # display(*list(sol.as_matrix()))
            sol_list += [
                sol.applyfunc(
                    lambda row: SimplifiedExpr(
                        row, ivar=self._t_list[0], parameters=self._t_list[1:]
                    ).full_expr
                ).doit()
            ]
            # sol_list += [sol]
        print("sol list")
        display(*sol_list)

        return sol_list

    @property
    def general_solution(self):

        order = self.order

        # print('my order',order)
        sol_list = []
        gen_sol = self._general_sol(order)

        # print('gen_sol')
        # display(*gen_sol)

        if self.eps in self.secular_eq:
            print("sec_eq")
            ode2check = self.secular_eq[self.eps]  # .as_first_ode_linear_system()
            print(type(ode2check))
            display(ode2check)
            # display(ode2check._as_fode())
            # display(ode2check.solution)
            # ode2check._as_fode()
            # ode2check.solution
            C_const_sol = ode2check.solution.as_explicit_dict()

            print("tut")
            display(C_const_sol)

        else:
            C_const_sol = {}
        # print('tut')
        for order, approx_sol in enumerate(gen_sol):

            # print('stale: ')
            # display(C_const_sol)
            # print('aproox bez subsa')

            # display(approx_sol.rhs.applyfunc(lambda obj: obj.expand().doit().subs(C_const_sol).subs(C_const_sol)
            #                                ))
            solution = approx_sol.applyfunc(
                lambda obj: obj.expand().doit().subs(C_const_sol).subs(C_const_sol)
            ).subs({self.t_list[1]: self.eps * self.ivar})

            # print('aproox ze subsa')
            # display(solution)

            sol_list += [
                (self.eps**order)
                * solution.subs(
                    {self.t_list[1]: self.eps * self.ivar, self.t_list[0]: self.ivar}
                ).rhs
            ]

        # display(*list(SimplifiedExpr._subs_container.values()))
        result = (sum(sol_list, Matrix(2 * len(self.dvars) * [0]))).applyfunc(
            lambda obj: obj.expand().doit()
        )

        new_res = PerturbationODESolution.from_vars_and_rhs(
            Matrix(list(self.dvars) + list(self.dvars.diff(self.ivar))), result
        )  # .set_small_parameter(self.eps)
        new_res.small_parameter = self.eps
        new_res.ivar = self.ivar

        return new_res

    @cached_property
    def steady_solution(self):

        result = Matrix(2 * len(self.dvars) * [0])

        new_res = PerturbationODESolution(
            Matrix(list(self.dvars) + list(self.dvars.diff(self.ivar))), result
        )  # .set_small_parameter(self.eps)
        new_res.small_parameter = self.eps
        new_res.ivar = self.ivar
        return new_res

    @cached_property
    def solution(self):

        return self.general_solution

    def compute_solution(
        self, t_span=None, ic_list=None, t_eval=None, params_values=None, method="RK45"
    ):

        #         print('compute_solution')
        #         display(params_values)
        #         display(ic_list)

        if ic_list:
            print("class ics has been taken")
            self.ics = ic_list

        return self.__call__(
            ivar=t_span, order=self._order, params_values=params_values
        )


class WeakNonlinearProblemSolution(LinearODESolution):

    def __init__(
        self,
        odes_system,
        ivar=Symbol("t"),
        dvars=[],
        eps=Symbol("varepsilon"),
        omega=None,
        t_span=[],
        params=[],
        params_values={},
        ic_point={},
        equation_type=None,
    ):

        super().__init__(
            odes_system=odes_system,
            ivar=ivar,
            dvars=dvars,
            t_span=t_span,
            params=params,
            params_values=params_values,
            ic_point=ic_point,
            equation_type=equation_type,
        )
        self.eps = eps

        #        display(omega)
        self.omega = 1
        if omega:
            self.omega = omega

    def __call__(self, time, order=1, params_values=None):

        if params_values:
            self.params_values = params_values

        if isinstance(time, Symbol):
            return (
                self.nth_order_solution(order)
                .subs(self.ivar, time)
                .subs(self.params_values)
            )
        else:
            nth_order_solution_fun = lambdify(
                self.ivar,
                self.nth_order_solution(order).subs(self.params_values),
                "numpy",
            )

            #            return nth_order_solution_fun(time)
            result = (
                TimeDataFrame(
                    index=time, data={"solution": (nth_order_solution_fun(time))}
                )
                .copy()
                .apply(np.real)
            )

            display(result)
            display(result.apply(np.real))

            return result

    def _format_solution(self, dvars, solution, dict=False, equation=False):

        if equation:
            solution = Matrix([Eq(lhs, rhs) for lhs, rhs in zip(dvars, list(solution))])

        if dict:
            solution = {lhs: rhs for lhs, rhs in zip(dvars, list(solution))}

        return solution

    def approximation_function(self, order, max_order=1):

        dvars_no = len(self.dvars)

        t_list = self.t_list

        aprrox_fun = [
            sym.Function("Y_{" + str(dvar_no) + str(order) + "}")(*t_list)
            for dvar_no, dvar in enumerate(self.dvars)
        ]

        return Matrix(aprrox_fun)

    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)

        solution = sum(
            (
                self.approximation_function(comp_ord) * self.eps**comp_ord
                for comp_ord in range(order + 1)
            ),
            sym.zeros(dvars_no, 1),
        )

        return self._format_solution(
            dvars=self.dvars, solution=solution, dict=dict, equation=equation
        )

    def eigenvalues(self):
        modes, eigs = (
            self.inertia_matrix().inv() * self.stiffness_matrix()
        ).diagonalize()

        return [eig for eig in eigs if eig != 0]

    def eigenvalue_approximation(
        self,
        order=1,
        eigenvalue=None,
        eigenvalue_no=None,
        eigenvalue_symbol=Symbol("omega"),
        series_coefficients=None,
    ):

        if not eigenvalue:
            eigenvalue = self.omega
        self.omega = eigenvalue

        if not series_coefficients:
            series_coefficients = [
                symbols("A_{}{}_1:{}".format(eom_no, no, order + 1))
                for eom_no, eom in enumerate(self.dvars)
                for no, coord in enumerate(self.dvars)
            ]

        # approximation={eigenvalue:eigenvalue_symbol**2 - sum(coeff**(eps_ord+1)  for eps_ord,coeff  in enumerate(series_coefficients) )}
        approximations = [
            eigenvalue**2
            - sum(
                [
                    coeff * self.eps ** (eps_ord + 1)
                    for eps_ord, coeff in enumerate(coeffs)
                ]
            )
            for coeffs in series_coefficients
        ]
        return approximations

    def eoms_approximation(self, order=3, odes_system=None):

        if not odes_system:
            odes_system = self.odes_system

        stiffness_mat = self.stiffness_matrix()

        eps_stiffness_mat = Matrix(
            *stiffness_mat.shape,
            [
                stiff_comp * approx
                for stiff_comp, approx in zip(
                    stiffness_mat, self.eigenvalue_approximation(order=order)
                )
            ],
        )

        odes_system = odes_system + (eps_stiffness_mat - stiffness_mat) * Matrix(
            self.dvars
        )

        eoms_approximated = Matrix(odes_system).subs(
            self.predicted_solution(order, dict=True)
        )

        return eoms_approximated

    def eoms_approximation_list(self, max_order=3, min_order=0, odes_system=None):

        eoms_approximated = self.eoms_approximation(
            order=max_order, odes_system=odes_system
        ).expand()

        return [
            eoms_approximated.applyfunc(
                lambda obj: obj.diff(self.eps, order) / factorial(order)
            )
            .subs(self.eps, 0)
            .doit()
            for order in range(min_order, max_order + 1)
        ]

    def nth_eoms_approximation(self, order=3):

        eoms_approximated = self.eoms_approximation_list(
            max_order=order, min_order=order
        )

        return eoms_approximated[0]

    def _determine_secular_terms(self, zeroth_approx):

        sol_zeroth = self._find_nth_solution(zeroth_approx, order=0, secular_comps={})

        # eig_min=sqrt(self.eigenvalues()[0])

        # secular_comps = {sin(eig_min*self.ivar),cos(eig_min*self.ivar)}
        secular_comps = sum(sol_zeroth).atoms(sin, cos)

        return secular_comps

    def nth_order_solution(self, order=3):

        eoms_list = self.eoms_approximation_list(max_order=order)

        stiffness_mat = self.stiffness_matrix()

        eps_stiffness_mat = Matrix(
            *stiffness_mat.shape,
            [
                stiff_comp * approx
                for stiff_comp, approx in zip(
                    stiffness_mat, self.eigenvalue_approximation(order=order)
                )
            ],
        )

        secular_components = {
            comp: 0
            for comp in self._determine_secular_terms(zeroth_approx=eoms_list[0])
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
                dict=True,
            )
            #            display(nth_solution)
            approx_dict.update(nth_solution)

        return Matrix(list(approx_dict.values()))

    def zeroth_approximation(self, dict=False, equation=False):

        # eoms = Matrix(self.odes_system).subs(self.eps, 0)

        eoms = self.nth_eoms_approximation(0)

        #        display(eoms, self.approximation_function(order=0))

        solution = LinearODESolution(
            eoms, ivar=self.ivar, dvars=self.approximation_function(order=0)
        ).solution()

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
            equation=equation,
        )

    def _find_nth_solution(
        self, eoms_nth, order, secular_comps, dict=False, equation=False
    ):

        eoms = eoms_nth
        eoms = (eoms.expand()).applyfunc(
            lambda eqn: TR10(TR8(TR10(eqn).expand()).expand()).expand()
        )

        #         print('='*100)
        #         display(*[row.coeff(comp)  for row in eoms  for comp in secular_comps])
        #         print('='*100)

        eoms = eoms.subs({comp: 0 for comp in secular_comps})

        #        display(eoms)

        solution = LinearODESolution(
            eoms, ivar=self.ivar, dvars=self.approximation_function(order=order)
        ).solution()

        #         print('Const const')
        #         print(LinearODESolution._const_list)

        return self._format_solution(
            dvars=self.approximation_function(order=order),
            solution=solution,
            dict=dict,
            equation=equation,
        )

    def nth_approximation(self, order=1, dict=False, equation=False):

        #         eoms = Matrix(self.odes_system).subs(
        #             self.predicted_solution(2, dict=True)).diff(self.eps,order).subs(
        #                 self.eps, 0).expand()

        eoms = self.nth_eoms_approximation(order)

        zeroth_approximation = self.zeroth_approximation(dict=True)

        secular_comps = sum(zeroth_approximation.values()).atoms(cos, sin)

        approx = [zeroth_approximation] + [
            self.nth_approximation(order=i, dict=True) for i in range(1, order)
        ]

        approx_dict = type({})(ChainMap(*approx))

        #        display('abc',approx_dict)

        return self._find_nth_solution(
            eoms.subs(approx_dict), order, secular_comps, dict=dict, equation=equation
        )

    def first_approximation(self, dict=False, equation=False):
        return self.nth_approximation(order=1, dict=dict, equation=equation)


class MultiTimeScaleMethod(LinearODESolution):

    _stored_solution = None
    _saved_solution = {}

    _saved_numerized = {}

    def __init__(
        self,
        odes_system,
        ivar=Symbol("t"),
        dvars=[],
        ics=None,
        eps=Symbol("varepsilon"),
        omega=None,
        order=1,
        t_span=[],
        params=[],
        params_values={},
        extra_params={},
        ic_point={},
        equation_type=None,
        label=None,
    ):

        super().__init__(
            odes_system=odes_system,
            ivar=ivar,
            dvars=dvars,
            t_span=t_span,
            params=params,
            params_values=params_values,
            ic_point=ic_point,
            equation_type=equation_type,
        )
        if label == None:
            label = self._label = (
                self.__class__.__name__ + " with " + str(len(self.dvars)) + " equations"
            )
        self._label = label

        self.eps = eps
        self._order = order

        self._stored_solution = None
        self._saved_solution = {}
        self.extra_params = extra_params

        print(extra_params)

        self.secular_eq = set()
        self.int_const = set()

        #        display(omega)
        self.omega = 1
        if omega:
            self.omega = omega

        self.ics = ics

        self._const_sol = {}

        self._eoms = self.odes_system

    def __str__(self):
        return self._label

    def __repr__(self):
        return self._label

    def set_solution_order(self, order=1):

        self._order = order

    @property
    def order(self):

        return self._order

    @property
    def ics_dvars(self):

        q = list(self.dvars)
        dq = list(Matrix(self.dvars).diff(self.ivar))

        return Matrix(q + dq)

    @order.setter
    def order(self, order):
        self._order = order

    @property
    def t_list(self):
        r"""
        Returns the list of time domanin slow varying function
        """

        self._t_list = [
            t_i(self.ivar)
            for t_i in symbols(f"t_0:{self._order+1}", cls=Function, real=True)
        ]
        return self._t_list

    def _numbers_dict(self, data):

        return {
            key: float(N(value))
            for key, value in data.items()
            if isinstance(N(value), Number)
        }

    def _notnumbers_dict(self, data):

        return {
            key: value
            for key, value in data.items()
            if not isinstance(N(value), Number)
        }

    def __call__(self, ivar, order=1, params_values=None):

        if params_values:
            self.params_values = params_values

        solution = self.nth_order_solution(order).rhs

        solution = solution.applyfunc(
            lambda x: x.subs(self.extra_params).subs(self.params_values)
        )

        if isinstance(ivar, Symbol):

            ics_dict = self._ic_from_sol(order=order, formula=True)
            display(ics_dict)

            return solution.subs(ics_dict).subs(self.extra_params).subs(self.ivar, ivar)
        else:

            ics_dict = self._ic_from_sol(order=order, formula=False)

            ics_symbols_dict = {
                sym: val for sym, val in zip(self.ics_symbols, self.ics)
            }

            ######################### old #############

            # nth_order_solution_fun = lambdify(self.ivar, (solution.subs(ics_dict).subs(self.extra_params).subs(ics_symbols_dict)).subs(self.params_values).n(), 'numpy')
            ######################### old #############

            notnumbers_dict = {
                **self._notnumbers_dict(self.extra_params),
                **self._notnumbers_dict(self.params_values),
            }
            numbers_dict = {
                **self._numbers_dict(self.extra_params),
                **self._numbers_dict(self.params_values),
            }

            nth_order_solution_expr = solution.subs(ics_dict).subs(notnumbers_dict)

            if (tuple(nth_order_solution_expr)) not in self.__class__._saved_numerized:
                self.__class__._saved_numerized[tuple(nth_order_solution_expr)] = (
                    lambdify(
                        [self.ivar]
                        + list(numbers_dict.keys())
                        + list(ics_symbols_dict.keys()),
                        nth_order_solution_expr.n(),
                    )
                )

            t_0 = time.time()

            nth_order_solution_fun = self.__class__._saved_numerized[
                tuple(nth_order_solution_expr)
            ]

            # return nth_order_solution_fun(ivar)
            solution = TimeDataFrame(
                data={
                    dvar: data[0]
                    for dvar, data in zip(
                        self.dvars,
                        nth_order_solution_fun(
                            ivar,
                            **{str(key): val for key, val in numbers_dict.items()},
                            **{str(key): val for key, val in ics_symbols_dict.items()},
                        ),
                    )
                    # for dvar, data in zip(self.dvars, nth_order_solution_fun(ivar))
                },
                index=ivar,
            )

            for dvar in self.dvars:
                solution[dvar.diff(self.ivar, 1)] = solution[dvar].gradient()

            for dvar in self.dvars:
                solution[dvar.diff(self.ivar, 2)] = solution[
                    dvar.diff(self.ivar, 1)
                ].gradient()
            t_e = time.time()

            t_d = t_e - t_0

            solution.index.name = self.ivar

            print("A_time" * 20, t_d)

            solution = solution.apply(np.real)

            solution._set_comp_time(t_d)

            print("B_time" * 20, solution._get_comp_time(t_d))

            return solution

    def _format_solution(self, dvars, solution, dict=False, equation=False):

        if equation:
            solution = Matrix([Eq(lhs, rhs) for lhs, rhs in zip(dvars, list(solution))])

        if dict:
            solution = {lhs: rhs for lhs, rhs in zip(dvars, list(solution))}

        return solution

    def approximation_function(self, order, max_order=1):

        dvars_no = len(self.dvars)

        t_list = self.t_list

        # Function('X')(t0(t),t1(t),t2(t)).diff(t,2).subs({t0(t).diff(t,1):1,t1(t).diff(t,1):eps,t2(t).diff(t,1):eps**2}).doit().expand().subs(Function('X')(t0(t),t1(t),t2(t)) ,Function('X0')(t0(t),t1(t),t2(t)) + Function('X1')(t0(t),t1(t),t2(t))*eps + Function('X2')(t0(t),t1(t),t2(t))  *eps**2 ).doit().expand().coeff(eps)

        # print(t_list)

        return Matrix(
            [
                sym.Function("Y_{" + str(dvar_no) + str(order) + "}")(*t_list)
                for dvar_no, dvar in enumerate(self.dvars)
            ]
        )

    def predicted_solution(self, order=1, dict=False, equation=False):

        dvars_no = len(self.dvars)

        solution = sum(
            (
                self.approximation_function(comp_ord, order) * self.eps**comp_ord
                for comp_ord in range(order + 1)
            ),
            sym.zeros(dvars_no, 1),
        )

        return self._format_solution(
            dvars=self.dvars, solution=solution, dict=dict, equation=equation
        )

    def eigenvalues(self):
        modes, eigs = (
            self.inertia_matrix().inv() * self.stiffness_matrix()
        ).diagonalize()

        return [eig for eig in eigs if eig != 0]

    def eigenvalue_approximation(
        self,
        order=1,
        eigenvalue=None,
        eigenvalue_no=None,
        eigenvalue_symbol=Symbol("omega"),
        series_coefficients=None,
    ):

        if not eigenvalue:
            eigenvalue = self.omega
        self.omega = eigenvalue

        if not series_coefficients:
            series_coefficients = [
                symbols("A_{}{}_1:{}".format(eom_no, no, order + 1))
                for eom_no, eom in enumerate(self.dvars)
                for no, coord in enumerate(self.dvars)
            ]

        # approximation={eigenvalue:eigenvalue_symbol**2 - sum(coeff**(eps_ord+1)  for eps_ord,coeff  in enumerate(series_coefficients) )}
        approximations = [
            eigenvalue**2
            - sum(
                [
                    coeff * self.eps ** (eps_ord + 1)
                    for eps_ord, coeff in enumerate(coeffs)
                ]
            )
            for coeffs in series_coefficients
        ]
        return approximations

    def eoms_approximation(self, order=3, odes_system=None):

        if not odes_system:
            odes_system = self.odes_system

        stiffness_mat = self.stiffness_matrix()

        eps_stiffness_mat = Matrix(
            *stiffness_mat.shape,
            [
                stiff_comp * approx
                for stiff_comp, approx in zip(
                    stiffness_mat, self.eigenvalue_approximation(order=order)
                )
            ],
        )

        first_ord_subs = {
            t_i.diff(self.ivar): self.eps**t_ord
            for t_ord, t_i in enumerate(self.t_list)
        }
        sec_ord_subs = {
            t_i.diff(self.ivar, 2): 0 for t_ord, t_i in enumerate(self.t_list)
        }

        odes_system = odes_system

        eoms_approximated = (
            Matrix(odes_system)
            .subs(self.predicted_solution(order, dict=True))
            .subs(sec_ord_subs)
            .doit()
            .subs(first_ord_subs)
        )

        t_fun_subs_dict = {
            expr_tmp: expr_tmp.subs(self.ivar, self.t_list[0])
            for expr_tmp in (
                eoms_approximated.atoms(Function)
                - {*self.t_list}
                - {
                    *sum(
                        [
                            list(self.approximation_function(ord_tmp, order))
                            for ord_tmp in range(order + 1)
                        ],
                        [],
                    )
                }
            )
        }

        return eoms_approximated.subs(t_fun_subs_dict)

    def eoms_approximation_list(self, max_order=3, min_order=0, odes_system=None):

        eoms_approximated = self.eoms_approximation(
            order=max_order, odes_system=odes_system
        ).expand()

        return [
            eoms_approximated.applyfunc(
                lambda obj: obj.diff(self.eps, order) / factorial(order)
            )
            .subs(self.eps, 0)
            .doit()
            for order in range(min_order, max_order + 1)
        ]

    def nth_eoms_approximation(self, order=3):

        eoms_approximated = self.eoms_approximation_list(
            max_order=order, min_order=order
        )

        return eoms_approximated[0]

    def _determine_secular_terms(self, zeroth_approx, order, ivar=None):

        #         print('zeroth approx')
        #         display(*zeroth_approx)

        if not ivar:
            ivar = self.t_list[0]

        #         print('===')
        sol_zeroth = self._find_nth_solution(
            zeroth_approx, order=0, secular_comps={}, ivar=ivar
        )

        #         print('+++ sol of zeroth')
        #         display('+'*100,sol_zeroth,'+'*100)
        #         print('+++ sol of zeroth')
        #         #eig_min=sqrt(self.eigenvalues()[0])

        #         #secular_comps = {sin(eig_min*self.ivar),cos(eig_min*self.ivar)}
        secular_comps = sum(sol_zeroth).atoms(sin, cos)

        self.int_const |= {
            eq.coeff(comp) for eq in sol_zeroth for comp in secular_comps
        }

        #         display('+'*100,secular_comps,'+'*100)

        return secular_comps

    def _nth_order_solution_with_C(self, order=1):

        t_list = self.t_list

        eoms_list = self.eoms_approximation_list(max_order=order)

        stiffness_mat = self.stiffness_matrix()

        eps_stiffness_mat = Matrix(
            *stiffness_mat.shape,
            [
                stiff_comp * approx
                for stiff_comp, approx in zip(
                    stiffness_mat, self.eigenvalue_approximation(order=order)
                )
            ],
        )

        secular_components = {
            comp: 0
            for comp in self._determine_secular_terms(
                zeroth_approx=eoms_list[0], order=order
            )
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
                dict=True,
                ivar=self.t_list[0],
            )
            #            display(nth_solution)
            approx_dict.update(nth_solution)

        sec_eq = Matrix(list(self.secular_eq - {S.Zero}))
        #         display(sec_eq)
        #         display(self.int_const)

        #         const_sol=FirstOrderODE((sec_eq),
        #                         self.t_list[1],
        #                         dvars=Matrix(list(self.int_const))
        #                         ).general_solution()
        #         display(const_sol)

        #         print('Const const')
        #         print(FirstOrderODE._const_list)

        t_back_subs = {
            t_i: self.ivar * self.eps**t_ord for t_ord, t_i in enumerate(self.t_list)
        }

        #         display(t_back_subs)

        general_form_dict = {
            var: eqn.subs(approx_dict)
            for var, eqn in self.predicted_solution(order=order, dict=True).items()
        }

        #         print('ics')
        #         display(*list({var: eqn.subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
        #         display(*list({var: eqn.diff(self.ivar).subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))

        #         display(const_vals)

        return general_form_dict

    def _ic_from_sol(self, order=1, formula=True):

        general_sol_obj = self.nth_order_solution(order=order)
        general_form_dict = general_sol_obj.as_dict()

        #         display(general_form_dict)

        ics_eqns = list(general_sol_obj.rhs.subs({self.ivar: 0})) + list(
            general_sol_obj.rhs.diff(self.ivar).subs({self.ivar: 0})
        )

        eqns = Matrix([eq.expand() for eq in ics_eqns])

        #         display(eqns)

        const_from_eqn = [
            var
            for var in list(eqns.atoms(Symbol, Function))
            if var in FirstOrderODE._const_list
        ]

        # display(Matrix(ics_eqns).jacobian(const_from_eqn))

        #         print('const from eqns')
        #         display(const_from_eqn)

        # display(ics_eqns)

        #         print('ics params check')
        #         display(self.params_values)

        if formula:

            const_vals_lin = (
                Matrix([eq.expand() for eq in ics_eqns])
                .jacobian(const_from_eqn)
                .subs({const: 0 for const in const_from_eqn})
                .subs(self.eps, 0)
            )  # .subs(self.params_values)#*Matrix(self.ics)
            const_vals_zth = Matrix([eq.expand() for eq in ics_eqns]).subs(
                {const: 0 for const in const_from_eqn}
            )  # .subs(self.params_values)

        else:

            const_vals_lin = (
                Matrix([eq.expand() for eq in ics_eqns])
                .jacobian(const_from_eqn)
                .subs({const: 0 for const in const_from_eqn})
                .subs(self.params_values)
            )  # *Matrix(self.ics)
            const_vals_zth = (
                Matrix([eq.expand() for eq in ics_eqns])
                .subs({const: 0 for const in const_from_eqn})
                .subs(self.params_values)
            )

        #         print('+++++++++  eqns for ics +++++++++++ ')
        #         display(const_vals_lin)
        #         display(const_vals_zth)

        self.ics_symbols = symbols(f"D0:{len(self.ics)}")

        self._ics_formula = {
            const: val
            for const, val in zip(
                const_from_eqn,
                const_vals_lin.inv() * (-const_vals_zth + Matrix(self.ics_symbols)),
            )
        }

        return self._ics_formula

    def _nth_order_solution(self, order=1):

        general_form_dict = self._nth_order_solution_with_C(order=order)

        sec_eq = Matrix(list(self.secular_eq - {S.Zero}))
        #         display(sec_eq)
        #         display(self.int_const)

        const_sol = FirstOrderODE(
            (sec_eq), self.t_list[1], dvars=Matrix(list(self.int_const))
        ).solution()

        #         print('trial for SS')
        steady_sol = FirstOrderODE(
            (sec_eq), self.t_list[1], dvars=Matrix(list(self.int_const))
        ).steady_solution()

        self._const_sol = {
            coord: const_sol[coord] + steady_sol[coord] for coord in const_sol.keys()
        }
        #         display(const_sol)

        #         print('Const const')
        #         print(FirstOrderODE._const_list)

        t_back_subs = {
            t_i: self.ivar * self.eps**t_ord for t_ord, t_i in enumerate(self.t_list)
        }

        #         display(t_back_subs)

        general_form_dict = {
            key: eq.subs(const_sol).subs(t_back_subs)
            for key, eq in general_form_dict.items()
        }

        #         print('ics')
        #         display(*list({var: eqn.subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))
        #         display(*list({var: eqn.diff(self.ivar).subs({self.ivar:0}).subs(self.eps,0)  for var, eqn in general_form_dict.items()}.values()))

        #         display(self._ics_formula)
        #         print('+++++++++  eqns for ics +++++++++++ ')

        return AnalyticalSolution.from_dict(general_form_dict)

    def nth_order_solution(self, order=1):

        if not (tuple(self._eoms), order) in self.__class__._saved_solution:

            self.__class__._saved_solution[(tuple(self._eoms), order)] = (
                self._nth_order_solution(order=order)
            )

        return self.__class__._saved_solution[(tuple(self._eoms), order)]

    def zeroth_approximation(self, dict=False, equation=False):

        # eoms = Matrix(self.odes_system).subs(self.eps, 0)

        eoms = self.nth_eoms_approximation(0)

        #        display(eoms, self.approximation_function(order=0))

        solution = LinearODESolution(
            eoms,
            ivar=self.ivar,
            dvars=self.approximation_function(order=len(self.t_list)),
        ).solution()

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
            equation=equation,
        )

    def _find_nth_solution(
        self, eoms_nth, order, secular_comps, dict=False, equation=False, ivar=None
    ):

        if not ivar:
            ivar = self.ivar

        eoms = eoms_nth

        # print('='*100,'eoms_nth for solve')
        # display(eoms )
        # print('='*100,'eoms_nth for solve')

        eoms = (eoms.expand()).applyfunc(
            lambda eqn: ((TR10(TR8(TR10(eqn).expand()).expand()).expand()))
        )

        # print('='*100,'eoms_nth for solve')
        # display(eoms )
        # print('='*100,'eoms_nth for solve')

        self.secular_eq |= {
            row.coeff(comp) for row in eoms for comp in secular_comps
        }  # test

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
            eoms, ivar=ivar, dvars=self.approximation_function(order=order)
        ).solution()

        #         print('=eoms and its linear sol'*100)
        #         display(eoms)
        #         display(solution)

        #         print('=eoms and its sol'*100)

        return self._format_solution(
            dvars=self.approximation_function(order=order),
            solution=solution,
            dict=dict,
            equation=equation,
        )

    def nth_approximation(self, order=1, dict=False, equation=False):

        #         eoms = Matrix(self.odes_system).subs(
        #             self.predicted_solution(2, dict=True)).diff(self.eps,order).subs(
        #                 self.eps, 0).expand()

        eoms = self.nth_eoms_approximation(order)

        zeroth_approximation = self.zeroth_approximation(dict=True, order=order)

        secular_comps = sum(zeroth_approximation.values()).atoms(cos, sin)

        approx = [zeroth_approximation] + [
            self.nth_approximation(order=i, dict=True) for i in range(1, order)
        ]

        approx_dict = type({})(ChainMap(*approx))

        #        display('abc',approx_dict)

        return self._find_nth_solution(
            eoms.subs(approx_dict), order, secular_comps, dict=dict, equation=equation
        )

    def first_approximation(self, dict=False, equation=False):
        return self.nth_approximation(order=1, dict=dict, equation=equation)

    def numerized(self, params_values={}, **kwargs):

        print("params values")
        print(params_values)

        if params_values == {}:
            return copy.copy(self)
        else:

            if "ics" in params_values:
                ics_list = params_values["ics"]
            else:
                ics_list = self.ics

            new_system = self.__class__(
                odes_system=self.governing_equations,
                ivar=self.ivar,
                dvars=self.dvars,
                ics=ics_list,
                eps=self.eps,
                omega=self.omega,
                order=self._order,
                t_span=[],
                params=[],
                params_values={**self.params_values, **params_values},
                ic_point=ics_list,
                equation_type=None,
                label=self._label,
            )

            new_system._stored_solution = copy.copy(self._stored_solution)
            new_system._saved_solution = copy.copy(self._saved_solution)

            return new_system

    def compute_solution(
        self, t_span=None, ic_list=None, t_eval=None, params_values=None, method="RK45"
    ):

        #         print('compute_solution')
        #         display(params_values)
        #         display(ic_list)

        if ic_list:
            print("class ics has been taken")
            self.ics = ic_list

        return self.__call__(
            ivar=t_span, order=self._order, params_values=params_values
        )
