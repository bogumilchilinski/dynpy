from ..mechanics import ReportComponent as BaseReportComponent
from ..mechanics import *


class ReportComponent(BaseReportComponent):

    @property
    def reported_object(self):

        from ....solvers.linear import ODESystem

        if isinstance(self._reported_object, ODESystem):
            return self._reported_object
        else:
            return self._reported_object._ode_system

    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object = obj

    @property
    def _system(self):
        print(
            "Kod do poprawienia #################### bo stosujesz starą zmienną self._system"
        )
        print("zamień se self._system na self.reported_object")

        return self.reported_object


class TitlePageComponent(Environment):

    latex_name = "titlepage"

    def __init__(
        self, system=None, options=None, arguments=None, start_arguments=None, **kwargs
    ):
        r"""
        Args
        ----
        options: str or list or  `~.Options`
            Options to be added to the ``\begin`` command
        arguments: str or list or `~.Arguments`
            Arguments to be added to the ``\begin`` command
        start_arguments: str or list or `~.Arguments`
            Arguments to be added before the options
        """

        self.system = system
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.system is not None:

            system = self.system

            self.append(NoEscape("\centering"))

            self.append(NoEscape("\\Huge DRGANIA MECHANICZNE \n \n"))

            self.append(Command("vspace", arguments="1cm"))

            if len(system.dvars) == 1:
                dof_str = "JEDNYM STOPNIU SWOBODY"
            else:
                dof_str = "WIELU STOPNIACH SWOBODY"

            if (
                system._dissipative_potential == 0
                or system._dissipative_potential is None
            ):
                damping_str = "NIETŁUMIONE"
            else:
                damping_str = "TŁUMIONE"

            self.append(NoEscape(f"\\Large {damping_str} UKŁADY O {dof_str} \n \n"))

            self.append(Command("vspace", arguments="1cm"))

            self.append(NoEscape(f"{system._label} \n \n"))

            self.append(Command("vspace", arguments="1cm"))

            self.append(Command("MyAuthor"))
            self.append(NoEscape(f"\\par"))
            self.append(Command("vspace", arguments="1cm"))

            # self.append(Command('vspace',arguments='1cm'))
            # self.append(NoEscape(f'\\protect\\par'))
            # self.append(NewLine())
            self.append(Command("MyDate"))


class ODESystemComponent(ReportComponent):

    title = "Differential equations"

    @property
    def header_text(self):

        return "The investigated system is described by differential equations being as follows:"

    @property
    def footer_text(self):

        return (
            "To solve the problem serve several methods depending on the equation type."
        )

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        for no, eq in enumerate(system):
            #             display(SympyFormula(Eq(system.lhs[no],system.rhs[no])))
            equation = Eq(system.lhs[no], system.rhs[no])
            display(SympyFormula(equation))
            display(SymbolsDescription(expr=equation))

        display(ReportText(self.footer_text))


class ODESystemCodeComponent(ReportComponent):

    title = "Differential equations"

    @property
    def header_text(self):

        return "The code that enables to obtain differential equations is as follows:"

    @property
    def footer_text(self):

        return "Presented method can change a type of the equation representation."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        display(
            Markdown(
                "\t from dynpy.solvers.linear import ODESystem  \n \t ode=ODESystem([]) \n \t display(ode)"
            )
        )

        display(ReportText(self.footer_text))


class VariablesComponent(ReportComponent):

    title = "System's variables"

    @property
    def header_text(self):

        return "As dynamic systems's behaviour is described by an ordinary differential equation, the variables of the equation are as follows:"

    @property
    def footer_text(self):

        return "The variables allow to study and analyse the system's dynamics varying with time."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        display(SympyFormula(system.ivar))
        display(SympyFormula(system.dvars[0]))
        display(SymbolsDescription(expr=system.ivar))

        display(ReportText(self.footer_text))


class GoverningEquationComponent(ReportComponent):

    title = "Equation of motion"

    @property
    def header_text(self):
        system = self.reported_object

        return f"The equation of motion was derived based on physical principles and governing laws regarding the considered system. The equation is described by the following expression ({AutoMarker(Eq(system.odes[0],0))})"

    @property
    def footer_text(self):
        return f"The determined equation constitutes a mathematical description of the dynamic properties of the system. Further analysis allows for an effective study on the modelled object's operation and determination of its mechanical parameters."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        for ode_expr in system.odes:
            display(SympyFormula(Eq(ode_expr, 0)))

        display(ReportText(self.footer_text))


class ApproximatedNonlinearGoverningEquationComponent(ReportComponent):

    title = "Approximated equation of motion"

    @property
    def header_text(self):

        return "The fully nonlinear form of the governing equation can be approximated by Taylor series expansion and then the system's motion equation is represented by the following:"

    @property
    def footer_text(self):

        return "Solving the nonlinear type of problems requires an involvement of dedicated methods of nonlinear theory of vibration."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        display(SympyFormula(Eq(system.approximated(3).simplify().expand(), 0)))

        display(ReportText(self.footer_text))


class FundamentalMatrixComponent(ReportComponent):

    title = "Determination of fundamental matrix"

    @property
    def header_text(self):

        return "The system physical properties of mass and stiffness are the basis to formulate the fundamental matrix of the system:"

    @property
    def body_text(self):

        system = self.reported_object

        Delta = Symbol("Delta")

        return f"The characteristic polynomial ${latex(Delta)}$ of the considered system derived based on its fundamental matrix is presented in {AutoMarker(Eq(Delta,system.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))}:"

    @property
    def footer_text(self):

        return "The system fundamental matrix and its characteristic equation allows for a determination of the system unsteady state and its eigenfrequencies related to free vibration."

    def append_elements(self):

        system = self.reported_object

        Delta = Symbol("Delta")

        display(ReportText(self.header_text))

        display(
            (SympyFormula(Eq(Symbol("A"), system.fundamental_matrix(), evaluate=False)))
        )

        display(ReportText(self.body_text))

        display(
            (
                SympyFormula(
                    Eq(
                        Delta,
                        system.fundamental_matrix().det().expand().simplify().expand(),
                        evaluate=False,
                    )
                )
            )
        )

        display(ReportText(self.footer_text))


class ODECharecteristicPolynomialComponent(ReportComponent):

    title = "Characteristic polynomial"

    @property
    def header_text(self):
        return "Characterisitc pylynomial for the considered system of equations of motion is the following:"

    @property
    def footer_text(self):
        return "Calculation of roots of this polynomial allows to find natural frequencies of the system."

    def append_elements(self):
        system = self.reported_object
        delta = Symbol("Delta")

        display(ReportText(self.header_text))
        polynomial = system.char_polynomial()
        display(SympyFormula(Eq(delta, polynomial)))

        display(ReportText(self.footer_text))


class ODECharecteristicPolynomialCodeComponent(ReportComponent):

    title = "Characteristic polynomial"

    @property
    def header_text(self):
        return "Presented piece of code generates the characteristic pylynomial of the desired system:"

    @property
    def footer_text(self):
        return "The following code is based on dynpy mechanical module."

    def append_elements(self):
        system = self.reported_object

        display(ReportText(self.header_text))

        display(
            Markdown(
                f"""

    from dynpy.models.mechanics import *
    from sympy import *
    from dynpy.solvers.linear import*

    dyn_sys = {system.__class__.__name__}()

    display( Eq(Symbol('Delta'),
                     dyn_sys._ode_system.char_polynomial()) )


"""
            )
        )

        display(ReportText(self.footer_text))


class ZerothOrderApproximatedEqComponent(ReportComponent):

    title = "Zeroth-order approximation"

    @property
    def header_text(self):

        from ....solvers.nonlinear import MultiTimeScaleSolution

        base_system = self.reported_object

        zeroth_ord_eq = Symbol("Omega")  # system.nth_eoms_approximation(0)
        eps = Symbol("\\varepsilon")

        return f"The ordering and separation of approximation equation in terms of the power of a small parameter ${eps}$ leads to obtaining a recursive sequence of linear equations of motion leading to the solution of the nonlinear equation. The zeroth-order approximate linear equation is given in the following paragraph where the dependency was assumed for the zeroth-order solution of the time variable:"

    @property
    def middle_text(self):

        from ....solvers.nonlinear import MultiTimeScaleSolution

        base_system = self.reported_object
        system = base_system.set_order(0)

        zeroth_ord_approx_eq = Symbol("Omega")

        return f"Or transformed to the second order ODE:"

    @property
    def footer_text(self):

        from ....solvers.nonlinear import MultiTimeScaleSolution

        base_system = self.reported_object
        system = system = base_system.set_order(0)

        t_list = base_system.t_list

        return f"Since ${latex(t_list[0])}$ and ${latex(t_list[1])}$ are treated as independent, the differential equation becomes a partial differential equation for a function of two variables ${latex(t_list[0])}$ and ${latex(t_list[1])}$."

    def append_elements(self):

        from ....solvers.nonlinear import MultiTimeScaleSolution

        base_system = self.reported_object
        system = base_system.set_order(0)

        t_list = system.t_list
        zeroth_ord_approx = system.eoms_approximation_list()[0]
        zeroth_ord_approx_eq = Eq(
            (zeroth_ord_approx.lhs - zeroth_ord_approx.rhs)[0], 0, evaluate=False
        )

        display(ReportText(self.header_text))

        for no, eq in enumerate(zeroth_ord_approx):

            eq_to_check = Eq(zeroth_ord_approx.lhs[no], zeroth_ord_approx.rhs[no])

            if isinstance(eq_to_check, Eq):
                display(SympyFormula(eq_to_check))

        display(ReportText(self.footer_text))


class FirstOrderApproximatedEqComponent(ReportComponent):

    title = "First-order approximation"

    @property
    def header_text(self):

        system = self.reported_object.set_order(1)
        first_ord_eq = system.eoms_approximation_list()[1]

        return f"The next component of a recursive sequence of linear equations of motion leading to the solution of the considered nonlinear one is given as follows:"  # {AutoMarker(first_ord_eq)}:"

    @property
    def middle_text(self):

        system = self.reported_object.set_order(1)
        first_ord_approx_eq = Eq(
            system.eoms_approximation_list()[0].lhs
            - system.eoms_approximation_list()[0].rhs,
            0,
        )

        return (
            f"Or transformed to the second order ODE {AutoMarker(first_ord_approx_eq)}:"
        )

    @property
    def footer_text(self):

        system = self.reported_object.set_order(1)
        t_list = system.t_list

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of ${latex(t_list[1])}$."

    def append_elements(self):

        system = self.reported_object.set_order(1)

        t_list = system.t_list
        first_ord_approx = system.eoms_approximation_list()[1]
        first_ord_approx_eq = Eq(
            (first_ord_approx.lhs[0] - first_ord_approx.rhs[0]), 0, evaluate=False
        )

        display(ReportText(self.header_text))

        for no, eq in enumerate(first_ord_approx):

            eq_to_check = Eq(first_ord_approx.lhs[no], first_ord_approx.rhs[no])

            if isinstance(eq_to_check, Eq):
                display(SympyFormula(eq_to_check))

        display(ReportText(self.footer_text))


class MSMCalculationsOrderComponent(ReportComponent):

    title = "Calculation algorithm for MSM"

    def append_elements(self):

        system = ode = self.reported_object

        display(
            ReportText(
                "First step is to establish the approximation function. This is done using the method"
            )
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation \node.approximation_function(ord)"
            )
        )

        display(
            ReportText("Next step is prediction of solution using following method")
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.predicted_solution(ord)"
            )
        )
        display(
            Markdown(
                "This method uses `approximation_function` method and multiplies it by small parameter to obtain the desired predicted solution."
            )
        )

        display(
            ReportText(
                "Later the derivative of predicted solution should be calculated. This is done using:"
            )
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.derivative_without_scale(ord)"
            )
        )
        display(
            Markdown(
                "This method uses `predicted_solution_without_scales` method and calculates a derivatives."
            )
        )

        display(
            ReportText("Derivatives of each approximation function is calculated by:")
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.part_derivative(ord)"
            )
        )
        display(
            ReportText("This method calculates derivatives of approximation function.")
        )

        display(
            ReportText("After that time scales and their derivatives are calculated.")
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node._scales_formula(ord)"
            )
        )

        display(ReportText("Substitution of time scales is done using the method:"))
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.part_derivative_subs(ord)"
            )
        )

        display(
            ReportText(
                "Next step is the substitution of approximation function to predicted solution."
            )
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.eoms_approximation(ord)"
            )
        )
        display(
            ReportText(
                f"This method returns object of {type(ode.eoms_approximation())} type that is used to create list of approximation odes."
            )
        )

        display(
            ReportText(
                f"Odes approximation allows to find list of new ode systems, it's handled by the following method:"
            )
        )
        display(
            ObjectCode(
                "ord = 2 # given order of approximation  \node.eoms_approximation_list(ord)"
            )
        )

        display(
            ReportText("General solution of the resulting equation can be found using")
        )
        display(ObjectCode("ord = 2 # given order of approximation  \node.general_sol"))
        display(Markdown("This method is using `eoms_approx_with_const` method."))

        display(ReportText("Resulting equation is solved"))
        display(ObjectCode("ord = 2 # given order of approximation  \node.solution"))


class PredictedSolutionComponent(ReportComponent):

    title = "Predicted solution equation"

    @property
    def header_text(self):

        return "Solving the given equation begins with assuming predicted solutions, which are determined based on the chosen order of approximation."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        for ord in range(system.order + 1):

            display(ReportText(f"\n Order of prediction: {ord}"))
            predicted_sol = Eq(
                system.dvars[0], system.predicted_solution(ord).rhs[0], evaluate=False
            )
            display(SympyFormula(predicted_sol))

        display(SymbolsDescription(expr=predicted_sol))


class RawPredictedSolutionComponent(ReportComponent):

    title = "Predicted solution equation in ODE System form"

    @property
    def header_text(self):

        return f"Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the consecutive approximations (depending on the accuracy assumed) have the following form:"

    @property
    def footer_text(self):

        return "Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth-order approximation {approx_fun} has the following form:"

    def append_elements(self):

        system = self.reported_object

        #         t_list = system.t_list

        display(ReportText(self.header_text))

        for ord in range(system.order + 1):

            display(ReportText(f"Order of prediction: {ord}"))
            display(SympyFormula(system.predicted_solution(ord)))


class DetailsOfPredictedSolutionComponent(ReportComponent):

    title = "Details of predicted solution with calls"

    @property
    def header_text(self):

        return f"Rough (first) approximation in time domain has the following form:\n"

    def append_elements(self):

        system = self.reported_object
        ivar = system.ivar

        #         t_list = system.t_list

        display(ReportText(self.header_text))
        display(
            SympyFormula(system.predicted_solution_without_scales().as_eq_list()[0])
        )
        display(SympyFormula(system.predicted_solution().as_eq_list()[0]))

        display(
            ReportText(
                "Comparison of Approximation with and without multiple scales has a following form \n"
            )
        )

        for ord in range(system.order + 1):

            display(ReportText(f"Order of prediction: {ord}"))

            display(
                SympyFormula(
                    Eq(
                        system.dvars[0],
                        system.predicted_solution_without_scales(ord).rhs[0],
                        evaluate=False,
                    )
                )
            )
            display(
                SympyFormula(
                    Eq(
                        system.dvars[0],
                        system.predicted_solution(ord).rhs[0],
                        evaluate=False,
                    )
                )
            )

            display(
                ObjectCode(
                    f"ord = {ord} # given order of approximation  \node.predicted_solution(ord)"
                )
            )

        display(
            ReportText(
                "Finally, formulas of predicted solution and its derivatives are as follows:"
            )
        )
        for degree in range(0, 3):
            display(
                SympyFormula(system.derivative_without_scale(degree).as_eq_list()[0])
            )

        display(ReportText("Also multiple time scales are introduced:"))

        predicted_sol_scales = system.predicted_solution(2).as_eq_list()[0]

        display(SympyFormula(predicted_sol_scales))
        display(
            SympyFormula(
                Eq(
                    predicted_sol_scales.lhs.diff(ivar),
                    predicted_sol_scales.rhs.diff(ivar),
                ).expand()
            )
        )
        display(
            SympyFormula(
                Eq(
                    predicted_sol_scales.lhs.diff(ivar, 2),
                    predicted_sol_scales.rhs.diff(ivar, 2),
                ).expand()
            )
        )

        display(ReportText("Time scales are given by the following functions:"))

        for tm_scale in system._scales_formula.as_eq_list():
            display(SympyFormula(tm_scale))

        display(
            ReportText(
                "First derivatives of the scales (needed for further calculations) are as follows:"
            )
        )

        for tm_scale in system._scales_formula.as_eq_list():
            display(SympyFormula(Eq(tm_scale.lhs.diff(ivar), tm_scale.rhs.diff(ivar))))

        display(ReportText("All second derivatives of the scales are zero:"))

        for tm_scale in system._scales_formula.as_eq_list():
            display(
                SympyFormula(Eq(tm_scale.lhs.diff(ivar, 2), tm_scale.rhs.diff(ivar, 2)))
            )


class GeneralSolutionComponent(ReportComponent):

    title = "Zeroth-order approximated solution"

    @property
    def header_text(self):

        system = self.reported_object

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of t. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

    @property
    def footer_text(self):

        return f"In order to determine the functions \(\operatorname{C_1}\left(t_{1}\right),\operatorname{C_2}\left(t_{1}\right)\) and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        for no, eq in enumerate(system._general_solution):
            display(
                SympyFormula(
                    Eq(
                        system._general_solution.lhs[no],
                        system._general_solution.rhs[no],
                    )
                )
            )


class FLODESolutionComponent(ReportComponent):

    title = "Solution test component"

    @property
    def header_text(self):

        system = self.reported_object

        return f"Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

    @property
    def footer_text(self):

        return "In order to determine the functions $\\operatorname{C_1}(t_{1})$, $\\operatorname{C_2}(t_{1})$ and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):

        system = self.reported_object

        #         display(ReportText(  self.header_text   ))

        from dynpy.solvers.linear import (
            FirstOrderLinearODESystem,
            FirstOrderLinearODESystemWithHarmonics,
        )

        display(
            FirstOrderLinearODESystemWithHarmonics.from_ode_system(
                self.reported_object
            ).solution
        )

        display(
            FirstOrderLinearODESystem.from_ode_system(self.reported_object).solution
        )


class ZerothOrderSolutionComponent(ReportComponent):

    title = "Zeroth-order approximated solution"

    @property
    def header_text(self):
        t_list = self.reported_object.t_list

        return f"Zeroth-order approximated solution is"

    @property
    def footer_text(self):

        return "Zeroth-order approximated solution has following consts"

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        ode_0th_approx = system.eoms_approximation_list()[0]
        display(
            SympyFormula(
                Eq(ode_0th_approx.solution.lhs[0], ode_0th_approx.solution.rhs[0])
            )
        )
        #         display(SympyFormula( Eq(ode_0th_approx.solution.lhs[1], ode_0th_approx.solution.rhs[1]  )))
        display(ReportText("Secular functions are as follows:"))

        for sec_fun in ode_0th_approx._secular_funcs:
            display(SympyFormula(sec_fun))


class SecularTermsEquationsComponent(ReportComponent):

    title = "First-order approximated solution and secular terms"

    @property
    def header_text(self):
        t_list = self.reported_object.t_list

        return f"Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

    @property
    def footer_text(self):

        return "In order to determine the functions $\\operatorname{C_1}(t_{1})$, $\\operatorname{C_2}(t_{1})$, the aforementioned first-order approximate equation has to be considered."

    def append_elements(self):

        system = self.reported_object
        system.order = 1

        display(ReportText(self.header_text))

        system.general_solution
        display(system.secular_eq)
        sec_eqns = (
            system.secular_eq[system.eps].as_first_ode_linear_system().lhs
            - system.secular_eq[system.eps].as_first_ode_linear_system().rhs
        )

        for no, eq in enumerate(sec_eqns):

            display(SympyFormula(Eq(eq, 0)))

        display(SymbolsDescription(expr=system.secular_eq[system.eps].dvars))

        display(ReportText(self.footer_text))


class GoverningEqnsWithConstFuncsComponent(ReportComponent):

    title = "Governing equations with constant functions"

    @property
    def header_text(self):
        return "The rearranged governing equations contaning the functions of $\\operatorname{C_1}(t_{1})$, $\\operatorname{C_2}(t_{1})$ are prepared in order to furhter spot and remove the secular terms from the equations. These are as follows:"

    @property
    def footer_text(self):
        return "At this stage it needs to be checked whether all the expected functions contained."

    def append_elements(self):

        system = self.reported_object
        approx_no = 1

        display(ReportText(self.header_text))

        for no, eq in enumerate(system.eoms_approximation_list(approx_no)):
            display(ReportText(f"The equation for the approximation number {no} is:"))
            display(
                SympyFormula(
                    Eq(
                        system.eoms_approx_with_const()[no].lhs,
                        system.eoms_approx_with_const()[no].rhs,
                    )
                )
            )

        display(ReportText(self.footer_text))


class SecularFuncsComponent(ReportComponent):

    title = "Secular functions"

    @property
    def header_text(self):
        return "Text"

    @property
    def footer_text(self):
        return "Text"

    def append_elements(self):

        system = self.reported_object
        approx_no = 1

        display(ReportText(self.header_text))

        for no, eq in enumerate(system.eoms_approximation_list(approx_no)):
            display(
                ReportText(
                    f"Set of secular functions for the approximated solution {no} is:"
                )
            )
            display(
                SympyFormula(
                    Eq(
                        Symbol(f"Set_{no}"),
                        system.eoms_approx_with_const()[no]._secular_funcs,
                        evaluate=False,
                    )
                )
            )

        display(ReportText(self.footer_text))


class SecularConditionsComponent(ReportComponent):

    title = "Secular conditions"

    @property
    def header_text(self):
        return "Text"

    @property
    def footer_text(self):
        return "Text"

    def append_elements(self):

        system = self.reported_object
        approx_no = 1

        display(ReportText(self.header_text))

        for no, eq in enumerate(system.eoms_approximation_list(approx_no)):
            display(
                ReportText(f"Secular conditions for the approximated solution {no} is:")
            )
            display(
                SympyFormula(
                    Eq(
                        Symbol(f"C_{no}"),
                        system.eoms_approx_with_const()[no]._secular_conditions,
                        evaluate=False,
                    )
                )
            )

        display(ReportText(self.footer_text))


class SolutionWithConstFuncsComponent(ReportComponent):

    title = "Approximated solutions with constant functions"

    @property
    def header_text(self):
        return "In order to investigate the early stage of particular solutions, the following equations are analysed and chcecked in terms of whether they contain the expected functions:"

    @property
    def footer_text(self):
        return "Having the solutions with constants $\\operatorname{C_1}(t_{1})$, $\\operatorname{C_2}(t_{1})$ as expected, the next step is to solve the secular equations."

    def append_elements(self):

        system = self.reported_object
        approx_no = 1

        display(ReportText(self.header_text))

        for no, eq in enumerate(system.eoms_approximation_list(approx_no)):
            display(ReportText(f"Solution for {no} approximation is then:"))
            display(
                SympyFormula(
                    Eq(
                        system.eoms_approx_with_const()[no].solution.lhs,
                        system.eoms_approx_with_const()[no].solution.rhs,
                    )
                )
            )
        #             display(system.eoms_approx_with_const()[no].solution)

        display(ReportText(self.footer_text))


class HomEquationComponent(ReportComponent):

    title = "Homogenous equation"

    @property
    def header_text(self):

        return "Homogenous equation is the equation in which all the external forces and excitations are assumed to be equal to zero."

    @property
    def footer_text(self):

        return "Homogenous equation is used to find the general solution of differential equation which allows to analyse free vibrations of the system."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        display(SympyFormula(system._hom_equation()))
        display(
            ReportText(
                "Equations presented above correspond to each generalized coordinate:"
            )
        )
        display(SympyFormula(system.dvars))

        display(ReportText(self.footer_text))


class HomEquationCodeComponent(ReportComponent):

    title = "Code responsible for calculating homogenous equation"

    @property
    def header_text(self):

        return "Code line responsible for extracting homogenous equation from ODESystsem is following."

    @property
    def footer_text(self):

        return "Homogenous equation is used to find the general solution of differential equation which allows to analyse free vibrations of the system."

    def append_elements(self):

        system = self.reported_object

        display(ReportText(self.header_text))

        display(
            Markdown(
                """
    from sympy import *
    from dynpy.solvers.linear import *

    display(SympyFormula( system._hom_equation()))
"""
            )
        )

        display(SympyFormula(system._hom_equation()))
        display(
            ReportText(
                "Equations presented above correspond to each generalized coordinate:"
            )
        )
        display(
            ReportText(
                "Code responsible for displaying generalized coordinates is following:"
            )
        )
        display(
            Markdown(
                """
        display(SympyFormula(system.dvars))
"""
            )
        )
        display(SympyFormula(system.dvars))

        display(ReportText(self.footer_text))


class ODEInitComponent(ReportComponent):

    title = "Exemplary system"

    @property
    def header_text(self):

        return "Exemplary ODE system can be created using following code."

    @property
    def footer_text(self):

        return "Code presented above shows the elements which are neccessary to build ODESystem."

    def append_elements(self):

        system = self.reported_object

        if len(system.dvars) < 2:
            display(
                ReportText("Differential equation to initialize class is as follows:")
            )
        else:
            display(
                ReportText("Differential equations to initialize class are as follows:")
            )
        for n in system.odes:
            display(SympyFormula(n))

        display(
            ReportText(
                "Depenedent variables of the presented case has a following form:"
            )
        )
        display(SympyFormula(system.dvars))

        display(ReportText(self.header_text))

        code_list = python(system.lhs).split("\n")
        var_name = "my_eq"
        code = "\n\n\t".join(code_list[:-1]) + "\n\n\t" + var_name + code_list[-1][1:]
        dvars = str(system.dvars)
        order = str(system.ode_order)
        if len(system.dvars) < 2:
            display(
                Markdown(
                    """Code presented below describes all of the symbols needed to construct ordinary differential equation. The `my_eq`, `funcs` and `ode_order` variables store the differential equation, unknown functions and equation order."""
                )
            )
        else:
            display(
                Markdown(
                    """Code presented below describes all of the symbols needed to construct ordinary differential equations. The `my_eq`, `funcs` and `ode_order` variables store the differential equations, unknown functions and equation order."""
                )
            )
        display(
            Markdown(
                f"""

    from sympy import *
    from dynpy.solvers.linear import ODESystem

    {code}

    order = {str(system.ode_order)}
    funcs = {str(system.dvars)}

    odesys = ODESystem(odes = my_eq, dvars = funcs, ode_order = order)
    display(odesys)
"""
            )
        )
        if len(system.dvars) < 2:
            display(
                ReportText(
                    "After creating our differential equation it can be analysed using ODESystem class. Additionally one has to provide dvars and ode_order arguments. Dvars is the generalized coordinate of the system and ode_order is the order of considered differential equation."
                )
            )
        else:
            display(
                ReportText(
                    "After creating our differential equations they can be analysed using ODESystem class. Additionally one has to provide dvars and ode_order arguments. Dvars is the generalized coordinate of the system and ode_order is the order of considered differential equations."
                )
            )
        display(SympyFormula(system))

        display(ReportText(self.footer_text))


# class ODENonHomogeneousEquationComponent(ReportComponent):

#     title="Differential equations"
#     @property
#     def header_text(self):

#         return "By separating the equation we obtain the form of equation with time dependent variables on one hand side and free terms on the other:"


#     @property
#     def footer_text(self):

#         return "To solve the problem serve several methods depending on the equation type."

#     def append_elements(self):

#         system = self._system

#         display(ReportText(  self.header_text   ))

#         display(SympyFormula( Eq(system._hom_equation().lhs, system._free_component())))

#         display(ReportText(  self.footer_text   ))

# class ODENonHomogeneousEquationCodeComponent(ReportComponent):

#     title="Differential equations"
#     @property
#     def header_text(self):

#         return "By separating the equation we obtain the form of equation with time dependent variables on one hand side and free terms on the other:"


# <<<<<<< HEAD
# =======
# class ODENonHomogeneousEquationCodeComponent(ReportComponent):

#     title="Differential equations"
#     @property
#     def header_text(self):

#         return "By separating the equation we obtain the form of equation with time dependent variables on one hand side and free terms on the other:"


# >>>>>>> f02b360 (Mechanics module maintenance and development)
#     @property
#     def footer_text(self):

#         return "To solve the problem serve several methods depending on the equation type."

#     def append_elements(self):

#         system = self._system

#         display(ReportText(  self.header_text   ))
# <<<<<<< HEAD

#         display(Markdown
# ('''
#     from sympy import *
#     from dynpy.solvers.linear import ODESystem

#     display(SympyFormula(Eq(system._hom_equation().lhs, system._free_component())))
# =======
# ''')
#         display(Markdown
# ('''
#     system = self._system
#     display(SympyFormula( Eq(system._hom_equation().lhs, system._free_component())))
# >>>>>>> f02b360 (Mechanics module maintenance and development)
# '''))

#         display(SympyFormula( Eq(system._hom_equation().lhs, system._free_component())))

# <<<<<<< HEAD
#         display(ReportText(  self.footer_text   ))


class SeparableODEIntroComponent(ReportComponent):
    title = (
        "The basic knowledge necessary to solve an equation with separated variables."
    )

    def append_elements(self):
        system = self.reported_object
        ode_sys = system

        display(
            ReportText(
                f"One of the fundamental types of differential equations is an equation with separated variables. The process of solving this type of equation involves separating the terms containing $y$ and $x$, moving them to different sides of the equation, and then applying integration. After integrating, its important not to forget to add the constant of integration to one side. Below is an example problem illustrating the process of solving this type of equation."
            )
        )


class EquationDefinitionComponent(ReportComponent):
    title = "We consider the following differential equation:"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system

        y_sym = system._variable_name()[2]
        display(SympyFormula(Eq(ode_sys.lhs[0], ode_sys.rhs[0])))
        display(ReportText("Where:"))
        display(SympyFormula(y_sym))
        display(ReportText("Is our encountered function."))


class VariablesSeparationComponent(ReportComponent):

    title = "Separation of variables"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._get_dvars_symbols()[fun.diff(ivar)]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._get_dvars_symbols()[fun]
        gx = Symbol("g(x)")
        hx = Symbol("h(x)")

        # ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]

        # sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        display(
            ReportText(
                f"If possible, we perform the commonly known separation of variables $x$ and $y$. By separating variables, we mean separating terms containing the variable $x$ from those containing the variable $y$. To accomplish this separation, we can assign these terms to functions $h(x)$ and $g(x)$, as shown below:"
            )
        )
        display(SympyFormula(Eq(gx, g_fun_ivar)))
        display(SympyFormula(Eq(hx, h_fun_dvar)))
        display(
            ReportText(
                f"Apart from these elements, our equation also contains the term ${latex(fun_str)}$, which can also be written as:"
            )
        )
        display(SympyFormula(Eq(fun_str, d_var / d_ivar)))
        display(
            ReportText(
                f"After substituting this expression, it can be observed that in subsequent steps, multiplying the equation on both sides by ${latex(d_ivar)}$ enables integrating both sides with respect to different variables."
            )
        )


class SeparatedVariablesIntegrationComponent(ReportComponent):

    # title="Integration of separated variables"
    title = "Integration of separated variables"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar
        gx = Symbol("g(x)")
        hx = Symbol("h(x)")

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]
        y_sym = system._variable_name()[2]
        fun_str = system._variable_name()[3]

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        #         display(ReportText(''))

        display(
            SympyFormula(Eq(Symbol(latex(d_var / d_ivar)), g_fun_ivar * h_fun_dvar))
        )
        display(
            ReportText(
                f"""To do this, we multiply both sides by ${latex(d_ivar)}$ and move ${latex(fun)}$ to the opposite side of the equation:"""
            )
        )

        display(
            SympyFormula(
                Symbol(
                    latex(d_var / h_fun_dvar) + "=" + latex(g_fun_ivar) + latex(d_ivar)
                )
            )
        )
        display(
            ReportText(
                """Then, to "go back" to the original function and obtain the result in the form of the f $y$, we integrate both sides of the obtained equation. """
            )
        )

        display(
            ReportText(
                """Importantly, when integrating, it's essential to remember to add the constant of integration in any form, typically denoted as $C$ or $D$. """
            )
        )
        display(
            SympyFormula(
                Symbol(
                    latex(Integral(1 / h_fun_dvar, y_sym))
                    + "= "
                    + latex(Integral(g_fun_ivar, ivar))
                )
            )
        )

        #         const = (system.solution._spot_constant())[0]
        const = Symbol("D")
        display(ReportText("After integration, we obtain the following result:"))
        display(
            SympyFormula(
                Symbol(
                    latex(Integral(1 / h_fun_dvar, y_sym).doit())
                    + "="
                    + latex(Integral(g_fun_ivar, ivar).doit())
                    + "+"
                    + latex(const)
                )
            )
        )


#         display(ReportText(  self.footer_text   ))


class SolutionComponent(ReportComponent):

    title = "Simplification of the solution obtained after integrating the equations"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]

        ode_rhs = solve(ode_sys.lhs[0], ode_sys.dvars[0].diff())[0]

        sep_vars_dict = separatevars(ode_rhs.subs(fun, y_sym), (ivar, y_sym), dict=True)

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]

        display(ReportText("After rearrangment a solution has the following form:"))
        display(SympyFormula(system.general_solution.as_eq_list()[0]))

        display(
            ReportText(
                "The obtained expression represents the solution to a differential equation with separated variables. The constant of integration can be determined when initial conditions are provided in the problem, i.e., the value of the function for a specific $x$."
            )
        )


class LinearODEIntroComponent(ReportComponent):
    title = "The basic knowledge necessary to solve a differential equation"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]
        Px = Symbol("P(x)")
        Qx = Symbol("Q(x)")
        n = Symbol("n")
        display(
            ReportText(
                "One of the most fundamental types of differential equations is a linear equation with constant coefficients. A linear equation consists only of linear terms, which may include constant coefficients. Transitioning to a homogeneous equation is required - we move all terms involving the function $y$ and its derivatives, then equate them to zero, obtaining a homogeneous equation. Next, it is necessary to treat the constant as a function dependent on $x$. We calculate the derivatives of both sides and utilize the results in the initial equation. The result is the sum of solutions from both parts."
            )
        )


class LinearTransformation(ReportComponent):
    title = "Transition to a homogeneous equation"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        fun_str = system._variable_name()[3]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy

        y_sym = system._variable_name()[2]  # y(x)
        const = Symbol("D")

        display(
            ReportText(
                "We move all terms involving the function $y$ and its derivatives to one side and equate it to zero, obtaining a homogeneous equation."
            )
        )
        # display(SympyFormula(Eq(system._hom_equation().lhs[0],0)))


class LinearToSeparable(ReportComponent):

    title = "The solution to the resulting separated variable equation"

    def append_elements(self):

        system = self.reported_object
        fun = system.dvars[0]

        ode_sys = system._hom_equation()
        # sep_ode_report=ode_sys.report #mniej więcej coś takiego tylko coś pomyliłem na 100%

        display(
            ReportText(
                "The obtained homogeneous equation is a separated variable equation and takes the form:"
            )
        )
        display(SympyFormula(Eq(ode_sys.lhs[0], 0)))

        # print("#"*200)
        sep_ode_report = ode_sys.report

        # print(list(sep_ode_report))
        # print('#'*3)
        # print(list(sep_ode_report)[3:])

        for elem in list(sep_ode_report)[3:]:
            self.append(elem)


class VariationOfConstant(ReportComponent):
    title = "Variation of constant"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        ode_sys2 = ode_sys
        x = ivar
        szur = ode_sys2.general_solution
        #         display(szur)
        ode_syso = Eq(ode_sys2.lhs[0], ode_sys2.rhs[0])
        Cx = Function(szur._spot_constant()[0])(x)
        Cx
        c_prim = Function("C'")(x)
        szur2 = szur.subs(szur._spot_constant()[0], Cx)
        # display(szur2)
        szur3 = Eq(szur2.lhs.diff(x), szur2.rhs.diff(x))
        # display(szur3)
        szur3_ = szur3.subs(Cx.diff(x), c_prim)
        szur3_
        # To implement, change left hand side to y_prim variable and subsitute variables with object variables

        ##### Tu sie teoretycznie zaczyna variation of constant ########################################
        sol_cx = Eq(Symbol("y"), szur2.rhs[0])
        display(
            ReportText(
                "The next step is to treat the constant as a function dependent on $x$:"
            )
        )
        display(SympyFormula(sol_cx))
        sol_diff3 = Eq(Symbol("y'"), szur3_.rhs[0])
        display(
            ReportText(
                "Next, you need to calculate the derivatives of both sides of the equation."
            )
        )
        ## dodać równanie
        display(ReportText("Upon obtaining the equation in the following form:"))
        display(SympyFormula(sol_diff3))  # gdzie? gdzie co?
        display(
            ReportText(
                "Thanks to this operation, we can now substitute $y$ and $y'$ into the original equation."
            )
        )
        display(SympyFormula(ode_syso))
        lolo = ode_syso.subs(y.diff(ivar), sol_diff3.rhs).subs(y, sol_cx.rhs)
        display(ReportText("And obtain:"))
        display(SympyFormula(lolo))
        c_p_sol = solve(lolo, c_prim)
        cp_eq = Eq(c_prim, c_p_sol[0])
        display(
            ReportText(
                "Then, if everything has been calculated correctly, the values of $C(x)$ should simplify to a form from which we can easily determine $C'(x)$, obtaining:"
            )
        )
        display(SympyFormula(cp_eq))
        cx_wyn = Eq(Cx, cp_eq.rhs.diff(ivar))
        display(
            ReportText(
                "After integrating both sides, we obtain the value of the constant:"
            )
        )
        display(SympyFormula(cx_wyn))
        FINALOWY_WYNIK = sol_cx.subs(Cx, cx_wyn.rhs)
        display(
            ReportText(
                "To obtain the final result, we need to substitute the obtained constant into the first equation where we introduced the constant as a variable:"
            )
        )
        display(SympyFormula(FINALOWY_WYNIK))


class BernoulliODEIntroComponent(ReportComponent):
    title = "The basic knowledge necessary to solve a differential equation:"

    def append_elements(self):

        system = self.reported_object
        y_sym = system.dvars[0]
        dane = system._get_dvars_symbols()
        Px = Symbol("P(x)")
        Qx = Symbol("Q(x)")
        n = Symbol("n")
        ynpow = Symbol("y") ** Symbol("n")
        display(
            ReportText(
                "The Bernoulli differential equation is a special type of first-order ordinary differential equation. Its general form is:"
            )
        )
        bernoulli_form = Eq(y_sym.diff() + Px * y_sym, Qx * y_sym**n)
        display(SympyFormula(bernoulli_form.subs(dane)))

        display(
            ReportText(
                f"Where $y$ is the unknown function dependent on the variable $x$, $P(x)$ and $Q(x)$ are given functions dependent on $x$, and $n$ is a constant real number, different from 0 and 1. Characteristic of Bernoulli equations is that they are nonlinear due to the presence of the term ${latex(ynpow)}$ on one side of the equation. Solutions to such equations can be challenging to obtain in general, but for certain values of $n$, certain techniques can be applied to simplify the solution. In the case of $n=0$, the Bernoulli equation becomes linear, and for $n=1$, it can be transformed into the form of a linear ordinary differential equation of the first order. In other cases, special techniques such as variable substitution or order reduction may be necessary to obtain the general solution."
            )
        )


class EquationDefinitionComponent(ReportComponent):
    title = "We consider the following differential equation:"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        y_sym = system.dvars[0]

        dane = system._get_dvars_symbols()
        ode_sub = ode_sys.subs(dane)
        display(SympyFormula(Eq(ode_sub.lhs[0], ode_sub.rhs[0])))
        display(ReportText("Where:"))
        display(SympyFormula(y_sym))
        display(ReportText("Is our encountered function."))


class BernoulliTransformation(ReportComponent):
    title = "Determining the new function $z$ dependent on $x$ in order to transition to a linear equation."

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar
        ###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy
        dane = system._get_dvars_symbols()
        bernoulli_sub = system._bernoulli_substitution()
        # bern_power = system

        n = Symbol("n")
        z = Function("z")(ivar)

        display(
            ReportText(
                "To make the substitution, it is necessary to define $n$ by finding the highest power of $y$, where $n$ is the exponent of this function."
            )
        )
        display(SympyFormula(Eq(Symbol("n"), n.subs(bernoulli_sub))))
        display(
            ReportText(
                "Next, you need to make the substitution and transition to the new function $z$, using the following formula:"
            )
        )
        display(SympyFormula(Eq(Symbol("z"), Symbol("y") ** (S.One - Symbol("n")))))
        display(ReportText("In our case:"))
        display(
            SympyFormula(
                Eq(Symbol("z"), Symbol("y") ** (S.One - n.subs(bernoulli_sub)))
            )
        )
        # display(SympyFormula(Eq(Symbol('z'),z.subs(bernoulli_sub))))
        display(
            ReportText(
                "The next step is to determine the first derivative of the function $z$:"
            )
        )
        display(SympyFormula(Eq(Symbol("z'"), bernoulli_sub["z_diff"])))
        display(ReportText("After substituting the $n$ value :"))
        display(
            SympyFormula(
                Eq(Symbol("z'"), bernoulli_sub["z_diff"].subs(n, bernoulli_sub[n]))
            )
        )

        display(
            ReportText(
                "Next, we need to determine $y$ and $y'$ from the above equations."
            )
        )
        display(SympyFormula(Eq(Symbol("y"), fun.subs(bernoulli_sub))))
        display(SympyFormula(Eq(Symbol("y'"), fun.diff().subs(bernoulli_sub))))


class BernoulliLinearTransformation(ReportComponent):
    title = " Transition to a linear equation"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar
        ###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy
        dane = system._get_dvars_symbols()

        rownanie = system._as_reduced_ode()

        display(
            ReportText(
                "The last step in transitioning to a linear equation is substituting the previously calculated $y$ and $y'$ into the original equation:"
            )
        )
        display(SympyFormula(Eq(rownanie.lhs[0], rownanie.rhs[0])))
        display(ReportText("We compute the obtained equation in the following manner:"))

        sep_ode_report = rownanie.report

        for elem in list(sep_ode_report)[3:]:
            self.append(elem)


class MultiScaleProcedureComponent(ReportComponent):
    title = "Algorithm needed to use MSM solver"

    def append_elements(self):
        system = self.reported_object
        display(
            ReportText(
                "Step by step instruction: \n\n Step 1: \n\n Step 2: \n\n Step 3: \n\n Step 4: \n\n Step 5:"
            )
        )


class MSMComparison(ReportComponent):
    title = "Plot Comparison: Numerically Simulated Results vs. Analytical Solutions"

    def append_elements(self):
        from dynpy.solvers.nonlinear import MultiTimeScaleSolution

        system = self.reported_object
        t_span = np.linspace(0.0, 10, 1000)
        eps = Symbol("\\varepsilon")
        delta = Symbol("delta", positive=True)
        omg = Symbol("omega", positive=True)
        sys_msm = MultiTimeScaleSolution(
            system.odes[0],
            system.dvars[0],
            ivar=system.ivar,
            omega=S.One,
            order=2,
            eps=eps,
        )
        sys_sol = sys_msm.general_solution
        sys_an = (
            sys_sol.subs({eps: 0.1, delta: 0.1, omg: 0.1})
            .with_ics([0.1, 0])
            .numerized()
            .compute_solution(t_span, [0.1, 0])
        )
        sys_num = (
            sys_msm.subs({eps: 0.1, delta: 0.1, omg: 0.1})
            .numerized()
            .compute_solution(t_span, [0.1, 0])
        )
        sys_diff = sys_num - sys_an
        display(ReportText("Analytical solution"))
        display(sys_an.plot())
        display(ReportText("Numerical solution"))
        display(sys_num.plot())
        display(ReportText("Difference between numerical and analytical solutions"))
        display(sys_diff.plot())


class NaturalFrequenciesAnalysisComponent(ReportComponent):

    title = "Determination of natural frequencies"

    @property
    def header_text(self):

        return """In order to describe the basic properties of the system, its inertia and stiffness matrices will serve to determine an influence of the system parameters on its natural vibration frequencies (circular frequencies). The system inertia matrix provides a complete description of the system components' behaviour under the action of accelerations on individual degrees of freedom. It describes an influence of inertia of the system's elements while the stiffness matrix describes the global response of the system to external forces. The individual components of the matrix result from elastic and geometrical properties of the whole structure.
        """

    @property
    def body_text(self):

        system = self.reported_object

        return f"The characteristic polynomial ${latex(Delta)}$ of the considered system derived based on its fundamental matrix is presented in {AutoMarker(Eq(Delta,system.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))}:"

    @property
    def footer_text(self):

        return """This allows for a more transparent notation of the vibration frequency and therefore easier physical analysis and tuning its TMD parameters."""

    def append_elements(self):

        system = self.reported_object
        t = system.ivar

        base_matrix = system._inertia_matrix().inv() * system._stiffness_matrix()
        #         display(base_matrix)

        nom = fraction(
            ((base_matrix.trace() / 2) ** 2 - base_matrix.det()).simplify().radsimp()
        )[0]
        denom = fraction(
            ((base_matrix.trace() / 2) ** 2 - base_matrix.det()).simplify().radsimp()
        )[1]
        delta_omg = sqrt(nom) / sqrt(denom)

        omg_mean = (base_matrix.trace() / 2).simplify()

        l1 = Symbol("l_alpha", positive=True)
        l2 = Symbol("l_beta", positive=True)
        g = Symbol("g", positive=True)
        omg_m = Symbol("\=\omega")
        delta = Symbol("\Delta")
        omg = Symbol("omega")
        omg_1 = Symbol("omega_1")
        omg_2 = Symbol("omega_2")

        delta_omg_simp = omg_m**2 + (g / l1) * simplify(base_matrix)[2]
        #         display(simplify(base_matrix)[1])

        display(ReportText(self.header_text))

        display(
            ReportText(
                f"""\n The inertia and stifness matrices are given in {AutoMarker(Eq(Symbol('M'),system._inertia_matrix(),evaluate=False))} and {AutoMarker(Eq(Symbol('M'),system._stiffness_matrix(),evaluate=False))}:
        """
            )
        )

        display(SympyFormula(Eq(Symbol("M"), system._inertia_matrix(), evaluate=False)))
        display(
            SympyFormula(Eq(Symbol("K"), system._stiffness_matrix(), evaluate=False))
        )

        r = Symbol("r")
        omg = Symbol("omega")

        display(
            ReportText(
                """The stiffness matrix makes possible to determine equilibrium positions of the system and together with the inertia matrix, to determine its natural frequencies creating its the fundamental matrix:
        """
            )
        )

        fund_mat = system.fundamental_matrix().subs({r: -omg})

        display(SympyFormula(Eq(Symbol("A"), fund_mat, evaluate=False)))

        display(
            ReportText(
                f"""By calculating a determinant of the fundamental matrix, the system natural frequencies can be determined.
        """
            )
        )

        display(SympyFormula(Eq(fund_mat.det().collect(omg**2), 0, evaluate=False)))

        display(
            ReportText(
                f"""Assuming certain dependencies, the formula {AutoMarker(Eq(Symbol('A'),fund_mat,evaluate=False))} can be rearranged and simplified to the following:
        """
            )
        )

        display(
            SympyFormula(Eq(Symbol("a"), fund_mat.det().collect(omg**2).coeff(omg**4)))
        )
        display(
            SympyFormula(Eq(Symbol("b"), fund_mat.det().collect(omg**2).coeff(omg**2)))
        )
        display(
            SympyFormula(
                Eq(
                    Symbol("c"),
                    fund_mat.det().collect(omg**2).subs({omg**4: 0, omg**2: 0}),
                )
            )
        )

        display(
            SympyFormula(
                Eq(Symbol("a") * omg**4 + Symbol("b") * omg**2 + Symbol("c"), 0)
            )
        )

        display(
            ReportText(
                f"""The purpose of the presented approach is to provide an opportunity to influence the eigenvalues and thus control the system resonant behaviour. The linear form of the equations of motion allows to obtain the system's natural frequencies ({AutoMarker(Eq(Symbol('omega_1'),sqrt(omg_mean-delta_omg)))},{AutoMarker(Eq(Symbol('omega_1'),sqrt(omg_mean+delta_omg)))}) utilising the fundamental matrix {AutoMarker(Eq(Symbol('A'),fund_mat,evaluate=False))}:
        """
            )
        )

        display(SympyFormula(Eq(omg_1, sqrt(omg_mean - delta_omg))))
        display(SympyFormula(Eq(omg_2, sqrt(omg_mean + delta_omg))))

        display(
            ReportText(
                """By creating an inverse matrix of matrix $M$ and multiplying it by matrix $K$, it is possible to identify and assign certain frequency relationships in the cosidered system:
        """
            )
        )

        display(
            SympyFormula(
                Eq(Symbol("A_{\omega^2}"), simplify(base_matrix), evaluate=False)
            )
        )

        display(
            SympyFormula(
                Eq(
                    Symbol("A_{\omega^2}"),
                    Matrix(
                        [
                            [Symbol("\omega_{MS}") ** 2, -Symbol("\omega_{C}") ** 2],
                            [-Symbol("\omega_{TMD}") ** 2, Symbol("\omega_{TMD}") ** 2],
                        ]
                    ),
                    evaluate=False,
                )
            )
        )

        display(
            ReportText(
                """where: $\omega_{MS}^2$ corresponds a frequency of a main system, $\omega_{C}^2$ corresponds a composed frequency and $\omega_{TMD}^2$ corresponds a frequency of a tuned mass damper.
        """
            )
        )

        display(
            ReportText(
                """Analyzing the natural frequency components, some dependencies can be noticed. An expression of "a base frequency" can be separated, which is changed by a $\Delta \omega$ value. It can be presented as follows:
        """
            )
        )

        display(SympyFormula(Eq(omg_1**2, omg_m**2 - delta * omg**2)))
        display(SympyFormula(Eq(omg_2**2, omg_m**2 + delta * omg**2)))

        display(
            ReportText(
                """And the individual components, where a base frequency is defined as:
        """
            )
        )

        display(SympyFormula(Eq(omg_m**2, omg_mean)))

        display(
            ReportText(
                """And the separated $\Delta \omega$ component containing also a base frequency term, after some rearangement can be simplified to the following expressions:
        """
            )
        )

        display(
            SympyFormula(
                Eq(Eq(delta * omg**2, delta_omg), sqrt(delta_omg_simp), evaluate=False)
            )
        )

        #         display(SympyFormula(Eq(delta*omg**2, sqrt(delta_omg_simp))))

        display(
            ReportText(
                f"""Finally, the natural frequencies can be simplified to the expressions presented in {AutoMarker(Eq(omg_1, omg_m**2 - sqrt(delta_omg_simp)))}, {AutoMarker(Eq(omg_2, omg_m**2 + sqrt(delta_omg_simp)))}
        """
            )
        )

        display(SympyFormula(Eq(omg_1, sqrt(omg_m**2 - sqrt(delta_omg_simp)))))
        display(SympyFormula(Eq(omg_2, sqrt(omg_m**2 + sqrt(delta_omg_simp)))))

        display(ReportText(self.footer_text))


class NormalisedNaturalFrequenciesAnalysisComponent(ReportComponent):

    title = "Determination of natural frequencies"

    @property
    def header_text(self):

        return """In order to describe the basic properties of the system, its inertia and stiffness matrices will serve to determine an influence of the system parameters on its natural vibration frequencies. The system inertia matrix provides a complete description of the system components' behaviour under the action of accelerations on individual degrees of freedom. It describes an influence of inertia of the system's elements while the stiffness matrix describes the global response of the system to external forces. The individual components of the matrix result from elastic and geometrical properties of the whole structure.
        """

    @property
    def footer_text(self):

        return """This allows for a more transparent notation of the system natural frequencies. Therefore an easier physical analysis and tuning the system TMD parameters can be done in a more comprehensive and controllable manner."""

    def append_elements(self):

        system = self.reported_object
        t = system.ivar

        base_matrix = system._inertia_matrix().inv() * system._stiffness_matrix()
        base = simplify(base_matrix)
        base_matrix = base / base[0]
        #         display(base_matrix)

        nom = fraction(
            ((base_matrix.trace() / 2) ** 2 - base_matrix.det()).simplify().radsimp()
        )[0]
        denom = fraction(
            ((base_matrix.trace() / 2) ** 2 - base_matrix.det()).simplify().radsimp()
        )[1]
        delta_omg = sqrt(nom) / sqrt(denom)

        omg_mean = (base_matrix.trace() / 2).simplify()

        l1 = Symbol("l_alpha", positive=True)
        l2 = Symbol("l_beta", positive=True)
        g = Symbol("g", positive=True)
        delta = Symbol("\Delta")
        omg = Symbol("omega_n")
        omg_1 = Symbol("\omega_{n_{1}}")
        omg_2 = Symbol("\omega_{n_{2}}")
        omg_n = Symbol("\\bar \omega_n")
        omg12_omg21 = Symbol("\omega_{12}") * Symbol("\omega_{21}")

        delta_omg_simp = omg_n**4 - simplify(base_matrix)[1] * simplify(base_matrix)[2]
        #         display(simplify(base_matrix)[1])

        display(ReportText(self.header_text))

        #         display(ReportText(f'''\n The inertia and stifness matrices are given in ({AutoMarker(Eq(Symbol('M'),system._inertia_matrix(),evaluate=False))}) and ({AutoMarker(Eq(Symbol('K'),system._stiffness_matrix(),evaluate=False))}):
        #         '''))

        #         display(SympyFormula(Eq(Symbol('M'),system._inertia_matrix(),evaluate=False)))
        #         display(SympyFormula(Eq(Symbol('K'),system._stiffness_matrix(),evaluate=False)))

        r = Symbol("r")
        omg = Symbol("omega")

        #         display(ReportText('''The stiffness matrix makes possible to determine equilibrium positions of the system and together with the inertia matrix, to determine its natural frequencies creating its the fundamental matrix:
        #         '''))

        fund_mat = system.fundamental_matrix().subs({r: -omg})

        display(SympyFormula(Eq(Symbol("A"), fund_mat, evaluate=False)))

        display(
            ReportText(
                f"""By calculating a determinant of the fundamental matrix, the system natural frequencies can be determined.
        """
            )
        )

        #         display(SympyFormula(Eq(fund_mat.det().collect(omg**2),0,evaluate=False)))

        display(
            ReportText(
                f"""Assuming certain dependencies, the formula ({AutoMarker(Eq(Symbol('A'),fund_mat,evaluate=False))}) can be rearranged into a biquadratic polynomial and simplified to the following, where its roots are the system eigenvalues:
        """
            )
        )

        display(
            SympyFormula(Eq(Symbol("a"), fund_mat.det().collect(omg**2).coeff(omg**4)))
        )
        display(
            SympyFormula(Eq(Symbol("b"), fund_mat.det().collect(omg**2).coeff(omg**2)))
        )
        display(
            SympyFormula(
                Eq(
                    Symbol("c"),
                    fund_mat.det().collect(omg**2).subs({omg**4: 0, omg**2: 0}),
                )
            )
        )

        display(
            SympyFormula(
                Eq(Symbol("a") * omg**4 + Symbol("b") * omg**2 + Symbol("c"), 0)
            )
        )

        display(
            ReportText(
                f"""The purpose of the presented approach is to provide an opportunity to influence the eigenvalues and thus control the system resonant behaviour. The linear form of the equations of motion allows to obtain the system's natural frequencies ({AutoMarker(Eq(omg_1,sqrt(omg_mean-delta_omg)))}),({AutoMarker(Eq(omg_2,sqrt(omg_mean+delta_omg)))}) utilising the fundamental matrix ({AutoMarker(Eq(Symbol('A'),fund_mat,evaluate=False))}):
        """
            )
        )

        #         display(SympyFormula(Eq(Symbol('\omega^2_{n_{1}}'),Symbol('\\frac{\omega^2_1}{\omega^2_m}'))))
        #         display(SympyFormula(Eq(Symbol('\omega^2_{n_{2}}'),Symbol('\\frac{\omega^2_2}{\omega^2_m}'))))

        display(SympyFormula(Eq(omg_1, sqrt(omg_mean - delta_omg))))
        display(SympyFormula(Eq(omg_2, sqrt(omg_mean + delta_omg))))

        omg11_str = r"$\omega_{11}$"
        Aomg_str = r"$\omega_{n^2}$"

        display(
            ReportText(
                f"""By creating an inverse matrix of matrix $M$ and multiplying it by matrix $K$, it is possible to identify and assign certain frequency relationships in the cosidered system. The matrix elements were normalised to the value of element {omg11_str} ({AutoMarker(Eq(Symbol('{omg11_str}'),simplify(system._inertia_matrix().inv()*system._stiffness_matrix()[0],evaluate=False)))}) corresponding to a frequency of the main system. This results in a normalised matrix presented in ({AutoMarker(Eq(Symbol('A_{Aomg_str}'),simplify(base_matrix),evaluate=False))}):
        """
            )
        )

        display(
            SympyFormula(
                Eq(
                    Symbol("\omega_11"),
                    simplify(
                        system._inertia_matrix().inv() * system._stiffness_matrix()
                    )[0],
                    evaluate=False,
                )
            )
        )

        display(
            SympyFormula(
                Eq(Symbol("A_{\omega_n^2}"), simplify(base_matrix), evaluate=False)
            )
        )

        display(
            ReportText(
                """Analyzing the natural frequency components, some dependencies can be noticed. An expression of "a base frequency" ($\omega_n^2$) can be separated. Its value is being changed by some $\Delta \omega^2$ value. It can be presented as follows:
        """
            )
        )

        display(SympyFormula(Eq(omg_1**2, omg_n**2 - delta * omg**2)))
        display(SympyFormula(Eq(omg_2**2, omg_n**2 + delta * omg**2)))

        display(
            ReportText(
                """And the individual components, where a base frequency $\omega_n^2$ is defined as:
        """
            )
        )

        display(SympyFormula(Eq(omg_n**2, omg_mean)))

        display(
            ReportText(
                """And the separated $\Delta \omega^2$ component containing also a base frequency term, after some rearrangement can be simplified to the following expressions:
        """
            )
        )

        display(
            SympyFormula(
                Eq(Eq(delta * omg**2, delta_omg), sqrt(delta_omg_simp), evaluate=False)
            )
        )

        #         display(SympyFormula(Eq(delta*omg**2, sqrt(omg_mean**2-simplify(base_matrix)[1]*simplify(base_matrix)[2]))))

        display(
            ReportText(
                f"""Finally, the natural frequencies can be simplified to the expressions presented in ({AutoMarker(Eq(omg_1, sqrt(omg_n**2 - sqrt(delta_omg_simp))))}), ({AutoMarker(Eq(omg_2, sqrt(omg_n**2 + sqrt(delta_omg_simp))))}):
        """
            )
        )

        display(SympyFormula(Eq(omg_1, sqrt(omg_n**2 - sqrt(delta_omg_simp)))))
        display(SympyFormula(Eq(omg_2, sqrt(omg_n**2 + sqrt(delta_omg_simp)))))

        display(
            ReportText(
                f"""By noting certain dependencies of the normalised matrix ({AutoMarker(Eq(Symbol('A_{omega_n^2}'),simplify(base_matrix),evaluate=False))}), one can simplify the notation of the natural frequencies into a formula containing the basic dependencies:
        """
            )
        )

        display(
            SympyFormula(
                Eq(omg_1, omg_n * sqrt(1 - sqrt(1 - (omg12_omg21) / omg_n**2)))
            )
        )
        display(
            SympyFormula(
                Eq(omg_2, omg_n * sqrt(1 + sqrt(1 - (omg12_omg21) / omg_n**2)))
            )
        )

        display(
            ReportText(
                """where $\omega_{12}$ and $\omega_{21}$ are the elements of the normalised matrix.
        """
            )
        )

        display(
            ReportText(
                """Having the eigenvalue formulas simplified, now the system properties can be put back into the equations resulting in final expressions:
        """
            )
        )

        from dynpy.solvers.linear import ODESystem

        eigs = ODESystem(
            odes=Matrix([system.lhs[0], system.lhs[1]]), dvars=Matrix(system.dvars)
        )._eigenvals_simplified()

        display(SympyFormula(Eq(omg_1, eigs[0])))
        display(SympyFormula(Eq(omg_2, eigs[3])))

        display(ReportText(self.footer_text))


class ODEInitCodeComponent(ReportComponent):

    title = "Creating a linear differential equation"

    def append_elements(self):

        from dynpy.solvers.tools import CodePrinter

        system = self.reported_object
        ode_sys = system

        _dvars = system.dvars
        _odes = system.odes

        _ode_order = system.ode_order

        if _ode_order == 1:
            cls_name = "FirstOrderODESystem"
        else:
            cls_name = "ODESystem"

        string_gc = CodePrinter(_odes[0])._generate_code()

        display(
            ReportText(
                f"Creating an instance of a class {cls_name} requires defining the following variables and equations"
            )
        )

        display(
            ObjectCode(
                f"""


from dynpy.solvers.linear import *
from dynpy.models.odes.linear import *
from dynpy.models.odes.numerical import *

from sympy import *

{string_gc}

dvars = {_dvars}
odes = {_odes}
ode_order = {_ode_order}

odesys = {cls_name}(odes = odes , dvars = dvars, ode_order = ode_order)

display(odesys)

        """
            )
        )

        display(
            ReportText(
                f"After executing the instruction, we have a ready instance of the class object {cls_name}"
            )
        )


class ODERootsComponent(ReportComponent):

    title = "Creating roots of the characteristic equation"

    def append_elements(self):
        system = self.reported_object
        #         _roots = system.fundamental_matrix
        _roots = system._as_fode()._fundamental_matrix

        display(ReportText(f"Creating a matrix of roots from a called out class :"))
        display(SympyFormula(_roots))


steady_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system.steady_solution
"""


class ODERootsCodeComponent(ReportComponent):

    title = "Creating roots of the characteristic equation"

    def append_elements(self):
        # Class was impelemnted in ODERootsComponent. There you can find call code and call
        pass


class ODESteadySolutionComponent(ReportComponent):

    title = "Wyznaczanie rozwiązania ogólnego"

    @property
    def import_text(self):

        return "Aby zaimportować przykładowy system  i wywołać rozwiązanie szczególne należy wywołać następujący kod:"

    @property
    def body_text(self):

        return "Rozwiązanie szczególne przyjmuje postać:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(steady_code))

        display(ReportText(self.body_text))

        display((SympyFormula(system.steady_solution)))


general_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system.general_solution
"""


class ODEGeneralSolutionComponent(ReportComponent):

    title = "Determining the general solution"

    @property
    def import_text(self):

        return "To import the sample system and call the generic solution, call the following code:"

    @property
    def body_text(self):

        return "The general solution takes the form:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(general_code))

        display(ReportText(self.body_text))

        display((SympyFormula(system.general_solution)))


import_code = """
from dynpy.models.odes.numerical import *

system=LinearFirstOrder.from_reference_data()
system.as_eq()  #returns equation
system.as_matrix()  #returns matrix
system.as_iterable()  #returns list
system.as_dict()  #returns dictionary
system.as_eq_list()  #returns an equation in the form of a list
system._as_fode()  #returns a first-order equation
system.lhs #returns left-hand side of equation
system.rhs #returns right-hand side of equation
"""


class ODESystemRepresentationComponent(ReportComponent):

    title = "Diverse representations of the system"

    @property
    def import_text(self):

        return "Methods that enable a variety of system representations should be called as follows:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(import_code))

        display(ReportText(r"Method .as_eq()"))
        display((system.as_eq()))

        display(ReportText(r"Method .as_matrix()"))
        display((system.as_matrix()))

        display(ReportText(r"Method .as_iterable()"))
        display((system.as_iterable()))

        display(ReportText(r"Method .as_dict()"))
        display((system.as_dict()))

        display(ReportText(r"Method .as_eq_list()"))
        display((system.as_eq_list()))

        display(ReportText(r"Method ._as_fode()"))
        display((system._as_fode()))

        display(ReportText(r"Method .lhs"))
        display((system.lhs))

        display(ReportText(r"Method .rhs"))
        display((system.rhs))


class ODESystemAdditionComponent(ODESystemOperations):

    title = "Adding two ODESystem Matrix"

    def append_elements(self):
        from sympy import Function, Matrix, Symbol

        systems = self.reported_object
        l = len(systems)
        matrix_add = self.get_zero_ode_sys()

        matrix_mat_lhs = Matrix([0])
        matrix_mat_rhs = Matrix([0])

        count = 1
        display(ReportText("#### Checking the sum of matrices:"))
        for ode_sys in systems:
            display(ReportText(f"Matrix {count}"))
            display(ode_sys)
            count += 1
            matrix_add += ode_sys

            matrix_mat_lhs += ode_sys.lhs
            matrix_mat_rhs += ode_sys.rhs

        display(ReportText("ODEsystem matrix sum"))
        display(matrix_add)
        display(ReportText("Sympy matrix sum"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemSubtractionComponent(ODESystemOperations):

    title = "Subtractiing two ODESystem Matrix"

    def append_elements(self):
        systems = self.reported_object
        l = len(systems)
        matrix_sub = systems[0]

        matrix_mat_lhs = systems[0].lhs
        matrix_mat_rhs = systems[0].rhs

        display(ReportText("#### Checking the difference of a matrix:"))
        display(ReportText(f"Matrix 1"))
        display(systems[0])
        count = 2
        for ode_sys in systems[1:]:
            display(ReportText(f"Matrix {count}"))
            display(ode_sys)
            count += 1
            matrix_sub -= ode_sys

            matrix_mat_lhs -= ode_sys.lhs
            matrix_mat_rhs -= ode_sys.rhs

        display(ReportText("ODEsystem matrix difference"))
        display(matrix_sub)
        display(ReportText("Sympy matrix difference"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemMultiplicationComponent(ODESystemOperations):

    title = "Multiplying two ODESystem Matrix"

    def append_elements(self):
        systems = self.reported_object
        l = len(systems)
        matrix_mul = 1
        count = 1

        matrix_mat_lhs = Matrix([1])
        matrix_mat_rhs = Matrix([1])

        display(ReportText("#### Checking the product of matrices:"))
        for ode_sys in systems:
            display(ReportText(f"Matrix {count}"))
            display(ode_sys)
            count += 1
            matrix_mul *= ode_sys

            matrix_mat_lhs *= ode_sys.lhs
            matrix_mat_rhs *= ode_sys.rhs

        display(ReportText("ODEsystem matrix product"))
        display(matrix_mul)
        display(ReportText("Sympy matrix product"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemExponentiationComponent(ReportComponent):
    title = "Exponentiation of ODESystem Matrix"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system

        display(ReportText("#### Checking the power of a matrix:"))
        display(ReportText("Matrix"))
        display(ode_sys)

        display(ReportText("ODEsystem Matrix power"))

        matrix_exp = ode_sys * ode_sys
        display(matrix_exp)

        matrix_mat_lhs = ode_sys.lhs * ode_sys.lhs
        matrix_mat_rhs = ode_sys.rhs * ode_sys.rhs

        display(ReportText("Sympy matrix power"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


excitation_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system.as_type(FirstOrderLinearODESystemWithHarmonics)._get_excitation_comps
"""


class PredictionOfSteadySolutionComponent(ReportComponent):

    title = "Determining the predicted solution of a special equation"

    @property
    def import_text(self):

        return "To import an example system and determine its predicted special solution, call the following code:"

    @property
    def body_text(self):

        return "The solution takes the form:"

    def append_elements(self):
        from dynpy.solvers.linear import FirstOrderLinearODESystemWithHarmonics

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars
        display(ReportText(self.import_text))
        display(ObjectCode(excitation_code))
        display(ReportText(self.body_text))
        system.as_type(FirstOrderLinearODESystemWithHarmonics)._get_excitation_comps
        ext_comps_list = system.as_type(
            FirstOrderLinearODESystemWithHarmonics
        )._get_excitation_comps
        display(ext_comps_list)
        for comp, amp in ext_comps_list:

            display(SympyFormula(comp))
            display(SympyFormula(amp))


class HomoPredictionIntroComponent(ReportComponent):

    title = "Second-order differential equation - prediction method"

    def append_elements(self):
        system = self.reported_object
        dvar = system.dvars[0]
        ivar = system.ivar
        display(
            ReportText(
                "One of the basic types of second-order differential equations is the homogeneous equation with constant coefficients. The process of solving this type of equation involves predicting a solution that takes the general form:"
            )
        )
        display(SympyFormula(Eq(dvar, Symbol("C") * exp(Symbol("r") * ivar))))
        display(
            ReportText(
                "Based on the obtained roots (solutions) of the predicted equation, the final general equation is formulated. Below is an example problem illustrating the process of solving this type of equation:"
            )
        )


class MainPredictionComponent(ReportComponent):

    title = "Predicting the homogeneous solution"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        C = Symbol("C", positive=True)
        r = Symbol("r")
        y = system.dvars[0]
        x = system.ivar
        sym_dict = system._get_dvars_symbols()
        sys_sub = Eq(system.lhs[0], system.rhs[0])

        prdct = Eq(y, C * exp(r * x))
        prdct_diff = Eq(prdct.lhs.diff(x), prdct.rhs.diff(x))
        prdct_ddiff = Eq(prdct.lhs.diff(x, 2), prdct.rhs.diff(x, 2))

        prdct_dict = {
            y: prdct.rhs,
            y.diff(x): prdct_diff.rhs,
            y.diff(x, 2): prdct_ddiff.rhs,
        }

        r_eq = Eq(system.fundamental_matrix()[0], 0)
        r0_eq = Eq(y.subs(sym_dict), Symbol("r^0"))
        r1_eq = Eq(y.diff(x).subs(sym_dict), Symbol("r^1"))
        r2_eq = Eq(y.diff(x, 2).subs(sym_dict), Symbol("r^2"))

        display(
            ReportText(
                "The first step in solving the given equation is to determine the derivatives of the predicted solution. Recalling the original equation, which we predict as our solution, has the following form:"
            )
        )
        display(SympyFormula(prdct.subs(sym_dict)))
        display(
            ReportText(
                "And then finding the necessary derivatives. In the case of second-order equations, as the name suggests, the highest derivative is of the second degree. Therefore, we need the first two derivatives of the predicted solution. The first derivative is:"
            )
        )
        display(SympyFormula(prdct_diff.subs(sym_dict)))
        display(ReportText("Second derivative:"))
        display(SympyFormula(prdct_ddiff.subs(sym_dict)))
        display(
            ReportText(
                "The next step is to substitute the respective predictions into the corresponding derivatives in the original equation:"
            )
        )
        display(SympyFormula(sys_sub.subs(prdct_dict)))
        display(
            ReportText(
                "One can notice recurring elements and factor them out before the parentheses:"
            )
        )
        display(SympyFormula(sys_sub.subs(prdct_dict).simplify()))
        display(
            ReportText(
                "The fundamental principle of mathematics states that if the product of two factors equals zero, then at least one of those factors must be zero. However, in the method of predictions, we do not consider two cases because when the constant of integration, and thus the exponential function, equals zero, the task becomes trivial and does not require further analysis. Therefore, we focus on the second factor, commonly referred to as the characteristic equation, which will provide us with the necessary roots to predict the final equation. The equation undergoing further analysis has the following form:"
            )
        )
        display(SympyFormula(r_eq))
        display(
            ReportText(
                "After solving a certain number of problems, one can observe a pattern arising from the above prediction, which allows skipping this step and substituting the appropriate values according to the following substitutions:"
            )
        )
        display(SympyFormula(r0_eq))
        display(SympyFormula(r1_eq))
        display(SympyFormula(r2_eq))


class RootsAnalysisComponent(ReportComponent):

    title = "Solution of the characteristic equation and general analysis of its solutions in the context of solving second-order differential equations"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        sym_dict = system._get_dvars_symbols()
        y = system.dvars[0].subs(sym_dict)
        x = system.ivar
        char = system.fundamental_matrix()[0]
        char_eq = Eq(char, 0)
        roots = system._eigenvalues()
        C1 = Symbol("C_1", positive=True)
        C2 = Symbol("C_2", positive=True)
        r0 = Symbol("r_0")
        r1 = Symbol("r_1")
        a = Symbol("alpha")
        b = Symbol("beta")
        i = Symbol("i")

        r1_img = Eq(r0, a + i * b)
        r2_img = Eq(r1, a - i * b)

        pred_real = Eq(y, C1 * exp(r0 * x) + C2 * exp(r1 * x))
        pred_one = Eq(y, C1 * exp(r0 * x) + C2 * x * exp(r0 * x))
        pred_imagine = Eq(y, exp(a * x) * (C1 * cos(b * x) + C2 * sin(b * x)))

        sol_eq = Eq(y, system.general_solution.rhs[0])

        display(
            ReportText(
                f"In second-order differential equations, the characteristic equation takes the form of a quadratic equation. Such an equation can yield three types of roots in its solution, and depending on the type obtained, we select the general solution of this equation. When the calculated discriminant is greater than zero, we obtain two real roots in the form ${latex(r0)}$ and ${latex(r1)}$. In this case, the predicted solution takes the form:"
            )
        )
        display(SympyFormula(pred_real))
        display(
            ReportText(
                f"When the calculated discriminant is equal to 0, we obtain one root in the form ${latex(r0)}$. In this case, the predicted solution takes the form:"
            )
        )
        display(SympyFormula(pred_one))
        display(
            ReportText(
                "When the calculated discriminant is less than 0, we obtain two complex roots in the form:"
            )
        )
        display(SympyFormula(r1_img))
        display(SympyFormula(r2_img))
        display(
            ReportText(
                f"Where ${latex(a)}$ is the real part of the root, and ${latex(b)}$ is the imaginary part. The predicted solution in this case takes the following form:"
            )
        )
        display(SympyFormula(pred_imagine))
        display(ReportText("The characteristic equation in the given example is:"))
        display(SympyFormula(char_eq))
        display(ReportText("The roots of this equation take the form:"))
        i = 0
        lista = system._roots_list()
        for elem in lista:
            if elem != 0:
                display(SympyFormula(Eq(Symbol(f"r_{i}"), elem), marker=None))
                i += 1
        display(
            ReportText(
                "After identifying the type of roots, we can express the solution as:"
            )
        )
        display(SympyFormula(sol_eq))


class SolutionComparisonComponent(ReportComponent):

    title = "Multiple time scale solution and ODESystem solution"

    def append_elements(self):

        from dynpy.solvers.linear import (
            FirstOrderLinearODESystem,
            FirstOrderLinearODESystemWithHarmonics,
            ODESystem,
        )
        from dynpy.solvers.nonlinear import MultiTimeScaleSolution

        system = self.reported_object

        display(Markdown("`as_type` method used:"))
        display(
            SympyFormula(
                Eq(
                    system.as_type(FirstOrderLinearODESystem).solution.lhs[0],
                    system.as_type(FirstOrderLinearODESystem).solution.rhs[0],
                )
            )
        )
        display(
            SympyFormula(
                Eq(
                    system.as_type(FirstOrderLinearODESystem).solution.lhs[1],
                    system.as_type(FirstOrderLinearODESystem).solution.rhs[1],
                )
            )
        )

        display(ReportText("Multiple time scales solution is following:"))
        display(
            SympyFormula(
                Eq(
                    system.as_type(FirstOrderLinearODESystemWithHarmonics).solution.lhs[
                        0
                    ],
                    system.as_type(FirstOrderLinearODESystemWithHarmonics).solution.rhs[
                        0
                    ],
                )
            )
        )
        display(
            SympyFormula(
                Eq(
                    system.as_type(FirstOrderLinearODESystemWithHarmonics).solution.lhs[
                        1
                    ],
                    system.as_type(FirstOrderLinearODESystemWithHarmonics).solution.rhs[
                        1
                    ],
                )
            )
        )


class RawSolutionComparisonComponent(ReportComponent):

    title = "Multiple time scale solution and ODESystem solution with ODESystem output"

    def append_elements(self):

        from dynpy.solvers.linear import (
            FirstOrderLinearODESystem,
            FirstOrderLinearODESystemWithHarmonics,
            ODESystem,
        )
        from dynpy.solvers.nonlinear import MultiTimeScaleSolution

        system = self.reported_object
        # if type(system) == ODESystem:
        display(SympyFormula(system.as_type(FirstOrderLinearODESystem).solution))
        display(ReportText("as_type method used:"))
        display(
            SympyFormula(
                system.as_type(FirstOrderLinearODESystemWithHarmonics).solution
            )
        )

        # elif type(system) == MultiTimeScaleSolution:
        display(ReportText("Multiple time scales solution is following:"))
        display(
            SympyFormula(
                system.as_type(FirstOrderLinearODESystemWithHarmonics).solution
            )
        )

    # else:


#       display(ReportText("Given system is neither ODESystem nor MultiTimeScaleSolution"))


class SolutionDiagnosticComponent(ReportComponent):

    title = "Multiple time scale solution and ODESystem solution"

    def append_elements(self):

        system = self.reported_object._as_fode()

        display(ReportText("_auxiliary_fundamental_matrix:"))
        display(SympyFormula(system._auxiliary_fundamental_matrix))

        display(ReportText("_reduced_fundamental_matrix:"))
        display(SympyFormula(system._reduced_fundamental_matrix))

        display(ReportText("_trig_fundamental_matrix:"))
        display(SympyFormula(system._trig_fundamental_matrix))

        display(ReportText("_reduced_eigenvalues:"))
        display(SympyFormula(system._reduced_eigenvalues))

        display(ReportText("_trig_eigenvalues:"))
        display(SympyFormula(system._trig_eigenvalues))

        display(ReportText("_reduced_modes:"))
        display(SympyFormula(system._reduced_modes))

        display(ReportText("_combined_modes:"))
        display(SympyFormula(system._combined_modes))

        display(ReportText("_trig_modes:"))
        display(SympyFormula(system._trig_modes))

        display(ReportText("eigenvalues:"))
        display(SympyFormula(system.eigenvalues))

        display(ReportText("_is_proportional_damping:"))
        display(SympyFormula(system._is_proportional_damping))

        display(ReportText("_general_solution:"))
        display(SympyFormula(system._general_solution))


class ParticularDerivativesComponent(ReportComponent):

    title = "Derivatives of approximation functions"

    def append_elements(self):

        from dynpy.solvers.linear import (
            FirstOrderLinearODESystem,
            FirstOrderLinearODESystemWithHarmonics,
            ODESystem,
        )
        from dynpy.solvers.perturbational import MultiTimeScaleSolution

        system = self.reported_object

        display(
            ReportText(
                "To clarify the solving process, the derivatives of the approximation functions are presented below, both with and without the substitution of time scales and their corresponding derivatives."
            )
        )
        for index, value in enumerate(system._scales_formula):
            display(
                ReportText(
                    f"For the function ${latex(system.approximation_function(index)[0])}$:"
                )
            )
            #             display(SympyFormula(system.approximation_function(index)[0]))
            display(SympyFormula(system.part_derivative(index).expand()))
            display(SympyFormula(system.sec_part_derivative(index).expand()))

        for index, value in enumerate(system._scales_formula):
            display(
                ReportText(
                    f"After the substitution of time scales and their derivatives for the function ${latex(system.approximation_function(index)[0])}$:"
                )
            )
            #             display(SympyFormula(system.approximation_function(index)[0]))
            display(SympyFormula(system.part_derivative_subs(index).expand()))
            display(SympyFormula(system.sec_part_derivative_subs(index).expand()))
