from  ..mechanics import *



class TitlePageComponent(Environment):
    
    latex_name='titlepage'
    
    def __init__(self, system=None, options=None, arguments=None, start_arguments=None,
                 **kwargs):
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

        
        
        super().__init__(options=options, arguments=arguments, start_arguments=start_arguments,**kwargs)
        
        if self.system is not None:

        
            system = self.system


            
            self.append(NoEscape('\centering'))

            self.append(NoEscape('\\Huge DRGANIA MECHANICZNE \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))

            
            
            if len(system.dvars)==1:
                dof_str = 'JEDNYM STOPNIU SWOBODY'
            else:
                dof_str = 'WIELU STOPNIACH SWOBODY'
                
            if system._dissipative_potential==0 or system._dissipative_potential is None:
                damping_str = 'NIETŁUMIONE'
            else:
                damping_str = 'TŁUMIONE'

           
            self.append(NoEscape(f'\\Large {damping_str} UKŁADY O {dof_str} \n \n'))    
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(NoEscape(f'{system._label} \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(Command('MyAuthor'))
            self.append(NoEscape(f'\\par'))
            self.append(Command('vspace',arguments='1cm'))
            
            #self.append(Command('vspace',arguments='1cm'))
            #self.append(NoEscape(f'\\protect\\par'))
            #self.append(NewLine())
            self.append(Command('MyDate'))


class ODESystemComponent(ReportComponent):
    
    title="Differential equations"
    @property
    def header_text(self):

        return "The investigated system is described by differential equations being as follows:"

        
    @property
    def footer_text(self):

        return "To solve the problem serve several methods depending on the equation type."

    def append_elements(self):

        system = self._system


        display(ReportText(  self.header_text   ))

        display(SympyFormula( system ))

        display(ReportText(  self.footer_text   ))

        
class ODESystemCodeComponent(ReportComponent):
    
    title="Differential equations"
    @property
    def header_text(self):

        return "The code that enables to obtain differential equations is as follows:"

        
    @property
    def footer_text(self):

        return "Presented method can change a type of the equation represantation."

    def append_elements(self):

        system = self._system


        display(ReportText(  self.header_text   ))

        display(Markdown('\t from dynpy.solvers.linear import ODESystem  \n \t ode=ODESystem([]) \n \t display(ode)'))

        display(ReportText(  self.footer_text   ))
        

class VariablesComponent(ReportComponent):
    
    title="System's variables"
    @property
    def header_text(self):

        return "As dynamic systems's behaviour is described by an ordinary differential equation, the variables of the equation are as follows:"

        
    @property
    def footer_text(self):

        return "The variables allow to study and analyse the system's dynamics varying with time."

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text   ))

        display(SympyFormula(  system.ivar))
        display(SympyFormula(  system.dvars))

        display(ReportText(  self.footer_text   ))



class GoverningEquationComponent(ReportComponent):
    
    title="Equation of motion"
    
    @property
    def header_text(self):
        system = self._system

        return f"The equation of motion was derived based on physical principles and governing laws regarding the considered system. The equation is described by the following expression ({AutoMarker(Eq(system.odes[0],0))})"

    @property
    def footer_text(self):
        return f"The determined equation constitutes a mathematical description of the dynamic properties of the system. Further analysis allows for an effective study on the modelled object's operation and determination of its mechanical parameters."
    
    def append_elements(self):
        
        system = self._system

        display(ReportText(self.header_text))

        display(SympyFormula(Eq(system.odes[0],0)))

        display(ReportText(self.footer_text))


class ApproximatedNonlinearGoverningEquationComponent(ReportComponent):
    
    title="Approximated equation of motion"
    @property
    def header_text(self):

        return "The fully nonlinear form of the governing equation can be approximated by Taylor series expansion and then the system's motion equation is represented by the following:"

    @property
    def footer_text(self):

        return "Solving the nonlinear type of problems requires an involvement of dedicated methods of nonlinear theory of vibration."

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text   ))

        display(SympyFormula( Eq(system.approximated(3).simplify().expand(),0) ))

        display(ReportText(  self.footer_text   ))
        
class FundamentalMatrixComponent(ReportComponent):
    
    title="Determination of fundamental matrix"
    
    @property
    def header_text(self):

        return "The system physical properties of mass and stiffness are the basis to formulate the fundamental matrix of the system:"

    @property
    def body_text(self):
        
        system = self._system
        system_lin=system.linearized()
        Delta = Symbol('Delta')

        return f"The characteristic polynomial ${latex(Delta)}$ of the considered system derived based on its fundamental matrix is presented in {AutoMarker(Eq(Delta,system_lin.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))}:"
    
    @property
    def footer_text(self):

        return "The system fundamental matrix and its characteristic equation allows for a determination of the system unsteady state and its eigenfrequencies related to free vibration."
    
    def append_elements(self):

        system = self._system
        system_lin=system.linearized()


        display(ReportText(self.header_text))

        display((SympyFormula(Eq(Symbol('A'),system_lin.fundamental_matrix(),evaluate=False))))

        display(ReportText(self.body_text))

        display((SympyFormula(Eq(Delta,system_lin.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))))

        display(ReportText(self.footer_text))


class ODECharecteristicPolynomialComponent(ReportComponent):

    title="Characteristic polynomial"
    
    @property
    def header_text(self):
        return "Characterisitc pylynomial for the considered system of equations of motion is the following:"
    @property
    def footer_text(self):
        return "Calculation of roots of this polynomial allows to find natural frequencies of the system."

    def append_elements(self):
        system = self._system
        delta=Symbol('Delta')

        display(ReportText(  self.header_text   ))
        polynomial=system._ode_system.char_polynomial()
        display(SympyFormula(Eq(delta,polynomial)))

        display(ReportText(  self.footer_text   ))

class ODECharecteristicPolynomialCodeComponent(ReportComponent):

    title="Characteristic polynomial"
    
    @property
    def header_text(self):
        return "Presented piece of code generates the characteristic pylynomial of the desired system:"
    @property
    def footer_text(self):
        return "The following code is based on dynpy mechanical module."

    def append_elements(self):
        system = self._system
        
        display(ReportText(  self.header_text   ))

        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *
    from dynpy.solvers.linear import*

    dyn_sys = {system.__class__.__name__}()

    display( Eq(Symbol('Delta'),
                     dyn_sys._ode_system.char_polynomial()) )


'''))


        display(ReportText(  self.footer_text   ))
        
class ZerothOrderApproximatedEqComponent(ReportComponent):
    
    title="Zeroth-order approximation"
    
    
    @property
    def header_text(self):
        
        system = self._system
        system.general_solution(1)
#         approx_sys = system.approximated(3)

        zeroth_ord_eq = Symbol('Omega')#system.nth_eoms_approximation(0)
        eps = Symbol('\\varepsilon')

        return "The ordering and separation of equation {AutoMarker(system.odes[0])} in terms of the power of a small parameter {eps} leads to obtaining a recursive sequence of linear equations of motion leading to the solution of the nonlinear equation. The zeroth-order approximate linear equation is given in {AutoMarker(zeroth_ord_eq)} where the dependency {approx_fun} was assumed for the zeroth-order solution of the time variable:"
    
    @property
    def middle_text(self):
        
        system = self._system
        system.general_solution(1)
        zeroth_ord_approx_eq = Symbol('Omega')#Eq(system.nth_eoms_approximation(0).lhs[1]-system.nth_eoms_approximation(0).rhs[1],0)

        return f"Or transformed to the second order ODE {AutoMarker(zeroth_ord_approx_eq)}:"

        
    @property
    def footer_text(self):
        
        system = self._system
        t_list = system.t_list

        return f"Since {t_list[0]} and {t_list[1]} are treated as independent, the differential equation becomes a partial differential equation for a function of two variables {t_list[0]} and {t_list[1]}. Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}."

    def append_elements(self):
        
        system = self._system
        system.general_solution(1)

        t_list = system.t_list
        zeroth_ord_approx = system.eoms_approximation_list()[0]
        zeroth_ord_approx_eq = Eq(zeroth_ord_approx.lhs-zeroth_ord_approx.rhs,0,evaluate=False)

        display(ReportText(  self.header_text   ))

        display(SympyFormula(zeroth_ord_approx))

        display(ReportText(  self.middle_text   ))

        display(SympyFormula(zeroth_ord_approx_eq))

        display(ReportText(  self.footer_text   ))
        
class FirstOrderApproximatedEqComponent(ReportComponent):
    
    title="First-order approximation"
    
    
    @property
    def header_text(self):
        
        system = self._system
        first_ord_eq = system.eoms_approximation_list()[1]

        return f"The next component of a recursive sequence of linear equations of motion leading to the solution of the considered nonlinear one is given in {AutoMarker(first_ord_eq)}:"
    
    @property
    def middle_text(self):
        
        system = self._system
        
        first_ord_approx_eq = Eq(system.eoms_approximation_list()[0].lhs-system.eoms_approximation_list()[0].rhs,0)

        return f"Or transformed to the second order ODE {AutoMarker(first_ord_approx_eq)}:"

    @property
    def footer_text(self):
        
        system = self._system
        t_list = system.t_list

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}."

    def append_elements(self):

        system = self._system

        t_list = system.t_list
        first_ord_approx = system.eoms_approximation_list()[1]
        first_ord_approx_eq = Eq(first_ord_approx.lhs-first_ord_approx.rhs,0,evaluate=False)

        display(ReportText(  self.header_text   ))

        display(SympyFormula(first_ord_approx))

        display(ReportText(  self.middle_text   ))

        display(SympyFormula(first_ord_approx_eq))

        display(ReportText(  self.footer_text   ))


class PredictedSolutionComponent(ReportComponent):
    
    title="Predicted solution equation"
    
    
    @property
    def header_text(self):

        return f"Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the consecutive approximations (depending on the accuracy assumed) have the following form:"

    @property
    def footer_text(self):

        return "Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth-order approximation {approx_fun} has the following form:"

    def append_elements(self):

        system = self._system

        t_list = system.t_list
        
        display(ReportText(  self.header_text   ))

        for ord in [0,1]:
            display(SympyFormula(Eq(system.dvars[0],system.predicted_solution(ord)[0])))

#         display(ReportText(  self.footer_text   ))


class GeneralSolutionComponent(ReportComponent):
    
    title="Zeroth-order approximated solution"
    
    
    @property
    def header_text(self):
        
        system = self._system

        t_list = system.t_list

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

        
    @property
    def footer_text(self):

        return f"In order to determine the functions \(\operatorname{C_1}\left(t_{1}\right),\operatorname{C_2}\left(t_{1}\right)\) and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):

        system = self._system
        
        display(ReportText(  self.header_text   ))

        for no,eq in enumerate(system._general_sol(1)):
            display(SympyFormula(Eq(system._general_sol(1)[no].lhs[0],system._general_sol(1)[no].rhs[0])))

#         display(ReportText(  self.footer_text   ))
        
class SecularTermsEquationsComponent(ReportComponent):
    
    title="Zeroth-order approximated solution"
    
    
    @property
    def header_text(self):
        t_list = self._system.t_list
        
        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

        
    @property
    def footer_text(self):

        return "In order to determine the functions \(\operatorname{C_1}\left(t_{1}\right),\operatorname{C_2}\left(t_{1}\right)\) and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):


        system = self._system
        
        display(ReportText(  self.header_text   ))
        sec_eqns = system.secular_eq[system.eps].as_first_ode_linear_system().lhs - system.secular_eq[system.eps].as_first_ode_linear_system().rhs

        for no,eq in enumerate(sec_eqns):
            
            display(SympyFormula(Eq(eq,0)))
            
        display(ReportText(  self.footer_text   ))
        
        
class HomEquationComponent(ReportComponent):
    
    title="Homogenous equation"
    @property
    def header_text(self):

        return "Homogenous equation is the equation in which all the external forces and excitations are assumed to be equal to zero."

        
    @property
    def footer_text(self):

        return "Homogenous equation is used to find the general solution of differential equation which allows to analyse free vibrations of the system."

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text   ))

        display(SympyFormula( system._hom_equation()))
        display(ReportText('Equations presented above correspond to each generalized coordinate:'))
        display(SympyFormula(system.dvars))

        display(ReportText(  self.footer_text   ))

class HomEquationCodeComponent(ReportComponent):
    
    title="Code responsible for calculating homogenous equation"
    @property
    def header_text(self):

        return "Code line responsible for extracting homogenous equation from ODESystsem is following."

        
    @property
    def footer_text(self):

        return "Homogenous equation is used to find the general solution of differential equation which allows to analyse free vibrations of the system."

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text   ))

        display(Markdown
('''
    def _hom_equation(self):
        free = (self.lhs - self.rhs).subs({var:0 for var in self.dvars}).doit()
        hom_eq=self.lhs - free
        return Eq(hom_eq, zeros(len(self.dvars),1))
 
    system = self._system
    display(SympyFormula( system._hom_equation()))
'''))

        display(SympyFormula( system._hom_equation()))
        display(ReportText('Equations presented above correspond to each generalized coordinate:'))
        display(ReportText('Code responsible for displaying generalized coordinates is following:'))
        display(Markdown
('''
        display(SympyFormula(system.dvars))
'''))
        display(SympyFormula(system.dvars))

        display(ReportText(  self.footer_text   ))
        
class ODENonHomogeneousEquationComponent(ReportComponent):
    
    title="Differential equations"
    @property
    def header_text(self):

        return "By separating the equation we obtain the form of equation with time dependent variables on one hand side and free terms on the other:"

        
    @property
    def footer_text(self):

        return "To solve the problem serve several methods depending on the equation type."

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text   ))

        display(SympyFormula( Eq(system._hom_equation().lhs, system._free_component())))

        display(ReportText(  self.footer_text   ))