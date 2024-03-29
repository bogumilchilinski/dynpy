from  ..mechanics import *



class ReportComponent(Subsection):

    latex_name = 'subsection'
    packages=[
              Package('standalone'),
              Package('siunitx')
             ]
    
    title='Report generic component'

    def __init__(self, reported_object, title=None, numbering=False, *, label=True, **kwargs):
        """
        Args
        ----
        title: str
            The section title.
        numbering: bool
            Add a number before the section title.
        label: Label or bool or str
            Can set a label manually or use a boolean to set
            preference between automatic or no label
        """



        self.reported_object = reported_object #it's forced by pylatex library
        
        
        if title is None:
            title = self.title
        
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)
        CurrentContainer(self)
        

        self.append_elements()
        
    def append_elements(self):
        pass

    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame

    @property
    def reported_object(self):

        from ....solvers.linear import ODESystem
        
        if isinstance(self._reported_object, ODESystem):
            return self._reported_object
        else:
            return self._reported_object._ode_system
            
    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj

    @property
    def _system(self):
        print('Kod do poprawienia #################### bo stosujesz starą zmienną self._system')
        print('zamień se self._system na self.reported_object')
        
        return self.reported_object


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

        system = self.reported_object


        display(ReportText(  self.header_text   ))

        for no,eq in enumerate(system):
            display(SympyFormula(Eq(system.lhs[no],system.rhs[no])))

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

        system = self.reported_object


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

        system = self.reported_object

        display(ReportText(  self.header_text   ))

        display(SympyFormula(  system.ivar))
        display(SympyFormula(  system.dvars[0]))

        display(ReportText(  self.footer_text   ))



class GoverningEquationComponent(ReportComponent):
    
    title="Equation of motion"
    
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
            display(SympyFormula(Eq(ode_expr,0)))

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

        system = self.reported_object

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
        
        system = self.reported_object

        Delta = Symbol('Delta')

        return f"The characteristic polynomial ${latex(Delta)}$ of the considered system derived based on its fundamental matrix is presented in {AutoMarker(Eq(Delta,system.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))}:"
    
    @property
    def footer_text(self):

        return "The system fundamental matrix and its characteristic equation allows for a determination of the system unsteady state and its eigenfrequencies related to free vibration."
    
    def append_elements(self):

        system = self.reported_object

        Delta = Symbol('Delta')
        
        display(ReportText(self.header_text))

        display((SympyFormula(Eq(Symbol('A'),system.fundamental_matrix(),evaluate=False))))

        display(ReportText(self.body_text))

        display((SympyFormula(Eq(Delta,system.fundamental_matrix().det().expand().simplify().expand(),evaluate=False))))

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
        system = self.reported_object
        delta=Symbol('Delta')

        display(ReportText(  self.header_text   ))
        polynomial=system.char_polynomial()
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
        system = self.reported_object
        
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
        
        
        system = self.reported_object.set_order(1)
#         approx_sys = system.approximated(3)

        zeroth_ord_eq = Symbol('Omega')#system.nth_eoms_approximation(0)
        eps = Symbol('\\varepsilon')

        return "The ordering and separation of equation {AutoMarker(system.odes[0])} in terms of the power of a small parameter {eps} leads to obtaining a recursive sequence of linear equations of motion leading to the solution of the nonlinear equation. The zeroth-order approximate linear equation is given in {AutoMarker(zeroth_ord_eq)} where the dependency {approx_fun} was assumed for the zeroth-order solution of the time variable:"
    
    @property
    def middle_text(self):
        
        system = self.reported_object.set_order(1)

        zeroth_ord_approx_eq = Symbol('Omega')#Eq(system.nth_eoms_approximation(0).lhs[1]-system.nth_eoms_approximation(0).rhs[1],0)

        return f"Or transformed to the second order ODE {AutoMarker(zeroth_ord_approx_eq)}:"

        
    @property
    def footer_text(self):
        
        system = self.reported_object
        t_list = system.t_list

        return f"Since {t_list[0]} and {t_list[1]} are treated as independent, the differential equation becomes a partial differential equation for a function of two variables {t_list[0]} and {t_list[1]}. Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}."

    def append_elements(self):
        
        system = self.reported_object.set_order(1)
        

        t_list = system.t_list
        zeroth_ord_approx = system.eoms_approximation_list()[0]
        zeroth_ord_approx_eq = Eq((zeroth_ord_approx.lhs-zeroth_ord_approx.rhs)[1],0,evaluate=False)

        display(ReportText(  self.header_text   ))

        for no,eq in enumerate(zeroth_ord_approx):
            
            eq_to_check=Eq(zeroth_ord_approx.lhs[no],zeroth_ord_approx.rhs[no])
            
            if isinstance(eq_to_check,Eq):
                display(SympyFormula(eq_to_check))

        display(ReportText(  self.middle_text   ))

        display(SympyFormula(zeroth_ord_approx_eq))

        display(ReportText(  self.footer_text   ))
        
class FirstOrderApproximatedEqComponent(ReportComponent):
    
    title="First-order approximation"
    
    
    @property
    def header_text(self):
        
        system = self.reported_object.set_order(1)
        first_ord_eq = system.eoms_approximation_list()[1]

        return f"The next component of a recursive sequence of linear equations of motion leading to the solution of the considered nonlinear one is given in {AutoMarker(first_ord_eq)}:"
    
    @property
    def middle_text(self):
        
        system = self.reported_object.set_order(1)
        
        first_ord_approx_eq = Eq(system.eoms_approximation_list()[0].lhs-system.eoms_approximation_list()[0].rhs,0)

        return f"Or transformed to the second order ODE {AutoMarker(first_ord_approx_eq)}:"

    @property
    def footer_text(self):
        
        system = self.reported_object.set_order(1)
        t_list = system.t_list

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}."

    def append_elements(self):

        system = self.reported_object.set_order(1)

        t_list = system.t_list
        first_ord_approx = system.eoms_approximation_list()[1]
        first_ord_approx_eq = Eq((first_ord_approx.lhs[1]-first_ord_approx.rhs[1]),0,evaluate=False)

        display(ReportText(  self.header_text   ))

        for no,eq in enumerate(first_ord_approx):
            
            eq_to_check=Eq(first_ord_approx.lhs[no],first_ord_approx.rhs[no])
            
            if isinstance(eq_to_check,Eq):
                display(SympyFormula(eq_to_check))

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

        system = self.reported_object

#         t_list = system.t_list
        
        display(ReportText(  self.header_text   ))

        for ord in [0,1]:
            display(SympyFormula(Eq(system.dvars[0],system.predicted_solution(ord)[0])))

#         display(ReportText(  self.footer_text   ))


class GeneralSolutionComponent(ReportComponent):
    
    title="Zeroth-order approximated solution"
    
    
    @property
    def header_text(self):
        
        system = self.reported_object

        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of t. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"


    @property
    def footer_text(self):

        return f"In order to determine the functions \(\operatorname{C_1}\left(t_{1}\right),\operatorname{C_2}\left(t_{1}\right)\) and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):

        system = self.reported_object
        
        display(ReportText(  self.header_text   ))
        
#         display(SympyFormula(system.general_solution))

        for no,eq in enumerate(system._general_solution):
            display(SympyFormula(Eq(system._general_solution.lhs[no],system._general_solution.rhs[no])))

#         display(ReportText(  self.footer_text   ))
        
class SecularTermsEquationsComponent(ReportComponent):
    
    title="Zeroth-order approximated solution"
    
    
    @property
    def header_text(self):
        t_list = self.reported_object.t_list
        
        return f"Therefore the general solution may be obtained from the general solution of the corresponding ordinary differential equation by the assumptions of the arbitrary constants becoming the arbitrary functions of {t_list[1]}. Thus solving the considered equation for the unformulated initial conditions, it can be assumed that the predicted solution for the zeroth and first order approximations are as follows:"

        
    @property
    def footer_text(self):

        return "In order to determine the functions \(\operatorname{C_1}\left(t_{1}\right),\operatorname{C_2}\left(t_{1}\right)\) and hence , the first-order approximate equation has to be considered:"

    def append_elements(self):


        system = self.reported_object
        system.order = 1
        
        display(ReportText(  self.header_text   ))
        
        
        system.general_solution
        #display(system.secular_eq)
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

        system = self.reported_object

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

        system = self.reported_object

        display(ReportText(  self.header_text   ))

        display(Markdown
('''
    from sympy import *
    from dynpy.solvers.linear import *
    
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
        
        
class ODEInitComponent(ReportComponent):

    title="Exemplary system"
    @property
    def header_text(self):

        return "Exemplary ODE system can be created using following code."


    @property
    def footer_text(self):

        return "Code presented above shows the elements which are neccessary to build ODESystem."

    def append_elements(self):

        system = self.reported_object
        
        if len(system.dvars)<2:
            display(ReportText('Differential equation to initialize class is as follows:'))
        else:
            display(ReportText('Differential equations to initialize class are as follows:'))
        for n in system.odes:
            display(SympyFormula( n ))

        display(ReportText('Depenedent variables of the presented case has a following form:'))
        display(SympyFormula( system.dvars ))


        display(ReportText(  self.header_text   ))



        code_list = python(system.lhs).split('\n')
        var_name = 'my_eq'
        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]
        dvars=str(system.dvars)
        order=str(system.ode_order)
        if len(system.dvars)<2:
            display(Markdown('''Code presented below describes all of the symbols needed to construct ordinary differential equation. The `my_eq`, `funcs` and `ode_order` variables store the differential equation, unknown functions and equation order.'''))
        else:display(Markdown('''Code presented below describes all of the symbols needed to construct ordinary differential equations. The `my_eq`, `funcs` and `ode_order` variables store the differential equations, unknown functions and equation order.'''))
        display(Markdown(
f'''

    from sympy import *
    from dynpy.solvers.linear import ODESystem

    {code}

    order = {str(system.ode_order)}
    funcs = {str(system.dvars)}

    odesys = ODESystem(odes = my_eq, dvars = funcs, ode_order = order)
    display(odesys)
'''))
        if len(system.dvars)<2:
            display(ReportText("After creating our differential equation it can be analysed using ODESystem class. Additionally one has to provide dvars and ode_order arguments. Dvars is the generalized coordinate of the system and ode_order is the order of considered differential equation."))
        else:
             display(ReportText("After creating our differential equations they can be analysed using ODESystem class. Additionally one has to provide dvars and ode_order arguments. Dvars is the generalized coordinate of the system and ode_order is the order of considered differential equations."))
        display(SympyFormula( system ))

        display(ReportText(  self.footer_text   ))

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
    title="The basic knowledge necessary to solve an equation with separated variables."
    def append_elements(self):
        system = self.reported_object
        ode_sys =system

        display(ReportText(f'One of the fundamental types of differential equations is an equation with separated variables. The process of solving this type of equation involves separating the terms containing $y$ and $x$, moving them to different sides of the equation, and then applying integration. After integrating, its important not to forget to add the constant of integration to one side. Below is an example problem illustrating the process of solving this type of equation.'))


class EquationDefinitionComponent(ReportComponent):
    title="We consider the following differential equation:"

    def append_elements(self):

        system = self.reported_object
        ode_sys =system

        y_sym = system._variable_name()[2]
        display(SympyFormula(Eq(ode_sys.lhs[0],ode_sys.rhs[0])))
        display(ReportText('Where:'))
        display(SympyFormula(y_sym))
        display(ReportText('Is our encountered function.'))

class VariablesSeparationComponent(ReportComponent):
    
    title="Separation of variables"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._get_dvars_symbols()[fun.diff(ivar)]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._get_dvars_symbols()[fun]
        gx=Symbol('g(x)')
        hx=Symbol('h(x)')

        
        
        #ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        
        #sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))
        
        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        display(ReportText(f'If possible, we perform the commonly known separation of variables $x$ and $y$. By separating variables, we mean separating terms containing the variable $x$ from those containing the variable $y$. To accomplish this separation, we can assign these terms to functions $h(x)$ and $g(x)$, as shown below:'))
        display(SympyFormula(Eq(gx,g_fun_ivar)))
        display(SympyFormula(Eq(hx,h_fun_dvar)))
        display(ReportText(f'Apart from these elements, our equation also contains the term ${latex(fun_str)}$, which can also be written as:'))
        display(SympyFormula(Eq(fun_str,d_var/d_ivar)))
        display(ReportText(f'After substituting this expression, it can be observed that in subsequent steps, multiplying the equation on both sides by ${latex(d_ivar)}$ enables integrating both sides with respect to different variables.'))

class SeparatedVariablesIntegrationComponent(ReportComponent):
    
    #title="Integration of separated variables"
    title="Integration of separated variables"
    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar
        gx=Symbol('g(x)')
        hx=Symbol('h(x)')
        

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]
        y_sym = system._variable_name()[2]
        fun_str = system._variable_name()[3]

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
#         display(ReportText(''))
        
        display(SympyFormula( Eq(Symbol(latex(d_var/d_ivar)),g_fun_ivar*h_fun_dvar)))
        display(ReportText(f'''To do this, we multiply both sides by ${latex(d_ivar)}$ and move ${latex(fun)}$ to the opposite side of the equation:''' ))
        
        display(SympyFormula( Symbol(latex(d_var/h_fun_dvar)+ '=' + latex(g_fun_ivar) + latex(d_ivar) )   ) )
        display(ReportText('''Then, to "go back" to the original function and obtain the result in the form of the f $y$, we integrate both sides of the obtained equation. '''))
        
        display(ReportText('''Importantly, when integrating, it's essential to remember to add the constant of integration in any form, typically denoted as $C$ or $D$. '''))
        display(SympyFormula(  Symbol(latex( Integral(1/h_fun_dvar,y_sym) )+ '= ' + latex(Integral(g_fun_ivar,ivar))  )   ) )

#         const = (system.solution._spot_constant())[0]
        const=Symbol('D')
        display(ReportText('After integration, we obtain the following result:'))
        display(SympyFormula(  Symbol(latex( Integral(1/h_fun_dvar,y_sym).doit() )+ '=' + latex(Integral(g_fun_ivar,ivar).doit()) + '+' + latex(const)   )   ) )

#         display(ReportText(  self.footer_text   ))

class SolutionComponent(ReportComponent):
    
    title="Simplification of the solution obtained after integrating the equations"

    def append_elements(self):

        system = self.reported_object
        ode_sys =system
        ############################ must have

#         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]


        
        ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]
        
        sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))
        
        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        
        display(ReportText('After rearrangment a solution has the following form:'  ))
        display(SympyFormula( system.general_solution.as_eq_list()[0] ))

        display(ReportText('The obtained expression represents the solution to a differential equation with separated variables. The constant of integration can be determined when initial conditions are provided in the problem, i.e., the value of the function for a specific $x$.'))
        
class LinearODEIntroComponent(ReportComponent):
    title="The basic knowledge necessary to solve a differential equation"
    
    def append_elements(self):

        system = self.reported_object
        ode_sys =system
        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]
        Px=Symbol('P(x)')
        Qx=Symbol('Q(x)')
        n=Symbol('n')
        display(ReportText('One of the most fundamental types of differential equations is a linear equation with constant coefficients. A linear equation consists only of linear terms, which may include constant coefficients. Transitioning to a homogeneous equation is required - we move all terms involving the function $y$ and its derivatives, then equate them to zero, obtaining a homogeneous equation. Next, it is necessary to treat the constant as a function dependent on $x$. We calculate the derivatives of both sides and utilize the results in the initial equation. The result is the sum of solutions from both parts.'))

class LinearTransformation(ReportComponent):
    title="Transition to a homogeneous equation"
    
    def append_elements(self):
        system = self.reported_object
        ode_sys =system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        fun_str = system._variable_name()[3] #y'

        d_ivar = system._variable_name()[0] #dx
        d_var = system._variable_name()[1] #dy

        y_sym = system._variable_name()[2] #y(x)
        const=Symbol('D')

        display(ReportText('We move all terms involving the function $y$ and its derivatives to one side and equate it to zero, obtaining a homogeneous equation.'))
        #display(SympyFormula(Eq(system._hom_equation().lhs[0],0)))
        
        
        


class LinearToSeparable(ReportComponent):
    
    title="The solution to the resulting separated variable equation"
    def append_elements(self):

        system = self.reported_object
        fun = system.dvars[0]
        
        ode_sys=system._hom_equation()
        #sep_ode_report=ode_sys.report #mniej więcej coś takiego tylko coś pomyliłem na 100%
        
        display(ReportText('The obtained homogeneous equation is a separated variable equation and takes the form:'))
        display(SympyFormula(Eq(ode_sys.lhs[0],0)))
        
        #print("#"*200)
        sep_ode_report=ode_sys.report
        
        
        #print(list(sep_ode_report))
        #print('#'*3)
        #print(list(sep_ode_report)[3:])
        
        for elem in list(sep_ode_report)[3:]:
            self.append(elem)



class VariationOfConstant(ReportComponent):
    title="Variation of constant"
    
    def append_elements(self):

        system = self.reported_object
        ode_sys =system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar
        

        ode_sys2 = ode_sys
        x=ivar
        szur=ode_sys2.general_solution
#         display(szur)
        ode_syso=Eq(ode_sys2.lhs[0],ode_sys2.rhs[0])
        Cx = Function(szur._spot_constant()[0])(x)
        Cx
        c_prim=Function("C'")(x)
        szur2=szur.subs(szur._spot_constant()[0], Cx)
        #display(szur2)
        szur3=Eq(szur2.lhs.diff(x), szur2.rhs.diff(x))
        #display(szur3)
        szur3_=szur3.subs(Cx.diff(x),c_prim)
        szur3_
        #To implement, change left hand side to y_prim variable and subsitute variables with object variables

        ##### Tu sie teoretycznie zaczyna variation of constant ########################################
        sol_cx=(Eq(Symbol('y'),szur2.rhs[0]))
        display(ReportText('The next step is to treat the constant as a function dependent on $x$:'))
        display(SympyFormula(sol_cx))
        sol_diff3=Eq(Symbol("y'"), szur3_.rhs[0])
        display(ReportText('Next, you need to calculate the derivatives of both sides of the equation.'))
        ## dodać równanie 
        display(ReportText('Upon obtaining the equation in the following form:'))
        display(SympyFormula(sol_diff3))#gdzie? gdzie co?
        display(ReportText("Thanks to this operation, we can now substitute $y$ and $y'$ into the original equation."))
        display(SympyFormula(ode_syso))
        lolo=ode_syso.subs(y.diff(ivar),sol_diff3.rhs).subs(y,sol_cx.rhs)
        display(ReportText('And obtain:'))
        display(SympyFormula(lolo))
        c_p_sol=solve(lolo,c_prim)
        cp_eq=Eq(c_prim,c_p_sol[0])
        display(ReportText("Then, if everything has been calculated correctly, the values of $C(x)$ should simplify to a form from which we can easily determine $C'(x)$, obtaining:"))
        display(SympyFormula(cp_eq))
        cx_wyn=Eq(Cx,cp_eq.rhs.diff(ivar))
        display(ReportText("After integrating both sides, we obtain the value of the constant:"))
        display(SympyFormula(cx_wyn))
        FINALOWY_WYNIK=sol_cx.subs(Cx,cx_wyn.rhs)
        display(ReportText("To obtain the final result, we need to substitute the obtained constant into the first equation where we introduced the constant as a variable:"))
        display(SympyFormula(FINALOWY_WYNIK))



class BernoulliODEIntroComponent(ReportComponent):
    title="The basic knowledge necessary to solve a differential equation:"
    
    def append_elements(self):

        system = self.reported_object
        y_sym = system.dvars[0]
        dane=system._get_dvars_symbols()
        Px=Symbol('P(x)')
        Qx=Symbol('Q(x)')
        n=Symbol('n')
        ynpow=Symbol('y')**Symbol('n')
        display(ReportText('The Bernoulli differential equation is a special type of first-order ordinary differential equation. Its general form is:'))
        bernoulli_form=Eq(y_sym.diff()+Px*y_sym,Qx*y_sym**n)
        display(SympyFormula(bernoulli_form.subs(dane)))
        
        display(ReportText(f'Where $y$ is the unknown function dependent on the variable $x$, $P(x)$ and $Q(x)$ are given functions dependent on $x$, and $n$ is a constant real number, different from 0 and 1. Characteristic of Bernoulli equations is that they are nonlinear due to the presence of the term ${latex(ynpow)}$ on one side of the equation. Solutions to such equations can be challenging to obtain in general, but for certain values of $n$, certain techniques can be applied to simplify the solution. In the case of $n=0$, the Bernoulli equation becomes linear, and for $n=1$, it can be transformed into the form of a linear ordinary differential equation of the first order. In other cases, special techniques such as variable substitution or order reduction may be necessary to obtain the general solution.'))

class EquationDefinitionComponent(ReportComponent):
    title="We consider the following differential equation:"

    def append_elements(self):

        system = self.reported_object
        ode_sys =system
        y_sym = system.dvars[0]

        dane=system._get_dvars_symbols()
        ode_sub=ode_sys.subs(dane)
        display(SympyFormula(Eq(ode_sub.lhs[0],ode_sub.rhs[0])))
        display(ReportText('Where:'))
        display(SympyFormula(y_sym))
        display(ReportText('Is our encountered function.'))
        

class BernoulliTransformation(ReportComponent):
    title="Determining the new function $z$ dependent on $x$ in order to transition to a linear equation."
    
    def append_elements(self):
        system = self.reported_object
        ode_sys =system
        fun = system.dvars[0]
        ivar = system.ivar
###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2] #y'

        d_ivar = system._variable_name()[0] #dx
        d_var = system._variable_name()[1] #dy
        dane=system._get_dvars_symbols()
        bernoulli_sub=system._bernoulli_substitution()
        #bern_power = system



        
        n=Symbol('n')
        z=Function('z')(ivar)

        display(ReportText('To make the substitution, it is necessary to define $n$ by finding the highest power of $y$, where $n$ is the exponent of this function.'))
        display(SympyFormula(Eq(Symbol('n'),n.subs(bernoulli_sub))))
        display(ReportText('Next, you need to make the substitution and transition to the new function $z$, using the following formula:'))
        display(SympyFormula(Eq(Symbol('z'),Symbol('y')**(S.One-Symbol('n')))))
        display(ReportText('In our case:'))
        display(SympyFormula(Eq(Symbol('z'),Symbol('y')**(S.One-n.subs(bernoulli_sub)))))
        #display(SympyFormula(Eq(Symbol('z'),z.subs(bernoulli_sub))))
        display(ReportText('The next step is to determine the first derivative of the function $z$:'))
        display(SympyFormula(Eq(Symbol("z'"), bernoulli_sub['z_diff'])))
        display(ReportText('After substituting the $n$ value :'))
        display(SympyFormula(Eq(Symbol("z'"), bernoulli_sub['z_diff'].subs(n,bernoulli_sub[n]))))
        
        display(ReportText("Next, we need to determine $y$ and $y'$ from the above equations."))
        display(SympyFormula(Eq(Symbol("y"),fun.subs(bernoulli_sub))))
        display(SympyFormula(Eq(Symbol("y'"),fun.diff().subs(bernoulli_sub))))

class BernoulliLinearTransformation(ReportComponent):
    title=" Transition to a linear equation"
    
    def append_elements(self):
        system = self.reported_object
        ode_sys =system
        fun = system.dvars[0]
        ivar = system.ivar
###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2] #y'

        d_ivar = system._variable_name()[0] #dx
        d_var = system._variable_name()[1] #dy
        dane=system._get_dvars_symbols()

        rownanie=system._as_reduced_ode()

        display(ReportText("The last step in transitioning to a linear equation is substituting the previously calculated $y$ and $y'$ into the original equation:"))
        display(SympyFormula(Eq(rownanie.lhs[0],rownanie.rhs[0])))
        display(ReportText('We compute the obtained equation in the following manner:'))

        sep_ode_report=rownanie.report


        for elem in list(sep_ode_report)[3:]:
            self.append(elem)
