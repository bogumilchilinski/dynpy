from  ..mechanics import *
from sympy.physics.mechanics import dynamicsymbols


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

            self.append(NoEscape('\\Huge Mechanical Vibration \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))

            
            
            if len(system.q)==1:
                dof_str = 'Single Degree of Freedom'
            else:
                dof_str = 'Multi Degree of Freedom'
                
            if system._dissipative_potential==0 or system._dissipative_potential is None:
                damping_str = 'UNDAMPED'
            else:
                damping_str = 'DAMPED'

           
            self.append(NoEscape(f'\\Large {dof_str} {damping_str} SYSTEM \n \n'))    
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(NoEscape(f'{system.__class__.__name__} \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(Command('MyAuthor'))
            self.append(NoEscape(f'\\par'))
            self.append(Command('vspace',arguments='1cm'))
            
            #self.append(Command('vspace',arguments='1cm'))
            #self.append(NoEscape(f'\\protect\\par'))
            #self.append(NewLine())
            self.append(Command('MyDate'))


            
######################### Initial point

    
    
# Damian
class ExemplaryPictureComponent(ReportComponent):
    title="Example of an real object"
    packages=[Package('float')] 
    @property
    def header_text(self):
        #"Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."
        return "The picture represents a real object."
    @property
    def footer_text(self):


        return "A real example is shown for better understanding of the schematic view and purpose of the system."
    def append_elements(self):
        
        system = self._system

        display(ReportText(self.header_text))
        

        display(Picture(system._real_example(),width='8cm'))


        display(ReportText(self.footer_text))

class ExemplaryPictureDynPyCodeComponent(ReportComponent):

    title="Dynpy Code for Exemplary Picture"

    @property
    def header_text(self):
        return "Code line responsible for extracting exemplary picture from dynamic system is as following:"


    @property
    def footer_text(self):
        return "Picture representing real object of considered dynamic system."

    def append_elements(self):

        system = self._system
        dyn_sys=system

        display(ReportText( self.header_text  ))


        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *
    from dynpy.utilities.report import Picture

    dyn_sys = {system.__class__.__name__}()

    display(Picture(dyn_sys._real_example()))


'''))

        display(ReportText( self.footer_text  ))
        
        
class ExemplaryPictureSymPyCodeComponent(ExemplaryPictureComponent):
    title="SymPy code for calling scheme of the system"

    @property
    def header_text(self):

        return "Code for calling an exemplary picture is as follows:"

    @property
    def footer_text(self):
        
        system = self._system

        return f"Analysis of the real example is needed to propose physical scheme."
    
    
    def append_elements(self):
        
        system = self.reported_object
        path=system._real_example()
        
        display(ReportText(  self.header_text ))
        
        display(Markdown(
f'''


    from dynpy.models.mechanics import *
    from sympy import *
    from dynpy.utilities.report import Picture
   
    path="{path}"
    
    display(Picture(path ,width='8cm'))



'''
            ))
        display(ReportText( self.footer_text ))
        
        
# Boogi        
class SchemeComponent(ExemplaryPictureComponent):
    title="Scheme of the system"

    @property
    def header_text(self):
        #"Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."
        return "Scheme of the real object is presented in the figure. It was obtained by the analysis of real mechanism."

        
    @property
    def footer_text(self):
        
        system = self._system
        
        #"Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}."
        return f"Analysis of the scheme allows to claim its number of degrees of freedom which is {len(system.q)}"
    
    
    def append_elements(self):
        
        system = self._system

        display(ReportText(  self.header_text ))

        display(Picture(system._scheme(),width='8cm'))
        


        display(ReportText( self.footer_text ))


class SchemeDynPyCodeComponent(ExemplaryPictureComponent):
    title="DynPy code for scheme of the system"

    @property
    def header_text(self):

        return "DynPy code for calling the scheme of the dynamic system is as following:"


    @property
    def footer_text(self):

        system = self._system

        return f"Analysis of the scheme allows to claim its number of degrees of freedom which is {len(system.q)}"

    def append_elements(self):

        system = self._system

        display(ReportText(  self.header_text ))

        display(Markdown(

f'''



    from dynpy.models.mechanics import *
    from dynpy.utilities.report import Picture
    from sympy import *

    dyn_sys = {system.__class__.__name__}()


    display(Picture(dyn_sys._scheme(), caption='Scheme of considered dynamic system.'))



'''))

        display(ReportText( self.footer_text ))


class SchemeSymPyCodeComponent(ExemplaryPictureComponent):
    title="SymPy code for scheme of the system"

    @property
    def header_text(self):
        return "SymPy code for calling the scheme of the dynamic system is as following:"

        
    @property
    def footer_text(self):
        
        system = self._system
        return f"Analysis of the scheme allows to claim its number of degrees of freedom which is {len(system.q)}"
    
    
    def append_elements(self):
        
        system = self._system
        path=system._default_folder_path + system.scheme_name

        display(ReportText(  self.header_text ))

        display(Markdown(
f'''


    from dynpy.models.mechanics import *
    from dynpy.utilities.report import Picture
    from sympy import *

    path="{path}"

    display(Picture(path ,width='8cm', caption='Scheme of considered system.'))



'''
            ))
        display(ReportText( self.footer_text ))


#TODO
class DamperPlot(ExemplaryPictureComponent):
    title="Siła tłumienia"
  

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys




class NumericalAnalysisComponent(ExemplaryPictureComponent):
    title="Symulacja numeryczna"

    _default_ics = None
    _default_parameters = {}
    _default_tspan = np.linspace(0,1,1000)
    
    @classmethod
    def set_default_ics(cls,default_ics=None):
        
        if default_ics is None:
            cls._default_ics = default_ics
        
        return cls

    @classmethod
    def set_default_parameters(cls,default_parameters=None):
        
        if default_parameters is None:
            cls._default_parameters = default_parameters
        
        return cls
    
    @classmethod
    def set_default_tspan(cls,default_tspan=None):
        
        if default_tspan is None:
            cls._default_tspan = default_tspan
        
        return cls
    

    def append_elements(self):
        
        system = self._system

        display(ReportText(f'''Dla Damiana :P
                            '''))
        if self._default_ics is None:
            ics = list(system.Y*0.0)
        else:
            ics = self._default_ics
            
        sym_res = system.subs(self._default_parameters).numerized().compute_solution(self._default_tspan,ics)

        LatexDataFrame.formatted(sym_res).plotted(preview=True)
        

        display(ReportText(f'''Dla Damiana :P
                            '''))

# Boogi
class KineticEnergyComponent(ReportComponent):
    
    title="Kinetic energy"

    @property
    def header_text(self):
        #"Energia kinetyczna układu wyrażona jest wzorem:"
        return "Kinetic energy of the system has a following form:"

        
    @property
    def footer_text(self):
        #"Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
        return "Determined formula specify energy of the system related to its inertial properties."
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys




        display(ReportText( self.header_text  ))
        

        display(SympyFormula( Eq(Symbol('T'),
                     dyn_sys_lin._kinetic_energy) , marker=None))
        
        
        
               
        display(ReportText( self.footer_text  ))
    
class KineticEnergyDynPyCodeComponent(ReportComponent):

    title="Dynpy Code for Kinetic Energy"

    @property
    def header_text(self):
        return "Code line responsible for extracting kinetic energy from dynamic system is as following:"

    @property
    def footer_text(self):
        return "Determined piece of code specify energy of the system related to its inertial properties."

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys



        display(ReportText( self.header_text  ))


        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *

    dyn_sys = {system.__class__.__name__}()

    display( Eq(Symbol('T'),
                     dyn_sys._kinetic_energy) )


'''))
        
        
class KineticEnergySymPyCodeComponent(ReportComponent):

    title="Sympy Code for Kinetic Energy"

    @property
    def header_text(self):
        return "Code line responsible for providing symbolic form of kinetic energy is as following:"


    @property
    def footer_text(self):
        return "Determined piece of code specify energy of the system related to its inertial properties."

    def append_elements(self):

        system = self._system
        dyn_sys=system


        display(ReportText( self.header_text  ))

#         Ek=dyn_sys._kinetic_energy

        code_list=python(system._kinetic_energy).split('\n')

        var_name = 'Ek'

        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]

        display(Markdown(
f'''




    from sympy import *

    {code}


    display(Eq(Symbol('T'),Ek))




'''))

        display(ReportText( self.footer_text  ))

    
class PotentialEnergyComponent(ReportComponent):#Jaś fasola
    
    title="Potential energy"
    @property
    def header_text(self):
        #"Energia potencjalna układu wyrażona jest wzorem:"
        return "Potential energy of the system has a following form:"

        
    @property
    def footer_text(self):
        #"Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."
        return "The presented relationship describes the interaction of potential force fields in which the object is located."

    def append_elements(self):

        system = self._system
        dyn_sys=system


        display(ReportText(  self.header_text   ))

        display(SympyFormula( Eq(Symbol('V'),
                     dyn_sys._potential_energy ), marker=None))

        display(ReportText(  self.footer_text   ))

class PotentialEnergyDynPyCodeComponent(ReportComponent):
    
    title="Dynpy Code for Potential Energy"

    @property
    def header_text(self):
        return "Code line responsible for extracting potential energy from dynamic system is as following:"

        
    @property
    def footer_text(self):
        return "Determined piece of code specify energy of the system, which is stored in object."
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system


        display(ReportText( self.header_text  ))


        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *
    
    dyn_sys = {system.__class__.__name__}()

    display( Eq(Symbol('Ep'),
                     dyn_sys._potential_energy) )

    
'''))

        display(ReportText( self.footer_text  ))
        
        
class PotentialEnergySymPyCodeComponent(ReportComponent):

    title="Sympy Code for Potential Energy"

    @property
    def header_text(self):
        return "Code line responsible for providing symbolic form of potential energy is as following:"


    @property
    def footer_text(self):
        return "Determined piece of code specify energy of the system, which is stored in object."

    def append_elements(self):

        system = self._system
        dyn_sys=system


        display(ReportText( self.header_text  ))

        Ep=dyn_sys._potential_energy

        code_list=python(system._potential_energy).split('\n')

        var_name = 'Ep'

        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]

        display(Markdown(
f'''




    from sympy import *

    {code}


    display(Eq(Symbol('E_p'),Ep))




'''))

        display(ReportText( self.footer_text  ))


# Marcel
class DissipationComponent(ReportComponent):
    
    #"Dyssypacyjna funkcja Rayleigh'a"
    title="Dissipative Rayleigh function"
    
    
    @property
    def header_text(self):
        #"Energia rozpraszana tłumieniem wyrażona jest wzorem:"
        return "The energy dissipated by attenuation is given by the formula:"

        
    @property
    def footer_text(self):
        #" Podana zależność stanowi potencjał dysynpacyjny Rayleigh'a, który poddany różniczkowaniu względem wektora prędkości uogólnionych pozwala na określenie sił wiskotycznego tłumienia."
        return "The given dependence is the Rayleigh dissipation potential, which when differentiated against the generalized velocity vector allows to determine the viscous damping forces."
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

        


        display(ReportText( self.header_text  ))
        

        display(SympyFormula( Eq(Symbol('D'),
             dyn_sys_lin._dissipative_potential) , marker=None ))

        display(ReportText( self.footer_text  ))

# Boogi          
class LagrangianComponent(ReportComponent):
    
    #title="Lagrangian (funkcja Lagrange'a)  układu"
    title="Lagrangian of the system (Lagranges function)"
        
    @property
    def header_text(self):
        #"Lagrangian systemu wyrażony jest wzorem ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):"
        
        system = self._system
        dyn_sys=system
        
        return f"System Lagrangian is described by the formula ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):"
        
    @property
    def body_text_1(self):
        #"Równania Eulera-Lagrange'a dla rozważanego przypadku są nastęujące: "   
        return "The Euler-Lagrange equations for the case under consideration are as follows: "
        
    @property
    def body_text_2(self):
        #"Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: "
        return "Subsequent derivatives obtained with the Euler-Lagrange equations are as follows: "
    
    @property
    def footer_text(self):
        #"Wyniki przedstawionych operacji wykorzystuje się wyznaczenia równań ruchu układu."
        return "The results of the presented operations are used to determine the equations of motion of the system."
        
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        
        print(system._scheme())
        
        mrk_lagrangian_nonlin=Marker('lagrangianNL',prefix='eq')

        #display(ReportText(f'''The following model is considered. The system's Lagrangian is described by the formula ({Ref(mrk_lagrangian_nonlin).dumps()}):
        #                    '''))
        display(ReportText(self.header_text ))

        display((SympyFormula(  Eq(Symbol('L'),dyn_sys.L.expand()[0])  , marker=mrk_lagrangian_nonlin )  ))
        
        q_sym =[ Symbol(f'{coord}'[0:-3]) for coord in dyn_sys.q]
        
        diffL_d=lambda coord: Symbol(f'\\frac{{ \\partial L}}{{  \\partial {vlatex(coord)}  }}')
        diffD_d=lambda coord: Symbol(f'\\frac{{ \\partial D}}{{  \\partial {vlatex(diff(coord))}  }}')
        d_dt_diffL_d=lambda coord: Symbol(f'\\frac{{ \\mathrm{{d}}  }}{{  \\mathrm{{d}} {vlatex(system.ivar)}  }} {vlatex(diffL_d(coord))} ')        

        display(ReportText(self.body_text_1))
        
        for coord in dyn_sys.q:
            display((SympyFormula(  Eq(d_dt_diffL_d(coord.diff(system.ivar)) - diffL_d(coord) + diffD_d(coord),Symbol(f'Q_{{ {vlatex(coord)} }}^N'))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        
        
        display(ReportText(self.body_text_2))
        
        for coord in dyn_sys.Y:
            display((SympyFormula(  Eq(diffL_d(coord),dyn_sys.L.expand()[0].diff(coord))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            


        for coord in dyn_sys.q.diff(system.ivar):
            display((SympyFormula(  Eq(d_dt_diffL_d(coord),dyn_sys.L.expand()[0].diff(coord).diff(system.ivar))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        for coord in dyn_sys.q:
            display((SympyFormula(  Eq(diffD_d(coord),(S.One/2 * diff(dyn_sys.q.transpose())* dyn_sys_lin.damping_matrix()* diff(dyn_sys.q))[0].diff(diff(coord)))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            #display(Markdown(f'\\begin{equation}    \\end{equation}').reported())
        #with doc_model.create(DMath()) as eq:
        #    eq.append(NoEscape(latex(Derivative(Symbol('L'),q_sym[0],evaluate=False))))
        #    eq.append(NoEscape('='))
        #    eq.append(NoEscape(vlatex(dyn_sys.L.expand()[0].diff(dyn_sys.q[0]))))

        mrk_gov_eq_nonlin=Marker('gov_eq_nonlin_sys',prefix='eq')

        #display(ReportText(f'''The governing equations of the system have a following form ({Ref(mrk_gov_eq_nonlin).dumps()}):
        #                    '''))

        display(ReportText(self.footer_text))


        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame

        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame


class LagrangianDynPyCodeComponent(ReportComponent):
    
    title="Dynpy Code for Lagrangian Component"

    @property
    def header_text(self):
        return "Code line responsible for providing symbolic form of Lagrange's function is as following:"

        
    @property
    def footer_text(self):
        return "Considered piece of code can be used to determine the equations of motion of a system."
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system


        display(ReportText( self.header_text  ))


        display(Markdown(
f'''



    from dynpy.models.mechanics import *
    from sympy import *

    dyn_sys = {system.__class__.__name__}()
    t = dyn_sys.ivar

    display( Eq(Symbol('L'),
                     dyn_sys.lagrangian()) )

    for q in dyn_sys.q:

        display(dyn_sys.lagrangian().diff(q))

        display(dyn_sys.lagrangian().diff(q.diff(t)).diff(t))



'''))

        display(ReportText( self.footer_text  ))


class LagrangianSymPyCodeComponent(ReportComponent):

    title="Sympy Code for Lagrangian Component"

    @property
    def header_text(self):
        return "Code line responsible for providing symbolic form of Lagrange's function is as following:"

        
    @property
    def footer_text(self):
        return "Considered piece of code can be used to determine the equations of motion of a system."
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system


        display(ReportText( self.header_text  ))


        L=dyn_sys.lagrangian()
        q=dyn_sys.q
        t=dyn_sys.ivar

        code_list=python(system.lagrangian()).split('\n')

        var_name = 'L'

        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]


        display(Markdown(
f'''




    from sympy import *

    {code}


    display(Eq(Symbol('L'),L))

    qs = {q}

    for q in qs:

        display(L.diff(q))

        display(L.diff(q.diff(t)).diff(t))



'''))

        display(ReportText( self.footer_text  ))

        
        
# Amadi & Damian
class GoverningEquationComponent(ReportComponent):
    
    #Równania ruchu
    title="Equation of motion"
    
    @property
    def entry_text_one(self):
        system = self._system
        dyn_sys=system
        return f"Using the calculated derivatives, the equation of motion is based on the appropriate formula. System equation of motion is described by the formula ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})"
    
    @property
    def entry_text_two(self):
        system = self._system
        dyn_sys=system
        #f"przedstawiają zależności ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"
        return f"Using the calculated derivatives, the equations of motions are based on the appropriate formulas. System equations of motions are described by the formulas ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-     ({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"

    @property
    def footer_text(self):
        #f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return f''' The determined equations constitute a mathematical dynamic description of the properties of the system.
                    Further analysis allows for an effective analysis of the modeled object's operation and determination of its mechanical parameters. '''
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system

        if len(system.q)==1:
            display(ReportText(self.entry_text_one))   
        else:
            display(ReportText(self.entry_text_two))   

        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=None))

        display(ReportText(self.footer_text))   

        
class LinearizedGoverningEquationComponent(ReportComponent):
    #Równania ruchu
    title="Equation of motion"
    
    @property
    def entry_text_one(self):
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()
        #f"przedstawia zależność ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})"
        return f"Using the calculated derivatives, the equation of motion is based on the appropriate formula. Linearized system equation of motion is described by the formula ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})"
    
    @property
    def entry_text_two(self):
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()
        #f"przedstawiają zależności ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"
        return f"Using the calculated derivatives, the equations of motions are based on the appropriate formulas. Linearized system equations of motions are described by the formulas ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys_lin._eoms[-1].simplify().expand(),0))})"
    
    @property
    def footer_text(self):
        #f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return f''' The determined equations constitute a mathematical dynamic description of the properties of the system.
                    Further analysis allows for an effective analysis of the modeled object's operation and determination of its mechanical parameters. '''
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()


        
        if len(system.q)==1:
            display(ReportText(self.entry_text_one))   
        else:
            display(ReportText(self.entry_text_two))   
        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=None))

        display(ReportText(self.footer_text))   
        
        
class GoverningEquationDynpyCodeComponent(ReportComponent):
    
    #Równania ruchu
    title="Equation of motion"
    
    @property
    def header_text(self):
      
      return 'Presented piece of code describes equations of motoion of the system distributed from dynamical system'

       

    @property
    def footer_text(self):
        #f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return 'Execution of this code returns equations of motion generated from dynamic system'
    
    def append_elements(self):
        system = self._system
        dyn_sys=system
        eoms=str(dyn_sys._eoms)
        
        display(ReportText(self.header_text))

        display(Markdown(
f'''
    
    


    from dynpy.models.mechanics import *
    from sympy import *
    
    dyn_sys = {system.__class__.__name__}()
    
    
    display(Eq(dyn_sys._eoms,0,evaluate=False))



    
'''))

        display(ReportText(self.footer_text))        
        
class GoverningEquationSympyCodeComponent(ReportComponent):
    
    #Równania ruchu
    title="Equation of motion"
    
    @property
    def header_text(self):
      
      return 'Presented piece of code describes equations of motoion of the system created using Sympy library'

       

    @property
    def footer_text(self):
        #f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return 'Execution of this code returns equations of motion in a form of Sympy outputs'
    
    def append_elements(self):
        system = self._system
        dyn_sys=system
        eoms=str(dyn_sys._eoms)
        display(ReportText(self.header_text))
        
        code_list=python(system._eoms).split('\n')
        
        var_name = 'eoms'
        
        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]

        

        display(Markdown(
f'''
    
    


    from dynpy.models.mechanics import *
    from sympy import *
    
    {code} 
    
    display(Eq(eoms,0,evaluate=False))

    
'''))

        display(ReportText(self.footer_text))

class LinearizationComponent(ReportComponent): # Szymon

    title="Linearization of equation of motion"

    @property
    def entry_text(self):
        return "Linearization of governing equations is about finding Taylor series with respect to generalized coordinates, velocities and accelerations in the neighbourhood of the equilibrium point. Following symbols have been introduced to make a simplification."

    @property
    def equilibrium_point_text(self):
        return "Equilibrium points of the system have following forms:"

    @property
    def eom_text(self):
        dyn_sys= self._system
        return '''Equation of motion for coordinate ${coord}$ can be presented as:'''

    @property
    def lagrange_text(self):
        dyn_sys= self._system
        return "Proper computings requires finding derivatives of generalized coordinates, which are components of Lagrange's equations"

    @property
    def derivative_text(self):
        return "The calculated derivatives have a following form:"

    @property
    def linearized_eq_text(self):
        return "The following equation (linearized) can be obtained after substitution of calculated derivatives."

    def append_elements(self):

        system = self._system
        ReportText.set_directory('./SDAresults')
        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        mrk_lagrangian_nonlin = Marker('lagrangLin',prefix='eq')
        mrk_lagrangian_lin = Marker('lagrangLin',prefix='eq')

        t=system.ivar

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()
        coords=tuple(list(dyn_sys.Y) + list(dyn_sys.q.diff(t,t)))
        op_point = {coord: 0 for coord in coords}
        op_point.update(dyn_sys._op_points(subs=True)[0])

        display(ReportText(self.entry_text))


        for coord in coords:
            display((SympyFormula(  Eq(Symbol(vlatex(coord)),Symbol(latex(coord))) , marker=None  )  )) 



        display(ReportText( self.equilibrium_point_text ))

        for eq_coord,val in op_point.items():
            display((SympyFormula(  Eq(eq_coord,val) , marker=mrk_lagrangian_lin  )  ))


        diffL_d=lambda coord: Symbol(latex(Derivative(Symbol('L'),Symbol(vlatex(coord))))  )



        for no,eom in enumerate(dyn_sys._eoms):




            eq_sym=Symbol(f'RR_{latex(dyn_sys.q[no])}')


            display(ReportText( self.eom_text.format(coord = latex(dyn_sys.q[no])) ))

            #display((SympyFormula(  Eq(eq_sym,eom,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))


            display(ReportText( self.lagrange_text ))


            display((SympyFormula(  Eq(MultivariableTaylorSeries(eq_sym,coords,n=1,x0=op_point)._symbolic_sum(),0) , marker=None,backend=latex  )  ))

            diff_list=MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).calculation_steps(expr_symbol=eq_sym)

            display(ReportText( self.derivative_text ))

            for diff_eq in diff_list:

                display((SympyFormula(  diff_eq , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display(ReportText( self.linearized_eq_text ))
            display((SympyFormula(  Eq(MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).doit().expand().simplify().expand(),0,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))



        AutoBreak.latex_backend = latex_store

class CriticalPointsComponent(ReportComponent):

    @property
    def header_text(self):
        return ('From the relvant kinematic dependencies the critical points of the system are calculated')
    
    @property
    def footer_text(self):
        return('Resultant values of generalized coordinates provide static equilibrium of the system')
    
    def append_elements(self):

        system = self._system
        dyn_sys = system

        display(ReportText(self.header_text))

        for point in dyn_sys._op_points():
          display([Eq(key,Val) for key, Val in point.items()])

        display(ReportText(self.footer_text))

class CriticalPointsDynPyCodeComponent(ReportComponent):

    @property
    def header_text(self):
        return ('The following piece of code allows calculation of critical points for the given system')
    
    @property
    def footer_text(self):
        return('Outputs of this code can be used for further calculation of non-linear systems')
    
    def append_elements(self):

        system = self._system
        dyn_sys = system

        display(ReportText(self.header_text))
        display(Markdown(
f'''
    
    from dynpy.models.mechanics import *
    from sympy import *
    
    dyn_sys = {system.__class__.__name__}()
    
    display(dyn_sys._op_points())
    
'''))
       
        display(ReportText(self.footer_text))

class CriticalPointsSymPyCodeComponent(ReportComponent):

    @property
    def header_text(self):
        return ('From the relvant kinematic dependencies the critical points of the system are calculated')
    
    @property
    def footer_text(self):
        return('Resultant values of generalized coordinates provide static equilibrium of the system')
    
    def append_elements(self):

        system = self._system
        dyn_sys=system

        critical_points=str(dyn_sys._op_points())

        display(ReportText(self.header_text))

        code_list=python(system._op_points()).split('\n')
        
        var_name = 'critical_points'
        
        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]

        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *
    
    {code} 
    
    display(critical_points)

'''))

        display(ReportText(self.footer_text))

# Marcel & Monika
#class FundamentalMatrixComponent(ReportComponent):
#    
#    #title="Wyznaczanie macierzy fundamentalnej"
#    title="Determining fundemental matrix component";
#    
#    @property
#    def header_text(self):
#        #"Z równań ruchu wyznaczono macierz mas i sztywności układu:"
#        # google tlumacz
#        return "The matrix of masses and stiffnesses of the system was determined from the equations of motion:"
#
#        
#    @property
#    def body_text(self):
#        #"Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu #${latex(Delta)}$, przedstawiają się następująco:"
#        # google tlumacz
#        return "The fundamental matrix, on the basis of which the characteristic equation of the considered #system ${latex (Delta)}$ was determined, is as follows:"
#    
#    @property
#    def footer_text(self):
#        #" Macierz fundamentalna pozwala określić rozwiązanie ustalone. Natomiast bazując na równaniu #charakterystycznym określa się częstości własne układu."
#        # google tlumacz
#        return "The fundamental matrix allows you to define a fixed solution. On the other hand, based on the #characteristic equation, the eigenfrequencies of the system are determined."
#    
#    
#    def append_elements(self):
#        
#        system = self._system
#        ReportText.set_directory('./SDAresults')
#
#        latex_store=AutoBreak.latex_backend
#        AutoBreak.latex_backend = latex_store
#        
#        t=system.ivar
#        
#
#        dyn_sys=system
#        dyn_sys_lin=dyn_sys.linearized()
#
#
#        display(ReportText(self.header_text))
#        
#        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker='a' )  ))
#
#        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker='a')  ))
#
#        Delta = Symbol('\Delta')
#
#        display(ReportText(self.body_text))
#
#        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker='a'  )  ))
 #       display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker='a',backend=latex  )  ))
#
#        display(ReportText(self.footer_text))
#
 #       AutoBreak.latex_backend = latex_store

class FundamentalMatrixDynPyCodeComponent(ReportComponent):
    
    #title="Wyznaczanie macierzy fundamentalnej"
    title="Determining fundemental matrix component";
    
    @property
    def header_text(self):
        #"Z równań ruchu wyznaczono macierz mas i sztywności układu:"
        # google tlumacz
        return "The matrix of masses and stiffnesses of the system was determined from the equations of motion:"

        
    @property
    def body_text(self):
        #"Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu ${latex(Delta)}$, przedstawiają się następująco:"
        # google tlumacz
        return "The fundamental matrix, on the basis of which the characteristic equation of the considered system ${latex (Delta)}$ was determined, is as follows:"
    
    @property
    def footer_text(self):
        #" Macierz fundamentalna pozwala określić rozwiązanie ustalone. Natomiast bazując na równaniu charakterystycznym określa się częstości własne układu."
        # google tlumacz
        return "The fundamental matrix allows you to define a fixed solution. On the other hand, based on the characteristic equation, the eigenfrequencies of the system are determined."
    
    
    def append_elements(self):
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(self.header_text))
        
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker='a' )  ))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker='a')  ))

        Delta = Symbol('\Delta')

        display(ReportText(self.body_text))

        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker='a'  )  ))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker='a',backend=latex  )  ))
        
        display(Markdown('''
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker='a' )  ))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker='a')  ))

        Delta = Symbol('\Delta')
        
        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker='a'  )  ))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker='a',backend=latex  )  ))
        
        
        '''))

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store    

class FundamentalMatrixSymPyCodeComponent(ReportComponent):
    
    title="Determining fundemental matrix component";
    
    @property
    def header_text(self):
        return "The mass and stiffness matrices of the system were determined from the equations of motion as follows:"

    @property
    def body_text(self):

        return "The fundamental matrix, based on which the characteristic equation of the considered system ${latex(Delta)}$ was determined, is represented as follows:"
    
    @property
    def footer_text(self):
        return "The fundamental matrix allows for determining the steady-state solution. Meanwhile, based on the characteristic equation, one often determines the eigenfrequencies of the system."
    
    def append_elements(self):
        
        system = self._system

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()

        display(ReportText(self.header_text))
        
        
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False))))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False))))

        Delta = Symbol('\Delta')

        display(ReportText(self.body_text))

        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False))))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False))))
        
        
        code_list=python(system.fundamental_matrix()).split('\n')

        var_name = 'Fund_mat'

        code = '\n\n\t'.join(code_list[:-1]) + '\n\n\t' + var_name +code_list[-1][1:]

        display(Markdown(
f'''

    from sympy import *

    {code}


    display(Eq(Symbol('Delta'),Fund_mat))

'''))
        
        display(ReportText(self.footer_text))
    
# Mateusz
class GeneralSolutionComponent(ReportComponent):
    #"Rozwiązanie ogólne"
    title="General solution"
    
    @property
    def header_text(self):
        #'Rozwiązanie ogólne przedstawia wyrażenie:'
        return "General solution is presented by expression:"
    @property
    def footer_text(self):
        #'Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.'
        return "General solution describes motion of the analised system - presents displacement i function of time - and is given by considerations about free vibrations of the system"
    
    def append_elements(self):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        #from dynpy import LagrangesDynamicSystem, HarmonicOscillator

        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        
        
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(self.header_text))
#         display((SympyFormula(  Eq(Symbol('X'),HarmonicOscillator(dyn_sys_lin).general_solution()[0].n(3),
#                                             evaluate=False) , marker='a',backend=latex  )  ))
#         display(dyn_sys_lin._eoms)
        display((SympyFormula(  Eq(Symbol('X_g'), dyn_sys_lin._ode_system.general_solution.n(3), evaluate=False) , marker='a',backend=latex  )  ))

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store
        
class GeneralSolutionDynpyCodeComponent(ReportComponent):
    #"Rozwiązanie ogólne"
    title="Dynpy code for general solution"
    
    @property
    def header_text(self):
        #'Rozwiązanie ogólne przedstawia wyrażenie:'
        return "Dynpy code for general solution is presented by expression:"
    @property
    def footer_text(self):
        #'Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.'
        return "General solution describes motion of the analised system - presents displacement i function of time - and is given by considerations about free vibrations of the system. This code presents a solution in dynpy code."
    
    def append_elements(self):
        system=self._system
        t=system.ivar
        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()
        display(ReportText(self.header_text))
        coord=dyn_sys.q
        
        ode=system._ode_system

        display(Markdown(
f'''

    from dynpy.models.mechanics import *
    from sympy import *

    dyn_sys = {system.__class__.__name__}()
    
    
    dyn_sys._ode_system.general_solution
    


'''))
        
        display(ReportText(self.footer_text))
        
class GeneralSolutionSympyCodeComponent(ReportComponent):
    #"Rozwiązanie ogólne"
    title="Sympy code for General Solution"
    
    @property
    def header_text(self):
        #'Rozwiązanie ogólne przedstawia wyrażenie:'
        return "Sympy code for General Solution is presented by expression:"
    @property
    def footer_text(self):
        #'Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.'
        return "General solution describes motion of the analised system - presents displacement i function of time - and is given by considerations about free vibrations of the system.  This code presents a solution in sympy code."
    
    def append_elements(self):
        system=self._system
        t=system.ivar
        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()
        ode=system._ode_system.lhs-system._ode_system.rhs
        cord=system.q

        
      
        display(ReportText(self.header_text))

        code_list_ode=python(ode).split('\n')

        var_name = 'ode'
        
        code_ode = '\n\n\t'.join(code_list_ode[:-1]) + '\n\n\t' + var_name +code_list_ode[-1][1:]

        display(ReportText(self.footer_text))

        code_list_dvars=python(system.q).split('\n')

        var_name = 'dvars'
        
        code_dvars = '\n\n\t'.join(code_list_dvars[:-1]) + '\n\n\t' + var_name +code_list_dvars[-1][1:]
        
        
        display(Markdown(
f'''
    from sympy import *
    init_printing()
    

    {code_ode}
    
    {code_dvars}
    
    dsolve(ode,dvars)


'''))

        
# Grześ
class FrequencyResponseFunctionComponent(ReportComponent):
    #"Charakterystyka Amplitudowo-Częstotliwościowa"
    title="amplitude-frequency characteristic"
    
    @property
    def header_text(self):
        #"funkcja odpowiedzi częstotliwościowej:"
        return "Frequency response is given by the formula:"
    @property
    def footer_text(self):
        #"jest to suma kwadratów amplitud pod pierwiastkiem:"
        return "it is the sum of the squared amplitudes under the root:"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('frf'),
                     dyn_sys.frequency_response_function() ), marker=None))

        display(ReportText(self.footer_text))
        
FRFComponent = FrequencyResponseFunctionComponent


#Pati
class FrequencyResponseFunctionComponentToSecond(ReportComponent):
    
    title="sens mocy układu"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

    @property
    def header_text(self):
        #"sens mocy układu:"
        return "system power sense:"

        display(SympyFormula( Eq(Symbol('frf^2'),
                     dyn_sys.frequency_response_function().doit()**2 ), marker=None))
    @property
    def footer_text(self):
        #"Charakterystyka Amplitudowo-Częstotliwościowa podniesiona do kwadratu"
        return "The amplitude-frequency response squared"
#Pati        
class SpringForce(ReportComponent):
    
    title="Spring force"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


    @property
    def header_text(self): 
        #"Siła od sprężyny wyrażona jest wzorem:""
        return "The spring force is given by the formula:"

    @property
    def middle_text(self): 
        #"zamiast x używam steady solution"
        return "steady solution was used instead of z"

    @property
    def footer_text(self):
           #"Siła od sprężyny, zwana również siłą naciągu, pojawia sie przy ściskaniu lub rozciaganiu. Siła, która działa jest przeciwnie skierowana do ruch i chce przywrócić do pierwotnego jej położenia. Zależy od sztywności sprężyny k oraz od tego o ile została rozciagnieta bądź skrócona x."
        return "Spring force, also known as pull force, occurs when compressed or stretched. The force that acts is opposite to the movement and wants to return it to its original position. It depends on the spring stiffness k and how much it is stretched or shortened z."

    def append_elements(self):


        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F'),-1*system._left_mount.stiffness*system.z,evaluate=False)))

        display(ReportText(self.middle_text))
        
        display(SympyFormula( Eq(Symbol('F'),-1*system._left_mount.stiffness*system._ode_system.steady_solution[0],evaluate=False)))
        
        display(ReportText(self.footer_text))
        
FRFComponent = FrequencyResponseFunctionComponent
#Pati
class DamperForce(ReportComponent):
    
    title="Siła tłumienia"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


    @property
    def header_text(self): 
        #"Siła tłumienia wyrażona jest wzorem:"
        return "The damping force is given by the formula:"

        display(SympyFormula( Eq(Symbol('F'),-1*system.c*system.v,evaluate=False)))

    @property
    def header_text(self): 
             #"zastępuje prędkość jako pochodna steady stolution po czasie:"
        return "replaces velocity as a derivative of the steady solution over time:"

        display(SympyFormula( Eq(Symbol('F'),-1*system.c*system.steady_solution().diff(system.ivar),evaluate=False)))

    @property
    def footer_text(self):
         #"Siła tłumienia zmniejsza amplitude drgań, ma zwrot przeciwny do prędkości. Zależy od współczynnika tłumienia b oraz od prędkości v."
        return "Siła tłumienia zmniejsza amplitude drgań, ma zwrot przeciwny do prędkości. Zależy od współczynnika tłumienia b oraz od prędkości v."
#Pati
class LogarithmicDecrement(ReportComponent):
    
    title="logarytmiczny dekrement tłumienia liczony z amplitud"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


    @property
    def header_text(self): 
        #"Logarytmiczny dekrement tłumienia wyrażona jest wzorem przy pomocy amplitud:"
        return "The logarithmic decrement of damping is given by the formula in terms of amplitudes."

        display(SympyFormula( Eq(Symbol('delta'),log(system.A_n1/system.A_n2))))

    @property
    def footer_text(self): 
            #"Zaprezentowana zależność opisuje bezwymaiarową wielkość charakteryzującą intensywność tłumienia drgań swobodnych w przypadku podkrytycznym. Jest to wielkość stała dla rozpaywanego                                  układu drgającego i nie zależy od warrunków początkowych. Zależy natomiast od maksymalnych wychyleń w chwilach różniących się o okres drań tłumionych."
            return "The presented dependence describes a dimensionless quantity characterizing the intensity of damping of free vibrations in the subcritical case. It is a constant value for the considered oscillating system and does not depend on the initial conditions. However, it depends on the maximum deflections at moments differing by the period of damped bastards."
    @property
    def header_text(self): 
           #"Logarytmiczny dekrement tłumienia u wyrażona jest wzorem przy pomocy okresu drgań tłumionych:"
        return "The logarithmic damping decrement u is given by the formula using the period of damped vibrations:"

        display(SympyFormula( Eq(Symbol('delta'),2*pi*system.damping_coefficient()*(system.natural_frequencies())**-1,evaluate=False)))
        
    @property
    def footer_text(self): 
           #"Zaprezentowana zależność opisuje bezwymaiarową wielkość charakteryzującą intensywność tłumienia drgań swobodnych w przypadku podkrytycznym. Jest to wielkość stała dla rozpaywanego układu drgającego i nie zależy od warrunków początkowych.Zależy natomiast od okresu drań tłumionych i współczynnika h."
        return "The presented dependence describes a dimensionless quantity characterizing the intensity of damping of free vibrations in the subcritical case. It is a constant value for the considered vibrating system and does not depend on the initial conditions, but it depends on the period of damped vibrations and the coefficient h."



        
# Marcel & Monika
class SteadySolutionComponent(ReportComponent):
    
        #"Rozwiązanie szczególne"
    title="Steady solution"
    _phi=False
    
    @property
    def header_text(self):
        #"Rozwiązanie szczególne przedstawia wyrażenie:"
        # google tlumacz
        return "The steady solution is given by the formula:"

        
    @property
    def footer_text(self):
        #" Rozwiązanie szczególne związane jest obecnością wielkości wymuszających ruch (drgania) analizowanego układu."
        # google tlumacz
        return "The specific solution is related to the presence of quantities that force motion (vibrations) of the analyzed system."
    
    def append_elements(self,phi=_phi):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()



        display(ReportText(self.header_text ))

#         display((SympyFormula(  Eq(Symbol('X_s'),
#                             HarmonicOscillator(dyn_sys_lin).steady_solution()[0].n(3),
#                             evaluate=False) , marker='b',backend=latex  )  ))
        
        display((SympyFormula(  Eq(Symbol('X_s'), dyn_sys_lin._ode_system.steady_solution.n(3), evaluate=False) , marker='b',backend=latex  )  ))

        AutoBreak.latex_backend = latex_store

        display(ReportText(self.footer_text ))
        
#### Amadi
class MaxStaticForce(ReportComponent):
    
    #Maksymalna siła statyczna
    title="Maximum static force"

    @property
    def header_text(self):
        #"Wartość maksymalna siły statycznej działającej na pojedynczy element mocujący:"
        return "Maximum value of the static force acting on a single element:"
    
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F_s'),
                     dyn_sys.max_static_force().doit()/2 ), marker=None))
        
#### Amadi
class MaxDynamicForce(ReportComponent):
    
    #Maksymalna siła statyczna
    title="Maximum dynamic force"

    @property
    def header_text(self):
        #"Wartość maksymalna siły dynamicznej działającej na pojedynczy element mocujący:"
        return "Maximum value of the dynamic force acting on a single element:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F_d'),
                     dyn_sys.max_dynamic_force().doit() ), marker=None))
        
        
#### Amadi
class DynamicPinDiameter(ReportComponent):
    
    title="Minimum diameter of pin due to dynamic force"

    @property
    def header_text(self):
        return "Minimum diameter of the pin due to dynamic force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('d'),
                     dyn_sys.dynamic_force_pin_diameter().doit() ), marker=None))
        
#### Amadi
class StaticPinDiameter(ReportComponent):
    
    title="Minimum diameter of pin due to static force"

    @property
    def header_text(self):
        return "Minimum diameter of the pin due to static force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('d'),
                     dyn_sys.static_force_pin_diameter().doit() ), marker=None))
        
#### Amadi
class TensionerForce(ReportComponent):
    
    title="Tensioner force"


    @property
    def header_text(self): 
        #"Siła od sprężyny wyrażona jest wzorem:""
        return "The tensioner force modelled as a spring is given by the formula:"

    @property
    def middle_text(self): 
        #"zamiast x używam steady solution"
        return "steady solution was used instead of z"

    @property
    def footer_text(self):
           #"Siła od sprężyny, zwana również siłą naciągu, pojawia sie przy ściskaniu lub rozciaganiu. Siła, która działa jest przeciwnie skierowana do ruch i chce przywrócić do pierwotnego jej położenia. Zależy od sztywności sprężyny k oraz od tego o ile została rozciagnieta bądź skrócona x."
        return "Spring force, also known as pull force, occurs when compressed or stretched. The force that acts is opposite to the movement and wants to return it to its original position. It depends on the spring stiffness k and how much it is stretched or shortened z."

    def append_elements(self):


        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F'),-1*system._tensioner.stiffness*(system.z0-system.z),evaluate=False)))

        display(ReportText(self.middle_text))
        
        display(SympyFormula( Eq(Symbol('F'),-1*system._tensioner.stiffness*system.steady_solution()[0],evaluate=False)))
        
        display(ReportText(self.footer_text))
        
#### Amadi
class TensionerDamperForce(ReportComponent):
    
    title="Tensioner damping force"


    @property
    def header_text(self): 
        #"Siła tłumienia wyrażona jest wzorem:"
        return "The damping force is given by the formula:"


    @property
    def middle_text(self): 
             #"zastępuje prędkość jako pochodna steady stolution po czasie:"
        return "replaces velocity as a derivative of the steady solution over time:"


    @property
    def footer_text(self):
         #"Siła tłumienia zmniejsza amplitude drgań, ma zwrot przeciwny do prędkości. Zależy od współczynnika tłumienia b oraz od prędkości v."
        return "The damping force, which reduces the amplitude of vibrations, has a direction opposite to the velocity. It depends on the damping coefficient c and the speed V."

    def append_elements(self):


        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F'),-1*system._tensioner_damping.c*Symbol('V'),evaluate=False)))

        display(ReportText(self.middle_text))
        
        display(SympyFormula( Eq(Symbol('F'),-1*system._tensioner_damping.c*system.steady_solution()[0].diff(system.ivar),evaluate=False)))
        
        display(ReportText(self.footer_text))
        

#### Amadi
class StaticKeyLength(ReportComponent):
    
    title="Minimum length of key due to static force"

    @property
    def header_text(self):
        return "Minimum length of the key due to static force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('l'),
                     dyn_sys.static_key_length().doit() ), marker=None))
        
#### Amadi
class DynamicKeyLength(ReportComponent):
    
    title="Minimum length of key due to dynamic force"

    @property
    def header_text(self):
        return "Minimum length of the key due to dynamic force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('l'),
                     dyn_sys.dynamic_key_length().doit() ), marker=None))

class StaticAndDynamicCableForceComponent(ReportComponent):
    
    title="Pendulum's cable longitudinal force value"

    @property
    def header_text(self):
        return "Pendulum's cable static force is represented by formula:"
    
    @property
    def footer_text(self):
        return "Pendulum's cable dynamic force is represented by formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F_s'),
                     dyn_sys.max_static_cable_force().doit() ), marker=None))
        
        display(ReportText(self.footer_text))

        display(SympyFormula( Eq(Symbol('T'),
                     dyn_sys.max_dynamic_cable_force().doit() ), marker=None))

        
class PendulumLongitudinalForce(StaticAndDynamicCableForceComponent):
    
    @property
    def middle_text(self):
        return "Pendulum's maximum angular velocity has following form:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.middle_text))

        display(SympyFormula( Eq(Eq(dynamicsymbols('varphi_max').diff(dyn_sys.ivar) , Symbol('A_omega'), evaluate=False),
                     dyn_sys.linearized().frequency_response_function() * dyn_sys.Omega, evaluate=False ), marker=None))
        

class MDoFGeneralSolutionComponent(ReportComponent):
    #"Rozwiązanie ogólne"
    title="General solution"
    
    @property
    def header_text(self):
        #'Rozwiązanie ogólne przedstawia wyrażenie:'
        return "General solution is presented by expression:"
    @property
    def footer_text(self):
        #'Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.'
        return "General solution describes motion of the analised system - presents displacement i function of time - and is given by considerations about free vibrations of the system"
    
    def append_elements(self):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        #from dynpy import LagrangesDynamicSystem, HarmonicOscillator

        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        
        
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(self.header_text))

        for gen_sol in HarmonicOscillator(dyn_sys_lin.linearized()).general_solution().n(3):
                 display((SympyFormula(Eq(Symbol('X'), gen_sol, evaluate=False) , marker='a',backend=latex)))

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store

class MDoFSteadySolutionComponent(ReportComponent):
    
        #"Rozwiązanie szczególne"
    title="Steady solution"
    _phi=False
    
    @property
    def header_text(self):
        #"Rozwiązanie szczególne przedstawia wyrażenie:"
        # google tlumacz
        return "The steady solution is given by the formula:"

        
    @property
    def footer_text(self):
        #" Rozwiązanie szczególne związane jest obecnością wielkości wymuszających ruch (drgania) analizowanego układu."
        # google tlumacz
        return "The specific solution is related to the presence of quantities that force motion (vibrations) of the analyzed system."
    
    def append_elements(self,phi=_phi):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()



        display(ReportText(self.header_text ))

        steady_sol=dyn_sys_lin._fodes_system.steady_solution

        for q in dyn_sys_lin.q:

            display((SympyFormula(  Eq(Symbol('X_s'), steady_sol.as_dict()[q], evaluate=False) , marker='b',backend=latex)))

        AutoBreak.latex_backend = latex_store

        display(ReportText(self.footer_text ))
        
        
        
        
# Grzes
class RightSpringForceComponent(ReportComponent):

    title="Force in top right spring of the system"

    @property
    def header_text(self):
        #"Energia kinetyczna układu wyrażona jest wzorem:"
        return "Force in the right spring is determined by the formula:"


    @property
    def footer_text(self):
        #"Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
        return "Determined formula specifies force in the spring"

    def append_elements(self):

        system = self._system
        kr=system.k_r
        x_2=system.x_2
        Fr=-kr*x_2
        
        display(ReportText( self.header_text))
        
        display(SympyFormula( Eq(Symbol('F_r'),
                     Fr) , marker=None))

        display(SympyFormula( Eq(Symbol('F_r'),
                     system.right_spring_force()) , marker=None))

        
        display(ReportText( self.footer_text))
        
class CentralSpringForceComponent(ReportComponent):

    title="Force in top centre spring of the system"

    @property
    def header_text(self):
        #"Energia kinetyczna układu wyrażona jest wzorem:"
        return "Force in the top centre spring is determined by the formula:"


    @property
    def footer_text(self):
        #"Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
        return "Determined formula specifies force in the spring"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        kc=system.k_c
        x_1=system.x_1
        x_2=system.x_2
        Fc=-kc*(x_2-x_1)
        
        display(ReportText( self.header_text))
        
        display(SympyFormula( Eq(Symbol('F_c'),
                     Fc) , marker=None))

        display(SympyFormula( Eq(Symbol('F_c'),
                     system.centre_spring_force()) , marker=None))


        display(ReportText( self.footer_text))
        
class FreeVibrationFrequencyComponent(ReportComponent):
    title = "Natural Frequency"


    @property
    def header_text(self):
        return " natural frequency is given by the formula:"

    @property
    def footer_text(self):
        return "The specific solution is related to the presence of quantities that force motion (vibrations) of the analyzed system."

    def append_elements(self):
        system = self._system
        t = system.ivar
        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        display(ReportText(self.header_text))

        # Calculate and display the steady solution
        char_poly = dyn_sys_lin.fundamental_matrix().det()
        display(SympyFormula(char_poly))

        omg = Symbol('omega',positive=True)

        # You may want to add more code here to display additional information or perform actions.
#         display(char_poly.coeff(omg**2))
#         display(char_poly.coeff(omg**1))
#         display(char_poly.subs(omg,0))


        display(ReportText(f"Natural Frequency in rad per sec (omage0): "))
#         display(SympyFormula(solve(char_poly,Symbol('omega',positive=True))[0]))
        eigs = dyn_sys_lin.natural_frequencies()

        display(SympyFormula( eigs ))
        
        
        display(SympyFormula( 2*pi*eigs ))
        
        
        display(SympyFormula( (2*pi*eigs)**-1 ))
        
        display(ReportText(self.footer_text))
        
class StaticForceComponent(ReportComponent):

    title="Force of the considered element"

    @property
    def header_text(self):
        return "Force of the considered element is following:"

    def get_force(self):
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        return dyn_sys.static_force().doit()

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

        force=self.get_force()

        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F'),
                         force), marker=None))

class DynamicForceComponent(StaticForceComponent):
    def get_force(self):
        system=self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        return dyn_sys.dynamic_force().doit()
    
class DynamicTensionForceComponent(StaticForceComponent):
    def get_force(self):
        system=self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        return dyn_sys.force_in_cable().doit()

class StaticCableDiameter(ReportComponent):
    
    title="Minimum diameter of the cable due to static force"

    @property
    def header_text(self):
        return "Minimum diameter of the cable due to static force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('d'),
                     dyn_sys.static_cable_diameter().doit() ), marker=None))

class DynamicCableDiameter(ReportComponent):
    
    title="Minimum diameter of the cable due to dynamic force"

    @property
    def header_text(self):
        return "Minimum diameter of the cable due to dynamic force formula:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('d'),
                     dyn_sys.dynamic_cable_diameter().doit() ), marker=None))

class DampedVibrationFrequency(ReportComponent):
    
    title="damped vibration frequency"
    @property
    def header_text(self):
        return "damped vibration frequency(omega'h):"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('omega_h'),
                     ((list(dyn_sys()._ode_system.general_solution.atoms(sin))[0].args[0])/dyn_sys.ivar)**2 ), marker=None))

class DampedEngineDynamicShearRatioComponent(ReportComponent):
    
    title="Ratio of undamped and damped dynamic shear force"

    @property
    def header_text(self):
        return "Ratio of undamped and damped dynamic shear force is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        from dynpy.models.mechanics.engine import EngineConstantVelocityVerticalSpringGravity
        
        undamped=EngineConstantVelocityVerticalSpringGravity()
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('phi'),
                     undamped.max_static_force()/dyn_sys.max_static_force()), marker=None))

        
class DampingMatrixComponent(ReportComponent):
    title="Damping matrix of the system"
    
    @property
    def header_text(self):
         return "Damping matrix of the system is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('C'),
                     dyn_sys.damping_matrix(), evaluate=False, marker=None)))

class TMDForceComponent(ReportComponent):
    title="Dynamic force in the system with TMD and without"
    
    @property
    def header_text(self):
         return "Dynamic force in the system with TMD and without is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F_tmd'),
                     dyn_sys.dynamic_force_TMD(), marker=None)))
        display(SympyFormula(Eq(Symbol('F'), dyn_sys.dynamic_force(), marker=None)))
        
class SmallParameterComponent(ReportComponent):
    title="Small parameter"
    
    @property
    def header_text(self):
         return "Small parameter of the system is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('epsilon'),
                     dyn_sys.small_parameter(), marker=None)))
        
class NonlinearSteadySolutionComponent(ReportComponent):
    title="Nonlinear Steady Solution"
    
    @property
    def header_text(self):
         return "Nonlinear steady solution of the system is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(dyn_sys.qs[0], solve(dyn_sys.frequency_response_function(),a)[2]*sin(dyn_sys.Omega*dyn_sys.ivar), marker=None)))
        
class AmplitudeAndFrequencyRelationComponent(ReportComponent):
    title="Nonlinear Steady Solution"
    
    @property
    def header_text(self):
         return "Nonlinear steady solution of the system is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('omega')**2, dyn_sys.freq_amp_rel(), marker=None)))
        
class DynamicTorqueComponent(ReportComponent):
    title="Dynamic Torque"
    
    @property
    def header_text(self):
         return "Dynamic torque of the system is following:"
    
    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys
        
        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('T_d'), dyn_sys.dynamic_torque(), marker=None)))
    
class GivenDataComponent(ReportComponent):

    title="Table with parameter values for calculations"

    @property
    def entry_text(self):

        return "The values of individual parameters adopted for the calculations are presented in the table {tabelkamarker}"
    
    @property
    def table_caption(self):
        return "Default parameter values"


    def append_elements(self):

        dyn_sys = self.reported_object
        
        given_data=dyn_sys._given_data

        if given_data == {}:
            params_dict = dyn_sys.get_reference_data()

        else: params_dict = given_data

        preview_mapper = lambda x: f'${latex(x)}$' if isinstance(x,(Symbol,Eq,Function,Mul)) else str(x)

        mapped_dict = {k: f'${latex(v)}$' if isinstance(v,(Symbol,Eq,Function,Mul)) else str(v) for k, v in params_dict.items()}

        params_dict_for_ldf={'Parameter':list(params_dict.keys()),'Value':list(mapped_dict.values())}
        params_ldf = pd.DataFrame(data=params_dict_for_ldf)
        params_ldf['Parameter'].apply(preview_mapper)
        params_ldf=params_ldf.set_index('Parameter')

        table_params=LatexDataFrame.formatted(params_ldf)

        display(ReportText(self.entry_text.format(tabelkamarker = AutoMarker(table_params))))

        display(table_params.reported(index=True, caption=self.table_caption))

# class SchemeDynPyCodeComponent(ExemplaryPictureComponent):
#     title="Scheme of the system"

#     @property
#     def header_text(self):
#         #"Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."
#         return "Scheme of the real object is presented in the figure. It was obtained by the analysis of real mechanism."


#     @property
#     def footer_text(self):

#         system = self._system

#         #"Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}."
#         return f"Analysis of the scheme allows to claim its number of degrees of freedom which is {len(system.q)}"

#     def append_elements(self):

#         system = self._system

#         display(ReportText(  self.header_text ))
#         display(Picture(system._scheme(),width='8cm'))

#         display(Markdown(

# '''The code for calling scheme of the system is following:

#     display(Picture(system._scheme(),width='8cm'))

# '''
#             ))

#         display(ReportText( self.footer_text ))
        
# class SchemeSymPyCodeComponent(ExemplaryPictureComponent):
#     title="Scheme of the system"

#     @property
#     def header_text(self):
#         #"Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."
#         return "Scheme of the real object is presented in the figure. It was obtained by the analysis of real mechanism."

        
#     @property
#     def footer_text(self):
        
#         system = self._system
        
#         #"Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}."
#         return f"Analysis of the scheme allows to claim its number of degrees of freedom which is {len(system.q)}"
    
    
#     def append_elements(self):
        
#         system = self._system
#         path=system._default_folder_path + system.scheme_name

#         display(ReportText(  self.header_text ))
#         display(Picture(path,width='8cm'))
        
#         display(Markdown(
# '''The code for calling scheme of the system is following:

#     path=self._default_folder_path + self.scheme_name
    
#     display(Picture(path ,width='8cm'))

# '''
#             ))
#         display(ReportText( self.footer_text ))
