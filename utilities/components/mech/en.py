from ..mechanics import *



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

            
            
            if len(system.q)==1:
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

        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example(),width='8cm')

        display(ReportText(self.footer_text))

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
          
        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._scheme(),width='8cm')

        display(ReportText( self.footer_text ))



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
        return "Kinetic energy of the system has a following form!:"

        
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
    
    
    
class PotentialEnergyComponent(ReportComponent):#Jaś fasola
    
    title="Potential energy"
    @property
    def header_text(self):
        #"Energia potencjalna układu wyrażona jest wzorem:"
        return "Potential energy of the system has a following form!:"

        
    @property
    def footer_text(self):
        #"Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."
        return "The presented relationship describes the interaction of potential force fields in which the object is located."

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(  self.header_text   ))

        display(SympyFormula( Eq(Symbol('V'),
                     dyn_sys_lin._potential_energy ), marker=None))

        display(ReportText(  self.footer_text   ))

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


# Marcel & Monika
class FundamentalMatrixComponent(ReportComponent):
    
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

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store

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
        display((SympyFormula(  Eq(Symbol('X'),HarmonicOscillator(dyn_sys_lin.linearized(
                                            )).general_solution().n(3),
                                            evaluate=False) , marker='a',backend=latex  )  ))

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store

        
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



class FrequencyResponseFunctionComponentToSecond(ReportComponent):
    
    title="sens mocy układu"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           sens mocy układu:
                           '''))

        display(SympyFormula( Eq(Symbol('frf^2'),
                     dyn_sys.frequency_response_function().doit()**2 ), marker=None))

        display(ReportText(f'''
                           Charakterystyka Amplitudowo-Częstotliwościowa podniesiona do kwadratu
                           '''))  
        


        
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

        display((SympyFormula(  Eq(Symbol('X_s'),
                            HarmonicOscillator(dyn_sys_lin.linearized(
                            )).steady_solution().n(3),
                            evaluate=False) , marker='b',backend=latex  )  ))

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
                     dyn_sys.max_static_force().doit() ), marker=None))
        
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