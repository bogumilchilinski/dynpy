from ..mechanics import *
from . import en
import dynpy.solvers.linear as comps


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
            
#             self.append(NoEscape(f'{system._label} \n \n'))

            if len(system._given_data) == 0:
                self.append(NoEscape(f'{system.readable_name} z {len(system.q)} St.S. \n \n'))
            else:
                data_eqns= [f'${latex(key)} = {latex(value)}$' for key, value in system._given_data.items()]

                self.append(NoEscape(f'''{system.readable_name} z {len(system.q)} St.S. dla {f', '.join(data_eqns)} \n \n'''))

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
    
    title="Przykład rzeczywistego obiektu"
    packages=[Package('float')] 
    @property
    def header_text(self):
        
        return "Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."
    @property
    def footer_text(self):


        return "Model dynamiczny układu określa się na podstawie analizy rozważanego przypadku. Należy pamiętać, że stopień odwzorowania (poziom abstrakcji) modelu zależy od tego do czego planuje się go używać."
    
    def append_elements(self):
        
        system = self.reported_object

        display(ReportText(self.header_text))

        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example(),width='8cm')

        display(ReportText(self.footer_text))

# Boogi        
class SchemeComponent(ExemplaryPictureComponent):
    title="Schemat systemu"

    @property
    def header_text(self):
        #
        return "Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."

        
    @property
    def footer_text(self):
        
        system = self._system

        return f"Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}."
    
    
    def append_elements(self):
        
        system = self.reported_object

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
        
        system = self.reported_object

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
class KineticEnergyComponent(en.KineticEnergyComponent):
    
    title="Energia kinetyczna"

    @property
    def header_text(self):
        #"Energia kinetyczna układu wyrażona jest wzorem:"
        return "Energia kinetyczna układu wyrażona jest wzorem:"

        
    @property
    def footer_text(self):
        #"Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
        return "Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
    

    
    
    
class PotentialEnergyComponent(ReportComponent):#Jaś fasola
    
    title="Energia potencjalna"
    @property
    def header_text(self):
        
        return "Energia potencjalna układu wyrażona jest wzorem:"

        
    @property
    def footer_text(self):
        
        return "Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."

    def append_elements(self):

        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(  self.header_text   ))

        display(SympyFormula( Eq(Symbol('V'),
                     dyn_sys_lin._potential_energy ), marker=None))

        display(ReportText(  self.footer_text   ))

# Marcel spolszczone
class DissipationComponent(ReportComponent):
    
    
    title="Dyssypacyjna funkcja Rayleigh'a"
    
    
    @property
    def header_text(self):
        
        return "Energia rozpraszana tłumieniem wyrażona jest wzorem:"

        
    @property
    def footer_text(self):
        
        return "Podana zależność stanowi potencjał dysynpacyjny Rayleigh'a, który poddany różniczkowaniu względem wektora prędkości uogólnionych pozwala na określenie sił wiskotycznego tłumienia."
    
    def append_elements(self):
        
        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText( self.header_text  ))
        

        display(SympyFormula( Eq(Symbol('D'),
             dyn_sys_lin._dissipative_potential) , marker=None ))

        display(ReportText( self.footer_text  ))

# Boogi          
class LagrangianComponent(en.LagrangianComponent):
    
    #title="Lagrangian (funkcja Lagrange'a)  układu"
    title="Lagrangian układu (Funkcja Lagrange'a)"
        
    @property
    def header_text(self):
        #"Lagrangian systemu wyrażony jest wzorem ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):"

        
        return f"Lagrangian układu dany jest następującym wyrażeniem (marker):"
        
    @property
    def body_text_1(self):
        #"Równania Eulera-Lagrange'a dla rozważanego przypadku są nastęujące: "   
        return "Równania Eulera Lagrange'a dla rozważanego przypadku są następujące: "
        
    @property
    def body_text_2(self):
        #"Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: "
        return "Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące:"
    
    @property
    def footer_text(self):
        #"Wyniki przedstawionych operacji wykorzystuje się wyznaczenia równań ruchu układu."
        return "Wyniki przedstawionych operacji wykorzystuje się wyznaczenia równań ruchu układu."
        
    def append_elements(self):
        
        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        
        print(system._scheme())
        
        mrk_lagrangian_nonlin=Marker('lagrangianNL',prefix='eq')

        #display(ReportText(f'''The following model is considered. The system's Lagrangian is described by the formula ({Ref(mrk_lagrangian_nonlin).dumps()}):
        #                    '''))
        eq_lagr = Eq(Symbol('L'),dyn_sys.L.expand()[0])
        
        display(ReportText(self.header_text.replace('marker', f'{AutoMarker(eq_lagr)}') ))

        display((SympyFormula(  eq_lagr)  ))
        
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
    title="Równanie ruchu"
    
    @property
    def entry_text_one(self):
        system = self.reported_object
        dyn_sys=system
        return f"Wykorzystując obliczone pochodne, wyznacza się równanie ruchu na podstawie odpowiedniego wzoru. Równanie ruchu układu przedstawia zależność: ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})"
    
    @property
    def entry_text_two(self):
        system = self.reported_object
        dyn_sys=system
        
        return f"Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru. Równania ruchu układu przedstawiają zależności: ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-     ({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"

    @property
    def footer_text(self):
        #f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return f'''  Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
    
    def append_elements(self):
        
        system = self.reported_object
        dyn_sys=system#.subs(system._given_data)

        if len(system.q)==1:
            display(ReportText(self.entry_text_one))   
        else:
            display(ReportText(self.entry_text_two))   

        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=None))

        display(ReportText(self.footer_text))   

class LinearizedGoverningEquationComponent(ReportComponent):
    #Równania ruchu
    title="Równanie ruchu"

    @property
    def entry_text_one(self):
        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()
        
        return f"Wykorzystując obliczone pochodne, wyznacza się równanie ruchu na podstawie odpowiedniego wzoru. Zlinearyzowane równanie ruchu układu przedstawia zależność: ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})"
    
    @property
    def entry_text_two(self):
        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()
        
        return f"Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru. Zlinearyzowane równania ruchu układu przedstawiają zależności: ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys_lin._eoms[-1].simplify().expand(),0))})"
    
    @property
    def footer_text(self):
        
        return f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
    
    def append_elements(self):
        
        system = self.reported_object
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

    title="Linearyzacja równań ruchu"

    @property
    def entry_text(self):
        return "Linearyzaja równań polega na znalezieniu ich rozwinięcia w szereg Taylora względem współrzędnych, prędkości i przyspieszeń uogólnionych w otoczeniu punktu równowagi. Celem uproszczenia wprowadzono następujące oznaczenia:"

    @property
    def equilibrium_point_text(self):
        return "Punkty równowagi rozważanego układu są następujące:"
    @property
    def eom_text(self):
        dyn_sys= self._system
#         return '''Równanie ruchu dla współrzędnej ${latex(dyn_sys.q[no])}$ można przestawić jako:'''
        return '''Równanie ruchu dla współrzędnej ${coord}$ można przestawić jako:'''
    @property
    def lagrange_text(self):
        dyn_sys= self._system
        return "Formalnie należy obliczyć pochodne cząstkowe wielkości uogólnionych ze składników równań Lagrange'a"

    @property
    def derivative_text(self):
        return "Poszczególne pochodne mają następującą postać:"

    @property
    def linearized_eq_text(self):
        return "Po podstawieniu obliczonych pochodnych, otrzumuje się następujące zlinearyzowane równanie:"

    def append_elements(self):

        system = self.reported_object
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

################################ NOT SURE ######################################################
            display(ReportText( self.eom_text.format(coord = latex(dyn_sys.q[no])) ))

            #display((SympyFormula(  Eq(eq_sym,eom,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))


            display(ReportText( self.lagrange_text ))


            display((SympyFormula(  Eq(comps.MultivariableTaylorSeries(eq_sym,coords,n=1,x0=op_point)._symbolic_sum(),0) , marker=None,backend=latex  )  ))

            diff_list=comps.MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).calculation_steps(expr_symbol=eq_sym)

            display(ReportText( self.derivative_text ))

            for diff_eq in diff_list:

                display((SympyFormula(  diff_eq , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display(ReportText( self.linearized_eq_text ))
            display((SympyFormula(  Eq(comps.MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).doit().expand().simplify().expand(),0,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))



        AutoBreak.latex_backend = latex_store


# Marcel & Monika spolszczone
class FundamentalMatrixComponent(ReportComponent):
    
    title="Wyznaczanie macierzy fundamentalnej";

    
    @property
    def header_text(self):
        return "Z równań ruchu wyznaczono macierz mas i sztywności układu::"

        
    @property
    def body_text(self):

        return "Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu ${delta}$, przedstawiają się następująco::"
    
    @property
    def footer_text(self):

        return "Macierz fundamentalna pozwala określić rozwiązanie ustalone. Natomiast bazując na równaniu charakterystycznym określa się częstości własne układu."
    
    
    def append_elements(self):
        
        system = self.reported_object
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()#.subs(system._given_data)


        display(ReportText(self.header_text))
        
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker='a' )  ))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker='a')  ))

        Delta = Symbol('\Delta')

        display(ReportText(self.body_text.format(delta = latex(Delta) )))

        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker='a'  )  ))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker='a',backend=latex  )  ))

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store

# Mateusz
class GeneralSolutionComponent(ReportComponent):
    title="Rozwiązanie ogólne"
    
    @property
    def header_text(self):
        
        return 'Rozwiązanie ogólne przedstawia wyrażenie:'
    @property
    def footer_text(self):
        
        return 'Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.'
    
    def append_elements(self):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        #from dynpy import LagrangesDynamicSystem, HarmonicOscillator

        system = self.reported_object
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        
        
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()#.subs(system._given_data)


        display(ReportText(self.header_text))


        for i,qs in enumerate(system.q):

            display((SympyFormula(  Eq(Symbol(f'X_g_-{qs}'),dyn_sys_lin._ode_system.general_solution.n(3)[i], evaluate=False) , marker='a',backend=latex  )  ))


        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store

        
# Grześ
class FrequencyResponseFunctionComponent(ReportComponent):
    #"amplitude-frequency characteristic"
    title="Charakterystyka Amplitudowo-Częstotliwościowa"
    
    @property
    def header_text(self):
        #"funkcja odpowiedzi częstotliwościowej:"
        return "funkcja odpowiedzi częstotliwościowej:"
    @property
    def footer_text(self):
        #"jest to suma kwadratów amplitud pod pierwiastkiem:"
        return "jest to suma kwadratów amplitud pod pierwiastkiem:"

    def append_elements(self):

        system = self.reported_object
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

        system = self.reported_object
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
    title="Rozwiązanie szczególne"
    _phi=False
    
    @property
    def header_text(self):
        #"Rozwiązanie szczególne przedstawia wyrażenie:"
        # google tlumacz
        return "Rozwiązanie szczególne dane jest następującym wyrażeniem:"

        
    @property
    def footer_text(self):
        #" Rozwiązanie szczególne związane jest obecnością wielkości wymuszających ruch (drgania) analizowanego układu."
        # google tlumacz
        return "Rozwiązanie szczególne układu przedstawia zależność położenia od czasu odpowiednią dla drgań wymuszonych"


    def append_elements(self,phi=_phi):

        from ....dynamics import LagrangesDynamicSystem, HarmonicOscillator
        
        system = self.reported_object
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()#.subs(system._given_data)



        display(ReportText(self.header_text ))


        for i,qs in enumerate(system.q):

            display((SympyFormula(  Eq(Symbol(f'X_s_-{qs}'),dyn_sys_lin._ode_system.steady_solution.n(3)[i], evaluate=False) , marker='b',backend=latex  )  ))


        AutoBreak.latex_backend = latex_store

        display(ReportText(self.footer_text ))
        
#### Amadi
class MaxStaticForce(ReportComponent):
    
    #Maksymalna siła statyczna
    title="Maximum static force"

    @property
    def header_text(self):
        
        return "Wartość maksymalna siły statycznej działającej na pojedynczy element mocujący:"
    
    
    def append_elements(self):

        system = self.reported_object
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
       
        return "Wartość maksymalna siły dynamicznej działającej na pojedynczy element mocujący:"
    
    def append_elements(self):

        system = self.reported_object
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(self.header_text))

        display(SympyFormula( Eq(Symbol('F_d'),
                     dyn_sys.max_dynamic_force().doit() ), marker=None))
        
class GivenDataComponent(ReportComponent):

    title="Tabelka z wartościami parametrów do obliczeń"

    @property
    def entry_text(self):

        return "Przyjęte do obliczeń wartości poszczególnych parametrów przedstawia tabelka {tabelkamarker}"
    
    @property
    def table_caption(self):
        return "Podstawowe wartości parametrów"


    def append_elements(self):

        dyn_sys = self.reported_object
        
        given_data=dyn_sys._given_data

        if given_data == {}:
            params_dict = dyn_sys.get_reference_data()

        else: params_dict = given_data

        preview_mapper = lambda x: f'${latex(x)}$' if isinstance(x,(Symbol,Eq,Function,Mul)) else str(x)

        mapped_dict = {k: f'${latex(v)}$' if isinstance(v,(Symbol,Eq,Function,Mul)) else str(v) for k, v in params_dict.items()}

        params_dict_for_ldf={'Parametr':list(params_dict.keys()),'Wartość':list(mapped_dict.values())}
        params_ldf = pd.DataFrame(data=params_dict_for_ldf)
        params_ldf['Parametr'].apply(preview_mapper)
        params_ldf=params_ldf.set_index('Parametr')

        table_params=LatexDataFrame.formatted(params_ldf)

        display(ReportText(self.entry_text.format(tabelkamarker = AutoMarker(table_params))))

        display(table_params.reported(index=True, caption=self.table_caption))
