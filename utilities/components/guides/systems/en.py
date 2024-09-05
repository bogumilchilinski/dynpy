from  ..mechanics import *
from  ..mechanics import ReportComponent as BaseReportComponent


from . import pl

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from ....solvers.linear import *
from ....dynamics import *

months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']

from ..mechanics import display
import datetime


miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}

df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)

class ReportComponent(BaseReportComponent):



    @property
    def reported_object(self):

        from ....solvers.linear import ODESystem
        
        if isinstance(self._reported_object, ODESystem):
            return self._reported_object

        else:
            return self._reported_object
            
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
    guide_title  = 'WPROWADZENIE DO OBSŁUGI ŚRODOWISKA COCALC'
    
    
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

            self.append(NoEscape('\\Huge MANUAL DO KRÓTKIEGO KURSU PYTHONA\'A \n \n'))
            self.append(NoEscape(f'\\Huge {self.guide_title} \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))

            
            
            if len(system.q)==1:
                dof_str = 'JEDNYM STOPNIU SWOBODY'
            else:
                dof_str = 'WIELU STOPNIACH SWOBODY'
                
            if system._dissipative_potential==0 or system._dissipative_potential is None:
                damping_str = 'NIETŁUMIONE'
            else:
                damping_str = 'TŁUMIONE'

           
            #self.append(NoEscape(f'\\Large {damping_str} UKŁADY O {dof_str} \n \n'))    
            
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


#obsluga kolaka

class CocalcLoginComponent(ReportComponent):
    
    title="Logowanie do CoCalc"


    def append_elements(self):

        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('''Link do CoCalca: [LINK](https://cocalc.com/) '''))

        pic1 = Picture('./dynpy/utilities/components/guides/images/ColacSignIn.png', caption = "Przed Zalogowaniem",width='9cm')

        display(pic1)

        display(ReportText('Teraz należy zalogować się'))

        pic2 = Picture('./dynpy/utilities/components/guides/images/MetLogow.png', caption = "Metody logowania",width='9cm')

        display(pic2)

        pic3 = Picture('./dynpy/utilities/components/guides/images/Zalogowany.png', caption = "Zalogowany",width='9cm')

        display(pic3)

        display(ReportText('''Poprzez tego linka można dołączyć do projektu Ongoing: [ONGOING](https://cocalc.com/app?project-invite=dS62jsbRJvcfj2Mu) '''))
        

#Zakladanie Jupytera        

dynpy_imports_code=(
'''

from sympy import* 
from sympy.physics.mechanics import dynamicsymbols, init_vprinting
from sympy.abc import* 
init_vprinting()
from sympy import S
from dynpy.solvers.linear import ODESystem
import pandas as pd
from dynpy.models.mechanics import Pendulum
from dynpy.models import mechanics

from pylatex import Document, Section, Subsection, Itemize, Package,  HorizontalSpace, Description, Marker, Command
from pylatex.section import Paragraph, Chapter
from pylatex.utils import italic, NoEscape
from dynpy.utilities.adaptable import *
from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase, EngeneeringDrawingGuide,DevelopmentGuide
from dynpy.utilities.templates.tikz import TikzCaSCStandalone
from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock
from dynpy.utilities.report import (SystemDynamicsAnalyzer, DataPlot, AccelerationComparison, FFTComparison, ReportEntry,             SimulationalBlock, ReportText, SimulationFFT, DataStorage, Markdown, SummaryTable, SympyFormula, SymbolsDescription, DescriptionsRegistry,             ObjectCode,CurrentContainer)
from dynpy.models import mechanics
import inspect

''')

class JupyterSetUpComponent(ReportComponent):
    
    title="Zakładanie Jupytera"


    def append_elements(self):

        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        display(ReportText('Posługując sie Jupyterem będziemy tworzyć cele i wykonywali w nich pewne zadania.  '))

        pic9 = Picture('./dynpy/utilities/components/guides/images/PoDodaniuCeli.png', width='9cm')

        display(pic9)

        display(ReportText('Aby dodac nową cele wystarczy kliknąć na "prostokąt" i pod nim pojawią się opcje takie jak Code czy Text'))


        pic10 = Picture('./dynpy/utilities/components/guides/images/Markdown.png', width='9cm')

        display(pic10)

        display(ReportText('Po wybraniu opcji Text będziemy mogli napisać opis zadania, które wykonujemy lub nadać sekcji nagłówek  '))

        pic11 = Picture('./dynpy/utilities/components/guides/images/Plotki-w-naglowk.png', width='9cm')


        display(ReportText('Aby dodać nagłowek wystarczy wpisac znak Hashtag/płotek, ważne jest to że każdy nagłówek będzie traktowany jako sekcja niższego stopnia (\# główna sekcja, \#\# podsekcja itd...)'))

        display(pic11)

        pic12 = Picture('./dynpy/utilities/components/guides/images//Nazywanie-naglowka.png', width='9cm')

        display(pic12)


        display(ReportText('Po wejsciu w zakładkę Edit ukażą nam się poniższe opcje, takie jak usuwanie cel.'))

        pic13 = Picture('./dynpy/utilities/components/guides/images/Del-inne-opcje.png', width='9cm')

        display(pic13)


        display(ReportText('Czasami gdy Jupyter się zawiesi można spróbować go zrestartować uzywając poniższego przycisku'))

        pic14 = Picture('./dynpy/utilities/components/guides/images/Restart.jpeg', width='9cm')

        display(pic14)

        pic15 = Picture('./dynpy/utilities/components/guides/images/Restart-kernel.png', width='9cm')

        display(pic15)

        display(ReportText('Aby mieć pewność że w naszym kodzie uwzględnione zostały wszystkie zaimportowane biblioteki warto uruchamiając jupytera użyć opcji Run all cells pokazanej na poniższych obrazkach'))

        pic16 = Picture('./dynpy/utilities/components/guides/images/RunAll.jpeg', width='9cm')

        display(pic16)

        pic17 = Picture('./dynpy/utilities/components/guides/images/Run_all.png', width='9cm')

        display(pic17)

        display(ReportText('Aby elementy z których korzystamy działały, należy zaimportować potrzebne biblioteki, takie jak na załączonym obrazku oraz w formie tekstu który można skopiowac'))

        pic18 = Picture('./dynpy/utilities/components/guides/images/ImportowanieBibliotek.png', width='9cm')

        # display(pic18)

        display(GuideCode(dynpy_imports_code))


class CocalcFolderComponent(ReportComponent):
    
    title="Tworzenie Folderu i Jupytera"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('Tworzenie Folderu i Jupytera'))
        
        pic4 = Picture('./dynpy/utilities/components/guides/images/OtwarcieOngo.png', width='9cm')
        
        display(pic4)
        
        display(ReportText('W wyszukiwarce szukamy "Ongoing" i otwieramy projekt  '))
        
        pic5 = Picture('./dynpy/utilities/components/guides/images/PoWejsciuWOngo.png', width='9cm')
        
        display(pic5)
        
        display(ReportText('Otworzenie folderu UCZESTNICY  '))
        
        
        display(ReportText('Nadajemy nazwę folderu od naszego imieniem i nazwiskiem jak na zalączonym obrazku i klikamy create folder:'))
        
        pic6 = Picture('./dynpy/utilities/components/guides/images/NazwanieFolderu.png', width='9cm')
        
        display(pic6)
        
        
        display(ReportText('Po utworzeniu folderu klikamy przycisk NEW i utworzymy pierwszego Jupytera. W tym celu po nazwaniu Jupytera wybieramy zakladke Jupyter Notebook.'))
        
        
        pic7 = Picture('./dynpy/utilities/components/guides/images/WyborKernela.png', width='9cm')
        
        display(pic7)
        
        display(ReportText('Wybieramy Python 3 (system - wide)  '))
        display(ReportText('Po wybraniu Kernela utworzony zostanie nasz Jupyter.'))
        
        pic8 = Picture('./dynpy/utilities/components/guides/images/JupekUtworzony.png', width='9cm')
        
        display(pic8)
        ####

#podstawy raportowania
        



        
code_dynsys_list_str = '''
#komendy importujące biblioteki
import sympy 
from sympy import Symbol

from dynpy.models.mechanics.pendulum import Pendulum

#komenda wywołująca preview
Pendulum().preview()

#komenda wywołująca pomocne informacje dotyczące wahadła
help(Pendulum)

'''



class CocalcDynSysListComponent(ReportComponent):
    
    title="Lista systemów dynamicznych"


    def append_elements(self):

        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        display(ReportText(' Możemy również wyświetlić podgląd takiego wahadła używając poniższej funkcji:  '))
        display(Picture('./dynpy/utilities/components/guides/images/preview.jpg'))
        display(ReportText(' Co więcej, używając poniższej komendy help, możemy otrzymać wszelkie informacje na temat konkretnego modelu:  '))
        display(Picture('./dynpy/utilities/components/guides/images/help.jpg'))

        display(ReportText('''Kod wywołania pomocy jest następujący:'''))
        
        display(GuideCode(   code_dynsys_list_str   ))


        display(ReportText('Poniżej znajduje się lista systemów mechanicznych, których można użyć:  '))
        from dynpy.models import mechanics


        moduly=inspect.getmembers(mechanics, inspect.ismodule)


        for nazwa,modul in moduly:

            mods_tmp = [name for name,cls in inspect.getmembers(modul, inspect.isclass) if issubclass(cls,mechanics.principles.ComposedSystem)]
            classes_str=',\n - '.join(mods_tmp)
            if mods_tmp != []:
                display(ReportText( '\\par' + f'-{nazwa}'+ '\\par'))
                display(ReportText(f'- {classes_str}'))


class ReportingBasicsComponent(ReportComponent):
    
    title="Podstawy raportowania"


    def append_elements(self):
        from dynpy.utilities.templates.document import Guide
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        display(ReportText('Pierwszy krokiem jest import "szablonu" dokumentu np. klasy `Guide` oraz wywołanie metody `base_setup`. Wynik wywołania jest następujący:'))

        display(GuideCode(
'''
from dynpy.utilities.templates.document import *
Guide.base_setup()
'''
        ))
        
        display(Guide.base_setup())

        display(ReportText('Metoda zwraca przykładowe cele gotowe do przekopiowania do notatnika.'))
                
        

        
#obsluga systemow dynamicznych


imports_code_str = '''
from sympy import *
import numpy as np
from pandas import *
from sympy.physics.mechanics import init_vprinting, dynamicsymbols
from pylatex import Document, Section, Subsection, Itemize, Package, HorizontalSpace, Description, Marker, Command 
from pylatex.section import Paragraph, Chapter 
from pylatex.utils import italic, NoEscape 
from dynpy.utilities.adaptable import *
from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase 
from dynpy.utilities.templates.tikz import TikzCaSCStandalone
from dynpy.utilities.report import (Markdown, Picture, SympyFormula, SymbolsDescription, DescriptionsRegistry, ObjectCode, CurrentContainer)
from dynpy.solvers.linear import ODESystem
init_vprinting()
from dynpy.models.mechanics import SDOFWinchSystem

'''



class DynamicSystemCallComponent(ReportComponent):
    
    title="Introduction"


    def append_elements(self):

        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText('Required library/class imports to create the report:'))


        display(GuideCode(imports_code_str.replace('SDOFWinchSystem',system_name)))


        display(ReportText('A class that allows you to solve a specific problem (in this case, a problem related to the description of a crane):'))

        display(GuideCode("""from dynpy.models.mechanics import SDOFWinchSystem""".replace('SDOFWinchSystem',system_name)      ))
        
        
        display(GuideCode(f'system=dyn_sys={system_name}()'  ))

        display(ReportText('Path to the sought-after class on the CoCalc platform:'))
        #display(Picture('./dynpy/utilities/components/guides/images/sciezka_w.jpg'))

        display(ReportText((system.__class__.__module__)))
        
        display(ReportText('The way to call out the preview of the class - the so-called preview:'))

        display(GuideCode(f'''{system_name}()._as_picture()'''  ))
        display(system._as_picture())
        
        
class DynamicSystemMethodsUsageComponent(ReportComponent):

    title="Wywolywanie rownan ruchu oraz innych metod"



    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        #from dynpy.models.mechanics import ForcedSpringMassSystem as SDOFWinchSystem
        
        eoms=system._eoms[0]

        display(ReportText('Proces wywołowania równania ruchu za pomoca metody eoms:'))
        display(GuideCode('''eoms=SDOFWinchSystem()._eoms[0]'''.replace('SDOFWinchSystem',system_name)))



        display(ReportText('Wynik jest następujący:'))
        #display(SympyFormula(eoms))


        eoms_eq=Eq(eoms,0)
        eoms_eq

        display(ReportText('Proces tworzenia równania z wcześniej wywołanego równania ruchu (równanie początkowo jest wywoływane jako ciąg znaków reprezentujący opis danego problemu, jednak nie jest on wówczas równaniem, które możemy zastosować do obliczeń):'))
        display(GuideCode('eoms_eq=Eq(eoms,0)'.replace('SDOFWinchSystem',system_name)))
        display(ReportText('Wynik jest następujący:'))
        display(SympyFormula(eoms))

        #display(Picture('./Image/eomseq_w.jpg'))
        
        slownik=system.get_random_parameters()
        
        slownik_numerical=system.get_numerical_parameters()
        
        display(ReportText('Słownik z losowymi wartościami parametrów symbolicznych systemu wygląda następująco:'))
        
        display(SympyFormula(Dict(slownik)))
        
        display(ReportText('Słownik z losowymi wartościami liczbowymi do analizy numerycznej systemu, wygląda następująco:'))
        
        display(SympyFormula(Dict(slownik_numerical)))

        display(ReportText('Tworzenie słownika (w tym przypadku z losowymi wartościami) umożliwiający nadanie wartości zmiennym:'))
        
        

        display(GuideCode('slownik=SDOFWinchSystem().get_random_parameters()'.replace('SDOFWinchSystem',system_name)))
        display(GuideCode('slownik_numerical=SDOFWinchSystem().get_numerical_parameters()'.replace('SDOFWinchSystem',system_name)))
        #display(Picture('./Image/slownik_w.jpg'))

        display(ReportText('Rozwiązanie ogólne równania:'))
        #display(Picture('./Image/ogolne_w.jpg'))

        solution=system._ode_system.solution[0]
        display(GuideCode('''solution=SDOFWinchSystem()._ode_system.solution[0]'''.replace('SDOFWinchSystem',system_name)))
        display(GuideCode('solution_sym=SDOFWinchSystem()._ode_system.solution'.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(solution))

        display(ReportText('Rozwiązanie szczególne równania:'))



        steady_solution=system._ode_system.steady_solution[0]
        display(GuideCode('''steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(steady_solution))
    
    
    
steady_sol_str=(
'''

slownik_numerical=system.get_numerical_parameters()

steady_solution=SDOFWinchSystem()._ode_system.subs(slownik_numerical).steady_solution
steady_solution_subs=steady_solution

''')



solution_table_str=(
'''
coord = system.q[0]
table=table_eq.numerized().compute_solution(t_span)
table_new=table[[coord]]

''')

report_table_define_str=(
'''

tabela = Section('Tabela')
CurrentContainer(tabela)
coord = system.q[0]
solution_sym=ForcedSpringMassSystem()._ode_system.subs(slownik_numerical).steady_solution

general_sol_matrix=solution_sym
display(general_sol_matrix)
table_eq=general_sol_matrix.numerized(slownik_numerical).compute_solution(t_span)
table=table_eq

table_for_report = table.iloc[0:10].to_latex_dataframe()

table_for_report.reported(caption=('Tabela'))


''')

class SimulationsComponent(ReportComponent):
    
    title="Creating graphs"

    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__


        display(ReportText('Definition of time vector:'))
        #display(Picture('./Image/tspan_w.jpg'))

        display(GuideCode('''t_span = np.linspace(0,100,200)'''.replace('SDOFWinchSystem',system_name)))
        t_span = np.linspace(0,100,200)


        display(ReportText('The way to create charts manually:'))
        display(ReportText('Substitution of data from a previously created dictionary:'))
        #display(Picture('./Image/stedisub_w.jpg'))

        
        
        display(GuideCode(steady_sol_str.replace('SDOFWinchSystem',system_name)))
        #dict=system.get_random_parameters()
        dict_numerical=system.get_numerical_parameters()
        
        steady_solution=system._ode_system.steady_solution.rhs[0]
        steady_solution_subs=steady_solution.subs(dict_numerical)
        display(steady_solution_subs)
        display(ReportText('Using the lambdify method to convert functions to anonymous form'))
        #display(Picture('./Image/lambdify_w.jpg'))


#         eq_lambdify=lambdify(system.ivar, list(steady_solution_subs.rhs))
        eq_lambdify=lambdify(system.ivar, steady_solution_subs)
        display(GuideCode('''eq_lambdify=lambdify(system.ivar, list(steady_solution_subs.rhs))'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(eq_lambdify))


        display(ReportText('Creating a graph:'))
        #display(Picture('./Image/plot2_w.jpg'))
        
        #pd.DataFrame(data=eq_lambdify(t_span),index=t_span).plot()
        #plt.show
        df = pd.DataFrame(data=eq_lambdify(t_span))
        display(df.plot())
        display(GuideCode('''df = pd.DataFrame(data=eq_lambdify(t_span))
df.T.set_index(t_span).plot()'''.replace('SDOFWinchSystem',system_name)))


        display(ReportText('A way to create charts directly from a previously called class:'))
        display(ReportText('Dictionary substitution:'))
        #display(Picture('./Image/table_in.jpg'))

        
        #zmiana slownik na slownik_numerical
        table_eq=system._ode_system.steady_solution.subs(dict_numerical)


        display(GuideCode('''dict_numerical=system.get_numerical_parameters()
table_eq=system._ode_system.steady_solution.subs(dict_numerical)
table_eq'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(table_eq))

        #jak tabele zrobic + ew jej output
        display(ReportText('The way to create a table:'))
        #display(Picture('./Image/table_in.jpg'))
        coord = system.q[0]
        table=table_eq.numerized().compute_solution(t_span)
        table_new=table[coord]

        display(GuideCode(solution_table_str.replace('SDOFWinchSystem',system_name)))
        display(table)
        display(table_new)

        display(ReportText('The way to create a visualization of the defined table:'))
        #display(Picture('./Image/table_in.jpg'))

##### nazwy image tez przetlumaczone
        display(ReportText('How to create the table displayed in the finished report:'))
        #display(Picture('./Image/creating_table_in.png'))

        display(GuideCode(report_table_define_str.replace('SDOFWinchSystem',system_name)))

        display(ReportText('Visualization of the graph:'))
        #display(Picture('./Image/tikzplot_w.jpg'))

        display(GuideCode('''table_new.to_pylatex_tikz().in_figure()'''.replace('SDOFWinchSystem',system_name)))    

section_development_str=(
'''
from dynpy.utilities.report import ReportText
intro=Section('Introduction')
CurrentContainer(intro)

#scheme=Picture('picture_winch.png')
#display(scheme)

picture = Picture(system._as_picture().__dict__['image'], caption = f'System diagram')
display(picture)
''')

forming_equation_str=(
'''
eom=Section('Equations of motion')
CurrentContainer(eom)
display(ReportText('Equations of motion: '))
eoms=system._eoms[0]

slownik=system.get_random_parameters()
slownik_numerical=system.get_numerical_parameters()

eoms_eq=Eq(eoms,0)
solution_subs = system._ode_system.steady_solution.subs(slownik_numerical)
steady_solution=system._ode_system.steady_solution[0]
steady_solution_subs=steady_solution.subs(slownik_numerical)
#zmiana slownik na slownik numerical

steady_solution_subs
display(solution_subs)

ode_system = system._ode_system #ESystem.from_dynamic_system(system.linearized())
param_dict = system.get_numerical_parameters()
t_span = np.linspace(1,10,1001)
ic_list = [0.0,0.0] ### Dla układu o jednym stopniu swobody
ode_solution = ode_system.subs(param_dict).steady_solution

display(ode_solution)

display(SympyFormula(eoms_eq))
display(SympyFormula(solution_subs))
display(SympyFormula(steady_solution_subs))
''')

figure_eoms=(
'''
wykres=Section('Graph of dynamic system ')
CurrentContainer(wykres)
display(ReportText('The graph shows the solution of the equations of motion for the dynamic system under consideration. The graph shows the dependence of the displacement (or other relevant quantity, depending on the type of system) as a function of time, obtained from the solution of the differential equations.'))

ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)
picture=ode_simulation.to_pylatex_tikz(subplots=True).in_figure(caption='Graph of the solution of the equations of motion for the dynamic system under consideration.')

display(picture)

display('Equations of motion:')

display(eoms_eq)
display(solution_subs)
display(steady_solution_subs)
''')

report_generation_libraries=(
'''
doc_final = MechanicalCase('./output/Document',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)

doc_final.append(intro)
doc_final.append(eom)
doc_final.append(wykres)

doc_final.generate_pdf()
''')


class SimulationReportComponent(ReportComponent):
    
    title="Raportowanie wykresów"

    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText('Sposób tworzenia sekcji przy użyciu metody CurrentContainer:'))

        #sekcja 1|
        display(GuideCode(section_development_str))


        display(ReportText('Poniżej przedstawiony jest sposób tworzenia równań. W tym przypadku przyrównamy wyrazy (przed przecinkiem) do 0 (po przecinku).'))
        display(GuideCode(forming_equation_str.replace('SDOFWinchSystem',system_name)))



        display(ReportText('poniżej przedstawiony jest sposób tworzenia wykresu dla podanych danych: '))
        
        eoms=system._eoms[0]
        
        slownik=system.get_random_parameters()
        slownik_numerical=system.get_numerical_parameters()
        
        eoms_eq=Eq(eoms,0)
        solution_subs = system._ode_system.steady_solution.subs(slownik_numerical)
        steady_solution=system._ode_system.steady_solution[0]
        steady_solution_subs=steady_solution.subs(slownik_numerical)
        #zmiana slownik na slownik numerical
        
        steady_solution_subs
#         eq_sol=Eq(solution_subs,0)
#         eq_steady_sol=Eq(steady_solution_subs,0)
        display(SympyFormula(solution_subs))

        ode_system = system._ode_system #ESystem.from_dynamic_system(system.linearized())
        param_dict = system.get_numerical_parameters()
        t_span = np.linspace(1,10,1001)
        ic_list = [0.0,0.0] ### Dla układu o jednym stopniu swobody
        ode_solution = ode_system.subs(param_dict).steady_solution
        display(SympyFormula(ode_solution))

        display(GuideCode(figure_eoms))
        
        ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)
        #de_simulation = steady_solution_subs.compute_solution(t_span, ic_list = ic_list)
        display(ReportText('Wykres: '))
        picture=ode_simulation.to_pylatex_tikz(subplots=True).in_figure(caption='Otrzymany wykres')

        display(picture)

        display(ReportText('równania ruchu: '))

        display(SympyFormula(eoms_eq))
        display(SympyFormula(solution_subs))
        display(SympyFormula(steady_solution_subs))




        display(GuideCode(report_generation_libraries))



        display(ReportText('Tak wygląda ostateczny raport:'))
        display(Picture('./dynpy/utilities/components/guides/images/sekcja3_w.jpg'))
        

class DynSysCodeComponent(ReportComponent):

    title="Calling dynamic system code"



    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText('''The code of the analyzed class is as follows'''))
        display(ObjectCode(system.__class__))
        display(ReportText(f'''The presented code of the {system_name} class was retrieved with the help of the ObjectCode class.  '''))
        
        
#pandas guide
data_code=(
'''
months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']
average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]
day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]
Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]
Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]
''')

data_dict_code=(
'''
data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}
''')
output_code=(
'''
import pandas as pd

df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)
df
''')

class PandasTableGenerationComponent(ReportComponent):
    
    title="Generating a table"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Creating a list of data:'))

        display(GuideCode(data_code))


        display(ReportText('Creating dictionary:'))

        display(GuideCode(data_dict_code))


        display(ReportText('Calling the table:'))
        display(GuideCode(output_code))
        display(df)
        

winch_pandas_code=(
'''

from dynpy.models.mechanics.tmac import SDOFWinchSystem
eoms= SDOFWinchSystem()._eoms[0]
eoms_eq_raw=Eq(eoms,0)
data=SDOFWinchSystem().get_random_parameters()
Ms0=Symbol('M_s0')
b=Symbol('b')

winch = SDOFWinchSystem()

Ms=winch.M_s
t = winch.ivar

dane={
    Ms0:4700,
    b:1500,
    SDOFWinchSystem().m0:10}



winch_up=SDOFWinchSystem().subs(Ms,Ms0-b*SDOFWinchSystem().phi.diff()).subs(data).subs(dane)
nowe_ode=winch_up._ode_system
solution=nowe_ode.solution
steady_sol=nowe_ode.steady_solution[0]
t_span = np.linspace(0,100,200)
phi=winch.q[0]



units_dict = {phi : ureg.radian,phi.diff(t):ureg.radian,t:ureg.second }
LatexDataFrame.set_default_units(units_dict)
''')        
        
class BasicSymComponent(ReportComponent):
    
    title="Symulacje - wariant podstawowy"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        SDOFWinchSystem = type(system)
#         from dynpy.models.mechanics.tmac import SDOFWinchSystem
        eoms= SDOFWinchSystem()._eoms[0]
        eoms_eq_raw=Eq(eoms,0)
        data=SDOFWinchSystem().get_random_parameters()
        Ms0=Symbol('M_s0')
        b=Symbol('b')

        winch = SDOFWinchSystem()

        Ms=winch.M_s
        t = winch.ivar

        dane={
            Ms0:4700,
            b:1500,
            SDOFWinchSystem().m0:10}



        winch_up=SDOFWinchSystem().subs(Ms,Ms0-b*SDOFWinchSystem().phi.diff()).subs(data).subs(dane)
        nowe_ode=winch_up._ode_system
        solution=nowe_ode.solution
        steady_sol=nowe_ode.steady_solution[0]
        t_span = np.linspace(0,100,200)
        phi=winch.q[0]



        units_dict = {phi : ureg.radian,phi.diff(t):ureg.radian,t:ureg.second }
        LatexDataFrame.set_default_units(units_dict)
        display(GuideCode(winch_pandas_code))


        display(ReportText('Analiza i symulacje klasy SDOFWinchSystem, tzw. wciągarki. Równania ruchu: '))
        display(SympyFormula(eoms_eq_raw))
        display(ReportText(f'Po podstawieniu danych do ogólnej postaci (do równania {AutoMarker(eoms_eq_raw)}) otrzymuje się: '))
        display(SympyFormula(Eq(winch_up._eoms[0],0)))
        display(ReportText('Rozwiązanie szczególne: '))
        display(SympyFormula(Eq(phi,steady_sol)))
        display(ReportText('Rozwiązanie ogólne: '))
        display(SympyFormula(Eq(phi,solution[0])))
        general_sol_matrix=nowe_ode.solution.with_ics([0,0])
        table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)
        table=table_eq[[phi.diff(t)]]
        table_for_report = table.iloc[0::15].to_latex_dataframe()
        display(GuideCode(f'Wyniki symulacji zostały zestawione w tabeli {AutoMarker(table_for_report)}'))
        display(table_for_report.reported(caption='Tabela'))
        sim_plot = table.to_pylatex_tikz().in_figure(caption='Wykres')
        display(GuideCode(f'Wyniki symulacji zostały przedstawione graficznie na wykresie {AutoMarker(sim_plot)}'))
        display(sim_plot)
        

class PandasMethodsComponent(ReportComponent):
    
    title="Pandas library methods"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        
        months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']

        average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

        Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

        Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

        day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

        data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}

        df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)

        display(ReportText('Iloc method:'))
        display(GuideCode('df.iloc[0:3]'))
        display(df.iloc[0:3])

        # display(Picture('./Image/iloc.jpg', caption = ""))

        display(ReportText('Loc method:'))
        display(GuideCode('df.loc["January"]'))
        display(df.loc["January"])
        # display(Picture('./Image/loc.jpg', caption = ""))

        display(ReportText('The set axis method works for both columns and rows. For rows:'))
        display(GuideCode("df.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'],axis = 'index')"))
        display(df.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'], axis = 'index'))

        # display(Picture('./Image/set_axis.jpg', caption = ""))

        display(ReportText('For columns:'))
        display(GuideCode("df.set_axis(['a','b','c'],axis = 'columns')"))
        display(df.set_axis(['a','b','c'],axis = 'columns'))
        # display(Picture('./Image/set_axis2.jpg', caption = ""))

        display(ReportText('Metoda rename - Zmienienie pojedynczej kolumny:'))
        display(GuideCode('''df.rename(columns = {"Day length in a month [h]":'AAA'}, index = {"January":'A'} )'''))
        display(df.rename(columns = {"Day length in a month [h]":'AAA'}, index = {"January":'A'} ))
        # display(Picture('./Image/rename.jpg', caption = ""))

        display(ReportText('Map method:'))
        display(GuideCode('df.map(lambda x:x+2)'))
#         display(df.map(lambda x:x+2)) - nowa wersja pandasa zmieniła nazwę "map" na "applymap" - Michał SZ
        display(df.map(lambda x:x+2))

        # display(ReportText('***Applymap mthod for colummn/row:***'))
        # display(Markdown('''
        # Remeber, that df[NAME] => creates series, a df[[NAME]] => creates new frame, which can be used with .applymap'''))
        # display(df[['Day length in a month']].applymap(lambda x:x+100))
        # pd_h = df[['Day length in a month']].applymap(lambda x:x+100)
        # list1 = pd_h['Day length in a month'].tolist()
        # df['Day length in a month'] = list1

        # display(df)

        # display(Picture('./Image/applymap.jpg', caption = ""))

        display(ReportText('Slicing:'))
        display(GuideCode('''df['Day length in a month [h]']'''))
        display(df['Day length in a month [h]'])
        # display(Picture('./Image/slice2.jpg', caption = ""))
        
        
steady_state_code=(
'''
steady_state = SDOFWinchSystem()._ode_system.solution[0]
steady_state
''')
    
list_code=(
'''
Model_number = [1,2,3,4] 
inertia_list= [1,2,3,4] 
stiffness_list = [1,2,3,4] 
damping_list = [1,2,3,4] 
Excitation_list = [1,2,3,4]

system_parameters = {'Inertia':inertia_list, 'Stiffness':stiffness_list,'Damping':damping_list,'Excitation':Excitation_list}
''')

df_code=(
'''
df = pd.DataFrame(index = Model_number, data = system_parameters)
df
''')

dict_code=(
'''
data_dict_1 = df.T[1].to_dict()
data_dict_2 = df.T[2].to_dict()
data_dict_3 = df.T[3].to_dict()
data_dict_4 = df.T[4].to_dict()
eq1 = steady_state.subs(data_dict_1).subs(SDOFWinchSystem().F,100)
eq2 = steady_state.subs(data_dict_2).subs(SDOFWinchSystem().F,100)
eq3 = steady_state.subs(data_dict_3).subs(SDOFWinchSystem().F,100)
eq4 = steady_state.subs(data_dict_4).subs(SDOFWinchSystem().F,100)
''')

lambdify_code_str=(
'''
function = lambdify(SDOFWinchSystem.ivar, steady_state)
t_span = np.linspace(0,100,200)
''')

sym_gen_sol=(
'''
table = SDOFWinchSystem().subs(data_dict_2).subs(SDOFWinchSystem.F,100).numerized().compute_solution(t_span,[0,0]).iloc[:,0] 
table.plot() 
''')

list_elem = (
'''
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
''')


solution_code=(
'''
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
pd.DataFrame(index = t_span, data = eq2_fun(t_span).plot()
''')



class DifferentSimulationsComponent(ReportComponent):
    
    title="Simulations - other variants"



    
    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        SDOFWinchSystem = type(system)

        display(ReportText('Generate a class particular solution:'))
        display(GuideCode(steady_state_code.replace('SDOFWinchSystem',system_name)))
        steady_state = SDOFWinchSystem.from_random_data()._ode_system.solution.with_ics([0,0])[0].subs(SDOFWinchSystem.m0,10)
        display(steady_state)

        # display(Picture('./Image/steady_sol.jpg', caption = ""))


        display(ReportText('Creating a list and a table'))
        display(GuideCode(list_code))
        Model_number = [1,2,3,4] 
        inertia_list= [1,2,3,4] 
        stiffness_list = [1,2,3,4] 
        damping_list = [1,2,3,4] 
        Excitation_list = [1,2,3,4]

        system_parameters = {'Inertia':inertia_list, 'Stiffness':stiffness_list,'Damping':damping_list,'Excitation':Excitation_list}


        display(GuideCode(df_code))

        df1 =  pd.DataFrame(index = Model_number, data = system_parameters)
        display(df1)

        # display(Picture('./Image/listatab.jpg', caption = ""))


        display(ReportText('Create a vector - a set of time values for our simulation: '))
        display(GuideCode(lambdify_code_str.replace('SDOFWinchSystem',system_name)))

        function = lambdify(SDOFWinchSystem.ivar, steady_state)
        t_span = np.linspace(0,100,200)

        # display(Picture('./Image/wektor1.jpg', caption = ""))


        display(ReportText('Creating data dictionaries and substituting the data in it into our equation so that it only depends on time: '))

        display(GuideCode(dict_code.replace('SDOFWinchSystem',system_name)))
        # display(Picture('./Image/slowniki_sym.jpg', caption = ""))


        display(ReportText('Create a list of elements by substituting time values one by one: '))
        display(GuideCode(list_elem.replace('SDOFWinchSystem',system_name)))
        display(Picture('./dynpy/utilities/components/guides/images/eq2_fun.jpg', caption = ""))


        display(ReportText('Generate a graph of a particular equation: '))
        display(GuideCode(solution_code.replace('SDOFWinchSystem',system_name)))
        display(Picture('./dynpy/utilities/components/guides/images/wykres_niefun.jpg', caption = ""))


        display(ReportText('Generate a graph of the general equation: '))
        display(GuideCode(sym_gen_sol.replace('SDOFWinchSystem',system_name)))
        display(Picture('./dynpy/utilities/components/guides/images/ogolne_tab.jpg', caption = ""))
        display(Picture('./dynpy/utilities/components/guides/images/ogolne_plot.jpg', caption = ""))

        
        
pandas_latax_df_code = '''

from dynpy.utilities.adaptable import LatexDataFrame
from sympy import symbols
from pint import UnitRegistry 
ureg = UnitRegistry()

a,b,c,d,e = symbols('a b c d e',positive = True)

units_dict = {a:ureg.meter,b:ureg.second, 'lubuskie':ureg.meter}
LatexDataFrame.set_default_units(units_dict)

LatexDataFrame.formatted(df.set_axis([a,b,c,d,e],axis='columns'))
'''


class BasicOperationsComponent(ReportComponent):

    title="Basic objects and operations"

    def append_elements(self):

        ss1 = pd.Series([2,34,5,-2,8,12])
        x1 = LatexDataFrame(ss1)

        ss2=ss1*10
        x2 = LatexDataFrame(ss2)
        ss3=ss1.abs()

        x3 = LatexDataFrame(ss3)
        ss4=ss1.describe()
        x4 = LatexDataFrame(ss4)

        s5=ss1
        s5.index =['first','second','third','fourth','fifth','sixth']
        x5 = LatexDataFrame(s5)

        voivodeships_list = ['Masovian','Greater Poland','Lublin','Warmian-Masurian','West Pomeranian','Podlaskie','Lower Silesian','Pomeranian','Łódź','Kuyavian-Pomeranian','Subcarpathian','Lesser Poland','Lubusz','Silesian','Holy Cross','Opole']
        area_list=[35558,29826,25122,24173,22892,20187,19947,18310,18219,17972, 17846,15183,13988,12333,11711,9412]
        population_list=[5349114,3475323,2139726,1439675,1710482,1188800,2904207,2307710,2493603,2086210, 2127657,3372618,1018075,4570849,1257179,996011]
        income_list=[4464,3371,3338,3447,3649,3552,3788,4104,3438,3508,3212, 3395,3187,3586,3268,3112]
        expenses_list=[787,596,623,597,767,697,742,1023,590,778,574,598,365,631,647,431]
        voivodeships_data={'area':area_list,'no. of people':population_list,'income':income_list,'expenses':expenses_list}
        
        df = pd.DataFrame(index=voivodeships_list,data=voivodeships_data)
        x6 = LatexDataFrame(df)

        area = Symbol('A',positive=True)
        df['population']=df['no. of people']/df['area']
        x7 = LatexDataFrame(df)
        
        df21=df.iloc[1:3,0:3]
        x8 = LatexDataFrame(df21)
        
        df22=df.loc['Masovian']
        x9 = LatexDataFrame(df22)
        
        df.loc['Warsaw',:]=[4000,600000,2500,300,120]
        x10 = LatexDataFrame(df)
        
        df30=df.set_axis([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],axis='index')
        x11 = LatexDataFrame(df30)
        
        a,b,c,d,e = symbols('a b c d e',positive = True)

        units_dict = {a:ureg.meter,b:ureg.second, 'Lubusz':ureg.meter}
        LatexDataFrame.set_default_units(units_dict)

        df31=df.set_axis([a,b,c,d,e],axis='columns')
        x12 = LatexDataFrame(df31)
        
        df40=df.rename(columns={'no. of people':'POPULATION'},index={'Lubusz':'LUBU...'})
        x21 = LatexDataFrame(df40)
        
        df50=df.income.apply(lambda x: x*2)
        x13 = LatexDataFrame(df50)
        df60=df.map(lambda x: x*2)
        x14 = LatexDataFrame(df60)
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        display(ReportText('Importing the Pandas library'))
        display(GuideCode('import pandas as pd'))

        display(ReportText('Create a series ( in analogy to Excel columns)'))
        display(GuideCode('s = pd.Series([2,34,5,-2,8,12])'))
        display(x1.reported())

        display(ReportText('Operations on series'))
        display(GuideCode('s*10'))
        display(x2.reported())

        display(ReportText('Abs() function - absolute value'))
        display(GuideCode('s.abs()'))
        display(x3.reported())

        display(ReportText('Describe() function - basic statistics'))
        display(GuideCode('s.describe()'))
        display(x4.reported())

        display(ReportText('Changing the names of indexes that start with 0 by default'))
        display(GuideCode('''s.index =['first','second','third','fourth','fifth','sixth']'''))
        display(x5.reported())
        
        display(ReportText('DataFrame - a DataFrame table that is created from a dictionary'))
        display(ReportText('Creating a list of voivodeships:'))
        display(GuideCode('''
voivodeships_list = ['Masovian','Greater Poland','Lublin','Warmian-Masurian','West Pomeranian','Podlaskie','Lower Silesian','Pomeranian','Łódź','Kuyavian-Pomeranian','Subcarpathian','Lesser Poland','Lubusz','Silesian','Holy Cross','Opole']
            '''))
        
        display(ReportText('Creating a list containing the area of each voivodeship in the ${km^2}$:'))
        display(GuideCode('area_list=[35558,29826,25122,24173,22892,20187,19947,18310,18219,17972,17846,15183,13988,12333,11711,9412]'))

        display(ReportText('Creating a list of population numbers for each voivodeship:'))
        display(GuideCode('population_list=[5349114,3475323,2139726,1439675,1710482,1188800,2904207,2307710,2493603,2086210,2127657,3372618,1018075,4570849,1257179,996011]'))
        
        display(ReportText('Creating a list of per capita incomes in a province:'))
        display(GuideCode('income_list=[4464,3371,3338,3447,3649,3552,3788,4104,3438,3508,3212,3395,3187,3586,3268,3112]'))

        display(ReportText('Creating a list of property/investment expenses per person in the province:'))
        display(GuideCode('expenses_list=[787,596,623,597,767,697,742,1023,590,778,574,598,365,631,647,431]'))

        display(ReportText('Creating a dictionary'))
        display(GuideCode('''voivodeships_data={'area':area_list,'no. of people':population_list,'income':income_list,'exp.':expenses_list}'''))

        display(ReportText('Calling a table'))
        display(GuideCode('''
df = pd.DataFrame(index=voivodeships_list,data=voivodeships_data)
df
        '''))
        display(x6.reported())

        display(ReportText('Operations on columns - adding a column with calculated population concentration'))
        display(GuideCode('''
df['popul.']=df['no. of people']/df['area']
df
        '''))
        display(x7.reported())

        display(ReportText('The iloc method - we use indexes to perform operations on rows and columns'))
        display(GuideCode('df.iloc[1:3,0:3]'))
        display(x8.reported())

        display(ReportText('Loc method - using labels, we perform operations on rows and columns'))
        display(GuideCode('''df.loc['Masovian']'''))
        display(x9.reported())

        display(ReportText('Adding a line using the loc method'))
        display(GuideCode('''
df.loc['Warsaw',:]=[4000,600000,2500,300,120]
df
        '''))
        display(x10.reported())

        display(ReportText('Set_axis method for rows'))
        display(GuideCode('''df.set_axis([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],axis='index')'''))
        display(x11.reported())

        display(ReportText('Set_axis method for columns'))
        display(GuideCode('''df.set_axis(['a','b','c','d','e'],axis='columns')'''))
        display(x12.reported())

        display(ReportText('Rename method - renames columns and rows'))
        display(GuideCode('''
df.rename(columns={'no. of people':'POPULATION'},index={'Lubusz':'LUBU...'})
        '''))
        display(x21.reported())

        display(ReportText('Apply method - implementation of a function on a selected table column'))
        display(GuideCode('df.income.apply(lambda x: x*2)'))
        display(x13.reported())

        display(ReportText('Map method - implementation of functions across the table'))
        display(GuideCode('df.map(lambda x: x*2)'))
        display(x14.reported())
        
class DynSysIntroComponent(ReportComponent):

    title="Introduction to using the DynSys module"
    
    def append_elements(self): 
        from dynpy.models.mechanics.tmac import SDOFWinchSystem
        
        display(ReportText('Required imports of libraries/individual classes to create the report:'))
        display(ObjectCode(
        '''
                from sympy import *
                import numpy as np
                from pandas import *
                from sympy.physics.mechanics import init_vprinting, dynamicsymbols
                from pylatex import Document, Section, Subsection, Itemize, Package, HorizontalSpace, Description, Marker, Command 
                from pylatex.section import Paragraph, Chapter 
                from pylatex.utils import italic, NoEscape 
                from dynpy.utilities.adaptable import *
                from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase 
                from dynpy.utilities.templates.tikz import TikzCaSCStandalone
                from dynpy.utilities.report import (SystemDynamicsAnalyzer, DataPlot, AccelerationComparison, FFTComparison, ReportEntry, SimulationalBlock,               ReportText, SimulationFFT, DataStorage, Markdown, SummaryTable, Picture, SympyFormula, SymbolsDescription, DescriptionsRegistry,
                ObjectCode, CurrentContainer)
                from dynpy.solvers.linear import ODESystem
                init_vprinting()
                from dynpy.models.mechanics.tmac import SDOFWinchSystem

        '''))
        display(ReportText('A class that allows you to solve a specific problem (in this case, a problem related to the description of a crane):'))
        display(Markdown('''
            from dynpy.models.mechanics.tmac import SDOFWinchSystem
            '''))

        display(ReportText('The way to call the class preview:'))

        display(ObjectCode('''
            SDOFWinchSystem().preview()
            '''))

  
        eoms=SDOFWinchSystem()._eoms[0]
        eoms

        display(ReportText('The process of calling the equation of motion using the eoms method:'))
        display(ObjectCode(
            '''
            eoms=SDOFWinchSystem()._eoms[0]
            eoms
            '''))

        display(ReportText('The result is:'))
        display(SympyFormula(eoms))


        eoms_eq=Eq(eoms,0)
        eoms_eq
        display(ReportText('The process of creating an equation from a previously called equation of motion (the equation is initially called as a string representing a description of the problem at hand, but it is not then an equation that we can use for calculations):'))
        display(ObjectCode('''
            eoms_eq=Eq(eoms,0)
            eoms_eq
            '''))
        display(ReportText('The result is:'))
        display(SympyFormula(eoms))

        slownik=SDOFWinchSystem().get_random_parameters()


        display(ReportText('Create a dictionary (in this case with random values) that allows you to assign values to variables:'))
        display(ObjectCode('''
            slownik=SDOFWinchSystem().get_random_parameters()
            slownik
            '''))
        display(slownik)

        display(ReportText('General solution of the equation:'))

        solution=SDOFWinchSystem()._ode_system.solution[0]
        display(ObjectCode(
            '''
                solution=SDOFWinchSystem()._ode_system.solution[0]
                solution

            '''))
        display(SympyFormula(solution))

        display(ReportText('Special solution of the equation:'))



        steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]
        display(ObjectCode(
            '''
                steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]
                steady_solution

            '''))
        display(SympyFormula(steady_solution))
        
##ODESystem wprowadzenie


obiekt_1=(
'''
from sympy import *
from dynpy.solvers.linear import ODESystem
from sympy.physics.mechanics import dynamicsymbols

omega = Symbol('\omega',positive = True)
Omega = Symbol('\Omega',positive = True)

t = Symbol('t')
x = Function('x')(t)
'''
)
obiekt_2=(
'''
m = Symbol('m',positive = True)
c = Symbol('c',positive = True)
k = Symbol('k',positive = True)
t=Symbol('t')
x=dynamicsymbols('x')

eq1 = Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0) \n

odesys1 = ODESystem(odes = Matrix([eq1.lhs-eq1.rhs]),dvars = Matrix([x], ode_order=2)) \n

odesys1.solution
''')

class BasicUsageOfODESystemComponent(ReportComponent):
    
    title="Podstawy używania klasy ODESystem"


    def append_elements(self):

        display(ReportText('Klasa ODESystem służy do rozwiązywania i analizowania równań różniczkowych. Tworzy ona obiekt w postaci równania różniczkowego zwyczajnego (ODE) w formie ogólnej i zapewnia metody podstawowych operacji, takich jak rozwiązanie ogólne czy szczególne.'))

        display(ReportText(
        '\n Do prawidłowego wykonania zadania należy przyjąć następujące kroki:  \n- Definiowanie potrzebnych symboli, które wykorzystamy do budowy ODESystem: '))

        display(ObjectCode(obiekt_1))
        display(ReportText('- Definiowanie równania, które wykorzystamy do budowy ODESystem: '))

        display(ObjectCode('''eq = Eq(x.diff(t, t) + omega**2 * x,sin(Omega*t))'''))


        display(ReportText('- Definiowanie ODESystem, gdzie " odes " to przeniesione na jedną stronę równanie, które będziemy rozwiązywać (najlepiej podawać w macierzy), a " dvars " to zmienne, po których rozwiązujemy nasze równanie (również podajemy jako macierz):'))

        display(ObjectCode('''odesys = ODESystem(odes=Matrix([eq.lhs-eq.rhs]),dvars=Matrix([x]), ode_order=2)'''))

        display(ReportText('- Rozwiązanie ogólne wcześniej zdefiniowanego równania metodą general solution, czyli rozwiązanie ogólne:'))

        display(ObjectCode('''odesys.general_solution[0]'''))

        display(ReportText('- Rozwiązanie szczególne zdefiniowanego równania metodą steady solution, rozwiązanie szczególne'))

        display(ObjectCode('''odesys.steady_solution'''))

        display(ReportText('- Przykład wykorzystania ODESystem dla oscylatora drgań:'))

        display(ObjectCode(obiekt_2))


        
harmonic_esc_str=(
'''
harmonic_oscilator_eq = Eq(m*x.diff(t,t) -F + c*x.diff(t) + k*x,0)
harmonic_oscilator_eq  
''')

harmonic_oscylator_ode = (
'''
harmonic_oscylator_ode = ODESystem(odes= Matrix([(harmonic_oscilator_eq.lhs/m).expand()]),dvars = Matrix([x]),ode_order=2 )
harmonic_oscylator_ode
''')

harmonic_oscylator_sol = (
'''
harmonic_oscylator_ode = ODESystem(odes= Matrix([(harmonic_oscilator_eq.lhs/m).expand()]),dvars = Matrix([x]),ode_order=2 )
harmonic_oscylator_ode
''')

harmonic_oscylator_sym = (
'''
harmonic_oscylator_sym = harmonic_oscylator_sol.subs(subs_dict).compute_solution(t_span)
harmonic_oscylator_sym
''')

ode_sys_P = (
'''
ode_sys_P = ODESystem.from_dynamic_system(dyn_sys_P.linearized())
ode_sys_P
'''
)

ode_sol_P = (
'''
ode_sol_P = ode_sys_P.solution.with_ics(ics_list)
ode_sol_P
''')

num_mod_P = (
'''
 = dyn_sys_P.linearized().subs(subs_dict_P).numerized().with_ics(ics_list)
num_mod_P
''')

num_sol_P = (
'''
num_sol_P = num_mod_P.compute_solution(t_span)
''')

num_sol_str = (
'''
display(num_sol_P.plot())
''')

num_wyn_SMS = (
'''
num_wyn_SMS = ode_sys_SMS.subs(subs_dict).numerized().with_ics(ics_list).compute_solution(t_span)
''')

ana_wyn_SMS = (
'''
ana_wyn_SMS = ode_sol_SMS.subs(subs_dict).with_ics(ics_list).compute_solution(t_span)
''')

num_wyn_SMS_plot = (
'''
display(num_wyn_SMS[x].plot())
''')

ana_wyn_SMS_plot = (
'''
display(ana_wyn_SMS[x].plot())
''')


ChapterContainer_str=(
'''
from pylatex import Document, Section, Package, Table, Subsubsection
from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, SummaryTable, BaseFrameFormatter, CurrentContainer, ObjectCode, DescriptionsRegistry, SymbolsDescription
dokument=Section(title='Stworzenie raportu  klasy ODESystem')
CurrentContainer(dokument)
''')

Chapter_SympyFormula=(
'''

from sympy import Symbol, Eq
a=Symbol('a')
b=Symbol('b')
nazwa_zdefiniowanego_rownania=Eq(a,b)
display(SympyFormula(nazwa_zdefiniowanego_rownania))
''')

ode_components_path_str=(
'''
MarekSkwarka/Shared_files/modules/dynpy/utilities/components/ode/en.py
''')



property_report_formula=(
"""
from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase, EngeneeringDrawingGuide,DevelopmentGuide
from pylatex.utils import italic, NoEscape

#doc_final = DevelopmentGuide('./Output/[PL][Guide][Programming][PL][Guide][ODESystem]Obsluga klasy ODESystem',documentclass=NoEscape('article'),document_options=['a4paper','eqn'],lmodern=False)

doc_final = DevelopmentGuide('./output/Obsluga klasy ODESystem',documentclass=NoEscape('article'),document_options=['a4paper','eqn'],lmodern=False)
doc_final.packages.append(Package('natbib',options=['numbers']))
doc_final.packages.append(Package('booktabs'))
doc_final.packages.append(Package('float'))
doc_final.packages.append(Package('siunitx'))

intro = Section(title='Intro - nazwa')
CurrentContainer(intro)
display(Markdown(f'''1
                tekst tekst'''))
dokument = Section(title='Dokument- nazwa')
CurrentContainer(dokument)
display(Markdown(f'''2
                tekst tekst'''))
komponent = Section(title='Komponent - nazwa')
CurrentContainer(komponent)
display(Markdown(f'''3
                tekst tekst'''))
symki = Section(title='Symki - nazwa')
CurrentContainer(symki)
display(Markdown(f'''4
                tekst tekst'''))

doc_final.append(intro)
doc_final.append(dokument)
doc_final.append(komponent)
doc_final.append(symki)

doc_final.generate_pdf()
result = False
result = doc_final.generate_pdf()
""")

class ODEReportComponent(ReportComponent):
    
    title="Stworzenie raportu"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Aby utworzyć dokument, standardowo należy utworzyć sekcje lub rozdziały w następujący sposób:'))
#        display(Picture('./dynpy/utilities/components/guides/images/dok.jpg'))
        display(ObjectCode(ChapterContainer_str))
        display(ReportText('Jeżeli chcemy żeby nasze równania ładnie wyswietlały się w dokumencie nie możemy zapominac o użyciu SympyFormula:'))
        display(ObjectCode(Chapter_SympyFormula))
#        display(Picture('./dynpy/utilities/components/guides/images/sympaj.jpg'))
        display(ReportText('Co więcej, ODESystem posiada ogromną ilosć komponentów raportujących, które tworzą sekcje raportu za nas. Scieżka do ich odnalezienia jest następująca:'))
        display(ObjectCode(ode_components_path_str)) 
#        display(Picture('./dynpy/utilities/components/guides/images/sciezka_ode.jpg'))
        display(ReportText('Co więcej, istnieje też property .report, o którym informacja znajduje się rozdział niżej. Aby wygenerować plik pdf z użyciem wybranego szablonu należy zdefiniować dokuemnt jako wybrany template. Przypomnienie: Sprawdź, czy wszystkie katalogi w ścieżce poniżej istnieją. Jeśli nie, utwórz brakujące foldery przed próbą wygenerowania pliku PDF. Dobrym nawykiem jest stworzenie katalogu na pliki tekstowe wewnątrz swojego projektu. Następnie dodaj wszystkie sekcje używając metody append'))
        display(ObjectCode(property_report_formula)) 
#        display(Picture('./dynpy/utilities/components/guides/images/appeend.jpg'))

sym_num_system_code_str = ('''

m,F,c,k,x=symbols('m,F,c,k,x', positive=True)
t=Symbol('t')
x=Function(x)(t)

subs_dict={m:100 ,
          F:200 ,
          c:50 ,
          k:500,
          }

ics_list = [0.0,0.0]

''')



sym_harmonic_oscilator_eq = ('''

harmonic_oscilator_eq = Eq(m*x.diff(t,t) -F + c*x.diff(t) + k*x,0)
harmonic_oscilator_eq

''')

sym_harmonic_oscylator_ode = ('''

harmonic_oscylator_ode = ODESystem(odes= Matrix([(harmonic_oscilator_eq.lhs/m).expand()]),dvars = Matrix([x]),ode_order=2 )
harmonic_oscylator_ode

''')


sym_harmonic_oscylator_sol = ('''

harmonic_oscylator_sol = harmonic_oscylator_ode.numerized(subs_dict)
display(harmonic_oscylator_sol)
harmonic_oscylator_sym = harmonic_oscylator_sol.compute_solution(t_span=t_span, ic_list=ics_list)
display(harmonic_oscylator_sym)
harmonic_oscylator_sym.plot()

''')




sym_Example_PulledPendulum =(
'''

from dynpy.models.mechanics.pendulum import  PulledPendulum
t_span=np.linspace(0, 2, 1000)

dyn_sys_P = PulledPendulum()
display(dyn_sys_P._eoms)

subs_dict_P = dyn_sys_P.get_numerical_parameters()
display(subs_dict_P)
subs_dict_P = { **subs_dict_P, dyn_sys_P.F: 10, dyn_sys_P.Omega: 3.14/5}
for key, item in subs_dict_P.items():
    print(key)
    if key == dyn_sys_P.l**2*dyn_sys_P.m:
        key_to_pop  = key
subs_dict_P.pop(key_to_pop)
display(subs_dict_P)
ics_list_P = [0.0,0.0] # 1 stopien swobody

''')

sym_ode_sol_P_plot = ('''

ode_sol_P = dyn_sys_P._ode_system.numerized(subs_dict_P)

display(ode_sol_P)

sym_wyn_P = ode_sol_P.compute_solution(t_span=t_span, ic_list=ics_list_P)
sym_wyn_P.plot()

''')

sym_ode_sys_P =( '''

ode_sys_P = ODESystem.from_dynamic_system(dyn_sys_P.linearized())
ode_sys_P


''')

sym_ode_sol_P =( '''

ode_sol_P = ode_sys_P.solution.with_ics(ics_list)
ode_sol_P

''')

sym_num_sol_P=(
'''
num_sol_P = num_mod_P.compute_solution(t_span=t_span, ic_list=ics_list_P)

display(num_sol_P.plot())
''')

sym_coords=(
'''
display(harmonic_oscylator_ode)

display(harmonic_oscylator_sol)

coords= ode_sol_P.lhs
display(coords)
display(ode_sys_P)
display(subs_dict_P)
nadf = NumericalAnalysisDataFrame.from_model(ode_sys_P,F,[10, 20, 30], reference_data=subs_dict_P,coordinates=list(coords), index=pd.Index(t_span, name=t),ics=None)

nadf_sol = nadf.perform_simulations()

nadf_sol.columns = nadf_sol.columns.droplevel(0)
nadf_sol

''')

sym_wywolanie_helpa=(
'''
help(NumericalAnalysisDataFrame)
''')

sym_nadf=(
'''
nadf_sol.plot()
''')

sym_sms=(
"""
Omega = Symbol('Omega', positive=True)

subs_dict={m:10,
    F:100,
    c:0.01,
    k:500,
    }

t_span = np.linspace(1,2000,2001)

ode_sys_SMS = harmonic_oscylator_ode.subs({F: F*sin(Omega*t)}).subs({Omega: sqrt(k/m)+0.01})

ode_sys_SMS
""")

sym_num_wyn_SMSe=(
'''
num_wyn_SMS = ode_sys_SMS.subs(subs_dict).numerized().with_ics(ics_list).compute_solution(t_span)
''')

sym_ana_wyn_SMS=(
'''
ana_wyn_SMS = ode_sol_SMS.subs(subs_dict).with_ics(ics_list).compute_solution(t_span)
''')

sym_sms_plot=(
'''
display(num_wyn_SMS[x].plot())
''')

sym_comparison=(
"""
sym_comparison = pd.DataFrame(num_wyn_SMS[x])
sym_comparison[Symbol('x(t)_a')] = ana_wyn_SMS[x]

sym_comparison = TimeDataFrame(sym_comparison)
""")

sym_comparison_plot=(
'''
sym_comparison.to_frequency_domain().single_sided_rms().truncate(0,0.2).plot()
''')



# sym_wyn_code=(
# '''
# sym_wyn_P = ode_sol_P.subs(subs_dict_P).compute_solution(t_span)
# sym_wyn_P.plot()
# ''')



class ODENumericalSimulationsComponent(ReportComponent):
    
    title="Symulacje numeryczne- różne techniki symulacji"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('Przeprowadzenie symulacji przy pomocy ODE'))
        display(ReportText('W pierwszej kolejności zdefiniowano parametry potrzebne do przeprowadzenia do obliczeń'))

        display(ObjectCode(sym_num_system_code_str))

        display(ReportText('Wyprowadzenie wzoru na oscylator harmoniczny'))
        display(ObjectCode(sym_harmonic_oscilator_eq))
        #display(ObjectCode(harmonic_esc_str))

        display(ReportText('Przeprowadzenie symulacji w ODE'))


        display(ObjectCode(sym_harmonic_oscylator_ode))


        display(ObjectCode(sym_harmonic_oscylator_sol))

        #display(ObjectCode)

        #display(ObjectCode())

        display(ReportText('Wprowadzenie systemu dynamicznego, zdefiniowano również pojawiającą się Omegę '))
        display(ObjectCode(sym_Example_PulledPendulum))

        display(ObjectCode(sym_ode_sol_P_plot))
        display(ObjectCode(sym_ode_sys_P))
        display(ObjectCode(sym_ode_sol_P))

        display(ReportText('Symulacje używając systemu dynamicznego zdefiniowanego przy użyciu modułu dynpy'))

        display(ReportText('Sprawdzenie za pomocą helpa czym tak na prawdę jest NumericalAnalysisDataFrame'))
        display(ObjectCode(sym_wywolanie_helpa))

        display(ReportText('Wyświetlenie poszczególnych wykresów, wyników oraz tablicy'))
        display(ObjectCode(sym_coords))
        display(ObjectCode(sym_nadf))
        display(ReportText('Analiza porównawcza symulacji rozwiązania analitycznego i numerycznego'))
        display(ObjectCode(sym_sms))
        display(ObjectCode(sym_num_wyn_SMSe))
        display(ObjectCode(sym_ana_wyn_SMS))
        display(ObjectCode(sym_sms_plot))
        display(ObjectCode(sym_comparison))

#         display(ObjectCode('''
#         display(harmonic_oscylator_ode)

#         display(harmonic_oscylator_sol)

#         coords= harmonic_oscylator_sol.lhs

#         nadf = NumericalAnalysisDataFrame.from_model(harmonic_oscylator_sol,F,[10, 20, 30], reference_data=subs_dict,coordinates=list(coords), index=pd.Index(t_span, name=t),ics=None)

#         nadf_sol = nadf.perform_simulations()

#         nadf_sol.columns = nadf_sol.columns.droplevel(0)

#         nadf_sol
#         '''))

#         display(ObjectCode('''nadf_sol.plot()'''))




#         display(ObjectCode('''
#         Omega = Symbol('Omega', positive=True)

#         subs_dict={m:10,
#                    F:100,
#                    c:0.01,
#                    k:500,
#                    }

#         t_span = np.linspace(1,2000,2001)

#         ode_sys_SMS = harmonic_oscylator_ode.subs({F: F*sin(Omega*t)}).subs({Omega: sqrt(k/m)+0.01})

#         ode_sys_SMS
#         '''))


#         display(ObjectCode('''
#         ode_sol_SMS = ode_sys_SMS.solution
#         ode_sol_SMS.with_ics(ics_list)

#         '''))

#         display(ObjectCode(num_wyn_SMS))


#         display(ObjectCode(ana_wyn_SMS))
#         display(ObjectCode(num_wyn_SMS_plot))
#         display(ObjectCode(ana_wyn_SMS_plot))

#         display(ObjectCode('''
#         sym_comparison = pd.DataFrame(num_wyn_SMS[x])
#         sym_comparison[Symbol('x(t)_a')] = ana_wyn_SMS[x]

#         sym_comparison = TimeDataFrame(sym_comparison)
#         '''))

#         display(ObjectCode('''
#         sym_comparison.to_frequency_domain().single_sided_rms().truncate(0,0.2).plot()
#         '''))

definicja_danych=(
'''
m=Symbol('m',positive=True)
g=Symbol('g')
alpha=Symbol('alpha')
v_0=Symbol('v_0')
t=Symbol("t")
x=Function("x")(t)
y=Function("y")(t)

dane = {m: 10, g: 9.81, alpha: pi/6, v_0:20}
''')

rownania_predkosci=(
'''
v_x=Eq(x.diff(t),v_0*cos(alpha))
v_y=Eq(y.diff(t),v_0*sin(alpha)-g*t)
display(v_x, v_y)
''')
definicja_ode=(
'''
ODE_x=ODESystem(odes=Matrix([v_x.lhs-v_x.rhs]),dvars=Matrix([x]),ode_order=1)
ODE_y=ODESystem(odes=Matrix([v_y.lhs-v_y.rhs]),dvars=Matrix([y]),ode_order=1)
display(ODE_x, ODE_y)
''')
ode_solution_code=(
'''
ODE_x_sol=ODE_x.subs(dane).solution.with_ics([20])
ODE_x_sol
''')
class ProjectileExampleComponent(ReportComponent):
    
    title="Przykład użycia ODESystem na podstawie rzutu ukośnego"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        from dynpy.solvers.linear import ODESystem
        m, g, v, alpha = symbols('m g v_0 alpha')

        t = Symbol('t')
        x = Function('x')(t)
        y = Function('y')(t)

        dane = {m: 10, g: 9.81, alpha: pi/6, v:20}
        predkosc_x=Eq(x.diff(t),v*cos(alpha))
        predkosc_y=Eq(y.diff(t),v*sin(alpha)-g*t)

        ode_x_x=ODESystem(odes=Matrix([predkosc_x.lhs-predkosc_x.rhs]),dvars=Matrix([x]),ode_order=1)
        ode_x_y=ODESystem(odes=Matrix([predkosc_y.lhs-predkosc_y.rhs]),dvars=Matrix([y]),ode_order=1)

        ODE_x_sol=ode_x_x.subs(dane).solution.with_ics([20])
        ODE_x_sol
        
        display(ReportText('Implementacja biblioteki ODESystem'))
        display(ObjectCode('''from sympy import *
from dynpy.solvers.linear import ODESystem
from sympy.physics.mechanics import dynamicsymbols'''))
  
        display(ReportText('Definiowanie zmiennych i utworzenie słownika'))
        display(ObjectCode(definicja_danych))
        
        display(ReportText('Równania ruchu dla osi x oraz y'))
        display(ObjectCode(rownania_predkosci))
        display(ReportText('Zapisanie równań za pomocą ODESystem'))
        display(ObjectCode(definicja_ode))

        display(ReportText('Wywołanie odpowiednio zmiennych i funkcji daje następujące rezultaty'))
        display(ObjectCode('''ODE_x'''))
        display(SympyFormula(ode_x_x))
        display(ObjectCode('''ODE_x.general_solution[0]'''))
        display(SympyFormula(ode_x_x.general_solution[0]))
        display(ObjectCode('''ODE_x.steady_solution[0]'''))
        display(SympyFormula(ode_x_x.steady_solution[0]))
        display(ObjectCode('''ODE_x.solution'''))
        display(SympyFormula(ode_x_x.solution[0]))
        display(ObjectCode('''ODE_x.solution.subs(dane)'''))

        display(SympyFormula(ode_x_x.solution.subs(dane)))#Fsimulations
        
        display(ObjectCode(ode_solution_code))

        display(SympyFormula(ODE_x_sol))
        display(ReportText('Ten sam wynik można osiągnąć stosując metoodę " dsolve ", jednak nie jest ona zalecana do równań różniczkowych'))
        display(ObjectCode('''dsolve(v_y.subs(dane))'''))
        display(SympyFormula(dsolve(predkosc_y.subs(dane))))

        
ode_import=(
'''
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from dynpy.solvers.linear import ODESystem
from sympy import Symbol, symbols, Function, Matrix
from dynpy.utilities.adaptable import *
''')

ode_variables=(
"""
t = Symbol('t')

U_z, R_w, L_w, E, U_Rw, U_Lw, k_e = symbols('U_z R_w L_w E U_Rw U_Lw k_e', positive=True)
M_s, B, J, M_obc, k_m = symbols('M_s B J M_obc k_m', positive=True)

i_w = Function('i_w')(t)
omega_s = Function('omega_s')(t)

dane = {
    U_z: 10,
    R_w: 2,
    L_w: 0.1,
    k_e: 0.1,
    k_m: 0.1,
    J: 0.1,
    B: 0.5,
    M_obc: 0
}

t_span = np.linspace(0, 1, 100)
""")

ode_objrctM=(
'''
ode = ODESystem(
    odes=Matrix([
        U_z - R_w * i_w - L_w * i_w.diff(t) - k_e * omega_s,
        k_m * i_w - J * omega_s.diff(t) - B * omega_s - M_obc
    ]),
    dvars=Matrix([i_w, omega_s]),
    ode_order=1
)
ode
''')

ode_met_1A=(
'''
ode.subs(dane).solution.with_ics([0,0])
''')

ode_tab_1A=(
'''
sym_an = ode.solution.with_ics([0,0]).subs(dane).numerized().compute_solution(t_span=t_span)
sym_an[[i_w, omega_s]]
''')

ode_met_1B=(
'''
sym_num = ode.subs(dane).numerized().compute_solution(t_span=t_span, ic_list=[0, 0])
sym_num[[i_w, omega_s]]
''')

ode_met_2=(
'''
X = [i_w, omega_s]

num_df = NumericalAnalysisDataFrame.from_model(
    ode.solution.with_ics([0, 0]),
    U_z,
    [10, 20, 50],
    reference_data=dane,
    coordinates=X,
    index=pd.Index(t_span, name=t)
)
num_df
''')

ode_met_PS=(
'''
sym_wyn = num_df.perform_simulations().droplevel(0,axis=1)

display(sym_wyn)
''')

ode_import_example=(
'''
from dynpy.models.mechanics.trolley import SpringDamperMassSystem

dyn_sys = SpringDamperMassSystem()

dyn_sys._ode_system
''')

ode_example2=(
'''
dane = dyn_sys.get_numerical_parameters()

data_dict = {**dane, dyn_sys.F: 1000}

data_dict
''')

ode_simulation=(
'''
simulation = dyn_sys._ode_system.subs(data_dict).numerized().compute_solution(
    np.linspace(0, 2, 1000), ic_list=[10, -5]
)

simulation
''')

ode_simulation_plot=(
'''
simulation.plot()
plt.grid()
plt.show()
''')




        
class ODESimulationComponent(ReportComponent):
    
    title="Przeprowadzanie symulacji"
    


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('To od czego standardowo powinno się zacząć to importy oraz definicja zmiennych. Potrzebne importy:'))
        display(ObjectCode(ode_import))
        display(ReportText('Definiowanie zmiennych:'))
        display(ObjectCode(ode_variables))
        display(ReportText('Kolejny niezbędny krok to tworzenie obiektu ODE. Należy go wykonać w następujący sposób:'))
        display(ObjectCode(ode_objrctM))
        display(ReportText('METODA 1A: Przeprowadzenie symulacji na gotowm rozwiazaniu analitycznym poprzez użycie property solution. Samo wywołanie                                 odpowiedniego rozwiązania analitycznego pokazuje nam poniższy kod:'))
        display(ObjectCode(ode_met_1A))
        display(ReportText('Następnie generujemy tabele niezbędną do symulacji w taki sposób:'))
        display(ObjectCode(ode_tab_1A))
        display(ReportText('METODA 1B: Innym sposobem jest przeprowadzenie symulacji numerycznych bez obliczeń analitycznych. Potrzebny jest do tego taki kod: '))
        display(ObjectCode(ode_met_1B))
        display(ReportText('METODA 2: Trzecim sposobem generowania symulacji jest użycie NumericalAnalysisDataFrame, kod wywołania prezentuje się następująco:'))
        display(ObjectCode(ode_met_2))
        display(ReportText('Aby przeprowadzić symulacje konieczne jest zastosowanie metody perform simulations. Bez zastosowania tej metody dostaniemy następujący wynik:'))
        display(ReportText('A z jej użyciem, wszystko powinno zadziałać i dostaniemy następujący output:'))
        display(ObjectCode(ode_met_PS))
        display(ReportText('Dla lepszego zrozumienia zagdnień poniżej przedstawione zostało zastosowanie symulacji na konkretnym przykładzie. Skorzystamy z                           prostej klasy SpringDamperMassSystem. Na dobry początek należy zaimportować tę klasę oraz przypisać ją do zmiennej. W przeprowadzeniu                                     symulacji pomoże nam wykorzystanie metody $ode system$ w celu pozyskania równania ruchu naszego systemu dynamicznego. Aby pozyskać dane w                                 prosty sposób korzystamy z metody $get numerical parameters$ oraz wrzucamy te dane w słownik. Kod do wykonania tych czynnosci prezentuje się                             następująco: '))
        display(ObjectCode(ode_import_example))
        display(ReportText('Następnie stosując jedną z wyżej pokazanych metod, możemy stworzyć tabelę:'))
        display(ObjectCode(ode_example2))
        display(ReportText('Ostatnim krokiem jest wizualizacja danych, czyli stworzenie wykresu. Aby taki wykres uzyskać wystarczy wpisać kod widoczny na                             obrazku:'))
        display(ObjectCode(ode_simulation))
        display(ObjectCode(ode_simulation_plot))

####
omega, omega2 = symbols('omega Omega')
t=Symbol('t')
x=Function('x')(t)
y=Function('y')(t)
ode=(
'''
omega, omega2 = symbols('omega Omega')
ODESystem(odes=Matrix([omega**2*x-sin(omega2*t)+x.diff(t,t)]),dvars=Matrix([x]),ode_order=1)
''')   
        
        
class ReportCompUseComponent(ReportComponent):
    
    title="Podstawy używania komponentu raportującego"


    def append_elements(self): 
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        komponent = Chapter('Podstawy używania komponentu raportującego')
        CurrentContainer(komponent)        
        display(ReportText('W celu użycia komponentu raportującego należy użyć $.report$. Przykładowo wpisując *odesys.report* otrzyma się następujący output:'))
        display(ReportText('The investigated system is described by differential equations being as follows:'))
        display(ObjectCode(ode))
        display(ReportText('To solve the problem serve several methods depending on the equation type. As dynamic systems s behaviour is described by an ordinary differential equation, the variables of the equation are as follows: t [x]'))
        display(ReportText('The variables allow to study and analyse the system s dynamics varying with time.'))
        
class NumericalAnalysisSimulationComponent(pl.NumericalAnalysisSimulationComponent):
    
    title="Numerical analysis simulation Component"

class AnalyticalSimulationComponent(pl.AnalyticalSimulationComponent):
    
    title="Analytical simulation Component"
    
class DynamicSystemCompletenessCheckComponent(pl.DynamicSystemCompletenessCheckComponent):
    
    title="Dynamic system completeness check component"
    
class IssuePreparationComponent(pl.IssuePreparationComponent):
    
    title="Issue preparation component"

    
    
class_call_str=(
'''
from dynpy.models.mechanics.principles import ComposedSystem
from sympy import *

class MyMaterialPointMovement(ComposedSystem):
    pass
''')

atributes_str=(
'''
m = Symbol('m', positive=True)
g = Symbol('g', positive=True)
c = Symbol('c', positive=True)
r = Symbol('r', positive=True)
phi = dynamicsymbols('phi')

c0 = Symbol('c0', positive=True)
r0 = Symbol('r0', positive=True)
phi0 = dynamicsymbols('phi0')
''')

innit_str=(
'''
def __init__(self,
             m=None,
             g=None,
cLatexDataFrame.set_default_units(unit_dict)

tabelka=LatexDataFrame.formatted(
    data=dane,
    index=['Value']).applymap(lambda x: f'${latex(x)}$').rename_axis(
        'Symbol',
        axis=1)

display(tabelka.reported(caption='Data given in the exercise'))

''')

component_str=(
'''
@property
def components(self):

    components = {}


    self._mass_x = MaterialPoint(self.m,
                                 pos1=self.r * sin(self.phi),
                                 qs=self.qs)
    self._mass_y = MaterialPoint(self.m,
                                 pos1=self.r * cos(self.phi),
                                 qs=self.qs)

    self._gravity_ = GravitationalForce(self.m,
                                        self.g,
                                        pos1=self.r * cos(self.phi),
                                        qs=self.qs)



    components['_mass_x']=self._mass_x
    components['_mass_y']=self._mass_y
    components['_gravity_']=self._gravity_


    return components
''')


symb_desc_str=(
'''
def symbols_description(self):
    self.sym_desc_dict = {
        self.m: r'Mass',
        self.g: r'Gravity constant',
        self.c: r'',
    }

    return self.sym_desc_dict
''')

def_data_str=(
'''
def get_default_data(self):

    m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

    default_data_dict = {
        self.m: [m0 * no for no in range(1, 8)],
        self.c: [c0 * no for no in range(1, 8)],
        self.r: [r0 * no for no in range(1, 8)],
        self.phi: [phi0 * no for no in range(1, 8)],
    }

    return default_data_dict
''')

num_data_str=(
'''
def get_numerical_data(self):

    m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

    default_data_dict = {
        self.m: [m0 * no for no in range(1, 8)],
        self.c: [c0 * no for no in range(1, 8)],
        self.r: [r0 * no for no in range(1, 8)],
        self.phi: [phi0 * no for no in range(1, 8)],
    }

    return default_data_dict
''')

dict_str=(
'''
def unit_dict(self):

    from pint import UnitRegistry
    ureg=UnitRegistry()

    unit_dict = {
        self.m: ureg.kilogram,
        self.g: ureg.meter/ureg.second/ureg.second,
        self.c: ureg.kilogram/ureg.second,
        self.r: ureg.meter,
        self.phi: ureg.radian,
        self.c0: ureg.kilogram/ureg.second,
        self.r0: ureg.meter,
        self.phi0:  ureg.radian
    }
    return unit_dict
''')

full_class_str=(
'''
from dynpy.models.mechanics.principles import ComposedSystem
from sympy import *

class MyMaterialPointMovement(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('phi')

    c0 = Symbol('c0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        self.ivar = ivar

        self.qs = [self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}


        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity_ = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)


      
        components['_mass_x']=self._mass_x
        components['_mass_y']=self._mass_y
        components['_gravity_']=self._gravity_
     

        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict

    def get_numerical_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict
        
    def unit_dict(self):

        from pint import UnitRegistry
        ureg=UnitRegistry()

        unit_dict = {
            self.m: ureg.kilogram,
            self.g: ureg.meter/ureg.second/ureg.second,
            self.c: ureg.kilogram/ureg.second,
            self.r: ureg.meter,
            self.phi: ureg.radian,
            self.c0: ureg.kilogram/ureg.second,
            self.r0: ureg.meter,
            self.phi0:  ureg.radian
        }
        return unit_dict        
''')

class DynSysImplementationComponent(ReportComponent):
    
    title="Implementacja systemow dynamicznych"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(ReportText("Tworzenie klasy należy zacząć od wybrania nazwy opisującej nasz system oraz zadeklarowania jej w odpowiedni sposób, który został przedstawiony poniżej:")) #itp
        display(GuideCode(class_call_str))
        display(ReportText('Następnie deklarujemy atrybuty/symbole, które są niezbędne do opisania danego obiektu'))
        display(GuideCode(atributes_str))
        display(ReportText('Do zdefiniowania inicjalizacji obiektu klasy z określonymi atrybutami służy poniższa struktura.'))
        display(GuideCode(innit_str))
        display(ReportText('W tym kroku korzystamy z predefiniowanych elementów takich jak punkty materialne, dyski, itp.'))
        display(GuideCode(component_str))
        display(ReportText('W tej części opisujemy parametry klasy, np. masę'))
        display(GuideCode(symb_desc_str))
        display(ReportText('Sekcja ta przedstawia inicjalizację słownika zawierającego dane symboliczne do równań analitycznych systemu:'))
        display(GuideCode(def_data_str))
        display(ReportText('Kolejnym krokiem jest przypisanie słowników zawierających dane do symulacji numerycznych:'))
        display(GuideCode(num_data_str))
        display(ReportText('Ta sekcja opisuje definiowanie jednostek używanych przez system dynamiczny:'))
        display(GuideCode(dict_str))
        
        display(ReportText('Ostatecznie kod implementacji całego systemu ma następującą postać:'))
        display(GuideCode(full_class_str))
        
## TABELKA InterimSchedule 

spotkanie_1='08.10.24'
spotkanie_2='15.10.24'
spotkanie_3='22.10.24'
spotkanie_4='29.10.24'
spotkanie_5='05.11.24'
spotkanie_6='12.11.24'
spotkanie_7='19.11.24'
spotkanie_8='26.11.24'
spotkanie_9='03.12.24'
spotkanie_10='10.12.24'


zadanie_1='Ustalenie tematu pracy przejściowej i wstępnego toku postępowania'
zadanie_2='Wyprowadzenie równań opisujących zjawisko'
zadanie_3='Zaimplementowanie równań w środowisku obliczeniowym'
zadanie_4='Dodanie opisu oraz listy symboli występujących w równaniach'
zadanie_5='Przeprowadzenie badań/symulacji, uzyskanie danych pomiarowych)'
zadanie_6='Wykonanie tabeli, przetworzenie jej na wykresy, opis zachodzących zależności'
zadanie_7='Opracowanie wniosków, wyników odnośnie przeprowadzonych badań/symulacji'
zadanie_8='Formowanie wstępu  oraz podsumowania pracy '
zadanie_9='Poprawa błędów'
zadanie_10='Poprawa błędów, oddanie pracy'

wl='Praca zasadnicza'
po='Poprawki'
        
class InterimScheduleComponent(ReportComponent):
    
    title="Implementacja systemow dynamicznych"


    @property
    def reported_object(self):

        default_data = {'date':datetime.datetime(2024,7,16),
                       'title':'Title of interim project',
                       'timedelta':datetime.timedelta(days=7),
                       'issue_no':359,
                       }

    
        
        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

#         elif isinstance(self._reported_object, str):
#             return {**default_data,'date':self._reported_object}
        elif isinstance(self._reported_object, str):
            return {**default_data, 'date': datetime.datetime.strptime(self._reported_object, "%Y-%m-%d")}
                    
        elif self._reported_object is None:
            return default_data

        else:
#             return self._reported_object
            return {**default_data, 'date': self._reported_object}


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj


    def append_elements(self):
        
        first_meeting_date = self.reported_object['date'] # it's useless in the case of permanent content - it's commented for future usage
        days_inc = self.reported_object['timedelta']


        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(ReportText(f"Pracę rozpoczynamy {first_meeting_date}")) #itp
        
        tabela_git = {
        'Tydzień prac': [1,2,3,4,5,6,7,8,9,10],
        'Data': [first_meeting_date+datetime.timedelta(days=days_inc.days*no) for no in  range(10)],
        'Temat konsultacji': [zadanie_1, zadanie_2, zadanie_3, zadanie_4, zadanie_5, zadanie_6, zadanie_7, zadanie_8, zadanie_9, zadanie_10],
        'Charakter spotkania': [wl,wl,wl,wl,wl,wl,wl,wl,po,po]
            
        
}

        tabelka_pandas=pd.DataFrame(tabela_git)
        tabelka_pandas.index += 1
        tabelka_latex=LatexDataFrame(tabelka_pandas).reported(caption='Harmonogram spotkań')
        display((display(tabelka_latex)))

class InterimTemplateComponent(InterimScheduleComponent):
    
    title="Implementacja systemow dynamicznych"


    
    
    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(ReportText("Tu kod dla podstawowej formy dokumentu")) #itp
    
    
    
project_issue_title_str = """
Maintenance of `{system_name}` class which is dynamic system representation        
"""
        
project_issue_desc_str = """

The following problems have to be checked or fixed, in order to ensure the correctness of implemented system:

- [ ] checking if a class is executable,

- [ ] validation of figures for schemes and real examples,

- [ ] validation of reference parameters,

- [ ] validation of random parameters,

- [ ] validation of  description of  description of parameters,

- [ ] validation of units.
"""
        
class InterimIssuesComponent(InterimScheduleComponent):
    
    title="Wzór na tworzenie issues w systemie GitHub dla przejsciowe"

    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText("Przykładowy tekst, który pozwoli na sprawne przygotowywanie issue ma następującą formę:"))

        display(ReportText("# Title: "))
        display(ObjectCode(project_issue_title_str.format(system_name=system_name)))

        display(ReportText("# Description: "))
        display(ObjectCode(project_issue_desc_str))
    
    
    
reportcomp_issue_title_str = """
Implementation of `{classname}` class that creates a part of a report
"""
        
reportcomp_issue_desc_str = """
Class is to implement in `{module}` module that creates a part (reporting component - formally child of `ReportComponent` class) of {field} for {target}.
"""
    
class ReportCompImplementationIssueComponent(ReportComponent):
    
    title="Issue na implementację komponentów raportujących"



    def append_elements(self):
        
        classname = self.reported_object['classname'] # it's useless in the case of permanent content - it's commented for future usage
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']


        display(ReportText("Przykładowy tekst, który pozwoli na sprawne przygotowywanie issue ma następującą formę:"))

        display(ReportText("# Title: "))
        display(ObjectCode(reportcomp_issue_title_str.format(classname=classname)))

        display(ReportText("# Description: "))
        display(ObjectCode(reportcomp_issue_desc_str.format(classname=classname,module=class_module,field=class_field,target=target)))
        
        
rep_comp_call = """
from dynpy.utilities.components.guides.en import ReportComponent
from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, GuideCode

class {classname}(ReportComponent):
"""

rep_comp_call_with_pass = rep_comp_call + '\n \t pass'


title_call_str="""
    title="Implementacja komponentów raportujących"
"""

append_elem_str ="""

    @property
    def reported_object(self):


        default_data = {'classname':'ReportingModuleIntroComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':359,
                       }


        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data

        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
"""

rep_obj_str ="""
    #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
"""




class ReportCompImplementationComponent(ReportComponent):
    
    title="Implementacja komponentów raportujących"


    @property
    def reported_object(self):

        
        default_data = {'classname':'ReportingModuleIntroComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':359,
                       }

        
        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}
        
        elif self._reported_object is None or not isinstance(self._reported_object, dict):
            return default_data

        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj
    
    
    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(ReportText("Tworzenie komponentu należy zacząć od wybrania nazwy opisującej nasz komponent raportujący oraz zadeklarowania go w odpowiedni sposób, który został przedstawiony poniżej:")) #itp
        display(GuideCode(rep_comp_call_with_pass.format(classname=classname)))
        display(ReportText('Następnie deklarujemy nazwę komponentu, która będzie wyświetlana jako nagłówek.'))
        display(GuideCode(title_call_str))
        
        display(ReportText('Kolejnym krokiem tworzenia komponentu, jest zadeklarowanie metody `append_elements`, która zawiera zasadniczą część komponentu.'))
        display(GuideCode(append_elem_str))
        display(ReportText('Dodatkowo na początku tej metody należy umieścić następujący kod, który w razie potrzeby pozwala na to aby komponent obsługiwał dowolny system dynamiczny.'))
        display(GuideCode(rep_obj_str))
        
        display(ReportText('Ostatecznie kod implementacji całego komponentu ma następującą postać:'))
        display(GuideCode(rep_comp_call.format(classname=classname)   + title_call_str + append_elem_str  ))
        
        
######### GUIDE REPORTING COMPONENTS
######### GUIDE REPORTING COMPONENTS
######### GUIDE REPORTING COMPONENTS

doc_gen_str =(
'''

doc_final = MechanicalCase('./output/nazwa_dokumentu',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)
doc_final.packages.append(Package('natbib', options=['numbers']))
doc_final.packages.append(Package('booktabs'))
doc_final.packages.append(Package('float'))
doc_final.packages.append(Package('siunitx'))


doc_final.append(sekcja)
doc_final.append(sekcja1)
doc_final.append(sekcja2)
doc_final.append(sekcja3)
doc_final.append(sekcja4)
doc_final.append(sekcja5)

doc_final.generate_pdf()

''')

class DocumentGenerationComponent(ReportComponent):

    #title="Generowanie dokumentu"
    title="Document generation"


    @property
    def reported_object(self):


        default_data = {'classname':'ReportingModuleIntroComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':359,
                       }


        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data

        elif self._reported_object is None or not isinstance(self._reported_object, dict):
            return default_data
        
        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        
       
        #implement reporting activieties here
        
        #display(ReportText('Ostatni krok to zaapendowanie sekcji i utworzenie dokumentu pdf :'))
        display(ReportText('The last step is to append a section and create a pdf document :'))
        display(GuideCode(f'{doc_gen_str}'))


unit_registry_str = """        


t=Symbol('t')
r=Symbol('r', positive=True)
a=Symbol('a', positive=True)
b=Symbol('b', positive=True)
h=Symbol('h', positive=True)
vphi = Function('varphi')(t)

dane={
    vphi.diff():25,
    r:0.15,
    a:0.2,
    b:0.6,
    h:0.4,
}

unit_dict = {
            t:ureg.second,
            r: ureg.meter,
            h: ureg.meter,
            a: ureg.meter,
            b:ureg.meter,
            vphi.diff(t):ureg.rad/ureg.sec,
   }

LatexDataFrame.set_default_units(unit_dict)

tabelka=LatexDataFrame.formatted(
data=dane,
index=['Value']).applymap(lambda x: f'${latex(x)}$').rename_axis(
    'Symbol',
    axis=1)

display(tabelka.reported(caption='Data given in the exercise'))

"""

class UnitRegistryIntroComponent(DocumentGenerationComponent):

    title="Rejestr jednostek"



    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        display(ReportText('To display the units of symbols used in the report in tables and charts, you must define a register of units:'))


        display(GuideCode(f'{unit_registry_str}'))

container_code='''
nazwa_sekcji=Section('Tytuł')
CurrentContainer(nazwa_sekcji)
'''
class CurrentContainerComponent(DocumentGenerationComponent):

    title="Klasa CurrentContainer i jej zastosowanie"



    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        display(ReportText('Wszystko co chcemy, aby znalazło się w naszym raporcie musi być wrzucone do sekcji, które potem zostaną dodane (*zaapendowane*) do naszego dokumentu. Aby utworzyć sekcję należy wywołać kod:'))
        display(ObjectCode(container_code))
        display(ReportText('Za pomocą **CurrentContainer** wszystko zostanie zebrane do naszej sekcji do momentu utworzenia nowej. (Zaleca się podział 1 sekcja - 1 cela w jupiterze).'))
        
        from dynpy.models.mechanics.engine import Engine
        
from dynpy.utilities.components.mech.en import KineticEnergyComponent, KineticEnergySymPyCodeComponent


SympyFormulaComponent_str='''
from sympy import *

M = Symbol('M')

t = Symbol('t')

m_e = Symbol('m_e')

e = Symbol('e')

z = Function('z')

varphi = Function('varphi')

Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2

display(Eq(Symbol('T'),Ek))
'''


class SympyFormulaComponent(DocumentGenerationComponent):

    title = "Klasa SympyFormula i jej zastosowanie"



    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        #implement reporting activieties here
        display(ReportText('Przykładowe zastosowanie komponentu do obliczeń symbolicznych, na przykładzie energii kinetycznej:'))
        display(GuideCode(f'{SympyFormulaComponent_str}'))



Predefined_code='''

rozdzial = KineticEnergyComponent(Engine())
rozdzial_kod = KineticEnergySymPyCodeComponent(Engine())
    
doc_final.append(rozdzial)
doc_final.append(rozdzial_kod)

'''


class PredefinedSectionComponent(DocumentGenerationComponent):

    #title="Dodawanie predefiniowanej sekcji"
    title="Adding a predefined section"



    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        
        #display(ReportText('Istnieje również możliwość wprowadzenia zpredefiniowanych sekcji. Możemy posłużyć się istniejącym komponentem dynamicznym na przykładzie silnika. Wykorzystamy zarówno opis teoretyczny jaki i reprezentację kodu służącą do wygenerowania równania. Aby użyć ich w dokumencie musimy wywołać te sekcje:'))
        display(ReportText('It is also possible to introduce predefined sections. We can use an existing dynamic component in the engine example. We will use both the theoretical description and the code representation used to generate the equation. To use them in the document we need to call these sections:'))
        from dynpy.models.mechanics import Engine
        display(ObjectCode(Predefined_code))
        display(KineticEnergyComponent(Engine()))
        display(KineticEnergySymPyCodeComponent(Engine()))
    
    

    
    

    

doc_comp_str = (
'''

doc = TikzCaSCStandalone(default_filepath='./output/Tytul_Dokumentu')
Picture._default_width = NoEscape('0.8\\columnwidth')

''')
class DocumentComponent(DocumentGenerationComponent):

    title="Tworzenie dokumentu"


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
        display(ReportText('Aby stworzyć dokument i nadać mu nazwę należy wywołać kod, jak w komórce nr 2 jupytera, tj. cell [2]:'))
        display(GuideCode(f'{doc_comp_str}'))
        



aut_indx_str1 = (
'''

x=Symbol('x')
y=Symbol('y')
eq=Eq(y,2*x)
display(ReportText('.. rozpatrujemy następujące równanie'))
display(SympyFormula(eq))
display(ReportText(f'Po wstawieniu danych do równania ({AutoMarker(eq)}) otrzymamy ...]'))

''')
aut_indx_str2 = (
'''

obrazek=Picture('./dynpy/models/images/engine.png')
display(ReportText('... rozpatrujemy następujący mechanizm'))
display(obrazek)
display(ReportText(f'Jak w wcześniej przedstawionym schemacie ({AutoMarker(obrazek)}) możemy zauważyć ...'))

''')
aut_indx_str3 = (
'''
phi=Function('phi')(t)
r=Symbol('r', positive=True)
a=Symbol('a', positive=True)
b=Symbol('b', positive=True)
h=Symbol('h', positive=True)
omg=Symbol('omega')
ob=Symbol('|OB|')
dane={
    phi:pi/2,
    omg:25,
    r:0.15,
    a:0.2,
    b:0.6,
    h:0.4,
    ob:0.42
}
tabelka=LatexDataFrame.formatted(
        data=dane,
        index=['Wartość']).applymap(lambda x: f'${latex(x)}$').rename_axis(
            'Symbol',
            axis=1).reported(caption='Dane z zadania')
display(ReportText('... rozpatrujemy następujące dane'))
display(tabelka)
display(ReportText(f'Podstawiając dane podane w zadaniu ({AutoMarker(tabelka)}) możemy otrzymać wynik ...'))

''')



class AutomarkerIntroComponent(DocumentGenerationComponent):

    title="Wstęp do wykorzystywania narzędzia AutoMarker"


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
        display(ReportText('Aby w raporcie referować do wcześniej wyświetlanych równań/obrazków/tabelek używamy narzędzia AutoMarker. W przypadku równania         przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy zdefiniować potrzebne symbole oraz równanie, a następnie               wyświetlić je. Ważne jest to, aby równanie zostało zdefiniowane przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonego           równania używamy AutoMarkera. Przkładowy kod pokazany został poniżej:'))
        display(GuideCode(f'{aut_indx_str1}'))
        display(ReportText(' W przypadku obrazka przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy zdefiniować zmienną           która zostanie przypisana do naszego obrazka stworzonego za pomocą klasy Picture (podając odpowiednią ścieżkę), a następnie wyświetlić ją. Ważne           jest to, aby zmienna została zdefiniowana przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonego obrazka ponownie używamy         AutoMarkera w ten sam sposób. Przykładowy kod pokazany został poniżej:'))
        display(GuideCode(f'{aut_indx_str2}'))
        display(ReportText(' W przypadku tabelki przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy utowrzyć zmienną która         zostanie przypisana do naszej tabelki(przykładowa tabelka pokazana sekcje wyżej), a następnie wyświetlić ją. Ważne jest to, aby zmienna została             zdefiniowana przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonej tabelki ponownie używamy AutoMarkera w ten sam sposób.         Przykładowy kod pokazany został poniżej:'))
        display(GuideCode(f'{aut_indx_str3}'))


class PictureComponent(DocumentGenerationComponent):

    title="Picture component implementation"



    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
        #display(ReportText('Picture is a class that represents a figure environment'))
        display(ReportText('''Picture is a class that represents a figure environment. To add picture to document simply call Picture() class. Possible argumnts in Picture class are:
        \n path - define path to the picture
        \n position - define picture positioning acording to latex 
        \n width - define picture width in cm
        \n height - define picture height in cm
        \n caption - caption of the picture
        \n Example picture class call: \n
        '''))
        display(ObjectCode('''Picture('pic.PNG', position = 'h', caption = 'Example text', height = '10', width = '20')'''))
        pic = Picture('./dynpy/utilities/components/guides/images/1.png', caption = 'Example text', width = '10')
        display(pic)



CodeEmbeddingComponent_str1 =(
'''

display(ObjectCode(GeometrySceneDG))

''')


CodeEmbeddingComponent_str2 =(
"""

class GeometrySceneDG:


    ax_2d = None
    ax_3d = None

    #def __init__(self,height=12,width=9,figsize=(12,9)):
    def __init__(self, init_3d=(30, 10), height=12, width=16, figsize=(12, 9)):
        
        
        
        plt.figure(figsize=figsize)
        ax_2d = plt.subplot(121)
        ax_2d.set(ylabel=(r'<-x | z ->'), xlabel='y')

        plt.xlim(0, width)
        plt.ylim(-height, height)
        plt.grid(True)

        ax_2d.set_yticks(range(-12, 12, 1))
        ax_2d.set_yticklabels(
            list(map(lambda tick: str(abs(tick)), range(-12, 12, 1))))

        ax_2d.set_xticks(range(0, 16, 1))
        ax_2d.set_xticklabels(
            list(map(lambda tick: str(abs(tick)), range(0, 16, 1))))

        ax_3d = plt.subplot(122, projection='3d')
        ax_3d.set(xlabel='x', ylabel='y', zlabel='z')

        plt.xlim(0, 16)
        plt.ylim(0, 16)

        ax_3d.set_zlim(0, 16)

        ax_3d.view_init(*init_3d)
        plt.tight_layout()

        self.__class__.ax_2d = ax_2d
        self.__class__.ax_3d = ax_3d
        self.__class__.size_3d=16

""")

class CodeEmbeddingComponent(DocumentGenerationComponent):

    title="Obsadzanie kodu"



    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
    
        display(ReportText('Klasa umożliwiająca dodanie sformatowanego fragmentu kodu pisanego programu. Jako przykład wykorzystana zostanie klasa **GeomerySceneDG** z biblioteki **dgeometry**.'))
        display(GuideCode(f'{CodeEmbeddingComponent_str1}'))
        display(GuideCode(f'{CodeEmbeddingComponent_str2}'))
        
ReportTextComponent_str=(
'''
display(ReportText('Tutaj wpisujemy tekst, który chcemy wyświetlić w naszym dokumencie'))
    
''')

class MarkdownComponent(DocumentGenerationComponent):
    
    title='Klasa Markdown'
    

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        
       
        #implement reporting activieties here



        display(ReportText('Modułem o rozszerzonej funkcjonalności jest funkcja **Markdown**. Funkcjonalność dostępna po wywołaniu kodu::'))
    #display(GuideCode(f'{MarkdownComponent_str}'))
        display(Markdown(
    '''

    - Tworzenie elementów wypisanych w punktach

    - Załączanie lącza internetowego

    - Wyswietlanie kodu

    '''
                ))
        display(ReportText('Aby wymienić w punktach należy napisac następujący kod:'))
        # markdown1=Picture('./Image/markdown1.jpg')

        display(Markdown(
        '''
            display(Markdown("""

                    - Tworzenie elementów wypisanych w punktach

                    - Załączanie lącza internetowego

                    - Wyswietlanie kodu

            """))
        '''))

        # display(markdown1)

        display(ReportText('Aby załączyć link, na przyklad do spotkania:'))

        display(Markdown(
        '''
            display(Markdown("""Link do spotkania: [SPOTKANIE]
            (https://wutwaw.sharepoint.com/sites/EfektywnyPythonbySiMR/_layouts/15/stream.aspx?id=%2Fsites%2EfektywnyPythonbySiMR%2FShared%20Documents%2FGeneral%2FRecordings%2FSpotkanie%20na%20kanale%20Ogólnym%2D20230628%5F175853%2DNagrywanie%20spotkania%2Emp4)
            """))
        '''))



        # markdown2=Picture('./Image/mrkdwn2.jpg')

        # display(markdown2)

        display(ReportText('Aby wyswietlic kod, przed apostrofami nalezy zapisać literę $f$. Przykładowo :'))

        display(Markdown(
        '''
            display(Markdown("""

                    - Tworzenie elementów wypisanych w punktach

                    - Załączanie lącza internetowego

                    - Wyswietlanie kodu

            """))
        '''))


        # markdown3=Picture('./Image/mrkdwn3.jpg')

        # display(markdown3)

        display(ReportText('Ostatni krok to zaapendowanie sekcji i utworzenie dokumentu pdf :'))
        display(Markdown(
        f'''


            doc_final = MechanicalCase('./output/nazwa_dokumentu',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)
            doc_final.packages.append(Package('natbib', options=['numbers']))
            doc_final.packages.append(Package('booktabs'))
            doc_final.packages.append(Package('float'))
            doc_final.packages.append(Package('siunitx'))


            doc_final.append(sekcja)
            doc_final.append(sekcja1)
            doc_final.append(sekcja2)
            doc_final.append(sekcja3)
            doc_final.append(sekcja4)
            doc_final.append(sekcja5)

            doc_final.generate_pdf()



        '''
                    ))


class ReportTextComponent(DocumentGenerationComponent):

    title="Klasa ReportText"


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        
       
        #implement reporting activieties here
        
        display(ReportText('Narzędziem o podstawowej funkcjonalności jest klasa **ReportText** - umożliwia dodanie i wyświetlanie tekstu, który ma zostać umieszczony w raporcie. Funkcjonalność ta jest dostępna po wywołaniu kodu:'))
        display(GuideCode(f'{ReportTextComponent_str}'))
        
TablesCreationComponent_str=(
'''
from dynpy.models.mechanics.tmac import SDOFWinchSystem
eoms= SDOFWinchSystem()._eoms[0]
eoms_eq_raw=Eq(eoms,0)
data=SDOFWinchSystem().get_random_parameters()
Ms0=Symbol('M_s0')
b=Symbol('b')
omg=Symbol('varphi')
ph1=Symbol('Phi')

winch = SDOFWinchSystem()

Ms=winch.M_s
t = winch.ivar

dane={
    Ms0:4700,
    b:1500,
    SDOFWinchSystem().m0:10}



winch_up=SDOFWinchSystem().subs(Ms,Ms0-b*SDOFWinchSystem().phi.diff()).subs(data).subs(dane)
nowe_ode=winch_up._ode_system
solution=nowe_ode.solution
steady_sol=nowe_ode.steady_solution[0]
t_span = np.linspace(0,100,200)
phi=winch.q[0]



units_dict = {phi : ureg.radian,phi.diff(t):ureg.radian,t:ureg.second }
LatexDataFrame.set_default_units(units_dict)

general_sol_matrix=nowe_ode.solution.with_ics([0,0])
table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)
table=table_eq[[phi.diff(t)]]
table_for_report = table[0::15].to_latex_dataframe()

display(Markdown(' Do tablicy `table_eq` zapisujemy równanie ogólne dla zmiennej czasowej `t_span`   '))
display(Markdown(
f"""
    
    table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)

"""))
display(Markdown(f' Do tablicy table przypisujemy wartości ${latex(omg)}$ powstałe w wyniku zróżniczkowania ${latex(ph1)}$ w tablicy `table_eq` '))
display(Markdown(
f"""
    table=table_eq[[phi.diff(t)]]
   
"""))
display(ReportText(' W następnym kroku wybieramy zakres wartości który chcemy umieścić w tabelce przy pomocy naiwasów kwadratowych. Nastepnie przekształcamy tabelę na zgodną z LateX przy pomocy metody: to\_latex\_dataframe()  '))
display(Markdown(
f"""

    table_for_report = table[0::15].to_latex_dataframe()

"""))
display(ReportText(' Naszą tabelkę opisujemy w tekście, pożądane jest użycie AutoMarkerów które zostały opisane w poniższej sekcji tego dokument '))
display(Markdown(
f"""

     display(Markdown(f'Wyniki symulacji zostały zestawione w tabeli {AutoMarker(table_for_report)}'))

"""))
display(ReportText(' Aby wyświetlic tabelkę korzystamy z  klasy display. Metoda .reported umożliwia nam nadanie tabelce podpisu. '))
display(Markdown(
f"""

    display(table_for_report.reported(caption='Tabela'))

"""))
    
''')

class TablesCreationComponent(DocumentGenerationComponent):

    title="Tworzenie tabelek"


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        
       
        #implement reporting activieties here
        
        display(ReportText('Aby utworzyć tabelkę najpierw musimy posiadać dane z obliczeń na interesującym nas systemie dynamicznym:'))
        display(GuideCode(f'{TablesCreationComponent_str}'))
        
        
kod_2='''
from dynpy.utilities.components.guides.en import ReportComponent
    from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, GuideCode, LatexDataFrame
    from dynpy.utilities.components.guides.en import ReportCompImplementationComponent,ReportCompImplementationIssueComponent
    from sympy import *
    from pint import UnitRegistry
    import pandas as pd
'''
class LibrariesImportComponent(DocumentGenerationComponent):

    #title="Klasa LibrariesImportComponent i jej zastosowanie"
    title="LibrariesImportComponent class and its use"

    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        #display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(ReportText('Before working in jupiter, import libraries'))
        display(ObjectCode(kod_2))
        #display(ReportText('Dzięki temu możemy korzystać z przygotowanych wcześniej klas i modułów'))
        display(ReportText('This allows us to use previously prepared classes and modules'))
        
markdown_reporting_str=("""
Koncepcję raportowania można przedstawić prosto w następujących punktach:  
​
* Biblioteka do raportowania ma budowę modułową. W kodowaniu obiektowo zorientowanym oznacza to, że składa się z szeregu klas. Poszczególne klasy to implementacje programistyczne standardowych elementów, z których składa się każdy raport. To znaczy zaimplementowano klasę do tworzenia dokumentu oraz klasy do kreowania typowych elementów dokumentu, jak rozdział, akapit (tj. blok tekstu), wzór, obraz, wykres.  
​
* Przyjęta została współbieżna zasada tworzenia raportów. Według tej zasady, aby zaraportować dany element najpierw musi on zostać wyświetlony w jupyterze, używając nomenklatury informatycznej - wydrukowany na wyjściu komórki jupytera, tzn. jako cell output.  
​
* Główna implementacja kodu, czyli używane klasy, znajdują się w module `dynpy.utilities.report`
​
* Do zbudowania standardowego raportu używa się następujących klas:  
    
    * class `Document`
    
    * class `CurrentContainer`
    
    * class `ReportText`
    
    * class `Markdown`
    
    * class `Picture`
    
    * class `SympyFormula`
​
* W efekcie tworzony jest dokument wynikowy w postacie pliku pdf, wygenrowany w języku zanczników LaTeX w oparciu o bazowy plik jupytera (*.ipynb), który zawiera przejrzysty podgląd dokumentu wynikowego.
​
""")

class ReportingModuleIntroComponent(DocumentGenerationComponent):

    title="Implementacja komponentów raportujących"



    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
        display(Markdown(markdown_reporting_str))
        
SympyFormulaComponent_str='''
from sympy import *

    M = Symbol('M')

    t = Symbol('t')

    m_e = Symbol('m_e')

    e = Symbol('e')

    z = Function('z')

    varphi = Function('varphi')

    Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2

    display(Eq(Symbol('T'),Ek))
'''


class SympyFormulaComponent(DocumentGenerationComponent):

    #title = "Klasa SympyFormula i jej zastosowanie"
    title = "SympyFormula class and its use"


    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        #implement reporting activieties here
        #display(ReportText('Przykładowe zastosowanie komponentu do obliczeń symbolicznych, na przykładzie energii kinetycznej:'))
        display(ReportText('Example application of the component for symbolic calculations, using kinetic energy as an example:'))
        display(GuideCode(f'{SympyFormulaComponent_str}'))
        
        
report_formating_str=("""
Należy zwrócić uwagę na następujęce aspekty:

- Cele powinny być poopisywane Markdownami żeby mozna bylo sie swobodnie poruszać po dokumencie przez Table of Contents.

- Appendowanie sekcji do dokumentu powinno znajdować się na końcu dokumentu, ponieważ aby zbudować dokument nie trzeba uruchamiać wszystkich cele po kolei.

- Nie można obsadzać nieformatowalnych kodów w 'ReportText'. 

- Jeżeli chcemy obsadzić jednoliniowy kod, musimy użyć ` (apostrof pod tyldą). np: display(Markdown(\"`help`\")). Ten sposób nie zapewnia kolorowania kodu.

- Obsadzamy w tekście zmienne w ten sposób: `${latex(omega_r)}$`, aby globalnie mieć do nich dostęp, możliwość podmiany symbolu bez manualnego szukania w tekście.

- Oznaczamy równania za pomocą klasy `AutoMarker`.

- Gdy mamy zbyt długie równania, szczególnie części pod pierwiastkiem, należy wyciągnąć je do zmiennej i pokazać niżej. Alternatywną metodą na rozwiązanie tego problemu jest wstawienie równania w system dynamiczny.

- Przechowywane wstępne obliczenia / symbole należy umieścić w pliku.py lub wykonać za pomocą klasy jeśli taka istnieje.

Cel pracy i metodyka powinny charakteryzować się nastęującymi cechami:

- stwierdzać fakty,

- być napisane w jednym czasie (przeszły lub teraźnieszy)

- nie uprzedzać faktów, które następują w dalszej częsci pracy,

- mieć charakter planu lub zbioru założeń/oczekiwań np. założono, że zostanie wykonana analiza wrażliwości; stwierdzono, że takie podejście zapewni poprawność modelu

Elementy, które wymagają usprawnienia w szablonie pracy:

- streszczenia w szablonie pracy dyplomowej,

- naprawienie obsługi komendy `\cite` w klasie `Markdown`

- stosowanie `{latex}`

- stosowanie klasy `AutoMarker`

- należy w szablonie PD przenieść Appendowanie na koniec 

""")

report_formating_en_str=('''
Pay attention to the following aspects:

- Objectives should be described in Markdown to be able to move freely through the document via Table of Contents.

- Appending sections to the document should be at the end of the document, because it is not necessary to run all the targets one by one to build the document.

- You can't cast unformatable codes in 'ReportText'. 

- If you want to cast a one-line code, you must use ` (apostrophe under the tilde). Ex: display(Markdown(`help`)). This way does not provide code coloring.

- We staff the variables in the text like this: `${latex(omega_r)}$` to globally access them, the ability to substitute the symbol without manually searching in the text.

- We mark equations using the `AutoMarker` class.

- When we have equations that are too long, especially the parts under the root, they should be extracted into a variable and shown below. An alternative method to solve this problem is to insert the equation in the dynamic system.

- The stored precalculations / symbols should be placed in a.py file or made using a class if one exists.

The purpose of the work and the methodology should have the following characteristics:

- state facts,

- be written in one tense (past or present)

- not anticipate the facts that follow in the later part of the work,

- have the character of a plan or a set of assumptions/expectations, e.g., it was assumed that a sensitivity analysis would be performed; it was stated that this approach would ensure the correctness of the model

Elements that need improvement in the thesis template:

- abstracts in the thesis template,

- fixing the handling of the `\cite` command in the `Markdown` class.

- use of `{latex}`

- application of `AutoMarker` class.

- move Append to the end in the PD template 

'''
)

class ReportFormattingGuidelinesComponent(DocumentGenerationComponent):

    #title="Wytyczne odnośnie formatowania raportu"
    title="Report formatting guidelines"

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']    
        
        #implement reporting activieties here
        display(Markdown(report_formating_en_str))

##########################
##SYNCHRONIZACJA GITHUB### 
##########################
        
class GitSynchroPanelAccessComponent(ReportComponent):

    title="Dostep do panelu synchronizacji Github"

    @property
    def reported_object(self):


        default_data = {'classname':'ReportingModuleIntroComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':359,
                       }


        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data

        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj
    
    
    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        display(ReportText(f''' Repozytoria Github'a zmieniają się dosyć dynamicznie, a w przypadku, gdy korzysta z niego więcej niż jedna osoba dobrze by było aby każdy miał dostęp do najnowszej wersji. W tym celu przeprowadza się regularnie synchronizację danego rezpozytorium. '''))

        display(ReportText('Jak się dostać do panelu synchronizacji?'))

        display(ReportText(f''' Wybierz projekt, który chcesz zsynchronizować.'''))
        display(Picture('./dynpy/utilities/components/guides/images/1.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('11cm'),caption='Wybór projektu z listy'))

        display(ReportText(f'''Przejdź do foledru, który chcesz zsynchronizować. W naszym przypadku zazwyczaj aktualizujemy moduły.'''))
        display(Picture('./dynpy/utilities/components/guides/images/2.png', position = 'H', height= NoEscape('18cm'), width = NoEscape('18cm'),caption='Wybór folderu do synchronizacji'))

        display(ReportText(f'''Wybierz ikonę VS Code - zostaniesz przeniesiony do nowego okna.'''))
        display(Picture('./dynpy/utilities/components/guides/images/3.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('11cm'),caption='Przycisk przenoszący do VS Code'))

        display(ReportText(f'''Będąc już w panelu VS Code przjedź do zakładki 'Source control' (trzecia od góry), to tu będą przeprowadzane wszystkie operacje.'''))
        display(Picture('./dynpy/utilities/components/guides/images/4.png', position = 'H', height= NoEscape('18cm'), width = NoEscape('18cm'),caption='Source control w panelu VS Code'))
        

example_commit_str=('''
Przykładowy Commit:

# Please enter the commit message for your changes. Lines starting
# with '#' will be ignored, and an empty message aborts the commit.
#
# On branch master
# Your branch is up to date with 'origin/master'.

# Changes to be committed:
#   modified: models/mechanics/disk.py
#   modified: utilitied/components/ode/en.py
#   modified: models/mechanics/bridge.py
''')

example_commented_commit_str=('''
# Please enter the commit message for your changes. Lines starting
# with '#' will be ignored, and an empty message aborts the commit.

Here place the title of the Commit

On branch master
Your branch is up to date with 'origin/master'.

Changes to be committed:
  modified: models/mechanics/disk.py - class added
  modified: utilitied/components/ode/en.py - class removed
  modified: models/mechanics/bridge.py - minor changes
''')

example_commit_comments_str=('''
class added
class fixed
class removed
minor changes
method added
method fixed
method removed
''')

class GitSynchroIntroComponent(GitSynchroPanelAccessComponent):

    title="Wprowadzenie do synchronizacji repozytorium Github"

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here

        display(ReportText('Synchronizacja repozytorium.'))

        display(ReportText('Przypadek 1. : "Commit" po lewej stronie NIE działa'))

        display(ReportText(f'''Jeśli "Commit" nie działa to kliknij od razu "Sync changes" lub kółeczka synchronizacji na górze'''))
        display(Picture('./dynpy/utilities/components/guides/images/10.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='Kółeczko synchronizacji'))
        display(Picture('./dynpy/utilities/components/guides/images/9.png', position = 'H', height= NoEscape('8cm'), width = NoEscape('6cm'),caption='Przycisk, na którym pojawia się "Commit" albo "Sync changes"'))

        display(ReportText('Przypadek 2. : "Commit" po lewej stronie wyświetla się normalnie'))

        display(ReportText(f'''Kliknij "Commit".'''))
        display(Picture('./dynpy/utilities/components/guides/images/5.png', position = 'H', height= NoEscape('8cm'), width = NoEscape('6cm'),caption='Przycisk "Commit"'))

        display(ReportText(f'''Otworzy się opis Commit'a - odkomentuj linijki z modyfikacjami, czyli wszystko od około 4. linii. '''))
        # display(Picture('./dynpy/utilities/components/guides/images/6.png', position = 'H', height= NoEscape('13cm'), width = NoEscape('13cm'),caption='''Opis Commit'a '''))

        display(ObjectCode(example_commit_str))

        display(ReportText(f'''Po odkomentowaniu dopisz tytuł synchronizacji oraz krótkie opisy poszczególnych zmian - jest to ważne w kwestii komunikacji między osobami pracującymi nad projektem. '''))
        # display(Picture('./dynpy/utilities/components/guides/images/7.png', position = 'H', height= NoEscape('13cm'), width = NoEscape('11cm'),caption=''' Przykladowy opis Commit'a '''))
        display(ReportText('Przykładowy opisany Commit:'))
        display(ObjectCode(example_commented_commit_str))
        display(ReportText('Spis przykładowych komenatrzy do Committów:'))
        display(ObjectCode(example_commit_comments_str))


        display(ReportText(f''' Dla pewności poprawności wprowadzonych zmian warto potwierdzić Commita na chatcie do synchronizacji na Slacku'''))
        display(Picture('./dynpy/utilities/components/guides/images/8_5.png', position = 'H', height= NoEscape('6cm'), width = NoEscape('6cm'),caption='''Kanał do synców '''))

        display(ReportText(f'''Zatwierdź Commit'a klikjąc tick w prawym górnym rogu. Jeśli przycisku nie ma, przełącz karty i spróbuj ponownie.'''))
        display(Picture('./dynpy/utilities/components/guides/images/8.png', position = 'H', height= NoEscape('20cm'), width = NoEscape('20cm'),caption='Zatwierdzanie zmian'))

        display(ReportText(f'''Po zatwierdzeniu, ze strony Github'a takie zmiany wyglądają następująco '''))
        display(Picture('./dynpy/utilities/components/guides/images/7_5.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='''Zatwierdzone zmiany z perspektywy Github'a'''))
        display(ReportText('Wystąpienie komunikatu merge'))

        display(ReportText(f'''Taki komunikat będzie oznacza błąd na wyższym szczeblu - nie ruszaj nic więcej i odezwij się do Pana Chilińskiego lub Pana Sierocińskiego. '''))
        display(Picture('./dynpy/utilities/components/guides/images/11.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='Error'))
        
class UsageOfGitHubInterfacesComponent(GitSynchroPanelAccessComponent):

    title="Wprowadzenie do GitHubInterface"

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here

        display(ReportText('Obsluga rytynowych czynności w dynpy jest wykonanywana przez klasę GitHubInterface'))

        display(ObjectCode(GitHubInterface))

        display(ReportText('Więcej informacji możesz uzyskać w helpie'))

        display(ObjectCode("help(GitHubInterface)"))
        display(ObjectCode(help(GitHubInterface)))

class AlgebraicExpressionComponent(ReportComponent):

    title = "Algebraic expression example"

    def append_elements(self):
        display(ReportText('Example Algebraic expression in Sympy library:'))
        display(GuideCode(f'{SympyFormulaComponent_str}'))
        
        from sympy import Symbol, Function

        M = Symbol('M')

        t = Symbol('t')

        m_e = Symbol('m_e')

        e = Symbol('e')

        z = Function('z')

        varphi = Function('varphi')

        Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2
        
        display(ReportText('After executing above code you can display kinetic energy algebraic expression:'))
        
        display(Eq(Symbol('T'),Ek))

ODEIcsComponent_dec_str = ('''
m, c, k, F, mu = symbols('m c k F mu', positive = True)
t = Symbol('t')
x = Function('x')(t)

subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))
display(ode_1)
''')
        
class ODEIcsComponent(ReportComponent):

    title = "ODESystem solution with set of initial conditions (ICS) component presentation"

    def append_elements(self):

        display(ReportText('First declare ODESystem component:'))
        
        display(ObjectCode(ODEIcsComponent_dec_str))
        
        from dynpy.solvers.linear import AnalyticalSolution, ODESolution, ODESystem
        import sympy
        
        m, c, k, F, mu = symbols('m c k F mu', positive = True)
        t = Symbol('t')
        x = Function('x')(t)

        subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

        ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))
        
        display(ReportText('After executing above code bellow equation will appear'))
        
        display(ode_1)

        display(ReportText('To present solution with initial set of conditions, first you need to calculate solution:'))
        
        display(ObjectCode('sol_sym_1_pre_ics = ode_1.subs(subs_data).solution'))
        
        sol_sym_1_pre_ics = ode_1.subs(subs_data).solution
        
        display(sol_sym_1_pre_ics)
        
        display(ReportText('To calculate solution with ICS, declare ics vector, and run with_ics method as shown bellow:'))
        
        display(ObjectCode('''
ics = [1, 2]
sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
display(sol_sym_1)'''))
        
        ics = [1, 2]
        sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
        display(sol_sym_1)
        
ODEIcsCodeComponent_str =('''
m, c, k, F, mu = symbols('m c k F mu', positive = True)
t = Symbol('t')
x = Function('x')(t)

subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))

sol_sym_1_pre_ics = ode_1.subs(subs_data).solution

ics = [1, 2]

sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
''')

class ODEIcsCodeComponent(ReportComponent):

    title = "ODESystem solution with set of initial conditions (ICS) component presentation"

    def append_elements(self):

        display(ReportText('For bellow ODESystem component example:'))
        
        display(ObjectCode(ODEIcsCodeComponent_str))
        
        from dynpy.solvers.linear import AnalyticalSolution, ODESolution, ODESystem
        import sympy
        
        m, c, k, F, mu = symbols('m c k F mu', positive = True)
        t = Symbol('t')
        x = Function('x')(t)

        subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

        ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))

        sol_sym_1_pre_ics = ode_1.subs(subs_data).solution

        ics = [1, 2]

        sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
        
        display(ReportText('To get call code for ICS, run bellow method:'))
        
        display(ObjectCode("sol_sym_1.get_ics()"))
        
        display(sol_sym_1.get_ics())
        
git_com_str=("""
from {module} import {classname}
{classname}({class_field});
""")


class IssueFeedbackComponent(ReportComponent):

    title="Komentarz zamykajacy issue z komponentami"

    @property
    def reported_object(self):
        from ....dynamics import LagrangesDynamicSystem

        default_data = {'classname':'Component',
                       'module':'dynpy.utilities.components.guides.en',
                       'field':'None',
                       'target':'`ReportComponent` class',
                       'issue_no':359,
                       }

        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        if isinstance(self._reported_object, ReportComponent) or isinstance(self._reported_object, LagrangesDynamicSystem):
            return {**default_data,'classname':self._reported_object.__class__.__name__,'module':self._reported_object.__class__.__module__}
        
        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data
        
        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        display(ReportText('The code was checked with the following call: '))
#        display(ObjectCode(git_com_str.format(classname=classname)))
        display(ObjectCode(git_com_str.format(classname=classname,module=class_module,class_field=class_field)))
    
old_class_code = '''class ExemplaryOldImplementedSystem(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('\\varphi')

    c0 = Symbol('c0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        self.ivar = ivar

        self.qs = [self.phi]

        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity_ = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)

        composed_system = self._mass_x + self._mass_y + self._gravity_

        super().__init__(composed_system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict
    
    def get_numerical_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict

    def max_static_force(self):
        return S.Zero

    def max_dynamic_force(self):
        return S.Zero'''

new_class_code = '''class MaterialPointMovement(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('phi')

    c0 = Symbol('c0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        self.ivar = ivar

        self.qs = [self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}


        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)



        components['_mass_x']=self._mass_x
        components['_mass_y']=self._mass_y
        components['_gravity']=self._gravity


        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict
    
    def get_numerical_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict    

    def max_static_force(self):
        return S.Zero

    def max_dynamic_force(self):
        return S.Zero'''


class CodeRefactorComponent(NumericalAnalysisSimulationComponent):

    title = "DynSys system code refactor component"

    def append_elements(self):

        from dynpy.models.mechanics import Engine
        
        system = self.reported_object
        
        display(ReportText('Refactoring is a systematic process of improving code without creating new functionality that can transforma mess into clean code and simple design.'))
        
        
        display(ReportText('Example of old DynSys component before refactoring:'))
        
        display(ObjectCode(old_class_code))
        
        display(ReportText('Example of current DynSys component after refactoring:'))
        
        display(ObjectCode(new_class_code))
        
        
        display(ReportText('After calling system_description() method for Engine model bellow result will be created:'))
        
        str_sys_des = system.system_description().replace('\n', '\n\n').replace('\t', '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')

        display(ReportText(str_sys_des))
        
        display(ReportText('\n\nCalling system_parameters() method will create bellow list:'))
        
        lst_sys_par = system.system_parameters()

        for element in lst_sys_par:

            display(ReportText(f"* {element}"))
            
class DynSysSummaryComponent(ReportComponent):

    title = "DynSys system code refactor component"

    def append_elements(self):

        from dynpy.models import mechanics


        display(ReportText('DynSys system consists of the following models:'))
        
        mods=inspect.getmembers(mechanics, inspect.ismodule)

        for nam,mod in mods:

            mods_tmp = [name for name,cls in inspect.getmembers(mod, inspect.isclass) if issubclass(cls,mechanics.principles.ComposedSystem)]
            classes_str=',\n - '.join(mods_tmp)
            if mods_tmp != []:
                display(ReportText(f' {nam}'))
                display(ReportText(f'- {classes_str}'))
                
        display(ReportText('All of this models should be refactored according to bellow example:'))
        CodeRefactorComponent(None)

ds_overview_str=("""
from dynpy.utilities.documents.document import DynSysOverviewReport
from dynpy.utilities.components.guides.en import *
from dynpy.models.mechanics.pendulum import Pendulum

DynSysOverviewReport(reported_object=Pendulum());
""")

dynsys_comp_check_str = '''from dynpy.models.mechanics import ForcedSpringMassSystem

comp = ForcedSpringMassSystem()

DynamicSystemCompletenessCheckComponent(comp)'''

basicsymcomp_str = '''from dynpy.models.mechanics.tmac import SDOFWinchSystem

comp = SDOFWinchSystem()

BasicSymComponent(comp)'''

coderefactorcomp_str = '''from dynpy.utilities.components.guides.en import CodeRefactorComponent

CodeRefactorComponent(None)'''

dynsysimpcomp_str = '''from dynpy.utilities.components.guides.pl import DynSysImplementationComponent

DynSysImplementationComponent(None)'''

diffsimcomp_str = '''from dynpy.utilities.components.guides.en import DifferentSimulationsComponent

comp = SDOFWinchSystem()

DifferentSimulationsComponent(comp)'''

dynintrocomp_str = '''from dynpy.utilities.components.guides.en import *
DynSysIntroComponent(None)'''

dynsumcomp_str = '''from dynpy.utilities.components.guides.en import *

DynSysSummaryComponent(None)'''

dynmetuscomp_str = '''from dynpy.utilities.components.guides.en import DynamicSystemMethodsUsageComponent

comp = SDOFWinchSystem()

DynamicSystemMethodsUsageComponent(comp)'''
dyncompov_str = '''Guide consist of bellow components which can be called separatelly with most of avaible dynsys models:
*DynamicSystemCallComponent
*NumericalAnalysisSimulationComponent
*AnalyticalSimulationComponent
*DynSysCodeComponent
*IssuePreparationComponent'''

class DynSysOverviewUsageComponent(ReportComponent):

    title = "Introduction to usage of DynSysOverviewReport"

    def append_elements(self):

        from dynpy.models import mechanics


        display(ReportText('This guide concers in-deepth analysis of exsisting dynamic system that are basic elements of `DynPy` library.'))
        display(ReportText('Basic call of it is as follows and runs default dynamic system which is `ForcedSpringMassSystem'))
        display(ObjectCode(ds_overview_str))
        
        display(ReportText(dyncompov_str))
        
        #display(ReportText('This guide concers in-deepth analysis of exsisting dynamic system that are basic elements of `DynPy` library.'))
        #display(ReportText('It can be run with any available dynamic system.'))
        '''
        #DynamicSystemCompletenessCheckComponent
        display(ReportText('This component calls and checks all key elements of the dynamic system'))
        display(ReportText('It can be run with any available dynamic systems.'))
        display(ObjectCode(dynsys_comp_check_str))
        
        #BasicSymComponent
        display(ReportText('This component shows how to do basic simulations for dynamic system.'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on below example call.'))
        display(ObjectCode(basicsymcomp_str))
        
        #CodeRefactorComponent
        display(ReportText('This component shows how to refactor DynSys component according to latest standard.'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(coderefactorcomp_str))
        
        #DifferentSimulationsComponent
        display(ReportText('This component explains how to implement dynamic system with ComposedSystem class'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(diffsimcomp_str))
        
        #DynSysImplementationComponent
        display(ReportText('This component explains how to implement dynamic system with ComposedSystem class'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on below example call.'))
        display(ObjectCode(dynsysimpcomp_str))
        
        #DynSysIntroComponent
        display(ReportText('This component is simple introduction to using the DynSys module'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on baeow example call.'))
        display(ObjectCode(dynintrocomp_str))
            
        #DynSysSummaryComponent
        display(ReportText('This component shows idea of code refactoring and additionally displays list of available dynamical systems'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(dynsumcomp_str))
        
        #DynamicSystemMethodsUsageComponent
        display(ReportText('This component shows equations of motion and other methods avaible on DynSys model'))
        display(ReportText('It can be run with any available dynamic systems.'))
        display(ObjectCode(dynmetuscomp_str))
        '''

        