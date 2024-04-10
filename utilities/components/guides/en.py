from  ..mechanics import *

import pandas as pd
import numpy as np
from sympy import lambdify


miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

długosc_dnia_w_miesiacach_godziny = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia temperatura':srednie_temp_list}

df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)

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
            return None
            
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

        display(ReportText('Aby elementy z których korzystamy działay należy zaimportować potrzebne biblioteki, takie jak na załączonym obrazku oraz w formie tekstu który można skopiowac'))

        pic18 = Picture('./dynpy/utilities/components/guides/images/ImportowanieBibliotek.png', width='9cm')

        # display(pic18)

        display(ObjectCode(dynpy_imports_code))


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
        
        display(ReportText(' Możemy również wyswietlić podgląd takiego wahadła używając poniższej funkcji:  '))
        display(Picture('./dynpy/utilities/components/guides/images/preview.jpg'))
        display(ReportText(' Co więcej, używając poniższej komendy help, możemy otrzymać wszelkie informacje na temat konkretnego modelu:  '))
        display(Picture('./dynpy/utilities/components/guides/images/help.jpg'))

        display(ReportText('''Kod wywołania pomocy jest następujący:'''))
        
        display(ObjectCode(   code_dynsys_list_str   ))


        display(ReportText('Poniżej znajduje się lista systemów mechanicznych których można użyć:  '))

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

        display(ObjectCode(
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
from dynpy.utilities.report import (SystemDynamicsAnalyzer, DataPlot, AccelerationComparison, FFTComparison, ReportEntry, SimulationalBlock, ReportText, SimulationFFT, DataStorage, Markdown, SummaryTable, Picture, SympyFormula, SymbolsDescription, DescriptionsRegistry, ObjectCode, CurrentContainer)
from dynpy.solvers.linear import ODESystem
init_vprinting()
from dynpy.models.mechanics.tmac import SDOFWinchSystem

'''



class DynamicSystemCallComponent(ReportComponent):
    
    title="Wprowadzenie"


    def append_elements(self):

        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        display(ReportText('Wymagane importy bibliotek/poszczególnych klas do stworzenia raportu:'))


        display(ObjectCode(imports_code_str))


        display(ReportText('Klasa umożliwiająca rozwiązanie konkretnej problematyki (w tym przypadku problematyki związanej z opisem dźwigu):'))

        display(ObjectCode(
        """
        from dynpy.models.mechanics.tmac import SDOFWinchSystem 
        """))

        display(ReportText('Ścieżka do poszukiwanej klasy na platformie CoCalc:'))
        display(Picture('./dynpy/utilities/components/guides/images/sciezka_w.jpg'))

        display(ReportText('Sposób wywowałania preview klasy - tzw. podglądu:'))

        display(ObjectCode('''
        SDOFWinchSystem().preview()
        '''))
        
        
class DynamicSystemMethodsUsageComponent(ReportComponent):

    title="Wywolywanie rownan ruchu oraz innych metod"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        #from dynpy.models.mechanics import ForcedSpringMassSystem as SDOFWinchSystem

        eoms=system._eoms[0]
        eoms

        display(ReportText('Proces wywołowania równania ruchu za pomoca metody eoms:'))
#         display(ObjectCode('''eoms=SDOFWinchSystem()._eoms[0]
#         eoms'''
#         ))



        display(ReportText('Wynik jest następujący:'))
        #display(SympyFormula(eoms))


        eoms_eq=Eq(eoms,0)
        eoms_eq

        display(ReportText('Proces tworzenia równania z wcześniej wywołanego równania ruchu (równanie początkowo jest wywoływane jako ciąg znaków reprezentujący opis danego problemu, jednak nie jest on wówczas równaniem, które możemy zastosować do obliczeń):'))
#         display(ObjectCode('''eoms_eq=Eq(eoms,0)
#         eoms_eq'''
#         ))
        display(ReportText('Wynik jest następujący:'))
        display(SympyFormula(eoms))

        #display(Picture('./Image/eomseq_w.jpg'))

        slownik=system.get_random_parameters()

        

        display(ReportText('Tworzenie słownika (w tym przypadku z losowymi wartościami) umożliwiający nadanie wartości zmiennym:'))

        display(Markdown(
"""

    slownik=SDOFWinchSystem().get_random_parameters()
    slownik

"""
        ))
        #display(Picture('./Image/slownik_w.jpg'))

        display(ReportText('Rozwiązanie ogólne równania:'))
        #display(Picture('./Image/ogolne_w.jpg'))

        solution=system._ode_system.solution[0]
#         display(ObjectCode('''solution=SDOFWinchSystem()._ode_system.solution[0]
#         solution'''
#         ))
        display(SympyFormula(solution))

        display(ReportText('Rozwiązanie szczególne równania:'))



        steady_solution=system._ode_system.steady_solution[0]
#         display(ObjectCode('''steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]
#         steady_solution'''))
        display(SympyFormula(steady_solution))
    
    
    
steady_sol_str=(
'''

steady_solution_subs=steady_solution.subs(slownik)
steady_solution_subsS

''')

solution_table_str=(
'''

table=table_eq.numerized().compute_solution(t_span)
table_new=table[[phi]]

''')

report_table_define_str=(
'''

tabela = Section('Tabela')
CurrentContainer(tabela)

table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)
table=table_eq[[phi]]

table_for_report = table[0::15].to_latex_dataframe()

#Wyniki symulacji zostały zestawione w tabeli

table_for_report.reported(caption=('Tabela'))

''')

class SimulationsComponent(ReportComponent):
    
    title="Tworzenie wykresów"

    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Definicja wektoru czasu:'))
        #display(Picture('./Image/tspan_w.jpg'))

        display(ObjectCode('''t_span = np.linspace(0,100,200)'''))
        t_span = np.linspace(0,100,200)


        display(ReportText('Sposób tworzenia wykresów z tak zwanego palca:'))
        display(ReportText('Podstawienie danych z wcześniej stworzonego słownika:'))
        #display(Picture('./Image/stedisub_w.jpg'))

        display(Markdown(steady_sol_str))
        slownik=system.get_random_parameters()
        steady_solution=system._ode_system.steady_solution[0]
        steady_solution_subs=steady_solution.subs(slownik)
        display(steady_solution_subs)
        display(ReportText('Zastosowanie metody lambdify umożliwiającej konwersje funkcji do postaci anonimowej'))
        #display(Picture('./Image/lambdify_w.jpg'))


        eq_lambdify=lambdify(system.ivar,steady_solution_subs.simplify())
        display(ObjectCode('''eq_lambdify=lambdify(system.ivar,steady_solution_subs.simplify())'''))
        display(SympyFormula(eq_lambdify))


        display(ReportText('Stworzenie wykresu:'))
        #display(Picture('./Image/plot2_w.jpg'))
        pd.DataFrame(index=t_span, data=eq_lambdify(t_span)).plot()
        display(ObjectCode('''pd.DataFrame(index=t_span, data=eq_lamdify(t_span)).plot()'''))


        display(ReportText('Sposób tworzenia wykresów bezpośrednio z wcześniej wywołanej klasy:'))
        display(ReportText('Podstawienie słownika:'))
        #display(Picture('./Image/tableq_w.jpg'))

        table_eq=system._ode_system.steady_solution.subs(slownik)


        display(Markdown('''table_eq=system._ode_system.steady_solution.subs(slownik)'''))
        display(SympyFormula(table_eq))

        #jak tabele zrobic + ew jej output
        display(ReportText('Sposób tworzenia tabeli:'))
        #display(Picture('./Image/tabela_w.jpg'))
        phi = system.phi
        table=table_eq.numerized().compute_solution(t_span)
        table_new=table[phi]

        display(ObjectCode(solution_table_str))
        display(table)
        display(table_new)

        display(ReportText('Sposób tworzenia wizualizacji zdefiniowanej tabeli:'))
        #display(Picture('./Image/tabelka_w.jpg'))

        display(ObjectCode('''table_new'''))

        display(ReportText('sposób tworzenia tabeli wyświetlanej w gotowym raporcie:'))
        #display(Picture('./Image/tworzenie_tabeli_w.png'))

        display(Markdown(report_table_define_str))

        display(ReportText('Wizualizacja wykresu::'))
        #display(Picture('./Image/tikzplot_w.jpg'))

        display(ObjectCode('''table_new.to_pylatex_tikz().in_figure()'''))    

section_development_str=(
'''
intro=Section('Wprowadzenie')
CurrentContainer(intro)
display(ReportText('schemat układu'))
schemat=Picture('picture_winch.png')
display(schemat)
''')

forming_equation_str=(
'''
solution_subs = SDOFWinchSystem()._ode_system.steady_solution.subs(slownik)
eq_sol=Eq(solution_subs,0)
eq_steady_sol=Eq(steady_solution_subs,0)

rownania_ruchu=Section('Równania ruchu')
CurrentContainer(rownania_ruchu)
display(ReportText('równania ruchu: '))

display(SympyFormula(eoms_eq))
display(SympyFormula(eq_sol))
display(SympyFormula(eq_steady_sol))
''')

report_generation_libraries=(
'''
doc_final = MechanicalCase('./Out/Document',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)

doc_final.append(intro)
doc_final.append(rownania_ruchu)
doc_final.append(wykres)

doc_final.generate_pdf()
''')

class SimulationReportComponent(ReportComponent):
    
    title="Raportowanie wykresów"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('Sposób tworzenia sekcji przy użyciu metody CurrentContainer:'))

        #sekcja 1|
        display(ObjectCode(section_development_str))


        display(ReportText('Poniżej przedstawiony jest sposób tworzenia równań. W tym przypadku przyrównamy wyrazy (przed przecinkiem) do 0 (po przecinku).'))
        display(ObjectCode(forming_equation_str))



        display(ReportText('poniżej przedstawiony jest sposób tworzenia wykresu dla podanych danych: '))
        eoms=system._eoms[0]
        slownik=system.get_random_parameters()
        eoms_eq=Eq(eoms,0)
        solution_subs = system._ode_system.steady_solution.subs(slownik)
        steady_solution=system._ode_system.steady_solution[0]
        steady_solution_subs=steady_solution.subs(slownik)
        steady_solution_subs
#         eq_sol=Eq(solution_subs,0)
#         eq_steady_sol=Eq(steady_solution_subs,0)
        display(SympyFormula(solution_subs))



        display(ReportText('równania ruchu: '))

        display(SympyFormula(eoms_eq))
        display(SympyFormula(solution_subs))
        display(SympyFormula(steady_solution_subs))




        display(ObjectCode(report_generation_libraries))



        display(ReportText('Tak wygląda ostateczny raport:'))
        display(Picture('./dynpy/utilities/components/guides/images/sekcja3_w.jpg'))
        
        
        
#pandas guide
data_code=(
'''
miesiace_list = ['styczeń','luty','marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']  
srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]  
Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]  
Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]
''')

dict_code=(
'''
data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii         [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia             temperatura':srednie_temp_list}
''')
output_code=(
'''
df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)
df
''')

class PandasTableGenerationComponent(ReportComponent):
    
    title="Tworzenie tabelki"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Stworzenie listy danych:'))

        display(ObjectCode(data_code))


        display(ReportText('Stworzenie słownika:'))

        display(ObjectCode(dict_code))


        display(ReportText('Wywolanie tabelki:'))
        display(ObjectCode(output_code))
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
        display(ObjectCode(winch_pandas_code))


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
        display(Markdown(f'Wyniki symulacji zostały zestawione w tabeli {AutoMarker(table_for_report)}'))
        display(table_for_report.reported(caption='Tabela'))
        sim_plot = table.to_pylatex_tikz().in_figure(caption='Wykres')
        display(Markdown(f'Wyniki symulacji zostały przedstawione graficznie na wykresie {AutoMarker(sim_plot)}'))
        display(sim_plot)
        
        
class PandasMethodsComponent(ReportComponent):
    
    title="Metody biblioteki Pandas"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        
        miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

        srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

        Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

        Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

        długosc_dnia_w_miesiacach_godziny = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

        data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia temperatura':srednie_temp_list}

        df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)

        display(ReportText('***Metoda iloc:***'))
        display(ObjectCode('df.iloc[0:3]'))
        display(df.iloc[0:3])

        # display(Picture('./Image/iloc.jpg', caption = ""))

        display(ReportText('***Metoda loc:***'))
        display(ObjectCode('df.loc["styczeń"]'))
        display(df.loc["styczeń"])
        # display(Picture('./Image/loc.jpg', caption = ""))

        display(ReportText('***Metoda set axis dziala za rowno dla kolumn jak i rzędów. Dla rzędów:***'))
        display(ObjectCode("def.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'],axis = 'index')"))
        display(df.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'], axis = 'index'))

        # display(Picture('./Image/set_axis.jpg', caption = ""))

        display(ReportText('***Dla kolumn:***'))
        display(ObjectCode("df.set_axis(['a','b','c'],axis = 'columns')"))
        display(df.set_axis(['a','b','c'],axis = 'columns'))
        # display(Picture('./Image/set_axis2.jpg', caption = ""))

        display(ReportText('***Metoda rename - Zmienienie pojedynczej kolumny:***'))
        display(ObjectCode('''df.rename(columns = {"Długość dnia w miesiącu":'AAA'}, index = {"styczeń":'A'} )'''))
        display(df.rename(columns = {"Długość dnia w miesiącu":'AAA'}, index = {"styczeń":'A'} ))
        # display(Picture('./Image/rename.jpg', caption = ""))

        display(ReportText('***Metoda applymap:***'))
        display(ObjectCode('     df.applymap(lambda x:x+2)'))
        display(df.applymap(lambda x:x+2))

        # display(ReportText('***Metoda applymap dla kolumny/wiersza:***'))
        # display(Markdown('''
        # Warto pamiętać, że df[NAZWA] => tworzy serię, a df[[Nazwa]] => tworzy nowy frame, który można wykorzystać do użycia .applymap'''))
        # display(df[['Długość dnia w miesiącu']].applymap(lambda x:x+100))
        # pd_h = df[['Długość dnia w miesiącu']].applymap(lambda x:x+100)
        # list1 = pd_h['Długość dnia w miesiącu'].tolist()
        # df['Długość dnia w miesiącu'] = list1

        # display(df)

        # display(Picture('./Image/applymap.jpg', caption = ""))

        display(ReportText('***Slicowanie:***'))
        display(ObjectCode('''df['Długość dnia w miesiącu']'''))
        display(df['Długość dnia w miesiącu'])
        # display(Picture('./Image/slice2.jpg', caption = ""))
        
        
steady_state_code=(
'''
steady_state = SDOFWinchSystem()._ode_system_solution()[0]
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
df = pd.DataFrame(index = Model_numbber, data = system_parameters)
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
    
    title="Symulacje - pozostałe warianty"



    
    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        SDOFWinchSystem = type(system)

        display(ReportText('Wygenerowanie rozwiazania szczegolnego klasy:'))
        display(ObjectCode(steady_state_code))
        steady_state = SDOFWinchSystem.from_random_data()._ode_system.solution.with_ics([0,0])[0].subs(SDOFWinchSystem.m0,10)
        display(steady_state)

        # display(Picture('./Image/steady_sol.jpg', caption = ""))


        display(ReportText('Stworzenie listy i tabelki'))
        display(ObjectCode(list_code))
        Model_number = [1,2,3,4] 
        inertia_list= [1,2,3,4] 
        stiffness_list = [1,2,3,4] 
        damping_list = [1,2,3,4] 
        Excitation_list = [1,2,3,4]

        system_parameters = {'Inertia':inertia_list, 'Stiffness':stiffness_list,'Damping':damping_list,'Excitation':Excitation_list}


        display(ObjectCode(df_code))

        df1 =  pd.DataFrame(index = Model_number, data = system_parameters)
        display(df1)

        # display(Picture('./Image/listatab.jpg', caption = ""))


        display(ReportText('Utworzenie wektora - zestaw wartosci czasu dla naszej symulacji '))
        display(Markdown('''
            function = lambdify(SDOFWinchSystem.ivar, steady_state)
            t_span = np.linspace(0,100,200)'''))

        function = lambdify(SDOFWinchSystem.ivar, steady_state)
        t_span = np.linspace(0,100,200)

        # display(Picture('./Image/wektor1.jpg', caption = ""))


        display(ReportText('Utworzenie słowników z danymi i podstawienie danych w nim zawartych do naszego rownania, tak aby bylo zalezne jedynie od czasu '))

        display(ObjectCode(dict_code))
        # display(Picture('./Image/slowniki_sym.jpg', caption = ""))


        display(ReportText('Utworzenie listy elementów podstawiajac po kolei wartosci czasu '))
        display(ObjectCode(list_elem))
        display(Picture('./dynpy/utilities/components/guides/images/eq2_fun.jpg', caption = ""))


        display(ReportText('Generowanie wykresu rownania szczegolnego '))
        display(ObjectCode(solution_code))
        display(Picture('./dynpy/utilities/components/guides/images/wykres_niefun.jpg', caption = ""))


        display(ReportText('Generowanie wykresu rownania ogolnego '))
        display(ObjectCode(sym_gen_sol))
        display(Picture('./dynpy/utilities/components/guides/images/ogolne_tab.jpg', caption = ""))
        display(Picture('./dynpy/utilities/components/guides/images/ogolne_plot.jpg', caption = ""))

class BasicOperationsComponent(ReportComponent):

    title="Podstawowe obekty i operacje"

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
        s5.index =['pierwszy','drugi','trzeci','czwarty','piąty','szósty']
        x5 = LatexDataFrame(s5)

        wojewodztwa_list = ['mazowieckie','wielkopolskie','lubelskie','warmińsko-mazurskie','zachodniopomorskie','podlaskie','dolnośląskie','pomorskie','łódzkie','kujawsko-pomorskie','podkarpackie','małopolskie','lubuskie','śląskie','świętokrzyskie','opolskie']
        powierzchnia_lista=[35558,29826,25122,24173,22892,20187,19947,18310,18219,17972, 17846,15183,13988,12333,11711,9412]
        ludnosc_lista=[5349114,3475323,2139726,1439675,1710482,1188800,2904207,2307710,2493603,2086210, 2127657,3372618,1018075,4570849,1257179,996011]
        dochody_lista=[4464,3371,3338,3447,3649,3552,3788,4104,3438,3508,3212, 3395,3187,3586,3268,3112]
        wydatki_lista=[787,596,623,597,767,697,742,1023,590,778,574,598,365,631,647,431]
        wojwodztwa_dane={'powierzchnia':powierzchnia_lista,'l.osob':ludnosc_lista,'dochody':dochody_lista,'wydatki':wydatki_lista}
        
        df = pd.DataFrame(index=wojewodztwa_list,data=wojwodztwa_dane)
        x6 = LatexDataFrame(df)

        area = Symbol('A',positive=True)
        df['zaludnienie']=df['l.osob']/df['powierzchnia']
        x7 = LatexDataFrame(df)
        
        df21=df.iloc[1:3,0:3]
        x8 = LatexDataFrame(df21)
        
        df22=df.loc['mazowieckie']
        x9 = LatexDataFrame(df22)
        
        df.loc['Warszawa',:]=[4000,600000,2500,300,120]
        x10 = LatexDataFrame(df)
        
        df30=df.set_axis([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],axis='index')
        x11 = LatexDataFrame(df30)
        
        a,b,c,d,e = symbols('a b c d e',positive = True)

        units_dict = {a:ureg.meter,b:ureg.second, 'lubuskie':ureg.meter}
        LatexDataFrame.set_default_units(units_dict)

        df31=df.set_axis([a,b,c,d,e],axis='columns')
        x12 = LatexDataFrame(df31)
        
        df40=df.rename(columns={'l.osob':'LUDNOŚĆ'},index={'lubuskie':'LUBU...'})
        x21 = LatexDataFrame(df40)
        
        df50=df.dochody.apply(lambda x: x*2)
        x13 = LatexDataFrame(df50)
        
        df60=df.applymap(lambda x: x*2)
        x14 = LatexDataFrame(df60)
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        display(ReportText('**Importowanie biblioteki Pandas**'))
        display(ObjectCode(
            '''
            import pandas as pd
            '''))

        display(ReportText('**Tworzenie serii (analogicznie do kolumn w Excelu)**'))
        display(x1.reported())

        display(ObjectCode(
            '''
            s = pd.Series([2,34,5,-2,8,12]) 
            '''))

        display(ReportText('**Operacje na seriach**s*10 '))
        display(x2.reported())

        display(ReportText('**Funkcja abs() zmiana wartrosci na dodtanie** s.abs()'))
        display(x3.reported())

        display(ReportText(' **Funkcja describe() - podstawowe statystyki** s.describe() '))
        display(x4.reported())

        display(Markdown(
        '''
        **Zmiana nazw indeskow domyslnie zaczynajacych sie od 0**
        s.index =['pierwszy','drugi','trzeci','czwarty','piąty','szósty']
        '''))
        display(x5.reported())
        
        display(ReportText('**DataFrame (tabele) - tworzenie DataFrame na bazie słownika** *tworzenie listy wojewodztw*'))
        display(ObjectCode(
            '''
            wojewodztwa_list = ['mazowieckie','wielkopolskie','lubelskie'
            'warmińsko-mazurskie','zachodniopomorskie','podlaskie'
            'dolnośląskie','pomorskie','łódzkie'
            'kujawsko-pomorskie','podkarpackie','małopolskie'
            'lubuskie','śląskie','świętokrzyskie','opolskie']
            '''))
        
        display(ReportText('*tworzenie listy powierzchni wojewodztw w km^2*'))
        display(ObjectCode(
            '''
            powierzchnia_lista=[35558,29826,25122,24173,22892,20187,19947,18310,18219,17972,17846,15183,13988,12333,11711,9412]
            '''))

        display(ReportText('*tworzenie listy liczby ludności*'))
        display(ObjectCode(
            '''
            ludnosc_lista=[5349114,3475323,2139726,1439675,1710482,1188800,2904207,2307710,2493603,2086210,2127657,3372618,1018075,4570849,1257179,996011]
            '''))
        
        display(ReportText('*tworzenie listy dochodow na osobe*'))
        display(ObjectCode(
            '''
            dochody_lista=[4464,3371,3338,3447,3649,3552,3788,4104,3438,3508,3212,3395,3187,3586,3268,3112]
            '''))

        display(ReportText('*tworzenie listy wydatków majątkowych/inwestycyjnych na osobe*'))
        display(ObjectCode(
            '''
            wydatki_lista=[787,596,623,597,767,697,742,1023,590,778,574,598,365,631,647,431]
            '''))

        display(ReportText('*Tworzenie słownika*'))
        display(ObjectCode(
            '''
           wojwodztwa_dane={'powierzchnia':powierzchnia,'l.osob':ludnosc_lista,'doch.':dochody_lista,'wyd.':wydatki_lista}
            '''))

        display(ReportText('*Wywołanie tabeli*'))
        display(ObjectCode(
            '''
            df = pd.DataFrame(index=wojewodztwa_list,data=wojwodztwa_dane)
            df
            '''))
        display(x6.reported())

        display(ReportText('**Operacje na kolumnach - dodanie kolumny z obliczoną gęstością zaludnienia**'))
        display(ObjectCode(
            '''
            df['zalud.']=df['l.os.']/df['pow.']
            df
            '''))
        display(x7.reported())
        
        
        display(ReportText('**Metoda iloc - za pomoca indeksów wykonujemy operacje na wierszach i kolumnach**'))
        display(ObjectCode(
            '''
            df.iloc[1:3,0:3]
            '''))
        display(x8.reported())
        
        display(ReportText('**Metoda loc**'))
        display(ObjectCode(
            '''
            df.loc['mazowieckie']
            '''))
        display(x9.reported())
        
        display(ReportText('**Dodanie wiersza za pomocą metody loc**'))
        display(ObjectCode(
            '''
            df.loc['Warszawa',:]=[4000,600000,2500,300,120]
            '''))
        display(x10.reported())

        display(ReportText('**Metoda set_axis dla wierszy - operacje na indeksach**'))
        display(ObjectCode(
            '''
            df.set_axis([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],axis='index')
            '''))
        display(x11.reported())

        display(ReportText('**Metoda set_axis dla kolumn**'))
        display(ObjectCode(
            '''
            a,b,c,d,e = symbols('a b c d e',positive = True)
            
            units_dict = {a:ureg.meter,b:ureg.second, 'lubuskie':ureg.meter}
            LatexDataFrame.set_default_units(units_dict)

            LatexDataFrame.formatted(df.set_axis([a,b,c,d,e],axis='columns'))
            '''))
        display(LatexDataFrame.formatted(x12).reported())
        
        display(ReportText('**Metoda rename - zmienia nazwy kolumni wierszy** '))
        display(ObjectCode(
            '''
            df.rename(columns={'l.osob':'LUDNOŚĆ'},index={'lubuskie':'LUBU...'})
            '''))
        display(x21.reported())

        display(ReportText('**Metoda apply - implemetacja funkcji na wybranej kolumnie tabeli**'))
        display(ObjectCode(
            '''
            df.dochody.apply(lambda x: x*2)
            '''))
        display(x13.reported())

        display(ReportText('**Metoda applymap - implemetacja funkcji na calej tabeli**'))
        display(ObjectCode(
            '''
            df.applymap(lambda x: x*2)
            '''))
        display(x14.reported())