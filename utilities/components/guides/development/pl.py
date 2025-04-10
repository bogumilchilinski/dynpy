from  ...mechanics import *
from  ...mechanics import ReportComponent as BaseReportComponent

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from .....solvers.linear import *
from .....dynamics import *

from ....report import display as ip_display


from  ...mechanics import display



miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

długosc_dnia_w_miesiacach_godziny = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia temperatura':srednie_temp_list}

df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)


# 
#Zakladanie Jupytera        
# 
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

# class JupyterSetUpComponent(en.JupyterSetUpComponent):
    
#     title="Zakładanie Jupytera"

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
        

    
class DynamicSystemCompletenessCheckComponent(ReportComponent):

    title="Wywoływanie i sprawdzanie wszystkich kluczowych elementów systemu dynamicznego"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        eoms=system._eoms[0]
        eoms

        display(ReportText('Wywołanie systemu w pierwszej kolejności można sprawdzić przez sprawdzenie równania ruchu metodą _eoms.'))

        display(ReportText('Wynik jest następujący:'))
        
        display(Markdown(
            '''
        eoms_eq=Eq(eoms,0)
        eoms_eq
            '''))
        
        eoms_eq=Eq(eoms,0)
        eoms_eq
        
        display(SympyFormula(eoms_eq))
        
        display(ReportText('Następnie wykonuje się analizę metody __init__ klasy:'))
        
        display(ObjectCode(system.__init__))
        
        display(ReportText('Kolejnym etapem jest sprawdzenie, czy schemat oraz zdjęcie reczywiste systemu są zdefiniowane:'))
        
        display(Markdown(
            '''
        system._as_picture()
        system.preview(example = 'detail_real_name')
            '''))
        
        system._as_picture()
        system.preview(example = 'detail_real_name')
        
        display(ReportText('Kolejnym etapem jest sprawdzenie, czy zdefiniowane są słowniki z parametrami:'))
        
        display(Markdown(
            '''
        system.get_random_parameters()
        system.get_numerical_parameters()
            '''))
        
        display(ReportText('Jeśli podczas wywołania pojawia się błąd o braku tej metody lub parametry są inne niż w metodzie __init__ należy zmodyfikować klasę analogicznie do pokazanego rozwiązania:'))
        
        display(Markdown(
            '''
        
        def get_default_data(self):

            m0, k0 = self.m0, self.k0

            default_data_dict = {

                self.m: [S.One * no * m0 /100 for no in range(80, 120)], # percentage
                self.k_r: [S.One * no * k0/100 for no in range(80, 120)], # percentage
                self.k_l: [S.One * no * k0/100 for no in range(70, 110)], # percentage
            }
            return default_data_dict

        def get_numerical_data(self):


            m0, k0 = 100, 1000


            default_data_dict = {
                self.m: [S.One * no * m0 /100 for no in range(80, 120)], # percentage
                self.k_r: [S.One * no * k0/100 for no in range(80, 120)], # percentage
                self.k_l: [S.One * no * k0/100 for no in range(70, 110)], # percentage

            }
            return default_data_dict
            '''))
        
        system.get_random_parameters()
        system.get_numerical_parameters()
        
        display(ReportText('Następnie konieczne jest sprawdzenie opisów symboli oraz jednostek:'))
        
        display(Markdown(
            '''
        system.symbols_description()
        system.unit_dict()
            '''))
        
        display(ReportText('Jeśli podczas wywołania tej metody pojawia się błąd, należy wykonać modyfikację w analogizny sposób do:'))
        
                
        display(Markdown(
            '''
        def symbols_description(self):
            self.sym_desc_dict = {
                self.m: r'system's mass',
                self.k_l: r'left spring stiffness',
                self.k_r: r'left spring stiffness',
            }
            return self.sym_desc_dict

        def unit_dict(self):
        
            from sympy.physics import units
            ureg=UnitRegistry()

            unit_dict = {
                self.m: ureg.kilogram,
                self.k_l: ureg.Newton/ureg.meter,
                self.k_r: ureg.Newton/ureg.meter,
            }
            return unit_dict

            '''))
        
        from sympy.physics import units
        ureg=UnitRegistry()
        
        system.symbols_description()
        system.unit_dict()

        
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


        
from dynpy.utilities.components.mech.en import KineticEnergyComponent, KineticEnergySymPyCodeComponent



kod_2='''
from dynpy.utilities.components.guides.en import ReportComponent
    from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, GuideCode, LatexDataFrame
    from dynpy.utilities.components.guides.en import ReportCompImplementationComponent,ReportCompImplementationIssueComponent
    from sympy import *
    from sympy.physics import units
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
        
class ModuleStructureComponent(ReportComponent):

    title="Introduction to ModuleStructure class"
    def append_elements(self):
        display(ReportText('ModuleStructure is a class that can be used to report dynsys class structure. There are few methods implemented. First method get_classes will return list of classes and modules. It can be called using bellow code:'))
        display(ObjectCode('''import dynpy
from dynpy.utilities.creators import ModuleStructure

ModuleStructure('dynpy.utilities.components').get_classes()
        '''))
        display(ReportText('You can also display this sstructure as tree using bellow method:'))
        display(ObjectCode('''import dynpy
from dynpy.utilities.creators import ModuleStructure

display(ObjectCode(ModuleStructure(dynpy.utilities.components).get_module_tree()))
        '''))
        display(ReportText('This class can also be used for updating init files in submodules. Bellow call will outout ready to use text for init file:'))
        display(ObjectCode('''from dynpy.utilities.creators import ModuleStructure

ModuleStructure.get_init_file_content()
        '''))
        display(ReportText('Another usefull feature might be geting import command for any class used in dynpy'))
        display(ObjectCode('''from dynpy.utilities.creators import ModuleStructure

ModuleStructure.get_import('ObjectCode')
        '''))