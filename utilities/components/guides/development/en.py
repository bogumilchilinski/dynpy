from  ...mechanics import *
from  ...mechanics import ReportComponent as BaseReportComponent


from . import pl

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from .....solvers.linear import *
from .....dynamics import *

months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']

from  ...mechanics import display
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

        from .....solvers.linear import ODESystem
        
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
        

    
class DynamicSystemCompletenessCheckComponent(pl.DynamicSystemCompletenessCheckComponent):
    
    title="Dynamic system completeness check component"

        
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

        
        
class ModellingInPythonGuidelinesComponent(ReportComponent):

    #title="Klasa LibrariesImportComponent i jej zastosowanie"
    title="Schedule of course on Advanced Modelling In Python"

    def append_elements(self):
        #classname = self.reported_object['classname']
        #class_module = self.reported_object['module']
        #class_field = self.reported_object['field']
        #target = self.reported_object['target']
        
        #display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(ReportText('Crucial information about the course on Advanced Modelling In Python'))
        display(ObjectCode(kod_2))
        #display(ReportText('Dzięki temu możemy korzystać z przygotowanych wcześniej klas i modułów'))
        display(ReportText('exemplary entry - fill as needed'))
        
        
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
        
class ModellingInPythonScheduleComponent(ReportComponent):

    #title="Klasa LibrariesImportComponent i jej zastosowanie"
    title="Schedule of course on Advanced Modelling In Python"

    def append_elements(self):
        #classname = self.reported_object['classname']
        #class_module = self.reported_object['module']
        #class_field = self.reported_object['field']
        #target = self.reported_object['target']
        
        #display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(ReportText('Crucial information about the course schedule....'))
        display(ObjectCode(kod_2))
        #display(ReportText('Dzięki temu możemy korzystać z przygotowanych wcześniej klas i modułów'))
        display(ReportText('exemplary entry - fill as needed'))        

class PythonBasicsScheduleComponent(ReportComponent):

    #title="Klasa LibrariesImportComponent i jej zastosowanie"
    title="Schedule of course on Python Basics"

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
        
        #display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(Markdown('''The course consists of following classes:

    -'zajecia_1': 'Wprowadzenie do przedmiotu, Przedstawienie guide z obsługi cocalca i zakładania jupytera, Syntax, typy zmiennych, przydzielenie tematów do prezentacji i rozpisanie issues - IntroToCocalc. (Tematy do prezentacji: Cocalc, Jupyter Notebook, import bibliotek - na przykłądnie sympy, numpy)'
    
    -'zajecia_2': "Podstawy gołego Pythona - podstawowe operacje na zmiennych, funkcje, klasy"
    
    -'zajecia_3': 'Biblioteki typu numpy, sympy itd. (studenkie propozycje przykladow zastoseowani z np prostej symulacji)',
    
    -'zajecia_4': 'Raportowanie zaawansowane (guide z raportowania)',
    
    -'zajecia_5': 'Odesystemy (przykładowe rozwiazanie równania różniczkowego i np. trajektoria)',
    
    -'zajecia_6': 'Pandas + symulacja a obiekcie pandas',
    
    -'zajecia_7': 'Wprowadzenia do modułu DynPy - helpy, guidy, tworzenie nowych systemów itd. gdzie szukać jak używać żeby pokazać co można finalnie stworzyć',
    
    -'zajecia_8': 'Budowanie komponentów i tutaj już jakiś zrobić np do systemu zrobionego wcześniej' ??,
    
    -'zajecia_9': 'Full raport na systemie i komponentach stworzony od zera (tutaj może już z tych trzech mini projekt w postaci zastowania pythona) - ???',
    
    -'zajecia_10': 'PROJEKT',
    
    -'zajecia_11': 'PROJEKT',
    
    -'zajecia_12': 'PROJEKT',
    
    -'zajecia_13': 'PROJEKT',
    
    -'zajecia_14': 'PROJEKT',
    
    -'zajecia_15': 'PROJEKT',

'''))

        display(ReportText('exemplary entry - fill as needed'))        

        
class PythonBasicsGuidelinesComponent(ReportComponent):

    #title="Klasa LibrariesImportComponent i jej zastosowanie"
    title="Main requirenents for Python Basics course"
        
    def append_elements(self):
#         classname = self.reported_object['classname']
#         class_module = self.reported_object['module']
#         class_field = self.reported_object['field']
#         target = self.reported_object['target']
        #display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(ReportText('Before working in jupiter, import libraries'))
        display(ObjectCode(kod_2))
        #display(ReportText('Dzięki temu możemy korzystać z przygotowanych wcześniej klas i modułów'))
        display(ReportText('This allows us to use previously prepared classes and modules'))
        zajecia = {}
        tematy = {}


        zajecia_list = [
            'Wprowadzenie do przedmiotu, cocalc biblioteki importy (na przykłądnie sympy nu,py) Kernel its (guide z obsługi cocalca i zakładania jupytera), raportowanie, przydzielenie prezenracja i rozpisanie issues - IntroToCocalc.',
            'Podstawy gołego pythona(prosty guide/raport z np importów i np min max, pętli.',
            'Biblioteki typu numpy sympy itd.(studenkie propozycje przyllakdow zastoseowani z np prostej symulacji).',
            'Raportowanie zaawansowane (guide z raportowanie ???)',
            'Odesystemy(przykładowe rozwiazanie rra i np trajektoria)',
            'Systemy dynamicczne tabelki(przedstawienie wyników obliczen z 6 i 7 w formie tabeli)',
            'Symulacje(pełny raport z symulacji)',
            'Przedstawienie obsługi dynpow typu helpy guidy its, gdzie szukać jak używać żeby pokazać co można finalnie stworzyć',
            'Budowanie systemów i tutaj już może coś zbudować',
            'Budowanie komponentów i tutaj już jakiś zrobić np do systemu zrobionego wcześniej',
            'Full raport na systemie i komponentach stworzony od zera#(tutaj może już z tych trzech mini projekt w postaci zastowania pythona'

    ]

        tematy_list = [
            'Celem pracy jest przygotowanie środowiska języka programowania Python do implementacji modeli związanych z mechaniką analityczną oraz klasyczną, z naciskiem na klasy związane z metodami małego parametru(metody pertubacyjne).',
            'Celem pracy jest przygotowanie środowiska języka programowania Python do implementacji modeli związanych z mechaniką analityczna oraz klasyczną, z naciskiem na opis różnych technik symulacji.',
            'Celem pracy jest zamodelowanie systemu autonomicznej jazdy w środowisku Python (aktywnego tempomatu) oraz zastosowanie Python’a jako narzędzia programistycznego do obliczeń inżynierskich.',
            'The objective of this paper is to examine the parametric vibration of a gearbox with variable mesh stiffness using the multiple time scale method. The study employs two linear models: a two-degree-of-freedom model for numerical analyses and a simplified one-degree-of-freedom gear model for deriving the analytical solution.',
            'Celem pracy jest rozpoznanie możliwości zastosowania języka programowania Python do analizy i badania układów dynamicznych oraz jako narzędzia dla inżynierów z wykorzystaniem biblioteki DynPy.',
            'Celem pracy jest przygotowanie klas wspomagających raportowanie w środowisku języka programowania Python. Stanowią one uzupełnienie do narzędzia umożliwiającego interaktywne rozwiązywanie zagadnień z zakresu mechaniki teoretycznej, z naciskiem na umożliwienie docelowemu użytkownikowi na uzyskanie szczegółowych informacji o metodach analizy.',
            'Celem pracy jest użycie języka programowania Python celem zbadania i analizy dynamiki ogniwa li-ion zastosowanego w modelu pojazdu elektrycznego.',
            'The aim of the thesis is to create model of theoretical battery in Python programming language, as module to Dynpy library.',
            'The thesis aims to conduct a comprehensive comparative analysis of various Battery Management System (BMS) models to assess their performance metrics using simulation techniques in Python. The study will focus on evaluating critical factors such as energy efficiency, state-of-charge estimation, battery lifespan maximization, and overall system reliability. It will also investigate the adaptability of these models to different types of battery chemistries and configurations. The goal is to identify the most effective BMS model that optimizes battery usage and extends the operational life of battery-powered systems.',
            'Celem pracy jest przygotowanie środowiska języka programowania Python do wspomagania prac projektowych z projektu Podstaw Konstrukcji Maszyn',
            'The aim of the paper is investigation of possiblities of utilisation of Python programming language in complex engineering calculations with an outcome in form of academic specialized text',
            'The  goal  of  the  paper  is  to  analyze  the  parametricvibration  of  a  gearbox  with  variable  mesh  stiffness  using  the  multiple  time  scale  method.  It  wasassumed  that  two  linear  models  will  be  used  for  the  study:  a  model  with  two  degrees  of  freedom(numerical  analyses)  and  a  reduced  gear  model  with  one  degree  of  freedom  (to  determine  the analytical solution).', 
            'Celem pracy jest użycie języka programowania Python do modelowania i analizy dynamiki pojazdów z uzgleniem wytrzymałości zmęczeniowej, zjawiska rezonansu oraz bilansu energetycznego'
    ]

        display(ReportText('''
        Python w zastosowaniach inżynierskich i naukowych

        - Przedmiot został stworzony na potrzeby nowej specjalności „Zaawansowane metody projektowania i rozwoju produktu w inżynierii mechanicznej” na studiach II stopnia (magisterskich). Odpowiada potrzebie opracowywania i zastosowania w praktyce innowacyjnych metod kształcenia opartych na pracy zespołowej i rozwiązywaniu problemów projektowych.

        - Idą było stworzenie i uruchamianie modułowych zajęć dydaktycznych o charakterze interdyscyplinarnym, które będą koncentrowały się na wykonywaniu zadań związanych z modelowaniem i rozwiązywaniem kompleksowych problemów inżynierskich oraz naukowych przy użyciu języka programowania Python.

        - Moduł zajęć dydaktyczny tworzą dwa powiązane przedmioty:
        1. Python w zastosowaniach inżynierskich i naukowych.
        2. Zaawansowane metody komputerowego modelowania maszyn i pojazdów.

        - Kompetencje zdobyte na przedmiocie Python w zastosowaniach inżynierskich i naukowych stanowią podstawę dla bardziej zaawansowanych zastosowań języka programowania Python, które realizowane są na przedmiocie Zaawansowane metody komputerowego modelowania maszyn i pojazdów.

        - Podczas zajęć studenci zapoznają się z metodykami programowania w paradygmacie obiektowym z wykorzystaniem mechanizmów języka Python. Opanowują składnię tego języka programowania w kontekście bogatego ekosystemu bibliotek i technologii, które przyczyniają się do jego dużej i wciąż rosnącej popularności. Jednocześnie zdobywają wiedzę z zakresu matematyki obliczeniowej: algorytmy, obliczenia symboliczne i numeryczne, przetwarzanie oraz analiza i wizualizacja danych.

        - Zadania będą realizowane zespołowo. Taka forma sprzyja rozwijaniu umiejętność pracy grupowej, tworzy przestrzeń do „burzy mózgów” i pozwala na efektywniejsze zdobywanie umiejętności projektowania i programowania obiektowego w kontekście rozważanego zagadnienia. Ułatwi to tworzenie adekwatnej architektury aplikacji z uwagi na zakres wymagań.

        - Tworząc przyjazną przestrzeń własnej pracy zespołowej, studenci nabędą umiejętności, które są przydatne podczas studiów, pracy naukowej, a przede wszystkim na współczesnym rynku pracy. Praktykują stosowanie technik i technologii programistycznych przedstawionych na zajęciach wykonując zadania związane z modelowaniem i rozwiązywaniem kompleksowych problemów.
        '''))

        display(ReportText('''
        Adaptacyjna forma zajęć:
        Jest to pierwsza edycja ww. zajęć, które będą poprowadzone na nowej specjalności studiów II st. magisterskich.
        - Wcześniej przeprowadziliśmy kilka względnie krótkich kursów „wdrożeniowych” dla studentów SiMR skupionych wokół grupy naukowo-dydaktycznej Effective Python, którą wspólnie założyliśmy i prowadzimy od 21.02.2022.
        - Tak więc, sami jesteśmy ciekawi jak to ostatecznie wyjdzie. Drodzy Studenci, wspólnie z nami będziecie kreować obraz zajęć w ramach tych dwóch przedmiotów modułowych. W dużej mierze to od Was, waszych chęci i zaangażowania będzie zależało to ile i jak uda nam się zrobić.
        - W ogólnym ujęciu zajęcia podzielone są na 2 fazy:
        I.	Szkolenie z podstaw programowania obiektowego z użyciem języka Python
        \(9 tygodni\)
        II.	Projekt (zespołowy) – praktyczne rozwiązanie zadanego/wybranego problemu
        \(6 tygodni\)
        '''))

        for i in range(11):
            zajecia[str(i + 1)] = zajecia_list[i]
            tematy[str(i + 1)] = tematy_list[i]

        df = pd.DataFrame(data={'Tematyka zajęć': zajecia, 'Guide': tematy})
        df1 = LatexDataFrame(df)


        display(df1.reported())
        display(Markdown('''
        FAZA I
        - Poszczególne zagadnienia przedmiotowe realizowane według harmonogramu zajęć.
        - Każde zajęcia podsumowuje raport konkretnego zastosowania Pythona.
        - Studenci prowadzą spotkania zajęciowe.
        Spokojnie, bez ciśnień, bez stresu. Chodzi o samokształcenie i swobodną formę zajęć, w ramach których na forum grupy każdy student przedstawia swoim kolegom/koleżankom wybrane zagadnienie podstawowe z harmonogramu, tzn. prezentuje i omawia dany guide, wyjaśnia prosto, jak korzystać. Tak jak to bywa w pracy na spotkaniach grupy projektowej czy programistycznej (sprinty Scrum).
        - Fazę I kończy spotkanie grupy, mające na celu podsumowanie dotychczasowych działań, wytyczenie kierunku dalszych prac, tzn. wydanie zadań do wykonania.
        '''))

        display(ReportText('''A przyświeca temu idea, którą określają trzy znane cytaty:

        “If you can't explain it simply, you don't understand it well enough”.
        (Jeśli nie potrafisz czegoś wyjaśnić prosto, nie rozumiesz tego wystarczająco dobrze.)
        Albert Einstein

        “What we don't understand we can make mean anything”.
        (To, czego nie rozumiemy, może mieć dowolne znaczenie.)
        Chuck Palahniuk
        A chodzi o to, żeby to miało właściwe znaczenie, żeby nie przeinaczasz sensu. ;-)

        “Knowing is not enough; we must apply. Willing is not enough; we must do”.
        (Wiedza nie wystarczy, musimy ją zastosować. Chęć nie wystarczy, musimy działać.
        Johann Wolfgang von Goethe'''))

        display(Markdown('''
        FAZA II
        - Planowanie i śledzenie postępu prac: zadania (issues) wyznaczone na GitHub. To jest hostingowy serwis internetowy przeznaczony do projektów programistycznych wspieranych przez system kontroli wersji Git.
        - Spotkania cotygodniowe poświęcone na weryfikację progresu prac (sprinty Scrum).
        - Wykonanie szablonowych raportów z zadań, według wytycznych.
        - Fazę II kończy spotkanie grupy, poświęcone podsumowaniu efektów prac: rozliczenie raportu i wykonanych zadań (issues).
        '''))
