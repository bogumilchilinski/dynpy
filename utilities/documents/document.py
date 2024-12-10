from pylatex import (Document, Package, Command, NewPage, Tabularx
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )
#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)
from ..report import Markdown, CurrentContainer, ReportText, IPMarkdown, ObjectCode,display
from ..components.mech import en as mech_compIntroDynPyProjectGuidlines
from ..components.guides import en as guide_comp
from ..components.ode import pl as ode_comp
from ..components.ode import en as ode_comp_en
from ..components.guides.reporting import en as reporting_comp
from ..components.guides.github import en as github_comp
from ..components.guides.systems import en as systems_comp
from ..components.guides.development import en as development_comp
from ..components.guides.pandas import en as pandas_comp
from .guides import (UsageOfDynamicSystemsGuide, Guide, EngeneeringDrawingGuide, DevelopmentGuide, IntroToPandasGuide, 
                    BasicsOfODESystemGuide, BasicsOfDynSysImplementationGuide, BasicsOfReportingGuide, ResearchProjectGuidelines,
                    InterimProjectGuidelines,IntroDynPyProjectGuidelines, BasicsOfReportComponentImplementationGuide, 
                    GithubSynchroGuide, IntroToCocalcGuide)
#from sympy import *
import datetime
import shutil
import os

class ReportMethods:

    _reported_object = None
    
    @classmethod
    def base_setup(cls):
        
        preliminary_str=(
"""

#Examplary setup is as follows:

## 
## Imports



    from dynpy.utilities.report import *
    from dynpy.utilities.templates.document import Guide

    doc = Guide('./output/report_name')
    

## CELL 2
## Text reporting
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_text = Section('Section that presents text reporting')
    CurrentContainer(sec_text)
    
    display(ReportText('Exemplary text'*100))
    
    #will not work in some projects, restrict usage 
    display(Markdown('Formatted text'*100)) 
    
    
## CELL 3
## Math 

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    from sympy import Eq, Symbol, symbols
    
    sec_formula = Section('Section that presents formulas reporting')
    CurrentContainer(sec_formula)
    
    display(ReportText('Mathematical formulas are reported with the support of sympy and it\\'s symbols.'))
    
    a,b = symbols('a b')
    display(SympyFormula(Eq(a,b)))


## CELL 4
## Picture
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_picture = Section('Section that presents pictures reporting')
    CurrentContainer(sec_picture)

    display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))



## CELL 5
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    guide_name = './output/report_name' #path for report file 

    doc = Guide(guide_name)
    doc.append(sec_text) # adding certain sections
    doc.append(sec_formula)
    doc.append(sec_picture)
    # Generating file
    doc.generate_pdf(clean_tex=True)
    


"""

)
        
        display(IPMarkdown(preliminary_str))    

        
        preliminary_str=(
"""
#Perform basic setup for document creation.

#This method initializes the document and prepares it for content addition.

#Example:

#To prepare a simple document with text and images:

#Good practice here is to allocate 1 section per 1 cell

## ############### CELL 1 ###########################
## Imports

from dynpy.utilities.report import *
from dynpy.utilities.templates.document import Guide

doc = Guide('./output/report_name')
    

## ############### CELL 2 ###########################
## Text reporting
    
sec_text = Section('Section that presents text reporting')
CurrentContainer(sec_text)

display(ReportText('Exemplary text'*100))

#will not work in some projects, restrict usage 
display(Markdown('Formatted text'*100)) 
    
    
## ############### CELL 3 ###########################
## Math 

from sympy import Eq, Symbol, symbols

sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

display(ReportText('Mathematical formulas are reported with the support of sympy and it\\'s symbols.'))

a,b = symbols('a b')
display(SympyFormula(Eq(a,b)))

## ############### CELL 4 ###########################
## Picture

sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))



## ############### CELL 5 ###########################
## Document

# Creating file
# Be sure *output* folder is in the current directory

guide_name = './output/report_name' #path for report file 

doc = Guide(guide_name)
doc.append(sec_text) # adding certain sections
doc.append(sec_formula)
doc.append(sec_picture)
# Generating file
doc.generate_pdf(clean_tex=True)

""")    
        return ObjectCode(preliminary_str) 
    
    
    
    @property
    def _report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        # mech_comp.SchemeComponent,
        # mech_comp.ExemplaryPictureComponent,
        # mech_comp.KineticEnergyComponent,
        # mech_comp.PotentialEnergyComponent,
        # mech_comp.LagrangianComponent,
        # mech_comp.GoverningEquationComponent,
        # #mech_comp.FundamentalMatrixComponent,
        # mech_comp.GeneralSolutionComponent,
        # #mech_comp.SteadySolutionComponent,
   
        ]
        
        return comp_list
    
    @property
    def default_reported_object(self):
        
        return None

    @property
    def reported_object(self):

        reported_obj=self._reported_object
        
        if reported_obj is None:
            return self.default_reported_object
        else:
            return reported_obj
        

    @reported_object.setter
    def reported_object(self,value):

        self._reported_object=value


    def append_components(self,reported_object=None):

        #self.reported_object = reported_object
        
        doc=self

        for comp in self._report_components:
            doc.append(comp(self.reported_object))

        return None



class CaseTemplate(Document,ReportMethods):

    latex_name = 'document'
    _documentclass = 'report'
    
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=['DynPy Team'],options=['R']),
                  Command('fancyhead', arguments=['Mechanical Vibration, 2021'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
                  #Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
                  #Command('newcommand{\dyplom}', arguments=['Inżynierska']),
                  #Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
                  #Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
                  #Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
                  #Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
                  #Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
                  #Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
                  #Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
                  #Command('newcommand{\\title}', arguments=['Wpisać tytuł pracy dyplomowej po angielsku']),
                  #Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
                  #Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
                  #Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
                  #Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
    ]


# \renewcommand{\headrulewidth}{1pt}
# \renewcommand{\footrulewidth}{1pt}
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 reported_object=None,
                 *,
                 documentclass=None,
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        if documentclass is not None: self._documentclass
        
        self._reported_object = reported_object
        
        super().__init__(
            default_filepath=default_filepath,
            documentclass=self._documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )

        self.append_components()





class ExampleTemplate(Guide):
    pass

        
class BPASTSPaper(Document):
    
    latex_name = 'document'
    packages = [
#                     Package('natbib', options=['numbers']),
                    Package('booktabs'),
                    Package('float'),
                    Package('standalone'),
                    Package('siunitx'),
                    Package('bpasts', options=['accepted']),
                    Package('t1enc'),
                    Package('amsmath'),
                    Package('amssymb'),
                    Package('amsfonts'),
                    Package('graphicx'),
                    Package('flushend'),
                    Package('hyperref',['colorlinks=true', 'allcolors=bpastsblue', NoEscape('pdfborder={0 0 0}')])

    ]
    
    abtitle='Paper for BPASTS'
    abauthor='Authors'
    title='Basic title'

    
    def __init__(self,
                 default_filepath='default_filepath',
                 title=None,
                 *,
                 documentclass='article',
                 document_options=['10pt','twoside','twocolumn','a4paper'], # for submission
                 fontenc=None,
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=False,
                 microtype=False,
                 page_numbers=True,
                 indent=None,
                 geometry_options=['inner=30mm', 'outer=20mm', 'bindingoffset=10mm', 'top=25mm', 'bottom=25mm'],#,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
                 data=None):

        if title is not None: self.title=title

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
        
        self.preamble.append(Command('abtitle',arguments=[self.abtitle]))
        self.preamble.append(NoEscape('\\title{'+f'{self.title}'+'}'))
        self.preamble.append(NoEscape('\\abauthor{'+f'{self.abauthor}'+'}'))
        
        self.preamble.append(NoEscape('%%%% EDITORS SECTION'))
        self.preamble.append(NoEscape('\\vol{XX} \\no{Y} \\year{2024}'))
        self.preamble.append(NoEscape('\\setcounter{page}{1}'))
        self.preamble.append(NoEscape('\\doi{10.24425/bpasts.yyyy.xxxxxx}'))
        self.preamble.append(NoEscape('%%%%%%%%%%%%%%%%%%%%'))
        self.append(Command('maketitle'))
        #self.append(NewPage())
        # tu implementować co tam potrzeba
        
    def authors(self,nameno,affiliation=None,coremail=None,corno=0):
        
        #nameno - dict; of author name (string) and affiliation number (0 to x integer)
        #affiliation - dcit; affiliation number (0 to x integer) and affiliation (string)
        #coremail - str; corresponding author email
        #corno - int; which of the authors in the dict is the corresponding author, count starting from 0
        
        #HOW TO USE?
        
        #Doc1 = BPASTSPaper(default_filepath='./output/method_test',title=NoEscape('''Your title'''))
        #author_names={'Damian Sierociński':1,'Bogumił Chiliński':1, 'Karolina Ziąbrowska':2, 'Jakub Szydłowski':2}
        #author_affiliations={1:'Institute of Machine Design Fundamentals, Faculty of Automotive and Construction Machinery Engineering, Warsaw University of Technology',2:'Faculty of Automotive and Construction Machinery Engineering, Warsaw University of Technology'}
        #correspondence='damian.sierocinski@pw.edu.pl'
        #Doc1.authors(nameno=author_names,affiliation=author_affiliations,coremail=correspondence)
        
        counter=0
        auth_string=''
        addr_string=''
        for name,no in nameno.items():
            auth_string=auth_string+name+'$^{'+f'{no}'+'}$'
            if counter==corno:
                auth_string=auth_string+'\email{'+coremail+'}'
            counter=counter+1
            
            if counter != len(nameno):
                auth_string=auth_string+', '
        
        counter_addr=1
        for no,addr in affiliation.items():
            addr_string=addr_string+'$^'+f'{no}'+'$'+addr
            if counter_addr!=len(affiliation):
                addr_string=addr_string+', '
            counter_addr=counter_addr+1

        self.preamble.append(Command('author',arguments=[NoEscape(auth_string)]))
        self.preamble.append(Command('Address',arguments=[NoEscape(addr_string)]))

    def abstract(self,container=None):
        for num, val in enumerate(container):
            self.preamble.append(Command('Abstract',arguments=[ReportText(val)]))
            
    def keywords(self,keywords=None):
        self.preamble.append(Command('Keywords',arguments=[keywords]))
        
class WutThesis(Document):
    
    latex_name = 'document'
    packages = [
                  Package('geometry',options=['lmargin=30mm', 'rmargin=30mm',  'top=25mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Package('graphicx'),
                  Package('indentfirst'),
                  Package('pdfpages'),
                  Package('amsmath'),

                  Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
                  Command('newcommand{\dyplom}', arguments=['Magisterska']),
                  Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
                  Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
                  Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
                  Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
                  Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
                  Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
                  Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
                  Command('newcommand{\\album}', arguments=['Wpisać numer albumu']),
                  Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
                  Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
                  Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
                  #Command('renewcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
                  Command('graphicspath{{../}}'),
                  Command('frenchspacing'),
                  Command('counterwithin{figure}{section}'),
                  Command('counterwithin{table}{section}'),
                  Command('fancypagestyle{headings}{\\fancyhead{} \\renewcommand{\headrulewidth}{1pt} \\fancyheadoffset{0cm} \\fancyhead[RO]{\\nouppercase{\\leftmark}} \\fancyhead[LE]{\\nouppercase{\\leftmark}} \\fancyfoot{} \\fancyfoot[LE,RO]{\\thepage}}'),
                  Command('fancypagestyle{plain}{\\fancyhf{} \\renewcommand{\\headrulewidth}{0pt} \\fancyfoot[LE,RO]{\\thepage}}'),
                  Command('numberwithin{equation}{section}'),
                  Command('renewcommand{\\familydefault}{\\sfdefault}'),
        #\renewcommand{\familydefault}{\sfdefault}
        
    ]
    
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 title='Basic title',
                 *,
                 documentclass='article',
                 document_options=['a4paper','11pt','twoside'],
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=['inner=30mm', 'outer=20mm', 'bindingoffset=10mm', 'top=25mm', 'bottom=25mm'],#,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
        self.title=title
        #self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
#         self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
#         self.packages.append(Command('newcommand{\praca}', arguments=['Praca dyplomowa']))
#         self.packages.append(Command('newcommand{\dyplom}', arguments=['Magisterska']))
#         self.packages.append(Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']))
#         self.packages.append(Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']))
#         self.packages.append(Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']))
#         self.packages.append(Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']))
#         self.packages.append(Command('newcommand{\promotor}', arguments=['Wpisać promotora']))
#         self.packages.append(Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']))
#         self.packages.append(Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']))
#         self.packages.append(Command('newcommand{\\album}', arguments=['303596']))
#         self.packages.append(Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']))
#         self.packages.append(Command('newcommand{\\rok}', arguments=['Rok składania pracy']))
#         self.packages.append(Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']))
#         #self.packages.append(Command('renewcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']))
        self.packages.append(Command('graphicspath{{../}}'))
#         self.append(Command('maketitle'))
        self.append(NoEscape('%%% New doc'))
        # tu implementować co tam potrzeba
        
    @classmethod
    def base_setup(cls, create_file = False):

        if create_file is True:
            return cls._create_base_setup_env()
        
        preliminary_str=(
"""

Examplary setup is as follows:

#References to guide imports needed to include - to see the list of available guides insert and run the code below:
```{python}
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```

## CELL 1
## Imports

```{python}
#Create file output
#Create file Images
#In file output create bibliography as .bib file (dynpy123)

from dynpy.utilities.report import *
from dynpy.utilities.documents.document import WutThesis
from dynpy.models.odes.linear import SpringMassEps
from sympy import symbols, Eq
from sympy.printing.latex import latex
from IPython.display import display, Math
import tikzplotlib
from dynpy.utilities.adaptable import TimeDataFrame
import pandas as pd


doc = WutThesis('./output/thesis_name')
```
    


## CELL 2
## Thesis introduction

```{python}    
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_intro = Section('Section that presents text reporting')
CurrentContainer(sec_intro)

sub_problem_outline = Subsection('Outline of the problem')
CurrentContainer(sub_problem_outline)
display(ReportText('This subsection provides information about investigated problem. '*100))

sub_obj_assum = Subsection('Objectives and assumptions')
CurrentContainer(sub_obj_assum)
display(ReportText('This subsection provides objectives and assumptions. '*100))

sub_SOT = Subsection('State of the art')
CurrentContainer(sub_SOT)
display(ReportText('This subsection provides state of the art. '*100))

sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology. '*100))
```
    
    
## CELL 3
## Math 

```{python}
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

from sympy import Eq, Symbol, symbols

sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

description = Subsection('Description of dynamic model')
CurrentContainer(description)
display(ReportText('This subsection provides description of the model. '*100))

from dynpy.models.mechanics.pendulum import Pendulum
dyn_sys = Pendulum()
dyn_sys.as_picture()

lagrangian = Subsection('Lagrangian and derivatives')
CurrentContainer(lagrangian)
display(ReportText('This subsection provides calculation of lagrangian and derivatives. '*100))

lagrangian = dyn_sys.L[0]
lagrangian

t=dyn_sys.ivar


for coord in dyn_sys.q:

    display(ReportText(f'Calculations of derivatives for ${vlatex(coord)}$'))


    vel=coord.diff(t)
    diff1 = lagrangian.diff(vel)
    diff2 = diff1.diff(t) 
    diff3 = lagrangian.diff(coord)


    diff1_sym = Symbol(f'\\frac{{\\partial L}}{{\\partial {vlatex(vel)}}}')
    diff2_sym = Symbol(f' \\frac{{d }}{{dt}}  {diff1_sym}'  )
    diff3_sym = Symbol(f'\\frac{{\\partial L}}{{\\partial {vlatex(coord)}}}')


    display(SympyFormula(  Eq(diff1_sym,diff1)  ))    
    display(SympyFormula(  Eq( diff2_sym ,diff2)  ))
    display(SympyFormula(  Eq( diff3_sym ,  diff3)  ))

display(ReportText('Outcomes of governing equations analysis. '*100))
    

equations = Subsection('Equations of motion')
CurrentContainer(equations)


display(ReportText('This subsection provides calculation of equations of motion and solution of the system. '*100))

ds1=dyn_sys.eoms

for eq1 in ds1.as_eq_list():
    display(SympyFormula(eq1.simplify()))
    
ds2=dyn_sys.linearized()._ode_system.general_solution

for eq2 in ds2.as_eq_list():
    display(SympyFormula(eq2.simplify()))
    
ds3=dyn_sys.linearized()._ode_system.steady_solution

for eq3 in ds3.as_eq_list():
    display(SympyFormula(eq3.simplify()))
    
ds4=dyn_sys.linearized().eoms.solution

for eq4 in ds4.as_eq_list():
    display(SympyFormula(eq4.simplify()))


```


## CELL 4
## Picture
    
    
```{python}
    
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))

```


## CELL 5
## Simulation
##
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    
     sec_simulation = Section('Section that contains simulation')
     CurrentContainer(sec_simulation)

     display(ReportText(' Basics of ODESystem based simulations are covered in guide to ODESystem '))
     display(ObjectCode(from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide))
     
     sec_ODESystem = Subsection('ODESystem simulation')
     CurrentContainer(sec_ODESystem)
     
     display(ReportText('Firstly create an ODESystem or import it from dynpy.modes.odes.linear.py'))
     display(ObjectCode(spring = SpringMassEps.from_reference_data()))
     
     display(ReportText('Secondly solve the equation using solution method:'))
     display(ObjectCode(spring_sol = spring.solution))
     
     display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
     spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()
     
     MSM_sim = Subsection('MSM simulation')
     CurrentContainer(MSM_sim)
     
     display(ReportText('Firstly create an MultiTimeScaleSolution system or import it from dynpy.modes.odes.linear.py'))
     display(ObjectCode(spring = SpringMassMSM.from_reference_data()))
     
     display(ReportText('Secondly solve the equation using solution method:'))
     display(ObjectCode(spring_sol = spring.solution))
     
     display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
     spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()))
     
     
     
## CELL 6
## Veryfication
 
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

```
from sympy import *
from pint import UnitRegistry
ureg = UnitRegistry()

import pandas as pd

t = Symbol('t')
power = Symbol('P')
work = Symbol('W')

unit_dict = {
t: ureg.second,
power: ureg.watt,
}

LatexDataFrame.set_default_units(unit_dict)



sec_verification = Section('Section that verificates research results')
CurrentContainer(sec_verification)

display(ReportText('Description of verifications concept. '*200))
     


### SELECT ONE CASE

## 1st CASE
## import of external data if path "./data_from_ml/test_1.csv" exists

#test_data = './data_from_ml/test_1.csv'
#df = pd.read_csv(filename, header = None)

## 2st CASE
df = pd.DataFrame({t:[0,1,2],'a':[2,3,4],'b':[8,7,6]}).set_index(t)

df_cols = df.set_axis([power,work],axis=1)

data_tab = TimeDataFrame(df_cols).to_latex_dataframe().reported(caption='Caption')
graph = TimeDataFrame(df_cols).to_pylatex_tikz().in_figure(caption='Caption')



display(data_tab)

display(ReportText('Obtained data analysis. '*200))

display(graph)

display(ReportText('Description of verifications outcomes. '*200))
```
    
## CELL 7
## Conclusion
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_conclusion = Section('Section that contains final conclusions')
    CurrentContainer(sec_conclusion)
    
    display(ReportText('Conclusions '*200))


## CELL 9
## Symbols description
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_symbols = Section('Section that contains all symbols descriptions')
    CurrentContainer(sec_symbols)
    
    E_K,F = symbols('E_K,F')
    descriptions = {
    E_K:r"Kinetic energy", 
    F:r"Force",  
    }
    syms_dict = descriptions
    DescriptionsRegistry().set_descriptions({**syms_dict})
    DescriptionsRegistry().reset_registry()
    SymbolsDescription.set_default_header('  ')

    display(SymbolsDescription({**syms_dict}))

## CELL 10
## Document

```{python}
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

# Creating file
# Be sure *output* folder is in the current directory

thesis_name = './output/report_name' #path for report file 

doc = WutThesis(thesis_name)
# Bibliography of quotations
### Select one bibligraphy managment system
####### BibLatex

doc.preamble.append(Package('biblatex',arguments=["backend=biber","sorting=none"]))
doc.preamble.append(Command('addbibresource','elementy_bibliagrafia.bib'))
####### Natbib
#doc.preamble.append(Package('natbib')
#doc.preamble.append(Command('bibliographystyle','unsrt'))

# TOC
doc.append(Command('tableofcontents')) #adds TOC

doc.append(sec_intro) # adding certain sections
doc.append(sub_problem_outline)
doc.append(sub_obj_assum)
doc.append(sub_SOT)
doc.append(sub_methodology)
doc.append(sec_formula)
doc.append(sec_picture)
doc.append(sec_simulation)
doc.append(sec_verification)
doc.append(sec_tables)
doc.append(sec_symbols)
doc.append(sec_conclusion)


### BibLatex
#doc.append(Command('printbibliography',arguments=["title={Bibliography}"])) - argument is to improve
doc.append(Command('printbibliography',options=[NoEscape("title={Bibliography}")]))

### Natbib
#doc.append(Command('bibliography',arguments=["references"])) # .bib file as "references"


#FIGURES LIST

doc.append(Command('addcontentsline{toc}{section}{List of figures}'))
doc.append(Command('listoffigures'))
doc.append(Command('pagestyle{plain}'))
doc.append(Command('newpage'))

#TABLES LIST

doc.append(Command('addcontentsline{toc}{section}{List of tables}'))
doc.append(Command('renewcommand{\listtablename}{List of tables}'))
doc.append(Command('listoftables'))
doc.append(Command('pagestyle{plain}'))
doc.append(Command('newpage'))

# Generating file
doc.generate_pdf(clean_tex=True)
```

"""
        )

        display(IPMarkdown(preliminary_str))    


        preliminary_str=(
"""

#Example:

#To prepare a simple document with text and images:

#Good practice here is to allocate 1 section per 1 cell



## CELL 1
## Imports



#Create file output
#In file output create bibliography as .bib file

from dynpy.utilities.report import *
from dynpy.utilities.templates.document import WutThesis

doc = WutThesis('./output/thesis_name')


## CELL 2
## Thesis introduction

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_intro = Section('Section that presents text reporting')
CurrentContainer(sec_intro)

sub_problem_outline = Subsection('Outline of the problem')
CurrentContainer(sub_problem_outline)
display(ReportText('This subsection provides information about investigated problem. '*10))

sub_obj_assum = Subsection('Objectives and assumptions')
CurrentContainer(sub_obj_assum)
display(ReportText('This subsection provides objectives and assumptions. '*10))

sub_SOT = Subsection('State of the art')
CurrentContainer(sub_SOT)
display(ReportText('This subsection provides state of the art. '*10))

sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology. '*10))



## CELL 3
## Math 

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

from sympy import Eq, Symbol, symbols

sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

sec_description = Subsection('Description of dynamic model')
CurrentContainer(sec_description)
display(ReportText('This subsection provides description of the model. '*10))

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys

dyn_sys = DynamicSys()
display(dyn_sys.as_picture())


sec_lagrangian = Subsection('Lagrangian and derivatives')
CurrentContainer(sec_lagrangian)
display(ReportText('This subsection provides calculation of lagrangian and derivatives. '*10))

lagrangian = dyn_sys.L[0]
lagrangian

t=dyn_sys.ivar


for coord in dyn_sys.q:

    display(ReportText(f'Calculations of derivatives for ${vlatex(coord)}$'))


    vel=coord.diff(t)
    diff1 = lagrangian.diff(vel)
    diff2 = diff1.diff(t) 
    diff3 = lagrangian.diff(coord)


    diff1_sym = Symbol(f'\\\\frac{{\\\\partial L}}{{\\\\partial {vlatex(vel)}}}')
    diff2_sym = Symbol(f' \\\\frac{{d }}{{dt}}  {diff1_sym}'  )
    diff3_sym = Symbol(f'\\\\frac{{\\\\partial L}}{{\\\\partial {vlatex(coord)}}}')


    display(SympyFormula(  Eq(diff1_sym,diff1)  ))    
    display(SympyFormula(  Eq( diff2_sym ,diff2)  ))
    display(SympyFormula(  Eq( diff3_sym ,  diff3)  ))

display(ReportText('Outcomes of governing equations analysis. '*10))
    

sec_equations = Subsection('Equations of motion')
CurrentContainer(sec_equations)


display(ReportText('This subsection provides calculation of equations of motion and solution of the system. '*10))

ds1=dyn_sys.eoms

for eq1 in ds1.as_eq_list():
    display(SympyFormula(eq1.simplify()))
    
ds2=dyn_sys.linearized()._ode_system.general_solution

for eq2 in ds2.as_eq_list():
    display(SympyFormula(eq2.simplify()))
    
ds3=dyn_sys.linearized()._ode_system.steady_solution

for eq3 in ds3.as_eq_list():
    display(SympyFormula(eq3.simplify()))
    
ds4=dyn_sys.linearized().eoms.solution

for eq4 in ds4.as_eq_list():
    display(SympyFormula(eq4.simplify()))




## CELL 4
## Picture

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))



## CELL 5
## Simulation

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_simulation = Section('Section that contains simulation')
CurrentContainer(sec_simulation)

display(ReportText('Simulation '*50))

## CELL 6 - Preparing of simulations - results

import numpy as np


from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
dyn_sys = DynamicSys()

param_dict = dyn_sys.get_numerical_parameters()
t_span = np.linspace(0,1,1001)

parameter = dyn_sys.system_parameters()[0] # 0-1 change a plot from const mass to const k
param_dict = {**param_dict,parameter:parameter}


display(param_dict)

na_df = dyn_sys.subs(param_dict).numerical_analysis(parameter=parameter,param_span=[10,20], t_span = t_span).with_ics([0.0,2.0])


na_df_sol_1 = na_df.compute_solution(t_span)
na_df_sol_1.iloc[:,[0,2]]#.plot() #it creates only preview in notebook uncomment if needed


#-----------------------------------------------------------

ic_list = [0.0,0.0] ### Dla ukladu o jednym stopniu swobody
parameter = dyn_sys.system_parameters()[0] # 0-1 change a plot from const mass to const k
param_dict = {**param_dict,parameter:parameter}
na_df = dyn_sys.subs(param_dict).eoms.numerical_analysis(parameter=parameter,param_span=[1,2,3], t_span = t_span)

na_df_sol_2 = na_df.with_ics([2.0,0.0]).compute_solution(t_span)
na_df_sol_2.iloc[:,[0,2]]#.plot() #it creates only preview in notebook uncomment if needed

# CELL 7 Results plotting

sec_plots = Subsection ('Plots')
CurrentContainer (sec_plots) 

plot_obj = na_df_sol_1.droplevel(0,axis=1)
my_wykres=plot_obj.to_pylatex_tikz(subplots=True).in_figure(caption='??')
display(my_wykres)

plot_obj = na_df_sol_2.droplevel(0,axis=1)
my_wykres=plot_obj.to_pylatex_tikz(subplots=True).in_figure(caption='??')
display(my_wykres)




## CELL 8
## Veryfication

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#


from sympy import *
from pint import UnitRegistry
ureg = UnitRegistry()

import pandas as pd

t = Symbol('t')
power = Symbol('P')
work = Symbol('W')

unit_dict = {
t: ureg.second,
power: ureg.watt,
}

LatexDataFrame.set_default_units(unit_dict)



sec_verification = Section('Section that verificates research results')
CurrentContainer(sec_verification)

display(ReportText('Description of verifications concept. '*20))



### SELECT ONE CASE

## 1st CASE
## import of external data if path "./data_from_ml/test_1.csv" exists

#test_data = './data_from_ml/test_1.csv'
#df = pd.read_csv(filename, header = None)

## 2st CASE
df = pd.DataFrame({t:[0,1,2],'a':[2,3,4],'b':[8,7,6]}).set_index(t)

df_cols = df.set_axis([power,work],axis=1)

data_tab = TimeDataFrame(df_cols).to_latex_dataframe().reported(caption='Caption')
graph = TimeDataFrame(df_cols).to_pylatex_tikz().in_figure(caption='Caption')



display(data_tab)

display(ReportText('Obtained data analysis. '*20))

display(graph)

display(ReportText('Description of verifications outcomes. '*20))

## CELL 9
## Conclusion

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_conclusion = Section('Section that contains final conclusions')
CurrentContainer(sec_conclusion)

display(ReportText('Conclusions '*20))


## CELL 10
## Symbols description

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

sec_symbols = Section('Section that contains all symbols descriptions')
CurrentContainer(sec_symbols)

E_K,F = symbols('E_K,F')
descriptions = {
E_K:r"Kinetic energy", 
F:r"Force",
}
syms_dict = descriptions
DescriptionsRegistry().set_descriptions({**syms_dict})
DescriptionsRegistry().reset_registry()
SymbolsDescription.set_default_header('  ')

display(SymbolsDescription({**syms_dict}))


## CELL 11
## Document

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

# Creating file
# Be sure *output* folder is in the current directory

thesis_name = './output/report_name' #path for report file 
doc = WutThesis(thesis_name)


# Bibliography of quotations
### Select one bibligraphy managment system
####### BibLatex

doc.preamble.append(Package('biblatex',["backend=biber","sorting=none"]))
doc.preamble.append(Command('addbibresource','elementy_bibliografia.bib'))
####### Natbib
#doc.preamble.append(Package('natbib')
#doc.preamble.append(Command('bibliographystyle','unsrt'))

# TOC
doc.append(Command('tableofcontents')) #adds TOC

doc.append(sec_intro) # adding certain sections
doc.append(sub_problem_outline)
doc.append(sub_obj_assum)
doc.append(sub_SOT)
doc.append(sub_methodology)
doc.append(sec_formula)
doc.append(sec_picture)

doc.append(sec_description)
doc.append(sec_lagrangian)
doc.append(sec_equations)

doc.append(sec_simulation)
doc.append(sec_plots)
doc.append(sec_verification)
#doc.append(sec_tables)
doc.append(sec_symbols)
doc.append(sec_conclusion)

### BibLatex
#doc.append(Command('printbibliography',arguments=["title={Bibliography}"])) - argument is to improve
doc.append(Command('printbibliography',options=[NoEscape("title={Bibliography}")]))

### Natbib
#doc.append(Command('bibliography',arguments=["references"])) # .bib file as "references"

#FIGURES LIST

doc.append(Command('addcontentsline{toc}{section}{List of figures}'))
doc.append(Command('listoffigures'))
doc.append(Command('pagestyle{plain}'))
doc.append(Command('newpage'))

#TABLES LIST

doc.append(Command('addcontentsline{toc}{section}{List of tables}'))
doc.append(Command('renewcommand{\listtablename}{List of tables}'))
doc.append(Command('listoftables'))
doc.append(Command('pagestyle{plain}'))
doc.append(Command('newpage'))

# Generating file
doc.generate_pdf(clean_tex=True)


"""
        )
        return ObjectCode(preliminary_str)

    Jupyter_file_content = """{
    "cells": [
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "d3000c",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": [
                "from dynpy.utilities.report import *",
                "\\n",
                "from dynpy.utilities.templates.document import WutThesis",
                "\\n",
                "",
                "doc = WutThesis('./output/thesis_name')"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "c3d54b",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": [
                "doc.base_setup()"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "11531d",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": []
        }
    ],
    "metadata": {
        "kernelspec": {
            "argv": [
                "/usr/bin/python3",
                "-m",
                "ipykernel",
                "--HistoryManager.enabled=False",
                "--matplotlib=inline",
                "-c",
                "%config InlineBackend.figure_formats = set(['retina'])import matplotlib; matplotlib.rcParams['figure.figsize'] = (12, 7)",
                "-f",
                "{connection_file}"
            ],
            "display_name": "Python 3 (system-wide)",
            "env": {},
            "language": "python",
            "metadata": {
                "cocalc": {
                    "description": "Python 3 programming language",
                    "priority": 100,
                    "url": "https://www.python.org/"
                }
            },
            "name": "python3",
            "resource_dir": "/ext/jupyter/kernels/python3"
        },
        "language_info": {
            "codemirror_mode": {
                "name": "ipython",
                "version": 3
            },
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10.12"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 4
}"""
    articles_bib = '''@INPROCEEDINGS{brief_history_of_simulations,
  author={Goldsman, David and Nance, Richard E. and Wilson, James R.},
  booktitle={Proceedings of the 2010 Winter Simulation Conference}, 
  title={A brief history of simulation revisited}, 
  year={2010},
  volume={},
  number={},
  pages={567-574},
  doi={10.1109/WSC.2010.5679129}
}
@article{eckhardt1987stan,
  title={Stan Ulam, John von Neumann},
  author={Eckhardt, Roger},
  journal={Los Alamos Science},
  volume={100},
  number={15},
  pages={131},
  year={1987},
  publisher={Los Alamos Scientific Laboratory}
}
@book{lutz2001programming,
  title={Programming python},
  author={Lutz, Mark},
  year={2001},
  publisher={" O'Reilly Media, Inc."}
}
@misc{NumPy, url={https://numpy.org/}, journal={NumPy}}
@misc{pandas, url={https://pandas.pydata.org/}, journal={pandas}}
@misc{RTPW,
  author = {mgr Tomasz Duda},
  title = {Rysunek Techniczny Podstawowe Wiadomości},
  journal = {},
  year = {},
  number = {},
  pages = {10},
  doi = {}
}
'''
    def _create_directory(path):
        '''
Method create directory from given path

Arguments:
path - (str), name of directory that method will create, if subdirectories shall be created, path should look like this: 'directory/subdirectory' 
        '''
        import os

        try:
            os.makedirs(path)
            print(f"Directory '{path}' created")

        except FileExistsError:
            print(f"Directory '{path}' exists.")

        except PermissionError:
            print(f"Permission denied, you don't have sudo permission")

        except Exception as e:
            print(f"An error occurred: {e}")
        
        
    def _create_file(name, content = None, path = None):
        '''
Method create file and write content to it

Arguments:
name - (str), name of file created by method,
content - (str), optional,  string with a content of file, that is going to be writen to it,
path - (str), optional, directory where file should be created
        '''
        import os

        if path is not None:
            if not os.path.exists(path):
                _create_directory(path)

            with open(os.path.join(path, name), 'w') as file: 
                file.write(content)

        else:
            with open(name, 'w') as file:
                file.write(content)

        file.close()
    @classmethod    
    def _create_base_setup_env(cls):
        '''
Method that create Jupyter notebook file with WUT thesis base setup and output directory with .bib file
        '''
        
        Jupyter_file_content = """{
    "cells": [
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "d3000c",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": [
                "from dynpy.utilities.report import *",
                "\\n",
                "from dynpy.utilities.templates.document import WutThesis",
                "\\n",
                "",
                "doc = WutThesis('./output/thesis_name')"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "c3d54b",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": [
                "doc.base_setup()"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": 0,
            "id": "11531d",
            "metadata": {
                "collapsed": false
            },
            "outputs": [],
            "source": []
        }
    ],
    "metadata": {
        "kernelspec": {
            "argv": [
                "/usr/bin/python3",
                "-m",
                "ipykernel",
                "--HistoryManager.enabled=False",
                "--matplotlib=inline",
                "-c",
                "%config InlineBackend.figure_formats = set(['retina'])import matplotlib; matplotlib.rcParams['figure.figsize'] = (12, 7)",
                "-f",
                "{connection_file}"
            ],
            "display_name": "Python 3 (system-wide)",
            "env": {},
            "language": "python",
            "metadata": {
                "cocalc": {
                    "description": "Python 3 programming language",
                    "priority": 100,
                    "url": "https://www.python.org/"
                }
            },
            "name": "python3",
            "resource_dir": "/ext/jupyter/kernels/python3"
        },
        "language_info": {
            "codemirror_mode": {
                "name": "ipython",
                "version": 3
            },
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10.12"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 4
}"""
        articles_bib = '''@INPROCEEDINGS{brief_history_of_simulations,
  author={Goldsman, David and Nance, Richard E. and Wilson, James R.},
  booktitle={Proceedings of the 2010 Winter Simulation Conference}, 
  title={A brief history of simulation revisited}, 
  year={2010},
  volume={},
  number={},
  pages={567-574},
  doi={10.1109/WSC.2010.5679129}
}
@article{eckhardt1987stan,
  title={Stan Ulam, John von Neumann},
  author={Eckhardt, Roger},
  journal={Los Alamos Science},
  volume={100},
  number={15},
  pages={131},
  year={1987},
  publisher={Los Alamos Scientific Laboratory}
}
@book{lutz2001programming,
  title={Programming python},
  author={Lutz, Mark},
  year={2001},
  publisher={" O'Reilly Media, Inc."}
}
@misc{NumPy, url={https://numpy.org/}, journal={NumPy}}
@misc{pandas, url={https://pandas.pydata.org/}, journal={pandas}}
@misc{RTPW,
  author = {mgr Tomasz Duda},
  title = {Rysunek Techniczny Podstawowe Wiadomości},
  journal = {},
  year = {},
  number = {},
  pages = {10},
  doi = {}
}
'''
        cls._create_file('WUT_Thesis_starter.ipynb', Jupyter_file_content)

        os.makedirs('output')
        os.makedirs('tikzplots')

        cls._create_file('articles.bib', articles_bib, 'output')
        
        

class ResearchProjectReport(WutThesis):

    @classmethod
    def base_setup(cls):

        
        preliminary_str=(
"""

Examplary setup is as follows:

## CELL 1
## Imports


    #Create file output
    #Create file Images
    
    #In file output create bibliography as .bib file (dynpy123)

    from dynpy.utilities.report import *
    from dynpy.utilities.templates.document import ResearchProjectReport

    doc = ResearchProjectReport('./output/thesis_name')
    


## CELL 2
## Thesis introduction
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_intro = Section('Section that presents text reporting')
    CurrentContainer(sec_intro)
    
    sub_problem_outline = Subsection('Outline of the problem')
    CurrentContainer(sub_problem_outline)
    display(ReportText('This subsection provides information about investigated problem. '*100))
    
    sub_obj_assum = Subsection('Objectives and assumptions')
    CurrentContainer(sub_obj_assum)
    display(ReportText('This subsection provides objectives and assumptions. '*100))
    
    sub_SOT = Subsection('State of the art')
    CurrentContainer(sub_SOT)
    display(ReportText('This subsection provides state of the art. '*100))
    
    sub_methodology = Subsection('Methodology')
    CurrentContainer(sub_methodology)
    display(ReportText('This subsection provides methodology. '*100))
    
    
    
## CELL 3
## Math 

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    from sympy import Eq, Symbol, symbols
    
    sec_formula = Section('Section that presents formulas reporting')
    CurrentContainer(sec_formula)
    
    display(ReportText('Mathematical formulas are reported with the support of sympy and it\\'s symbols.'))
    
    a,b = symbols('a b')
    display(SympyFormula(Eq(a,b)))


## CELL 4
## Picture
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_picture = Section('Section that presents pictures reporting')
    CurrentContainer(sec_picture)

    display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))




## CELL 5
## Simulation
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_simulation = Section('Section that contains simulation')
    CurrentContainer(sec_simulation)
    
     display(ReportText('Simulation '*200))
    
## CELL 6
## Veryfication
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_verification = Section('Section that verificates research results')
    CurrentContainer(sec_verification)
    
     display(ReportText('Verifications '*200))
    
## CELL 7
## Conclusion
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_conclusion = Section('Section that contains final conclusions')
    CurrentContainer(sec_conclusion)
    
    display(ReportText('Conclusions '*200))

## CELL 8
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    thesis_name = './output/report_name' #path for report file 

    doc = ResearchProjectReport(thesis_name)
    # Bibliography of quotations
    doc.preamble.append(NoEscape(r'\\usepackage[backend=bibtex, sorting=none]{biblatex}'))
    doc.preamble.append(NoEscape(r'\addbibresource{elementy_bibliagrafia.bib}'))
    doc.append(sec_intro) # adding certain sections
    doc.append(sub_problem_outline)
    doc.append(sub_obj_assum)
    doc.append(sub_SOT)
    doc.append(sub_methodology)
    doc.append(sec_formula)
    doc.append(sec_picture)
    doc.append(sec_simulation)
    doc.append(sec_verification)
    doc.append(sec_conclusion)
    # Generating file
    doc.generate_pdf(clean_tex=True)
    
    
    
    

"""
        )

        display(IPMarkdown(preliminary_str))    


        preliminary_str=(
"""

Examplary setup is as follows:

## CELL 1
## Imports



    #Create file output
    #In file output create bibliography as .bib file

    from dynpy.utilities.report import *
    from dynpy.utilities.templates.document import ResearchProjectReport

    doc = ResearchProjectReport('./output/thesis_name')
    

## CELL 2
## Thesis introduction
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_intro = Section('Section that presents text reporting')
    CurrentContainer(sec_intro)
    
    sub_problem_outline = Subsection('Outline of the problem')
    CurrentContainer(sub_problem_outline)
    display(ReportText('This subsection provides information about investigated problem. '*100))
    
    sub_obj_assum = Subsection('Objectives and assumptions')
    CurrentContainer(sub_obj_assum)
    display(ReportText('This subsection provides objectives and assumptions. '*100))
    
    sub_SOT = Subsection('State of the art')
    CurrentContainer(sub_SOT)
    display(ReportText('This subsection provides state of the art. '*100))
    
    sub_methodology = Subsection('Methodology')
    CurrentContainer(sub_methodology)
    display(ReportText('This subsection provides methodology. '*100))
    
    
    
## CELL 3
## Math 

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    from sympy import Eq, Symbol, symbols
    
    sec_formula = Section('Section that presents formulas reporting')
    CurrentContainer(sec_formula)
    
    display(ReportText('Mathematical formulas are reported with the support of sympy and it\\'s symbols.'))
    
    a,b = symbols('a b')
    display(SympyFormula(Eq(a,b)))


## CELL 4
## Picture
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_picture = Section('Section that presents pictures reporting')
    CurrentContainer(sec_picture)

    display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))



## CELL 5
## Simulation
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_simulation = Section('Section that contains simulation')
    CurrentContainer(sec_simulation)
    
     display(ReportText('Simulation '*200))
    
## CELL 6
## Veryfication
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_verification = Section('Section that verificates research results')
    CurrentContainer(sec_verification)
    
     display(ReportText('Verifications '*200))
    
## CELL 7
## Conclusion
 
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_conclusion = Section('Section that contains final conclusions')
    CurrentContainer(sec_conclusion)
    
    display(ReportText('Conclusions '*200))

    
    

## CELL 8
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    thesis_name = './output/report_name' #path for report file 

    # Bibliography of quotations
    doc.preamble.append(NoEscape(r'\\usepackage[backend=bibtex, sorting=none]{biblatex}'))
    doc.preamble.append(NoEscape(r'\addbibresource{elementy_bibliagrafia.bib}'))
    doc = ResearchProjectReport(thesis_name)
    doc.append(sec_intro) # adding certain sections
    doc.append(sub_problem_outline)
    doc.append(sub_obj_assum)
    doc.append(sub_SOT)
    doc.append(sub_methodology)
    doc.append(sec_formula)
    doc.append(sec_picture)
    doc.append(sec_simulation)
    doc.append(sec_verification)
    doc.append(sec_conclusion)
    doc.append(NoEscape(r'\printbibliography[title={Bibliografia}]'))
    # Generating file
    doc.generate_pdf(clean_tex=True)
    
    
    
    

""")
        return ObjectCode(preliminary_str)
    
    
        
        
class StateOfTheArtReport(WutThesis):

    @classmethod
    def base_setup(cls):

        
        preliminary_str=(
"""

Examplary setup is as follows:

## CELL 1
## Imports



    #Create file output
    #In file output create bibliography as .bib file

    from dynpy.utilities.report import *
    from dynpy.utilities.templates.document import ResearchProjectReport

    doc = ResearchProjectReport('./output/thesis_name')
    

## CELL 2
## Thesis introduction
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_intro = Section('Section that presents text reporting')
    CurrentContainer(sec_intro)
    
    sub_problem_outline = Subsection('Outline of the problem')
    CurrentContainer(sub_problem_outline)
    display(ReportText('This subsection provides information about investigated problem. '*100))

    sub_SOT = Subsection('State of the art')
    CurrentContainer(sub_SOT)
    display(ReportText('This subsection provides state of the art. '*100))
    
    sub_conclusions = Subsection('Conclusions')
    CurrentContainer(sub_conclusions)
    display(ReportText('This subsection provides conclusions. '*100))

## CELL 3
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    thesis_name = './output/report_name' #path for report file 

    # Bibliography of quotations
    doc.preamble.append(NoEscape(r'\\usepackage[backend=bibtex, sorting=none]{biblatex}'))
    doc.preamble.append(NoEscape(r'\addbibresource{elementy_bibliagrafia.bib}'))
    doc = ResearchProjectReport(thesis_name)
    doc.append(sec_intro) # adding certain sections
    doc.append(sub_problem_outline)
    doc.append(sub_SOT)
    doc.append(sub_conclusions)
    doc.append(sec_conclusion)
    doc.append(NoEscape(r'\printbibliography[title={Bibliografia}]'))
    # Generating file
    doc.generate_pdf(clean_tex=True)
    
    
    

"""
        )

        display(IPMarkdown(preliminary_str))    


        preliminary_str=(
"""

Examplary setup is as follows:

## CELL 1
## Imports



    #Create file output
    #In file output create bibliography as .bib file

    from dynpy.utilities.report import *
    from dynpy.utilities.templates.document import ResearchProjectReport

    doc = ResearchProjectReport('./output/thesis_name')
    

## CELL 2
## Thesis introduction
    
    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#
    
    sec_intro = Section('Section that presents text reporting')
    CurrentContainer(sec_intro)
    
    sub_problem_outline = Subsection('Outline of the problem')
    CurrentContainer(sub_problem_outline)
    display(ReportText('This subsection provides information about investigated problem. '*100))

    sub_SOT = Subsection('State of the art')
    CurrentContainer(sub_SOT)
    display(ReportText('This subsection provides state of the art. '*100))
    
    sub_conclusions = Subsection('Conclusions')
    CurrentContainer(sub_conclusions)
    display(ReportText('This subsection provides conclusions. '*100))

## CELL 3
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    thesis_name = './output/report_name' #path for report file 

    # Bibliography of quotations
    doc.preamble.append(NoEscape(r'\\usepackage[backend=bibtex, sorting=none]{biblatex}'))
    doc.preamble.append(NoEscape(r'\addbibresource{elementy_bibliagrafia.bib}'))
    doc = ResearchProjectReport(thesis_name)
    doc.append(sec_intro) # adding certain sections
    doc.append(sub_problem_outline)
    doc.append(sub_SOT)
    doc.append(sub_conclusions)
    doc.append(sec_conclusion)
    doc.append(NoEscape(r'\printbibliography[title={Bibliografia}]'))
    # Generating file
    doc.generate_pdf(clean_tex=True)
    
    
    
    

""")
        return ObjectCode(preliminary_str)        
class ThesisTemplate(WutThesis):
    pass
        


    
class MechanicalCase(Guide):
    
    latex_name = 'document'
    packages = [
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=['DynPy Team'],options=['R']),
                  Command('fancyhead', arguments=['Mechanics, 2023'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
        ]

    
    
class CaseStudy(Guide):
    
    latex_name = 'document'
    _documentclass = 'article'
    packages = [
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=['B. Chiliński'],options=['R']),
                  Command('fancyhead', arguments=['Studium przypadku, 2023'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
        ]
    
    
class TechThriveMechanicalCase(Guide):
    
    latex_name = 'document'
    _documentclass = 'article'
    packages = [
                  Package('float'),
                  Package('graphicx'),
                  Command('graphicspath{{../}}'),
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Package('svg'),        
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=[NoEscape(r'\includegraphics[height=0.8cm]{./dynpy/models/images/TT.png}')],options=['C']), #
                  Command('fancyhead', arguments=['Mechanika Ogólna'],options=['L']),
                  Command('fancyhead', arguments=[NoEscape('\\today')],options=['R']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
                  Command('fancyfoot',  arguments=['Copyright TechThrive 2023'],options=['L']),
                  Command('fancyfoot',  arguments=['Powered by DynPy'],options=['R']),
                  Command('graphicspath{{../}}'),
        ]
    
class TechThriveODECase(TechThriveMechanicalCase):
    packages = [
                  Command('fancyhead', arguments=['Równania różniczkowe'],options=['L']),
        ]
    
class TechThriveMVCase(TechThriveMechanicalCase):
    

    
    packages = [
                  Command('fancyhead', arguments=['Drgania mechaniczne'],options=['L']),
        ]


    @property
    def default_reported_object(self):

        from ...models.mechanics import ForcedDampedTrolleysWithSprings,ForcedSpringMassSystem
        
        return ForcedSpringMassSystem()   

    @property
    def _report_components(self):
        
        comp_list=[
#         mech_comp.TitlePageComponent,
#         #mech_comp.SchemeComponent,
#         #mech_comp.ExemplaryPictureComponent,
#         mech_comp.KineticEnergyComponent,
#         mech_comp.PotentialEnergyComponent,
#         mech_comp.LagrangianComponent,
#         mech_comp.GoverningEquationComponent,
#         #mech_comp.FundamentalMatrixComponent,
#         mech_comp.GeneralSolutionComponent,
#         #mech_comp.SteadySolutionComponent,
   
        ]
        
        return comp_list



        
class QuotationTemplate(Document):
    packages = [
                  Package('geometry',options=['left=1in','top=1in','right=1in','bottom=1in']),
                  Package('inputenc',options=['utf8']),
                  Package('hyperref',options=['colorlinks']),
                  Package('graphicx'),
                  Package('tabularx'),
                  Package('multirow'),
                  Package('ragged2e'),
                  Package('hhline'),
                  Package('array'),
                  Package('booktabs'),
                  Package('polski',options=['MeX']),]

    def __init__(self,
                 default_filepath='default_filepath',
                 title='Wycena zlecenia',

                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
#         self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
# #         self.append(Command('maketitle'))
#         self.append(NewPage())
        # tu implementować co tam potrzeba
        
class QuotationTemplate(Document):
    packages = [
                  Package('geometry',options=['left=1cm','top=1cm','right=1cm','bottom=1cm']),
                  Package('inputenc',options=['utf8']),
                  Package('hyperref',options=['colorlinks']),
                  Package('graphicx'),
                  Package('tabularx'),
                  Package('multirow'),
                  Package('ragged2e'),
                  Package('hhline'),
                  Package('array'),
                  Package('booktabs'),
                  Package('polski',options=['MeX']),]

    def __init__(self,

                 default_filepath='default_filepath',
                 title='Wycena zlecenia',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
#         self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
# #         self.append(Command('maketitle'))
#         self.append(NewPage())
        # tu implementować co tam potrzeba
        
class ScheduleTemplate(QuotationTemplate):
    pass
    
class LandscapeScheduleTemplate(ScheduleTemplate):
    def __init__(self,
                 default_filepath='default_filepath',
                 title='Wycena zlecenia',
                 *,
                 documentclass='report',
                 document_options='landscape',
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
        
class BeamerPresentation(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  #Package('titlesec'),
                  #Package('fancyhdr'),
                  #Command('pagestyle', arguments=['fancy']),
                  Command('author', arguments=['Szymon Kozłowski & Bogumił Chiliński']),
                  #Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
                  #Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
                  #Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
                  #Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']), 
                  Command('usetheme', arguments=['Madrid']),
                  Command('graphicspath', arguments=[NoEscape('{../}')])

            ]
    @classmethod
    def base_setup(cls):
        
        
        preliminary_str=(
"""

Examplary setup is as follows:

from dynpy.utilities.documents.document import BeamerPresentation
from dynpy.utilities.report import*


#CELL_1
frame_1=Frame(title='Introduction',options=[])
CurrentContainer(frame_1)
display(ReportText('abc '*100))


#CELL_2
frame_2=Frame(title='Physical model',options=['allowframebreaks'])
CurrentContainer(frame_2)
display(ReportText('abc '*100))


#CELL_3
frame_3=CodeFrame(title='Simulation results',options=[])
CurrentContainer(frame_3)
display(ReportText('abc '*100))

#CELL_4
doc = BeamerPresentation('Example_beamer',title='Exemplary presentation')
doc.append(frame_1)
doc.append(frame_2)
doc.append(frame_3)
doc.generate_pdf(clean_tex=False)
    

"""
)
        return ObjectCode(preliminary_str)

    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 title=None,
                 author=None,
                 initials=None,
                 documentclass='beamer',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='footnotesize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data, 
            )
        self.title = title
        self.author = author
        self.initials = initials
        
        if self.title is not None:
            self.packages.append(Command('title', arguments=[self.title], options=['']))
            
        if self.author is not None:
            if self.initials is not None:
                self.packages.append(Command('author', arguments=[self.author], options=[self.initials]))
            else:
                self.packages.append(Command('author', arguments=[self.author], options=['']))
            
        self.append(Command('frame', arguments=[NoEscape(r'\titlepage')]))
        
        
class BeamerTemplate(BeamerPresentation):

    pass

class DGBeamer(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  #Package('english',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  #Package('titlesec'),
                  #Package('fancyhdr'),
                  #Command('pagestyle', arguments=['fancy']),
                  Command('author', arguments=['Bogumił Chiliński & Krzysztof Twardoch']),
                  #Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
                  #Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
                  #Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
                  #Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']), 
                  Command('usetheme', arguments=['Madrid']),
                  Command('graphicspath', arguments=[NoEscape('{../}')])

            ]

    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 title=None,
                 documentclass='beamer',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data, 
            )
        self.title = title
        
        if self.title is not None:
            self.packages.append(Command('title', arguments=[self.title]))
        self.append(Command('frame', arguments=[NoEscape(r'\titlepage')]))
        
class PosterTemplate(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('siunitx'),
                  #Package('titlesec'),
                  #Package('fancyhdr'),
                  #Command('pagestyle', arguments=['fancy']),
                  Command('author', arguments=['Anna Mackojć & Bogumił Chiliński']),
                  #Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
                  #Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
                  #Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
                  #Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']), 
                  Command('usetheme', arguments=['Simple']),
                  Command('institute', arguments=['Institute of Machine Design Fundamentals, Warsaw University Of Technology']),        
                  Command('graphicspath', arguments=[NoEscape('{../}')])

            ]
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 title=None,
                 documentclass='tikzposter',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data, 
            )
        self.title = title
        title_with_box = NoEscape(r'\parbox{0.9\linewidth}{\centering ' +self.title+ ' }')
        
        if self.title is not None:
            self.packages.append(Command('title', arguments=[title_with_box]))

        #self.append(Command('frame', arguments=[NoEscape(r'\titlepage')]))
        self.append(Command('maketitle'))

  

    
class DynSysOverviewReport(UsageOfDynamicSystemsGuide):

    
    @property
    def default_reported_object(self):

        #from ...models.mechanics.tmac import SDOFWinchSystem
        from ...models.mechanics import ForcedSpringMassSystem as DynamicSystem

        return DynamicSystem()
    
    
    @property
    def _report_components(self):

        comp_list=[
            systems_comp.DynSysOverviewUsageComponent,
            systems_comp.DynamicSystemCallComponent,
            pandas_comp.NumericalAnalysisSimulationComponent,
            pandas_comp.AnalyticalSimulationComponent,
            #systems_comp.SimulationsComponent,
            #reporting_comp.SimulationReportComponent,
            systems_comp.DynSysCodeComponent,
            github_comp.IssuePreparationComponent,
            #systems_comp.DynamicSystemCheckerComponent,

        ]

        return comp_list


    
    
    
class ODESystemOverviewReport(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[
            systems_comp.ODESystemOverviewUsageComponent,
            systems_comp.ODEInitCodeComponent,
            systems_comp.ODEGeneralSolutionComponent,
            systems_comp.ODESteadySolutionComponent,
            systems_comp.ODESystemRepresentationComponent,
            ode_comp_en.ODESystemExponentiationComponent,
            ode_comp_en.PredictionOfSteadySolutionComponent,
            ode_comp_en.HomoPredictionIntroComponent,
            ode_comp_en.MainPredictionComponent,
            ode_comp_en.RootsAnalysisComponent,

        ]

        return comp_list

    @property
    def default_reported_object(self):

        #from ...models.mechanics.tmac import SDOFWinchSystem
        from ...models.odes.linear import LinearFirstOrder
        
        return LinearFirstOrder.from_reference_data()
    

    

class MSMethodOverviewReport(UsageOfDynamicSystemsGuide):
    @property
    def _report_components(self):
        comp_list=[
            ode_comp.ODESystemComponent,
            ode_comp.VariablesComponent,
            ode_comp.ODESystemCodeComponent,
            ode_comp.MSMCalculationsOrderComponent,
            ode_comp.PredictedSolutionComponent,
            ode_comp.ParticularDerivativesComponent,
            ode_comp.ZerothOrderApproximatedEqComponent,
            ode_comp.FirstOrderApproximatedEqComponent,
#             ode_comp.GeneralSolutionComponent,
        ]
        return comp_list

    @property
    def default_reported_object(self):
        from sympy import Symbol, sin, S
        from dynpy.solvers.perturbational import MultiTimeScaleSolution
        from dynpy.models.mechanics.gears import EquivalentSDOFGearModel
        sys = EquivalentSDOFGearModel()
        t=ivar=Symbol('t')
        z = sys.z
        delta = Symbol('delta',positive=True)
        eps=sys.eps
        nonlin_ode = MultiTimeScaleSolution(z.diff(t,2)+z*(1+delta*eps+eps*sin(t)), z, ivar=ivar, omega=S.One, order=3,eps=eps)
        return nonlin_ode
    
class MDPIPaper(Document):
    
    latex_name = 'document'
    packages = []

    def __init__(self,
                 default_filepath='default_filepath',
                 title=None,
                 *,
                 documentclass=None,
                 journal=None,
                 document_options=None, # for submission
                 fontenc=None,
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=False,
                 microtype=False,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
                 data=None):


        if document_options is None:
            document_options=[journal, 'article','submit','pdftex','moreauthors']

        super().__init__(
            default_filepath=default_filepath,
            documentclass=['Definitions/mdpi'],
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
        


        self.preamble.append(Command('firstpage','1'))
        self.preamble.append(Command('makeatletter'))
        self.preamble.append(NoEscape('\\setcounter{page}{\@firstpage}'))
        self.preamble.append(Command('makeatother'))
        self.preamble.append(Command('pubvolume','1'))
        self.preamble.append(Command('issuenum','1'))
        self.preamble.append(Command('articlenumber','0'))
        self.preamble.append(Command('pubyear','2024'))
        self.preamble.append(Command('copyrightyear','2024'))
        self.preamble.append(Command('datereceived',' '))
        self.preamble.append(Command('daterevised',' '))
        self.preamble.append(Command('dateaccepted',' '))
        self.preamble.append(Command('datepublished',' '))
        self.preamble.append(Command('hreflink','https://doi.org/'))
        self.preamble.append(Command('Title',title))
        self.preamble.append(Command('TitleCitation',title))
        self.preamble.append(Command('firstnote','Current address: Affiliation'))
        self.preamble.append(Command('secondnote','These authors contributed equally to this work.'))
        #self.append(NewPage())
        # tu implementować co tam potrzeba
        
        cwd=os.getcwd()
        
        source_path=f'/home/user/Shared files/modules/dynpy/utilities/documents/Definitions'
        path1=f'{cwd}/Definitions'
        path2=f'{cwd}/output/Definitions'
        
        if os.path.exists(path1)==False:
            shutil.copytree(source_path, path1)
        if os.path.exists(path2)==False:
            shutil.copytree(source_path, path2)
            
    def authors(self,orcid_dict,corr_no=1):
        
        ######## should be edited further if more affiliations!
        #for orcid_dict provide: consecutive capital latin alphabet letters : [orcid_number,firstnameandlastname]; {A:[00001,'Damian Sierociński'],B:[000002,'Bogumił Chiliński'],...}
        #corr_no: which of the authors in provided list is the corresponding author; int; count starts from 1
        affiliation_no=1
        author_no=1
        authors_no=len(orcid_dict)
        author_string='\\Author{'
        authornames_string='\\AuthorNames{'
        authorcitation_string='\\AuthorCitation{'
        self.od=orcid_dict
        for key, val in orcid_dict.items():
            display(key)
            self.preamble.append(NoEscape('\\newcommand{\\orcidauthor'+key+'}{'+val[0]+'}'))
            author_string+=val[1]+' $^{'+str(affiliation_no)+',\\dagger,\\ddagger'
            authornames_string+=val[1]
            if author_no == corr_no:
                author_string+=',*'
            author_string+='}$\\orcid'+key+'{}'
            if author_no != authors_no and author_no!=authors_no-1:
                author_string+=', '
                authornames_string+=', '
            elif author_no==authors_no-1:
                author_string+=' and '
                authornames_string+=' and '
                
            author_no+=1

        self.preamble.append(NoEscape(author_string+'}'))
        self.preamble.append(NoEscape(authornames_string+'}'))
        self.preamble.append(NoEscape(authorcitation_string+self._citation_format()+'}'))
        

        
        
    def _citation_format(self):
        display(self.od)
        cit_string=''
        for key,val in self.od.items():
            firstname=val[1].split(' ')[0]
            firstname=firstname[0]+'.'
            lastname=val[1].split(' ')[1]
            display(firstname)
            display(lastname)
            cit_string+=lastname+' '+firstname+'; '
            
        return cit_string[:-2]
    def address(self):
        ######## should be edited further if more affiliations!
        self.preamble.append(NoEscape('\\address{$^{1}$ \\quad Department of Computer Techniques, Institute of Machine Design Fundamentals, Faculty of Automotive and Construction Machinery Engineering, Warsaw University of Technology; 84 Ludwika Narbutta Street, 02-524 Warsaw, Poland}'))
        
    def correspondence(self,email):
        self.preamble.append(NoEscape('\\corres{Correspondence: '+email+'}'))
    def abstract(self,text):
        self.preamble.append(Command('abstract',text))
    def keywords(self,keywords_list):
        
        kwl_string='\\keyword{'
        key_num=0
        for num,val in enumerate(keywords_list):
            kwl_string+=val
            if key_num != len(keywords_list)-1:
                kwl_string+='; '
            key_num+=1
        self.preamble.append(NoEscape(kwl_string+'}'))


#             \AuthorCitation{Lastname, F.; Lastname, F.; Lastname, F.}
        
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        