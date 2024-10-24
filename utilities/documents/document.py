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

## CELL 1
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
                  Command('newcommand{\\album}', arguments=['297537']),
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
    def base_setup(cls):
        
        
        preliminary_str=(
"""

Examplary setup is as follows:

## CELL 1
## Imports



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
## Document

    #!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
    #!!!       BECAUSE OF NEEDED IMPORTS    !!!#

    # Creating file
    # Be sure *output* folder is in the current directory

    thesis_name = './output/report_name' #path for report file 

    doc = WutThesis(thesis_name)
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
            systems_comp.DynamicSystemCheckerComponent,

        ]

        return comp_list


    
    
    
class ODESystemOverviewReport(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[
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
            ode_comp.ZerothOrderApproximatedEqComponent,
            ode_comp.FirstOrderApproximatedEqComponent,
        ]
        return comp_list

    @property
    def default_reported_object(self):
        from sympy import Symbol, sin, S
        from dynpy.solvers.nonlinear import MultiTimeScaleSolution
        from dynpy.models.mechanics.gears import EquivalentSDOFGearModel
        sys = EquivalentSDOFGearModel()
        t=ivar=Symbol('t')
        z = sys.z
        delta = Symbol('delta',positive=True)
        eps=sys.eps
        nonlin_ode = MultiTimeScaleSolution(z.diff(t,2)+z*(1+delta*eps+eps*sin(t)), z, ivar=ivar, omega=S.One, order=3,eps=eps)
        return nonlin_ode