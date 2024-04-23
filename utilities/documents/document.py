from pylatex import (Document, Package, Command, NewPage, Tabularx
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )

#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)
from ..report import Markdown, CurrentContainer, ReportText, IPMarkdown, ObjectCode,display
from ..components.mech import en as mech_comp
from ..components.guides import en as guide_comp

class ReportMethods:

    _reported_object = None
    
    @classmethod
    def base_setup(cls):
        
        preliminary_str=(
"""

Examplary setup is as follows:

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

    oc = Guide(guide_name)
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

oc = Guide(guide_name)
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



        
class Guide(Document,ReportMethods):

    _documentclass = 'article'
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
                  Command('fancyhead', arguments=['Mechanical vibration, 2023'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
                  ]
        


    
    def __init__(self,
                 default_filepath='default_filepath',
                 title='Basic title',
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

#         label=self.label
        self.title='Mechanical vibration'
        #self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
        #self.packages.append(Command('author', arguments=['DynPy Team']))
        self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
        #self.append(Command('maketitle'))
        self.append(NewPage())
        # tu implementować co tam potrzeba
        self.append_components()


class ExampleTemplate(Guide):
    pass

        
class BPASTSPaper(Document):
    
    latex_name = 'document'
    packages = [
                    Package('natbib', options=['numbers']),
                    Package('booktabs'),
                    Package('float'),
                    Package('standalone'),
                    Package('siunitx'),
                    Package('bpasts', options=['accepted']),

    ]


    def __init__(self,
                 default_filepath='default_filepath',
                 title='Basic title',
                 *,
                 documentclass='article',
                 document_options=['10pt','twoside','twocolumn','a4paper'],
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
        self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
        self.packages.append(Command('author', arguments=['Author']))
        self.packages.append(Command('abauthor', arguments=['Author']))
        self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
        self.append(Command('maketitle'))
        #self.append(NewPage())
        # tu implementować co tam potrzeba
        
        
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
        mech_comp.TitlePageComponent,
        #mech_comp.SchemeComponent,
        #mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        #mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        #mech_comp.SteadySolutionComponent,
   
        ]
        
        return comp_list


class EngeneeringDrawingGuide(Guide):
    
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
                  Command('fancyhead', arguments=['Engeneering Drawing, 2023'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
    ]
    
    
class DevelopmentGuide(Guide):
    
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
                  Command('fancyhead', arguments=['DynPy development guide, 2023'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
        ]
        
        
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
        
class BeamerTemplate(Document):

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
        
        
class IntroToCocalcGuide(Guide):

    @property
    def _report_components(self):
        
        comp_list=[

            guide_comp.CocalcLoginComponent,
            guide_comp.JupyterSetUpComponent,
            guide_comp.CocalcFolderComponent,

            guide_comp.CocalcDynSysListComponent,
            guide_comp.ReportingBasicsComponent,
            

        ]

        return comp_list

class UsageOfDynamicSystemsGuide(Guide):

    @property
    def _report_components(self):

        comp_list=[

            guide_comp.DynamicSystemCallComponent,
            guide_comp.DynamicSystemMethodsUsageComponent,
            guide_comp.SimulationsComponent,
            guide_comp.SimulationReportComponent

        ]

        return comp_list
    
    @property
    def default_reported_object(self):
        
        #from ...models.mechanics.tmac import SDOFWinchSystem
        from ...models.mechanics import ForcedSpringMassSystem as DynamicSystem
        
        return DynamicSystem()

class IntroToPandasGuide(Guide):

    @property
    def _report_components(self):

        comp_list=[

            guide_comp.PandasTableGenerationComponent,
            
            guide_comp.PandasMethodsComponent,
            guide_comp.SimulationsComponent, #Common with *UsageOfDynamicSystemsGuide* class
            guide_comp.DifferentSimulationsComponent,
            guide_comp.BasicOperationsComponent

        ]

        return comp_list

    @property
    def default_reported_object(self):

        #from ...models.mechanics.tmac import SDOFWinchSystem
        from ...models.mechanics import ForcedSpringMassSystem as DynamicSystem

        return DynamicSystem()    

class BasicsOfODESystemGuide(Guide):

    @property
    def _report_components(self):

        comp_list=[

            guide_comp.BasicUsageOfODESystemComponent,
            guide_comp.ODEReportComponent,
            guide_comp.ReportCompUseComponent,
            guide_comp.ProjectileExampleComponent,
            guide_comp.ODESimulationComponent,
            guide_comp.ODENumericalSimulationsComponent

        ]

        return comp_list    
    
class DynSysOverviewReport(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[

            guide_comp.DynamicSystemCallComponent,
            guide_comp.DynamicSystemMethodsUsageComponent,
            guide_comp.SimulationsComponent,
            guide_comp.SimulationReportComponent,
            guide_comp.DynSysCodeComponent


        ]

        return comp_list



        