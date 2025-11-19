# from .guides import (UsageOfDynamicSystemsGuide, Guide, EngeneeringDrawingGuide, DevelopmentGuide, IntroToPandasGuide,
#                     BasicsOfODESystemGuide, BasicsOfDynSysImplementationGuide, BasicsOfReportingGuide, ResearchProjectGuidelines,
#                     InterimProjectGuidelines,IntroDynPyProjectGuidelines, BasicsOfReportComponentImplementationGuide,
#                     GithubSynchroGuide, IntroToCocalcGuide)
# from sympy import *
import datetime
import os
import shutil
import csv
from typing import List, Optional, Union
import pandas as pd

from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    NewPage,
    Package,
    Tabularx,
)
from pylatex.base_classes import Environment

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape  # italic,

from ..report import (
    CurrentContainer,
    IPMarkdown,
    Markdown,
    ObjectCode,
    ReportText,
    display,
)

imports_str = """
#Create file output
#Create file Images
#In file output create bibliography as .bib file (biblio.bib)

'''# File content begin


@misc{DynPi,
author={GitHub},
title="bogumilchilinski/dynpy",
url={https://github.com/bogumilchilinski/dynpy}",
note="Accessed:2024-05-04"
}


@book{lutz2001programming,
  title={Programming python},
  author={Lutz, Mark},
  year={2001},
  publisher={" O'Reilly Media, Inc."}
}
@misc{NumPy, url={https://numpy.org/}, journal={NumPy}}
@misc{pandas, url={https://pandas.pydata.org/}, journal={pandas}}
'''# File content end

from dynpy.utilities.report import *
from dynpy.utilities.documents.document import WutThesis
#from dynpy.models.odes.linear import SpringMassDamperMSM
from sympy import symbols, Eq
from sympy.printing.latex import latex
from IPython.display import display, Math

from dynpy.utilities.adaptable import TimeDataFrame
import pandas as pd

global_width=('8cm')#setting picture default size
global_height=('6cm')
doc = WutThesis('./output/thesis_name')
#doc.preamble.append(NoEscape(r'\\usepackage[MeX]{polski}')) #to set polish as main language"""

thesis_introduction_str = """
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
display(ReportText('This subsection provides state of the art with multiple reference for \\cite{DynPi}. '*50))
display(Markdown('Remeber about `biblio.bib` file with referencens. Place it in working directory and subfolder with output files. '*50))


sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology. '*100))"""

thesis_introduction_str_research = """
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

sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology / model description / parameters. '*100))"""

math_str = """
from sympy import Eq, Symbol, symbols
from dynpy.utilities.report import *


sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

sec_description = Subsection('Description of dynamic model')
CurrentContainer(sec_description)
display(ReportText('This subsection provides description of the model. '*10))

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys

dyn_sys = DynamicSys()
display(dyn_sys.as_picture())

display(ReportText('Summary of the subsection highlighting its major achievements. '*10))

sub_model_code = Subsection('Virtual model of the object')
CurrentContainer(sub_model_code)

from dynpy.utilities.report import ObjectCode

display(ReportText('This subsection provides code of the model. '*10))
display(ObjectCode(DynamicSys))
display(ReportText('Summary of the subsection highlighting its major achievements. '*10))


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
ds4_eqns = ds4.as_eq_list()

for eq4 in ds4_eqns:
    display(SympyFormula(eq4.simplify()))


display(ReportText(f'Outcomes of governing equations {AutoMarker(ds4_eqns[0])} {AutoMarker(ds4_eqns[1])}) analysis.'))

sec_math_desc = Subsection('Section that contains all symbols descriptions')
CurrentContainer(sec_math_desc)



E_K,F = symbols('E_K,F')
descriptions = {
E_K:r"Kinetic energy",
F:r"Force",
}
syms_dict = descriptions
DescriptionsRegistry().set_descriptions({**syms_dict})
DescriptionsRegistry().reset_registry()
SymbolsDescription.set_default_header('where:')

display(ReportText('Symbols description is as follows: '))


display(SymbolsDescription(expr=E_K*F))"""

picture_str = """
sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture',width=global_width,height=global_height))"""

simulaion_str = """
from dynpy.utilities.report import *
from sympy import *


sec_simulation = Section('Section that contains simulation')
CurrentContainer(sec_simulation)

# Basics of ODESystem based simulations are covered in guide to ODESystem, use the followin call
# from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
# BasicsOfODESystemGuide()

sec_ODESystem = Subsection('ODESystem simulation')
CurrentContainer(sec_ODESystem)

display(ReportText('Firstly create an ODESystem or import it from dynpy.modes.odes.linear.py'))
display(ObjectCode('spring = SpringMassEps.from_reference_data()''))

display(ReportText('Secondly solve the equation using solution method:'))
display(ObjectCode('spring_sol = spring.solution'))

display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()

# MSM_sim = Subsection('MSM simulation')
# CurrentContainer(MSM_sim)

# display(ReportText('Firstly create an MultiTimeScaleSolution system or import it from dynpy.modes.odes.linear.py'))
# display(ObjectCode('spring = SpringMassMSM.from_reference_data()'))

# display(ReportText('Secondly solve the equation using solution method:'))
# display(ObjectCode('spring_sol = spring.solution'))

# display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
# spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()"""

veryfication_str = """
from sympy import *
from dynpy.utilities.adaptable import *
from dynpy.utilities.report import CurrentContainer, ReportText
from sympy import *
from sympy.physics import units




from dynpy.utilities.adaptable import *
from dynpy.utilities.report import CurrentContainer, ReportText
import pandas as pd

t = Symbol('t')
power = Symbol('P')
work = Symbol('W')

unit_dict = {
t: units.second,
power: units.watt,
work: units.joule,
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

display(ReportText('Description of verifications outcomes. '*200))"""

conclusion_str = """
sec_conclusion = Section('Section that contains final conclusions')
CurrentContainer(sec_conclusion)

display(ReportText('Conclusions '*200))"""

symbols_description_str = """
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

display(SymbolsDescription({**syms_dict}))"""

document_str = """
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

## Units
doc.preamble.append(Package('siunitx'))

# TOC

doc.append(Command('includepdf{./Images/Front_Page.pdf}')) #includes front page##path shouldn't be changed, it must look like ./Images/Your_file.pdf
doc.append(Command('cleardoublepage'))
doc.append(Command('includepdf{./Images/oswiadczenie_autora_pracy.pdf}'))#path shouldn't be changed, it must look like ./Images/Your_file.pdf
doc.append(Command('pagestyle{plain}'))
doc.append(Command('cleardoublepage'))
doc.append(Command('includepdf{./output/oswiadczenie_biblioteczne.pdf}'))#path shouldn't be changed, it must look like ./Images/Your_file.pdf
doc.append(Command('pagestyle{plain}'))
doc.append(Command('cleardoublepage'))

doc.append(Command('tableofcontents')) #adds TOC

doc.append(sec_intro) # adding certain sections
doc.append(sub_problem_outline)
doc.append(sub_obj_assum)
doc.append(sub_SOT)
doc.append(sub_methodology)
doc.append(sec_formula)

doc.append(sub_model_code)
doc.append(sec_math_desc)
doc.append(sec_simulation)

doc.append(sec_picture)
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
doc.generate_pdf(clean_tex=True)"""


class ReportMethods:

    _reported_object = None

    @classmethod
    def base_setup(cls):

        preliminary_str = """

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

        display(IPMarkdown(preliminary_str))

        preliminary_str = """
#Code below is the exact same, as the one presented in the beginning of the guide. Duplicated for easier copying and pasting.

#Method presented below initializes a simple document and prepares it for content addition.

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

"""
        return ObjectCode(preliminary_str)

    @property
    def _report_components(self):

        from ..components.guides import en as guide_comp

        comp_list = [
            # mech_comp.TitlePageComponent,
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

        reported_obj = self._reported_object

        if reported_obj is None:
            return self.default_reported_object
        else:
            return reported_obj

    @reported_object.setter
    def reported_object(self, value):

        self._reported_object = value

    def append_components(self, reported_object=None):

        # self.reported_object = reported_object

        doc = self

        for comp in self._report_components:
            doc.append(comp(self.reported_object))

        return None


class Guide(Document, ReportMethods):
    """
    A class to generate a structured LaTeX document with pre-configured packages and settings.

    Inherits from `Document` and provides additional methods to simplify LaTeX report creation.

    Attributes:
        _documentclass (str): The LaTeX document class (default is 'article').
        latex_name (str): Name of the document.
        packages (List[Union[Package, Command]]): List of LaTeX packages and commands to include in the document.
        _reported_object (Optional[object]): The object being reported on, if any.

    Exemplary Usage:
        >>> from dynpy.utilities.report import *
        >>> from dynpy.utilities.documents.guides import Guide

        >>> doc = Guide('./output/sample_report', title="Sample Report")

        >>> section = Section('Exemplary section name')
        >>> CurrentContainer(section)

        >>> display(Markdown(''' Exemplary Markdown text in the section '''))
        >>> display(ReportText(' Exemplary text appended into section '))

        >>> doc.append(section)

        >>> doc.generate_pdf(clean_tex=True)



    Example of Customization of the document properties:
        >>> from dynpy.utilities.report import *
        >>> from dynpy.utilities.documents.guides import Guide

        >>> custom_geometry = ['lmargin=20mm', 'rmargin=20mm', 'top=25mm', 'bmargin=25mm']

        >>> doc = Guide(
        >>>     default_filepath='./output/custom_report',
        >>>     title='Custom Report',
        >>>     geometry_options=custom_geometry
        >>> )

        >>> section = Section('Exemplary Custom Section Name')

        >>> doc.append(section)

        >>> doc.generate_pdf(clean_tex=True)
    """

    _documentclass = "article"
    latex_name = "document"
    packages = [
        Package(
            "geometry",
            options=[
                "lmargin=25mm",
                "rmargin=25mm",
                "top=30mm",
                "bmargin=25mm",
                "headheight=50mm",
            ],
        ),
        #Package("microtype"),
        Package("authoraftertitle"),
        Package("polski", options=["MeX"]),
        # Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
        Package("listings"),
        Package("titlesec"),
        Package("fancyhdr"),
        Command("pagestyle", arguments=["fancy"]),
        Command("fancyhf", arguments=[""]),
        Command("fancyhead", arguments=["DynPy Team"], options=["R"]),
        Command("fancyhead", arguments=["Mechanical vibration, 2023"], options=["L"]),
        Command("fancyfoot", arguments=[NoEscape("\\thepage")], options=["C"]),
    ]

    def __init__(
        self,
        default_filepath: str = "default_filepath",
        title: str = "Basic title",
        reported_object: Optional[object] = None,
        *,
        documentclass: Optional[str] = None,
        document_options: Optional[List[str]] = None,
        fontenc: str = "T1",
        inputenc: str = "utf8",
        font_size: str = "normalsize",
        lmodern: bool = False,
        textcomp: bool = True,
        microtype: bool = True,
        page_numbers: bool = True,
        indent: Optional[Union[str, int]] = None,
        geometry_options: Optional[List[str]] = None,
        data: Optional[dict] = None,
    ):
        """
        Initialize the Guide class with optional customization.

        Args:
            default_filepath (str): Path for the generated document.
            title (str): Title of the document.
            reported_object (Optional[object]): Object being reported on, if any.
            documentclass (Optional[str]): LaTeX document class (e.g., 'article').
            document_options (Optional[List[str]]): Options for the document class.
            fontenc (str): Font encoding (default 'T1').
            inputenc (str): Input encoding (default 'utf8').
            font_size (str): Font size (default 'normalsize').
            lmodern (bool): Whether to use Latin Modern fonts.
            textcomp (bool): Whether to use the `textcomp` package.
            microtype (bool): Whether to use the `microtype` package.
            page_numbers (bool): Whether to include page numbers.
            indent (Optional[Union[str, int]]): Indentation settings.
            geometry_options (Optional[List[str]]): Geometry options for page layout.
            data (Optional[dict]): Additional data for customization.
        """

        if documentclass is not None:
            self._documentclass

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
        self.title = "Mechanical vibration"
        # self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
        # self.packages.append(Command('author', arguments=['DynPy Team']))
        self.packages.append(Command("date", arguments=[NoEscape("\\today")]))
        # self.append(Command('maketitle'))
        self.append(NewPage())
        # tu implementować co tam potrzeba
        self.append_components()


class CaseTemplate(Document, ReportMethods):

    latex_name = "document"
    _documentclass = "report"

    packages = [
        Package("microtype"),
        Package("polski", options=["MeX"]),
        Package(
            "geometry",
            options=[
                "lmargin=25mm",
                "rmargin=25mm",
                "top=30mm",
                "bmargin=25mm",
                "headheight=50mm",
            ],
        ),
        Package("listings"),
        Package("titlesec"),
        Package("fancyhdr"),
        Command("pagestyle", arguments=["fancy"]),
        Command("fancyhf", arguments=[""]),
        Command("fancyhead", arguments=["DynPy Team"], options=["R"]),
        Command("fancyhead", arguments=["Mechanical Vibration, 2021"], options=["L"]),
        Command("fancyfoot", arguments=[NoEscape("\\thepage")], options=["C"]),
        # Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
        # Command('newcommand{\dyplom}', arguments=['Inżynierska']),
        # Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
        # Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
        # Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
        # Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
        # Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
        # Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
        # Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
        # Command('newcommand{\\title}', arguments=['Wpisać tytuł pracy dyplomowej po angielsku']),
        # Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
        # Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
        # Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
        # Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
    ]

    # \renewcommand{\headrulewidth}{1pt}
    # \renewcommand{\footrulewidth}{1pt}

    def __init__(
        self,
        default_filepath="default_filepath",
        reported_object=None,
        *,
        documentclass=None,
        document_options=None,
        fontenc="T1",
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=True,
        microtype=True,
        page_numbers=True,
        indent=None,
        geometry_options=None,  # ['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
        data=None,
    ):

        if documentclass is not None:
            self._documentclass

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

    latex_name = "document"
    packages = [
        #                     Package('natbib', options=['numbers']),
        Package("booktabs"),
        Package("float"),
        Package("standalone"),
        Package("siunitx"),
        #                     Package('bpasts', options=['accepted']),
        Package("bpasts"),
        Package("t1enc"),
        Package("amsmath"),
        Package("amssymb"),
        Package("amsfonts"),
        Package("graphicx"),
        Package("flushend"),
        Package("textcomp"),
        Package("xcolor"),
        Package(
            "hyperref",
            ["colorlinks=true", "allcolors=bpastsblue", NoEscape("pdfborder={0 0 0}")],
        ),
    ]

    abtitle = "Paper for BPASTS"
    abauthor = "Authors"
    title = "Basic title"

    def __init__(
        self,
        default_filepath="default_filepath",
        title=None,
        *,
        documentclass="article",
        document_options=["10pt", "twoside", "twocolumn", "a4paper"],  # for submission
        fontenc=None,
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=False,
        microtype=False,
        page_numbers=True,
        indent=None,
        geometry_options=[
            "inner=30mm",
            "outer=20mm",
            "bindingoffset=10mm",
            "top=25mm",
            "bottom=25mm",
        ],  # ,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
        data=None,
    ):

        if title is not None:
            self.title = title

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

        self.preamble.append(Command("abtitle", arguments=[self.abtitle]))
        self.preamble.append(NoEscape("\\title{" + f"{self.title}" + "}"))
        self.preamble.append(NoEscape("\\abauthor{" + f"{self.abauthor}" + "}"))

        self.preamble.append(NoEscape("%%%% EDITORS SECTION"))
        self.preamble.append(NoEscape("\\vol{XX} \\no{Y} \\year{2024}"))
        self.preamble.append(NoEscape("\\setcounter{page}{1}"))
        self.preamble.append(NoEscape("\\doi{10.24425/bpasts.yyyy.xxxxxx}"))
        self.preamble.append(NoEscape("%%%%%%%%%%%%%%%%%%%%"))
        self.append(Command("maketitle"))
        # self.append(NewPage())
        # tu implementować co tam potrzeba

    def cwd_setup(self):
        cwd = os.getcwd()

        source_path = (
            f"/home/user/Shared files/modules/dynpy/utilities/documents/Definitions"
        )
        path1 = f"{cwd}/bpast.sty"
        path2 = f"{cwd}/output/bpast.sty"

        if os.path.exists(path1) == False:
            shutil.copytree(source_path, path1)
        if os.path.exists(path2) == False:
            shutil.copytree(source_path, path2)

    def authors(self, nameno, affiliation=None, coremail=None, corno=0):

        # nameno - dict; of author name (string) and affiliation number (0 to x integer)
        # affiliation - dcit; affiliation number (0 to x integer) and affiliation (string)
        # coremail - str; corresponding author email
        # corno - int; which of the authors in the dict is the corresponding author, count starting from 0

        # HOW TO USE?

        # Doc1 = BPASTSPaper(default_filepath='./output/method_test',title=NoEscape('''Your title'''))
        # author_names={'Damian Sierociński':1,'Bogumił Chiliński':1, 'Karolina Ziąbrowska':2, 'Jakub Szydłowski':2}
        # author_affiliations={1:'Institute of Machine Design Fundamentals, Faculty of Automotive and Construction Machinery Engineering, Warsaw University of Technology',2:'Faculty of Automotive and Construction Machinery Engineering, Warsaw University of Technology'}
        # correspondence='damian.sierocinski@pw.edu.pl'
        # Doc1.authors(nameno=author_names,affiliation=author_affiliations,coremail=correspondence)

        counter = 0
        auth_string = ""
        addr_string = ""
        for name, no in nameno.items():
            auth_string = auth_string + name + "$^{" + f"{no}" + "}$"
            if counter == corno:
                auth_string = auth_string + "\email{" + coremail + "}"
            counter = counter + 1

            if counter != len(nameno):
                auth_string = auth_string + ", "

        counter_addr = 1
        for no, addr in affiliation.items():
            addr_string = addr_string + "$^" + f"{no}" + "$" + addr
            if counter_addr != len(affiliation):
                addr_string = addr_string + ", "
            counter_addr = counter_addr + 1

        self.preamble.append(Command("author", arguments=[NoEscape(auth_string)]))
        self.preamble.append(Command("Address", arguments=[NoEscape(addr_string)]))

    def abstract(self, container=None):
        for num, val in enumerate(container):
            self.preamble.append(Command("Abstract", arguments=[ReportText(val)]))

    def keywords(self, keywords=None):
        self.preamble.append(Command("Keywords", arguments=[keywords]))

    @classmethod
    def base_setup(cls, create_file: bool = False):
        """
        Sets up the base template for the BPASTSPaper document.

        This method provides an easy way to initialize a thesis document with the necessary sections, structure, and configuration.

        Args:
            create_file (bool): If True, creates a base setup environment including necessary directories and files.

        Returns:
            Optional[str]: Returns setup content if `create_file` is False, otherwise creates the environment.

        Usage Example:
            >>> from dynpy.utilities.documents.document import BPASTSPaper
            >>> BPASTSPaper.base_setup()

        When to use:
            - Use this method when starting a new BPASTSPaper document.
            - If `create_file=True`, it will generate all required directories and files for a structured thesis.
            - If `create_file=False`, it returns a string with setup instructions that can be used manually.
        """

        if create_file is True:
            return cls._create_base_setup_env()

        # doc = WutThesis('./output/thesis_name')
        # doc.preamble.append(NoEscape(r'\\usepackage[MeX]{polski}')) #to set polish as main language"""

        thesis_introduction_str = """#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
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
display(ReportText('This subsection provides state of the art \\cite{pandas}. '*100))

sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology. '*100))"""

        math_str = '''from sympy import Eq, Symbol, symbols
from dynpy.utilities.report import *


sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

sec_description = Subsection('Description of dynamic model')
CurrentContainer(sec_description)
display(ReportText('This subsection provides description of the model. '*10))

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys

dyn_sys = DynamicSys()
display(dyn_sys.as_picture())

display(ReportText('Summary of the subsection highlighting its major achievements. '*10))

sub_model_code = Subsection('Virtual model of the object')
CurrentContainer(sub_model_code)

from dynpy.utilities.report import ObjectCode

display(ReportText('This subsection provides code of the model. '*10))
display(ObjectCode(DynamicSys))
display(ReportText('Summary of the subsection highlighting its major achievements. '*10))


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
ds1_eqns=ds1.as_eq_list()

for eq1 in ds1.as_eq_list():
    display(SympyFormula(eq1.simplify()))

if ds1.is_solvable():

    ds2=dyn_sys.linearized()._ode_system.general_solution

    for eq2 in ds2.as_eq_list():
        display(SympyFormula(eq2.simplify()))

    ds3=dyn_sys.linearized()._ode_system.steady_solution

    for eq3 in ds3.as_eq_list():
        display(SympyFormula(eq3.simplify()))

    ds4=dyn_sys.linearized().eoms.solution
    ds4_eqns = ds4.as_eq_list()

    for eq4 in ds4_eqns:
        display(SympyFormula(eq4.simplify()))


display(ReportText(f'Outcomes of governing equations {AutoMarker(ds1_eqns[0])} {AutoMarker(ds1_eqns[-1])}) analysis.'))

sec_math_desc = Subsection('Section that contains all symbols descriptions')
CurrentContainer(sec_math_desc)



E_K,F = symbols('E_K,F')
descriptions = {
E_K:r"Kinetic energy",
F:r"Force",
}
syms_dict = descriptions
DescriptionsRegistry().set_descriptions({**syms_dict})
DescriptionsRegistry().reset_registry()
SymbolsDescription.set_default_header('where:')

display(ReportText('Symbols description is as follows: '))


display(SymbolsDescription(expr=E_K*F))


display(ObjectCode("""
# Guide
from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
BasicsOfODESystemGuide();

"""))


# Basics of DynSys usage based simulations are covered in guide to DynSys usage, use the followin call
# from dynpy.utilities.documents.guides import UsageOfDynamicSystemsGuide
# UsageOfDynamicSystemsGuide()
'''

        picture_str = """sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))"""

        simulaion_str = """from dynpy.utilities.report import *
from sympy import *


sec_simulation = Section('Section that contains simulation')
CurrentContainer(sec_simulation)

# Basics of ODESystem based simulations are covered in guide to ODESystem, use the followin call
# from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
# BasicsOfODESystemGuide()

sec_ODESystem = Subsection('ODESystem simulation')
CurrentContainer(sec_ODESystem)

display(ReportText('Firstly create an ODESystem or import it from dynpy.modes.odes.linear.py'))
display(ObjectCode('spring = SpringMassEps.from_reference_data()'))

display(ReportText('Secondly solve the equation using solution method:'))
display(ObjectCode('spring_sol = spring.solution'))

display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
data_spring = dyn_sys.get_numerical_parameters()
#ds1.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()

# MSM_sim = Subsection('MSM simulation')
# CurrentContainer(MSM_sim)

# display(ReportText('Firstly create an MultiTimeScaleSolution system or import it from dynpy.modes.odes.linear.py'))
# display(ObjectCode('spring = SpringMassMSM.from_reference_data()'))

# display(ReportText('Secondly solve the equation using solution method:'))
# display(ObjectCode('spring_sol = spring.solution'))

# display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
# spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()"""

        veryfication_str = """from sympy import *

from sympy import *
from sympy.physics import units




from dynpy.utilities.adaptable import *
from dynpy.utilities.report import CurrentContainer, ReportText
import pandas as pd

t = Symbol('t')
power = Symbol('P')
work = Symbol('W')

unit_dict = {
t: units.second,
power: units.watt,
work: units.joule,
}

LatexDataFrame.set_default_units(unit_dict)



sec_verification = Section('Section that verificates research results')
CurrentContainer(sec_verification)

display(ReportText('Description of verifications concept. '*200))



### SELECT ONE CASE
## 1st CASE
## import of external data if path "./data_from_ml/test_1.csv" exists

#test_data = './data_from_ml/test_1.csv'
#df = pd.read_csv(test_data, header = None)

#data_tab = TimeDataFrame(df).to_latex_dataframe().reported(caption='Caption')
#graph=data_tab.iloc[:,[1]]

#display(data_tab)
display(ReportText('Obtained data analysis. '*200))
#display(graph)

## 2st CASE
## Creating graph from a
df = pd.DataFrame({t:[0,1,2],'a':[2,3,4],'b':[8,7,6]}).set_index(t)

df_cols = df.set_axis([power,work],axis=1)

data_tab = TimeDataFrame(df_cols).to_latex_dataframe().reported(caption='Caption')
graph = TimeDataFrame(df_cols).to_pylatex_tikz().in_figure(caption='Caption')

display(data_tab)
display(ReportText('Obtained data analysis. '*200))

display(graph)

## 3rd CASE
## dictionary data presentation in DataFrame
#dictionary ={
#     power:300000,
#     work:5000
# }
#display(
#     LatexDataFrame.formatted(
#         data=dictionary, #import data from a dictionary
#         index=['Value'] #set rows title
#     ).map(lambda x: f'${latex(x)}$')
#     .rename_axis('Parameter', axis=1) #set name for a column
#     .transpose()  # This swaps rows and columns, making the DataFrame vertical
#     .reported(caption='A caption for your data frame')
# )
#display(ReportText('Obtained data analysis. '*200))


## 4th CASE
## creating a graph based on simulation
#from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
#F=DynamicSys.F
#g=DynamicSys.g
#m=DynamicSys.m
#k=DynamicSys.k
#slownik ={
#F:2,
#g:10,
#m:5,
#k:100,
#}
#t_span = np.linspace(0.0,1,30)
#sym=DynamicSys().eoms.subs(slownik)
#wynik= sym.numerized(backend='numpy').compute_solution(t_span,[0.0,0.0])#.plot()
#wynik_tab = TimeDataFrame(wynik).to_latex_dataframe().reported(caption='Caption')
#display(wynik_tab)
#display(ReportText('Obtained data analysis. '*200))
#wynik.plot()


display(ReportText('Description of verifications outcomes. '*200))"""

        conclusion_str = """sec_conclusion = Section('Section that contains final conclusions')
CurrentContainer(sec_conclusion)

display(ReportText('Conclusions '*200))"""

        symbols_description_str = """sec_symbols = Section('Section that contains all symbols descriptions')
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

display(SymbolsDescription({**syms_dict}))"""

        document_str = """# Creating file
# Be sure *output* folder is in the current directory

paper_name = './output/report_name' #path for report file
global_width=('8cm')#setting picture default size
global_height=('6cm')
doc = BPASTSPaper(paper_name)
# Bibliography of quotations
### Select one bibligraphy managment system
####### BibLatex

doc.preamble.append(Package('biblatex',["backend=biber","sorting=none"]))
doc.preamble.append(Command('addbibresource','elementy_bibliagrafia.bib'))
####### Natbib
#doc.preamble.append(Package('natbib')
#doc.preamble.append(Command('bibliographystyle','unsrt'))

## Units
doc.preamble.append(Package('siunitx'))




doc.append(sec_intro) # adding certain sections
doc.append(sub_problem_outline)
doc.append(sub_obj_assum)
doc.append(sub_SOT)
doc.append(sub_methodology)
doc.append(sec_formula)

doc.append(sub_model_code)
doc.append(sec_math_desc)
doc.append(sec_simulation)

doc.append(sec_picture)
doc.append(sec_verification)
#doc.append(sec_tables)
doc.append(sec_symbols)
doc.append(sec_conclusion)




### Natbib
doc.append(Command('bibliography',arguments=["references"])) # .bib file as "references"


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
doc.generate_pdf(clean_tex=False)"""

        #         preliminary_str=(f"""
        # =======
        preliminary_str = f"""

Examplary setup is as follows:

#References to guide imports needed to include - to see the list of available guides insert and run the code below:
```python
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```

## CELL 1
## Imports

```python
{imports_str}


```



## CELL 2
## Thesis introduction

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{thesis_introduction_str}
```



## CELL 3
## Math

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{math_str}
```


## CELL 4
## Picture


```python

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{picture_str}

```


## CELL 5
## Simulation

Basics of ODESystem based simulations are covered in guide to ODESystem, use the following call:

```python
from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
BasicsOfODESystemGuide();
```


```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#


{simulaion_str}

```


## CELL 6
## Veryfication

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{veryfication_str}

```

## CELL 7
## Conclusion

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{conclusion_str}

```

## CELL 9
## Symbols description

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{symbols_description_str}

```

## CELL 10
## Document

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{document_str}
```

"""

        display(IPMarkdown(preliminary_str))

        preliminary_str = f"""
#Example:

#To prepare a simple document with text and images:

#Good practice here is to allocate 1 section per 1 cell

## CELL 1
## Imports

{imports_str}

## CELL 2
## Thesis introduction

{thesis_introduction_str}

## CELL 3
## Math

{math_str}

## CELL 4
## Picture

{picture_str}

## CELL 5
## Simulation

{simulaion_str}

## CELL 6
## Veryfication

{veryfication_str}

## CELL 7
## Conclusion

{conclusion_str}

## CELL 9
## Symbols description

{symbols_description_str}

## CELL 10
## Document

{document_str}

"""

        return ObjectCode(preliminary_str)


class APPAPaper(Document):

    latex_name = "document"
    packages = [
        #                     Package('natbib', options=['numbers']),
        Package("booktabs"),
        Package("float"),
        Package("standalone"),
        Package("siunitx"),
        #                     Package('bpasts', options=['accepted']),
        #Package("bpasts"),
        Package("t1enc"),
        Package("amsmath"),
        Package("amssymb"),
        Package("amsfonts"),
        Package("graphicx"),
        Package("flushend"),
        Package("textcomp"),
        Package("xcolor"),
#         Package(
#             "hyperref",
#             ["colorlinks=true", "allcolors=bpastsblue", NoEscape("pdfborder={0 0 0}")],
#         ),
    ]

    abtitle = "Paper for BPASTS"
    abauthor = "Authors"
    title = "Basic title"
    author = r"""
Author's Name Surname  - Co-author's Name Surname
"""
    abstract=r"""
The Abstract should not exceed 250 words. The Abstract should state the principal objectives and the scope of the investigation, as well as the methodology employed. It should summarize the results and state the principal conclusions. An effective abstract stands on its own — it can be understood fully even when made available without the full paper. To this end, avoid referring to figures or the bibliography in the abstract. Please introduce any acronyms the first time you use them in the abstract (if needed), and do so again in the full paper. About 4 to 6 significant key words should follow the abstract to aid indexing.
"""

    keywords = "key word, key word, key word, key word, key word"

    highlights=r"""
\item Highlights (4 to 6) are a short collection of bullet points that convey the core findings and provide readers with a quick textual overview of the article.
\item These four to six bullet points should describe the essence of the research (e.g. results or conclusions) and highlight what is distinctive about it.
\item See latest SV-JME papers for examples.
"""


    def __init__(
        self,
        default_filepath="default_filepath",
        title=None,
        *,
        documentclass="article",
        document_options=["10pt",
                          #"twoside", 
                        "twocolumn", "a4paper",
                         ],  # for submission
        fontenc=None,
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=False,
        microtype=False,
        page_numbers=True,
        indent=None,
        #geometry_options=[
        #    "inner=30mm",
        #    "outer=20mm",
        #    "bindingoffset=10mm",
        #    "top=25mm",
        #    "bottom=25mm",
        #],  # ,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
        geometry_options=None,
        data=None,
    ):

        if title is not None:
            self.title = title

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



#         self.preamble.append(NoEscape(r"""
# % to be set by Editor:
# % ------------------------------------------------------------------------------------------------------------------
# %\renewcommand{\svCorr}{*Corr. Author's Address: Name of institution, Address, City, Country, \email{xxx.yyy@wwww.zzz}}
# \renewcommand{\svSV}{Strojniški vestnik - Journal of Mechanical Engineering 63(2017)3,  XXX-4}
# \renewcommand{\svAuthors}{Author-A, Author-B}
# \renewcommand{\svCyear}{2017}
# \renewcommand{\svDOI}{DOI: 10.5545/sv-jme.2017.4027}
# \svDates{2016-11-04}{2017-01-14}{2017-02-12}
# \renewcommand{\svType}{Review Paper}
# \setcounter{page}{1} % first page of the paper
# % ------------------------------------------------------------------------------------------------------------------
# """))


        self.preamble.append(Command("title", arguments=[self.title]))
        #self.preamble.append(NoEscape("\\title{" + f"{self.title}" + "}"))
        #self.preamble.append(NoEscape("\\abauthor{" + f"{self.abauthor}" + "}"))


        self.preamble.append(Command("author", arguments=[NoEscape(self.author)]))
        #self.preamble.append(Command('affil', 'Politechnika Warszawska'))

#         self.preamble.append(NoEscape(r"""

# %\date{\today\  -- \clock}
# \date{}

# \crop[cross,axes]
# %\crop[frame,axes]
# %\crop[off]"""))

        self.append(NoEscape(r"""
\twocolumn[
\maketitle
\begin{abstract}
"""))

        self.append(self.abstract)

        self.append(NoEscape(r"""
\end{abstract}]
"""))
#         self.append(Command("svKeywords", arguments=[NoEscape(self.keywords)]))

#         self.append(NoEscape(r"""
# \begin{svHigh}
# """))
#         self.append(NoEscape(self.highlights))

#         self.append(NoEscape(r"""

# \end{svHigh}

# \end{svHead}]
#"""))


        # self.append(NewPage())
        # tu implementować co tam potrzeba

    # def cwd_setup(self):
    #     cwd = os.getcwd()

    #     source_path = (
    #         f"/home/user/Shared files/modules/dynpy/utilities/documents/Definitions"
    #     )
    #     path1 = f"{cwd}/bpast.sty"
    #     path2 = f"{cwd}/output/bpast.sty"

    #     if os.path.exists(path1) == False:
    #         shutil.copytree(source_path, path1)
    #     if os.path.exists(path2) == False:
    #         shutil.copytree(source_path, path2)

class JoMEPaper(Document):

    latex_name = "document"
    packages = [
        #                     Package('natbib', options=['numbers']),
        Package("booktabs"),
        Package("float"),
        Package("standalone"),
        Package("siunitx"),
        #                     Package('bpasts', options=['accepted']),
        #Package("bpasts"),
        Package("t1enc"),
        Package("amsmath"),
        Package("amssymb"),
        Package("amsfonts"),
        Package("graphicx"),
        Package("flushend"),
        Package("textcomp"),
        Package("xcolor"),
#         Package(
#             "hyperref",
#             ["colorlinks=true", "allcolors=bpastsblue", NoEscape("pdfborder={0 0 0}")],
#         ),
    ]

    abtitle = "Paper for BPASTS"
    abauthor = "Authors"
    title = "Basic title"
    author = r"""
Author's Name Surname\svMark{1,*}  - Co-author's Name Surname\svMark{2} \\[1mm]
\svAffil{1}{Author's institution} \\
\svAffil{2}{Co-Author's institution}
"""
    abstract=r"""
The Abstract should not exceed 250 words. The Abstract should state the principal objectives and the scope of the investigation, as well as the methodology employed. It should summarize the results and state the principal conclusions. An effective abstract stands on its own — it can be understood fully even when made available without the full paper. To this end, avoid referring to figures or the bibliography in the abstract. Please introduce any acronyms the first time you use them in the abstract (if needed), and do so again in the full paper. About 4 to 6 significant key words should follow the abstract to aid indexing.
"""

    keywords = "key word, key word, key word, key word, key word"

    highlights=r"""
\item Highlights (4 to 6) are a short collection of bullet points that convey the core findings and provide readers with a quick textual overview of the article.
\item These four to six bullet points should describe the essence of the research (e.g. results or conclusions) and highlight what is distinctive about it.
\item See latest SV-JME papers for examples.
"""


    def __init__(
        self,
        default_filepath="default_filepath",
        title=None,
        *,
        documentclass="JoME",
        document_options=["10pt",
                          #"twoside", "twocolumn", "a4paper"
                         ],  # for submission
        fontenc=None,
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=False,
        microtype=False,
        page_numbers=True,
        indent=None,
        #geometry_options=[
        #    "inner=30mm",
        #    "outer=20mm",
        #    "bindingoffset=10mm",
        #    "top=25mm",
        #    "bottom=25mm",
        #],  # ,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
        geometry_options=None,
        data=None,
    ):

        if title is not None:
            self.title = title

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



        self.preamble.append(NoEscape(r"""
% to be set by Editor:
% ------------------------------------------------------------------------------------------------------------------
%\renewcommand{\svCorr}{*Corr. Author's Address: Name of institution, Address, City, Country, \email{xxx.yyy@wwww.zzz}}
\renewcommand{\svSV}{Strojniški vestnik - Journal of Mechanical Engineering 63(2017)3,  XXX-4}
\renewcommand{\svAuthors}{Author-A, Author-B}
\renewcommand{\svCyear}{2017}
\renewcommand{\svDOI}{DOI: 10.5545/sv-jme.2017.4027}
\svDates{2016-11-04}{2017-01-14}{2017-02-12}
\renewcommand{\svType}{Review Paper}
\setcounter{page}{1} % first page of the paper
% ------------------------------------------------------------------------------------------------------------------
"""))


        self.preamble.append(Command("svTitle", arguments=[self.title]))
        self.preamble.append(NoEscape("\\title{" + f"{self.title}" + "}"))
        #self.preamble.append(NoEscape("\\abauthor{" + f"{self.abauthor}" + "}"))


        self.preamble.append(Command("author", arguments=[NoEscape(self.author)]))

        self.preamble.append(NoEscape(r"""

%\date{\today\  -- \clock}
\date{}

\crop[cross,axes]
%\crop[frame,axes]
%\crop[off]"""))

        self.append(NoEscape(r"""
\twocolumn[\begin{svHead}

\begin{svAbstract}
"""))

        self.append(self.abstract)

        self.append(NoEscape(r"""
\end{svAbstract}
"""))
        self.append(Command("svKeywords", arguments=[NoEscape(self.keywords)]))

        self.append(NoEscape(r"""
\begin{svHigh}
"""))
        self.append(NoEscape(self.highlights))

        self.append(NoEscape(r"""

\end{svHigh}

\end{svHead}]
"""))


        # self.append(NewPage())
        # tu implementować co tam potrzeba

    def cwd_setup(self):
        cwd = os.getcwd()

        source_path = (
            f"/home/user/Shared files/modules/dynpy/utilities/documents/Definitions"
        )
        path1 = f"{cwd}/bpast.sty"
        path2 = f"{cwd}/output/bpast.sty"

        if os.path.exists(path1) == False:
            shutil.copytree(source_path, path1)
        if os.path.exists(path2) == False:
            shutil.copytree(source_path, path2)



class Keywords(Environment):
    latex_name = "keyword"


class Abstract(Environment):
    latex_name = "abstract"


class ElsevierPaper(Document):
    """
    Doc = ElsevierTemplate(default_filepath='./output/paper')
    Doc.journal('Advances in Engineering Software')

    Doc.begin_frontmatter

    Doc.corresponding_author('Damian Sierociński','damian.sierocinski@pw.edu.pl')
    Doc.authors(['Bogumił Chiliński','Franciszek Gawiński','Piotr Przybyłowicz','Amadeusz Radomski'])
    Doc.abstract([NoEscape('\\textit{DynPy} is an open-source library implemented in \\textit{Python} programming language which aims to provide a versatile set of functionalities, that enable the user to model, solve, simulate, and report an analysis of a dynamic system in a single environment. In the paper examples for obtaining analytical and numerical solutions of the systems described with ordinary differential equations were presented. The assessment of solver accuracy was conducted utilising a model of a direct current motor and \\textit{MATLAB/Simulink} was used as a reference tool. The model was solved in \\textit{DynPy} with the hybrid analytical-numerical method and fully analytically, while in \\textit{MATLAB/Simulink} strictly numerical simulations were run. The comparison of the results obtained from both tools not only proved the credibility of the developed library but also showed its superiority in specific conditions. Moreover, the versatility of the tool was confirmed, as the whole process – from the definition of the model to the output of the formatted paper – was handled in the scope of a single \\textit{Jupyter} notebook.')])
    Doc.keywords(['engineering software','numerical simulations','analytical solution','electrical circuit',NoEscape('\\textit{Python} programming language')])
    Doc.title(NoEscape("An assessment of \\textit{DynPy} solver accuracy"))
    Doc.end_frontmatter
    """

    latex_name = "document"
    packages = [
        #                   Package('geometry',options=['lmargin=30mm', 'rmargin=30mm',  'top=25mm', 'bmargin=25mm', 'headheight=50mm']),
        #                   Package('microtype'),
        #                   Package('authoraftertitle'),
        #                   Package('listings'),
        #                   Package('titlesec'),
        #                   Package('fancyhdr'),
        #                   Package('graphicx'),
        #                   Package('indentfirst'),
        #                   Package('pdfpages'),
        Package("fontspec"),
        Package("amsmath"),
        Package("amssymb"),
        Package("epsfig"),
        #                   Command('graphicspath{{../}}'),
        #                   Command('frenchspacing'),
        #                   Command('counterwithin{figure}{section}'),
        #                   Command('counterwithin{table}{section}'),
        #                   Command('fancypagestyle{headings}{\\fancyhead{} \\renewcommand{\headrulewidth}{1pt} \\fancyheadoffset{0cm} \\fancyhead[RO]{\\nouppercase{\\leftmark}} \\fancyhead[LE]{\\nouppercase{\\leftmark}} \\fancyfoot{} \\fancyfoot[LE,RO]{\\thepage}}'),
        #                   Command('fancypagestyle{plain}{\\fancyhf{} \\renewcommand{\\headrulewidth}{0pt} \\fancyfoot[LE,RO]{\\thepage}}'),
        #                   Command('numberwithin{equation}{section}'),
        #                   Command('renewcommand{\\familydefault}{\\sfdefault}'),
        # \renewcommand{\familydefault}{\sfdefault}
    ]

    def __init__(
        self,
        default_filepath="default_filepath",
        title="Basic title",
        *,
        documentclass="elsarticle",
        document_options=["final", "times"],  # for submission
        #                  document_options=['final', '5p', 'times' ,'twocolumn'], for preview
        fontenc=None,
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=True,
        microtype=True,
        page_numbers=True,
        indent=None,
        geometry_options=[
            "inner=30mm",
            "outer=20mm",
            "bindingoffset=10mm",
            "top=25mm",
            "bottom=25mm",
        ],  # ,inner=20mm, outer=20mm, bindingoffset=10mm, top=25mm, bottom=25mm
        data=None,
    ):

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

    @property
    def begin_frontmatter(self):
        self.append(ReportText(NoEscape("\\begin{frontmatter}")))

    @property
    def end_frontmatter(self):
        self.append(ReportText(NoEscape("\\end{frontmatter}")))

    def corresponding_author(self, author=None, email=None):

        self.append(
            Command("author[label1]", arguments=[NoEscape(f"{author}\\corref{{cor1}}")])
        )
        self.append(Command("ead", arguments=[f"{email}"]))
        #         self.append(Command('fntext[label2]'))
        self.append(Command("cortext[cor1]{Corresponding author.}"))
        self.append(
            Command(
                "affiliation[label1]",
                arguments=[
                    NoEscape(
                        "organization={Warsaw Univeristy of Technology, Faculty of Automotive and Constrction Machinery Engineering},addressline={Ludwika Narbutta 84},city={Warsaw},postcode={02-524},state={Masovian Voivodeship},country={Poland}"
                    )
                ],
            )
        )

    def keywords(self, keywords_list=None):
        kwl = len(keywords_list)
        with self.create(Keywords()):
            for num, val in enumerate(keywords_list):
                self.append(val)
                if num < kwl - 1:
                    self.append(NoEscape("\sep"))

    def authors(self, authors_list=None):
        self.authors_list = authors_list

        for num, val in enumerate(self.authors_list):
            self.append(Command("author[label1]", arguments=[val]))

    def journal(self, journal_name=None):
        self.preamble.append(Command("journal", arguments=[journal_name]))

    def title(self, title=None):
        self.preamble.append(Command("title", arguments=[title]))

    def abstract(self, container=None):
        with self.create(Abstract()):
            for num, val in enumerate(container):
                self.append(ReportText(val))


class WutThesis(Document):
    """
    A class for creating a Warsaw University of Technology (WUT) thesis document.

    This class extends the `Document` class to include all necessary LaTeX packages, commands, and configurations required for WUT thesis formatting. It provides an easy way to create and manage thesis documents with pre-configured settings.

    Attributes:
        latex_name (str): Name of the document type ('document').
        packages (list): List of LaTeX packages and commands preloaded for thesis formatting.

    Example Usage:
        >>> from dynpy.utilities.documents.document import WutThesis
        >>> from pylatex import Section
        >>> doc = WutThesis('./output/thesis_name')
        >>> section = Section('Introduction')
        >>> doc.append(section)
        >>> doc.generate_pdf(clean_tex=True)

    Example of Full Thesis Workflow:
        >>> from dynpy.utilities.documents.document import WutThesis
        >>> from pylatex import Section, Command
        >>> doc = WutThesis('./output/thesis_name')
        >>> doc.append(Command('title', arguments=['Thesis Title']))
        >>> doc.append(Command('author', arguments=['Student Name']))
        >>> doc.append(Command('date', arguments=['\today']))
        >>> doc.append(Command('maketitle'))
        >>> section = Section('Introduction')
        >>> doc.append(section)
        >>> doc.generate_pdf(clean_tex=True)
    """

    latex_name = "document"
    packages = [
        Package(
            "geometry",
            options=[
                "lmargin=30mm",
                "rmargin=30mm",
                "top=25mm",
                "bmargin=25mm",
                "headheight=50mm",
            ],
        ),
        Package("microtype"),
        Package("authoraftertitle"),
        # Package('polski',options=['MeX']),
        # Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
        Package("listings"),
        Package("titlesec"),
        Package("fancyhdr"),
        Package("graphicx"),
        Package("indentfirst"),
        Package("pdfpages"),
        Package("amsmath"),
        Package("markdown"),
        Package("hyperref"),
        Command("newcommand{\praca}", arguments=["Praca dyplomowa"]),
        Command("newcommand{\dyplom}", arguments=["Magisterska"]),
        Command("newcommand{\kierunek}", arguments=["Wpisać kierunek"]),
        Command("newcommand{\specjalnosc}", arguments=["Wpisać specjalność"]),
        Command("newcommand{\\autor}", arguments=["Imię i nazwisko autora"]),
        Command("newcommand{\opiekun}", arguments=["Wpisać opiekuna"]),
        Command("newcommand{\promotor}", arguments=["Wpisać promotora"]),
        Command("newcommand{\konsultant}", arguments=["Wpisać konsultanta"]),
        Command(
            "newcommand{\\tytul}", arguments=["Wpisać tytuł pracy dyplomowej po polsku"]
        ),
        Command("newcommand{\\album}", arguments=["Wpisać numer albumu"]),
        Command("newcommand{\supervisor}", arguments=["dr inż. Bogumił Chiliński"]),
        Command("newcommand{\\rok}", arguments=["Rok składania pracy"]),
        Command(
            "newcommand{\kluczowe}",
            arguments=["Słowa kluczowe: Wpisać słowa kluczowe po polsku"],
        ),
        # Command('renewcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
        Command("graphicspath{{../}}"),
        Command("frenchspacing"),
        Command("counterwithin{figure}{section}"),
        Command("counterwithin{table}{section}"),
        Command(
            "fancypagestyle{headings}{\\fancyhead{} \\renewcommand{\headrulewidth}{1pt} \\fancyheadoffset{0cm} \\fancyhead[RO]{\\nouppercase{\\leftmark}} \\fancyhead[LE]{\\nouppercase{\\leftmark}} \\fancyfoot{} \\fancyfoot[LE,RO]{\\thepage}}"
        ),
        Command(
            "fancypagestyle{plain}{\\fancyhf{} \\renewcommand{\\headrulewidth}{0pt} \\fancyfoot[LE,RO]{\\thepage}}"
        ),
        Command("numberwithin{equation}{section}"),
        Command("renewcommand{\\familydefault}{\\sfdefault}"),
        # \renewcommand{\familydefault}{\sfdefault}
    ]

    def __init__(
        self,
        default_filepath: str = "default_filepath",
        title: str = "Basic title",
        *,
        documentclass: str = "article",
        document_options: list = ["a4paper", "11pt", "twoside"],
        fontenc: str = "T1",
        inputenc: str = "utf8",
        font_size: str = "normalsize",
        lmodern: bool = False,
        textcomp: bool = True,
        microtype: bool = True,
        page_numbers: bool = True,
        indent: Optional[Union[str, int]] = None,
        geometry_options: Optional[list] = [
            "inner=30mm",
            "outer=20mm",
            "bindingoffset=10mm",
            "top=25mm",
            "bottom=25mm",
        ],
        data: Optional[dict] = None,
    ):
        """
        Initialize the WutThesis class with optional customization.

        Args:
            default_filepath (str): Path for the generated document.
            title (str): Title of the document.
            documentclass (str): LaTeX document class (e.g., 'article').
            document_options (list): Options for the document class.
            fontenc (str): Font encoding (default 'T1').
            inputenc (str): Input encoding (default 'utf8').
            font_size (str): Font size (default 'normalsize').
            lmodern (bool): Whether to use Latin Modern fonts.
            textcomp (bool): Whether to use the `textcomp` package.
            microtype (bool): Whether to use the `microtype` package.
            page_numbers (bool): Whether to include page numbers.
            indent (Optional[Union[str, int]]): Indentation settings.
            geometry_options (Optional[list]): Geometry options for page layout.
            data (Optional[dict]): Additional data for customization.
        """

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
        self.title = title
        # self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
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
        self.packages.append(Command("graphicspath{{../}}"))
        #         self.append(Command('maketitle'))
        self.append(NoEscape("%%% New doc"))
        # tu implementować co tam potrzeba

    @classmethod
    def base_setup(cls, create_file: bool = False):
        """
        Sets up the base template for the WUT thesis document.

        This method provides an easy way to initialize a thesis document with the necessary sections, structure, and configuration.

        Args:
            create_file (bool): If True, creates a base setup environment including necessary directories and files.

        Returns:
            Optional[str]: Returns setup content if `create_file` is False, otherwise creates the environment.

        Usage Example:
            >>> from dynpy.utilities.documents.document import WutThesis
            >>> WutThesis.base_setup()

        When to use:
            - Use this method when starting a new WUT thesis document.
            - If `create_file=True`, it will generate all required directories and files for a structured thesis.
            - If `create_file=False`, it returns a string with setup instructions that can be used manually.
        """

        if create_file is True:
            return cls._create_base_setup_env()

        # doc = WutThesis('./output/thesis_name')
        # doc.preamble.append(NoEscape(r'\\usepackage[MeX]{polski}')) #to set polish as main language"""

        thesis_introduction_str = """#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
from dynpy.utilities.report import *

sec_intro = Section('Section that presents text reporting')
CurrentContainer(sec_intro)

sub_problem_outline = Subsection('Outline of the problem')
CurrentContainer(sub_problem_outline)
display(ReportText('This subsection provides information about investigated problem. '*1))

sub_obj_assum = Subsection('Objectives and assumptions')
CurrentContainer(sub_obj_assum)
display(ReportText('This subsection provides objectives and assumptions. '*1))

sub_SOT = Subsection('State of the art')
CurrentContainer(sub_SOT)
display(ReportText('This subsection provides state of the art. '*1))

sub_methodology = Subsection('Methodology')
CurrentContainer(sub_methodology)
display(ReportText('This subsection provides methodology. '*1))"""

        math_str = '''#!!! NO NEEDED RUN PREVIOUS CELLS !!!#
from sympy import Eq, Symbol, symbols
from dynpy.utilities.report import *


sec_formula = Section('Section that presents formulas reporting')
CurrentContainer(sec_formula)

sec_description = Subsection('Description of dynamic model')
CurrentContainer(sec_description)
display(ReportText('This subsection provides description of the model. '*1))

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys

dyn_sys = DynamicSys()
display(dyn_sys.as_picture())

display(ReportText('Summary of the subsection highlighting its major achievements. '*1))

sub_model_code = Subsection('Virtual model of the object')
CurrentContainer(sub_model_code)

from dynpy.utilities.report import ObjectCode

display(ReportText('This subsection provides code of the model. '*1))
display(ObjectCode(DynamicSys))
display(ReportText('Summary of the subsection highlighting its major achievements. '*1))


sec_lagrangian = Subsection('Lagrangian and derivatives')
CurrentContainer(sec_lagrangian)
display(ReportText('This subsection provides calculation of lagrangian and derivatives. '*1))

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

display(ReportText('Outcomes of governing equations analysis. '*10))


sec_equations = Subsection('Equations of motion')
CurrentContainer(sec_equations)


display(ReportText('This subsection provides calculation of equations of motion and solution of the system. '*1))

ds1=dyn_sys.eoms
ds1_eqns=ds1.as_eq_list()

for eq1 in ds1.as_eq_list():
    display(SympyFormula(eq1.simplify()))

if ds1.is_solvable():

    ds2=dyn_sys.linearized()._ode_system.general_solution

    for eq2 in ds2.as_eq_list():
        display(SympyFormula(eq2.simplify()))

    ds3=dyn_sys.linearized()._ode_system.steady_solution

    for eq3 in ds3.as_eq_list():
        display(SympyFormula(eq3.simplify()))

    ds4=dyn_sys.linearized().eoms.solution
    ds4_eqns = ds4.as_eq_list()

    for eq4 in ds4_eqns:
        display(SympyFormula(eq4.simplify()))


display(ReportText(f'Outcomes of governing equations {AutoMarker(ds1_eqns[0])} {AutoMarker(ds1_eqns[-1])}) analysis.'))

sec_math_desc = Subsection('Section that contains all symbols descriptions')
CurrentContainer(sec_math_desc)



E_K,F = symbols('E_K,F')
descriptions = {
E_K:r"Kinetic energy",
F:r"Force",
}
syms_dict = descriptions
DescriptionsRegistry().set_descriptions({**syms_dict})
DescriptionsRegistry().reset_registry()
SymbolsDescription.set_default_header('where:')

display(ReportText('Symbols description is as follows: '))


display(SymbolsDescription(expr=E_K*F))


display(ObjectCode("""
# Guide
from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
BasicsOfODESystemGuide();

"""))


# Basics of DynSys usage based simulations are covered in guide to DynSys usage, use the followin call
# from dynpy.utilities.documents.guides import UsageOfDynamicSystemsGuide
# UsageOfDynamicSystemsGuide()
'''

        picture_str = """
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
from dynpy.utilities.report import *

sec_picture = Section('Section that presents pictures reporting')
CurrentContainer(sec_picture)

display(Picture('./dynpy/models/images/taipei101.png',caption = 'Caption of picture'))"""

        simulaion_str = """#!!! NO NEEDED RUNNING PREVIOUS CELLS !!!#


from dynpy.utilities.report import *
from sympy import *
import numpy as np

sec_simulation = Section('Section that contains simulation')
CurrentContainer(sec_simulation)

# Basics of ODESystem based simulations are covered in guide to ODESystem, use the followin call
# from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
# BasicsOfODESystemGuide()

sec_ODESystem = Subsection('ODESystem simulation')
CurrentContainer(sec_ODESystem)

display(ReportText('Firstly create an ODESystem or import it from dynpy.modes.odes.linear.py'))
display(ObjectCode('spring = SpringMassEps.from_reference_data()'))

display(ReportText('Secondly solve the equation using solution method:'))
display(ObjectCode('spring_sol = spring.solution'))

display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
dyn_sys = DynamicSys()

data_spring = dyn_sys.get_numerical_parameters()
#ds1.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()

from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
F=DynamicSys.F
g=DynamicSys.g
m=DynamicSys.m
k=DynamicSys.k

slownik ={
F:2,
g:10,
m:5,
k:100,
}

t_span = np.linspace(0.0,1,30)
sym=DynamicSys().eoms.subs(slownik)
wynik= sym.numerized(backend='numpy').compute_solution(t_span,[0.0,0.0])#.plot()
wynik_tab = TimeDataFrame(wynik).to_latex_dataframe().reported(caption='Caption')
display(wynik_tab)
wynik.plot()


# MSM_sim = Subsection('MSM simulation')
# CurrentContainer(MSM_sim)

# display(ReportText('Firstly create an MultiTimeScaleSolution system or import it from dynpy.modes.odes.linear.py'))
# display(ObjectCode('spring = SpringMassMSM.from_reference_data()'))

# display(ReportText('Secondly solve the equation using solution method:'))
# display(ObjectCode('spring_sol = spring.solution'))

# display(ReportText('Last step is data substitution and using numerized method. After that use $compute_solution$ to finalize the simulation:'))
# spring_sol.subs(data_spring).with_ics([10,0]).numerized().compute_solution(t_span).plot()"""

        veryfication_str = """#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

from sympy import *

from sympy.physics import units



from dynpy.utilities.adaptable import *
from dynpy.utilities.report import CurrentContainer, ReportText
import pandas as pd

t = Symbol('t')
power = Symbol('P')
work = Symbol('W')

unit_dict = {
t: units.second,
power: units.watt,
work: units.joule,
}



LatexDataFrame.set_default_units(unit_dict)



sec_verification = Section('Section that verificates research results')
CurrentContainer(sec_verification)

display(ReportText('Description of verifications concept. '*200))



### SELECT ONE CASE
## 1st CASE
##import of external data if path "./data_from_ml/test_1.csv" exists

# test_data = './data_from_ml/test_1.csv'
# df = pd.read_csv(test_data, header = None)

# data_tab = TimeDataFrame(df).to_latex_dataframe().reported(caption='Caption')
# graph=data_tab.iloc[:,[1]]

## 2st CASE
## Creating graph from a
# df = pd.DataFrame({t:[0,1,2],'a':[2,3,4],'b':[8,7,6]}).set_index(t)

# df_cols = df.set_axis([power,work],axis=1)

# data_tab = TimeDataFrame(df_cols).to_latex_dataframe().reported(caption='Caption')
# graph = TimeDataFrame(df_cols).to_pylatex_tikz().in_figure(caption='Caption')

# display(data_tab)

# display(ReportText('Obtained data analysis. '*200))

# display(graph)

# display(ReportText('Description of verifications outcomes. '*200))

## 3rd CASE
## dictionary data presentation in DataFrame
# dictionary ={
#     power:300000,
#     work:5000
# }
# display(
#     LatexDataFrame.formatted(
#         data=dictionary, #import data from a dictionary
#         index=['Value'] #set rows title
#     ).map(lambda x: f'${latex(x)}$')
#     .rename_axis('Parameter', axis=1) #set name for a column
#     .transpose()  # This swaps rows and columns, making the DataFrame vertical
#     .reported(caption='A caption for your data frame')
# )

## 4th CASE
## creating a graph based on simulation
# from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
# F=DynamicSys.F
# g=DynamicSys.g
# m=DynamicSys.m
# k=DynamicSys.k
# slownik ={
# F:2,
# g:10,
# m:5,
# k:100,
# }
# t_span = np.linspace(0.0,1,30)
# sym=DynamicSys().eoms.subs(slownik)
# wynik= sym.numerized(backend='numpy').compute_solution(t_span,[0.0,0.0])#.plot()
# wynik_tab = TimeDataFrame(wynik).to_latex_dataframe().reported(caption='Caption')
# display(wynik_tab)
# wynik.plot()

"""

        conclusion_str = """### NO NEEDED RUN OF ALL PREVIOUS CELLS ###

from dynpy.utilities.report import *

sec_conclusion = Section('Section that contains final conclusions')
CurrentContainer(sec_conclusion)

display(ReportText('Conclusions '*200))"""

        symbols_description_str = """#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

from dynpy.utilities.report import *

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

display(SymbolsDescription({**syms_dict}))"""

        document_str = """#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
from dynpy.utilities.documents.document import WutThesis
from dynpy.utilities.report import *

# Creating file
# Be sure *output* folder is in the current directory

thesis_name = './output/report_name' #path for report file

doc = WutThesis(thesis_name)
# Bibliography of quotations
### Select one bibligraphy managment system
####### BibLatex

doc.preamble.append(Package('biblatex',["backend=biber","sorting=none"]))
doc.preamble.append(Command('addbibresource','elementy_bibliagrafia.bib'))
####### Natbib
#doc.preamble.append(Package('natbib')
#doc.preamble.append(Command('bibliographystyle','unsrt'))

## Units
doc.preamble.append(Package('siunitx'))

# TOC

###############################################################################
######### UNCOMMENT IT IF YOU NEED ADD WUT THESIS TITLE PAGES

# doc.append(Command('includepdf{./Images/Front_Page.pdf}')) #includes front page
# doc.append(Command('pagestyle{plain}'))
# doc.append(Command('cleardoublepage'))
# doc.append(Command('includepdf{./Images/oswiadczenie_autora_pracy.pdf}'))
# doc.append(Command('pagestyle{plain}'))
# doc.append(Command('cleardoublepage'))
# doc.append(Command('includepdf{./output/oswiadczenie_biblioteczne.pdf}'))
# doc.append(Command('pagestyle{plain}'))
# doc.append(Command('cleardoublepage'))

###############################################################################
###############################################################################


doc.append(Command('tableofcontents')) #adds TOC

doc.append(sec_intro) # adding certain sections
doc.append(sub_problem_outline)
doc.append(sub_obj_assum)
doc.append(sub_SOT)
doc.append(sub_methodology)
doc.append(sec_formula)

doc.append(sub_model_code)
doc.append(sec_math_desc)
doc.append(sec_simulation)

doc.append(sec_picture)
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
doc.generate_pdf(clean_tex=False)
"""

        #         preliminary_str=(f"""
        # =======
        preliminary_str = f"""

Examplary setup is as follows:

#References to guide imports needed to include - to see the list of available guides insert and run the code below:
```python
from dynpy.utilities.creators import list_of_guides
list_of_guides()
```

## CELL 1
## Imports

```python
{imports_str}
```



## CELL 2
## Thesis introduction

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{thesis_introduction_str}
```



## CELL 3
## Math

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{math_str}
```


## CELL 4
## Picture


```python

#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#
{picture_str}

```


## CELL 5
## Simulation

Basics of ODESystem based simulations are covered in guide to ODESystem, use the following call:

```python
from dynpy.utilities.documents.guides import BasicsOfODESystemGuide,UsageOfDynamicSystemsGuide
BasicsOfODESystemGuide();
```


```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#


{simulaion_str}

```


## CELL 6
## Veryfication

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{veryfication_str}

```

## CELL 7
## Conclusion

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{conclusion_str}

```

## CELL 9
## Symbols description

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{symbols_description_str}

```

## CELL 10
## Document

```python
#!!! BE SURE ALL PREVIOUS CELLS ARE RUN !!!#
#!!!       BECAUSE OF NEEDED IMPORTS    !!!#

{document_str}
```

"""

        display(IPMarkdown(preliminary_str))

        preliminary_str = f"""
#Example:

#To prepare a simple document with text and images:

#Good practice here is to allocate 1 section per 1 cell

## CELL 1
## Imports

{imports_str}

## CELL 2
## Thesis introduction

{thesis_introduction_str}

## CELL 3
## Math

{math_str}

## CELL 4
## Picture

{picture_str}

## CELL 5
## Simulation

{simulaion_str}

## CELL 6
## Veryfication

{veryfication_str}

## CELL 7
## Conclusion

{conclusion_str}

## CELL 9
## Symbols description

{symbols_description_str}

## CELL 10
## Document

{document_str}

"""

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
    articles_bib = """@INPROCEEDINGS{brief_history_of_simulations,
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
"""

    def _create_directory(path: str) -> None:
        # """
        # Creates a directory at the specified path.

        # Args:
        #     path (str): The directory path to be created. If subdirectories are needed     #     the path should be structured like 'directory/subdirectory'.

        # """

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

    def _create_file(name: str, content: str = None, path: str = None) -> None:
        # """
        #     Creates a file and writes content to it.

        #     Args:
        #         name (str): The name of the file to be created.
        #         content (str, optional): The content to write into the file.
        #         path (str, optional): The directory where the file should be created.
        #                               If specified, it will ensure the directory exists.
        # """
        import os

        if path is not None:
            if not os.path.exists(path):
                _create_directory(path)

            with open(os.path.join(path, name), "w") as file:
                file.write(content)

        else:
            with open(name, "w") as file:
                file.write(content)

        file.close()






class BeamerPresentation(Document):

    latex_name = "document"
    packages = [
        #Package("microtype"),
        # Package('polski',options=['MeX']),
        # Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
        Package("listings"),
        # Package('titlesec'),
        # Package('fancyhdr'),
        # Command('pagestyle', arguments=['fancy']),
        Command(
            "author", arguments=["Bogumił Chiliński"]
        ),  # ['Szymon Kozłowski & Bogumił Chiliński']),
        # Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
        # Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
        # Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
        # Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
        Command("usetheme", arguments=["Madrid"]),
        Command("graphicspath", arguments=[NoEscape("{../}")]),
    ]

    @classmethod
    def base_setup(cls):

        preliminary_str = """

# Examplary setup is as follows:

from dynpy.utilities.documents.document import BeamerPresentation
from dynpy.utilities.report import*
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display


#CELL_1
frame_1=Frame(title='Introduction',options=[])
CurrentContainer(frame_1)
display(ReportText('abc '*100))


#CELL_2
frame_2=Frame(title='Physical model',options=['allowframebreaks'])
CurrentContainer(frame_2)
display(ReportText('abc '*100))


#CELL_3
from sympy import *

frame_3=CodeFrame(title='Simulation results',options=[])
CurrentContainer(frame_3)

display(ObjectCode('from dynpy.utilities.report import *'))

#CELL_4
from sympy import *

frame_4=Frame(title='Simulation results',options=[])
CurrentContainer(frame_4)

display(SympyFormula(Symbol('a')**2))

#CELL_5
doc = BeamerPresentation('Example_beamer',title='Exemplary presentation')
doc.append(frame_1)
doc.append(frame_2)
doc.append(frame_3)
doc.append(frame_4)
doc.generate_pdf(clean_tex=False)


"""
        return ObjectCode(preliminary_str)

    def __init__(
        self,
        default_filepath="default_filepath",
        *,
        title=None,
        author=None,
        initials=None,
        documentclass="beamer",
        document_options=None,
        fontenc="T1",
        inputenc="utf8",
        font_size="footnotesize",
        lmodern=False,
        textcomp=True,
        microtype=True,
        page_numbers=True,
        indent=None,
        geometry_options=None,
        data=None,
    ):

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
            self.packages.append(Command("title", arguments=[self.title], options=[""]))

        if self.author is not None:
            if self.initials is not None:
                self.packages.append(
                    Command("author", arguments=[self.author], options=[self.initials])
                )
            else:
                self.packages.append(
                    Command("author", arguments=[self.author], options=[""])
                )

        self.append(Command("frame", arguments=[NoEscape(r"\titlepage")]))


class BeamerTemplate(BeamerPresentation):

    pass



class DynPyBasicsPresentation(BeamerPresentation):

    @classmethod
    def base_setup(cls):

        preliminary_str = r"""

### CELL 1 - IMPORTS

from dynpy.utilities.documents.document import BeamerPresentation
from dynpy.utilities.report import *
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display

# Fix Sympy formula formatting globally
SympyFormula._break_mode = 'eq'  # necessary for correct equation display

### CELL 2 - LIBRARY
slide1 = CodeFrame(title='Library Import',options=[])

CurrentContainer(slide1)
display(ReportText('''
DynPy is an advanced Python library. To use it, you need to import the modules as shown in the example below. Thanks to this operation, the user will be able to use the definitions of basic classes and methods available in the library. In this lab, the available methods for text reporting will be discussed. Finally, a sample report will be generated in the form of a PDF file.'''))

display(ObjectCode('''import dynpy
import inspect
from datetime import datetime

from dynpy.solvers.linear import ODESystem, AnalyticalSolution

from dynpy.utilities.creators import ModuleStructure, ClassLister
from dynpy.utilities.report import *
from dynpy.utilities.components.guides.github.en import ReportComponent

from pylatex import Subsubsection'''))

### CELL 3 - SECTIONS AND FRAMES (DOMCUMENTS AND PRESENTATIONS)
slide2 = CodeFrame(title='Sections and Containers',options=[])

CurrentContainer(slide2)
display(ReportText('Before you can report anything, you need to create a section and a container to hold the items within that section. You can do this as shown in the example below.'))
display(ReportText('All content placed under CurrentContainer(sec1) will automatically be included in section sec1.'))


display(ObjectCode('''sec_1 = Section(title='My first section', numbering=False)
CurrentContainer(sec_1)'''))

display(ReportText('In Python, Sections are used for creating structured reports, while Frames and CodeFrames are designed for building interactive presentations.'))

display(ObjectCode('''frame_1 = Frame(title='My first frame')
CurrentContainer(frame_1)'''))

display(ObjectCode('''codeframe_1 = CodeFrame(title='My first code frame')
CurrentContainer(codeframe_1)'''))

### CELL 4 - REPORTTEXT
slide3 = CodeFrame(title='Text Reporting Methods: ReportText', options=[])

CurrentContainer(slide3)
display(ReportText('''
Depending on the user's preferences and needs, there are several classes available for text reporting. The main purpose of all these classes is to allow placing text in the final file according to the user's specific preferences and requirements.
'''))



### CELL 5 - REPORTTEXT PREVIEW

slide31 = CodeFrame(title='Text Reporting Methods: ReportText Preview', options=[])

CurrentContainer(slide31)
display(ReportText(r'\hspace*{-1cm}'))

display(Picture('./sample_reporttext.png', width='13cm'))


### CELL 6 - REPORTTEXT RESULT

slide32 = Frame(title='Text Reporting Methods: ReportText Result', options=[])

CurrentContainer(slide32)
display(ReportText(
'''

This is a basic and also the most commonly used method. It allows the user to format text within a declared string without the need to use special commands, as is done in LaTeX or Markdown environments. Below is an example of how this class can be used:

Test Text:

- aaaaa

- bbbbb

- ccccc
'''))


### CELL 7 - MARKDOWN

slide4 = CodeFrame(title='Text Reporting Methods: Markdown', options=[])

CurrentContainer(slide4)

#display(ReportText('''Another method of reporting text is by using the Markdown class. It allows text to be reported using the syntax of the markup language of the same name. A major advantage of this class is its high level of text formatting capabilities. Markdown is a very popular language that continues to evolve, which means the possibilities of this class are nearly endless. The biggest drawback, however, is the need to understand the syntax of this language, which can be problematic and unintuitive for many users unlike the previously described ReportText class. Below is an example of how this class can be used:
#'''))
display(Markdown(
'''
Another method of reporting text is by using the Markdown class. \n
It allows text to be reported using the syntax of the markup language of the same name. \n

\+ high level of text formatting capabilities

\+ continues to evolve, which means the possibilities of this class are nearly endless

\- need to understand the syntax of this language (problematic and unintuitive in compto ReportText)

Below is an example of how this class can be used:
'''))

### CELL 8 - MARKDOWN CONT.

slide41 = CodeFrame(title='Text Reporting Methods: Markdown cont.', options=[])

CurrentContainer(slide41)

display(ReportText('''
abc abc
'''))



### CELL 9 - MARKDOWN PREVIEW

slide42 = CodeFrame(title='Text Reporting Methods: Markdown Preview', options=[])

CurrentContainer(slide42)

display(Picture('./sample_markdown.png', width='11cm'))

### CELL 10 - MARKDOWN RESULT

slide43 = Frame(title='Text Reporting Methods: Markdown Results', options=[])

CurrentContainer(slide43)

display(Markdown('''
Lorem Ipsum is simply dummy text of the printing and typesetting industry.

Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown...

Test Text:

- aaaaa

- bbbbb

- ccccc

Now is the time for all good men to come to
the aid of their country. This is just a
regular paragraph.

> This is a blockquote.
>
> This is the second paragraph in the blockquote.'''))


### CELL 11 - OBJECTCODE

slide5 = CodeFrame(title='Code reporting', options=[])

CurrentContainer(slide5)


display(ReportText('''
abc
'''))

### CELL 12 - MATH BASICS

from dynpy.utilities.report import *
from sympy import symbols, Eq, solve

slide51 = Frame(title='Code reporting: result', options=[])

CurrentContainer(slide51)

a = 10
b = 20

sum_ab = a + b

display(ReportText('Sum of a and b is equal to'))
display(SympyFormula(sum_ab))

display(ReportText('Result in the of the expression is as follows:'))
display(SympyFormula(Eq(Symbol('\\Sigma'),sum_ab)))

### CELL 13 - SYMPY

slide6 = CodeFrame(title='SympyFormula', options=[])

CurrentContainer(slide6)


display(ObjectCode('''sec_5 = Section(title='SympyFormula', numbering=False)

CurrentContainer(sec_5)

from sympy import symbols, Eq, solve

x = symbols('x')
equation = Eq(2*x + 3, 7)

display(SympyFormula(equation))'''))

from sympy import symbols, Eq, solve

x = symbols('x')
equation = Eq(2*x + 3, 7)

display(SympyFormula(equation))

### CELL 14 - IMAGES

slide7 = CodeFrame(title='Images', options=[])

CurrentContainer(slide7)


display(ObjectCode('''sec_6 = Section(title="Picture", numbering=False)

CurrentContainer(sec_6)

display(Picture('./input/logo.png', caption='Caption of the image.'))'''))
display(Picture('./logo.png', caption='Caption of the image.', width='5cm'))

### CELL 15 - TABLES

slide8 = CodeFrame(title='Tables', options=[])

CurrentContainer(slide8)


display(ObjectCode('''sec_7 = Section(title="Table", numbering=False)

CurrentContainer(sec_7)

import pandas as pd

# Sample data
df = pd.DataFrame({
    "Students": ["Anna", "Peter", "Mark"],
    "Age": [28, 34, 22],
    "City": ["Warsaw", "Cracow", "Gdansk"],
    "Grade": ["4", "5", "3.5"]
})

display(LatexDataFrame.formatted(df).reported(caption="Students Grades"))'''))

### CELL 16 - TABLES PREVIEW

slide81 = CodeFrame(title='Tables: Preview', options=[])

CurrentContainer(slide81)
display(Picture('./sample_table.png', width='10cm'))

### CELL 17 - TABLES RESULT

slide82 = CodeFrame(title='Tables: Result', options=[])

CurrentContainer(slide82)

import pandas as pd

# Sample data
df = pd.DataFrame({
    "Students": ["Anna", "Peter", "Mark"],
    "Age": [28, 34, 22],
    "City": ["Warsaw", "Cracow", "Gdansk"],
    "Grade": ["4", "5", "3.5"]
})

display(LatexDataFrame.formatted(df).reported(caption="Students Grades"))

### CELL 18 - PLOT

slide9 = CodeFrame(title='Plots', options=[])

CurrentContainer(slide9)

display(ObjectCode('''import matplotlib.pyplot as plt
from dynpy.utilities.report import PltPlot

sec_8 = Section(title="Plots", numbering=False)

CurrentContainer(sec_7)

plt.plot([0, 1, 2], [0, 1, 4], label='linear')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Sample Plot')
plt.legend()

display(PltPlot(caption='Sample plot'))'''))

### CELL 19 - PLOT RESULT

slide91 = Frame(title='Plots Result', options=['allowframebreaks'])

CurrentContainer(slide91)
import matplotlib.pyplot as plt
from dynpy.utilities.report import PltPlot

plt.plot([0, 1, 2], [0, 1, 4], label='linear')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Sample Plot')
plt.legend()

display(PltPlot(caption='Sample plot'))


### CELL 20 - DOC FORMAT

slide10 = CodeFrame(title='Document Fromatting', options=[])

CurrentContainer(slide10)

display(ReportText('''An important element of document creation is proper formatting. For this purpose, the DynPy library uses a structure compatible with the LaTeX format. According to this structure, sections are defined in the code, which correspond to chapters and subsections, to which user-defined elements are added in the form of code created using the previously described classes. The container definition is demonstrated in the example below.'''))

display(ObjectCode('''
sec_1 = Section(title='Document Formatting', numbering=False)
CurrentContainer(sec_9)
sec_91 = Subsection(title='Second-level Chapter', numbering=False) # Subchapter definition
CurrentContainer(sec_91)
display(ReportText('This is second-level section.'))

sec_911 = Subsubsection(title='Third-level Chapter', numbering=False) # Sub-subchapter definition
CurrentContainer(sec_911)
display(ReportText('This is t-level section.'))'''))

### CELL 21 - DOC FORMAT RESULT

slide101 = Frame(title='Document Formating', options=[])

CurrentContainer(slide101)
display(Picture('./sample_doc.png', width='11cm'))

### CELL 22 - LAST SLIDE

slide11 = CodeFrame(title='', options=[])

CurrentContainer(slide11)

display(ReportText(r'''
\vbox{}
{\usebeamercolor[fg]{titlegraphic}\inserttitlegraphic\par}

\begin{center}
    \Huge Thank you for your attention!
\end{center}

\vfill

'''))

### CELL 23 - PRESENTATION INFO

doc = BeamerPresentation('Example_dynpy',title='Presentation of basics of DynPy')
doc.preamble.append(Command('author', arguments=['XYZ']))
doc.preamble.append(Command('date', arguments=['Month XX, XXXX']))
doc.append(slide1)
doc.append(slide2)
doc.append(slide3)
doc.append(slide31)
doc.append(slide32)
doc.append(slide4)
doc.append(slide41)
doc.append(slide42)
doc.append(slide43)
doc.append(slide5)
doc.append(slide51)
doc.append(slide6)
doc.append(slide7)
doc.append(slide8)
doc.append(slide81)
doc.append(slide82)
doc.append(slide9)
doc.append(slide91)
doc.append(slide10)
doc.append(slide101)
doc.append(slide11)

doc.generate_pdf(clean_tex=False)


"""
        return ObjectCode(preliminary_str)

    
    
class ExternalDataBasics(BeamerPresentation):
    
    csv_data = [
                ["x", "f(x) = x^2", "Area"],
                [0, 0, 0],
                [0.1, 0.01, 0.0005],
                [0.2, 0.04, 0.0025],
                [0.3, 0.09, 0.0065],
                [0.4, 0.16, 0.0125],
                [0.5, 0.25, 0.0205],
                [0.6, 0.36, 0.0305],
                [0.7, 0.49, 0.0425],
                [0.8, 0.64, 0.0565],
                [0.9, 0.81, 0.0725],
                [1.0, 1.0, 0.0905]
                ]
    
    @classmethod
    def create_input_file(cls, data):
        headers = data[0]
        rows = data[1:]
        pd.DataFrame(rows, columns=headers).to_csv('./Integral.csv')

    @classmethod
    def base_setup(cls):
        cls.create_input_file(cls.csv_data)
        
        preliminary_str = r"""
### Woring with external data
# One of the fundamental functionalities of any library is the import of external data. In the DynPy library, the Pandas library is used for working with tabular data. It is one of the most popular libraries for data analysis and manipulation, serving as an alternative to tools such as Excel.

# To use it, it must be imported into the document, as shown in the example below:
    
    
### CELL 1 - Library Import
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display

slide1 = CodeFrame(title='Library Import',options=[])
CurrentContainer(slide1)

display(ObjectCode('''
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display
'''))

# CELL 2 - import the file
df = pd.read_csv('./Integral.csv')

display(df)

slide2 = CodeFrame(title='Reading csv file',options=[])
CurrentContainer(slide2)

display(ReportText('To read a csv file use read_csv function from Pandas library, as shown below.'))
display(ObjectCode('''
df = pd.read_csv('./Integral.csv')
'''))
display(LatexDataFrame.formatted(df).reported(caption="Data from exported file"))


### CELL 3 - loc function
#To display only a part of the table, simply use the loc() function.
display(df.loc[[0, 10]])


### CELL 4 - iloc function
#If only one column of data is being analyzed, the reference to it looks as follows:

### CELL 5 - Adding new column
# To add a new column, you can use the assign() function.
display(df.assign(new_column=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]))


### CELL 6 - Basic Operations
#Basic functions of spreadsheets include calculating the mean, sum, and finding the minimum and maximum values within a given range.
#Pandas allows you to easily obtain these statistics.
display(df.mean())
display(df.max())
display(df.min())
display(df.sum())

# Pandas allows you to combine described operations
display(df.iloc[:,1].max())


### CELL 7 
# Using these methods, we can calculate the values of a definite integral and add it as a new column. To do this, we first need to create partial sums:
sum = []
temp_sum = 0

for val in df.iloc[:,2]:
    temp_sum += val
    sum.append(temp_sum)
    
   
### CELL 8 - Adding sum column into DataFrame
df = df.assign(Sum = sum)
display(df)


### CELL 9 - Plots
display(df.plot())


### CELL 10 - Saving results
df.to_csv('./dataframe.csv', index=False)


### CELL 11 - Creating presentation
doc = BeamerPresentation('ExternalDataBasics_example',title='Working with external files', author='Author')
doc.preamble.append(Command('date', arguments=['Month XX, XXXX']))

doc.append(slide1)
doc.append(slide2)

doc.generate_pdf(clean_tex=False)
"""
        return ObjectCode(preliminary_str)
    
    
class ExternalDataBasics(BeamerPresentation):
    
    csv_data = [
                ["x", "$f(x) = x^2$", "Area"],
                [0, 0, 0],
                [0.1, 0.01, 0.0005],
                [0.2, 0.04, 0.0025],
                [0.3, 0.09, 0.0065],
                [0.4, 0.16, 0.0125],
                [0.5, 0.25, 0.0205],
                [0.6, 0.36, 0.0305],
                [0.7, 0.49, 0.0425],
                [0.8, 0.64, 0.0565],
                [0.9, 0.81, 0.0725],
                [1.0, 1.0, 0.0905]
                ]
    
    @classmethod
    def create_input_file(cls, data):
        headers = data[0]
        rows = data[1:]
        pd.DataFrame(rows, columns=headers).to_csv('./Integral.csv')

    @classmethod
    def base_setup(cls):
        cls.create_input_file(cls.csv_data)
        
        preliminary_str = r"""
### Woring with external data
# One of the fundamental functionalities of any library is the import of external data. In the DynPy library, the Pandas library is used for working with tabular data. It is one of the most popular libraries for data analysis and manipulation, serving as an alternative to tools such as Excel.

# To use it, it must be imported into the document, as shown in the example below:
    
    
### CELL 1 - Library Import
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display

slide1 = CodeFrame(title='Library Import',options=[])
CurrentContainer(slide1)

display(ObjectCode('''
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display
'''))

# CELL 2 - import the file
df = pd.read_csv('./Integral.csv', index_col=0)

display(df)

slide2 = CodeFrame(title='Reading csv file',options=[])
CurrentContainer(slide2)

display(Markdown('To read a csv file use `read_csv` function from Pandas library, as shown below.'))
display(ObjectCode('''
df = pd.read_csv('./Integral.csv', index_col=0)
'''))
display(LatexDataFrame.formatted(df).reported(caption="Data from exported file", index=False))


### CELL 3 - loc function
#To display only a part of the table, simply use the loc() function.
display(df.loc[[0, 10]])


### CELL 4 - iloc function
#If only one column of data is being analyzed, the reference to it looks as follows:

### CELL 5 - Adding new column
# To add a new column, you can use the assign() function.
display(df.assign(new_column=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]))


### CELL 6 - Basic Operations
#Basic functions of spreadsheets include calculating the mean, sum, and finding the minimum and maximum values within a given range.
#Pandas allows you to easily obtain these statistics.
display(df.mean())
display(df.max())
display(df.min())
display(df.sum())

# Pandas allows you to combine described operations
display(df.iloc[:,1].max())


### CELL 7 
# Using these methods, we can calculate the values of a definite integral and add it as a new column. To do this, we first need to create partial sums:
sum = []
temp_sum = 0

for val in df.iloc[:,2]:
    temp_sum += val
    sum.append(temp_sum)
    
   
### CELL 8 - Adding sum column into DataFrame
df = df.assign(Sum = sum)
display(df)


### CELL 9 - Plots
display(df.plot())


### CELL 10 - Saving results
df.to_csv('./dataframe.csv', index=False)


### CELL 11 - Creating presentation
doc = BeamerPresentation('ExternalDataBasics_example',title='Working with external files', author='Author')
doc.preamble.append(Command('date', arguments=['Month XX, XXXX']))

doc.append(slide1)
doc.append(slide2)

doc.generate_pdf(clean_tex=False)
"""
        return ObjectCode(preliminary_str)
    
    
    
class SympyIntegration(BeamerPresentation):
    
    csv_data = [
                ["$x$", "$f(x) = x^2$", "$Area$", "$Sum$"],
                [0.0, 0.0, 0.0, 0.0],
                [0.1, 0.01, 0.0005, 0.0005],
                [0.2, 0.04, 0.0025, 0.003],
                [0.3, 0.09, 0.0065, 0.0095],
                [0.4, 0.16, 0.0125, 0.022],
                [0.5, 0.25, 0.0205, 0.042499999999999996],
                [0.6, 0.36, 0.0305, 0.073],
                [0.7, 0.49, 0.0425, 0.11549999999999999],
                [0.8, 0.64, 0.0565, 0.172],
                [0.9, 0.81, 0.0725, 0.2445],
                [1.0, 1.0, 0.0905, 0.33499999999999996]
                ]
    
    @classmethod
    def create_input_file(cls, data):
        headers = data[0]
        rows = data[1:]
        pd.DataFrame(rows, columns=headers).to_csv('./dataframe.csv')

    @classmethod
    def base_setup(cls):
        cls.create_input_file(cls.csv_data)
        
        preliminary_str = r"""
### CELL 1 - Import Libraries
import sympy as sp
import numpy as np
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display

slide1 = CodeFrame(title='Import Libraries',options=[])
CurrentContainer(slide1)

display(ObjectCode('''
import sympy as sp
import numpy as np
import pandas as pd
from dynpy.utilities.report import *
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display
'''))


### CELL 2 - Define Symbols and Limits
x = sp.symbols('x')
f = x**2

a, b = 0, 1

slide2 = CodeFrame(title='Define Symbols and Limits',options=[])
CurrentContainer(slide2)

display(ObjectCode('''
x = sp.symbols('x')
f = x**2

a, b = 0, 1
'''))


### CELL 3 - Symbolic Integration
integral_symbolic = sp.integrate(f, (x, a, b))
print(f"Symbolic integral result = {integral_symbolic} or {integral_symbolic.evalf()}")


### CELL 4 - Numerical Integration with NumPy (Trapezoidal Rule)
num_points = 100
x_values = np.linspace(a, b, num_points)
f_values = x_values**2

# Results comparison
integral_numeric = np.trapz(f_values, x_values)
print(f"Numerical integral (trapezoidal) result = {integral_numeric}")


### CELL 5 - NumPy Error Calculation
difference_numpy = abs(float(integral_symbolic) - integral_numeric)
print(f"Difference between symbolic and numerical (Numpy): {difference_numpy}")


### CELL 6 - Load Results from Excel (Previous Class)
df = pd.read_csv('./dataframe.csv')
display(df)
# Getting value of integral from 0 to 1
integral_excel = df['$Sum$'][10]

### CELL 7 - Excel Error Calculation
difference_excel = abs(float(integral_symbolic) - integral_excel)
print(f"Difference between symbolic and numerical (Excel): {difference_excel}")


### CELL 8 - Comparison of Results
if difference_numpy < difference_excel:
    print('Numpy wins!')
elif difference_numpy > difference_excel:
    print('Excel forever!')
else:
    print("Woow! It's a draw")
    
    
### CELL 9 - Save Results
solutions_df = pd.DataFrame(
    {'Symbolic Solution': [integral_symbolic.evalf()],
     'Numpy Solution': [integral_numeric],
     'Excel Solution': [integral_excel],
    })

display(solutions_df)
solutions_df.to_csv('./results.csv', index=False)


### CELL 10 - Creating presentation
doc = BeamerPresentation('SymPyIntegration_example',title='Integration with SymPy', author='Author')
doc.preamble.append(Command('date', arguments=['Month XX, XXXX']))

doc.append(slide1)
doc.append(slide2)

doc.generate_pdf(clean_tex=False)
"""
        return ObjectCode(preliminary_str)
    
    
    
class SimulationBasics(BeamerPresentation):
    
    @classmethod
    def base_setup(cls):
        
        preliminary_str = r"""
### CELL 1 - Import Libraries
from sympy import Eq, Symbol, symbols, Function
from dynpy.utilities.report import *
from dynpy.solvers.linear import ODESystem
from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display

slide1 = CodeFrame(title='Import Libraries',options=[])
CurrentContainer(slide1)

display(ObjectCode('''
from sympy import Eq, Symbol, symbols, Function
from dynpy.utilities.report import *
from dynpy.solvers.linear import ODESystem
from dynpy.models.mechanics import ForcedSpringMassSystem as DynamicSys
from dynpy.utilities.documents.document import BeamerPresentation
SympyFormula._break_mode = 'eq' # neccesarry for correct equation display
'''))


### CELL 2 - Analyzed object
dyn_sys = DynamicSys()
dyn_sys.as_picture()

slide2 = CodeFrame(title='Import Libraries',options=[])
CurrentContainer(slide2)

display(ObjectCode('''
dyn_sys = DynamicSys()
dyn_sys.as_picture()
'''))
display(dyn_sys.as_picture(width='3cm'))


### CELL 3 - Defining symbols
k, g, m, F = symbols('k g m F', positive=True)
t = Symbol('t')

z = Function('z')(t)


### CELL 4 - Defining equations of motion
# Once all symbols and functions are defined, we can proceed to create the equations of motion. We will use the ODESystem class from the DynPy library for this purpose.
# To define the equations of motion, you need to provide at least two arguments:
# - odes - ordinary differential equations (ODEs)
# - dvars - dependent variables with respect to which differentiation is performed
# You can get additional information about the ODESystem class using the help() method.
ode = ODESystem(odes=[m*z.diff(t, t) - F + g*m + k*z], dvars=[z])
display(ode)


### CELL 5 - Solution of the defined equation
# Defining data
data = {
F:2,
g:10,
m:5,
k:100,
}
t_span = np.linspace(0.0,1,30)

# Symbolic solution
sol_sym = ode.solution()
display(sol_sym)
sol_sym = sol_sym.subs(data)
display(sol_sym)

# Numerical solution
ode_sub = ode.subs(data)
display(ode_sub)
sol = ode_sub.numerized().compute_solution(t_span, [0.0, 0.0])
display(sol)


### CELL 6 - Using predefined objects
ode_2 = dyn_sys.eoms
display(ode_2)

# You can use function solution as well to get symbolic solution
sym_sol_2 = dyn_sys.eoms.solution()
display(sym_sol_2)
display(sym_sol_2.subs(data))

# But if you want to solve this solution numerically you have to define symbols in data in a way presented below:
F=DynamicSys.F
g=DynamicSys.g
m=DynamicSys.m
k=DynamicSys.k

data ={
F:2,
g:10,
m:5,
k:100,
}

# Numerical solution
sol_num_2 = dyn_sys.eoms.subs(data).numerized().compute_solution(t_span, [0.0, 0.0])
display(sol_num_2)


### CELL 7 - Creating presentation
doc = BeamerPresentation('Simulations_example',title='Simuations in DynPy', author='Author')
doc.preamble.append(Command('date', arguments=['Month XX, XXXX']))

doc.append(slide1)
doc.append(slide2)

doc.generate_pdf(clean_tex=False)
"""
        return ObjectCode(preliminary_str)
    
    
    
class AdvancedReporting(BeamerPresentation):
    
    @classmethod
    def base_setup(cls):
        
        preliminary_str = r"""
### CELL 1 - Import Libraries
from sympy import*
from dynpy.utilities.report import*
from dynpy.utilities.templates.document import WutThesis
import pandas as pd


### CELL 2 - AutoMarker
automarker=Section("AutoMarker class")
CurrentContainer(automarker)

#
display(ReportText("Here are some figures:"))

# Please add an image named ‘pic’ to your current folder.
pic=Picture('./pic.png', caption="Sample Image", width="5cm")

a=Symbol("a")
b=Symbol("b")

dane={a:[10,20,30,40,50],b:[1,4,6,8,10]}

table=(
    LatexDataFrame.formatted(
        data=dane, index=[1,2,3,4,5]).rename_axis(
            't',
            axis=1)).reported(caption='Table with values')

plot=table.to_pylatex_tikz().in_figure(caption="Simple plot",width="10cm")

display(pic)
display(table)
display(plot)
display(ReportText(f"And now I will refer firstly to the picture ({AutoMarker(pic)}) then I need to mention the table ({AutoMarker(table)}) as well. In this short section the plot ({AutoMarker(plot)}) was done as well. "))

### CELL 3 - Biblography
biblio_mngr=BibliographyManager(arguments='biblio.bib')

biblio_entries=NoEscape(
r'''
@article{Radomski24,
author = {Amadeusz Radomski, Damian Sierociński and Bogumił Chiliński},
year = {2024},
month = {03},
pages = {1-10},
title = {Proposition of a structural health monitoring model for a concept of an innovative variable mass pendular tuned mass damper},
volume = {25},
journal = {Diagnostyka},
doi = {10.29354/diag/185458}
}
''')

biblio_mngr.append(biblio_entries)


### CELL 4 - Citing articles from the bibliography
cite_output=Section("Citation output")
CurrentContainer(cite_output)

display(ReportText("Many researches were done and some revealed that ... \cite{Radomski24}. It mentions..."))


### CELL 5 - Appending bibliography and table of contents
doc = WutThesis('AdvancedReporting')

# TO ATTACH BIBLIOGRAPHY
doc.preamble.append(Package('biblatex',["backend=biber","sorting=none"]))
doc.preamble.append(Command('addbibresource','biblio.bib'))

doc.packages.append(Package('siunitx'))

# TO ATTACH TABLE OF CONTENTS 
doc.append(Command('tableofcontents'))

doc.append(automarker)
doc.append(cite_output)

#TO DISPLAY BIBLIOGRAPHY
doc.append(biblio_mngr)
doc.append(Command('printbibliography',options=[NoEscape("title={Bibliography}")]))

doc.generate_pdf(clean_tex=False)
"""
        return ObjectCode(preliminary_str)