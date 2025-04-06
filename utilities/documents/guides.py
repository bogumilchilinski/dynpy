from pylatex import (Document, Package, Command, NewPage, Tabularx
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )

#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)
from ..report import Markdown, CurrentContainer, ReportText, IPMarkdown, ObjectCode,display
from ..components.mech import en as mech_comp
from ..components.guides import en as guide_comp
from ..components.ode import pl as ode_comp

from ..components.guides.reporting import en as reporting_comp
from ..components.guides.github import en as github_comp
from ..components.guides.systems import en as systems_comp
from ..components.guides.development import en as development_comp
from ..components.guides.pandas import en as pandas_comp
from typing import Optional, List, Union

#from sympy import *
import datetime


from .document import Guide,Document,ReportMethods

class ExampleTemplate(Guide):
    pass



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

        
        
class IntroToCocalcGuideV2(Guide):

    @property
    def _report_components(self):
        
        comp_list=[

#             github_comp.CocalcLoginComponent,
#              development_comp.JupyterSetUpComponent,
#             github_comp.CocalcFolderComponent,

#             github_comp.CocalcDynSysListComponent,
            reporting_comp.ReportingBasicsComponent,

        ]

        return comp_list
        
        
class IntroToCocalcGuide(Guide):

    @property
    def _report_components(self):
        
        comp_list=[

            #github_comp.CocalcLoginComponent,
            github_comp.CocalcUsageComponent,
            development_comp.JupyterSetUpComponent,
            github_comp.CocalcFolderComponent,

            github_comp.CocalcDynSysListComponent,
            reporting_comp.ReportingBasicsComponent,
            

        ]

        return comp_list

class UsageOfDynamicSystemsGuide(Guide):

    @property
    def _report_components(self):

        comp_list=[
            systems_comp.DynamicSystemsUsageIntroComponent,
            systems_comp.DynamicSystemCallComponent,
            systems_comp.DynamicSystemMethodsUsageComponent,

            pandas_comp.NumericalAnalysisSimulationComponent,
            
            systems_comp.SimulationsComponent, # entire component to rewrite
            reporting_comp.SimulationReportComponent,

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
            pandas_comp.IntroToPandasUsageComponent,
            pandas_comp.PandasTableGenerationComponent,
            pandas_comp.PandasMethodsComponent,
            reporting_comp.BasicOperationsComponent,
            systems_comp.DynamicSystemCallComponent,
            systems_comp.SimulationsComponent, #Common with *UsageOfDynamicSystemsGuide* class
            reporting_comp.DifferentSimulationsComponent,


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
            systems_comp.ODESystemsUsageIntroComponent,
            reporting_comp.BasicUsageOfODESystemComponent,
            #reporting_comp.ODEReportComponent,
            reporting_comp.ReportCompUseComponent,
            reporting_comp.ProjectileExampleComponent,
            systems_comp.ODESimulationComponent,
            systems_comp.ODENumericalSimulationsComponent

        ]

        return comp_list



class BasicsOfDynSysImplementationGuide(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[
            systems_comp.BasicsOfDynSysImplementationIntroComponent,
            systems_comp.DynSysImplementationComponent,
            systems_comp.DynamicSystemCallComponent,
            systems_comp.DynamicSystemMethodsUsageComponent,
            systems_comp.SimulationsComponent,
            systems_comp.DynSysCodeComponent,

        ]

        return comp_list


class BasicsOfReportingGuide(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[
            reporting_comp.BasicsOfReportingIntroComponent,
            reporting_comp.ReportingBasicsComponent,
            systems_comp.DynamicSystemCallComponent,
            systems_comp.SimulationsComponent,
            reporting_comp.SimulationReportComponent,
            reporting_comp.ReportingModuleIntroComponent,
            reporting_comp.LibrariesImportComponent,
            reporting_comp.DocumentComponent,
            reporting_comp.CurrentContainerComponent,
            reporting_comp.ReportTextComponent,
#             guide_comp.
            reporting_comp.PictureComponent,
            reporting_comp.SympyFormulaComponent,
            reporting_comp.DocumentGenerationComponent,
            reporting_comp.PredefinedSectionComponent,
#             guide_comp.
            pandas_comp.TablesCreationComponent,
            reporting_comp.AutomarkerIntroComponent,
            reporting_comp.CodeEmbeddingComponent,

            reporting_comp.UnitRegistryIntroComponent,
            reporting_comp.ReportFormattingGuidelinesComponent,
#             guide_comp.

        ]

        return comp_list

   
class ResearchProjectGuidelines(BasicsOfReportingGuide):

    @property
    def _report_components(self):

        comp_list=[


           development_comp.InterimTemplateComponent,
           development_comp.ModellingInPythonGuidelinesComponent,
            # systems_comp.DynamicSystemCallComponent,
            #systems_comp.SimulationsComponent,
            #reporting_comp.SimulationReportComponent,

        ]

        return comp_list 
    
class ThesisGuidelines(ResearchProjectGuidelines):

    @property
    def _report_components(self):

        comp_list=[


           development_comp.InterimTemplateComponent,
           development_comp.ModellingInPythonGuidelinesComponent,
            # systems_comp.DynamicSystemCallComponent,
            #systems_comp.SimulationsComponent,
            #reporting_comp.SimulationReportComponent,

        ]

        return comp_list


    
    
class InterimProjectGuidelines(ResearchProjectGuidelines):

    @property
    def _report_components(self):

        comp_list=[

            development_comp.InterimScheduleComponent,
            development_comp.InterimTemplateComponent,
#             systems_comp.DynamicSystemCallComponent,
#             systems_comp.SimulationsComponent,
#             reporting_comp.SimulationReportComponent,

        ]

        return comp_list
    
    @property
    def default_reported_object(self):

        #from ...models.mechanics.tmac import SDOFWinchSystem
        #from ...models.odes.linear import LinearFirstOrder
        
        return datetime.datetime(2024,7,13)
    
    

    
class IntroDynPyProjectGuidelines(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[

            development_comp.PythonBasicsGuidelinesComponent,
            development_comp.InterimTemplateComponent,
#             guide_comp.DynamicSystemCallComponent,
#             guide_comp.SimulationsComponent,
#             guide_comp.SimulationReportComponent,

        ]

        return comp_list
    
    
class BasicsOfReportComponentImplementationGuide(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[

            reporting_comp.ReportCompImplementationComponent,
            reporting_comp.ReportingCompsUsageComponent,
#             reporting_comp.ReportCompImplementationIssueComponent, #Obecnie jest problem z argumentem reported_object, dokładniej classname i sypie błędem
            reporting_comp.ReportingComponentsList

        ]

        return comp_list
    
    @property
    def default_reported_object(self):

        return None
    
class GithubSynchroGuide(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list=[
            github_comp.GithubSynchroIntroComponent,
            github_comp.GitSynchroPanelAccessComponent,
            github_comp.GitSynchroIntroComponent,
            github_comp.UsageOfGitHubInterfacesComponent,
            github_comp.UsageOfMeetingCreatorComponent,
            #github_comp.GithubIssueReportComponent, #komponent do listowania issue, chyba nie jest tutaj potrzebny do pokazania

        ]

        return comp_list
    
    @property
    def default_reported_object(self):

        default_data = {'classname':'GitSynchroPanelAccessComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':123,
                       }


        return default_data
