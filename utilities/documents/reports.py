from pylatex import (Document, Package, Command, NewPage, Tabularx
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )
from pylatex.base_classes import Environment
#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)

from ..components.mech import en as mech_compIntroDynPyProjectGuidlines
from ..components.guides import en as guide_comp
from ..components.ode import pl as ode_comp
from ..components.ode import en as ode_comp_en
from ..components.guides.reporting import en as reporting_comp
from ..components.guides.github import en as github_comp
from ..components.guides.systems import en as systems_comp
from ..components.guides.development import en as development_comp
from ..components.guides.pandas import en as pandas_comp

from .guides import Guide,UsageOfDynamicSystemsGuide

# from .guides import (UsageOfDynamicSystemsGuide, Guide, EngeneeringDrawingGuide, DevelopmentGuide, IntroToPandasGuide, 
#                     BasicsOfODESystemGuide, BasicsOfDynSysImplementationGuide, BasicsOfReportingGuide, ResearchProjectGuidelines,
#                     InterimProjectGuidelines,IntroDynPyProjectGuidelines, BasicsOfReportComponentImplementationGuide, 
#                     GithubSynchroGuide, IntroToCocalcGuide)
#from sympy import *
import datetime
import shutil
import os
from typing import Optional, Union




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
        from ...modf.odes.linear import LinearFirstOrder
        
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
                ode_comp.DetailsOfPredictedSolutionComponent,
                ode_comp.ParticularDerivativesComponent,
                ode_comp_en.GoverningEqnsWithConstFuncsComponent,
                ode_comp_en.SecularFuncsComponent,
                ode_comp_en.SecularConditionsComponent,
                ode_comp_en.SolutionWithConstFuncsComponent,
                ode_comp.ZerothOrderApproximatedEqComponent,
                ode_comp.FirstOrderApproximatedEqComponent,
                ode_comp.SecularTermsEquationsComponent,
                ode_comp.ZerothOrderSolutionComponent
        ]
        return comp_list

    @property
    def default_reported_object(self):
        from sympy import Symbol, sin, S
        from dynpy.solvers.perturbational import MultiTimeScaleSolution
        from dynpy.models.mechanics.gears import EquivalentGearModel
        sys = EquivalentGearModel()
        t=ivar=Symbol('t')
        z = sys.z
        delta = Symbol('delta',positive=True)
        eps=sys.eps
        nonlin_ode = MultiTimeScaleSolution(z.diff(t,2)+z*(1+delta*eps+eps*sin(t)), z, ivar=ivar, omega=S.One, order=3,eps=eps)
        return nonlin_ode