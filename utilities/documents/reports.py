# from .guides import (UsageOfDynamicSystemsGuide, Guide, EngeneeringDrawingGuide, DevelopmentGuide, IntroToPandasGuide,
#                     BasicsOfODESystemGuide, BasicsOfDynSysImplementationGuide, BasicsOfReportingGuide, ResearchProjectGuidelines,
#                     InterimProjectGuidelines,IntroDynPyProjectGuidelines, BasicsOfReportComponentImplementationGuide,
#                     GithubSynchroGuide, IntroToCocalcGuide)
# from sympy import *
import datetime
import os
import shutil
from typing import Optional, Union
from pandas import DataFrame

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

from ..components.guides import en as guide_comp
from ..components.guides.development import en as development_comp
from ..components.guides.github import en as github_comp
from ..components.guides.pandas import en as pandas_comp
from ..components.guides.reporting import en as reporting_comp
from ..components.guides.systems import en as systems_comp
from ..components.mech import en as mech_compIntroDynPyProjectGuidlines
from ..components.ode import en as ode_comp_en
from ..components.ode import pl as ode_comp
from .guides import Guide, UsageOfDynamicSystemsGuide


class DynSysOverviewReport(UsageOfDynamicSystemsGuide):

    @property
    def default_reported_object(self):

        # from ...models.mechanics.tmac import SDOFWinchSystem
        from ...models.mechanics import ForcedSpringMassSystem as DynamicSystem

        return DynamicSystem()

    @property
    def _report_components(self):

        comp_list = [
            systems_comp.DynSysOverviewUsageComponent,
            systems_comp.DynamicSystemCallComponent,
            pandas_comp.NumericalAnalysisSimulationComponent,
            pandas_comp.AnalyticalSimulationComponent,
            # systems_comp.SimulationsComponent,
            # reporting_comp.SimulationReportComponent,
            systems_comp.DynSysCodeComponent,
            github_comp.IssuePreparationComponent,
            # systems_comp.DynamicSystemCheckerComponent,
        ]

        return comp_list


class ODESystemOverviewReport(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list = [
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

        # from ...models.mechanics.tmac import SDOFWinchSystem
        from ...modf.odes.linear import LinearFirstOrder

        return LinearFirstOrder.from_reference_data()


class MSMethodOverviewReport(UsageOfDynamicSystemsGuide):
    @property
    def _report_components(self):
        comp_list = [
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
            ode_comp.ZerothOrderSolutionComponent,
        ]
        return comp_list

    @property
    def default_reported_object(self):
        from sympy import S, Symbol, sin

        from dynpy.models.mechanics.gears import EquivalentGearModel
        from dynpy.solvers.perturbational import MultiTimeScaleSolution

        sys = EquivalentGearModel()
        t = ivar = Symbol("t")
        z = sys.z
        delta = Symbol("delta", positive=True)
        eps = sys.eps
        nonlin_ode = MultiTimeScaleSolution(
            z.diff(t, 2) + z * (1 + delta * eps + eps * sin(t)),
            z,
            ivar=ivar,
            omega=S.One,
            order=3,
            eps=eps,
        )
        return nonlin_ode
class DataFrameAnalyzer:
    def __init__(self, df: DataFrame):
        self.df = df

    def generate(self):
        # Wybór kolumn numerycznych
        num_cols = self.df.select_dtypes(include=[np.number])
        
        if num_cols.empty:
            display(ReportText("Zbiór nie zawiera danych numerycznych."))
            return

        # Lista, do której będziemy zbierać zdania
        zdania = []

        # Zdanie wstępne
        zdania.append(f"Analizowany zestaw danych obejmuje {len(num_cols.columns)} zmiennych numerycznych.")

        for col in num_cols.columns:
            seria = num_cols[col]
            
            # Obliczenia
            mn = seria.min()
            mx = seria.max()
            avg = seria.mean()
            med = seria.median()
            std = seria.std()
            
            # Konstrukcja zdania dla danej kolumny (ciągła narracja)
            # Używamy formatu .2f aby liczby nie miały za dużo miejsc po przecinku
            zdanie = (
                f"Dla zmiennej {col} odnotowano wartości w przedziale od {mn:.2f} do {mx:.2f}, "
                f"przy czym średnia arytmetyczna wynosi {avg:.2f}, a mediana {med:.2f} "
                f"(odchylenie standardowe {std:.2f})."
            )
            zdania.append(zdanie)

        # Łączenie wszystkich zdań spacjami w jeden akapit
        pelny_tekst = " ".join(zdania)
        
        # Wyświetlenie wyniku
        display(ReportText(pelny_tekst))






