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
from nutree import Tree
from collections import defaultdict

class TreeFromTextGenerator:
    def __init__(self):
        self._name_counts = defaultdict(lambda: defaultdict(int))

    def _clean_name(self, line):
        # Czyści tekst ze śmieci typu "-Name composed of: ..."
        name = line.strip()
        if name.startswith('-'): name = name[1:].strip()
        name = name.replace("composed of:", "").strip()
        name = name.rstrip('.,')
        return name

    def _get_unique_name(self, parent, base_name):
        # Obsługa duplikatów dla nutree 1.1.0 (używa adresu pamięci rodzica)
        pid = id(parent)
        self._name_counts[pid][base_name] += 1
        count = self._name_counts[pid][base_name]
        return f"{base_name} ({count})" if count > 1 else base_name

    def generate(self, text):
        # Reset liczników
        self._name_counts.clear()
        
        # Dzielimy tekst na linie
        lines = text.strip().split('\n')

        tree = None
        # STOS: Przechowuje pary (wcięcie, obiekt_węzła)
        # Dzięki temu zawsze wiemy, kto jest "aktywnym" rodzicem
        stack = []

        for line in lines:
            # Pomiń puste linie
            if not line.strip(): continue

            # --- KLUCZOWE DLA TWOJEGO STRINGA ---
            # Liczymy TYLKO tabulatory. Ignorujemy spacje.
            indent = line.count('\t')
            # ------------------------------------
            
            name = self._clean_name(line)

            if tree is None:
                # Tworzymy korzeń (zakładamy, że pierwsza linia to zawsze korzeń)
                tree = Tree(name)
                stack.append((indent, tree))
            else:
                # LOGIKA STOSU (Backtracking):
                # Jeśli aktualne wcięcie jest mniejsze lub równe temu na górze stosu,
                # to znaczy, że zamknęliśmy poprzednią gałąź.
                # Zdejmujemy elementy ze stosu, aż znajdziemy "prawdziwego" rodzica (który ma mniejsze wcięcie).
                while stack and stack[-1][0] >= indent:
                    stack.pop()

                if stack:
                    # Rodzic to ten, kto został na szczycie stosu
                    parent = stack[-1][1]
                    
                    unique_name = self._get_unique_name(parent, name)
                    new_node = parent.add(unique_name)
                    
                    # Dodajemy nowy węzeł na stos (bo może on mieć własne dzieci)
                    stack.append((indent, new_node))

        return tree





