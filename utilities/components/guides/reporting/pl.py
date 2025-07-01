import numpy as np
import pandas as pd
from pandas import *
from sympy import *
from sympy import lambdify
from sympy.physics.mechanics import dynamicsymbols

from .....dynamics import *
from .....solvers.linear import *
from ....report import display as ip_display
from ...mechanics import ReportComponent as BaseReportComponent
from ...mechanics import *
from ...mechanics import display

miesiace_list = [
    "styczeń",
    "luty",
    "marzec",
    "kwiecień",
    "maj",
    "czerwiec",
    "lipiec",
    "sierpień",
    "wrzesień",
    "październik",
    "listopad",
    "grudzień",
]

srednie_temp_list = [-1.9, -0.8, 3.2, 9.3, 14.6, 18, 20.1, 19.5, 14.7, 9.3, 4.8, 0.5]

Eg_dzienne_list_watogodziny_na_metr2 = [
    600,
    1000,
    3000,
    3800,
    4800,
    5400,
    5300,
    4900,
    3300,
    1700,
    700,
    500,
]

Eg_dzienne_kilowatogodziny_na_metr2 = [
    0.6,
    1,
    3,
    3.8,
    4.8,
    5.3,
    4.9,
    3.3,
    1.7,
    0.7,
    0.5,
]

długosc_dnia_w_miesiacach_godziny = [
    8.3,
    10.0,
    11.8,
    13.9,
    15.7,
    16.7,
    16.3,
    14.7,
    12.7,
    10.7,
    8.8,
    7.8,
]

data_warunki_atmosferyczne = {
    "Długość dnia w miesiącu": długosc_dnia_w_miesiacach_godziny,
    "Dzienne natezenie energii [kWh/m^2]": Eg_dzienne_list_watogodziny_na_metr2,
    "Średnia temperatura": srednie_temp_list,
}

df = pd.DataFrame(index=miesiace_list, data=data_warunki_atmosferyczne)


class TitlePageComponent(Environment):

    latex_name = "titlepage"
    guide_title = "WPROWADZENIE DO OBSŁUGI ŚRODOWISKA COCALC"

    def __init__(
        self, system=None, options=None, arguments=None, start_arguments=None, **kwargs
    ):
        r"""
        Args
        ----
        options: str or list or  `~.Options`
            Options to be added to the ``\begin`` command
        arguments: str or list or `~.Arguments`
            Arguments to be added to the ``\begin`` command
        start_arguments: str or list or `~.Arguments`
            Arguments to be added before the options
        """

        self.system = system
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.system is not None:

            system = self.system

            self.append(NoEscape("\centering"))

            self.append(NoEscape("\\Huge MANUAL DO KRÓTKIEGO KURSU PYTHONA'A \n \n"))
            self.append(NoEscape(f"\\Huge {self.guide_title} \n \n"))

            self.append(Command("vspace", arguments="1cm"))

            if len(system.q) == 1:
                dof_str = "JEDNYM STOPNIU SWOBODY"
            else:
                dof_str = "WIELU STOPNIACH SWOBODY"

            if (
                system._dissipative_potential == 0
                or system._dissipative_potential is None
            ):
                damping_str = "NIETŁUMIONE"
            else:
                damping_str = "TŁUMIONE"

            # self.append(NoEscape(f'\\Large {damping_str} UKŁADY O {dof_str} \n \n'))

            self.append(Command("vspace", arguments="1cm"))

            self.append(NoEscape(f"{system._label} \n \n"))

            self.append(Command("vspace", arguments="1cm"))

            self.append(Command("MyAuthor"))
            self.append(NoEscape(f"\\par"))
            self.append(Command("vspace", arguments="1cm"))

            # self.append(Command('vspace',arguments='1cm'))
            # self.append(NoEscape(f'\\protect\\par'))
            # self.append(NewLine())
            self.append(Command("MyDate"))


dynpy_imports_code = """

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

"""


class JupyterSetUpComponent(ReportComponent):

    title = "Zakładanie Jupytera"

    def append_elements(self):

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        display(
            ReportText(
                "Posługując sie Jupyterem będziemy tworzyć cele i wykonywali w nich pewne zadania.  "
            )
        )

        pic9 = Picture(
            "./dynpy/utilities/components/guides/images/PoDodaniuCeli.png", width="9cm"
        )

        display(pic9)

        display(
            ReportText(
                'Aby dodac nową cele wystarczy kliknąć na "prostokąt" i pod nim pojawią się opcje takie jak Code czy Text'
            )
        )

        pic10 = Picture(
            "./dynpy/utilities/components/guides/images/Markdown.png", width="9cm"
        )

        display(pic10)

        display(
            ReportText(
                "Po wybraniu opcji Text będziemy mogli napisać opis zadania, które wykonujemy lub nadać sekcji nagłówek  "
            )
        )

        pic11 = Picture(
            "./dynpy/utilities/components/guides/images/Plotki-w-naglowk.png",
            width="9cm",
        )

        display(
            ReportText(
                "Aby dodać nagłowek wystarczy wpisac znak Hashtag/płotek, ważne jest to że każdy nagłówek będzie traktowany jako sekcja niższego stopnia (\# główna sekcja, \#\# podsekcja itd...)"
            )
        )

        display(pic11)

        pic12 = Picture(
            "./dynpy/utilities/components/guides/images//Nazywanie-naglowka.png",
            width="9cm",
        )

        display(pic12)

        display(
            ReportText(
                "Po wejsciu w zakładkę Edit ukażą nam się poniższe opcje, takie jak usuwanie cel."
            )
        )

        pic13 = Picture(
            "./dynpy/utilities/components/guides/images/Del-inne-opcje.png", width="9cm"
        )

        display(pic13)

        display(
            ReportText(
                "Czasami gdy Jupyter się zawiesi można spróbować go zrestartować uzywając poniższego przycisku"
            )
        )

        pic14 = Picture(
            "./dynpy/utilities/components/guides/images/Restart.jpeg", width="9cm"
        )

        display(pic14)

        pic15 = Picture(
            "./dynpy/utilities/components/guides/images/Restart-kernel.png", width="9cm"
        )

        display(pic15)

        display(
            ReportText(
                "Aby mieć pewność że w naszym kodzie uwzględnione zostały wszystkie zaimportowane biblioteki warto uruchamiając jupytera użyć opcji Run all cells pokazanej na poniższych obrazkach"
            )
        )

        pic16 = Picture(
            "./dynpy/utilities/components/guides/images/RunAll.jpeg", width="9cm"
        )

        display(pic16)

        pic17 = Picture(
            "./dynpy/utilities/components/guides/images/Run_all.png", width="9cm"
        )

        display(pic17)

        display(
            ReportText(
                "Aby elementy z których korzystamy działały, należy zaimportować potrzebne biblioteki, takie jak na załączonym obrazku oraz w formie tekstu który można skopiowac"
            )
        )

        pic18 = Picture(
            "./dynpy/utilities/components/guides/images/ImportowanieBibliotek.png",
            width="9cm",
        )

        # display(pic18)

        display(GuideCode(dynpy_imports_code))


# podstawy raportowania


class ReportingBasicsComponent(ReportComponent):

    title = "Podstawy raportowania"

    def append_elements(self):
        from dynpy.utilities.templates.document import Guide

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(
            ReportText(
                'Pierwszy krokiem jest import "szablonu" dokumentu np. klasy `Guide` oraz wywołanie metody `base_setup`. Wynik wywołania jest następujący:'
            )
        )

        display(
            ObjectCode(
                """
from dynpy.utilities.templates.document import *
Guide.base_setup()
"""
            )
        )

        display(Guide.base_setup())

        display(
            ReportText(
                "Metoda zwraca przykładowe cele gotowe do przekopiowania do notatnika."
            )
        )


# obsluga systemow dynamicznych


section_development_str = """
intro=Section('Wprowadzenie')
CurrentContainer(intro)
display(ReportText('schemat układu'))
schemat=Picture('picture_winch.png')
display(schemat)
"""

forming_equation_str = """
solution_subs = SDOFWinchSystem()._ode_system.steady_solution.subs(slownik)
eq_sol=Eq(solution_subs,0)
eq_steady_sol=Eq(steady_solution_subs,0)

rownania_ruchu=Section('Równania ruchu')
CurrentContainer(rownania_ruchu)
display(ReportText('równania ruchu: '))

display(SympyFormula(eoms_eq))
display(SympyFormula(eq_sol))
display(SympyFormula(eq_steady_sol))
"""

plot_generation_str = """

"""

report_generation_libraries = """
doc_final = MechanicalCase('./Out/Document',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)

doc_final.append(intro)
doc_final.append(rownania_ruchu)
doc_final.append(wykres)

doc_final.generate_pdf()
"""


class SimulationReportComponent(ReportComponent):

    title = "Raportowanie wykresów"

    def append_elements(self):

        system = (
            self.reported_object
        )  # it's useless in the case of permanent content - it's commented for future usage

        display(
            ReportText("Sposób tworzenia sekcji przy użyciu metody CurrentContainer:")
        )

        # sekcja 1|
        display(ObjectCode(section_development_str))

        display(
            ReportText(
                "Poniżej przedstawiony jest sposób tworzenia równań. W tym przypadku przyrównamy wyrazy (przed przecinkiem) do 0 (po przecinku)."
            )
        )
        display(ObjectCode(forming_equation_str))

        display(
            ReportText(
                "poniżej przedstawiony jest sposób tworzenia wykresu dla podanych danych: "
            )
        )
        eoms = system._eoms[0]
        slownik = system.get_random_parameters()
        eoms_eq = Eq(eoms, 0)
        solution_subs = system._ode_system.steady_solution.subs(slownik)
        steady_solution = system._ode_system.steady_solution[0]
        steady_solution_subs = steady_solution.subs(slownik)
        steady_solution_subs
        #         eq_sol=Eq(solution_subs,0)
        #         eq_steady_sol=Eq(steady_solution_subs,0)
        display(SympyFormula(solution_subs))

        #         ode_system = ODESystem.from_dynamic_system(system.linearized())
        #         param_dict = system.get_numerical_parameters()
        t_span = np.linspace(1, 100, 101)
        ic_list = [0.0, 0.0]  ### Dla układu o jednym stopniu swobody
        #         ode_solution = ode_system.subs(param_dict).steady_solution
        #         display(SympyFormula(ode_solution))

        ode_simulation = steady_solution_subs.compute_solution(
            t_span, ic_list=ic_list
        )  # ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)

        display(ReportText("Wykres: "))
        picture = ode_simulation.to_pylatex_tikz(subplots=True).in_figure(
            caption="Otrzymany wykres"
        )

        display(picture)

        display(ReportText("równania ruchu: "))

        display(SympyFormula(eoms_eq))
        display(SympyFormula(solution_subs))
        display(SympyFormula(steady_solution_subs))

        display(ObjectCode(report_generation_libraries))

        display(ReportText("Tak wygląda ostateczny raport:"))
        display(Picture("./dynpy/utilities/components/guides/images/sekcja3_w.jpg"))


winch_pandas_code = """

from dynpy.models.mechanics.tmac import SDOFWinchSystem
eoms= SDOFWinchSystem()._eoms[0]
eoms_eq_raw=Eq(eoms,0)
data=SDOFWinchSystem().get_random_parameters()
Ms0=Symbol('M_s0')
b=Symbol('b')

winch = SDOFWinchSystem()

Ms=winch.M_s
t = winch.ivar

dane={
    Ms0:4700,
    b:1500,
    SDOFWinchSystem().m0:10}



winch_up=SDOFWinchSystem().subs(Ms,Ms0-b*SDOFWinchSystem().phi.diff()).subs(data).subs(dane)
nowe_ode=winch_up._ode_system
solution=nowe_ode.solution
steady_sol=nowe_ode.steady_solution[0]
t_span = np.linspace(0,100,200)
phi=winch.q[0]



units_dict = {phi : ureg.radian,phi.diff(t):ureg.radian,t:ureg.second }
LatexDataFrame.set_default_units(units_dict)
"""


class BasicSymComponent(ReportComponent):

    title = "Symulacje - wariant podstawowy"

    def append_elements(self):

        system = (
            self.reported_object
        )  # it's useless in the case of permanent content - it's commented for future usage

        # przekopiowac tu fragment kodu z guide pamietajac o indentach
        SDOFWinchSystem = type(system)
        #         from dynpy.models.mechanics.tmac import SDOFWinchSystem
        eoms = SDOFWinchSystem()._eoms[0]
        eoms_eq_raw = Eq(eoms, 0)
        data = SDOFWinchSystem().get_random_parameters()
        Ms0 = Symbol("M_s0")
        b = Symbol("b")

        winch = SDOFWinchSystem()

        Ms = winch.M_s
        t = winch.ivar

        dane = {Ms0: 4700, b: 1500, SDOFWinchSystem().m0: 10}

        winch_up = (
            SDOFWinchSystem()
            .subs(Ms, Ms0 - b * SDOFWinchSystem().phi.diff())
            .subs(data)
            .subs(dane)
        )
        nowe_ode = winch_up._ode_system
        solution = nowe_ode.solution
        steady_sol = nowe_ode.steady_solution[0]
        t_span = np.linspace(0, 100, 200)
        phi = winch.q[0]

        units_dict = {phi: ureg.radian, phi.diff(t): ureg.radian, t: ureg.second}
        LatexDataFrame.set_default_units(units_dict)
        display(ObjectCode(winch_pandas_code))

        display(
            ReportText(
                "Analiza i symulacje klasy SDOFWinchSystem, tzw. wciągarki. Równania ruchu: "
            )
        )
        display(SympyFormula(eoms_eq_raw))
        display(
            ReportText(
                f"Po podstawieniu danych do ogólnej postaci (do równania {AutoMarker(eoms_eq_raw)}) otrzymuje się: "
            )
        )
        display(SympyFormula(Eq(winch_up._eoms[0], 0)))
        display(ReportText("Rozwiązanie szczególne: "))
        display(SympyFormula(Eq(phi, steady_sol)))
        display(ReportText("Rozwiązanie ogólne: "))
        display(SympyFormula(Eq(phi, solution[0])))
        general_sol_matrix = nowe_ode.solution.with_ics([0, 0])
        table_eq = general_sol_matrix.n().numerized().compute_solution(t_span)
        table = table_eq[[phi.diff(t)]]
        table_for_report = table.iloc[0::15].to_latex_dataframe()
        display(
            Markdown(
                f"Wyniki symulacji zostały zestawione w tabeli {AutoMarker(table_for_report)}"
            )
        )
        display(table_for_report.reported(caption="Tabela"))
        sim_plot = table.to_pylatex_tikz().in_figure(caption="Wykres")
        display(
            Markdown(
                f"Wyniki symulacji zostały przedstawione graficznie na wykresie {AutoMarker(sim_plot)}"
            )
        )
        display(sim_plot)


steady_state_code = """
steady_state = SDOFWinchSystem().eoms.solution[0]
steady_state
"""

list_code = """
Model_number = [1,2,3,4] 
inertia_list= [1,2,3,4] 
stiffness_list = [1,2,3,4] 
damping_list = [1,2,3,4] 
Excitation_list = [1,2,3,4]

system_parameters = {'Inertia':inertia_list, 'Stiffness':stiffness_list,'Damping':damping_list,'Excitation':Excitation_list}
"""

df_code = """
df = pd.DataFrame(index = Model_numbber, data = system_parameters)
df
"""

dict_code = """
data_dict_1 = df.T[1].to_dict()
data_dict_2 = df.T[2].to_dict()
data_dict_3 = df.T[3].to_dict()
data_dict_4 = df.T[4].to_dict()
eq1 = steady_state.subs(data_dict_1).subs(SDOFWinchSystem().F,100)
eq2 = steady_state.subs(data_dict_2).subs(SDOFWinchSystem().F,100)
eq3 = steady_state.subs(data_dict_3).subs(SDOFWinchSystem().F,100)
eq4 = steady_state.subs(data_dict_4).subs(SDOFWinchSystem().F,100)
"""


sym_gen_sol = """
    table = SDOFWinchSystem().subs(data_dict_2).subs(SDOFWinchSystem.F,100).numerized().compute_solution(t_span,[0,0]).iloc[:,0] 
    table.plot() 
"""

list_elem = """
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
"""


solution_code = """
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
pd.DataFrame(index = t_span, data = eq2_fun(t_span).plot()
"""


class DifferentSimulationsComponent(ReportComponent):

    title = "Symulacje - pozostałe warianty"

    def append_elements(self):

        system = (
            self.reported_object
        )  # it's useless in the case of permanent content - it's commented for future usage

        SDOFWinchSystem = type(system)

        display(ReportText("Wygenerowanie rozwiazania szczegolnego klasy:"))

        display(
            ObjectCode(
                steady_state_code.replace("SDOFWinchSystem", system.__class__.__name__)
            )
        )

        steady_state = (
            SDOFWinchSystem.from_random_data()
            ._ode_system.solution.with_ics([0, 0])[0]
            .subs(SDOFWinchSystem.m0, 10)
        )
        display(steady_state)

        # display(Picture('./Image/steady_sol.jpg', caption = ""))

        display(ReportText("Stworzenie listy i tabelki"))
        display(ObjectCode(list_code))
        Model_number = [1, 2, 3, 4]
        inertia_list = [1, 2, 3, 4]
        stiffness_list = [1, 2, 3, 4]
        damping_list = [1, 2, 3, 4]
        Excitation_list = [1, 2, 3, 4]

        system_parameters = {
            "Inertia": inertia_list,
            "Stiffness": stiffness_list,
            "Damping": damping_list,
            "Excitation": Excitation_list,
        }

        display(ObjectCode(df_code))

        df1 = pd.DataFrame(index=Model_number, data=system_parameters)
        display(df1)

        # display(Picture('./Image/listatab.jpg', caption = ""))

        display(
            ReportText(
                "Utworzenie wektora - zestaw wartosci czasu dla naszej symulacji "
            )
        )
        display(
            Markdown(
                """
            function = lambdify(SDOFWinchSystem.ivar, steady_state)
            t_span = np.linspace(0,100,200)"""
            )
        )

        function = lambdify(SDOFWinchSystem.ivar, steady_state)
        t_span = np.linspace(0, 100, 200)

        # display(Picture('./Image/wektor1.jpg', caption = ""))

        display(
            ReportText(
                "Utworzenie słowników z danymi i podstawienie danych w nim zawartych do naszego rownania, tak aby bylo zalezne jedynie od czasu "
            )
        )

        display(ObjectCode(dict_code))
        # display(Picture('./Image/slowniki_sym.jpg', caption = ""))

        display(
            ReportText(
                "Utworzenie listy elementów podstawiajac po kolei wartosci czasu "
            )
        )
        display(ObjectCode(list_elem))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/eq2_fun.jpg", caption=""
            )
        )

        display(ReportText("Generowanie wykresu rownania szczegolnego "))
        display(ObjectCode(solution_code))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/wykres_niefun.jpg",
                caption="",
            )
        )

        display(ReportText("Generowanie wykresu rownania ogolnego "))
        display(ObjectCode(sym_gen_sol))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/ogolne_tab.jpg", caption=""
            )
        )
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/ogolne_plot.jpg", caption=""
            )
        )


class BasicOperationsComponent(ReportComponent):

    title = "Podstawowe obekty i operacje"

    def append_elements(self):

        ss1 = pd.Series([2, 34, 5, -2, 8, 12])
        x1 = LatexDataFrame(ss1)

        ss2 = ss1 * 10
        x2 = LatexDataFrame(ss2)
        ss3 = ss1.abs()

        x3 = LatexDataFrame(ss3)
        ss4 = ss1.describe()
        x4 = LatexDataFrame(ss4)

        s5 = ss1
        s5.index = ["pierwszy", "drugi", "trzeci", "czwarty", "piąty", "szósty"]
        x5 = LatexDataFrame(s5)

        wojewodztwa_list = [
            "mazowieckie",
            "wielkopolskie",
            "lubelskie",
            "warmińsko-mazurskie",
            "zachodniopomorskie",
            "podlaskie",
            "dolnośląskie",
            "pomorskie",
            "łódzkie",
            "kujawsko-pomorskie",
            "podkarpackie",
            "małopolskie",
            "lubuskie",
            "śląskie",
            "świętokrzyskie",
            "opolskie",
        ]
        powierzchnia_lista = [
            35558,
            29826,
            25122,
            24173,
            22892,
            20187,
            19947,
            18310,
            18219,
            17972,
            17846,
            15183,
            13988,
            12333,
            11711,
            9412,
        ]
        ludnosc_lista = [
            5349114,
            3475323,
            2139726,
            1439675,
            1710482,
            1188800,
            2904207,
            2307710,
            2493603,
            2086210,
            2127657,
            3372618,
            1018075,
            4570849,
            1257179,
            996011,
        ]
        dochody_lista = [
            4464,
            3371,
            3338,
            3447,
            3649,
            3552,
            3788,
            4104,
            3438,
            3508,
            3212,
            3395,
            3187,
            3586,
            3268,
            3112,
        ]
        wydatki_lista = [
            787,
            596,
            623,
            597,
            767,
            697,
            742,
            1023,
            590,
            778,
            574,
            598,
            365,
            631,
            647,
            431,
        ]
        wojwodztwa_dane = {
            "powierzchnia": powierzchnia_lista,
            "l.osob": ludnosc_lista,
            "dochody": dochody_lista,
            "wydatki": wydatki_lista,
        }

        df = pd.DataFrame(index=wojewodztwa_list, data=wojwodztwa_dane)
        x6 = LatexDataFrame(df)

        area = Symbol("A", positive=True)
        df["zaludnienie"] = df["l.osob"] / df["powierzchnia"]
        x7 = LatexDataFrame(df)

        df21 = df.iloc[1:3, 0:3]
        x8 = LatexDataFrame(df21)

        df22 = df.loc["mazowieckie"]
        x9 = LatexDataFrame(df22)

        df.loc["Warszawa", :] = [4000, 600000, 2500, 300, 120]
        x10 = LatexDataFrame(df)

        df30 = df.set_axis(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], axis="index"
        )
        x11 = LatexDataFrame(df30)

        a, b, c, d, e = symbols("a b c d e", positive=True)

        units_dict = {a: ureg.meter, b: ureg.second, "lubuskie": ureg.meter}
        LatexDataFrame.set_default_units(units_dict)

        df31 = df.set_axis([a, b, c, d, e], axis="columns")
        x12 = LatexDataFrame(df31)

        df40 = df.rename(columns={"l.osob": "LUDNOŚĆ"}, index={"lubuskie": "LUBU..."})
        x21 = LatexDataFrame(df40)

        df50 = df.dochody.apply(lambda x: x * 2)
        x13 = LatexDataFrame(df50)

        df60 = df.applymap(lambda x: x * 2)
        x14 = LatexDataFrame(df60)

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText("**Importowanie biblioteki Pandas**"))
        display(
            ObjectCode(
                """
            import pandas as pd
            """
            )
        )

        display(ReportText("**Tworzenie serii (analogicznie do kolumn w Excelu)**"))
        display(x1.reported())

        display(
            ObjectCode(
                """
            s = pd.Series([2,34,5,-2,8,12]) 
            """
            )
        )

        display(ReportText("**Operacje na seriach**s*10 "))
        display(x2.reported())

        display(ReportText("**Funkcja abs() zmiana wartrosci na dodtanie** s.abs()"))
        display(x3.reported())

        display(
            ReportText(" **Funkcja describe() - podstawowe statystyki** s.describe() ")
        )
        display(x4.reported())

        display(
            Markdown(
                """
        **Zmiana nazw indeskow domyslnie zaczynajacych sie od 0**
        s.index =['pierwszy','drugi','trzeci','czwarty','piąty','szósty']
        """
            )
        )
        display(x5.reported())

        display(
            ReportText(
                "**DataFrame (tabele) - tworzenie DataFrame na bazie słownika** *tworzenie listy wojewodztw*"
            )
        )
        display(
            ObjectCode(
                """
            wojewodztwa_list = ['mazowieckie','wielkopolskie','lubelskie'
            'warmińsko-mazurskie','zachodniopomorskie','podlaskie'
            'dolnośląskie','pomorskie','łódzkie'
            'kujawsko-pomorskie','podkarpackie','małopolskie'
            'lubuskie','śląskie','świętokrzyskie','opolskie']
            """
            )
        )

        display(ReportText("*tworzenie listy powierzchni wojewodztw w km^2*"))
        display(
            ObjectCode(
                """
            powierzchnia_lista=[35558,29826,25122,24173,22892,20187,19947,18310,18219,17972,17846,15183,13988,12333,11711,9412]
            """
            )
        )

        display(ReportText("*tworzenie listy liczby ludności*"))
        display(
            ObjectCode(
                """
            ludnosc_lista=[5349114,3475323,2139726,1439675,1710482,1188800,2904207,2307710,2493603,2086210,2127657,3372618,1018075,4570849,1257179,996011]
            """
            )
        )

        display(ReportText("*tworzenie listy dochodow na osobe*"))
        display(
            ObjectCode(
                """
            dochody_lista=[4464,3371,3338,3447,3649,3552,3788,4104,3438,3508,3212,3395,3187,3586,3268,3112]
            """
            )
        )

        display(
            ReportText("*tworzenie listy wydatków majątkowych/inwestycyjnych na osobe*")
        )
        display(
            ObjectCode(
                """
            wydatki_lista=[787,596,623,597,767,697,742,1023,590,778,574,598,365,631,647,431]
            """
            )
        )

        display(ReportText("*Tworzenie słownika*"))
        display(
            ObjectCode(
                """
           wojwodztwa_dane={'powierzchnia':powierzchnia,'l.osob':ludnosc_lista,'doch.':dochody_lista,'wyd.':wydatki_lista}
            """
            )
        )

        display(ReportText("*Wywołanie tabeli*"))
        display(
            ObjectCode(
                """
            df = pd.DataFrame(index=wojewodztwa_list,data=wojwodztwa_dane)
            df
            """
            )
        )
        display(x6.reported())

        display(
            ReportText(
                "**Operacje na kolumnach - dodanie kolumny z obliczoną gęstością zaludnienia**"
            )
        )
        display(
            ObjectCode(
                """
            df['zalud.']=df['l.os.']/df['pow.']
            df
            """
            )
        )
        display(x7.reported())

        display(
            ReportText(
                "**Metoda iloc - za pomoca indeksów wykonujemy operacje na wierszach i kolumnach**"
            )
        )
        display(
            ObjectCode(
                """
            df.iloc[1:3,0:3]
            """
            )
        )
        display(x8.reported())

        display(ReportText("**Metoda loc**"))
        display(
            ObjectCode(
                """
            df.loc['mazowieckie']
            """
            )
        )
        display(x9.reported())

        display(ReportText("**Dodanie wiersza za pomocą metody loc**"))
        display(
            ObjectCode(
                """
            df.loc['Warszawa',:]=[4000,600000,2500,300,120]
            """
            )
        )
        display(x10.reported())

        display(ReportText("**Metoda set_axis dla wierszy - operacje na indeksach**"))
        display(
            ObjectCode(
                """
            df.set_axis([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],axis='index')
            """
            )
        )
        display(x11.reported())

        display(ReportText("**Metoda set_axis dla kolumn**"))
        display(
            ObjectCode(
                """
            a,b,c,d,e = symbols('a b c d e',positive = True)
            
            units_dict = {a:ureg.meter,b:ureg.second, 'lubuskie':ureg.meter}
            LatexDataFrame.set_default_units(units_dict)

            LatexDataFrame.formatted(df.set_axis([a,b,c,d,e],axis='columns'))
            """
            )
        )
        display(LatexDataFrame.formatted(x12).reported())

        display(ReportText("**Metoda rename - zmienia nazwy kolumni wierszy** "))
        display(
            ObjectCode(
                """
            df.rename(columns={'l.osob':'LUDNOŚĆ'},index={'lubuskie':'LUBU...'})
            """
            )
        )
        display(x21.reported())

        display(
            ReportText(
                "**Metoda apply - implemetacja funkcji na wybranej kolumnie tabeli**"
            )
        )
        display(
            ObjectCode(
                """
            df.dochody.apply(lambda x: x*2)
            """
            )
        )
        display(x13.reported())

        display(
            ReportText("**Metoda applymap - implemetacja funkcji na calej tabeli**")
        )
        display(
            ObjectCode(
                """
            df.applymap(lambda x: x*2)
            """
            )
        )
        display(x14.reported())


# ODESystem
obiekt_1 = """
x = dynamicsymbols('x')
omega = Symbol('\\omega',positive = True)
Omega = Symbol('\\Omega',positive = True)
t = Symbol('t')
"""
obiekt_2 = """
m = Symbol('m',positive = True)
c = Symbol('c',positive = True)
k = Symbol('k',positive = True)
t=Symbol('t')
x=dynamicsymbols('x')

eq1 = Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0) \n

odesys1 = ODESystem(odes = Matrix([eq1.lhs-eq1.rhs]),dvars = Matrix([x], ode_order=2)) \n

odesys1.solution
"""


class BasicUsageOfODESystemComponent(ReportComponent):

    title = "Podstawy używania klasy ODESystem"

    def append_elements(self):

        display(
            ReportText(
                "Klasa ODESystem służy do rozwiązywania i analizowania równań różniczkowych. Tworzy ona obiekt w postaci równania różniczkowego zwyczajnego (ODE) w formie ogólnej i zapewnia metody podstawowych operacji, takich jak rozwiązanie ogólne czy szczególne."
            )
        )

        display(
            ReportText(
                """\n Do prawidłowego wykonania zadania należy przyjąć następujące kroki:

        - Definiowanie potrzebnych symboli, które wykorzystamy do budowy ODESystem: 
        """
            )
        )

        display(ObjectCode(obiekt_1))
        display(
            ReportText(
                "- Definiowanie równania, które wykorzystamy do budowy ODESystem: "
            )
        )

        display(ObjectCode("""eq = Eq(x.diff(t, t) + omega**2 * x,sin(Omega*t))"""))

        display(
            ReportText(
                '- Definiowanie ODESystem, gdzie " odes " to przeniesione na jedną stronę równanie, które będziemy rozwiązywać (najlepiej podawać w macierzy), a " dvars " to zmienne, po których rozwiązujemy nasze równanie (równierz podajemy jako macierz):'
            )
        )

        display(
            ObjectCode(
                """odesys = ODESystem(odes=Matrix([eq.lhs-eq.rhs]),dvars=Matrix([x]), ode_order=2)"""
            )
        )

        display(
            ReportText(
                "- Rozwiązanie ogólne wcześniej zdefiniowanego równania metodą general solution, czyli rozwiązanie ogólne:"
            )
        )

        display(ObjectCode("""odesys.general_solution[0]"""))

        display(
            ReportText(
                "- Rozwiązanie szczególne zdefiniowanego równania metodą steady solution, rozwiązanie szczególne"
            )
        )

        display(ObjectCode("""odesys.steady_solution"""))

        display(ReportText("- Przykład wykorzystania ODESystem dla oscylatora drgań:"))

        display(ObjectCode(obiekt_2))


class ProjectileExampleComponent(ReportComponent):

    title = "Przykład użycia ODESystem na podstawie rzutu ukośnego"

    def append_elements(self):
        from dynpy.solvers.linear import ODESystem

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        m, g, v, alpha = symbols("m g v \\alpha")

        t = Symbol("t")
        x = Function("x")(t)
        y = Function("y")(t)

        dane = {m: 10, g: 9.81, alpha: pi / 6, v: 20}
        predkosc_x = Eq(x.diff(t), v * cos(alpha))
        predkosc_y = Eq(y.diff(t), v * sin(alpha) - g * t)

        ode_x_x = ODESystem(
            odes=Matrix([predkosc_x.lhs - predkosc_x.rhs]),
            dvars=Matrix([x]),
            ode_order=1,
        )
        ode_x_y = ODESystem(
            odes=Matrix([predkosc_y.lhs - predkosc_y.rhs]),
            dvars=Matrix([y]),
            ode_order=1,
        )

        ODE_x_sol = ode_x_x.subs(dane).solution.with_ics([20])
        ODE_x_sol

        display(ReportText("Implementacja biblioteki ODESystem"))
        display(
            ObjectCode(
                """
            from dynpy.solvers.linear import ODESystem
            """
            )
        )

        display(ReportText("Definiowanie zmiennych i utworzenie słownika"))
        display(
            ObjectCode(
                """
            m=Symbol('m',positive=True)
            g=Symbol('g')
            alfa=Symbol('\\alpha')
            v_0=Symbol('v_0')
            t=Symbol("t")
            x=Function("x")(t)
            y=Function("y")(t)

            dane={m:10,g:10,v_0:20,alfa:pi/4}
            """
            )
        )

        display(ReportText("Równania ruchu dla osi x oraz y"))
        display(
            ObjectCode(
                """
            v_x=Eq(x.diff(t),v_0*cos(alfa))
            v_y=Eq(y.diff(t),v_0*sin(alfa)-g*t)
            """
            )
        )
        display(ReportText("Zapisanie równań za pomocą ODESystem"))
        display(
            ObjectCode(
                """
            ODE_x=ODESystem(odes=Matrix([v_x.lhs-v_x.rhs]),dvars=Matrix([x]),ode_order=1)
            ODE_y=ODESystem(odes=Matrix([v_y.lhs-v_y.rhs]),dvars=Matrix([y]),ode_order=1)
            """
            )
        )

        display(
            ReportText(
                "Wywołanie odpowiednio zmiennych i funkcji daje następujące rezultaty"
            )
        )
        display(
            ObjectCode(
                """
            ODE_x
            """
            )
        )
        display(SympyFormula(ode_x_x))
        display(
            ObjectCode(
                """
            ODE_x.general_solution[0]
            """
            )
        )
        display(SympyFormula(ode_x_x.general_solution[0]))
        display(
            ObjectCode(
                """
            ODE_x.steady_solution[0]
            """
            )
        )
        display(SympyFormula(ode_x_x.steady_solution[0]))
        display(
            ObjectCode(
                """
            ODE_x.solution
            """
            )
        )
        display(SympyFormula(ode_x_x.solution[0]))
        display(
            ObjectCode(
                """
            ODE_x.solution.subs(dane)
            """
            )
        )

        display(SympyFormula(ode_x_x.solution.subs(dane)))
        display(
            ObjectCode(
                """
            ODE_x_sol=ODE_x.subs(dane).solution.with_ics([20])
            ODE_x_sol
            """
            )
        )

        display(SympyFormula(ODE_x_sol))
        display(
            ReportText(
                'Ten sam wynik można osiągnąć stosując metoodę " dsolve ", jednak nie jest ona zalecana do równań różniczkowych'
            )
        )
        display(
            ObjectCode(
                """
            dsolve(v_y.subs(dane))
            """
            )
        )
        display(SympyFormula(dsolve(predkosc_y.subs(dane))))


####
omega, omega2 = symbols("omega Omega")
t = Symbol("t")
x = Function("x")(t)
y = Function("y")(t)
ode = """
ODESystem(odes=Matrix([omega**2*x-sin(omega2*t)+x.diff(t,t)]),dvars=Matrix([x]),ode_order=1)
"""


class ReportCompUseComponent(ReportComponent):

    title = "Podstawy używania komponentu raportującego"

    def append_elements(self):
        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        komponent = Chapter("Podstawy używania komponentu raportującego")
        CurrentContainer(komponent)

        display(
            ReportText(
                "W celu użycia komponentu raportującego należy użyć $.report$. Przykładowo wpisując *odesys.report* otrzyma się następujący output:"
            )
        )
        display(
            ReportText(
                "The investigated system is described by differential equations being as follows:"
            )
        )
        display(ObjectCode(ode))
        display(
            ReportText(
                "To solve the problem serve several methods depending on the equation type. As dynamic systems s behaviour is described by an ordinary differential equation, the variables of the equation are as follows: t [x]"
            )
        )
        display(
            ReportText(
                "The variables allow to study and analyse the system s dynamics varying with time."
            )
        )


reportcomp_issue_title_str = """
Implementation of `{classname}` class that creates a part of a report
"""

reportcomp_issue_desc_str = """
Class is to implement in `{module}` module that creates a part (reporting component - formally child of `ReportComponent` class) of {field} for {target}.
"""


class ReportCompImplementationIssueComponent(ReportComponent):

    title = "Issue na implementację komponentów raportujących"

    def append_elements(self):

        classname = self.reported_object[
            "classname"
        ]  # it's useless in the case of permanent content - it's commented for future usage
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        display(
            ReportText(
                "Przykładowy tekst, który pozwoli na sprawne przygotowywanie issue ma następującą formę:"
            )
        )

        display(ReportText("# Title: "))
        display(ObjectCode(reportcomp_issue_title_str.format(classname=classname)))

        display(ReportText("# Description: "))
        display(
            ObjectCode(
                reportcomp_issue_desc_str.format(
                    classname=classname,
                    module=class_module,
                    field=class_field,
                    target=target,
                )
            )
        )


rep_comp_call = """
from dynpy.utilities.components.guides.en import ReportComponent
from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, GuideCode

class {classname}(ReportComponent):
"""

rep_comp_call_with_pass = rep_comp_call + "\n \t pass"


title_call_str = """
    title="Implementacja komponentów raportujących"
"""

append_elem_str = """

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
"""

rep_obj_str = """
    #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
"""


class ReportCompImplementationComponent(ReportComponent):

    title = "Implementacja komponentów raportujących"

    @property
    def reported_object(self):

        default_data = {
            "classname": "ReportingModuleIntroComponent",
            "module": "guide.en.py",
            "field": "guide or report",
            "target": "`ODESystem` class",
            "issue_no": 359,
        }

        if isinstance(self._reported_object, dict):
            return {**default_data, **self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data, "classname": self._reported_object}

        elif self._reported_object is None or not isinstance(
            self._reported_object, dict
        ):
            return default_data

        else:
            return self._reported_object

    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object = obj

    def append_elements(self):

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(
            ReportText(
                "Tworzenie komponentu należy zacząć od wybrania nazwy opisującej nasz komponent raportujący oraz zadeklarowania go w odpowiedni sposób, który został przedstawiony poniżej:"
            )
        )  # itp
        display(GuideCode(rep_comp_call_with_pass.format(classname=classname)))
        display(
            ReportText(
                "Następnie deklarujemy nazwę komponentu, która będzie wyświetlana jako nagłówek."
            )
        )
        display(GuideCode(title_call_str))

        display(
            ReportText(
                "Kolejnym krokiem tworzenia komponentu, jest zadeklarowanie metody `append_elements`, która zawiera zasadniczą część komponentu."
            )
        )
        display(GuideCode(append_elem_str))
        display(
            ReportText(
                "Dodatkowo na początku tej metody należy umieścić następujący kod, który w razie potrzeby pozwala na to aby komponent obsługiwał dowolny system dynamiczny."
            )
        )
        display(GuideCode(rep_obj_str))

        display(
            ReportText(
                "Ostatecznie kod implementacji całego komponentu ma następującą postać:"
            )
        )
        display(
            GuideCode(
                rep_comp_call.format(classname=classname)
                + title_call_str
                + append_elem_str
            )
        )


######### GUIDE REPORTING COMPONENTS
######### GUIDE REPORTING COMPONENTS
######### GUIDE REPORTING COMPONENTS

doc_gen_str = """

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

"""


class DocumentGenerationComponent(ReportComponent):

    # title="Generowanie dokumentu"
    title = "Document generation"

    @property
    def reported_object(self):

        default_data = {
            "classname": "ReportingModuleIntroComponent",
            "module": "guide.en.py",
            "field": "guide or report",
            "target": "`ODESystem` class",
            "issue_no": 359,
        }

        if isinstance(self._reported_object, dict):
            return {**default_data, **self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data, "classname": self._reported_object}

        elif self._reported_object is None:
            return default_data

        elif self._reported_object is None or not isinstance(
            self._reported_object, dict
        ):
            return default_data

        else:
            return self._reported_object

    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object = obj

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here

        # display(ReportText('Ostatni krok to zaapendowanie sekcji i utworzenie dokumentu pdf :'))
        display(
            ReportText(
                "The last step is to append a section and create a pdf document :"
            )
        )
        display(GuideCode(f"{doc_gen_str}"))


doc_comp_str = """

doc = TikzCaSCStandalone(default_filepath='./output/Tytul_Dokumentu')
Picture._default_width = NoEscape('0.8\\columnwidth')

"""


class DocumentComponent(DocumentGenerationComponent):

    title = "Tworzenie dokumentu"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        display(
            ReportText(
                "Aby stworzyć dokument i nadać mu nazwę należy wywołać kod, jak w komórce nr 2 jupytera, tj. cell [2]:"
            )
        )
        display(GuideCode(f"{doc_comp_str}"))


sim_code = """
ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)
ode_simulation.plot()
"""


class DynamicSystemCompletenessCheckComponent(ReportComponent):

    title = (
        "Wywoływanie i sprawdzanie wszystkich kluczowych elementów systemu dynamicznego"
    )

    def append_elements(self):

        system = (
            self.reported_object
        )  # it's useless in the case of permanent content - it's commented for future usage

        eoms = system._eoms[0]
        eoms

        display(
            ReportText(
                "Wywołanie systemu w pierwszej kolejności można sprawdzić przez sprawdzenie równania ruchu metodą _eoms."
            )
        )

        display(ReportText("Wynik jest następujący:"))

        display(
            Markdown(
                """
        eoms_eq=Eq(eoms,0)
        eoms_eq
            """
            )
        )

        eoms_eq = Eq(eoms, 0)
        eoms_eq

        display(SympyFormula(eoms_eq))

        display(ReportText("Następnie wykonuje się analizę metody __init__ klasy:"))

        display(ObjectCode(system.__init__))

        display(
            ReportText(
                "Kolejnym etapem jest sprawdzenie, czy schemat oraz zdjęcie reczywiste systemu są zdefiniowane:"
            )
        )

        display(
            Markdown(
                """
        system._as_picture()
        system.preview(example = 'detail_real_name')
            """
            )
        )

        system._as_picture()
        system.preview(example="detail_real_name")

        display(
            ReportText(
                "Kolejnym etapem jest sprawdzenie, czy zdefiniowane są słowniki z parametrami:"
            )
        )

        display(
            Markdown(
                """
        system.get_random_parameters()
        system.get_numerical_parameters()
            """
            )
        )

        display(
            ReportText(
                "Jeśli podczas wywołania pojawia się błąd o braku tej metody lub parametry są inne niż w metodzie __init__ należy zmodyfikować klasę analogicznie do pokazanego rozwiązania:"
            )
        )

        display(
            Markdown(
                """
        
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
            """
            )
        )

        system.get_random_parameters()
        system.get_numerical_parameters()

        display(
            ReportText(
                "Następnie konieczne jest sprawdzenie opisów symboli oraz jednostek:"
            )
        )

        display(
            Markdown(
                """
        system.symbols_description()
        system.unit_dict()
            """
            )
        )

        display(
            ReportText(
                "Jeśli podczas wywołania tej metody pojawia się błąd, należy wykonać modyfikację w analogizny sposób do:"
            )
        )

        display(
            Markdown(
                """
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

            """
            )
        )

        from sympy.physics import units

        ureg = UnitRegistry()

        system.symbols_description()
        system.unit_dict()


class AutomarkerIntroComponent(DocumentGenerationComponent):

    title = "Wstęp do wykorzystywania narzędzia AutoMarker"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        display(
            ReportText(
                "Aby w raporcie referować do wcześniej wyświetlanych równań/obrazków/tabelek używamy narzędzia AutoMarker. W przypadku równania         przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy zdefiniować potrzebne symbole oraz równanie, a następnie               wyświetlić je. Ważne jest to, aby równanie zostało zdefiniowane przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonego           równania używamy AutoMarkera. Przkładowy kod pokazany został poniżej:"
            )
        )
        display(GuideCode(f"{aut_indx_str1}"))
        display(
            ReportText(
                " W przypadku obrazka przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy zdefiniować zmienną           która zostanie przypisana do naszego obrazka stworzonego za pomocą klasy Picture (podając odpowiednią ścieżkę), a następnie wyświetlić ją. Ważne           jest to, aby zmienna została zdefiniowana przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonego obrazka ponownie używamy         AutoMarkera w ten sam sposób. Przykładowy kod pokazany został poniżej:"
            )
        )
        display(GuideCode(f"{aut_indx_str2}"))
        display(
            ReportText(
                " W przypadku tabelki przykładowy kod z użyciem tego narzędzia tworzy się następująco: Na początek należy utowrzyć zmienną która         zostanie przypisana do naszej tabelki(przykładowa tabelka pokazana sekcje wyżej), a następnie wyświetlić ją. Ważne jest to, aby zmienna została             zdefiniowana przed użyciem odniesienia. Następnie, aby zareferować do wcześniej wyświetlonej tabelki ponownie używamy AutoMarkera w ten sam sposób.         Przykładowy kod pokazany został poniżej:"
            )
        )
        display(GuideCode(f"{aut_indx_str3}"))


class PictureComponent(DocumentGenerationComponent):

    title = "Picture component implementation"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        # display(ReportText('Picture is a class that represents a figure environment'))
        display(
            ReportText(
                """Picture is a class that represents a figure environment. To add picture to document simply call Picture() class. Possible argumnts in Picture class are:
        \n path - define path to the picture
        \n position - define picture positioning acording to latex 
        \n width - define picture width in cm
        \n height - define picture height in cm
        \n caption - caption of the picture
        \n Example picture class call: \n
        """
            )
        )
        display(
            ObjectCode(
                """Picture('pic.PNG', position = 'h', caption = 'Example text', height = '10', width = '20')"""
            )
        )
        pic = Picture(
            "./dynpy/utilities/components/guides/images/1.png",
            caption="Example text",
            width="10",
        )
        display(pic)


class BasicsOfReportComponentImplementationGuide(UsageOfDynamicSystemsGuide):

    @property
    def _report_components(self):

        comp_list = [
            guide_comp.ReportCompImplementationComponent,
            guide_comp.ReportingCompsUsageComponent,
            #             guide_comp.ReportCompImplementationIssueComponent, #Obecnie jest problem z argumentem reported_object, dokładniej classname i sypie błędem
            guide_comp.ReportingComponentsList,
        ]

        return comp_list

    @property
    def default_reported_object(self):

        return None


CodeEmbeddingComponent_str1 = """

display(ObjectCode(GeometrySceneDG))

"""


CodeEmbeddingComponent_str2 = """

class GeometrySceneDG:


    ax_2d = None
    ax_3d = None

    #def __init__(self,height=12,width=9,figsize=(12,9)):
    def __init__(self, init_3d=(30, 10), height=12, width=16, figsize=(12, 9)):
        
        
        
        plt.figure(figsize=figsize)
        ax_2d = plt.subplot(121)
        ax_2d.set(ylabel=(r'<-x | z ->'), xlabel='y')

        plt.xlim(0, width)
        plt.ylim(-height, height)
        plt.grid(True)

        ax_2d.set_yticks(range(-12, 12, 1))
        ax_2d.set_yticklabels(
            list(map(lambda tick: str(abs(tick)), range(-12, 12, 1))))

        ax_2d.set_xticks(range(0, 16, 1))
        ax_2d.set_xticklabels(
            list(map(lambda tick: str(abs(tick)), range(0, 16, 1))))

        ax_3d = plt.subplot(122, projection='3d')
        ax_3d.set(xlabel='x', ylabel='y', zlabel='z')

        plt.xlim(0, 16)
        plt.ylim(0, 16)

        ax_3d.set_zlim(0, 16)

        ax_3d.view_init(*init_3d)
        plt.tight_layout()

        self.__class__.ax_2d = ax_2d
        self.__class__.ax_3d = ax_3d
        self.__class__.size_3d=16

"""


class CodeEmbeddingComponent(DocumentGenerationComponent):

    title = "Obsadzanie kodu"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        display(
            ReportText(
                "Klasa umożliwiająca dodanie sformatowanego fragmentu kodu pisanego programu. Jako przykład wykorzystana zostanie klasa **GeomerySceneDG** z biblioteki **dgeometry**."
            )
        )
        display(GuideCode(f"{CodeEmbeddingComponent_str1}"))
        display(GuideCode(f"{CodeEmbeddingComponent_str2}"))


ReportTextComponent_str = """
display(ReportText('Tutaj wpisujemy tekst, który chcemy wyświetlić w naszym dokumencie'))
    
"""


class MarkdownComponent(DocumentGenerationComponent):

    title = "Klasa Markdown"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here

        display(
            ReportText(
                "Modułem o rozszerzonej funkcjonalności jest funkcja **Markdown**. Funkcjonalność dostępna po wywołaniu kodu::"
            )
        )
        # display(GuideCode(f'{MarkdownComponent_str}'))
        display(
            Markdown(
                """

    - Tworzenie elementów wypisanych w punktach

    - Załączanie lącza internetowego

    - Wyswietlanie kodu

    """
            )
        )
        display(ReportText("Aby wymienić w punktach należy napisac następujący kod:"))
        # markdown1=Picture('./Image/markdown1.jpg')

        display(
            Markdown(
                '''
            display(Markdown("""

                    - Tworzenie elementów wypisanych w punktach

                    - Załączanie lącza internetowego

                    - Wyswietlanie kodu

            """))
        '''
            )
        )

        # display(markdown1)

        display(ReportText("Aby załączyć link, na przyklad do spotkania:"))

        display(
            Markdown(
                '''
            display(Markdown("""Link do spotkania: [SPOTKANIE]
            (https://wutwaw.sharepoint.com/sites/EfektywnyPythonbySiMR/_layouts/15/stream.aspx?id=%2Fsites%2EfektywnyPythonbySiMR%2FShared%20Documents%2FGeneral%2FRecordings%2FSpotkanie%20na%20kanale%20Ogólnym%2D20230628%5F175853%2DNagrywanie%20spotkania%2Emp4)
            """))
        '''
            )
        )

        # markdown2=Picture('./Image/mrkdwn2.jpg')

        # display(markdown2)

        display(
            ReportText(
                "Aby wyswietlic kod, przed apostrofami nalezy zapisać literę $f$. Przykładowo :"
            )
        )

        display(
            Markdown(
                '''
            display(Markdown("""

                    - Tworzenie elementów wypisanych w punktach

                    - Załączanie lącza internetowego

                    - Wyswietlanie kodu

            """))
        '''
            )
        )

        # markdown3=Picture('./Image/mrkdwn3.jpg')

        # display(markdown3)

        display(
            ReportText(
                "Ostatni krok to zaapendowanie sekcji i utworzenie dokumentu pdf :"
            )
        )
        display(
            Markdown(
                f"""


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



        """
            )
        )


class ReportTextComponent(DocumentGenerationComponent):

    title = "Klasa ReportText"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here

        display(
            ReportText(
                "Narzędziem o podstawowej funkcjonalności jest klasa **ReportText** - umożliwia dodanie i wyświetlanie tekstu, który ma zostać umieszczony w raporcie. Funkcjonalność ta jest dostępna po wywołaniu kodu:"
            )
        )
        display(GuideCode(f"{ReportTextComponent_str}"))


kod_2 = """
from dynpy.utilities.components.guides.en import ReportComponent
    from dynpy.utilities.report import ReportText, Markdown, Picture, SympyFormula, Frame, ObjectCode, Block, AlertBlock, ExampleBlock, GuideCode, LatexDataFrame
    from dynpy.utilities.components.guides.en import ReportCompImplementationComponent,ReportCompImplementationIssueComponent
    from sympy import *
    from sympy.physics import units
    import pandas as pd
"""


class LibrariesImportComponent(DocumentGenerationComponent):

    # title="Klasa LibrariesImportComponent i jej zastosowanie"
    title = "LibrariesImportComponent class and its use"

    def append_elements(self):
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]
        # display(ReportText('Przed pracą w jupiterze należy zaimportować biblioteki'))
        display(ReportText("Before working in jupiter, import libraries"))
        display(ObjectCode(kod_2))
        # display(ReportText('Dzięki temu możemy korzystać z przygotowanych wcześniej klas i modułów'))
        display(
            ReportText("This allows us to use previously prepared classes and modules")
        )


markdown_reporting_str = """
Koncepcję raportowania można przedstawić prosto w następujących punktach:  
​
* Biblioteka do raportowania ma budowę modułową. W kodowaniu obiektowo zorientowanym oznacza to, że składa się z szeregu klas. Poszczególne klasy to implementacje programistyczne standardowych elementów, z których składa się każdy raport. To znaczy zaimplementowano klasę do tworzenia dokumentu oraz klasy do kreowania typowych elementów dokumentu, jak rozdział, akapit (tj. blok tekstu), wzór, obraz, wykres.  
​
* Przyjęta została współbieżna zasada tworzenia raportów. Według tej zasady, aby zaraportować dany element najpierw musi on zostać wyświetlony w jupyterze, używając nomenklatury informatycznej - wydrukowany na wyjściu komórki jupytera, tzn. jako cell output.  
​
* Główna implementacja kodu, czyli używane klasy, znajdują się w module `dynpy.utilities.report`
​
* Do zbudowania standardowego raportu używa się następujących klas:  
    
    * class `Document`
    
    * class `CurrentContainer`
    
    * class `ReportText`
    
    * class `Markdown`
    
    * class `Picture`
    
    * class `SympyFormula`
​
* W efekcie tworzony jest dokument wynikowy w postacie pliku pdf, wygenrowany w języku zanczników LaTeX w oparciu o bazowy plik jupytera (*.ipynb), który zawiera przejrzysty podgląd dokumentu wynikowego.
​
"""


class ReportingModuleIntroComponent(DocumentGenerationComponent):

    title = "Implementacja komponentów raportujących"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        display(Markdown(markdown_reporting_str))


SympyFormulaComponent_str = """
from sympy import *

    M = Symbol('M')

    t = Symbol('t')

    m_e = Symbol('m_e')

    e = Symbol('e')

    z = Function('z')

    varphi = Function('varphi')

    Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2

    display(Eq(Symbol('T'),Ek))
"""


class SympyFormulaComponent(DocumentGenerationComponent):

    # title = "Klasa SympyFormula i jej zastosowanie"
    title = "SympyFormula class and its use"

    def append_elements(self):
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        # display(ReportText('Przykładowe zastosowanie komponentu do obliczeń symbolicznych, na przykładzie energii kinetycznej:'))
        display(
            ReportText(
                "Example application of the component for symbolic calculations, using kinetic energy as an example:"
            )
        )
        display(GuideCode(f"{SympyFormulaComponent_str}"))


report_formating_str = """
Należy zwrócić uwagę na następujęce aspekty:

- Cele powinny być poopisywane Markdownami żeby mozna bylo sie swobodnie poruszać po dokumencie przez Table of Contents.

- Appendowanie sekcji do dokumentu powinno znajdować się na końcu dokumentu, ponieważ aby zbudować dokument nie trzeba uruchamiać wszystkich cele po kolei.

- Nie można obsadzać nieformatowalnych kodów w 'ReportText'. 

- Jeżeli chcemy obsadzić jednoliniowy kod, musimy użyć ` (apostrof pod tyldą). np: display(Markdown(\"`help`\")). Ten sposób nie zapewnia kolorowania kodu.

- Obsadzamy w tekście zmienne w ten sposób: `${latex(omega_r)}$`, aby globalnie mieć do nich dostęp, możliwość podmiany symbolu bez manualnego szukania w tekście.

- Oznaczamy równania za pomocą klasy `AutoMarker`.

- Gdy mamy zbyt długie równania, szczególnie części pod pierwiastkiem, należy wyciągnąć je do zmiennej i pokazać niżej. Alternatywną metodą na rozwiązanie tego problemu jest wstawienie równania w system dynamiczny.

- Przechowywane wstępne obliczenia / symbole należy umieścić w pliku.py lub wykonać za pomocą klasy jeśli taka istnieje.

Cel pracy i metodyka powinny charakteryzować się nastęującymi cechami:

- stwierdzać fakty,

- być napisane w jednym czasie (przeszły lub teraźnieszy)

- nie uprzedzać faktów, które następują w dalszej częsci pracy,

- mieć charakter planu lub zbioru założeń/oczekiwań np. założono, że zostanie wykonana analiza wrażliwości; stwierdzono, że takie podejście zapewni poprawność modelu

Elementy, które wymagają usprawnienia w szablonie pracy:

- streszczenia w szablonie pracy dyplomowej,

- naprawienie obsługi komendy `\cite` w klasie `Markdown`

- stosowanie `{latex}`

- stosowanie klasy `AutoMarker`

- należy w szablonie PD przenieść Appendowanie na koniec 

"""

report_formating_en_str = """
Pay attention to the following aspects:

- Objectives should be described in Markdown to be able to move freely through the document via Table of Contents.

- Appending sections to the document should be at the end of the document, because it is not necessary to run all the targets one by one to build the document.

- You can't cast unformatable codes in 'ReportText'. 

- If you want to cast a one-line code, you must use ` (apostrophe under the tilde). Ex: display(Markdown(`help`)). This way does not provide code coloring.

- We staff the variables in the text like this: `${latex(omega_r)}$` to globally access them, the ability to substitute the symbol without manually searching in the text.

- We mark equations using the `AutoMarker` class.

- When we have equations that are too long, especially the parts under the root, they should be extracted into a variable and shown below. An alternative method to solve this problem is to insert the equation in the dynamic system.

- The stored precalculations / symbols should be placed in a.py file or made using a class if one exists.

The purpose of the work and the methodology should have the following characteristics:

- state facts,

- be written in one tense (past or present)

- not anticipate the facts that follow in the later part of the work,

- have the character of a plan or a set of assumptions/expectations, e.g., it was assumed that a sensitivity analysis would be performed; it was stated that this approach would ensure the correctness of the model

Elements that need improvement in the thesis template:

- abstracts in the thesis template,

- fixing the handling of the `\cite` command in the `Markdown` class.

- use of `{latex}`

- application of `AutoMarker` class.

- move Append to the end in the PD template 

"""


class ReportFormattingGuidelinesComponent(DocumentGenerationComponent):

    # title="Wytyczne odnośnie formatowania raportu"
    title = "Report formatting guidelines"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here
        display(Markdown(report_formating_en_str))


old_class_code = """class ExemplaryOldImplementedSystem(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('\\varphi')

    c0 = Symbol('c0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        self.ivar = ivar

        self.qs = [self.phi]

        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity_ = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)

        composed_system = self._mass_x + self._mass_y + self._gravity_

        super().__init__(composed_system, **kwargs)

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict
    
    def get_numerical_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict

    def max_static_force(self):
        return S.Zero

    def max_dynamic_force(self):
        return S.Zero"""
#
new_class_code = """class MaterialPointMovement(ComposedSystem):

    m = Symbol('m', positive=True)
    g = Symbol('g', positive=True)
    c = Symbol('c', positive=True)
    r = Symbol('r', positive=True)
    phi = dynamicsymbols('phi')

    c0 = Symbol('c0', positive=True)
    r0 = Symbol('r0', positive=True)
    phi0 = dynamicsymbols('phi0')

    def __init__(self,
                 m=None,
                 g=None,
                 c=None,
                 r=None,
                 phi=None,
                 ivar=Symbol('t'),
                 **kwargs):

        if m is not None: self.m = m
        if g is not None: self.g = g
        if c is not None: self.c = c
        if r is not None: self.r = r
        if phi is not None: self.phi = phi
        self.ivar = ivar

        self.qs = [self.phi]

        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}


        self._mass_x = MaterialPoint(self.m,
                                     pos1=self.r * sin(self.phi),
                                     qs=self.qs)
        self._mass_y = MaterialPoint(self.m,
                                     pos1=self.r * cos(self.phi),
                                     qs=self.qs)

        self._gravity = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)



        components['_mass_x']=self._mass_x
        components['_mass_y']=self._mass_y
        components['_gravity']=self._gravity


        return components
        
    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'Mass',
            self.g: r'Gravity constant',
            self.c: r'',
        }

        return self.sym_desc_dict

    def get_default_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict
    
    def get_numerical_data(self):

        m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

        default_data_dict = {
            self.m: [m0 * no for no in range(1, 8)],
            self.c: [c0 * no for no in range(1, 8)],
            self.r: [r0 * no for no in range(1, 8)],
            self.phi: [phi0 * no for no in range(1, 8)],
        }

        return default_data_dict    

    def max_static_force(self):
        return S.Zero

    def max_dynamic_force(self):
        return S.Zero"""


class CodeRefactorComponent(ReportComponent):

    title = "DynSys system code refactor component"

    def append_elements(self):

        from dynpy.models.mechanics import Engine

        system = self.reported_object

        display(
            ReportText(
                "Refactoring is a systematic process of improving code without creating new functionality that can transforma mess into clean code and simple design."
            )
        )

        display(ReportText("Example of old DynSys component before refactoring:"))

        display(ObjectCode(old_class_code))

        display(ReportText("Example of current DynSys component after refactoring:"))

        display(ObjectCode(new_class_code))

        display(
            ReportText(
                "After calling system_description() method for Engine model bellow result will be created:"
            )
        )

        str_sys_des = (
            system.system_description()
            .replace("\n", "\n\n")
            .replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;")
        )

        display(ReportText(str_sys_des))

        display(
            ReportText(
                "\n\nCalling system_parameters() method will create bellow list:"
            )
        )

        lst_sys_par = system.system_parameters()

        for element in lst_sys_par:

            display(ReportText(f"* {element}"))


container_code = """
nazwa_sekcji=Section('Tytuł')
CurrentContainer(nazwa_sekcji)
"""

unit_registry_str = """        


t=Symbol('t')
r=Symbol('r', positive=True)
a=Symbol('a', positive=True)
b=Symbol('b', positive=True)
h=Symbol('h', positive=True)
vphi = Function('varphi')(t)

dane={
    vphi.diff():25,
    r:0.15,
    a:0.2,
    b:0.6,
    h:0.4,
}

unit_dict = {
            t:ureg.second,
            r: ureg.meter,
            h: ureg.meter,
            a: ureg.meter,
            b:ureg.meter,
            vphi.diff(t):ureg.rad/ureg.sec,
   }

LatexDataFrame.set_default_units(unit_dict)

tabelka=LatexDataFrame.formatted(
data=dane,
index=['Value']).applymap(lambda x: f'${latex(x)}$').rename_axis(
    'Symbol',
    axis=1)

display(tabelka.reported(caption='Data given in the exercise'))

"""


class UnitRegistryIntroComponent(DocumentGenerationComponent):

    title = "Rejestr jednostek"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        display(
            ReportText(
                "To display the units of symbols used in the report in tables and charts, you must define a register of units:"
            )
        )

        display(GuideCode(f"{unit_registry_str}"))


container_code = """
nazwa_sekcji=Section('Tytuł')
CurrentContainer(nazwa_sekcji)
"""


class CurrentContainerComponent(DocumentGenerationComponent):

    title = "Klasa CurrentContainer i jej zastosowanie"

    def append_elements(self):
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]
        display(
            ReportText(
                "Wszystko co chcemy, aby znalazło się w naszym raporcie musi być wrzucone do sekcji, które potem zostaną dodane (*zaapendowane*) do naszego dokumentu. Aby utworzyć sekcję należy wywołać kod:"
            )
        )
        display(ObjectCode(container_code))
        display(
            ReportText(
                "Za pomocą **CurrentContainer** wszystko zostanie zebrane do naszej sekcji do momentu utworzenia nowej. (Zaleca się podział 1 sekcja - 1 cela w jupiterze)."
            )
        )

        from dynpy.models.mechanics.engine import Engine


Predefined_code = """

rozdzial = KineticEnergyComponent(Engine())
rozdzial_kod = KineticEnergySymPyCodeComponent(Engine())
    
doc_final.append(rozdzial)
doc_final.append(rozdzial_kod)

"""


class PredefinedSectionComponent(DocumentGenerationComponent):

    # title="Dodawanie predefiniowanej sekcji"
    title = "Adding a predefined section"

    def append_elements(self):
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # display(ReportText('Istnieje również możliwość wprowadzenia zpredefiniowanych sekcji. Możemy posłużyć się istniejącym komponentem dynamicznym na przykładzie silnika. Wykorzystamy zarówno opis teoretyczny jaki i reprezentację kodu służącą do wygenerowania równania. Aby użyć ich w dokumencie musimy wywołać te sekcje:'))
        display(
            ReportText(
                "It is also possible to introduce predefined sections. We can use an existing dynamic component in the engine example. We will use both the theoretical description and the code representation used to generate the equation. To use them in the document we need to call these sections:"
            )
        )
        from dynpy.models.mechanics import Engine

        display(ObjectCode(Predefined_code))
        display(KineticEnergyComponent(Engine()))
        display(KineticEnergySymPyCodeComponent(Engine()))


steady_state_code = """
steady_state = SDOFWinchSystem()._ode_system.solution[0]
steady_state
"""

list_code = """
Model_number = [1,2,3,4] 
inertia_list= [1,2,3,4] 
stiffness_list = [1,2,3,4] 
damping_list = [1,2,3,4] 
Excitation_list = [1,2,3,4]

system_parameters = {'Inertia':inertia_list, 'Stiffness':stiffness_list,'Damping':damping_list,'Excitation':Excitation_list}
"""

df_code = """
df = pd.DataFrame(index = Model_number, data = system_parameters)
df
"""

dict_code = """
data_dict_1 = df.T[1].to_dict()
data_dict_2 = df.T[2].to_dict()
data_dict_3 = df.T[3].to_dict()
data_dict_4 = df.T[4].to_dict()
eq1 = steady_state.subs(data_dict_1).subs(SDOFWinchSystem().F,100)
eq2 = steady_state.subs(data_dict_2).subs(SDOFWinchSystem().F,100)
eq3 = steady_state.subs(data_dict_3).subs(SDOFWinchSystem().F,100)
eq4 = steady_state.subs(data_dict_4).subs(SDOFWinchSystem().F,100)
"""

lambdify_code_str = """
function = lambdify(SDOFWinchSystem.ivar, steady_state)
t_span = np.linspace(0,100,200)
"""

sym_gen_sol = """
table = SDOFWinchSystem().subs(data_dict_2).subs(SDOFWinchSystem.F,100).numerized().compute_solution(t_span,[0,0]).iloc[:,0] 
table.plot() 
"""

list_elem = """
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
"""


solution_code = """
eq2_fun = lambdify(SDOFWinchSystem().ivar,eq2.simplify())
eq2_fun(t_span)
pd.DataFrame(index = t_span, data = eq2_fun(t_span).plot()
"""


class DifferentSimulationsComponent(ReportComponent):

    title = "Simulations - other variants"

    def append_elements(self):

        system = (
            self.reported_object
        )  # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__

        SDOFWinchSystem = type(system)

        display(ReportText("Generate a class particular solution:"))
        display(GuideCode(steady_state_code.replace("SDOFWinchSystem", system_name)))
        steady_state = (
            SDOFWinchSystem.from_random_data()
            ._ode_system.solution.with_ics([0, 0])[0]
            .subs(SDOFWinchSystem.m0, 10)
        )
        display(steady_state)

        # display(Picture('./Image/steady_sol.jpg', caption = ""))

        display(ReportText("Creating a list and a table"))
        display(GuideCode(list_code))
        Model_number = [1, 2, 3, 4]
        inertia_list = [1, 2, 3, 4]
        stiffness_list = [1, 2, 3, 4]
        damping_list = [1, 2, 3, 4]
        Excitation_list = [1, 2, 3, 4]

        system_parameters = {
            "Inertia": inertia_list,
            "Stiffness": stiffness_list,
            "Damping": damping_list,
            "Excitation": Excitation_list,
        }

        display(GuideCode(df_code))

        df1 = pd.DataFrame(index=Model_number, data=system_parameters)
        display(df1)

        # display(Picture('./Image/listatab.jpg', caption = ""))

        display(
            ReportText("Create a vector - a set of time values for our simulation: ")
        )
        display(GuideCode(lambdify_code_str.replace("SDOFWinchSystem", system_name)))

        function = lambdify(SDOFWinchSystem.ivar, steady_state)
        t_span = np.linspace(0, 100, 200)

        # display(Picture('./Image/wektor1.jpg', caption = ""))

        display(
            ReportText(
                "Creating data dictionaries and substituting the data in it into our equation so that it only depends on time: "
            )
        )

        display(GuideCode(dict_code.replace("SDOFWinchSystem", system_name)))
        # display(Picture('./Image/slowniki_sym.jpg', caption = ""))

        display(
            ReportText(
                "Create a list of elements by substituting time values one by one: "
            )
        )
        display(GuideCode(list_elem.replace("SDOFWinchSystem", system_name)))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/eq2_fun.jpg", caption=""
            )
        )

        display(ReportText("Generate a graph of a particular equation: "))
        display(GuideCode(solution_code.replace("SDOFWinchSystem", system_name)))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/wykres_niefun.jpg",
                caption="",
            )
        )

        display(ReportText("Generate a graph of the general equation: "))
        display(GuideCode(sym_gen_sol.replace("SDOFWinchSystem", system_name)))
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/ogolne_tab.jpg", caption=""
            )
        )
        display(
            Picture(
                "./dynpy/utilities/components/guides/images/ogolne_plot.jpg", caption=""
            )
        )


class ReportingComponentsList(ReportComponent):
    title = "Lista komponentów raportujących"

    def append_elements(self):

        system = self.reported_object
        import dynpy

        from .....utilities.creators import ModuleStructure

        display(
            ObjectCode(ModuleStructure(dynpy.utilities.components).get_module_tree())
        )


issue_title_str = """
Maintenance of `{system_name}` class which is dynamic system representation        
"""

issue_desc_str = """

The following problems have to be checked or fixed, in order to ensure the correctness of implemented system:

- [ ] checking if a class is executable,

- [ ] validation of figures for schemes and real examples,

- [ ] validation of reference parameters,

- [ ] validation of random parameters,

- [ ] validation of  description of  description of parameters,

- [ ] validation of units.
"""

rep_comp_call = """
    from dynpy.utilities.components.guides.en import GithubIssueReportComponent
    GithubIssueReportComponent(None)
"""

comp_output_str = """
    Details of GitHub issue
    Issue title: TEST Issue number: 999



    Issue description: Przykladowy opis - nie został podany żaden argument


"""
comp_var_str = """
    def append_elements(self):
        
        issue = self.reported_object
"""

rep_obj_str = """
    class GithubIssueReportComponent(ReportComponent):

    title="Details of GitHub issue"

    def append_elements(self):
        
        issue = self.reported_object #Tutaj przypisujemy do zmiennej argument
        if issue==None:
            from github.Issue import Issue
            issue = Issue(requester = 'lsikor', headers = {'title' : 'TEST'}, attributes = {'body':'Przykladowy opis - nie został podany żaden argument', 'title':'TEST', 'number': 999}, completed = False)
        display(Markdown((f'Issue title: {issue.title} Issue number: {issue.number}')))
        display(ReportText('\\newline'))
        if issue.body is None:
            display(Markdown("No issue description"))
        else:
            display(Markdown("Issue description: " + issue.body))
        display(ReportText('\\newline'))
        display(ReportText('-'*130))
        display(ReportText('\\newline'))
"""

rep_jup_str = """
    from dynpy.utilities.creators import GitHubInterface
    from dynpy.utilities.components.guides.en import GithubIssueReportComponent
    from datetime import datetime
    user = '*' # Podaj nazwę uzytkownika github, podanie * oznacza wylistowanie wszystkich
    since_date = "2024-09-05" # Zmień datę, od której ma listować issues
    # token = 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 # Token autoryzujący dostep do repo, każdy powinien wygenerować własny, można googlowac
    client = GitHubInterface()

    date_string = since_date + "T00:00:00Z" # Dodanie do podanej wyżej daty odpowiedniego formatu, który przyjmuje atrybut Pygithub
    date_object = datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%SZ") # Zamiana daty na dateobject z odpowiednim formatowaniem

    list_open=[] #pusta lista
    list_open = client.get_issues_list(repo_name='bogumilchilinski/dynpy', state='open', assignee=user, sort='created', since=date_object) # pozyskanie listy issues z repo

    for single_issue in list_open: # pętla na wylistowanie issues 
        GithubIssueReportComponent(single_issue) # wykorzystanie zaimplementowanego komponentu, argument jaki przekazujemy to pojedyńcze issue z listy
"""


class ReportingCompsUsageComponent(ReportComponent):

    title = "Działanie komponentów raportujących"

    @property
    def reported_object(self):

        default_data = {
            "classname": "ReportingCompsUsageComponent",
            "module": "guide.en.py",
            "field": "guide or report",
            "target": "`ODESystem` class",
            "issue_no": 359,
        }

        if isinstance(self._reported_object, dict):
            return {**default_data, **self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data, "classname": self._reported_object}

        elif self._reported_object is None or not isinstance(
            self._reported_object, dict
        ):
            return default_data

        else:
            return self._reported_object

    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object = obj

    def append_elements(self):

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        display(
            ReportText(
                """Po stworzeniu komponentu możesz go wywołać poprzez zaimportowanie odpowiedniej klasy we własnym 
                           jupku,poniżej zostanie opisany konkretny przykład, aby lepiej zrozumieć działanie komponentu"""
            )
        )
        display(ObjectCode(rep_comp_call))
        display(
            ReportText(
                "Aby dokładniej zrozumieć metodę `append_elements` należy przyjrzeć się temu co mamy w (), domyślnie każdy komponent sprawdzamy wpisując tam `None`"
            )
        )
        display(
            ReportText(
                """`None` - jest to argument jaki przekazujemy klasie, w tym przypadku oznacza to, że nie przekazujemy żadnego konkretnego argumentu, ponieważ None = Nic
                           Jeżeli komponent został zabezpieczny to podanie tego argumentu pokaże tylko przykład zaimplementowany w metodzie
                           Uruchamiając poprzedni przykład w jupku otrzymamy poniższy output:"""
            )
        )
        display(ObjectCode(comp_output_str))
        display(
            ReportText(
                "W metodzie `append_elements` zostało przypisane, aby komponent przekazał argument do zmiennej."
            )
        )
        display(ObjectCode(comp_var_str))
        display(
            ReportText(
                "Przypisany w () argument zostaje przekazany do klasy ReportComponent, dokładniej do metody _init_ a następnie przypisany do `self.reported_object`"
            )
        )
        display(
            ReportText(
                "Wygląda to następująco: JUPYTER > WYWOŁANIE KOMPONENTU > _INIT_(ReportComponent) > SELF.REPORTED_OBJECT"
            )
        )
        rep_pic1 = Picture(
            "./dynpy/utilities/components/guides/images/OpOfComp.png",
            caption="Droga podanego argumentu",
            width="9cm",
        )
        display(rep_pic1)
        display(ReportText("Poniżej pełny kod tego przykładowego komponentu:"))
        display(ObjectCode(rep_obj_str))
        display(
            ReportText("A tutaj przykładowy jupek z wykorzystaniem tego komponentu:")
        )
        display(ObjectCode(rep_jup_str))
