import numpy as np
import pandas as pd
from pandas import *
from sympy import *
from sympy import lambdify
from sympy.physics.mechanics import dynamicsymbols

from .....dynamics import *
from .....solvers.linear import *
from ...mechanics import ReportComponent as BaseReportComponent
from ...mechanics import *
from . import pl

months_list = [
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
]

import datetime

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

average_temp_list = [-1.9, -0.8, 3.2, 9.3, 14.6, 18, 20.1, 19.5, 14.7, 9.3, 4.8, 0.5]

Eg_daily_list_Wh_per_meter2 = [
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

Eg_daily_list_kWh_per_meter2 = [0.6, 1, 3, 3.8, 4.8, 5.3, 4.9, 3.3, 1.7, 0.7, 0.5]

day_length_in_months_hours = [
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

data_atmospheric_conditions = {
    "Day length in a month [h]": day_length_in_months_hours,
    "Daily energy intensity [${kWh/m^2}$]": Eg_daily_list_Wh_per_meter2,
    "Average temperature [$^{\circ}$C]": average_temp_list,
}

df = pd.DataFrame(index=months_list, data=data_atmospheric_conditions)


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
        self._reported_object = obj

    @property
    def _system(self):
        print(
            "Kod do poprawienia #################### bo stosujesz starą zmienną self._system"
        )
        print("zamień se self._system na self.reported_object")

        return self.reported_object


# pandas guide
data_code = """
months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']
average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]
day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]
Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]
Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]
"""

data_dict_code = """
data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}
"""
output_code = """
import pandas as pd

df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)
df
"""


class PandasTableGenerationComponent(ReportComponent):

    title = "Generating a table"

    def append_elements(self):

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText("Creating a list of data:"))

        display(GuideCode(data_code))

        display(ReportText("Creating dictionary:"))

        display(GuideCode(data_dict_code))

        display(ReportText("Calling the table:"))
        display(GuideCode(output_code))
        display(df)


class PandasMethodsComponent(ReportComponent):

    title = "Pandas library methods"

    def append_elements(self):

        # system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        months_list = [
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]

        average_temp_list = [
            -1.9,
            -0.8,
            3.2,
            9.3,
            14.6,
            18,
            20.1,
            19.5,
            14.7,
            9.3,
            4.8,
            0.5,
        ]

        Eg_daily_list_Wh_per_meter2 = [
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

        Eg_daily_list_kWh_per_meter2 = [
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

        day_length_in_months_hours = [
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

        data_atmospheric_conditions = {
            "Day length in a month [h]": day_length_in_months_hours,
            "Daily energy intensity [${kWh/m^2}$]": Eg_daily_list_Wh_per_meter2,
            "Average temperature [$^{\circ}$C]": average_temp_list,
        }

        df = pd.DataFrame(index=months_list, data=data_atmospheric_conditions)

        display(ReportText("Iloc method:"))
        display(GuideCode("df.iloc[0:3]"))
        display(df.iloc[0:3])

        # display(Picture('./Image/iloc.jpg', caption = ""))

        display(ReportText("Loc method:"))
        display(GuideCode('df.loc["January"]'))
        display(df.loc["January"])
        # display(Picture('./Image/loc.jpg', caption = ""))

        display(
            ReportText("The set axis method works for both columns and rows. For rows:")
        )
        display(
            GuideCode(
                "df.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'],axis = 'index')"
            )
        )
        display(
            df.set_axis(
                ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
                axis="index",
            )
        )

        # display(Picture('./Image/set_axis.jpg', caption = ""))

        display(ReportText("For columns:"))
        display(GuideCode("df.set_axis(['a','b','c'],axis = 'columns')"))
        display(df.set_axis(["a", "b", "c"], axis="columns"))
        # display(Picture('./Image/set_axis2.jpg', caption = ""))

        display(ReportText("Metoda rename - Zmienienie pojedynczej kolumny:"))
        display(
            GuideCode(
                """df.rename(columns = {"Day length in a month [h]":'AAA'}, index = {"January":'A'} )"""
            )
        )
        display(
            df.rename(
                columns={"Day length in a month [h]": "AAA"}, index={"January": "A"}
            )
        )
        # display(Picture('./Image/rename.jpg', caption = ""))

        display(ReportText("Map method:"))
        display(GuideCode("df.map(lambda x:x+2)"))
        #         display(df.map(lambda x:x+2)) - nowa wersja pandasa zmieniła nazwę "map" na "applymap" - Michał SZ
        display(df.map(lambda x: x + 2))

        # display(ReportText('***Applymap mthod for colummn/row:***'))
        # display(Markdown('''
        # Remeber, that df[NAME] => creates series, a df[[NAME]] => creates new frame, which can be used with .applymap'''))
        # display(df[['Day length in a month']].applymap(lambda x:x+100))
        # pd_h = df[['Day length in a month']].applymap(lambda x:x+100)
        # list1 = pd_h['Day length in a month'].tolist()
        # df['Day length in a month'] = list1

        # display(df)

        # display(Picture('./Image/applymap.jpg', caption = ""))

        display(ReportText("Slicing:"))
        display(GuideCode("""df['Day length in a month [h]']"""))
        display(df["Day length in a month [h]"])
        # display(Picture('./Image/slice2.jpg', caption = ""))


class NumericalAnalysisSimulationComponent(pl.NumericalAnalysisSimulationComponent):

    title = "Numerical analysis simulation Component"


class AnalyticalSimulationComponent(pl.AnalyticalSimulationComponent):

    title = "Analytical simulation Component"


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


TablesCreationComponent_str = '''
from dynpy.models.mechanics.tmac import SDOFWinchSystem
eoms= SDOFWinchSystem()._eoms[0]
eoms_eq_raw=Eq(eoms,0)
data=SDOFWinchSystem().get_random_parameters()
Ms0=Symbol('M_s0')
b=Symbol('b')
omg=Symbol('varphi')
ph1=Symbol('Phi')

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

general_sol_matrix=nowe_ode.solution.with_ics([0,0])
table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)
table=table_eq[[phi.diff(t)]]
table_for_report = table[0::15].to_latex_dataframe()

display(Markdown(' Do tablicy `table_eq` zapisujemy równanie ogólne dla zmiennej czasowej `t_span`   '))
display(Markdown(
f"""
    
    table_eq=general_sol_matrix.n().numerized().compute_solution(t_span)

"""))
display(Markdown(f' Do tablicy table przypisujemy wartości ${latex(omg)}$ powstałe w wyniku zróżniczkowania ${latex(ph1)}$ w tablicy `table_eq` '))
display(Markdown(
f"""
    table=table_eq[[phi.diff(t)]]
   
"""))
display(ReportText(' W następnym kroku wybieramy zakres wartości który chcemy umieścić w tabelce przy pomocy naiwasów kwadratowych. Nastepnie przekształcamy tabelę na zgodną z LateX przy pomocy metody: to\_latex\_dataframe()  '))
display(Markdown(
f"""

    table_for_report = table[0::15].to_latex_dataframe()

"""))
display(ReportText(' Naszą tabelkę opisujemy w tekście, pożądane jest użycie AutoMarkerów które zostały opisane w poniższej sekcji tego dokument '))
display(Markdown(
f"""

     display(Markdown(f'Wyniki symulacji zostały zestawione w tabeli {AutoMarker(table_for_report)}'))

"""))
display(ReportText(' Aby wyświetlic tabelkę korzystamy z  klasy display. Metoda .reported umożliwia nam nadanie tabelce podpisu. '))
display(Markdown(
f"""

    display(table_for_report.reported(caption='Tabela'))

"""))
    
'''


class TablesCreationComponent(DocumentGenerationComponent):

    title = "Tworzenie tabelek"

    def append_elements(self):
        # variables provided by `reported_object` arg
        classname = self.reported_object["classname"]
        class_module = self.reported_object["module"]
        class_field = self.reported_object["field"]
        target = self.reported_object["target"]

        # implement reporting activieties here

        display(
            ReportText(
                "Aby utworzyć tabelkę najpierw musimy posiadać dane z obliczeń na interesującym nas systemie dynamicznym:"
            )
        )
        display(GuideCode(f"{TablesCreationComponent_str}"))


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


class AlgebraicExpressionComponent(ReportComponent):

    title = "Algebraic expression example"

    def append_elements(self):
        display(ReportText("Example Algebraic expression in Sympy library:"))
        display(GuideCode(f"{SympyFormulaComponent_str}"))

        from sympy import Function, Symbol

        M = Symbol("M")

        t = Symbol("t")

        m_e = Symbol("m_e")

        e = Symbol("e")

        z = Function("z")

        varphi = Function("varphi")

        Ek = (
            M * Derivative(z(t), t) ** 2 / 2
            + m_e
            * (e * sin(varphi(t)) * Derivative(varphi(t), t) - Derivative(z(t), t)) ** 2
            / 2
        )

        display(
            ReportText(
                "After executing above code you can display kinetic energy algebraic expression:"
            )
        )

        display(Eq(Symbol("T"), Ek))


pandas_intro_call_str = """
from dynpy.utilities.documents.guides import IntroToPandasGuide
from dynpy.utilities.components.guides.en import *
from dynpy.models.mechanics.pendulum import ForcedSpringMassSystem

IntroToPandasGuide(reported_object=ForcedSpringMassSystem());

"""


pandas_intro_str = """Guide consist of bellow components which can be called separatelly:
*IntroToPandasUsageComponent
*PandasTableGenerationComponent
*PandasMethodsComponent
*BasicOperationsComponent
*DynamicSystemCallComponent
*SimulationsComponent
*DifferentSimulationsComponent
"""


class IntroToPandasUsageComponent(ReportComponent):

    title = "Introduction to basics of pandas library"

    def append_elements(self):

        from dynpy.models import mechanics

        display(
            ReportText(
                "This guide concers basic usage of pandas library with dynamical systems simulations in `DynPy` library."
            )
        )
        display(
            ReportText(
                "Basic call of it is as follows and runs default dynamic system which is `ForcedSpringMassSystem"
            )
        )
        display(ObjectCode(pandas_intro_call_str))

        display(ReportText(pandas_intro_str))
