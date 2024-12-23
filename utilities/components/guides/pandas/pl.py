from  ...mechanics import *
from  ...mechanics import ReportComponent as BaseReportComponent

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from .....solvers.linear import *
from .....dynamics import *

from ....report import display as ip_display


from  ...mechanics import display



miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

długosc_dnia_w_miesiacach_godziny = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia temperatura':srednie_temp_list}

df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)


        
#pandas guide
data_code=(
'''
miesiace_list = ['styczeń','luty','marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']  
srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]  
Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]  
Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]
''')

dict_code=(
'''
data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii         [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia             temperatura':srednie_temp_list}
''')
output_code=(
'''
df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)
df
''')

class PandasTableGenerationComponent(ReportComponent):
    
    title="Tworzenie tabelki"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Stworzenie listy danych:'))

        display(ObjectCode(data_code))


        display(ReportText('Stworzenie słownika:'))

        display(ObjectCode(dict_code))


        display(ReportText('Wywolanie tabelki:'))
        display(ObjectCode(output_code))
        display(df)
        

        
class PandasMethodsComponent(ReportComponent):
    
    title="Metody biblioteki Pandas"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        
        miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

        srednie_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

        Eg_dzienne_list_watogodziny_na_metr2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

        Eg_dzienne_kilowatogodziny_na_metr2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

        długosc_dnia_w_miesiacach_godziny = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

        data_warunki_atmosferyczne = {'Długość dnia w miesiącu':długosc_dnia_w_miesiacach_godziny,'Dzienne natezenie energii [kWh/m^2]':Eg_dzienne_list_watogodziny_na_metr2,'Średnia temperatura':srednie_temp_list}

        df = pd.DataFrame(index = miesiace_list,data = data_warunki_atmosferyczne)

        display(ReportText('***Metoda iloc:***'))
        display(ObjectCode('df.iloc[0:3]'))
        display(df.iloc[0:3])

        # display(Picture('./Image/iloc.jpg', caption = ""))

        display(ReportText('***Metoda loc:***'))
        display(ObjectCode('df.loc["styczeń"]'))
        display(df.loc["styczeń"])
        # display(Picture('./Image/loc.jpg', caption = ""))

        display(ReportText('***Metoda set axis dziala za rowno dla kolumn jak i rzędów. Dla rzędów:***'))
        display(ObjectCode("def.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'],axis = 'index')"))
        display(df.set_axis(['a','b','c','d','e','f','g','h','i','j','k','l'], axis = 'index'))

        # display(Picture('./Image/set_axis.jpg', caption = ""))

        display(ReportText('***Dla kolumn:***'))
        display(ObjectCode("df.set_axis(['a','b','c'],axis = 'columns')"))
        display(df.set_axis(['a','b','c'],axis = 'columns'))
        # display(Picture('./Image/set_axis2.jpg', caption = ""))

        display(ReportText('***Metoda rename - Zmienienie pojedynczej kolumny:***'))
        display(ObjectCode('''df.rename(columns = {"Długość dnia w miesiącu":'AAA'}, index = {"styczeń":'A'} )'''))
        display(df.rename(columns = {"Długość dnia w miesiącu":'AAA'}, index = {"styczeń":'A'} ))
        # display(Picture('./Image/rename.jpg', caption = ""))

        display(ReportText('***Metoda applymap:***'))
        display(ObjectCode('     df.applymap(lambda x:x+2)'))
        display(df.applymap(lambda x:x+2))

        # display(ReportText('***Metoda applymap dla kolumny/wiersza:***'))
        # display(Markdown('''
        # Warto pamiętać, że df[NAZWA] => tworzy serię, a df[[Nazwa]] => tworzy nowy frame, który można wykorzystać do użycia .applymap'''))
        # display(df[['Długość dnia w miesiącu']].applymap(lambda x:x+100))
        # pd_h = df[['Długość dnia w miesiącu']].applymap(lambda x:x+100)
        # list1 = pd_h['Długość dnia w miesiącu'].tolist()
        # df['Długość dnia w miesiącu'] = list1

        # display(df)

        # display(Picture('./Image/applymap.jpg', caption = ""))

        display(ReportText('***Slicowanie:***'))
        display(ObjectCode('''df['Długość dnia w miesiącu']'''))
        display(df['Długość dnia w miesiącu'])
        # display(Picture('./Image/slice2.jpg', caption = ""))
        

eoms_code=(
'''
eoms_eq=Eq(eoms,0)
eoms_eq''')
dict_code=(
'''
param_dict = dyn_sys.get_numerical_parameters()
t_span = np.linspace(0,100,1001)
ic_list = [0.0,0.0] ### Dla układu o jednym stopniu swobody
#param_dict = {**param_dict,dyn_sys.m:dyn_sys.m}
''')
params_code=(
'''
parameter = dyn_sys.system_parameters()[0]
param_dict = {**param_dict,parameter:parameter}
''')

num_sys_df_code=(
'''
num_sys=dyn_sys.subs(param_dict)('Num')
'''
)

na_df_code=(
'''
na_df = num_sys.numerical_analysis(parameter=parameter,param_span=[1,2,3], t_span = t_span)
na_df
''')
na_df_code_eoms=(
'''
na_df = num_sys.eoms.numerical_analysis(parameter=parameter,param_span=[1,2,3], t_span = t_span)
na_df
''')

na_df_sol_code=(
'''
na_df_sol = na_df.with_ics([2.0,0.0]).compute_solution(t_span)
na_df_sol.plot()
''')

full_example_na_df_sol_code=(
f'''
{eoms_code}

{dict_code}

{params_code}

{num_sys_df_code}

{na_df_code}

{na_df_code_eoms}

{na_df_sol_code}

''')

class NumericalAnalysisSimulationComponent(ReportComponent):

    title="Wykonanie symulacji z wykorzystaniem metody numerical_analysis"


    def append_elements(self):

        from .....solvers.linear import ODESystem
        
        dyn_sys = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        eoms=dyn_sys._eoms[0]
        eoms

        display(ReportText('Proces wywołowania równania ruchu za pomoca metody eoms:'))

        display(ReportText('Wynik jest następujący:'))
        
        display(ObjectCode(eoms_code))
        
        eoms_eq=Eq(eoms,0)
        eoms_eq
        
        display(SympyFormula(eoms_eq))
        
        display(ReportText('Zaprezentowane równanie pozwala sprawdzić, czy rozważany układ jest poprawnie zaimportowany do przestrzeni roboczej.'))
        
        display(ReportText('Kolejnym kluczowym elementem jest zdefiniowanie parametrów układu, wektora czasu i warunków początkowych.'))
        
        param_dict = dyn_sys.get_numerical_parameters()
        t_span = np.linspace(0,100,1001)
        ic_list = [1.0,0.0] ### Dla układu o jednym stopniu swobody
        param = dyn_sys.system_parameters()[0]
        
        
        display(ObjectCode(dict_code))
        
        display(ReportText('Następnie należy określić, który parametr ma zostać poddany analizie. Parametr ten należy wówczas zmienić w słowniku parametrów z wartości liczbowej na symbol.'))
        
        parameter = dyn_sys.system_parameters()[0]
        param_dict = {**param_dict,parameter:parameter}
        
        display(ObjectCode(params_code))
        display(ObjectCode(num_sys_df_code))
        
        
        display(ReportText('Ostatnim etapem jest wywołanie metody numerical_analysis dla rozważanego systemu. Możemy to zrobić przy użyciu property eoms:'))
        display(ObjectCode(na_df_code_eoms))
        display(ReportText('Lub bez niej:'))
        display(ObjectCode(na_df_code))
        num_sys=dyn_sys.subs(param_dict)('Num')
        na_df = dyn_sys.subs(param_dict).numerical_analysis(parameter=dyn_sys.system_parameters()[0],param_span=[1,2,3], t_span = t_span).with_ics([2.0,0.0])
        #na_df = system.subs(param_dict).numerized()
        na_df
        
        display(ReportText('Symulacje numeryczną wykonuje się w analogiczny sposób do pozostałych symulacji w następujący sposób:'))

        display(ObjectCode(na_df_sol_code))
        
        na_df_sol = na_df.compute_solution(t_span = t_span, ic_list = [2.0,0.0])
        (na_df_sol.plot())
        
        import matplotlib.pyplot as plt
        plt.show()
        
        display(ReportText('Pełny kod potrzebny do uzyskania symulacji wygląda następująco:'))
        display(ObjectCode(full_example_na_df_sol_code))


sim_code=(
'''
ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)
ode_simulation.plot()
''')

class AnalyticalSimulationComponent(ReportComponent):

    title="Wykonanie symulacji z wykorzystaniem rozwiązania analitycznego"

    
    
    def append_elements(self):

        from .....solvers.linear import ODESystem
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        eoms=system._eoms[0]
        eoms

        display(ReportText('Proces wywołowania równania ruchu za pomoca metody eoms:'))

        display(ReportText('Wynik jest następujący:'))
        
        display(ObjectCode(eoms_code))
        
        eoms_eq=Eq(eoms,0)
        eoms_eq
        
        display(SympyFormula(eoms_eq))
        
        ode_system = ODESystem.from_dynamic_system(system.linearized())

        display(ReportText('Podane równanie w tym przypadku nie będzie odpowiednim wyborem. Do wykonania symulacji analitycznej należy zastosować system dynamiczny w postaci równań różniczkowych zwyczajnych, który można uzyskać za pomocą klasy ODESystem:'))
        
        display(ReportText('Wynik jest następujący:'))

        display(SympyFormula(ode_system))
        
        param_dict = system.get_numerical_parameters()
        t_span = np.linspace(1,100,101)
        ic_list = [1.0,0.0] ### Dla układu o jednym stopniu swobody
        
        
        #param_dict = {**param_dict,system.m:system.m}
        
        
        display(ObjectCode(dict_code))
        
        
        display(ReportText('Aby przeprowadzić symulację analityczną w pierwszej kolejności konieczne jest uzyskanie rozwiązania równań ruchu dla podstawionych danych subs(param_dict), co można wykonać za pomocą metod `solution` i `steady_solution` w zależności od skomplikowania systemu.'))
        
        ode_solution = ode_system.subs(param_dict).steady_solution
        
        display(ReportText('Wynik jest następujący:'))
        
        display(SympyFormula(ode_solution))
        
        display(ReportText('Do przeprowadzenia symulacji konieczne jest zdefiniowanie parametrów układu, wektora czasu oraz warunków początkowych. Można je uzyskać w następujący sposób:'))
        

        
        display(ObjectCode(sim_code))
        
        display(ReportText('Dla zdefiniowanych parametrów symulacja analityczna wygląda w następujący sposób:'))
        
        ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)

        (ode_simulation.plot())
        import matplotlib.pyplot as plt
        plt.show()

doc_gen_str =(
'''

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

''')

class DocumentGenerationComponent(ReportComponent):

    #title="Generowanie dokumentu"
    title="Document generation"


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

        elif self._reported_object is None or not isinstance(self._reported_object, dict):
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
        
        #display(ReportText('Ostatni krok to zaapendowanie sekcji i utworzenie dokumentu pdf :'))
        display(ReportText('The last step is to append a section and create a pdf document :'))
        display(GuideCode(f'{doc_gen_str}'))


        
TablesCreationComponent_str=(
'''
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
    
''')

class TablesCreationComponent(DocumentGenerationComponent):

    title="Tworzenie tabelek"


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
        
       
        #implement reporting activieties here
        
        display(ReportText('Aby utworzyć tabelkę najpierw musimy posiadać dane z obliczeń na interesującym nas systemie dynamicznym:'))
        display(GuideCode(f'{TablesCreationComponent_str}'))



class AlgebraicExpressionComponent(ReportComponent):

    title = "Algebraic expression example"

    def append_elements(self):
        display(ReportText('Example Algebraic expression in Sympy library:'))
        display(GuideCode(f'{SympyFormulaComponent_str}'))
        
        from sympy import Symbol, Function

        M = Symbol('M')

        t = Symbol('t')

        m_e = Symbol('m_e')

        e = Symbol('e')

        z = Function('z')

        varphi = Function('varphi')

        Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2
        
        display(ReportText('After executing above code you can display kinetic energy algebraic expression:'))
        
        display(Eq(Symbol('T'),Ek))