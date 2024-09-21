from  ...mechanics import *
from  ...mechanics import ReportComponent as BaseReportComponent

from . import pl

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from .....solvers.linear import *
from .....dynamics import *

months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']

from  ...mechanics import display
import datetime


miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}

df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)

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
        self._reported_object=obj

    @property
    def _system(self):
        print('Kod do poprawienia #################### bo stosujesz starą zmienną self._system')
        print('zamień se self._system na self.reported_object')
        
        return self.reported_object



        
# #obsluga systemow dynamicznych


imports_code_str = '''
from sympy import *
import numpy as np
from pandas import *
from sympy.physics.mechanics import init_vprinting, dynamicsymbols
from pylatex import Document, Section, Subsection, Itemize, Package, HorizontalSpace, Description, Marker, Command 
from pylatex.section import Paragraph, Chapter 
from pylatex.utils import italic, NoEscape 
from dynpy.utilities.adaptable import *
from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase 
from dynpy.utilities.templates.tikz import TikzCaSCStandalone
from dynpy.utilities.report import (Markdown, Picture, SympyFormula, SymbolsDescription, DescriptionsRegistry, ObjectCode, CurrentContainer)
from dynpy.solvers.linear import ODESystem
init_vprinting()
from dynpy.models.mechanics import SDOFWinchSystem

'''



class DynamicSystemCallComponent(ReportComponent):
    
    title="Introduction"


    def append_elements(self):

        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText('Required library/class imports to create the report:'))


        display(GuideCode(imports_code_str.replace('SDOFWinchSystem',system_name)))


        display(ReportText('A class that allows you to solve a specific problem (in this case, a problem related to the description of a crane):'))

        display(GuideCode("""from dynpy.models.mechanics import SDOFWinchSystem""".replace('SDOFWinchSystem',system_name)      ))
        
        
        display(GuideCode(f'system=dyn_sys={system_name}()'  ))

        display(ReportText('Path to the sought-after class on the CoCalc platform:'))
        #display(Picture('./dynpy/utilities/components/guides/images/sciezka_w.jpg'))

        display(ReportText((system.__class__.__module__)))
        
        display(ReportText('The way to call out the preview of the class - the so-called preview:'))

        display(GuideCode(f'''{system_name}().as_picture()'''  ))
        display(system.as_picture())
        
        
class DynamicSystemMethodsUsageComponent(ReportComponent):

    title="Wywolywanie rownan ruchu oraz innych metod"



    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        #from dynpy.models.mechanics import ForcedSpringMassSystem as SDOFWinchSystem
        
        eoms=system._eoms[0]

        display(ReportText('Proces wywołowania równania ruchu za pomoca metody eoms:'))
        display(GuideCode('''eoms=SDOFWinchSystem()._eoms[0]'''.replace('SDOFWinchSystem',system_name)))



        display(ReportText('Wynik jest następujący:'))
        #display(SympyFormula(eoms))


        eoms_eq=Eq(eoms,0)
        eoms_eq

        display(ReportText('Proces tworzenia równania z wcześniej wywołanego równania ruchu (równanie początkowo jest wywoływane jako ciąg znaków reprezentujący opis danego problemu, jednak nie jest on wówczas równaniem, które możemy zastosować do obliczeń):'))
        display(GuideCode('eoms_eq=Eq(eoms,0)'.replace('SDOFWinchSystem',system_name)))
        display(ReportText('Wynik jest następujący:'))
        display(SympyFormula(eoms))

        #display(Picture('./Image/eomseq_w.jpg'))
        
        slownik=system.get_random_parameters()
        
        slownik_numerical=system.get_numerical_parameters()
        
        display(ReportText('Słownik z losowymi wartościami parametrów symbolicznych systemu wygląda następująco:'))
        
        display(SympyFormula(Dict(slownik)))
        
        display(ReportText('Słownik z losowymi wartościami liczbowymi do analizy numerycznej systemu, wygląda następująco:'))
        
        display(SympyFormula(Dict(slownik_numerical)))

        display(ReportText('Tworzenie słownika (w tym przypadku z losowymi wartościami) umożliwiający nadanie wartości zmiennym:'))
        
        

        display(GuideCode('slownik=SDOFWinchSystem().get_random_parameters()'.replace('SDOFWinchSystem',system_name)))
        display(GuideCode('slownik_numerical=SDOFWinchSystem().get_numerical_parameters()'.replace('SDOFWinchSystem',system_name)))
        #display(Picture('./Image/slownik_w.jpg'))

        display(ReportText('Rozwiązanie ogólne równania:'))
        #display(Picture('./Image/ogolne_w.jpg'))

        solution=system._ode_system.solution[0]
        display(GuideCode('''solution=SDOFWinchSystem()._ode_system.solution[0]'''.replace('SDOFWinchSystem',system_name)))
        display(GuideCode('solution_sym=SDOFWinchSystem()._ode_system.solution'.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(solution))

        display(ReportText('Rozwiązanie szczególne równania:'))



        steady_solution=system._ode_system.steady_solution[0]
        display(GuideCode('''steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(steady_solution))
    
    


class SimulationsComponent(ReportComponent):
    
    title="Creating graphs"

    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__


        display(ReportText('Definition of time vector:'))
        #display(Picture('./Image/tspan_w.jpg'))

        display(GuideCode('''t_span = np.linspace(0,100,200)'''.replace('SDOFWinchSystem',system_name)))
        t_span = np.linspace(0,100,200)


        display(ReportText('The way to create charts manually:'))
        display(ReportText('Substitution of data from a previously created dictionary:'))
        #display(Picture('./Image/stedisub_w.jpg'))

        
        
        display(GuideCode(steady_sol_str.replace('SDOFWinchSystem',system_name)))
        #dict=system.get_random_parameters()
        dict_numerical=system.get_numerical_parameters()
        
        steady_solution=system._ode_system.steady_solution.rhs[0]
        steady_solution_subs=steady_solution.subs(dict_numerical)
        display(steady_solution_subs)
        display(ReportText('Using the lambdify method to convert functions to anonymous form'))
        #display(Picture('./Image/lambdify_w.jpg'))


#         eq_lambdify=lambdify(system.ivar, list(steady_solution_subs.rhs))
        eq_lambdify=lambdify(system.ivar, steady_solution_subs)
        display(GuideCode('''eq_lambdify=lambdify(system.ivar, list(steady_solution_subs.rhs))'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(eq_lambdify))


        display(ReportText('Creating a graph:'))
        #display(Picture('./Image/plot2_w.jpg'))
        
        #pd.DataFrame(data=eq_lambdify(t_span),index=t_span).plot()
        #plt.show
        df = pd.DataFrame(data=eq_lambdify(t_span))
        display(df.plot())
        display(GuideCode('''df = pd.DataFrame(data=eq_lambdify(t_span))
df.T.set_index(t_span).plot()'''.replace('SDOFWinchSystem',system_name)))


        display(ReportText('A way to create charts directly from a previously called class:'))
        display(ReportText('Dictionary substitution:'))
        #display(Picture('./Image/table_in.jpg'))

        
        #zmiana slownik na slownik_numerical
        table_eq=system._ode_system.steady_solution.subs(dict_numerical)


        display(GuideCode('''dict_numerical=system.get_numerical_parameters()
table_eq=system._ode_system.steady_solution.subs(dict_numerical)
table_eq'''.replace('SDOFWinchSystem',system_name)))
        display(SympyFormula(table_eq))

        #jak tabele zrobic + ew jej output
        display(ReportText('The way to create a table:'))
        #display(Picture('./Image/table_in.jpg'))
        coord = system.q[0]
        table=table_eq.numerized().compute_solution(t_span)
        table_new=table[coord]

        display(GuideCode(solution_table_str.replace('SDOFWinchSystem',system_name)))
        display(table)
        display(table_new)

        display(ReportText('The way to create a visualization of the defined table:'))
        #display(Picture('./Image/table_in.jpg'))

##### nazwy image tez przetlumaczone
        display(ReportText('How to create the table displayed in the finished report:'))
        #display(Picture('./Image/creating_table_in.png'))

        display(GuideCode(report_table_define_str.replace('SDOFWinchSystem',system_name)))

        display(ReportText('Visualization of the graph:'))
        #display(Picture('./Image/tikzplot_w.jpg'))

        display(GuideCode('''table_new.to_pylatex_tikz().in_figure()'''.replace('SDOFWinchSystem',system_name)))    

section_development_str=(
'''
from dynpy.utilities.report import ReportText
intro=Section('Introduction')
CurrentContainer(intro)

#scheme=Picture('picture_winch.png')
#display(scheme)

picture = Picture(system._as_picture().__dict__['image'], caption = f'System diagram')
display(picture)
''')

forming_equation_str=(
'''
eom=Section('Equations of motion')
CurrentContainer(eom)
display(ReportText('Equations of motion: '))
eoms=system._eoms[0]

slownik=system.get_random_parameters()
slownik_numerical=system.get_numerical_parameters()

eoms_eq=Eq(eoms,0)
solution_subs = system._ode_system.steady_solution.subs(slownik_numerical)
steady_solution=system._ode_system.steady_solution[0]
steady_solution_subs=steady_solution.subs(slownik_numerical)
#zmiana slownik na slownik numerical

steady_solution_subs
display(solution_subs)

ode_system = system._ode_system #ESystem.from_dynamic_system(system.linearized())
param_dict = system.get_numerical_parameters()
t_span = np.linspace(1,10,1001)
ic_list = [0.0,0.0] ### Dla układu o jednym stopniu swobody
ode_solution = ode_system.subs(param_dict).steady_solution

display(ode_solution)

display(SympyFormula(eoms_eq))
display(SympyFormula(solution_subs))
display(SympyFormula(steady_solution_subs))
''')

figure_eoms=(
'''
wykres=Section('Graph of dynamic system ')
CurrentContainer(wykres)
display(ReportText('The graph shows the solution of the equations of motion for the dynamic system under consideration. The graph shows the dependence of the displacement (or other relevant quantity, depending on the type of system) as a function of time, obtained from the solution of the differential equations.'))

ode_simulation = ode_solution.subs(param_dict).compute_solution(t_span, ic_list = ic_list)
picture=ode_simulation.to_pylatex_tikz(subplots=True).in_figure(caption='Graph of the solution of the equations of motion for the dynamic system under consideration.')

display(picture)

display('Equations of motion:')

display(eoms_eq)
display(solution_subs)
display(steady_solution_subs)
''')

report_generation_libraries=(
'''
doc_final = MechanicalCase('./output/Document',documentclass=NoEscape('article'),document_options=['a4paper','fleqn'],lmodern=False)

doc_final.append(intro)
doc_final.append(eom)
doc_final.append(wykres)

doc_final.generate_pdf()
''')




class DynSysCodeComponent(ReportComponent):

    title="Calling dynamic system code"



    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        system_name = system.__class__.__name__
        
        display(ReportText('''The code of the analyzed class is as follows'''))
        display(ObjectCode(system.__class__))
        display(ReportText(f'''The presented code of the {system_name} class was retrieved with the help of the ObjectCode class.  '''))
        
        

        
class DynSysIntroComponent(ReportComponent):

    title="Introduction to using the DynSys module"
    
    def append_elements(self): 
        from dynpy.models.mechanics.tmac import SDOFWinchSystem
        
        display(ReportText('Required imports of libraries/individual classes to create the report:'))
        display(ObjectCode(
        '''
                from sympy import *
                import numpy as np
                from pandas import *
                from sympy.physics.mechanics import init_vprinting, dynamicsymbols
                from pylatex import Document, Section, Subsection, Itemize, Package, HorizontalSpace, Description, Marker, Command 
                from pylatex.section import Paragraph, Chapter 
                from pylatex.utils import italic, NoEscape 
                from dynpy.utilities.adaptable import *
                from dynpy.utilities.templates.document import BeamerTemplate, MechanicalCase 
                from dynpy.utilities.templates.tikz import TikzCaSCStandalone
                from dynpy.utilities.report import (SystemDynamicsAnalyzer, DataPlot, AccelerationComparison, FFTComparison, ReportEntry, SimulationalBlock,               ReportText, SimulationFFT, DataStorage, Markdown, SummaryTable, Picture, SympyFormula, SymbolsDescription, DescriptionsRegistry,
                ObjectCode, CurrentContainer)
                from dynpy.solvers.linear import ODESystem
                init_vprinting()
                from dynpy.models.mechanics.tmac import SDOFWinchSystem

        '''))
        display(ReportText('A class that allows you to solve a specific problem (in this case, a problem related to the description of a crane):'))
        display(Markdown('''
            from dynpy.models.mechanics.tmac import SDOFWinchSystem
            '''))

        display(ReportText('The way to call the class preview:'))

        display(ObjectCode('''
            SDOFWinchSystem().preview()
            '''))

  
        eoms=SDOFWinchSystem()._eoms[0]
        eoms

        display(ReportText('The process of calling the equation of motion using the eoms method:'))
        display(ObjectCode(
            '''
            eoms=SDOFWinchSystem()._eoms[0]
            eoms
            '''))

        display(ReportText('The result is:'))
        display(SympyFormula(eoms))


        eoms_eq=Eq(eoms,0)
        eoms_eq
        display(ReportText('The process of creating an equation from a previously called equation of motion (the equation is initially called as a string representing a description of the problem at hand, but it is not then an equation that we can use for calculations):'))
        display(ObjectCode('''
            eoms_eq=Eq(eoms,0)
            eoms_eq
            '''))
        display(ReportText('The result is:'))
        display(SympyFormula(eoms))

        slownik=SDOFWinchSystem().get_random_parameters()


        display(ReportText('Create a dictionary (in this case with random values) that allows you to assign values to variables:'))
        display(ObjectCode('''
            slownik=SDOFWinchSystem().get_random_parameters()
            slownik
            '''))
        display(slownik)

        display(ReportText('General solution of the equation:'))

        solution=SDOFWinchSystem()._ode_system.solution[0]
        display(ObjectCode(
            '''
                solution=SDOFWinchSystem()._ode_system.solution[0]
                solution

            '''))
        display(SympyFormula(solution))

        display(ReportText('Special solution of the equation:'))



        steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]
        display(ObjectCode(
            '''
                steady_solution=SDOFWinchSystem()._ode_system.steady_solution[0]
                steady_solution

            '''))
        display(SympyFormula(steady_solution))


class ODEReportComponent(ReportComponent):
    
    title="Stworzenie raportu"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('Aby utworzyć dokument, standardowo należy utworzyć sekcje lub rozdziały w następujący sposób:'))
#        display(Picture('./dynpy/utilities/components/guides/images/dok.jpg'))
        display(ObjectCode(ChapterContainer_str))
        display(ReportText('Jeżeli chcemy żeby nasze równania ładnie wyswietlały się w dokumencie nie możemy zapominac o użyciu SympyFormula:'))
        display(ObjectCode(Chapter_SympyFormula))
#        display(Picture('./dynpy/utilities/components/guides/images/sympaj.jpg'))
        display(ReportText('Co więcej, ODESystem posiada ogromną ilosć komponentów raportujących, które tworzą sekcje raportu za nas. Scieżka do ich odnalezienia jest następująca:'))
        display(ObjectCode(ode_components_path_str)) 
#        display(Picture('./dynpy/utilities/components/guides/images/sciezka_ode.jpg'))
        display(ReportText('Co więcej, istnieje też property .report, o którym informacja znajduje się rozdział niżej. Aby wygenerować plik pdf z użyciem wybranego szablonu należy zdefiniować dokuemnt jako wybrany template. Przypomnienie: Sprawdź, czy wszystkie katalogi w ścieżce poniżej istnieją. Jeśli nie, utwórz brakujące foldery przed próbą wygenerowania pliku PDF. Dobrym nawykiem jest stworzenie katalogu na pliki tekstowe wewnątrz swojego projektu. Następnie dodaj wszystkie sekcje używając metody append'))
        display(ObjectCode(property_report_formula)) 
#        display(Picture('./dynpy/utilities/components/guides/images/appeend.jpg'))

sym_num_system_code_str = ('''

m,F,c,k,x=symbols('m,F,c,k,x', positive=True)
t=Symbol('t')
x=Function(x)(t)

subs_dict={m:100 ,
          F:200 ,
          c:50 ,
          k:500,
          }

ics_list = [0.0,0.0]

''')



sym_harmonic_oscilator_eq = ('''

harmonic_oscilator_eq = Eq(m*x.diff(t,t) -F + c*x.diff(t) + k*x,0)
harmonic_oscilator_eq

''')

sym_harmonic_oscylator_ode = ('''

harmonic_oscylator_ode = ODESystem(odes= Matrix([(harmonic_oscilator_eq.lhs/m).expand()]),dvars = Matrix([x]),ode_order=2 )
harmonic_oscylator_ode

''')


sym_harmonic_oscylator_sol = ('''

harmonic_oscylator_sol = harmonic_oscylator_ode.numerized(subs_dict)
display(harmonic_oscylator_sol)
harmonic_oscylator_sym = harmonic_oscylator_sol.compute_solution(t_span=t_span, ic_list=ics_list)
display(harmonic_oscylator_sym)
harmonic_oscylator_sym.plot()

''')




sym_Example_PulledPendulum =(
'''

from dynpy.models.mechanics.pendulum import  PulledPendulum
t_span=np.linspace(0, 2, 1000)

dyn_sys_P = PulledPendulum()
display(dyn_sys_P._eoms)

subs_dict_P = dyn_sys_P.get_numerical_parameters()
display(subs_dict_P)
subs_dict_P = { **subs_dict_P, dyn_sys_P.F: 10, dyn_sys_P.Omega: 3.14/5}
for key, item in subs_dict_P.items():
    print(key)
    if key == dyn_sys_P.l**2*dyn_sys_P.m:
        key_to_pop  = key
subs_dict_P.pop(key_to_pop)
display(subs_dict_P)
ics_list_P = [0.0,0.0] # 1 stopien swobody

''')

sym_ode_sol_P_plot = ('''

ode_sol_P = dyn_sys_P._ode_system.numerized(subs_dict_P)

display(ode_sol_P)

sym_wyn_P = ode_sol_P.compute_solution(t_span=t_span, ic_list=ics_list_P)
sym_wyn_P.plot()

''')

sym_ode_sys_P =( '''

ode_sys_P = ODESystem.from_dynamic_system(dyn_sys_P.linearized())
ode_sys_P


''')

sym_ode_sol_P =( '''

ode_sol_P = ode_sys_P.solution.with_ics(ics_list)
ode_sol_P

''')

sym_num_sol_P=(
'''
num_sol_P = num_mod_P.compute_solution(t_span=t_span, ic_list=ics_list_P)

display(num_sol_P.plot())
''')

sym_coords=(
'''
display(harmonic_oscylator_ode)

display(harmonic_oscylator_sol)

coords= ode_sol_P.lhs
display(coords)
display(ode_sys_P)
display(subs_dict_P)
nadf = NumericalAnalysisDataFrame.from_model(ode_sys_P,F,[10, 20, 30], reference_data=subs_dict_P,coordinates=list(coords), index=pd.Index(t_span, name=t),ics=None)

nadf_sol = nadf.perform_simulations()

nadf_sol.columns = nadf_sol.columns.droplevel(0)
nadf_sol

''')

sym_wywolanie_helpa=(
'''
help(NumericalAnalysisDataFrame)
''')

sym_nadf=(
'''
nadf_sol.plot()
''')

sym_sms=(
"""
Omega = Symbol('Omega', positive=True)

subs_dict={m:10,
    F:100,
    c:0.01,
    k:500,
    }

t_span = np.linspace(1,2000,2001)

ode_sys_SMS = harmonic_oscylator_ode.subs({F: F*sin(Omega*t)}).subs({Omega: sqrt(k/m)+0.01})

ode_sys_SMS
""")

sym_num_wyn_SMSe=(
'''
num_wyn_SMS = ode_sys_SMS.subs(subs_dict).numerized().with_ics(ics_list).compute_solution(t_span)
''')

sym_ana_wyn_SMS=(
'''
ana_wyn_SMS = ode_sol_SMS.subs(subs_dict).with_ics(ics_list).compute_solution(t_span)
''')

sym_sms_plot=(
'''
display(num_wyn_SMS[x].plot())
''')

sym_comparison=(
"""
sym_comparison = pd.DataFrame(num_wyn_SMS[x])
sym_comparison[Symbol('x(t)_a')] = ana_wyn_SMS[x]

sym_comparison = TimeDataFrame(sym_comparison)
""")

sym_comparison_plot=(
'''
sym_comparison.to_frequency_domain().single_sided_rms().truncate(0,0.2).plot()
''')



sym_wyn_code=(
'''
sym_wyn_P = ode_sol_P.subs(subs_dict_P).compute_solution(t_span)
sym_wyn_P.plot()
''')



class ODENumericalSimulationsComponent(ReportComponent):
    
    title="Symulacje numeryczne- różne techniki symulacji"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('Przeprowadzenie symulacji przy pomocy ODE'))
        display(ReportText('W pierwszej kolejności zdefiniowano parametry potrzebne do przeprowadzenia do obliczeń'))

        display(ObjectCode(sym_num_system_code_str))

        display(ReportText('Wyprowadzenie wzoru na oscylator harmoniczny'))
        display(ObjectCode(sym_harmonic_oscilator_eq))
        #display(ObjectCode(harmonic_esc_str))

        display(ReportText('Przeprowadzenie symulacji w ODE'))


        display(ObjectCode(sym_harmonic_oscylator_ode))


        display(ObjectCode(sym_harmonic_oscylator_sol))

        #display(ObjectCode)

        #display(ObjectCode())

        display(ReportText('Wprowadzenie systemu dynamicznego, zdefiniowano również pojawiającą się Omegę '))
        display(ObjectCode(sym_Example_PulledPendulum))

        display(ObjectCode(sym_ode_sol_P_plot))
        display(ObjectCode(sym_ode_sys_P))
        display(ObjectCode(sym_ode_sol_P))

        display(ReportText('Symulacje używając systemu dynamicznego zdefiniowanego przy użyciu modułu dynpy'))

        display(ReportText('Sprawdzenie za pomocą helpa czym tak na prawdę jest NumericalAnalysisDataFrame'))
        display(ObjectCode(sym_wywolanie_helpa))

        display(ReportText('Wyświetlenie poszczególnych wykresów, wyników oraz tablicy'))
        display(ObjectCode(sym_coords))
        display(ObjectCode(sym_nadf))
        display(ReportText('Analiza porównawcza symulacji rozwiązania analitycznego i numerycznego'))
        display(ObjectCode(sym_sms))
        display(ObjectCode(sym_num_wyn_SMSe))
        display(ObjectCode(sym_ana_wyn_SMS))
        display(ObjectCode(sym_sms_plot))
        display(ObjectCode(sym_comparison))

        display(ObjectCode('''
        display(harmonic_oscylator_ode)

        display(harmonic_oscylator_sol)

        coords= harmonic_oscylator_sol.lhs

        nadf = NumericalAnalysisDataFrame.from_model(harmonic_oscylator_sol,F,[10, 20, 30], reference_data=subs_dict,coordinates=list(coords), index=pd.Index(t_span, name=t),ics=None)

        nadf_sol = nadf.perform_simulations()

        nadf_sol.columns = nadf_sol.columns.droplevel(0)

        nadf_sol
        '''))

        display(ObjectCode('''nadf_sol.plot()'''))




        display(ObjectCode('''
        Omega = Symbol('Omega', positive=True)

        subs_dict={m:10,
                   F:100,
                   c:0.01,
                   k:500,
                   }

        t_span = np.linspace(1,2000,2001)

        ode_sys_SMS = harmonic_oscylator_ode.subs({F: F*sin(Omega*t)}).subs({Omega: sqrt(k/m)+0.01})

        ode_sys_SMS
        '''))


        display(ObjectCode('''
        ode_sol_SMS = ode_sys_SMS.solution
        ode_sol_SMS.with_ics(ics_list)

        '''))

        display(ObjectCode('''num_wyn_SMS'''))


        display(ObjectCode('''ana_wyn_SMS'''))
        display(ObjectCode('''num_wyn_SMS_plot'''))
        display(ObjectCode('''ana_wyn_SMS_plot'''))

        display(ObjectCode('''
        sym_comparison = pd.DataFrame(num_wyn_SMS[x])
        sym_comparison[Symbol('x(t)_a')] = ana_wyn_SMS[x]

        sym_comparison = TimeDataFrame(sym_comparison)
        '''))

        display(ObjectCode('''
        sym_comparison.to_frequency_domain().single_sided_rms().truncate(0,0.2).plot()
        '''))

definicja_danych=(
'''
m=Symbol('m',positive=True)
g=Symbol('g')
alpha=Symbol('alpha')
v_0=Symbol('v_0')
t=Symbol("t")
x=Function("x")(t)
y=Function("y")(t)

dane = {m: 10, g: 9.81, alpha: pi/6, v_0:20}
''')

rownania_predkosci=(
'''
v_x=Eq(x.diff(t),v_0*cos(alpha))
v_y=Eq(y.diff(t),v_0*sin(alpha)-g*t)
display(v_x, v_y)
''')
definicja_ode=(
'''
ODE_x=ODESystem(odes=Matrix([v_x.lhs-v_x.rhs]),dvars=Matrix([x]),ode_order=1)
ODE_y=ODESystem(odes=Matrix([v_y.lhs-v_y.rhs]),dvars=Matrix([y]),ode_order=1)
display(ODE_x, ODE_y)
''')
ode_solution_code=(
'''
ODE_x_sol=ODE_x.subs(dane).solution.with_ics([20])
ODE_x_sol
''')
class ProjectileExampleComponent(ReportComponent):
    
    title="Przykład użycia ODESystem na podstawie rzutu ukośnego"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        from dynpy.solvers.linear import ODESystem
        m, g, v, alpha = symbols('m g v_0 alpha')

        t = Symbol('t')
        x = Function('x')(t)
        y = Function('y')(t)

        dane = {m: 10, g: 9.81, alpha: pi/6, v:20}
        predkosc_x=Eq(x.diff(t),v*cos(alpha))
        predkosc_y=Eq(y.diff(t),v*sin(alpha)-g*t)

        ode_x_x=ODESystem(odes=Matrix([predkosc_x.lhs-predkosc_x.rhs]),dvars=Matrix([x]),ode_order=1)
        ode_x_y=ODESystem(odes=Matrix([predkosc_y.lhs-predkosc_y.rhs]),dvars=Matrix([y]),ode_order=1)

        ODE_x_sol=ode_x_x.subs(dane).solution.with_ics([20])
        ODE_x_sol
        
        display(ReportText('Implementacja biblioteki ODESystem'))
        display(ObjectCode('''from sympy import *
from dynpy.solvers.linear import ODESystem
from sympy.physics.mechanics import dynamicsymbols'''))
  
        display(ReportText('Definiowanie zmiennych i utworzenie słownika'))
        display(ObjectCode(definicja_danych))
        
        display(ReportText('Równania ruchu dla osi x oraz y'))
        display(ObjectCode(rownania_predkosci))
        display(ReportText('Zapisanie równań za pomocą ODESystem'))
        display(ObjectCode(definicja_ode))

        display(ReportText('Wywołanie odpowiednio zmiennych i funkcji daje następujące rezultaty'))
        display(ObjectCode('''ODE_x'''))
        display(SympyFormula(ode_x_x))
        display(ObjectCode('''ODE_x.general_solution[0]'''))
        display(SympyFormula(ode_x_x.general_solution[0]))
        display(ObjectCode('''ODE_x.steady_solution[0]'''))
        display(SympyFormula(ode_x_x.steady_solution[0]))
        display(ObjectCode('''ODE_x.solution'''))
        display(SympyFormula(ode_x_x.solution[0]))
        display(ObjectCode('''ODE_x.solution.subs(dane)'''))

        display(SympyFormula(ode_x_x.solution.subs(dane)))#Fsimulations
        
        display(ObjectCode(ode_solution_code))

        display(SympyFormula(ODE_x_sol))
        display(ReportText('Ten sam wynik można osiągnąć stosując metoodę " dsolve ", jednak nie jest ona zalecana do równań różniczkowych'))
        display(ObjectCode('''dsolve(v_y.subs(dane))'''))
        display(SympyFormula(dsolve(predkosc_y.subs(dane))))

        
ode_import=(
'''
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from dynpy.solvers.linear import ODESystem
from sympy import Symbol, symbols, Function, Matrix
from dynpy.utilities.adaptable import *
''')

ode_variables=(
"""
t = Symbol('t')

U_z, R_w, L_w, E, U_Rw, U_Lw, k_e = symbols('U_z R_w L_w E U_Rw U_Lw k_e', positive=True)
M_s, B, J, M_obc, k_m = symbols('M_s B J M_obc k_m', positive=True)

i_w = Function('i_w')(t)
omega_s = Function('omega_s')(t)

dane = {
    U_z: 10,
    R_w: 2,
    L_w: 0.1,
    k_e: 0.1,
    k_m: 0.1,
    J: 0.1,
    B: 0.5,
    M_obc: 0
}

t_span = np.linspace(0, 1, 100)
""")

ode_objrctM=(
'''
ode = ODESystem(
    odes=Matrix([
        U_z - R_w * i_w - L_w * i_w.diff(t) - k_e * omega_s,
        k_m * i_w - J * omega_s.diff(t) - B * omega_s - M_obc
    ]),
    dvars=Matrix([i_w, omega_s]),
    ode_order=1
)
ode
''')

ode_met_1A=(
'''
ode.subs(dane).solution.with_ics([0,0])
''')

ode_tab_1A=(
'''
sym_an = ode.solution.with_ics([0,0]).subs(dane).numerized().compute_solution(t_span=t_span)
sym_an[[i_w, omega_s]]
''')

ode_met_1B=(
'''
sym_num = ode.subs(dane).numerized().compute_solution(t_span=t_span, ic_list=[0, 0])
sym_num[[i_w, omega_s]]
''')

ode_met_2=(
'''
X = [i_w, omega_s]

num_df = NumericalAnalysisDataFrame.from_model(
    ode.solution.with_ics([0, 0]),
    U_z,
    [10, 20, 50],
    reference_data=dane,
    coordinates=X,
    index=pd.Index(t_span, name=t)
)
num_df
''')

ode_met_PS=(
'''
sym_wyn = num_df.perform_simulations().droplevel(0,axis=1)

display(sym_wyn)
''')

ode_import_example=(
'''
from dynpy.models.mechanics.trolley import SpringDamperMassSystem

dyn_sys = SpringDamperMassSystem()

dyn_sys._ode_system
''')

ode_example2=(
'''
dane = dyn_sys.get_numerical_parameters()

data_dict = {**dane, dyn_sys.F: 1000}

data_dict
''')

ode_simulation=(
'''
simulation = dyn_sys._ode_system.subs(data_dict).numerized().compute_solution(
    np.linspace(0, 2, 1000), ic_list=[10, -5]
)

simulation
''')

ode_simulation_plot=(
'''
simulation.plot()
plt.grid()
plt.show()
''')




        
class ODESimulationComponent(ReportComponent):
    
    title="Przeprowadzanie symulacji"
    


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('To od czego standardowo powinno się zacząć to importy oraz definicja zmiennych. Potrzebne importy:'))
        display(ObjectCode(ode_import))
        display(ReportText('Definiowanie zmiennych:'))
        display(ObjectCode(ode_variables))
        display(ReportText('Kolejny niezbędny krok to tworzenie obiektu ODE. Należy go wykonać w następujący sposób:'))
        display(ObjectCode(ode_objrctM))
        display(ReportText('METODA 1A: Przeprowadzenie symulacji na gotowm rozwiazaniu analitycznym poprzez użycie property solution. Samo wywołanie                                 odpowiedniego rozwiązania analitycznego pokazuje nam poniższy kod:'))
        display(ObjectCode(ode_met_1A))
        display(ReportText('Następnie generujemy tabele niezbędną do symulacji w taki sposób:'))
        display(ObjectCode(ode_tab_1A))
        display(ReportText('METODA 1B: Innym sposobem jest przeprowadzenie symulacji numerycznych bez obliczeń analitycznych. Potrzebny jest do tego taki kod: '))
        display(ObjectCode(ode_met_1B))
        display(ReportText('METODA 2: Trzecim sposobem generowania symulacji jest użycie NumericalAnalysisDataFrame, kod wywołania prezentuje się następująco:'))
        display(ObjectCode(ode_met_2))
        display(ReportText('Aby przeprowadzić symulacje konieczne jest zastosowanie metody perform simulations. Bez zastosowania tej metody dostaniemy następujący wynik:'))
        display(ReportText('A z jej użyciem, wszystko powinno zadziałać i dostaniemy następujący output:'))
        display(ObjectCode(ode_met_PS))
        display(ReportText('Dla lepszego zrozumienia zagdnień poniżej przedstawione zostało zastosowanie symulacji na konkretnym przykładzie. Skorzystamy z                           prostej klasy SpringDamperMassSystem. Na dobry początek należy zaimportować tę klasę oraz przypisać ją do zmiennej. W przeprowadzeniu                                     symulacji pomoże nam wykorzystanie metody $ode system$ w celu pozyskania równania ruchu naszego systemu dynamicznego. Aby pozyskać dane w                                 prosty sposób korzystamy z metody $get numerical parameters$ oraz wrzucamy te dane w słownik. Kod do wykonania tych czynnosci prezentuje się                             następująco: '))
        display(ObjectCode(ode_import_example))
        display(ReportText('Następnie stosując jedną z wyżej pokazanych metod, możemy stworzyć tabelę:'))
        display(ObjectCode(ode_example2))
        display(ReportText('Ostatnim krokiem jest wizualizacja danych, czyli stworzenie wykresu. Aby taki wykres uzyskać wystarczy wpisać kod widoczny na                             obrazku:'))
        display(ObjectCode(ode_simulation))
        display(ObjectCode(ode_simulation_plot))


class_call_str=(
'''
from dynpy.models.mechanics.principles import ComposedSystem
from sympy import *

class MyMaterialPointMovement(ComposedSystem):
    pass
''')

atributes_str=(
'''
m = Symbol('m', positive=True)
g = Symbol('g', positive=True)
c = Symbol('c', positive=True)
r = Symbol('r', positive=True)
phi = dynamicsymbols('phi')

c0 = Symbol('c0', positive=True)
r0 = Symbol('r0', positive=True)
phi0 = dynamicsymbols('phi0')
''')

innit_str=(
'''
def __init__(self,
             m=None,
             g=None,
cLatexDataFrame.set_default_units(unit_dict)

tabelka=LatexDataFrame.formatted(
    data=dane,
    index=['Value']).applymap(lambda x: f'${latex(x)}$').rename_axis(
        'Symbol',
        axis=1)

display(tabelka.reported(caption='Data given in the exercise'))

''')

component_str=(
'''
@property
def components(self):

    components = {}


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



    components['_mass_x']=self._mass_x
    components['_mass_y']=self._mass_y
    components['_gravity_']=self._gravity_


    return components
''')


symb_desc_str=(
'''
def symbols_description(self):
    self.sym_desc_dict = {
        self.m: r'Mass',
        self.g: r'Gravity constant',
        self.c: r'',
    }

    return self.sym_desc_dict
''')

def_data_str=(
'''
def get_default_data(self):

    m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

    default_data_dict = {
        self.m: [m0 * no for no in range(1, 8)],
        self.c: [c0 * no for no in range(1, 8)],
        self.r: [r0 * no for no in range(1, 8)],
        self.phi: [phi0 * no for no in range(1, 8)],
    }

    return default_data_dict
''')

num_data_str=(
'''
def get_numerical_data(self):

    m0, c0, r0, phi0 = self.m0, self.c0, self.r0, self.phi0

    default_data_dict = {
        self.m: [m0 * no for no in range(1, 8)],
        self.c: [c0 * no for no in range(1, 8)],
        self.r: [r0 * no for no in range(1, 8)],
        self.phi: [phi0 * no for no in range(1, 8)],
    }

    return default_data_dict
''')

dict_str=(
'''
def unit_dict(self):

    from pint import UnitRegistry
    ureg=UnitRegistry()

    unit_dict = {
        self.m: ureg.kilogram,
        self.g: ureg.meter/ureg.second/ureg.second,
        self.c: ureg.kilogram/ureg.second,
        self.r: ureg.meter,
        self.phi: ureg.radian,
        self.c0: ureg.kilogram/ureg.second,
        self.r0: ureg.meter,
        self.phi0:  ureg.radian
    }
    return unit_dict
''')

full_class_str=(
'''
from dynpy.models.mechanics.principles import ComposedSystem
from sympy import *

class MyMaterialPointMovement(ComposedSystem):

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

        self._gravity_ = GravitationalForce(self.m,
                                            self.g,
                                            pos1=self.r * cos(self.phi),
                                            qs=self.qs)


      
        components['_mass_x']=self._mass_x
        components['_mass_y']=self._mass_y
        components['_gravity_']=self._gravity_
     

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
        
    def unit_dict(self):

        from pint import UnitRegistry
        ureg=UnitRegistry()

        unit_dict = {
            self.m: ureg.kilogram,
            self.g: ureg.meter/ureg.second/ureg.second,
            self.c: ureg.kilogram/ureg.second,
            self.r: ureg.meter,
            self.phi: ureg.radian,
            self.c0: ureg.kilogram/ureg.second,
            self.r0: ureg.meter,
            self.phi0:  ureg.radian
        }
        return unit_dict        
''')

class DynSysImplementationComponent(ReportComponent):
    
    title="Implementacja systemow dynamicznych"


    def append_elements(self):
        
        #system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        #przekopiowac tu fragment kodu z guide pamietajac o indentach
        display(ReportText("Tworzenie klasy należy zacząć od wybrania nazwy opisującej nasz system oraz zadeklarowania jej w odpowiedni sposób, który został przedstawiony poniżej:")) #itp
        display(GuideCode(class_call_str))
        display(ReportText('Następnie deklarujemy atrybuty/symbole, które są niezbędne do opisania danego obiektu'))
        display(GuideCode(atributes_str))
        display(ReportText('Do zdefiniowania inicjalizacji obiektu klasy z określonymi atrybutami służy poniższa struktura.'))
        display(GuideCode(innit_str))
        display(ReportText('W tym kroku korzystamy z predefiniowanych elementów takich jak punkty materialne, dyski, itp.'))
        display(GuideCode(component_str))
        display(ReportText('W tej części opisujemy parametry klasy, np. masę'))
        display(GuideCode(symb_desc_str))
        display(ReportText('Sekcja ta przedstawia inicjalizację słownika zawierającego dane symboliczne do równań analitycznych systemu:'))
        display(GuideCode(def_data_str))
        display(ReportText('Kolejnym krokiem jest przypisanie słowników zawierających dane do symulacji numerycznych:'))
        display(GuideCode(num_data_str))
        display(ReportText('Ta sekcja opisuje definiowanie jednostek używanych przez system dynamiczny:'))
        display(GuideCode(dict_str))
        
        display(ReportText('Ostatecznie kod implementacji całego systemu ma następującą postać:'))
        display(GuideCode(full_class_str))


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




SympyFormulaComponent_str=(
'''
from sympy import *

M = Symbol('M')

t = Symbol('t')

m_e = Symbol('m_e')

e = Symbol('e')

z = Function('z')

varphi = Function('varphi')

Ek = M*Derivative(z(t), t)**2/2 + m_e*(e*sin(varphi(t))*Derivative(varphi(t), t) - Derivative(z(t), t))**2/2

display(Eq(Symbol('T'),Ek))
''')


class SympyFormulaComponent(DocumentGenerationComponent):

    title = "Klasa SympyFormula i jej zastosowanie"



    def append_elements(self):
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        #implement reporting activieties here
        display(ReportText('Przykładowe zastosowanie komponentu do obliczeń symbolicznych, na przykładzie energii kinetycznej:'))
        display(GuideCode(f'{SympyFormulaComponent_str}'))


ODEIcsComponent_dec_str = ('''
m, c, k, F, mu = symbols('m c k F mu', positive = True)
t = Symbol('t')
x = Function('x')(t)

subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))
display(ode_1)
''')
        
class ODEIcsComponent(ReportComponent):

    title = "ODESystem solution with set of initial conditions (ICS) component presentation"

    def append_elements(self):

        display(ReportText('First declare ODESystem component:'))
        
        display(ObjectCode(ODEIcsComponent_dec_str))
        
        from dynpy.solvers.linear import AnalyticalSolution, ODESolution, ODESystem
        import sympy
        
        m, c, k, F, mu = symbols('m c k F mu', positive = True)
        t = Symbol('t')
        x = Function('x')(t)

        subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

        ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))
        
        display(ReportText('After executing above code bellow equation will appear'))
        
        display(ode_1)

        display(ReportText('To present solution with initial set of conditions, first you need to calculate solution:'))
        
        display(ObjectCode('sol_sym_1_pre_ics = ode_1.subs(subs_data).solution'))
        
        sol_sym_1_pre_ics = ode_1.subs(subs_data).solution
        
        display(sol_sym_1_pre_ics)
        
        display(ReportText('To calculate solution with ICS, declare ics vector, and run with_ics method as shown bellow:'))
        
        display(ObjectCode('''
ics = [1, 2]
sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
display(sol_sym_1)'''))
        
        ics = [1, 2]
        sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
        display(sol_sym_1)
        
ODEIcsCodeComponent_str =('''
m, c, k, F, mu = symbols('m c k F mu', positive = True)
t = Symbol('t')
x = Function('x')(t)

subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))

sol_sym_1_pre_ics = ode_1.subs(subs_data).solution

ics = [1, 2]

sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
''')

class ODEIcsCodeComponent(ReportComponent):

    title = "ODESystem solution with set of initial conditions (ICS) component presentation"

    def append_elements(self):

        display(ReportText('For bellow ODESystem component example:'))
        
        display(ObjectCode(ODEIcsCodeComponent_str))
        
        from dynpy.solvers.linear import AnalyticalSolution, ODESolution, ODESystem
        import sympy
        
        m, c, k, F, mu = symbols('m c k F mu', positive = True)
        t = Symbol('t')
        x = Function('x')(t)

        subs_data = {m: 10, c: 10, k: 10, F: 5, mu: 2}

        ode_1 = ODESystem(odes = Matrix([m*x.diff(t,t) + c *x.diff(t)  + k * x]), odes_rhs = Matrix([F*sin(mu * t)]), dvars = Matrix([x]))

        sol_sym_1_pre_ics = ode_1.subs(subs_data).solution

        ics = [1, 2]

        sol_sym_1 = sol_sym_1_pre_ics.with_ics(ics)
        
        display(ReportText('To get call code for ICS, run bellow method:'))
        
        display(ObjectCode("sol_sym_1.get_ics()"))
        
        display(sol_sym_1.get_ics())
        


ds_overview_str=("""
from dynpy.utilities.documents.document import DynSysOverviewReport
from dynpy.utilities.components.guides.en import *
from dynpy.models.mechanics.pendulum import Pendulum

DynSysOverviewReport(reported_object=Pendulum());
""")

dynsys_comp_check_str = '''from dynpy.models.mechanics import ForcedSpringMassSystem

comp = ForcedSpringMassSystem()

DynamicSystemCompletenessCheckComponent(comp)'''

basicsymcomp_str = '''from dynpy.models.mechanics.tmac import SDOFWinchSystem

comp = SDOFWinchSystem()

BasicSymComponent(comp)'''

coderefactorcomp_str = '''from dynpy.utilities.components.guides.en import CodeRefactorComponent

CodeRefactorComponent(None)'''

dynsysimpcomp_str = '''from dynpy.utilities.components.guides.pl import DynSysImplementationComponent

DynSysImplementationComponent(None)'''

diffsimcomp_str = '''from dynpy.utilities.components.guides.en import DifferentSimulationsComponent

comp = SDOFWinchSystem()

DifferentSimulationsComponent(comp)'''

dynintrocomp_str = '''from dynpy.utilities.components.guides.en import *
DynSysIntroComponent(None)'''

dynsumcomp_str = '''from dynpy.utilities.components.guides.en import *

DynSysSummaryComponent(None)'''

dynmetuscomp_str = '''from dynpy.utilities.components.guides.en import DynamicSystemMethodsUsageComponent

comp = SDOFWinchSystem()

DynamicSystemMethodsUsageComponent(comp)'''
dyncompov_str = '''Guide consist of bellow components which can be called separatelly with most of avaible dynsys models:
*DynamicSystemCallComponent
*NumericalAnalysisSimulationComponent
*AnalyticalSimulationComponent
*DynSysCodeComponent
*IssuePreparationComponent'''

class DynSysOverviewUsageComponent(ReportComponent):

    title = "Introduction to usage of DynSysOverviewReport"

    def append_elements(self):

        from dynpy.models import mechanics


        display(ReportText('This guide concers in-deepth analysis of exsisting dynamic system that are basic elements of `DynPy` library.'))
        display(ReportText('Basic call of it is as follows and runs default dynamic system which is `ForcedSpringMassSystem'))
        display(ObjectCode(ds_overview_str))
        
        display(ReportText(dyncompov_str))
        
        #display(ReportText('This guide concers in-deepth analysis of exsisting dynamic system that are basic elements of `DynPy` library.'))
        #display(ReportText('It can be run with any available dynamic system.'))
        '''
        #DynamicSystemCompletenessCheckComponent
        display(ReportText('This component calls and checks all key elements of the dynamic system'))
        display(ReportText('It can be run with any available dynamic systems.'))
        display(ObjectCode(dynsys_comp_check_str))
        
        #BasicSymComponent
        display(ReportText('This component shows how to do basic simulations for dynamic system.'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on below example call.'))
        display(ObjectCode(basicsymcomp_str))
        
        #CodeRefactorComponent
        display(ReportText('This component shows how to refactor DynSys component according to latest standard.'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(coderefactorcomp_str))
        
        #DifferentSimulationsComponent
        display(ReportText('This component explains how to implement dynamic system with ComposedSystem class'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(diffsimcomp_str))
        
        #DynSysImplementationComponent
        display(ReportText('This component explains how to implement dynamic system with ComposedSystem class'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on below example call.'))
        display(ObjectCode(dynsysimpcomp_str))
        
        #DynSysIntroComponent
        display(ReportText('This component is simple introduction to using the DynSys module'))
        display(ReportText('It can be run for SDOFWinchSystem, as shown on baeow example call.'))
        display(ObjectCode(dynintrocomp_str))
            
        #DynSysSummaryComponent
        display(ReportText('This component shows idea of code refactoring and additionally displays list of available dynamical systems'))
        display(ReportText('It is called without any argument.'))
        display(ObjectCode(dynsumcomp_str))
        
        #DynamicSystemMethodsUsageComponent
        display(ReportText('This component shows equations of motion and other methods avaible on DynSys model'))
        display(ReportText('It can be run with any available dynamic systems.'))
        display(ObjectCode(dynmetuscomp_str))
        '''


class DynamicSystemCompletenessCheckComponent(ReportComponent):

    title="Wywoływanie i sprawdzanie wszystkich kluczowych elementów systemu dynamicznego"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        eoms=system._eoms[0]
        eoms

        display(ReportText('Wywołanie systemu w pierwszej kolejności można sprawdzić przez sprawdzenie równania ruchu metodą _eoms.'))

        display(ReportText('Wynik jest następujący:'))
        
        display(Markdown(
            '''
        eoms_eq=Eq(eoms,0)
        eoms_eq
            '''))
        
        eoms_eq=Eq(eoms,0)
        eoms_eq
        
        display(SympyFormula(eoms_eq))
        
        display(ReportText('Następnie wykonuje się analizę metody __init__ klasy:'))
        
        display(ObjectCode(system.__init__))
        
        display(ReportText('Kolejnym etapem jest sprawdzenie, czy schemat oraz zdjęcie reczywiste systemu są zdefiniowane:'))
        
        display(Markdown(
            '''
        system._as_picture()
        system.preview(example = 'detail_real_name')
            '''))
        
        system._as_picture()
        system.preview(example = 'detail_real_name')
        
        display(ReportText('Kolejnym etapem jest sprawdzenie, czy zdefiniowane są słowniki z parametrami:'))
        
        display(Markdown(
            '''
        system.get_random_parameters()
        system.get_numerical_parameters()
            '''))
        
        display(ReportText('Jeśli podczas wywołania pojawia się błąd o braku tej metody lub parametry są inne niż w metodzie __init__ należy zmodyfikować klasę analogicznie do pokazanego rozwiązania:'))
        
        display(Markdown(
            '''
        
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
            '''))
        
        system.get_random_parameters()
        system.get_numerical_parameters()
        
        display(ReportText('Następnie konieczne jest sprawdzenie opisów symboli oraz jednostek:'))
        
        display(Markdown(
            '''
        system.symbols_description()
        system.unit_dict()
            '''))
        
        display(ReportText('Jeśli podczas wywołania tej metody pojawia się błąd, należy wykonać modyfikację w analogizny sposób do:'))
        
                
        display(Markdown(
            '''
        def symbols_description(self):
            self.sym_desc_dict = {
                self.m: r'system's mass',
                self.k_l: r'left spring stiffness',
                self.k_r: r'left spring stiffness',
            }
            return self.sym_desc_dict

        def unit_dict(self):
        
            from pint import UnitRegistry
            ureg=UnitRegistry()

            unit_dict = {
                self.m: ureg.kilogram,
                self.k_l: ureg.Newton/ureg.meter,
                self.k_r: ureg.Newton/ureg.meter,
            }
            return unit_dict

            '''))
        
        from pint import UnitRegistry
        ureg=UnitRegistry()
        
        system.symbols_description()
        system.unit_dict()

