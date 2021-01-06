import numpy as np
import pandas as pd

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Alignat, Package, Quantity, Command, Label, Ref, Marker, NewPage, Eqref, Table
from pylatex.section import Chapter
from pylatex.utils import italic, NoEscape

from pylatex.base_classes import Environment
from pylatex.package import Package

import os

from sympy.physics.vector.printing import vpprint, vlatex
import sympy.physics.mechanics as me

import matplotlib.pyplot as plt
import random
from sympy import *


class Equation(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package('mathtools')]
    escape = False
    content_separator = "\n"


class DMath(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package('breqn'), Package('flexisym')]
    escape = False
    content_separator = " "

# class EqRef(Environment):
#         """A class to wrap LaTeX's alltt environment."""

#     packages = [Package('breqn'), Package('flexisym')]
#     escape = False
#     content_separator = " "

class ReportModelAnalysis:
    
    chapter_name = 'Model analityczny układ o swobody'
    
    def __init__(self,model=None):
        geometry_options = {
            "margin": "1cm",
        }
        doc = Document(documentclass='report',
                       geometry_options=geometry_options)
        #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

        doc.packages.append(Package('fontenc'))
        #doc.packages.append(Package('fontspec'))
        doc.packages.append(Package('standalone'))
        doc.packages.append(Package('siunitx'))
        doc.packages.append(Package('amsmath'))
        doc.packages.append(Package('float'))
        doc.packages.append(Package('tikz'))
        doc.packages.append(Package('pgfplots'))        
        #doc.packages.append(Package('preamble'))
        doc.packages.append(Package('polski', options=["MeX"]))
        doc.packages.append(Command('usepgfplotslibrary', arguments='units'))
#         doc.packages.append(Package('preamble'))
        self.doc = doc
        
        self.model = model
        
        self.chapter_name= type(self).chapter_name

    def lagranges(self):

        doc = self.doc

        model = self.model
        
        #musisz też dopisać o jednym stopniu - to będzie ten indeks
        with doc.create(Chapter(self.chapter_name)):

            with doc.create(DMath()):
                doc.append('T=')

                doc.append(vlatex(model.lagrangian() +model.lagrangian().subs({vel:0  for vel  in model.u})))
            with doc.create(DMath()):
                doc.append('V=')

                doc.append(vlatex(-self.model.lagrangian().subs({vel:0  for vel  in model.u})))

    

            with doc.create(DMath()):
                doc.append('L=')

                doc.append(vlatex(self.model.lagrangian()))
        return doc
    def stat_and_dyn(self):
        doc = self.doc

        model = self.model

        with doc.create(Chapter(self.chapter_name)):

            for i in range(len(model.q)):
                with doc.create(DMath()):
                    doc.append(vlatex(self.model.equilibrium_equation().doit()[i]))
            for lhs,rhs in zip(symbols('k_1:{ndof}_1:{ndof}'.format(ndof=len(model.q)+1)),(list(model.linearized().stiffness_matrix()))):
                with doc.create(DMath()):
                    doc.append(vlatex([Eq(lhs,rhs)]))  
            for i in range(len(model.q)):
                with doc.create(DMath()):
                    doc.append(vlatex(self.model._eoms[i]))
            for lhs,rhs in zip(symbols('m_1:{ndof}_1:{ndof}'.format(ndof=len(model.q)+1)),(list(model.inertia_matrix()))):
                with doc.create(DMath()):
                    doc.append(vlatex([Eq(lhs,rhs)]))  
        return doc


    def doit(self):

        doc = self.doc

        model = self.model
        

        with doc.create(Chapter(self.chapter_name)):
            doc.append(
                '''Model przyjęty do weryfikacji z danymi doświadaczalnymi stanowił złożenie modelu o trzech stopniach oraz układu napędowego o jednym stopniu swobody. Funkcja Lagrangea tego układu miała następującą postać:'''
            )

            with doc.create(DMath()):
                doc.append('L=')

                doc.append(vlatex(model.lagrangian()))

            doc.append(
                '''Bazując na wyznaczonej funkcji Lagrange'a określono równania równowagi (specjalny przypadek równań Lagrange'a drugiego rodzaju")'''
            )

            for no, balance_eq in enumerate(model.equilibrium_equation()):
                with doc.create(DMath()):

                    #doc.append(vlatex(chair_dict['linearized_ramp'].q[no])+ ':~~~~'  +  vlatex(Eq(balance_eq.doit(),0) )   )
                    doc.append(vlatex(Eq(balance_eq.doit(), 0)))

            doc.append(
                '''Bazując na wyznaczonej funkcji Lagrange'a określono równania równowagi (specjalny przypadek równań Lagrange'a drugiego rodzaju)'''
            )

            for no, eom in enumerate(model._eoms):
                with doc.create(DMath()):

                    #doc.append(vlatex(chair_dict['linearized_ramp'].q[no])+ ':~~~~'  +  vlatex(Eq(balance_eq.doit(),0) )   )
                    doc.append(vlatex(Eq(eom.doit(), 0)))

        return doc


    
class RandomDescription:
    

    def __init__(self,*args,random_or_order=None,**kwargs):
        
        self.pool_list=args
        self.random_or_order=random_or_order
        

    def to_string(self):
        
        if hasattr(self.random_or_order, '__contains__'):
            self.order=list(random_or_order)
        else:
            self.selected_sentences=[random.choice(pool)  for pool in self.pool_list]
        
        
        
        
        return ' '.join(self.selected_sentences)
    
    def __str__(self):
        
        
        return self.to_string()
    
    def __repr__(self):
        return self.__class__.__name__+': '+self.to_string()

    
    
    
time_domain_summary_1 = [
'Przedstawione wyniki pokazują wpływ zmian badanego parametru \({param_name}\) na dynamikę wózka. ',
'Opracowane wykresy porównawcze pozwoliły na oszacowanie wpływu zmiennej \({param_name}\) na dynamikę rozpatrywanego układu. ',
'Zaprezentowane wyniki posłużyły ocenie wrażliwości modelowanego obiektu na badaną zmienną \({param_name}\). ',
'Przygotowane przebiegi służą zaprezentowaniu wyników symulacji numerycznych badanej zmiennej w funkcji wybranego parametru wpływowego \({param_name}\). ',
'Analizowany model wykazuje wrażliwość na zmianę rozważanego parametru \({param_name}\). ',
]

time_domain_summary_2 = [
'Zaobserwowano wpływ badanego parametru \({x(t)}\) na drgania wózka, których maksymalna wartość amplitudy wyniosła {x(t)_max}. ',
'Analizowany przebieg parametru \({x(t)}\) osiąga wartość maksymalną równą {x(t)_max}. '
'Minimalna wartość dla przebiegu \({x(t)}\) jest równa {x(t)_min}. '
'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({x(t)}\) nie występują amplitudy o wartości mniejszej niż {x(t)_min} oraz większej niż {x(t)_max}. '
'Dokonując analizy wyników symulacji dla parametru \({x(t)}\) stwierdza się, że amplitudy drgań nie przekraczają wartości {x(t)_max}, a wartość mininalna wynosi {x(t)_min}. ',
]


time_domain_summary_3 = [ 
'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi(t)}\) wynosi {varphi(t)_max}. ',
'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi(t)_max} dla współrzędnej \({varphi(t)}\). '
'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {varphi(t)_min}. ',
'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
]


time_domain_summary_4 = [ 
'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({z(t)}\) wynosi {z(t)_max}. ',
'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {z(t)_max} dla współrzędnej \({z(t)}\). '
'Wartość minimalna przebiegu współrzędnej \({z(t)}\) jest równa {z(t)_min}. ',
'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({z(t)}\) nie przekracza {z(t)_max} i jest ograniczona z dołu przez {z(t)_min}. ',
'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({z(t)}\) nie przekraczaa wartości {z(t)_max}, a wartość mininalna wynosi {z(t)_min}. ',
]


time_domain_summary_5 = [ 
'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi_RC(t)}\) wynosi {varphi_RC(t)_max}. ',
'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi_RC(t)_max} dla współrzędnej \({varphi_RC(t)}\).'
'Wartość minimalna przebiegu współrzędnej \({varphi_RC(t)}\) jest równa {varphi_RC(t)_min}. ',
'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi_RC(t)}\) nie przekracza {varphi_RC(t)_max} i jest ograniczona z dołu przez {varphi_RC(t)_min}. ',
'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi_RC(t)}\) nie przekraczaa wartości {varphi_RC(t)_max}, a wartość mininalna wynosi {varphi_RC(t)_min}. ',
]


time_domain_summary_7 = [
'Zaobserwowane odczyty wynikające z wpływu \({param_name}\) występują w konketrnym czasie ruchu {t_val}. ',
'Odczytana maksymalna wartość amplitudy {max_val} wynika z działania \({param_name}\) wystąpiła konkretnie w czasie {t_val}. '
'Minimalny zaobserwowany wynik amplitudy {min_val} miał miejsce w {t_val} czasu trwania pomiaru. ',
]

time_domain_summary_8 = [ 
'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionych są następujące: {max_val_list}. ',
'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się następującymi wartościami maksymalnymi: {max_val_list}. '
'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {min_val}. ',
'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
]



freq_domain_summary =[
'Ponadto sprawdzono własności spektralne otrzymanych wyników symulacyjnych.',
'Według przeprowadzonych badań modelowych struktura częstotliwościowa nie zawiera komonentów większych niż {max_val_spec}.',
'Amplitudy drgań badanego parametru \({q_name}\) przyjmują maksymalną wartość równą {max_val_spec}.',
'Na podstawie analizy widmowej minimalna wartość częstotliwości rozpatrywanego parametru wynosi {min_val_spec}.',


]


summary_bank_composed_model=[str(RandomDescription(
        time_domain_summary_1,
        time_domain_summary_2,
        time_domain_summary_3,
        time_domain_summary_4,
        time_domain_summary_5,
            )) for obj in range(30)  ]    



sec_introduction_bank_1=[
'Wykorzystując przygotowany model dynamiczny, przeprowadzono serię symulacji mających na celu ocenę wrażliwości modelu na zmiany wybranych parametrów. ',
'Przeprowadzenie analizy numerycznej wymaga znajomości modelu dynamicznego oraz zdeterminowania parametrów oddziałujących na zachowanie układu wskutek ich zmian. ',
'Wykorzystując przygotowane środowisko obliczeniowe, wykonano szereg symulacji, opracowano dane numeryczne i przedstawiono je w postaci następujących wykresów. ',
'Celem zbadania i oceny dynamiki rozpatrywanego systemu przeprowadzono symulacje numeryczne przedstawijące przebiegi zmiennych układu w funkcji czasu. '
'W oparciu o wybrany model dynamiczny wózka, wykonano analizy numeryczne zaprezentowane na kolejnych wykresach. '
]

sec_introduction_bank_2=[
'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
]


sec_introduction_bank_3=[
'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
]


sec_intro_composed_model=[str(RandomDescription(sec_introduction_bank_1,sec_introduction_bank_2,sec_introduction_bank_3)) for obj in range(30)  ]  

introduction_bank_1=[
'Na wykresie {nr_rys} przedstawiano zmiany wielkości dynamicznych charakteryzujących ruch wózka. Po przeanalizowaniu można zanotować wzajemną zależność poszczególnych wielkości dynamicznych.',
'Charakter przebiegu wartości dynamicznych wózka został przedstawiony na rysunku {nr_rys}.',
'Wykres {nr_rys} pokazuje zmienność parametrów dynamicznych wózka w trakcie symulownego przejazdu.',
'Rysunek {nr_rys} prezentuje wykres opisujący zmienność parametrów dynamicznych wózka w trakcie jazdy po zamodelowanym torze.',
]

introduction_bank_2=[
'Symulacja została przeprowadzona dla następujących danych: {given_data}.',
'Zaprezentowane wyniki numeryczne otrzymano dla danych równych:{given_data}.',
'Eksperyment numeryczny został przeprowadzona dla następujących danych:{given_data}.',
'Wyniki stmulacji numerycznych otrzymano dla:{given_data}.',
]



introduction_bank_3=[
'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu stwierdza się oscylacyjny charakter odpowiedzi układu. Dynamika systemu ma stabilny charakter',
'Charakter zmian przebiegów poszczególnch współrzędnych jest stabilny. Otrzymane sygnały mieszczą się w zakresie dopuszczalnych uwarunkowaniami fizycznymi.',
'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu mają drganiowy charakter i są zgodne z oczekiwaniami.',
'Bazując na otrzymanych wynikach można stwierdzić wzajemną zależność poszczególnych wielkości dynamicznych. Dodatkowo wielkości modelu stabilizują się w czasie. '
]



intro_bank_composed_model=[str(RandomDescription(introduction_bank_1,introduction_bank_2,introduction_bank_3)) for obj in range(30)  ]


class PlottedData(Figure):
    
    _latex_name='figure'
    
    def __init__(self,numerical_data,fig_name,*, preview=False, position=None, **kwargs):
        super().__init__(position=position, **kwargs)

        self._numerical_data=numerical_data
        self.fig_name=str(fig_name)
        #self._latex_name='figure' #super()._latex_name
        self.preview=preview


        
        
    def add_data_plot(self,numerical_data=None):
        numerical_data=self._numerical_data

        ax = numerical_data.plot(subplots=True, figsize=(10, 7))
        #ax=solution_tmp.plot(subplots=True)
        ([
            ax_tmp.legend(['$' + vlatex(sym) + '$'])
            for ax_tmp, sym in zip(ax, numerical_data.columns)
        ])
        plt.savefig(self.fig_name+'.png')
        self.add_image(self.fig_name,width=NoEscape('15cm'))

        if self.preview==True:
            plt.show()

        plt.close()
        



class ReportSection(Chapter):
    
    _latex_name='section'
    
    def __init__(self,
                 results_frame,
                 analysis_key,
                 analysis_name='dynamics',
                 results_col_name='simulations',
                 preview=True):
        
        self.analysis_key = analysis_key
        self.analysis_name = analysis_name
        self.sims_name = results_col_name
        self.preview=preview


class ReportSimulationsAnalysis:

    chapter_name='Analiza wrażliwości dla parametru \({param_name}\)'
    summary_sec_name='Spostrzeżenia i wnioski z analizy wpływu parametru  \({param_name}\)'
    
    sensitivity_level={0.02:'znikomy wpływ',0.1:'mały wpływ',0.5:'średni wpływ',0.7:'średni wpływ',0.9:'istotny wpływ'}
    
    units_dict={}
    
    sec_intro_bank= [str(RandomDescription(
                ['Section description 1'+str(i) for i in range(5)],
                ['Section description 2'+str(i) for i in range(5)],
                ['Section description 3'+str(i) for i in range(5)],
                ['Section description 4'+str(i) for i in range(5)],
                ['Section description 5'+str(i) for i in range(5)],
                    )) for obj in range(30)  ]



    intro_bank=[str(RandomDescription(
                ['Intro 1'+str(i) for i in range(5)],
                ['Intro 2'+str(i) for i in range(5)],
                ['Intro 3'+str(i) for i in range(5)],
                ['Intro 4'+str(i) for i in range(5)],
                ['Intro 5'+str(i) for i in range(5)],
                    )) for obj in range(30)  ]


    summary_bank=[str(RandomDescription(
                ['Summary 1'+str(i) for i in range(5)],
                ['Summary 2'+str(i) for i in range(5)],
                ['Summary 3'+str(i) for i in range(5)],
                ['Summary 4'+str(i) for i in range(5)],
                ['Summary 5'+str(i) for i in range(5)],
                    )) for obj in range(30)  ]

    
    conclusion_bank=[str(RandomDescription(
#                 ['{varp(t)_max} dla argumentu {x(t)_idmax} 1'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 2'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 3'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 4'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 5'+str(i) for i in range(5)],
                ['Conclusion 1'+str(i) for i in range(5)],
                ['Conclusion 2'+str(i) for i in range(5)],
                ['Concluison 3'+str(i) for i in range(5)],
                ['Conclusion 4'+str(i) for i in range(5)],
                ['Conclusion 5'+str(i) for i in range(5)],
                    )) for obj in range(30)  ]
    
    optimization_desc_bank=[str(RandomDescription(
#                 ['{varp(t)_max} dla argumentu {x(t)_idmax} 1'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 2'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 3'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 4'+str(i) for i in range(5)],
#                 ['{x(t)_max} dla argumentu {x(t)_idmax} 5'+str(i) for i in range(5)],
                ['Opti 1'+str(i) for i in range(5)],
                ['Opti 2'+str(i) for i in range(5)],
                ['Opti 3'+str(i) for i in range(5)],
                ['Opti 4'+str(i) for i in range(5)],
                ['Opti 5'+str(i) for i in range(5)],
                    )) for obj in range(30)  ]    
    
    
    def __init__(self,
                 results_frame,
                 analysis_key,
                 analysis_name='dynamics',
                 results_col_name='simulations',
                 preview=True):
        geometry_options = {
            "margin": "1cm",
        }
        doc = Document(documentclass='report',
                       geometry_options=geometry_options)
        #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

        doc.packages.append(Package('fontenc'))
        #doc.packages.append(Package('fontspec'))
        doc.packages.append(Package('standalone'))
        doc.packages.append(Package('siunitx'))
        doc.packages.append(Package('amsmath'))
        doc.packages.append(Package('float'))
        doc.packages.append(Package('tikz'))
        doc.packages.append(Package('booktabs'))
        doc.packages.append(Package('pgfplots'))
        doc.packages.append(Package('polski', options=["MeX"]))
        doc.packages.append(Command('usepgfplotslibrary', arguments='units'))

        self.doc = doc
        self.analysis_key = analysis_key
        self.analysis_name = analysis_name
        self.sims_name = results_col_name
        self.preview=preview
        
        #self.model= model

        self.sec_intro_bank= type(self).sec_intro_bank
        

        
        self.intro_bank=type(self).intro_bank
        
       
        self.summary_bank=type(self).summary_bank
        self.conclusion_bank=type(self).conclusion_bank

        self.sensitivity_level=type(self).sensitivity_level

        


        
        self.caption='Przebiegi czasowe modelu dla danych: {given_data}'
        # 'Wartość największa wynosi: {max_val}'.format(**{'max_val':10})

        self.results_frame = results_frame
        
        self.key_dict={'given_data': 'listed system parameters', #to override in the particular step
                       'nr_rys': 'NR', #to override in the particular steps
                       'param_name':vlatex(self.analysis_key),

            
                      }
        
        self.step_key_dict= {}
        
        self.preview=preview
        
        self.chapter_name=type(self).chapter_name
        self.summary_sec_name=type(self).summary_sec_name
        self.optimization_desc_bank=type(self). optimization_desc_bank
        self.units_dict=type(self).units_dict

    def content_preview(self,content=','):
        
        if self.preview:
            print('\n',content,'\n')
    
    #tmp_dict = {key:row[1].to_dict()[key]  for key in  model.system_parameters()  }
        
    def add_described_result(self,numerical_result,model,marker,format_dict={},caption=None):
        
        doc=self.doc

        ref_name=str(marker)
        fig_name='./plots/'+ref_name

        column_names=numerical_result.columns



        step_key_max_dict={str(name)+'_max':numerical_result[name].round(1).max() for name   in column_names}
        step_key_min_dict={str(name)+'_min':numerical_result[name].round(1).min() for name   in column_names}
        
        self.step_key_dict['nr_rys']=Ref(marker).dumps()
        self.step_key_dict={**self.step_key_dict,**step_key_min_dict,**step_key_max_dict,**format_dict}

        

        random_intro=random.choice(self.intro_bank)
        self.content_preview(random_intro)

        doc.append(
            NoEscape(
                 random_intro.format(
                     **{**self.key_dict,**self.step_key_dict})+'\n'*2
            ))


        with doc.create(Figure(position='H')) as fig:


            #                     var_filename='./tikz_plots/graph'+str(no)
            #                     #solution_tmp.columns=[vlatex(name) for name in solution_tmp.columns]
            #                     solution_tmp.to_tikz_plot(var_filename)

            #                     fig.append(Command('input',arguments=NoEscape(var_filename)))


            ax = numerical_result.plot(subplots=True, figsize=(10, 7))
            #ax=solution_tmp.plot(subplots=True)
            ([
                ax_tmp.legend(['$' + vlatex(sym) + '$'])
                for ax_tmp, sym in zip(ax, numerical_result.columns)
            ])

            plt.savefig(fig_name+'.png')
            fig.add_image(fig_name,width=NoEscape('15cm'))
            
            if self.preview:
                plt.show()
                
            plt.close()
            
            
            fig.add_caption(
                NoEscape(self.caption.format(**{**self.key_dict,**self.step_key_dict}))  )
            fig.append(Label(marker))
            


        random_summary=random.choice(self.summary_bank)
        self.content_preview(random_summary)

        doc.append(
            NoEscape(
                random_summary.format(
                     **{**self.key_dict,**self.step_key_dict})+'\n'*2
                ))

        doc.append(NewPage())
        
        return doc

    def influeance_level(self,data):
        

        indicator=round(2*float(
                            abs((data.abs().max()-data.abs().min())/(0.001+data.abs().max()+data.abs().min()))
                            ),1)/2
            
        select=self.sensitivity_level

        a=np.asarray(list(select.keys()))

        print('indicator: ',indicator)
        print(select[a.flat[np.abs(a - indicator).argmin()]])
        print(select[a.flat[np.abs(a - indicator).argmin()]])

        return select[a.flat[np.abs(a - indicator).argmin()]]
    
    def cost_val(self, result_frame,key='cost_val'):
        
        cost_val=result_frame[key]
        
        return cost_val
        
    def doit(self):

        doc = self.doc
        simulation_results_frame = self.results_frame
        
        generalized_coords=(next((simulation_results_frame['simulations']).items())[1].columns)
        
        description=next((simulation_results_frame['description']).items())[1]
        
        summary_frame_max = pd.DataFrame(
                data=[
                    a.round(1).max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=simulation_results_frame[self.analysis_key]).round(2)

        summary_frame_min = pd.DataFrame(
                data=[
                    a.round(1).min()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=simulation_results_frame[self.analysis_key]).round(2)
        
        
        strings_dict={'param_name':vlatex(self.analysis_key) ,'max_val':10,'min_val':10,'q_name':vlatex(self.analysis_key),'t_val':100}
        
        strings_dict={**strings_dict,**{str(coord):vlatex(coord)   for coord in generalized_coords}}
        
        
        doc = self.doc
        simulation_results_frame = self.results_frame
        doc.append(
            Chapter(NoEscape(
                self.chapter_name.format(
                    **strings_dict
                    )
            )))

        with doc.create(Section('Wyniki symulacji')):

            simulation_results_frame = self.results_frame

            random_sec_intro=random.choice(self.sec_intro_bank)
            print('\n',random_sec_intro,'\n')
            
            
            
            doc.append(
                NoEscape(
                     random_sec_intro.format(
                         **{**self.key_dict})+'\n'*2
                ))
            
            symbols_df=pd.DataFrame(data=[
                        {'Zmienna':'$ {eq}  $'.format(eq=vlatex(key)),
                         'Liczba przypadków':len(simulation_results_frame[key].unique()),
                         'Wartość minimalna':simulation_results_frame[key].min(),
                         'Wartość maksymalna':simulation_results_frame[key].max(),
                         'Jednostka': '$ {unit}  $'.format(unit=(self.units_dict[key]).dumps() ),
                        } for key in  simulation_results_frame.columns if isinstance(key,Symbol)  ]).round(2)
                         
            doc.append(NoEscape( symbols_df.to_latex(index=False,escape=False))  )
            
            doc.append(NewPage())
            
            for no, row in enumerate(simulation_results_frame.iterrows()):
                    

        
                ref_name=self.analysis_name+'_'+ str(self.analysis_key)+'_'+ str(no)
                mrk = Marker(ref_name, prefix='fig')
            
                num_data=row[1][self.sims_name]
                model=row[1]['model']
                
                coordinates_names_dict={str(coord):vlatex(coord)   for coord in model.Y}
                
                step_vals_dict = {key:row[1].to_dict()[key]  for key in  model.system_parameters()  }
                given_data_str='$' +'$, $'.join([vlatex(Eq(lhs, np.round(rhs,1),evaluate=False))+self.units_dict[lhs].dumps() for lhs, rhs in step_vals_dict.items()]) + '$.'
                
                self.step_key_dict={**self.step_key_dict,'given_data':given_data_str,**step_vals_dict,**coordinates_names_dict}
                
                
                self.add_described_result(numerical_result=num_data,model=model,marker=mrk,format_dict={},caption=None)
                
                doc.append(NewPage())
            

        with doc.create(
                Section(NoEscape(
                    self.summary_sec_name.format(**{**self.key_dict,**self.step_key_dict}))
                 )):

            summary_frame = pd.DataFrame(
                data=[
                    a.max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=[str(elem) for elem in simulation_results_frame[self.analysis_key].round(2)]
                )
            
            column_names=summary_frame.columns
            

            step_key_max_dict={str(name)+'_max':summary_frame[name].round(1).abs().max() for name   in column_names}
            step_key_min_dict={str(name)+'_min':summary_frame[name].round(1).abs().min() for name   in column_names}
            step_key_idmax_dict={str(name)+'_idmax':summary_frame[name].round(1).abs().idxmax() for name   in column_names}
            step_key_idmin_dict={str(name)+'_idmin':summary_frame[name].round(1).abs().idxmin() for name   in column_names}
            
#             def i
#                         indicator=round(2*float(
#                             ((summary_frame[name].abs().max()-summary_frame[name].abs().min())/(summary_frame[name].abs().max()+summary_frame[name].abs().min()).abs())
#                             ),1)/2
            
#             select=self.sensitivity_level

#             a=np.asarray(list(select.keys()))

#             print(select[a.flat[np.abs(a - indicator).argmin()]])
#             print(select[a.flat[np.abs(a - indicator).argmin()]])
#             print(select[a.flat[np.abs(a - indicator).argmin()]])
            
            step_key_influence_dict={str(name)+'_inf':self.influeance_level(summary_frame[name]) for name   in column_names}
            
            print(step_key_influence_dict)
            self.step_key_dict={**self.step_key_dict,**step_key_max_dict,**step_key_min_dict,**step_key_idmax_dict,**step_key_idmin_dict,**step_key_influence_dict}

            step_key_max_dict={str(name)+'_max':summary_frame[name].round(1).max() for name   in column_names}
            step_key_min_dict={str(name)+'_min':summary_frame[name].round(1).min() for name   in column_names}
            step_key_idmax_dict={str(name)+'_idmax':summary_frame[name].round(1).idxmax() for name   in column_names}
            step_key_idmin_dict={str(name)+'_idmin':summary_frame[name].round(1).idxmin() for name   in column_names}
            
            self.step_key_dict={**self.step_key_dict,**step_key_max_dict,**step_key_min_dict,**step_key_idmax_dict,**step_key_idmin_dict}

            
            doc.append(
                NoEscape(
                'Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty symulacji numerycznych. Dla uzyskanych danych symulacyjnych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego układu. '
                'Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\). '
                .format(**{**self.key_dict,**self.step_key_dict}))
            )
            
            ref_summary_name=str(self.analysis_key)+'.png'
            summary_fig_name='./plots/summary_'+ref_summary_name
            mrk_summary = Marker(ref_summary_name, prefix='fig')
           
            with doc.create(Figure(position='H')) as fig:

                ax = summary_frame.plot(subplots=True, figsize=(10, 7))
                #ax=solution_tmp.plot(subplots=True)
                ([
                    ax_tmp.legend(['$' + vlatex(sym) + '$'])
                    for ax_tmp, sym in zip(ax, summary_frame.columns)
                ])

                
                
                
                plt.savefig(summary_fig_name)
                fig.add_image(summary_fig_name,width=NoEscape('15cm'))
                plt.show()
                fig.add_caption(
                    NoEscape(
                        'Zestawienie zmienności parametru \({param_name}\)'
                        .format(**{**self.key_dict,**self.step_key_dict})))
                fig.append(Label(mrk_summary))


            


            
            
            self.step_key_dict={**self.step_key_dict,**step_key_max_dict,**step_key_min_dict}
            
            random_conclusions=random.choice(self.conclusion_bank)
            print('\n',random_conclusions,'\n')
            
            summary_description=next((simulation_results_frame['description']).items())[1]
            #print({**self.key_dict,**self.step_key_dict})
            doc.append(NoEscape(
                random_conclusions.format(**{**self.key_dict,**self.step_key_dict})
                )
            )
        return doc
    
class DataTable(Table):
    _latex_name='table'
    def __init__(self,numerical_data,position=None):
        super().__init__(position=position)
        print(numerical_data)
        self._numerical_data = numerical_data
        self.position = position
        
    def add_table(self,numerical_data=None):
        
        if numerical_data!=None:
            self._numerical_data=numerical_data
        
        tab =  self._numerical_data
        self.append(NoEscape(tab.to_latex(index=False,escape=False)))