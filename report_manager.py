import numpy as np
import pandas as pd

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Alignat, Package, Quantity, Command, Label, Ref, Marker, NewPage
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


class ReportModelAnalysis:
    def __init__(self, model=None):
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
        doc.packages.append(Package('preamble'))

        self.doc = doc

        self.model = model

    def doit(self):

        doc = self.doc

        model = self.model

        with doc.create(Chapter('Weryfikacja z danymi doświadczalnymi')):
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
    
class ReportSimulationsAnalysis:
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
        
        
        
        
       
        
        self.sec_introduction_bank_1=[
            'Wykorzystując przygotowany model dynamiczny, przeprowadzono serię symulacji mających na celu ocenę wrażliwości modelu na zmiany wybranych parametrów. ',
            'Przeprowadzenie analizy numerycznej wymaga znajomości modelu dynamicznego oraz zdeterminowania parametrów oddziałujących na zachowanie układu wskutek ich zmian. ',
            'Wykorzystując przygotowane środowisko obliczeniowe, wykonano szereg symulacji, opracowano dane numeryczne i przedstawiono je w postaci następujących wykresów. ',
            'Celem zbadania i oceny dynamiki rozpatrywanego systemu przeprowadzono symulacje numeryczne przedstawijące przebiegi zmiennych układu w funkcji czasu. '
            'W oparciu o wybrany model dynamiczny wózka, wykonano analizy numeryczne zaprezentowane na kolejnych wykresach. '
            ]
        
        self.sec_introduction_bank_2=[
            'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
            'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
            'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
            'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
        ]
        
        
        self.sec_introduction_bank_3=[
            'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
            'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
            'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
            'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
        ]        
        
        self.sec_intro_bank=[str(RandomDescription(self.sec_introduction_bank_1,self.sec_introduction_bank_2,self.sec_introduction_bank_3)) for obj in range(30)  ]        
        
        self.introduction_bank_1=[
            'Na wykresie {nr_rys} przedstawiano zmiany wielkości dynamicznych charakteryzujących ruch wózka. Po przeanalizowaniu można zanotować wzajemną zależność poszczególnych wielkości dynamicznych.',
            'Charakter przebiegu wartości dynamicznych wózka został przedstawiony na rysunku {nr_rys}.',
            'Wykres {nr_rys} pokazuje zmienność parametrów dynamicznych wózka w trakcie symulownego przejazdu.',
            'Rysunek {nr_rys} prezentuje wykres opisujący zmienność parametrów dynamicznych wózka w trakcie jazdy po zamodelowanym torze.',
        ]
        
        self.introduction_bank_2=[
            'Symulacja została przeprowadzona dla następujących danych: {given_data}.',
            'Zaprezentowane wyniki numeryczne otrzymano dla danych równych:{given_data}.',
            'Eksperyment numeryczny został przeprowadzona dla następujących danych:{given_data}.',
            'Wyniki stmulacji numerycznych otrzymano dla:{given_data}.',
        ]
        

        
        self.introduction_bank_3=[
            'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu stwierdza się oscylacyjny charakter odpowiedzi układu. Dynamika systemu ma stabilny charakter',
            'Charakter zmian przebiegów poszczególnch współrzędnych jest stabilny. Otrzymane sygnały mieszczą się w zakresie dopuszczalnych uwarunkowaniami fizycznymi.',
            'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu mają drganiowy charakter i są zgodne z oczekiwaniami.',
            'Bazując na otrzymanych wynikach można stwierdzić wzajemną zależność poszczególnych wielkości dynamicznych. Dodatkowo wielkości modelu stabilizują się w czasie. '
        ]
        
        
        self.intro_bank=[str(RandomDescription(self.introduction_bank_1,self.introduction_bank_2,self.introduction_bank_3)) for obj in range(30)  ]
        
        
        
        self.time_domain_summary_1 = [
            'Przedstawione wyniki pokazują wpływ zmian badanego parametru \({param_name}\) na dynamikę wózka. ',
            'Opracowane wykresy porównawcze pozwoliły na oszacowanie wpływu zmiennej \({param_name}\) na dynamikę rozpatrywanego układu. ',
            'Zaprezentowane wyniki posłużyły ocenie wrażliwości modelowanego obiektu na badaną zmienną \({param_name}\). ',
            'Przygotowane przebiegi służą zaprezentowaniu wyników symulacji numerycznych badanej zmiennej w funkcji wybranego parametru wpływowego \({param_name}\). ',
            'Analizowany model wykazuje wrażliwość na zmianę rozważanego parametru \({param_name}\). ',
        ]
        
        self.time_domain_summary_2 = [
            'Zaobserwowano wpływ badanego parametru \({x(t)}\) na drgania wózka, których maksymalna wartość amplitudy wyniosła {x(t)_max}. ',
            'Analizowany przebieg parametru \({x(t)}\) osiąga wartość maksymalną równą {x(t)_max}. '
            'Minimalna wartość dla przebiegu \({x(t)}\) jest równa {min_val}. '
            'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({x(t)}\) nie występują amplitudy o wartości mniejszej niż {x(t)_min} oraz większej niż {x(t)_max}. '
            'Dokonując analizy wyników symulacji dla parametru \({x(t)}\) stwierdza się, że amplitudy drgań nie przekraczają wartości {x(t)_max}, a wartość mininalna wynosi {x(t)_min}. ',
        ]
        
        
        self.time_domain_summary_3 = [ 
            'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi(t)}\) wynosi {varphi(t)_max}. ',
            'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi(t)_max} dla współrzędnej \({varphi(t)}\). '
            'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {varphi(t)_min}. ',
            'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
            'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
        ]
        

        self.time_domain_summary_4 = [ 
            'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({z(t)}\) wynosi {z(t)_max}. ',
            'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {z(t)_max} dla współrzędnej \({z(t)}\). '
            'Wartość minimalna przebiegu współrzędnej \({z(t)}\) jest równa {z(t)_min}. ',
            'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({z(t)}\) nie przekracza {z(t)_max} i jest ograniczona z dołu przez {z(t)_min}. ',
            'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({z(t)}\) nie przekraczaa wartości {z(t)_max}, a wartość mininalna wynosi {z(t)_min}. ',
        ]        
        
        
        self.time_domain_summary_5 = [ 
            'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi_RC(t)}\) wynosi {varphi_RC(t)_max}. ',
            'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi_RC(t)_max} dla współrzędnej \({varphi_RC(t)}\).'
            'Wartość minimalna przebiegu współrzędnej \({varphi_RC(t)}\) jest równa {varphi_RC(t)_min}. ',
            'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi_RC(t)}\) nie przekracza {varphi_RC(t)_max} i jest ograniczona z dołu przez {varphi_RC(t)_min}. ',
            'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi_RC(t)}\) nie przekraczaa wartości {varphi_RC(t)_max}, a wartość mininalna wynosi {varphi_RC(t)_min}. ',
        ]        
        
        
        self.summary_bank=[str(RandomDescription(
                    self.time_domain_summary_1,
                    self.time_domain_summary_2,
                    self.time_domain_summary_3,
                    self.time_domain_summary_4,
                    self.time_domain_summary_5,
                        )) for obj in range(30)  ]
        
        self.time_domain_summary_7 = [
            'Zaobserwowane odczyty wynikające z wpływu \({param_name}\) występują w konketrnym czasie ruchu {t_val}. ',
            'Odczytana maksymalna wartość amplitudy {max_val} wynika z działania \({param_name}\) wystąpiła konkretnie w czasie {t_val}. '
            'Minimalny zaobserwowany wynik amplitudy {min_val} miał miejsce w {t_val} czasu trwania pomiaru. ',
        ]

        self.time_domain_summary_8 = [ 
            'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionych są następujące: {max_val_list}. ',
            'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się następującymi wartościami maksymalnymi: {max_val_list}. '
            'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {min_val}. ',
            'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
            'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
        ]
        
        
       
        self.freq_domain_summary =[
           'Ponadto sprawdzono własności spektralne otrzymanych wyników symulacyjnych.',
           'Według przeprowadzonych badań modelowych struktura częstotliwościowa nie zawiera komonentów większych niż {max_val_spec}.',
           'Amplitudy drgań badanego parametru \({q_name}\) przyjmują maksymalną wartość równą {max_val_spec}.',
           'Na podstawie analizy widmowej minimalna wartość częstotliwości rozpatrywanego parametru wynosi {min_val_spec}.',
            
                                   
        ]
        # 'Wartość największa wynosi: {max_val}'.format(**{'max_val':10})

        self.results_frame = results_frame

    def doit(self):

        doc = self.doc
        simulation_results_frame = self.results_frame
        
        generalized_coords=(next((simulation_results_frame['simulations']).items())[1].columns)
        
        summary_frame_max = pd.DataFrame(
                data=[
                    a.max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=simulation_results_frame[self.analysis_key]).round(2)

        summary_frame_min = pd.DataFrame(
                data=[
                    a.min()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=simulation_results_frame[self.analysis_key]).round(2)
        
        
        strings_dict={'param_name':vlatex(self.analysis_key) ,'max_val':10,'min_val':10,'q_name':vlatex(self.analysis_key),'t_val':100}
        
        strings_dict={**strings_dict,**{str(coord):vlatex(coord)   for coord in generalized_coords}}
        
        
        doc = self.doc
        simulation_results_frame = self.results_frame
        doc.append(
            Chapter(NoEscape(
                'Analiza wrażliwości dla parametru \({param_name}\)'.format(
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
                         **strings_dict)+'\n'*2
                ))
            
            symbols_df=pd.DataFrame(data=[
                        {'Zmienna':'$ {eq}  $'.format(eq=vlatex(key)),
                         'Liczba przypadków':len(simulation_results_frame[key].unique()),
                         'Wartość minimalna':simulation_results_frame[key].min(),
                         'Wartość maksymalna':simulation_results_frame[key].max() } for key in  simulation_results_frame.columns if isinstance(key,Symbol)  ]).round(2)
            doc.append(NoEscape( symbols_df.to_latex(index=False,escape=False))  )
            
            doc.append(NewPage())
            
            for no, row in enumerate(simulation_results_frame.iterrows()):
                
                model=row[1]['model']
                
                ref_name=self.analysis_name+'_'+ str(self.analysis_key)+'_'+ str(no)
                fig_name='./plots/'+ref_name
                column_names=(next((simulation_results_frame['simulations']).items())[1].columns)
                
                mrk = Marker(ref_name, prefix='fig')
                
                strings_dict_max={str(name)+'_max':summary_frame_max[name][row[1][self.analysis_key]] for name   in column_names}
                strings_dict_min={str(name)+'_min':summary_frame_min[name][row[1][self.analysis_key]] for name   in column_names}
                strings_dict_max['nr_rys']=Ref(mrk).dumps()
             
                strings_dict_max['max_val_list']=', '.join([ '\( {coord} \) - \( {coord_max} \)'.format(**{'coord':vlatex(coord),'coord_max':str(summary_frame_max[coord][row[1][self.analysis_key]])})   for coord in column_names])
            
            
                
            
                strings_dict_max={**strings_dict_max,**strings_dict_min}
            
            
                
            
            
                #print([ '\( {coord} \) - \( {coord_max} \)'.format(**{'coord':str(coord),'coord_max':str(summary_frame_max[coord][row[1][self.analysis_key]])})   for coord in column_names])
                
            
                #display(strings_dict_max)
                #display(simulation_results_frame)

                solution_tmp = row[1]['simulations']

                tmp_dict = {key:row[1].to_dict()[key]  for key in  model.system_parameters()  }
                strings_dict_max['given_data']='$' +'$, $'.join([vlatex(Eq(lhs, rhs))for lhs, rhs in tmp_dict.items()]) + '$.'

#                 tmp_dict.pop('simulations')
#                 model=tmp_dict.pop('model')
#                 tmp_dict.pop('num_model')

                n = random.randint(0, 99)

                #mrk = Marker(self.analysis_name + str(no), prefix='fig')

        
                random_intro=random.choice(self.intro_bank)
                print('\n',random_intro,'\n')
        
                doc.append(
                    NoEscape(
                         random_intro.format(
                             **strings_dict,**strings_dict_max)+'\n'*2
                    ))
                
            
                with doc.create(Figure(position='H')) as fig:
                    #                 kitten_pic.add_image(image_filename, width='120px')
                    #                 kitten_pic.add_caption('Look it\'s on its back')

                    #solution_tmp.columns=[vlatex(sym) for sym in (solution_tmp.columns)]

                    #                     var_filename='./tikz_plots/graph'+str(no)
                    #                     #solution_tmp.columns=[vlatex(name) for name in solution_tmp.columns]
                    #                     solution_tmp.to_tikz_plot(var_filename)

                    #                     fig.append(Command('input',arguments=NoEscape(var_filename)))

                    #$'+vlatex(coord) + '$ dla parametru $'+vlatex(param_label)+'$'))

                    ax = solution_tmp.plot(subplots=True, figsize=(10, 7))
                    #ax=solution_tmp.plot(subplots=True)
                    ([
                        ax_tmp.legend(['$' + vlatex(sym) + '$'])
                        for ax_tmp, sym in zip(ax, solution_tmp.columns)
                    ])
                    #plt.figure(figsize=(15,10))
                    #fig.add_plot(width=NoEscape('0.9\\textwidth'))
                    
                    plt.savefig(fig_name+'.png')
                    fig.add_image(fig_name,width=NoEscape('20cm'))
                    plt.show()
                    fig.add_caption(
                        NoEscape('Przebiegi czasowe modelu dla danych: {given_data}'.format(**strings_dict,**strings_dict_max))  )
                    fig.append(Label(mrk))
                    #fig.add_caption(NoEscape('Przebiegi czasowe zmiennej zależnej')) #$'+vlatex(coord) + '$ dla parametru $'+vlatex(param_label)+'$'))

                n = random.randint(0, 99)
                
                random_summary=random.choice(self.summary_bank)
                print('\n',random_summary,'\n')
                
                doc.append(
                    NoEscape(
                        random_summary.format(
                             **strings_dict,**strings_dict_max)+'\n'*2
                        ))
                
                doc.append(NewPage())
            

        with doc.create(
                Section(NoEscape(
                    'Spostrzeżenia i wnioski z analizy wpływu parametru  \({param_name}\)'
                    .format(**strings_dict,**strings_dict_max))
                 )):

            summary_frame = pd.DataFrame(
                data=[
                    a.max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=simulation_results_frame[self.analysis_key])
            doc.append(
                NoEscape(
                'Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie wartości maksymalnych w funkcji parametru \({param_name}\)'
                .format(**strings_dict))
            )
            
            ref_summary_name=str(self.analysis_key)+'.png'
            summary_fig_name='./plots/summary_'+ref_summary_name
            mrk_summary = Marker(ref_summary_name, prefix='fig')
           
            with doc.create(Figure(position='H')) as fig:

                ax = summary_frame.plot(subplots=True, figsize=(10, 7))
                #ax=solution_tmp.plot(subplots=True)
                ([
                    ax_tmp.legend(['$' + vlatex(sym) + '$'])
                    for ax_tmp, sym in zip(ax, solution_tmp.columns)
                ])

                
                
                
                plt.savefig(summary_fig_name)
                fig.add_image(summary_fig_name,width=NoEscape('20cm'))
                plt.show()
                fig.add_caption(
                    NoEscape(
                        'Zestawienie zmienności parametru \({param_name}\)'
                        .format(**strings_dict)))
                fig.append(Label(mrk_summary))

            doc.append(NoEscape(
                'Przedstawione zestawienie wskazuje, że parametr \({param_name}\) ...'
                .format(**strings_dict))
            )
        return doc
