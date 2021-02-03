import numpy as np
import pandas as pd

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Alignat, Package, Quantity, Command, Label, Ref, Marker, NewPage, NewLine, Eqref, Table
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
import pint

import datetime as dtime


class Equation(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package('mathtools')]
    escape = False
    content_separator = "\n"


class DMath(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package('breqn'), Package('flexisym')]
    escape = False
    content_separator = "\n"


# class EqRef(Environment):
#         """A class to wrap LaTeX's alltt environment."""

#     packages = [Package('breqn'), Package('flexisym')]
#     escape = False
#     content_separator = " "


class ReportModelAnalysis:

    chapter_name = 'Model analityczny układ o swobody'

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
        #         doc.packages.append(Package('preamble'))
        self.doc = doc

        self.model = model

        self.chapter_name = type(self).chapter_name

    def lagranges(self):

        doc = self.doc

        model = self.model

        #musisz też dopisać o jednym stopniu - to będzie ten indeks
        with doc.create(Chapter(self.chapter_name)):

            with doc.create(DMath()):
                doc.append('T=')

                doc.append(
                    vlatex(model.lagrangian() +
                           model.lagrangian().subs({vel: 0
                                                    for vel in model.u})))
            with doc.create(DMath()):
                doc.append('V=')

                doc.append(
                    vlatex(-self.model.lagrangian().subs(
                        {vel: 0
                         for vel in model.u})))

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
                    doc.append(
                        vlatex(self.model.equilibrium_equation().doit()[i]))
            for lhs, rhs in zip(
                    symbols('k_1:{ndof}_1:{ndof}'.format(ndof=len(model.q) +
                                                         1)),
                (list(model.linearized().stiffness_matrix()))):
                with doc.create(DMath()):
                    doc.append(vlatex([Eq(lhs, rhs)]))
            for i in range(len(model.q)):
                with doc.create(DMath()):
                    doc.append(vlatex(self.model._eoms[i]))
            for lhs, rhs in zip(
                    symbols('m_1:{ndof}_1:{ndof}'.format(ndof=len(model.q) +
                                                         1)),
                (list(model.inertia_matrix()))):
                with doc.create(DMath()):
                    doc.append(vlatex([Eq(lhs, rhs)]))
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
    def __init__(self, *args, random_or_order=None, **kwargs):

        self.pool_list = args
        self.random_or_order = random_or_order

    def to_string(self):

        if hasattr(self.random_or_order, '__contains__'):
            self.order = list(random_or_order)
        else:
            self.selected_sentences = [
                random.choice(pool) for pool in self.pool_list
            ]

        return ' '.join(self.selected_sentences)

    def __str__(self):

        return self.to_string()

    def __repr__(self):
        return self.__class__.__name__ + ': ' + self.to_string()


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

freq_domain_summary = [
    'Ponadto sprawdzono własności spektralne otrzymanych wyników symulacyjnych.',
    'Według przeprowadzonych badań modelowych struktura częstotliwościowa nie zawiera komonentów większych niż {max_val_spec}.',
    'Amplitudy drgań badanego parametru \({q_name}\) przyjmują maksymalną wartość równą {max_val_spec}.',
    'Na podstawie analizy widmowej minimalna wartość częstotliwości rozpatrywanego parametru wynosi {min_val_spec}.',
]

summary_bank_composed_model_gen = RandomDescription(
    time_domain_summary_1,
    time_domain_summary_2,
    time_domain_summary_3,
    time_domain_summary_4,
    time_domain_summary_5,
)


summary_bank_chair_model_gen = RandomDescription(
    time_domain_summary_1,
    time_domain_summary_2,
    time_domain_summary_3,
    time_domain_summary_4,
#    time_domain_summary_5,
)


summary_bank_chair_model_gen = RandomDescription(
    time_domain_summary_1,
    time_domain_summary_2,
    time_domain_summary_3,
    time_domain_summary_4,
#    time_domain_summary_5,
)


summary_bank_drive_model_gen = RandomDescription(
    time_domain_summary_1,
#     time_domain_summary_2,
#     time_domain_summary_3,
#     time_domain_summary_4,
    time_domain_summary_5,
)




summary_bank_chair_move_model_gen = RandomDescription(
    time_domain_summary_1,
    time_domain_summary_2,
    time_domain_summary_3,
#     time_domain_summary_4,
 #   time_domain_summary_5,
)

summary_bank_drive_model = [
    str(
        RandomDescription(
            time_domain_summary_1,
#             time_domain_summary_2,
#             time_domain_summary_3,
#             time_domain_summary_4,
            time_domain_summary_5,
        )) for obj in range(30)
]

measurement_summary_1 = [
    'Przedstawione wyniki pokazują wpływ zmian badanego parametru \({param_name}\) na dynamikę wózka określonego z różnic porcentowych przyśpieszeń akcelerometru. ',
    'Opracowane wykresy porównawcze pozwoliły na określenie wpływu zmiennej \({param_name}\) na dynamikę rozpatrywanego układu poprzez porówanie wykresów z wykonanych pomiarów. ',
    'Zaprezentowane wyniki posłużyły ocenie wrażliwości modelowanego obiektu na badaną zmienną \({param_name}\). ',
    'Przygotowane przebiegi służą zaprezentowaniu wyników pomiarów fizycznych z wykorzystaniem akcelerometru badanej w funkcji wybranego parametru wpływowego \({param_name}\). ',
    'Analizowany model wykazuje wrażliwość na zmianę rozważanego parametru \({param_name}\). ',
]

measurement_summary_2 = [
    'Zaobserwowano wpływ badanego parametru \({a_ox}\) na drgania wózka, których maksymalna wartość amplitudy procentu (%) przyśpieszenia ziemskiego wyniosła {a_ox_max}. ',
    'Analizowany przebieg parametru \({a_ox}\) osiąga wartość maksymalną równą {a_ox_max}. ',
    'Minimalna wartość dla przebiegu \({a_ox}\) jest równa {a_ox_min}. ',
    'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_ox}\) nie występują amplitudy o wartości mniejszej niż {a_ox_min} oraz większej niż {a_ox_max}. ',
    'Dokonując analizy wyników symulacji dla parametru \({a_ox}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_ox_max}, a wartość mininalna wynosi {a_ox_min}.',
    'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_ox}\) nie przekracza {a_ox_idmax} i jest ograniczona z dołu przez {a_ox_idmin}. ',
    'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_ox}\) nie przekraczaa wartości {a_ox_idmax}, a wartość mininalna wynosi {a_ox_idmin}. ',
]

measurement_summary_3 = [
    'Zaobserwowano wpływ badanego parametru \({a_oz}\) na drgania wózka, których maksymalna wartość amplitudy porecntu (%) przyśpieszenia ziemskiego wyniosła {a_oz_max}. ',
    'Analizowany przebieg parametru \({a_oz}\) osiąga wartość maksymalną równą {a_oz_max}. ',
    'Minimalna wartość dla przebiegu \({a_oz}\) jest równa {a_oz_min}. ',
    'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_oz}\) nie występują amplitudy o wartości mniejszej niż {a_oz_min} oraz większej niż {a_oz_max}. ',
    'Dokonując analizy wyników symulacji dla parametru \({a_oz}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_oz_max}, a wartość mininalna wynosi {a_oz_min}.',
    'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_oz}\) nie przekracza {a_oz_idmax} i jest ograniczona z dołu przez {a_oz_idmin}. ',
    'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_oz}\) nie przekraczaa wartości {a_oz_idmax}, a wartość mininalna wynosi {a_oz_idmin}. ',
]

measurement_summary_4 = [
    'Zaobserwowano wpływ badanego parametru \({a_rx}\) na drgania wózka, których maksymalna wartość amplitudy porecntu (%) przyśpieszenia ziemskiego wyniosła {a_rx_max}. ',
    'Analizowany przebieg parametru \({a_rx}\) osiąga wartość maksymalną równą {a_rx_max}. ',
    'Minimalna wartość dla przebiegu \({a_rx}\) jest równa {a_rx_min}. ',
    'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_rx}\) nie występują amplitudy o wartości mniejszej niż {a_rx_min} oraz większej niż {a_rx_max}. ',
    'Dokonując analizy wyników symulacji dla parametru \({a_rx}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_rx_max}, a wartość mininalna wynosi {a_rx_min}.',
    'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_rx}\) nie przekracza {a_rx_idmax} i jest ograniczona z dołu przez {a_rx_idmin}. ',
    'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_rx}\) nie przekraczaa wartości {a_rx_idmax}, a wartość mininalna wynosi {a_rx_idmin}. ',
]

measurement_summary_5 = [
    'Zaobserwowano wpływ badanego parametru \({a_rz}\) na drgania wózka, których maksymalna wartość przyśpieszenia wyniosła {a_rz_max}. ',
    'Analizowany przebieg parametru \({a_rz}\) osiąga wartość maksymalną równą {a_rz_max}. ',
    'Minimalna wartość dla przebiegu \({a_rz}\) jest równa {a_ox_min}. ',
    'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_rz}\) nie występują amplitudy o wartości mniejszej niż {a_rz_min} oraz większej niż {a_rz_max}. ',
    'Dokonując analizy wyników symulacji dla parametru \({a_rz}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_rz_max}, a wartość mininalna wynosi {a_rz_min}.',
    'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_rz}\) nie przekracza {a_rz_idmax} i jest ograniczona z dołu przez {a_rz_idmin}. ',
    'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_rz}\) nie przekraczaa wartości {a_rz_idmax}, a wartość mininalna wynosi {a_rz_idmin}. ',
]

summary_bank_composed_model_measurements = [
    str(
        RandomDescription(
            measurement_summary_1,
            measurement_summary_2,
            measurement_summary_3,
            measurement_summary_4,
            measurement_summary_5,
        )) for obj in range(30)
]

sec_introduction_bank_0={'Zaprezentowane symulacje wykonano dla następujących danych: {given_data}.','Do wykonania symulacji użyto danych jak następuje: {given_data}','Zakres danych użytych do procesu symulacji jest następujący: {given_data}.','Do opracowania wyników symulacyjnych posłużono się niniejszymi danymi: {given_data}.'}
sec_introduction_bank_1 = [
    'Wykorzystując przygotowany model dynamiczny, przeprowadzono serię symulacji mających na celu ocenę wrażliwości modelu na zmiany wybranych parametrów. ',
    'Przeprowadzenie analizy numerycznej wymaga znajomości modelu dynamicznego oraz zdeterminowania parametrów oddziałujących na zachowanie układu wskutek ich zmian. ',
    'Wykorzystując przygotowane środowisko obliczeniowe, wykonano szereg symulacji, opracowano dane numeryczne i przedstawiono je w postaci następujących wykresów. ',
    'Celem zbadania i oceny dynamiki rozpatrywanego systemu przeprowadzono symulacje numeryczne przedstawijące przebiegi zmiennych układu w funkcji czasu. '
    'W oparciu o wybrany model dynamiczny obiektu, wykonano analizy numeryczne zaprezentowane na kolejnych wykresach. '
]

sec_introduction_bank_2 = [
    'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
    'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
    'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
    'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
]

sec_introduction_bank_3 = [
    'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
    'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
    'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
    'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
]

sec_intro_composed_model_gen = RandomDescription(sec_introduction_bank_1,
                                                 sec_introduction_bank_2,
                                                 sec_introduction_bank_3,
                                                sec_introduction_bank_0)

sec_intro_composed_model = [
    str(
        RandomDescription(sec_introduction_bank_1, sec_introduction_bank_2,
                          sec_introduction_bank_3)) for obj in range(30)
]

introduction_bank_meas_1 = [  #poprawione
    'Na wykresie {nr_rys} przedstawiano zmiany wielkości dynamicznych charakteryzujących ruch wózka. Po przeanalizowaniu można zanotować wzajemną zależność poszczególnych wielkości dynamicznych.',
    'Charakter przebiegu wartości dynamicznych wózka został przedstawiony na rysunku {nr_rys}.',
    'Wykres {nr_rys} pokazuje zmienność parametrów dynamicznych wózka w trakcie przejazdu pomiarowego.',
    'Rysunek {nr_rys} prezentuje wykres opisujący zmienność parametrów dynamicznych wózka w trakcie jazdy pomiarowej po torze.',
]

introduction_bank_meas_2 = [  #poprawione
    'Pomiar został przeprowadzony dla następujących danych: {given_data}.',
    'Zaprezentowane wyniki pomiarowe otrzymano dla danych równych:{given_data}.',
    'Badanie pomiarowe zostało przeprowadzone dla następujących danych:{given_data}.',
    'Wyniki pomiaru otrzymano dla:{given_data}.',
]

introduction_bank_meas_3 = [  #poprawione (nie trzeba zmieniać zdań)
    'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu stwierdza się oscylacyjny charakter odpowiedzi układu. Dynamika systemu ma stabilny charakter',
    'Charakter zmian przebiegów poszczególnch współrzędnych jest stabilny. Otrzymane sygnały mieszczą się w zakresie dopuszczalnych uwarunkowaniami fizycznymi.',
    'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu mają drganiowy charakter i są zgodne z oczekiwaniami.',
    'Bazując na otrzymanych wynikach można stwierdzić wzajemną zależność poszczególnych wielkości dynamicznych. Dodatkowo wielkości modelu stabilizują się w czasie. '
]

introduction_bank_1 = [
    'Na wykresie {nr_rys} przedstawiano zmiany wielkości dynamicznych charakteryzujących ruch obiektu. Po przeanalizowaniu można zanotować wzajemną zależność poszczególnych wielkości dynamicznych.',
    'Charakter przebiegu wartości dynamicznych układu został przedstawiony na rysunku {nr_rys}.',
    'Wykres {nr_rys} pokazuje zmienność parametrów dynamicznych wózka w trakcie symulownego przejazdu.',
    'Rysunek {nr_rys} prezentuje wykres opisujący zmienność parametrów dynamicznych obiektu w trakcie jazdy po zamodelowanym torze.',
]

introduction_bank_2 = [
    'Symulacja została przeprowadzona dla następujących danych: {given_data}.',
    'Zaprezentowane wyniki numeryczne otrzymano dla danych równych:{given_data}.',
    'Eksperyment numeryczny został przeprowadzona dla następujących danych:{given_data}.',
    'Wyniki stmulacji numerycznych otrzymano dla:{given_data}.',
]

introduction_bank_3 = [
    'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu stwierdza się oscylacyjny charakter odpowiedzi układu. Dynamika systemu ma stabilny charakter',
    'Charakter zmian przebiegów poszczególnch współrzędnych jest stabilny. Otrzymane sygnały mieszczą się w zakresie dopuszczalnych uwarunkowaniami fizycznymi.',
    'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu mają drganiowy charakter i są zgodne z oczekiwaniami.',
    'Bazując na otrzymanych wynikach można stwierdzić wzajemną zależność poszczególnych wielkości dynamicznych. Dodatkowo wielkości modelu stabilizują się w czasie. '
]

intro_bank_composed_model_gen = RandomDescription(introduction_bank_1,
                                                  introduction_bank_2,
                                                  introduction_bank_3)

intro_bank_composed_model = [
    str(
        RandomDescription(introduction_bank_1, introduction_bank_2,
                          introduction_bank_3)) for obj in range(30)
]

conclusion_bank_x = [
    'Zauważa się {x(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({x(t)}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {x(t)_max} oraz {Derivative(x(t), t)_max} dla {x(t)_idmax} i {Derivative(x(t), t)_idmax}. ',
    'Zaobserwowano {x(t)_inf} parametru \({param_name}\) na maksymalną wartość pokonanej drogi oraz osiąganą wartość prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {x(t)_max} oraz {Derivative(x(t), t)_max} dla maksymalnej wartości badanego parametru równej {x(t)_idmax}. ',
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({x(t)}\), gdzie największa wartość pokonanej drogi to {x(t)_max}. W konsekwencji zaobserwowano {x(t)_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(x(t), t)}\), dla której minimalna wartość wynosi {Derivative(x(t), t)_min}, a największą osiąganą wartością jest {Derivative(x(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(x(t), t)_idmin} oraz {Derivative(x(t), t)_idmax}. ',
]

measurement_conclusion_bank_x = [
    'Zauważa się {a_ox_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_ox}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ mierzonego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {a_ox_max} dla {a_ox_idmax}. ',
    'Zaobserwowano {a_ox_inf} parametru \({param_name}\) na maksymalną wartość pokonanej drogi oraz osiąganą wartość prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {a_ox_max} dla maksymalnej wartości mierzonego parametru równej {a_ox_idmax}. ',
    #'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_ox}\), gdzie największa wartość zmierzonego przyspieszenia to {a_ox_max}. W konsekwencji zaobserwowano {a_ox_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(a_ox, t)}\), dla której minimalna wartość wynosi {Derivative(a_ox, t)_min}, a największą osiąganą wartością jest {Derivative(a_ox, t)_max} odpowiednio dla wartości parmametru: {Derivative(a_ox, t)_idmin} oraz {Derivative(a_ox, t)_idmax}. ',
]

measurement_conclusion_bank_x_rc = [
    'Zauważa się {a_rx_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_rx}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ mierzonego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {a_rx_max} dla {a_rx_idmax}. ',
    'Zaobserwowano {a_rx_inf} parametru \({param_name}\) na maksymalną wartość pokonanej drogi oraz osiąganą wartość prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {a_rx_max} dla maksymalnej wartości mierzonego parametru równej {a_rx_idmax}. ',
    #'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rx}\), gdzie największa wartość pokonanej drogi to {a_rx_max}. W konsekwencji zaobserwowano {a_rx_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(a_rx, t)}\), dla której minimalna wartość wynosi {Derivative(a_rx, t)_min}, a największą osiąganą wartością jest {Derivative(a_rx, t)_max} odpowiednio dla wartości parmametru: {Derivative(a_rx, t)_idmin} oraz {Derivative(a_rx, t)_idmax}. ',
]

conclusion_bank_varphi_rc = [
    'Zaobserwowano {varphi_RC(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań i prędkość kątową napędu RC. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi_RC(t)_min} oraz większej niż {varphi_RC(t)_max} odpowiednio dla wartości parametru: {varphi_RC(t)_idmin} oraz {varphi_RC(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_max}.',
    'Zmiana \({param_name}\) ma {varphi_RC(t)_inf} na wartość drgań i prędkość kątową napędu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi_RC(t)_min} oraz większej niż {varphi_RC(t)_max} odpowiednio dla wartości parametru: {varphi_RC(t)_idmin} oraz {varphi_RC(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_min}. ',
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({varphi_RC(t)}\), gdzie największa wartość amplitudy to {varphi_RC(t)_max}, a najmniejsza {varphi_RC(t)_min}. W konsekwencji zaobserwowano {varphi_RC(t)_inf} analizowanej zmiennej na wartość prędkości kątową napędu RC. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(varphi_RC(t), t)_idmin} oraz {Derivative(varphi_RC(t), t)_idmax}. ',
]

conclusion_bank_varphi_zwrc = [
    'Zaobserwowano {z_wrc(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań pionowych i prędkość napędu RC. Przemieszczenia pionowe nie przyjmują wartości mniejszej niż {z_wrc(t)_min} oraz większej niż {z_wrc(t)_max} odpowiednio dla wartości parametru: {z_wrc(t)_idmin} oraz {z_wrc(t)_idmax}. Dla prędkości drgań pionowych napędu minimalna wartość amplitudy to {Derivative(z_wrc(t), t)_min}, a największą osiąganą wartością jest {Derivative(z_wrc(t), t)_max}. ',
    'Zmiana \({param_name}\) ma {z_wrc(t)_inf} na wartość drgań i prędkość kątową analizowanego układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {z_wrc(t)_min} oraz większej niż {z_wrc(t)_max} odpowiednio dla wartości parametru: {z_wrc(z_wrct)_idmin} oraz {z_wrc(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(z_wrc(t), t)_min}, a największą osiąganą wartością jest {Derivative(z_wrc(t), t)_min}. ',
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({varphi(t)}\), gdzie największa wartość amplitudy to {varphi_RC(t)_max}, a najmniejsza {varphi(t)_min}. W konsekwencji zaobserwowano {varphi(t)_inf} analizowanej zmiennej na wartość prędkości drgań kątowych badanego układu, dla którego minimalna wartość amplitudy wynosi {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(varphi(t), t)_idmin} oraz {Derivative(varphi(t), t)_idmax}. '
]

conclusion_bank_phi = [
    'Zmiana \({param_name}\) ma {varphi(t)_inf} na wartość drgań i prędkość kątową układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi(t)_min} oraz większej niż {varphi(t)_max} odpowiednio dla wartości parametru: {varphi(t)_idmin} oraz {varphi(t)_idmax}. Dla prędkości kątowej wózka minimalna wartość amplitudy to {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_min}. ',
    'Zaobserwowano {varphi(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań i prędkość kątową układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi(t)_min} oraz większej niż {varphi(t)_max} odpowiednio dla wartości parametru: {varphi(t)_idmin} oraz {varphi(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_max}. ',
    'Zauważa się {varphi(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({varphi(t)}\) oraz odpowiadającej temu przemieszczniu kątowemu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {varphi(t)_max} oraz {Derivative(varphi(t), t)_max} dla {varphi(t)_idmax} i {Derivative(varphi(t), t)_idmax}. ',
]

conclusion_bank_z = [
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({z(t)}\), gdzie największa wartość amplitudy to {z(t)_max}, a najmniejsza {z(t)_min}. W konsekwencji zaobserwowano {z(t)_inf} analizowanej zmiennej na wartość drgań prędkości drgań pionowych, dla których minimalna wartość amplitudy to {Derivative(z(t), t)_min}, a największą osiąganą wartością jest {Derivative(z(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(z(t), t)_idmin} oraz {Derivative(z(t), t)_idmax}. ',
    'Zauważa się {z(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({z(t)}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {z(t)_max} oraz {Derivative(z(t), t)_max} dla {z(t)_idmax} i {Derivative(z(t), t)_idmax}. ',
    'Zaobserwowano {z(t)_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {z(t)_max} oraz {Derivative(z(t), t)_max} dla maksymalnej wartości badanego parametru równej {z(t)_idmax}. ',
]

measurement_conclusion_bank_z = [
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_oz}\), gdzie największa wartość amplitudy to {a_oz_max}, a najmniejsza {a_oz_min}.',
    'Zauważa się {a_oz_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_oz}\) oraz odpowiadającej temu przemieszczniu prędkości.',
    'Zaobserwowano {a_oz_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości.'
]

measurement_conclusion_bank_z_rc = [
    'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rz}\), gdzie największa wartość amplitudy to {a_rz_max}, a najmniejsza {a_rz_min}.',
    'Zauważa się {a_rz_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_rz}\) oraz odpowiadającej temu przemieszczniu prędkości.',
    'Zaobserwowano {a_rz_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości.',
]

conclusion_bank_no_impact = [
    'Nie odmotowano natomiast istotnego wpływu tego parametru na pozostałe stopnie swobody układu. ',
    'Nie zauważono kluczowych zmian pozostałych zmiennych układu pod wpływem badanych wartości rozpatrywanego parametru. ',
    'Pozostałe zmienne układu nie wykazują wpływu oddziałującego w sposób istoty na dynamikę rozpatrywanego systemu. ',
    'Badany system dynamiczny nie wykazuje szczególnej wrażlowości na zmiany innych analizowanych parametrów. ',
]

measurement_conclusion_bank_no_impact = [
    'Nie odmotowano natomiast istotnego wpływu tego parametru na pozostałe stopnie swobody wózka. ',
    'Nie zauważono kluczowych zmian pozostałych zmiennych wózka pod wpływem badanych wartości rozpatrywanego parametru. ',
    'Pozostałe zmienne układu nie wykazują wpływu oddziałującego w sposób istoty na dynamikę rozpatrywanego wózka. ',
    'Badany wózek nie wykazuje szczególnej wrażlowości na zmiany innych analizowanych parametrów. ',
]

measurement_conclusion_bank_summary = [
    'Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty pomiarów. Dla uzyskanych danych liczbowych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego wózka. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).'
]

conclusion_bank_composed_model_gen = RandomDescription(
    conclusion_bank_varphi_rc, conclusion_bank_x, conclusion_bank_z,
    conclusion_bank_phi, conclusion_bank_no_impact)

conclusion_bank_composed_model = [
    str(
        RandomDescription(conclusion_bank_varphi_rc, conclusion_bank_x,
                          conclusion_bank_z, conclusion_bank_phi,
                          conclusion_bank_no_impact)) for obj in range(30)
]

measurement_conclusion_bank_composed_model = [
    str(
        RandomDescription(
            measurement_conclusion_bank_x_rc,
            measurement_conclusion_bank_x,
            measurement_conclusion_bank_z,
            measurement_conclusion_bank_z_rc,
            measurement_conclusion_bank_no_impact,
        )) for obj in range(30)
]

conclusion_bank_chair_model_gen = RandomDescription(conclusion_bank_x,
                                                    conclusion_bank_z,
                                                    conclusion_bank_phi,
                                                    conclusion_bank_no_impact)


conclusion_bank_chair_move_model_gen = RandomDescription(conclusion_bank_x,
                                                    conclusion_bank_phi,
                                                    conclusion_bank_no_impact)

conclusion_bank_chair_model = [
    str(
        RandomDescription(conclusion_bank_x, conclusion_bank_z,
                          conclusion_bank_phi, conclusion_bank_no_impact))
    for obj in range(30)
]

measurement_conclusion_bank_chair_model = [
    str(
        RandomDescription(measurement_conclusion_bank_x,
                          measurement_conclusion_bank_x_rc,
                          measurement_conclusion_bank_z,
                          measurement_conclusion_bank_z_rc,
                          measurement_conclusion_bank_no_impact))
    for obj in range(30)
]

conclusion_bank_drive_model_gen = RandomDescription(conclusion_bank_varphi_rc,
                                                    conclusion_bank_no_impact)

conclusion_bank_drive_model = [
    str(RandomDescription(conclusion_bank_varphi_rc,
                          conclusion_bank_no_impact)) for obj in range(30)
]

measurement_conclusion_bank_drive_model = [
    str(
        RandomDescription(measurement_conclusion_bank_x_rc,
                          measurement_conclusion_bank_no_impact))
    for obj in range(30)
]

analysis_intro_ending = '''Tabela {nr_rys} przedstawia zakres parametrów przyjętych do wykonania symulacji numerycznych. Oceniono, że przyjęty zakres badnia odpowiada możliwym do uzyskania w praktyce wartością i przeprowadzone symulacje będę dobrze reprezentować dynamikę układu.
'''

simulations_summary_str = ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty symulacji numerycznych. Dla uzyskanych danych symulacyjnych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego układu. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''


class PlottedData(Figure):

    _latex_name = 'figure'

    def __init__(self,
                 numerical_data,
                 fig_name,
                 *,
                 preview=False,
                 position=None,
                 **kwargs):
        super().__init__(position=position, **kwargs)

        self._numerical_data = numerical_data
        self.fig_name = str(fig_name)
        #self._latex_name='figure' #super()._latex_name
        self.preview = preview

    def add_data_plot(self, numerical_data=None):
        numerical_data = self._numerical_data

        ax = numerical_data.plot(subplots=True, figsize=(10, 7))
        #ax=solution_tmp.plot(subplots=True)
        ([
            ax_tmp.legend(['$' + vlatex(sym) + '$'])
            for ax_tmp, sym in zip(ax, numerical_data.columns)
        ])
        plt.savefig(self.fig_name + '.png')
        self.add_image(self.fig_name, width=NoEscape('15cm'))

        if self.preview == True:
            plt.show()

        plt.close()


class DataTable(Table):
    _latex_name = 'table'

    def __init__(self, numerical_data, position=None):
        super().__init__(position=position)
        #print(numerical_data)
        self._numerical_data = numerical_data
        self.position = position

    def add_table(self, numerical_data=None):
        self.append(NoEscape('%%%%%%%%%%%%%% Table %%%%%%%%%%%%%%%'))
        #         if numerical_data!=None:
        #             self._numerical_data=numerical_data

        tab = self._numerical_data
        self.append(NoEscape(tab.to_latex(index=False, escape=False)))


class ReportSection(Section):

    _latex_name = 'section'

    def __init__(self,
                 title,
                 results_frame=None,
                 analysis_key='',
                 analysis_name='dynamics',
                 results_col_name='simulations',
                 preview=True):

        self.analysis_key = analysis_key
        self.analysis_name = analysis_name
        self.sims_name = results_col_name
        self.preview = preview
        super().__init__(title)

    def influeance_level(self, data):

        sensitivity_level = {
            0.02: 'znikomy wpływ',
            0.1: 'mały wpływ',
            0.5: 'średni wpływ',
            0.7: 'znaczny wpływ',
            0.9: 'istotny wpływ'
        }

        indicator = round(
            2 * float(
                abs((data.abs().max() - data.abs().min()) /
                    (0.001 + data.abs().max() + data.abs().min()))), 1) / 2

        select = sensitivity_level

        a = np.asarray(list(select.keys()))

        print('indicator: ', indicator)
        print(select[a.flat[np.abs(a - indicator).argmin()]])
        print(select[a.flat[np.abs(a - indicator).argmin()]])

        return select[a.flat[np.abs(a - indicator).argmin()]]

    def match_unit_to_data(self, parameter, value, units_dict=None):

        if parameter in units_dict:
            return value * units_dict[parameter]
        else:
            return value * pint.Quantity(1)

    def get_feature_dict(self,
                         numerical_data=None,
                         given_data_dict=None,
                         units_dict={},
                         marker=None):

        feature_dict = {
            'param_name': vlatex(self.analysis_key),
        }
        if marker:
            feature_dict.update({'nr_rys': Ref(marker).dumps()})

        if type(numerical_data) != type(None):

            column_names = numerical_data.columns
            feature_dict.update(
                {str(name): vlatex(name)
                 for name in column_names})
            feature_dict.update({
                str(name) + '_max':
                '{val:~Lx}'.format(val=self.match_unit_to_data(
                    name,
                    numerical_data[name].round(1).max(),
                    units_dict=units_dict))
                for name in column_names
            })
            feature_dict.update({
                str(name) + '_min':
                '{val:~Lx}'.format(val=self.match_unit_to_data(
                    name,
                    numerical_data[name].round(1).min(),
                    units_dict=units_dict))
                for name in column_names
            })
            feature_dict.update({
                str(name) + '_idmax':
                '{val:~Lx}'.format(val=self.match_unit_to_data(
                    name,
                    numerical_data[name].round(1).idxmax(),
                    units_dict=units_dict))
                for name in column_names
            })
            feature_dict.update({
                str(name) + '_idmin':
                '{val:~Lx}'.format(val=self.match_unit_to_data(
                    name,
                    numerical_data[name].round(1).idxmin(),
                    units_dict=units_dict))
                for name in column_names
            })

        if type(given_data_dict) != type(None):
            feature_dict.update({
                'given_data':
                '\(' + '\), \('.join([
                    vlatex(lhs) +
                    '={rhs:~Lx.2f}'.format(rhs=self.match_unit_to_data(
                        lhs, rhs, units_dict=units_dict))
                    for lhs, rhs in given_data_dict.items()
                ]) + '\)'
            })

        return feature_dict

    def add_introduction(
            self,
            title='Zakres prowadzonej analizy',
            numerical_data=None,
            units_dict={},
            initial_description='Initial description',  #random.choice(sec_intro_composed_model)
            ending_summary='Ending summary',  # random.choice(sec_intro_composed_model)
    ):

        with self.create(Subsection(NoEscape(title))) as subsec:

            tab_mrk = Marker('given_data_' + str(self.analysis_name) + '_' +
                             str(self.analysis_key),
                             prefix='tab')

            symbols_df = pd.DataFrame(data=[
                {
                    'Zmienna': '$ {eq}  $'.format(eq=vlatex(key)),
                    'Liczba przypadków': len(column.unique()),
                    'Wartość minimalna': column.round(2).min(),
                    'Wartość maksymalna': column.round(2).max(),
                    'Jednostka': '$ {unit:~Lx}  $'.format(
                        unit=units_dict[key]),
                } for key, column in numerical_data.items()
            ])

            format_dict = (self.get_feature_dict(numerical_data=None,
                                                 given_data_dict=None,
                                                 units_dict=units_dict,
                                                 marker=tab_mrk))

            subsec.append(NoEscape(initial_description.format(**format_dict)))

            with subsec.create(DataTable(symbols_df, position='H')) as table:

                table.add_caption(
                    NoEscape(
                        'Zestawienie parametrów modelu'.format(**format_dict)))
                table.add_table(numerical_data.T)

            subsec.append(NoEscape(ending_summary.format(**format_dict)))

            #subsec.append(NewPage())
            #subsec.append(NoEscape('\\'))
            subsec.append(NoEscape('\par'))
            
    def add_simulational_results(
            self,
            title='Wyniki symulacji numerycznych',
            numerical_data=None,
            units_dict={},
            initial_description='Initial description',  #random.choice(intro_bank_composed_model)
            ending_summary='Ending summary',  #random.choice(summary_bank_composed_model)
    ):

        with self.create(Subsection(title)) as subsec:

            simulation_results_frame = numerical_data

            for key, row in simulation_results_frame.iterrows():

                data_with_units = {
                    parameter: value
                    for parameter, value in row.items()
                    if isinstance(parameter, Symbol)
                }

                current_time = dtime.datetime.now().timestamp()
                current_fig_mrk = Marker('data_plot_' +
                                         str(self.analysis_key) + '_' +
                                         str(current_time),
                                         prefix='fig')

                format_dict = {
                    **(self.get_feature_dict(numerical_data=row['simulations'],
                                             given_data_dict=data_with_units,
                                             units_dict=units_dict,
                                             marker=current_fig_mrk)),
                    #**{str(name):vlatex(name) for name in row['simulations'].columns}
                }

                subsec.append(
                    NoEscape(initial_description.format(**format_dict)))

                with subsec.create(
                        PlottedData(row['simulations'],
                                    './plots/fig_' + str(current_time),
                                    position='H',
                                    preview=self.preview)) as fig:
                    fig.add_data_plot(row['simulations'])
                    fig.add_caption(
                        NoEscape(
                            #'Przebiegi czasowe modelu dla danych: {given_data}.'.format(**format_dict)
                            'Przebiegi czasowe modelu dla rozważanego zakresu danych.'
                            .format(**format_dict)))
                    fig.append(Label(current_fig_mrk))

                subsec.append(NoEscape(ending_summary.format(**format_dict)))

                #subsec.append(NewPage())
                #subsec.append(NoEscape('\\'))
                subsec.append(NoEscape('\par'))

    def prepare_summary_data(self,simulation_results_frame):
        summary_frame = pd.DataFrame(
                data=[
                    a.max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=[
                    str(elem) for elem in simulation_results_frame[
                        self.analysis_key].round(2)
                ])
        return summary_frame

    
    def add_summary(
            self,
            title='Analiza otrzymanych wyników',
            numerical_data=None,
            units_dict={},
            initial_description='Initial description',
            ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
    ):
        ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty symulacji numerycznych. Dla uzyskanych danych symulacyjnych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego układu. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''

        summary_frame = numerical_data
        with self.create(Subsection(title)) as subsec:

            current_time = dtime.datetime.now().timestamp()
            summary_mrk = Marker('summary_' + str(self.analysis_name) + '_' +
                                 str(self.analysis_key),
                                 prefix='fig')
            
            

#             summary_frame = pd.DataFrame(
#                 data=[
#                     a.max()
#                     for a in simulation_results_frame['simulations'].values
#                 ],
#                 index=[
#                     str(elem) for elem in simulation_results_frame[
#                         self.analysis_key].round(2)
#                 ])

            format_dict = (self.get_feature_dict(numerical_data=summary_frame,
                                                 given_data_dict=None,
                                                 units_dict=units_dict,
                                                 marker=summary_mrk))

            print(format_dict)

            subsec.append(NoEscape(initial_description.format(**format_dict)))

            step_key_influence_dict = {
                str(name) + '_inf': self.influeance_level(summary_frame[name])
                for name in summary_frame.columns
            }

            print(step_key_influence_dict)

            with subsec.create(
                    PlottedData(summary_frame,
                                './plots/fig_summary_' + str(current_time),
                                position='H',
                                preview=True)) as fig:
                fig.add_data_plot(summary_frame)
                fig.add_caption(
                    NoEscape(
                        'Zestawienie wyników przeprowadzonej analizy.'
                        .format(**format_dict)))
                fig.append(Label(summary_mrk))

            subsec.append(
                NoEscape(
                    ending_summary.format(
                        **{
                            **format_dict,
                            **step_key_influence_dict,
                            #**{str(name):vlatex(name) for name in row['simulations'].columns}
                        })))

            #subsec.append(NewPage())
            #subsec.append(NoEscape('\\'))
            subsec.append(NoEscape('\par'))

    def add_optimization_result(
            self,
            title='Wyniki optymalizacji jednokryterialnej',
            numerical_data=None,
            units_dict={},
            initial_description='Initial description',
            ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
    ):
        
        
        
        self.add_summary(
            title=title,
            numerical_data=numerical_data,
            units_dict=units_dict,
            initial_description=initial_description,
            ending_summary=
            'Ending summary',  # random.choice(conclusion_bank_composed_model)
        )

    def add_measurement_summary(
            self,
            title='Analiza otrzymanych wyników pomiarów',
            numerical_data=None,
            units_dict={},
            initial_description='Initial description',
            ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
    ):
        ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty pomiarów. Dla uzyskanych danych liczbowych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego wózka. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''

        simulation_results_frame = numerical_data
        with self.create(Subsection(title)) as subsec:

            current_time = dtime.datetime.now().timestamp()
            summary_mrk = Marker('measurment_summary_' +
                                 str(self.analysis_name) +
                                 str(self.analysis_key),
                                 prefix='fig')

            summary_frame = pd.DataFrame(
                data=[
                    a.max()
                    for a in simulation_results_frame['simulations'].values
                ],
                index=[
                    str(elem) for elem in simulation_results_frame[
                        self.analysis_key].round(2)
                ])

            format_dict = (self.get_feature_dict(numerical_data=None,
                                                 given_data_dict=None,
                                                 units_dict=units_dict,
                                                 marker=summary_mrk))

            subsec.append(NoEscape(initial_description.format(**format_dict)))

            format_dict = (self.get_feature_dict(numerical_data=summary_frame,
                                                 given_data_dict=None,
                                                 units_dict=units_dict,
                                                 marker=summary_mrk))

            step_key_influence_dict = {
                str(name) + '_inf': self.influeance_level(summary_frame[name])
                for name in summary_frame.columns
            }

            with subsec.create(
                    PlottedData(summary_frame,
                                './plots/fig_summary_' + str(current_time),
                                position='H',
                                preview=True)) as fig:
                fig.add_data_plot(summary_frame)
                fig.add_caption(
                    NoEscape(
                        'Wpływ zmiany rozważanego parametru na dynamikę układu.'
                        .format(**format_dict)))
                fig.append(Label(summary_mrk))

            subsec.append(
                NoEscape(
                    ending_summary.format(
                        **{
                            **format_dict,
                            **step_key_influence_dict,
                            #**{str(name):vlatex(name) for name in row['simulations'].columns}
                        })))

            subsec.append(NewPage())
