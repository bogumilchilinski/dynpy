from .en import *


class TitlePageComponent(Environment):
    """
    - university: np. "Politechnika Warszawska"
    - faculty:    np. "Wydział Samochodów i Maszyn Roboczych"
    - thesis_type: np. "Praca Przejściowa"
    - thesis_title: np. "Dopracowanie biblioteki dynpy"
    - author:     np. "Mike Wazowski"
    - supervisor: np. "Jeremy Clarkson"
    - place:      np. "Warszawa"
    - year:       np. 2025
    """

    latex_name = "titlepage"

    # example use:
    # doc.append(TitlePageComponent(
    #     university='Politechnika Warszawska',
    #     faculty='Wydział Samochodów i Maszyn Roboczych',
    #     thesis_type='EXAMPLE',
    #     thesis_title='EXAMPLE',
    #     author='Richard Feynman',
    #     supervisor='Enrico Fermi',
    #     place='Warszawa',
    #     year=2025
    # ))

    def __init__(
        self,
        university: str,
        faculty: str,
        thesis_type: str,
        thesis_title: str,
        author: str,
        supervisor: str,
        place: str,
        year: int,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        # wyrównanie do środka
        self.append(NoEscape(r"\centering"))
        # uczelnia
        self.append(NoEscape(r"\Large " + university + r"\\[0.5cm]"))
        # wydział
        self.append(NoEscape(r"\large " + faculty + r"\\[1cm]"))
        # rodzaj pracy
        self.append(NoEscape(r"\textbf{" + thesis_type + r"}\\[1.5cm]"))
        # tytuł pracy
        self.append(NoEscape(r'\LARGE "' + thesis_title + r'"\\[2cm]'))
        # autor i opiekun
        self.append(NoEscape(r"\normalsize Autor: " + author + r"\\"))
        self.append(NoEscape(r"Opiekun: " + supervisor + r"\\[1.5cm]"))
        # miejsce i rok
        self.append(NoEscape(place + ", " + str(year)))


######################### Initial point


# Damian
class ExemplaryPictureComponent(ReportComponent):

    title = "Przykład rzeczywistego obiektu"
    packages = [Package("float")]

    @property
    def header_text(self):

        return "Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."

    @property
    def footer_text(self):

        return "Model dynamiczny układu określa się na podstawie analizy rozważanego przypadku. Należy pamiętać, że stopień odwzorowania (poziom abstrakcji) modelu zależy od tego do czego planuje się go używać."

    def append_elements(self):

        system = self._system

        display(ReportText(self.header_text))

        with self.create(Figure(position="H")) as fig:
            fig.add_image(system._real_example(), width="8cm")

        display(ReportText(self.footer_text))


# Boogi
class SchemeComponent(ExemplaryPictureComponent):
    title = "Schemat systemu"

    @property
    def header_text(self):
        #
        return "Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu."

    @property
    def footer_text(self):

        system = self._system

        return f"Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}."

    def append_elements(self):

        system = self._system

        display(ReportText(self.header_text))

        with self.create(Figure(position="H")) as fig:
            fig.add_image(system._scheme(), width="8cm")

        display(ReportText(self.footer_text))


class NumericalAnalysisComponent(ExemplaryPictureComponent):
    title = "Symulacja numeryczna"

    _default_ics = None
    _default_parameters = {}
    _default_tspan = np.linspace(0, 1, 1000)

    @classmethod
    def set_default_ics(cls, default_ics=None):

        if default_ics is None:
            cls._default_ics = default_ics

        return cls

    @classmethod
    def set_default_parameters(cls, default_parameters=None):

        if default_parameters is None:
            cls._default_parameters = default_parameters

        return cls

    @classmethod
    def set_default_tspan(cls, default_tspan=None):

        if default_tspan is None:
            cls._default_tspan = default_tspan

        return cls

    def append_elements(self):

        system = self._system

        display(
            ReportText(
                f"""Dla Damiana :P
                            """
            )
        )
        if self._default_ics is None:
            ics = list(system.Y * 0.0)
        else:
            ics = self._default_ics

        sym_res = (
            system.subs(self._default_parameters)
            .numerized()
            .compute_solution(self._default_tspan, ics)
        )

        LatexDataFrame.formatted(sym_res).plotted(preview=True)

        display(
            ReportText(
                f"""Dla Damiana :P
                            """
            )
        )


# Boogi
class KineticEnergyComponent(ReportComponent):

    title = "Energia kinetyczna"

    @property
    def header_text(self):
        # "Energia kinetyczna układu wyrażona jest wzorem:"
        return "Energia kinetyczna układu wyrażona jest wzorem:"

    @property
    def footer_text(self):
        # "Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych)."
        return "Determined formula specify energy of the system related to its inertial properties."

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(SympyFormula(Eq(Symbol("T"), dyn_sys_lin._kinetic_energy), marker=None))

        display(ReportText(self.footer_text))


class PotentialEnergyComponent(ReportComponent):  # Jaś fasola

    title = "Energia potencjalna"

    @property
    def header_text(self):

        return "Energia potencjalna układu wyrażona jest wzorem:"

    @property
    def footer_text(self):

        return "Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(
            SympyFormula(Eq(Symbol("V"), dyn_sys_lin._potential_energy), marker=None)
        )

        display(ReportText(self.footer_text))


# Marcel spolszczone
class DissipationComponent(ReportComponent):

    title = "Dyssypacyjna funkcja Rayleigh'a"

    @property
    def header_text(self):

        return "Energia rozpraszana tłumieniem wyrażona jest wzorem:"

    @property
    def footer_text(self):

        return "Podana zależność stanowi potencjał dysynpacyjny Rayleigh'a, który poddany różniczkowaniu względem wektora prędkości uogólnionych pozwala na określenie sił wiskotycznego tłumienia."

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(
            SympyFormula(
                Eq(Symbol("D"), dyn_sys_lin._dissipative_potential), marker=None
            )
        )

        display(ReportText(self.footer_text))


# Boogi
class LagrangianComponent(ReportComponent):

    # title="Lagrangian (funkcja Lagrange'a)  układu"
    title = "Lagrangian of the system (Lagranges function)"

    @property
    def header_text(self):
        # "Lagrangian systemu wyrażony jest wzorem ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):"

        system = self._system
        dyn_sys = system

        return f"System Lagrangian is described by the formula ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):"

    @property
    def body_text_1(self):
        # "Równania Eulera-Lagrange'a dla rozważanego przypadku są nastęujące: "
        return "The Euler-Lagrange equations for the case under consideration are as follows: "

    @property
    def body_text_2(self):
        # "Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: "
        return "Subsequent derivatives obtained with the Euler-Lagrange equations are as follows: "

    @property
    def footer_text(self):
        # "Wyniki przedstawionych operacji wykorzystuje się wyznaczenia równań ruchu układu."
        return "The results of the presented operations are used to determine the equations of motion of the system."

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        print(system._scheme())

        mrk_lagrangian_nonlin = Marker("lagrangianNL", prefix="eq")

        # display(ReportText(f'''The following model is considered. The system's Lagrangian is described by the formula ({Ref(mrk_lagrangian_nonlin).dumps()}):
        #                    '''))
        display(ReportText(self.header_text))

        display(
            (
                SympyFormula(
                    Eq(Symbol("L"), dyn_sys.L.expand()[0]), marker=mrk_lagrangian_nonlin
                )
            )
        )

        q_sym = [Symbol(f"{coord}"[0:-3]) for coord in dyn_sys.q]

        diffL_d = lambda coord: Symbol(
            f"\\frac{{ \\partial L}}{{  \\partial {vlatex(coord)}  }}"
        )
        diffD_d = lambda coord: Symbol(
            f"\\frac{{ \\partial D}}{{  \\partial {vlatex(diff(coord))}  }}"
        )
        d_dt_diffL_d = lambda coord: Symbol(
            f"\\frac{{ \\mathrm{{d}}  }}{{  \\mathrm{{d}} {vlatex(system.ivar)}  }} {vlatex(diffL_d(coord))} "
        )

        display(ReportText(self.body_text_1))

        for coord in dyn_sys.q:
            display(
                (
                    SympyFormula(
                        Eq(
                            d_dt_diffL_d(coord.diff(system.ivar))
                            - diffL_d(coord)
                            + diffD_d(coord),
                            Symbol(f"Q_{{ {vlatex(coord)} }}^N"),
                        ),
                        marker=mrk_lagrangian_nonlin,
                        backend=vlatex,
                    )
                )
            )

        display(ReportText(self.body_text_2))

        for coord in dyn_sys.Y:
            display(
                (
                    SympyFormula(
                        Eq(diffL_d(coord), dyn_sys.L.expand()[0].diff(coord)),
                        marker=mrk_lagrangian_nonlin,
                        backend=vlatex,
                    )
                )
            )

        for coord in dyn_sys.q.diff(system.ivar):
            display(
                (
                    SympyFormula(
                        Eq(
                            d_dt_diffL_d(coord),
                            dyn_sys.L.expand()[0].diff(coord).diff(system.ivar),
                        ),
                        marker=mrk_lagrangian_nonlin,
                        backend=vlatex,
                    )
                )
            )
        for coord in dyn_sys.q:
            display(
                (
                    SympyFormula(
                        Eq(
                            diffD_d(coord),
                            (
                                S.One
                                / 2
                                * diff(dyn_sys.q.transpose())
                                * dyn_sys_lin.damping_matrix()
                                * diff(dyn_sys.q)
                            )[0].diff(diff(coord)),
                        ),
                        marker=mrk_lagrangian_nonlin,
                        backend=vlatex,
                    )
                )
            )
            # display(Markdown(f'\\begin{equation}    \\end{equation}').reported())
        # with doc_model.create(DMath()) as eq:
        #    eq.append(NoEscape(latex(Derivative(Symbol('L'),q_sym[0],evaluate=False))))
        #    eq.append(NoEscape('='))
        #    eq.append(NoEscape(vlatex(dyn_sys.L.expand()[0].diff(dyn_sys.q[0]))))

        mrk_gov_eq_nonlin = Marker("gov_eq_nonlin_sys", prefix="eq")

        # display(ReportText(f'''The governing equations of the system have a following form ({Ref(mrk_gov_eq_nonlin).dumps()}):
        #                    '''))

        display(ReportText(self.footer_text))

    def as_frame(self):
        frame = Frame(title=self.title, options=["allowframebreaks"])
        # frame.packages +(self.packages)
        frame += list(self)
        return frame

    def as_frame(self):
        frame = Frame(title=self.title, options=["allowframebreaks"])
        # frame.packages +(self.packages)
        frame += list(self)
        return frame


# Amadi & Damian
class GoverningEquationComponent(ReportComponent):

    # Równania ruchu
    title = "Equation of motion"

    @property
    def entry_text_one(self):
        system = self._system
        dyn_sys = system
        return f"Wykorzystując obliczone pochodne, wyznacza się równanie ruchu na podstawie odpowiedniego wzoru. Równanie ruchu układu przedstawia zależność: ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})"

    @property
    def entry_text_two(self):
        system = self._system
        dyn_sys = system

        return f"Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru. Równania ruchu układu przedstawiają zależności: ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-     ({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"

    @property
    def footer_text(self):
        # f''' Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. '''
        return f"""  Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. """

    def append_elements(self):

        system = self._system
        dyn_sys = system

        if len(system.q) == 1:
            display(ReportText(self.entry_text_one))
        else:
            display(ReportText(self.entry_text_two))

        for eq in dyn_sys._eoms:
            display(SympyFormula(Eq(eq.simplify().expand(), 0), marker=None))

        display(ReportText(self.footer_text))


class LinearizedGoverningEquationComponent(ReportComponent):
    # Równania ruchu
    title = "Equation of motion"

    @property
    def entry_text_one(self):
        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        return f"Wykorzystując obliczone pochodne, wyznacza się równanie ruchu na podstawie odpowiedniego wzoru. Zlinearyzowane równanie ruchu układu przedstawia zależność: ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})"

    @property
    def entry_text_two(self):
        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        return f"Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru. Zlinearyzowane równania ruchu układu przedstawiają zależności: ({AutoMarker(Eq(dyn_sys_lin._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys_lin._eoms[-1].simplify().expand(),0))})"

    @property
    def footer_text(self):

        return f""" Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu. Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych. """

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        if len(system.q) == 1:
            display(ReportText(self.entry_text_one))
        else:
            display(ReportText(self.entry_text_two))

        for eq in dyn_sys._eoms:
            display(SympyFormula(Eq(eq.simplify().expand(), 0), marker=None))

        display(ReportText(self.footer_text))


class LinearizationComponent(ReportComponent):  # Szymon

    title = "Linearyzacja równań ruchu"

    @property
    def entry_text(self):
        return "Linearyzaja równań polega na znalezieniu ich rozwinięcia w szereg Taylora względem współrzędnych, prędkości i przyspieszeń uogólnionych w otoczeniu punktu równowagi. Celem uproszczenia wprowadzono następujące oznaczenia:"

    @property
    def equilibrium_point_text(self):
        return "Punkty równowagi rozważanego układu są następujące:"

    @property
    def eom_text(self):
        dyn_sys = self._system
        return """Równanie ruchu dla współrzędnej ${latex(dyn_sys.q[no])}$ można przestawić jako:"""

    @property
    def lagrange_text(self):
        dyn_sys = self._system
        return "Formalnie należy obliczyć pochodne cząstkowe wielkości uogólnionych ze składników równań Lagrange'a"

    @property
    def derivative_text(self):
        return "Poszczególne pochodne mają następującą postać:"

    @property
    def linearized_eq_text(self):
        return "Po podstawieniu obliczonych pochodnych, otrzumuje się następujące zlinearyzowane równanie:"

    def append_elements(self):

        system = self._system
        ReportText.set_directory("./SDAresults")
        latex_store = AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        mrk_lagrangian_nonlin = Marker("lagrangLin", prefix="eq")
        mrk_lagrangian_lin = Marker("lagrangLin", prefix="eq")

        t = system.ivar

        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()
        coords = tuple(list(dyn_sys.Y) + list(dyn_sys.q.diff(t, t)))
        op_point = {coord: 0 for coord in coords}
        op_point.update(dyn_sys._op_points(subs=True)[0])

        display(ReportText(self.entry_text))

        for coord in coords:
            display(
                (
                    SympyFormula(
                        Eq(Symbol(vlatex(coord)), Symbol(latex(coord))), marker=None
                    )
                )
            )

        display(ReportText(self.equilibrium_point_text))

        for eq_coord, val in op_point.items():
            display((SympyFormula(Eq(eq_coord, val), marker=mrk_lagrangian_lin)))

        diffL_d = lambda coord: Symbol(
            latex(Derivative(Symbol("L"), Symbol(vlatex(coord))))
        )

        for no, eom in enumerate(dyn_sys._eoms):

            eq_sym = Symbol(f"RR_{latex(dyn_sys.q[no])}")

            display(ReportText(self.eom_text.format(coord=latex(dyn_sys.q[no]))))

            # display((SympyFormula(  Eq(eq_sym,eom,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display(ReportText(self.lagrange_text))

            display(
                (
                    SympyFormula(
                        Eq(
                            comps.MultivariableTaylorSeries(
                                eq_sym, coords, n=1, x0=op_point
                            )._symbolic_sum(),
                            0,
                        ),
                        marker=None,
                        backend=latex,
                    )
                )
            )

            diff_list = comps.MultivariableTaylorSeries(
                eom, coords, n=1, x0=op_point
            ).calculation_steps(expr_symbol=eq_sym)

            display(ReportText(self.derivative_text))

            for diff_eq in diff_list:

                display(
                    (SympyFormula(diff_eq, marker=mrk_lagrangian_lin, backend=latex))
                )

            display(ReportText(self.linearized_eq_text))
            display(
                (
                    SympyFormula(
                        Eq(
                            comps.MultivariableTaylorSeries(
                                eom, coords, n=1, x0=op_point
                            )
                            .doit()
                            .expand()
                            .simplify()
                            .expand(),
                            0,
                            evaluate=False,
                        ),
                        marker=mrk_lagrangian_lin,
                        backend=latex,
                    )
                )
            )

        AutoBreak.latex_backend = latex_store


# Marcel & Monika spolszczone
class FundamentalMatrixComponent(ReportComponent):

    title = "Wyznaczanie macierzy fundamentalnej"

    @property
    def header_text(self):
        return "Z równań ruchu wyznaczono macierz mas i sztywności układu::"

    @property
    def body_text(self):

        return "Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu ${latex(Delta)}$, przedstawiają się następująco::"

    @property
    def footer_text(self):

        return "Macierz fundamentalna pozwala określić rozwiązanie ustalone. Natomiast bazując na równaniu charakterystycznym określa się częstości własne układu."

    def append_elements(self):

        system = self._system
        ReportText.set_directory("./SDAresults")

        latex_store = AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store

        t = system.ivar

        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        display(ReportText(self.header_text))

        display(
            (
                SympyFormula(
                    Eq(Symbol("M"), dyn_sys_lin.inertia_matrix(), evaluate=False),
                    marker="a",
                )
            )
        )

        display(
            (
                SympyFormula(
                    Eq(Symbol("K"), dyn_sys_lin.stiffness_matrix(), evaluate=False),
                    marker="a",
                )
            )
        )

        Delta = Symbol("\Delta")

        display(ReportText(self.body_text))

        display(
            (
                SympyFormula(
                    Eq(Symbol("A"), dyn_sys_lin.fundamental_matrix(), evaluate=False),
                    marker="a",
                )
            )
        )
        display(
            (
                SympyFormula(
                    Eq(
                        Delta,
                        dyn_sys_lin.fundamental_matrix()
                        .det()
                        .expand()
                        .simplify()
                        .simplify()
                        .expand(),
                        evaluate=False,
                    ),
                    marker="a",
                    backend=latex,
                )
            )
        )

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store


# Mateusz
class GeneralSolutionComponent(ReportComponent):
    title = "Rozwiązanie ogólne"

    @property
    def header_text(self):

        return "Rozwiązanie ogólne przedstawia wyrażenie:"

    @property
    def footer_text(self):

        return "Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu."

    def append_elements(self):

        from ....dynamics import HarmonicOscillator, LagrangesDynamicSystem

        # from dynpy import LagrangesDynamicSystem, HarmonicOscillator

        system = self._system
        ReportText.set_directory("./SDAresults")

        latex_store = AutoBreak.latex_backend

        AutoBreak.latex_backend = latex_store

        t = system.ivar

        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        display(ReportText(self.header_text))
        display(
            (
                SympyFormula(
                    Eq(
                        Symbol("X"),
                        HarmonicOscillator(dyn_sys_lin.linearized())
                        .general_solution()
                        .n(3),
                        evaluate=False,
                    ),
                    marker="a",
                    backend=latex,
                )
            )
        )

        display(ReportText(self.footer_text))

        AutoBreak.latex_backend = latex_store


# Grześ
class FrequencyResponseFunctionComponent(ReportComponent):
    # "amplitude-frequency characteristic"
    title = "Charakterystyka Amplitudowo-Częstotliwościowa"

    @property
    def header_text(self):
        # "funkcja odpowiedzi częstotliwościowej:"
        return "funkcja odpowiedzi częstotliwościowej:"

    @property
    def footer_text(self):
        # "jest to suma kwadratów amplitud pod pierwiastkiem:"
        return "jest to suma kwadratów amplitud pod pierwiastkiem:"

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(
            SympyFormula(
                Eq(Symbol("frf"), dyn_sys.frequency_response_function()), marker=None
            )
        )

        display(ReportText(self.footer_text))


FRFComponent = FrequencyResponseFunctionComponent


class FrequencyResponseFunctionComponentToSecond(ReportComponent):

    title = "sens mocy układu"

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(
            ReportText(
                f"""
                           sens mocy układu:
                           """
            )
        )

        display(
            SympyFormula(
                Eq(Symbol("frf^2"), dyn_sys.frequency_response_function().doit() ** 2),
                marker=None,
            )
        )

        display(
            ReportText(
                f"""
                           Charakterystyka Amplitudowo-Częstotliwościowa podniesiona do kwadratu
                           """
            )
        )


# Marcel & Monika
class SteadySolutionComponent(ReportComponent):

    # "Rozwiązanie szczególne"
    title = "Steady solution"
    _phi = False

    @property
    def header_text(self):
        # "Rozwiązanie szczególne przedstawia wyrażenie:"
        # google tlumacz
        return "The steady solution is given by the formula:"

    @property
    def footer_text(self):
        # " Rozwiązanie szczególne związane jest obecnością wielkości wymuszających ruch (drgania) analizowanego układu."
        # google tlumacz
        return "The specific solution is related to the presence of quantities that force motion (vibrations) of the analyzed system."

    def append_elements(self, phi=_phi):

        from ....dynamics import HarmonicOscillator, LagrangesDynamicSystem

        system = self._system
        ReportText.set_directory("./SDAresults")

        latex_store = AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store

        t = system.ivar

        dyn_sys = system
        dyn_sys_lin = dyn_sys.linearized()

        display(ReportText(self.header_text))

        display(
            (
                SympyFormula(
                    Eq(
                        Symbol("X_s"),
                        HarmonicOscillator(dyn_sys_lin.linearized())
                        .steady_solution()
                        .n(3),
                        evaluate=False,
                    ),
                    marker="b",
                    backend=latex,
                )
            )
        )

        AutoBreak.latex_backend = latex_store

        display(ReportText(self.footer_text))


#### Amadi
class MaxStaticForce(ReportComponent):

    # Maksymalna siła statyczna
    title = "Maximum static force"

    @property
    def header_text(self):

        return "Wartość maksymalna siły statycznej działającej na pojedynczy element mocujący:"

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(
            SympyFormula(
                Eq(Symbol("F_s"), dyn_sys.max_static_force().doit()), marker=None
            )
        )


#### Amadi
class MaxDynamicForce(ReportComponent):

    # Maksymalna siła statyczna
    title = "Maximum dynamic force"

    @property
    def header_text(self):

        return "Wartość maksymalna siły dynamicznej działającej na pojedynczy element mocujący:"

    def append_elements(self):

        system = self._system
        dyn_sys = system
        dyn_sys_lin = dyn_sys

        display(ReportText(self.header_text))

        display(
            SympyFormula(
                Eq(Symbol("F_d"), dyn_sys.max_dynamic_force().doit()), marker=None
            )
        )


##### Karolina + Kamil + Bogi rr
class SeparableODEIntroComponent(ReportComponent):
    title = (
        "Podstawowa wiedza niezbędna do rozwiązania równania o zmiennych rozdzielonych"
    )

    def append_elements(self):
        system = self.reported_object
        ode_sys = system

        display(
            ReportText(
                "Jednym z fundamentalnych typów równań różniczkowych jest równanie o zmiennych rozdzielonych. Proces rozwiązania tego typu równania polega na separacji wyrazów zawierających $y$ i $x$, przeniesieniu ich na różne strony równania, a następnie zastosowaniu operacji całkowania. Po przeprowadzeniu całkowania nie można zapomnieć o dodaniu do jednej ze stron stałej całkowania. Poniżej znajduje się przykład zadania ilustrującego proces rozwiązania tego rodzaju równania.:"
            )
        )


class EquationDefinitionComponent(ReportComponent):
    title = "Rozpatrujemy następujące równanie różniczkowe: "

    def append_elements(self):

        system = self.reported_object
        ode_sys = system

        y_sym = system.dvars[0]
        display(
            SympyFormula(
                Eq(ode_sys.lhs[0], ode_sys.rhs[0]).subs(ode_sys._get_dvars_symbols())
            )
        )
        display(ReportText("Gdzie:"))
        display(SympyFormula(y_sym))
        display(ReportText("Jest naszą szukaną funkcją."))


class VariablesSeparationComponent(ReportComponent):

    # title="Separation of variables"
    title = "Rodzielanie zmiennych"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._get_dvars_symbols()[fun.diff(ivar)]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._get_dvars_symbols()[fun]
        gx = Symbol("g(x)")
        hx = Symbol("h(x)")

        # ode_rhs = solve(ode_sys.lhs[0],ode_sys.dvars[0].diff())[0]

        # sep_vars_dict = (separatevars(ode_rhs.subs(fun,y_sym),(ivar,y_sym),dict=True))

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        # display(ReportText('If it is possible variables $x$ and $y$ should be separated like shown below:'))
        display(
            ReportText(
                "Jeżeli jest to możliwe, dokonujemy powszechnie znanego rozdzielenia zmiennych $x$ i $y$. Poprzez rozdzielenie zmiennych rozumie się oddzielenie wyrazów zawierających zmienną $x$ od tych, które zawierają zmienną $y$. W celu dokonania tego podziału, możemy przyporządkować te wyrazy do funkcji $h(x)$ i $g(x)$, co pokazano poniżej:"
            )
        )
        display(SympyFormula(Eq(gx, g_fun_ivar)))
        display(SympyFormula(Eq(hx, h_fun_dvar)))
        display(
            ReportText(
                f"Poza tymi elementami nasze równanie zawiera również wyraz ${latex(fun_str)}$, który można również zapisać jako:"
            )
        )
        display(SympyFormula(Eq(fun_str, d_var / d_ivar)))
        display(
            ReportText(
                f"Po wstawieniu tego wyrażenia można zauważyć, że w pózniejszych krokach pomnożenie równania obustronnie przez ${latex(d_ivar)}$ umożliwia scałkowanie obu stron względem innych zmiennych."
            )
        )


class SeparatedVariablesIntegrationComponent(ReportComponent):

    # title="Integration of separated variables"
    title = "Całkowanie zmiennych rozdzielonych"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar
        gx = Symbol("g(x)")
        hx = Symbol("h(x)")

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]
        y_sym = system._variable_name()[2]
        fun_str = system._variable_name()[3]

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]
        # display(ReportText(f'We treat our ${latex(hx)}$ function as a ${latex(y_sym)}$'))
        display(
            ReportText("Poniższe równanie należy zmodyfikować:")
        )  ###nie jestem pewien co tutaj wsyawic jako ten opis

        display(
            SympyFormula(Eq(Symbol(latex(d_var / d_ivar)), g_fun_ivar * h_fun_dvar))
        )

        # display(ReportText('Two functions of seprarate variables are to define. These are as follows'  ))
        display(
            ReportText(
                f"W tym celu mnożymy obustronnie przez ${latex(d_ivar)}$ oraz przenosimy ${latex(fun)}$ na przeciwną stronę równania:"
            )
        )

        display(
            SympyFormula(
                Symbol(
                    latex(d_var / h_fun_dvar) + "=" + latex(g_fun_ivar) + latex(d_ivar)
                )
            )
        )
        # display(ReportText('Integration on both sides'))## uzupełnić opis
        display(
            ReportText(
                """Następnie, aby "wrócić" do funkcji pierwotnej oraz otrzymać wynik w postaci funkcji $y$ całkujemy otrzymane równanie obustronnie: """
            )
        )

        # display(SympyFormula(  Symbol('\\int ' + latex(d_var/y_sym) + '=  \\int ' + latex(g_fun_ivar) + latex(d_ivar) )   ) )
        display(
            ReportText(
                """Co ważne, przy całkowaniu należy pamiętać o dodaniu stałej całkowania w dowolnej postaci naogół zapisywanej jako $C$ bądź $D$. """
            )
        )
        display(
            SympyFormula(
                Symbol(
                    latex(Integral(1 / h_fun_dvar, y_sym))
                    + "= "
                    + latex(Integral(g_fun_ivar, ivar))
                )
            )
        )

        #         const = (system.solution._spot_constant())[0]
        const = Symbol("D")
        display(ReportText(" Po scałkowaniu otrzymujemy następujący wynik:"))
        display(
            SympyFormula(
                Symbol(
                    latex(Integral(1 / h_fun_dvar, y_sym).doit())
                    + "="
                    + latex(Integral(g_fun_ivar, ivar).doit())
                    + "+"
                    + latex(const)
                )
            )
        )


#         display(ReportText('Equations presented above correspond to each generalized coordinate:'))
#         display(SympyFormula(system.dvars))

#         display(ReportText(  self.footer_text   ))


class SolutionComponent(ReportComponent):

    title = "Uproszczenie rozwiązania otrzymanego po scałkowaniu równań"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        ############################ must have

        #         display(ReportText(  self.header_text   ))

        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]

        ode_rhs = solve(ode_sys.lhs[0], ode_sys.dvars[0].diff())[0]

        sep_vars_dict = separatevars(ode_rhs.subs(fun, y_sym), (ivar, y_sym), dict=True)

        g_fun_ivar = system._function_separation()[0]
        h_fun_dvar = system._function_separation()[1]

        # display(ReportText('After rearrangment a solution has the following form:'  ))
        display(
            ReportText(
                "Po uporządkowaniu oraz przeobrażeniu stałej wyrażenie ma następującą postać:"
            )
        )  ## trzeba rozwinąć hasło przeobrażenie
        display(SympyFormula(system.general_solution.as_eq_list()[0]))

        # display(ReportText('Obtained expresions enable to rewrite equation under consideration.'  ))
        display(
            ReportText(
                "Otrzymane wyrażenie stanowi rozwiązania równania różniczkowego o zmiennych rozdzielonych. Stała całkowania jest możliwa do policzenia w momencie w którym w zadaniu zostały podane warunki początkowe, czyli wartość funkcji dla konkretnego $x$."
            )
        )


class ClassifyODEComponent(ReportComponent):
    title = "Klasyfikacja równania różniczkowego"

    def append_elements(self):

        system = self.reported_object
        types = system._classify_ode()
        display(
            ReportText("To równanie różniczkowe możemy określić jako równanie typu:")
        )
        for feature in types:
            display(ReportText(f"\n\n - { feature } ".replace("_", " ")))


class LinearODEIntroComponent(ReportComponent):
    title = "Podstawowa wiedza niezbędna do rozwiązania równania różniczkowego"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar

        fun_str = system._variable_name()[3]

        d_ivar = system._variable_name()[0]
        d_var = system._variable_name()[1]

        y_sym = system._variable_name()[2]
        Px = Symbol("P(x)")
        Qx = Symbol("Q(x)")
        n = Symbol("n")
        display(
            ReportText(
                "Jednym z najbardziej podstawowych typów równań różniczkowych jest równanie liniowe o stałych współczynnikach. Równanie liniowe składa się jedynie z liniowych wyrazów, które mogą zawierać stałe współczynniki. Wymagane jest przejście na rówanie jednorodne- przenosimy wszystkie wyrazy z funkcją $y$ i jej pochodnymi, po czym przyrównujemy do zera otrzymując równanie jednorodne. Następnie wymagane jest uzmiennienie stałej i traktowanie jej jako funkcji zależnej od $x$. Liczymy pochodne obu stron i wykorzystujemy wyniki w początkowym równaniu. Wynikiem jest suma rozwiązań z obu części"
            )
        )  ###dodać: mające postać a*y'+b*y+R(x),gdzie r(x) to funckja x opisana jawnie gdzie a=... b=... r(x)=... i tutaj wpisać dla naszego przykładu z tego spot consant)


class LinearTransformation(ReportComponent):
    title = " Przejscie na równanie jednorodne"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        fun_str = system._variable_name()[3]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy

        y_sym = system._variable_name()[2]  # y(x)
        const = Symbol("D")

        display(
            ReportText(
                "Przenosimy wszystkie elementy z funkcją $y$ i jej pochodnymi na jedną stronę i przyrównujemy do zera, otrzymując równanie jednorodne."
            )
        )
        # display(SympyFormula(Eq(system._hom_equation().lhs[0],0)))


class LinearToSeparable(ReportComponent):

    title = "Rozwiązanie powstałego równania o zmiennych rozdzielonych"

    def append_elements(self):

        system = self.reported_object
        fun = system.dvars[0]

        ode_sys = system._hom_equation()
        # sep_ode_report=ode_sys.report #mniej więcej coś takiego tylko coś pomyliłem na 100%

        display(
            ReportText(
                "Uzyskane równanie jednorodne jest równaniem o zmiennych rozdzielonych i przyjmuje postać:"
            )
        )
        display(SympyFormula(Eq(ode_sys.lhs[0], 0)))

        # print("#"*200)
        sep_ode_report = ode_sys.report

        # print(list(sep_ode_report))
        # print('#'*3)
        # print(list(sep_ode_report)[3:])

        for elem in list(sep_ode_report)[3:]:
            self.append(elem)


class VariationOfConstant(ReportComponent):
    title = "Uzmiennienie stałej"

    def append_elements(self):

        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        y = fun
        ivar = system.ivar

        ode_sys2 = ode_sys
        x = ivar
        szur = ode_sys2.general_solution
        #         display(szur)
        ode_syso = Eq(ode_sys2.lhs[0], ode_sys2.rhs[0])
        Cx = Function(szur._spot_constant()[0])(x)
        Cx
        c_prim = Function("C'")(x)
        szur2 = szur.subs(szur._spot_constant()[0], Cx)
        # display(szur2)
        szur3 = Eq(szur2.lhs.diff(x), szur2.rhs.diff(x))
        # display(szur3)
        szur3_ = szur3.subs(Cx.diff(x), c_prim)
        szur3_
        # To implement, change left hand side to y_prim variable and subsitute variables with object variables

        ##### Tu sie teoretycznie zaczyna variation of constant ########################################
        sol_cx = Eq(Symbol("y"), szur2.rhs[0])
        display(
            ReportText(
                "Kolejnym krokiem jest potraktowanie stałej jako funkcji zależnej od $x$:"
            )
        )
        display(SympyFormula(sol_cx))
        sol_diff3 = Eq(Symbol("y'"), szur3_.rhs[0])
        display(ReportText("Następnie należy policzyć pochodne obu stron równania"))
        ## dodać równanie
        display(ReportText("Otrzymując równanie następującej postaci"))
        ##dodać spochodnione równanie chyba soldiff3 bedzie git

        display(SympyFormula(sol_diff3))  # gdzie? gdzie co?
        display(
            ReportText(
                "Dzięki tej operacji możemy teraz podstawić $y$ i $y'$ do oryginalnego równania"
            )
        )
        display(SympyFormula(ode_syso))
        lolo = ode_syso.subs(y.diff(ivar), sol_diff3.rhs).subs(y, sol_cx.rhs)
        display(ReportText("I otrzymamy:"))
        display(SympyFormula(lolo))
        c_p_sol = solve(lolo, c_prim)
        cp_eq = Eq(c_prim, c_p_sol[0])
        display(
            ReportText(
                "Następnie, jeżeli wszystko zostało policzne poprawnie to wartosci $C(x)$ powinny się skrócić i uprościć do postaci z której w łatwy sposób możemy wyznaczyć $C'(x)$ otrzymując:"
            )
        )
        display(SympyFormula(cp_eq))
        cx_wyn = Eq(Cx, cp_eq.rhs.diff(ivar))
        display(ReportText("Po obustronnym scałkowaniu otrzymamy wartosć stałej:"))
        display(SympyFormula(cx_wyn))
        FINALOWY_WYNIK = sol_cx.subs(Cx, cx_wyn.rhs)
        display(
            ReportText(
                "Aby otrzymać ostateczny wynik musimy podstawić otrzymaną stałą do pierwszego równania w którym uzmiennilismy stałą:"
            )
        )
        display(SympyFormula(FINALOWY_WYNIK))


class BernoulliODEIntroComponent(ReportComponent):
    title = "Podstawowa wiedza niezbędna do rozwiązania równania różniczkowego"

    def append_elements(self):

        system = self.reported_object
        y_sym = system.dvars[0]
        dane = system._get_dvars_symbols()
        Px = Symbol("P(x)")
        Qx = Symbol("Q(x)")
        n = Symbol("n")
        ynpow = Symbol("y") ** Symbol("n")
        display(
            ReportText(
                "Równanie różniczkowe typu Bernoulliego to szczególny rodzaj równania różniczkowego zwyczajnego pierwszego rzędu. Ma postać ogólną:"
            )
        )
        bernoulli_form = Eq(y_sym.diff() + Px * y_sym, Qx * y_sym**n)
        display(SympyFormula(bernoulli_form.subs(dane)))

        display(
            ReportText(
                f"Gdzie $y$ to nieznana funkcja zależna od zmiennej $x$, $P(x)$ i $Q(x)$ o dane funkcje zależne od $x$, a $n$ to stała liczba rzeczywista, różna od 0 i 1. Charakterystyczne dla równań Bernoulliego jest to, że są one nieliniowe z powodu obecności wyrazu ${latex(ynpow)}$ po jednej stronie równania. Rozwiązania tego typu równań mogą być trudne do uzyskania w ogólności, ale dla pewnych wartości $n$ można zastosować pewne techniki upraszczające rozwiązanie. W przypadku $n=0$ równanie Bernoulliego staje się liniowe, a dla $n=1$ można je przekształcić do postaci liniowego równania różniczkowego zwyczajnego pierwszego rzędu. W innych przypadkach konieczne może być zastosowanie specjalnych technik, takich jak zamiana zmiennych czy redukcja rzędu, aby uzyskać ogólne rozwiązanie."
            )
        )


class BernoulliTransformation(ReportComponent):
    title = " Wyznaczenie nowej funkcji z zależnej od x w celu przejścia na równanie liniowe"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar
        ###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy
        dane = system._get_dvars_symbols()
        bernoulli_sub = system._bernoulli_substitution()
        # bern_power = system

        n = Symbol("n")
        z = Function("z")(ivar)

        display(
            ReportText(
                "Aby dokonać podstawienia, konieczne jest zdefiniowanie $n$ poprzez znalezienie najwyższej potęgi $y$, gdzie $n$ to wykładnik tej funkcji. "
            )
        )
        display(SympyFormula(Eq(Symbol("n"), n.subs(bernoulli_sub))))
        display(
            ReportText(
                "Następnie należy dokonać podstawienia i przejść na nową funkcję $z$, używając poniższego wzoru:"
            )
        )
        display(SympyFormula(Eq(Symbol("z"), Symbol("y") ** (S.One - Symbol("n")))))
        display(ReportText("W naszym przypadku:"))
        display(
            SympyFormula(
                Eq(Symbol("z"), Symbol("y") ** (S.One - n.subs(bernoulli_sub)))
            )
        )
        # display(SympyFormula(Eq(Symbol('z'),z.subs(bernoulli_sub))))
        display(
            ReportText(
                "Następnym krokiem jest wyznaczenie pierwszej pochodnej funkcji $z$:"
            )
        )
        display(SympyFormula(Eq(Symbol("z'"), bernoulli_sub["z_diff"])))
        display(ReportText("Po podstawieniu wartosci n :"))
        display(
            SympyFormula(
                Eq(Symbol("z'"), bernoulli_sub["z_diff"].subs(n, bernoulli_sub[n]))
            )
        )

        display(
            ReportText("Następnie należy wyznaczyć z powyższych równań $y$ oraz $y'$")
        )
        display(SympyFormula(Eq(Symbol("y"), fun.subs(bernoulli_sub))))
        display(SympyFormula(Eq(Symbol("y'"), fun.diff().subs(bernoulli_sub))))


class BernoulliLinearTransformation(ReportComponent):
    title = " Przejście na równanie liniowe"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system
        fun = system.dvars[0]
        ivar = system.ivar
        ###########################DO ZMIANYYYYYYYY########################
        fun_str = system._variable_name()[2]  # y'

        d_ivar = system._variable_name()[0]  # dx
        d_var = system._variable_name()[1]  # dy
        dane = system._get_dvars_symbols()

        rownanie = system._as_reduced_ode()

        display(
            ReportText(
                "Ostatnim krokiem do przejscia na równanie liniowe jest podstawienie wyliczonych wyżej $y$ i $y'$ do pierwotnego równania:"
            )
        )
        display(SympyFormula(Eq(rownanie.lhs[0], rownanie.rhs[0])))
        display(ReportText("Otrzymane równanie liczymy w następujący sposób:"))

        sep_ode_report = rownanie.report

        for elem in list(sep_ode_report)[3:]:
            self.append(elem)


# Karolina i Karol done

general_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system.general_solution
"""


class ODEGeneralSolutionComponent(ReportComponent):

    title = "Wyznaczanie rozwiązania ogólnego"

    @property
    def import_text(self):

        return "Aby zaimportować przykładowy system  i wywołać rozwiązanie ogólne należy wywołać następujący kod:"

    @property
    def body_text(self):

        return "Rozwiązanie ogólne przyjmuje postać:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(general_code))

        display(ReportText(self.body_text))

        display((SympyFormula(system.general_solution)))


steady_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system.steady_solution
"""


class ODESteadySolutionComponent(ReportComponent):

    title = "Wyznaczanie rozwiązania ogólnego"

    @property
    def import_text(self):

        return "Aby zaimportować przykładowy system  i wywołać rozwiązanie szczególne należy wywołać następujący kod:"

    @property
    def body_text(self):

        return "Rozwiązanie szczególne przyjmuje postać:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(steady_code))

        display(ReportText(self.body_text))

        display((SympyFormula(system.steady_solution)))


import_code = """
from dynpy.models.odes.numerical import *

system=LinearFirstOrder.from_reference_data()
system.as_eq()  #zwraca równanie
system.as_matrix()  #zwraca macierz
system.as_iterable()  #zwraca listę
system.as_dict()  #zwraca słownik
system.as_eq_list()  #zwraca równanie w postaci listy
system._as_fode()  #zwraca równanie pierwszego rzędu

"""


class ODESystemRepresentationComponent(ReportComponent):

    title = "Różnorodne reprezentacje systemu"

    @property
    def import_text(self):

        return "Metody umożliwiające różnorodne reprezentacje systemu powinny być wywoływane w następujący sposób:"

    def append_elements(self):

        system = self.reported_object
        t = system.ivar
        dvars = system.dvars

        display(ReportText(self.import_text))
        display(ObjectCode(import_code))

        display(ReportText(r"Metoda .as_eq()"))
        display((system.as_eq()))

        display(ReportText(r"Metoda .as_matrix()"))
        display((system.as_matrix()))

        display(ReportText(r"Metoda .as_iterable()"))
        display((system.as_iterable()))

        display(ReportText(r"Metoda .as_dict()"))
        display((system.as_dict()))

        display(ReportText(r"Metoda .as_eq_list()"))
        display((system.as_eq_list()))

        display(ReportText(r"Metoda ._as_fode()"))
        display((system._as_fode()))


class ODEInitCodeComponent(ReportComponent):

    title = "Stworzenie liniowego równania różniczkowego"

    def append_elements(self):

        from dynpy.solvers.tools import CodePrinter

        system = self.reported_object
        ode_sys = system

        _dvars = system.dvars
        _odes = system.odes

        _ode_order = system.ode_order

        if _ode_order == 1:
            cls_name = "FirstOrderODESystem"
        else:
            cls_name = "ODESystem"

        string_gc = CodePrinter(_odes[0])._generate_code()

        display(
            ReportText(
                f"Stworzenie instacji klasy {cls_name} wymaga zdefiniowania następujących zmiennych oraz równań:"
            )
        )

        display(
            ObjectCode(
                f"""


from dynpy.solvers.linear import *
from dynpy.models.odes.linear import *
from sympy import *

{string_gc}

dvars = {_dvars}
odes = {_odes}
ode_order = {_ode_order}

odesys = {cls_name}(odes = odes , dvars = dvars, ode_order = ode_order)

display(odesys)

        """
            )
        )

        display(
            ReportText(
                f"Po wykonaniu instrukcji mamy gotową instancję obiektu klasy {cls_name}"
            )
        )


class ODESystemAdditionComponent(ODESystemOperations):

    title = "Adding two ODESystem Matrix"

    def append_elements(self):
        from sympy import Function, Matrix, Symbol

        systems = self.reported_object
        l = len(systems)
        matrix_add = self.get_zero_ode_sys()

        matrix_mat_lhs = Matrix([0])
        matrix_mat_rhs = Matrix([0])

        count = 1
        display(ReportText("#### Sprawdzanie sumy macierzy:"))
        for ode_sys in systems:
            display(ReportText(f"Macierz {count}"))
            display(ode_sys)
            count += 1
            matrix_add += ode_sys

            matrix_mat_lhs += ode_sys.lhs
            matrix_mat_rhs += ode_sys.rhs

        display(ReportText("Suma macierzy ODEsystem"))
        display(matrix_add)
        display(ReportText("Suma macierzy Sympy"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemSubtractionComponent(ODESystemOperations):

    title = "Subtractiing two ODESystem Matrix"

    def append_elements(self):
        systems = self.reported_object
        l = len(systems)
        matrix_sub = systems[0]

        matrix_mat_lhs = systems[0].lhs
        matrix_mat_rhs = systems[0].rhs

        display(ReportText("#### Sprawdzanie różnicy macierzy:"))
        display(ReportText(f"Macierz 1"))
        display(systems[0])
        count = 2
        for ode_sys in systems[1:]:
            display(ReportText(f"Macierz {count}"))
            display(ode_sys)
            count += 1
            matrix_sub -= ode_sys

            matrix_mat_lhs -= ode_sys.lhs
            matrix_mat_rhs -= ode_sys.rhs

        display(ReportText("Różnica macierzy ODEsystem"))
        display(matrix_sub)
        display(ReportText("Różnica macierzy Sympy"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemMultiplicationComponent(ODESystemOperations):

    title = "Multiplying two ODESystem Matrix"

    def append_elements(self):
        systems = self.reported_object
        l = len(systems)
        matrix_mul = 1
        count = 1

        matrix_mat_lhs = Matrix([1])
        matrix_mat_rhs = Matrix([1])

        display(ReportText("#### Sprawdzanie iloczynu macierzy:"))
        for ode_sys in systems:
            display(ReportText(f"Macierz {count}"))
            display(ode_sys)
            count += 1
            matrix_mul *= ode_sys

            matrix_mat_lhs *= ode_sys.lhs
            matrix_mat_rhs *= ode_sys.rhs

        display(ReportText("Iloczyn macierzy ODEsystem"))
        display(matrix_mul)
        display(ReportText("Iloczyn macierzy Sympy"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


class ODESystemExponentiationComponent(ReportComponent):
    title = "Exponentiation of ODESystem Matrix"

    def append_elements(self):
        system = self.reported_object
        ode_sys = system

        display(ReportText("#### Sprawdzanie potęgi macierzy:"))
        display(ReportText("Macierz"))
        display(ode_sys)

        display(ReportText("Potęga macierzy ODEsystem"))

        matrix_exp = ode_sys * ode_sys
        display(matrix_exp)

        matrix_mat_lhs = ode_sys.lhs * ode_sys.lhs
        matrix_mat_rhs = ode_sys.rhs * ode_sys.rhs

        display(ReportText("Potęga macierzy Sympy"))
        display(Eq(matrix_mat_lhs, matrix_mat_rhs, evaluate=False))


excitation_code = """
from dynpy.models.odes.numerical import *
system=LinearFirstOrder.from_reference_data()
system._as_fode()._get_excitation_comps
"""


class PredictionOfSteadySolutionComponent(ReportComponent):

    title = "Wyznaczanie przewidywanego rozwiązania szczególnego równania"

    @property
    def import_text(self):

        return "Aby zaimportować przykładowy system  i wyznaczyć jego przewidywane rozwiązanie szczególne należy wywołać następujący kod:"

    @property
    def body_text(self):

        return "Rozwiązanie przyjmuje postać:"

    def append_elements(self):
        system = self.reported_object
        t = system.ivar
        dvars = system.dvars
        display(ReportText(self.import_text))
        display(ObjectCode(excitation_code))
        display(ReportText(self.body_text))
        system._as_fode()._get_excitation_comps

        ext_comps_list = system._as_fode()._get_excitation_comps

        for comp, amp in ext_comps_list:

            display(SympyFormula(comp))
            display(SympyFormula(amp))


class HomoPredictionIntroComponent(ReportComponent):

    title = "Równanie różniczkowe drugiego rzędu - metoda przewidywań"

    def append_elements(self):
        system = self.reported_object
        dvar = system.dvars[0]
        ivar = system.ivar
        display(
            ReportText(
                "Jednym z podstawowych typów równań różniczkowych drugiego rzędu jest równanie jednorodne o stałych współczynnikach . Proces rozwiązania tego typu równania polega na przewidzeniu rozwiązania, które przyjmuje ogólną postać:"
            )
        )
        display(SympyFormula(Eq(dvar, Symbol("C") * exp(Symbol("r") * ivar))))
        display(
            ReportText(
                "Na podstawie otrzymanych pierwiastków (rozwiązań) przewidywanego równania dobiera się docelowe równanie ogólne. Poniżej znajduje się przykład zadania ilustrującego proces rozwiązania tego rodzaju równania:"
            )
        )


class MainPredictionComponent(ReportComponent):

    title = "Przewidywanie rozwiązania równania jednorodnego"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        C = Symbol("C", positive=True)
        r = Symbol("r")
        y = system.dvars[0]
        x = system.ivar
        sym_dict = system._get_dvars_symbols()
        sys_sub = Eq(system.lhs[0], system.rhs[0])

        prdct = Eq(y, C * exp(r * x))
        prdct_diff = Eq(prdct.lhs.diff(x), prdct.rhs.diff(x))
        prdct_ddiff = Eq(prdct.lhs.diff(x, 2), prdct.rhs.diff(x, 2))

        prdct_dict = {
            y: prdct.rhs,
            y.diff(x): prdct_diff.rhs,
            y.diff(x, 2): prdct_ddiff.rhs,
        }

        r_eq = Eq(system.fundamental_matrix()[0], 0)
        r0_eq = Eq(y.subs(sym_dict), Symbol("r^0"))
        r1_eq = Eq(y.diff(x).subs(sym_dict), Symbol("r^1"))
        r2_eq = Eq(y.diff(x, 2).subs(sym_dict), Symbol("r^2"))

        display(
            ReportText(
                "Pierwszym krokiem do rozwiązania danego równania jest wyznaczenie pochodnych przewidywanego rozwiązania. Przypominając pierwotne równanie które przewidujemy jako nasze rozwiązanie ma następującą postać:"
            )
        )
        display(SympyFormula(prdct.subs(sym_dict)))
        display(
            ReportText(
                "A następnie znalezienie potrzebnych pochodnych. W przypadku równań drugiego rzędu jak nazwa wskazuje najwyższa pochodna jest drugiego stopnia stąd potrzebne nam dwie pierwsze pochodne przewidywanego rozwiązania. Pierwsza pochodna:"
            )
        )
        display(SympyFormula(prdct_diff.subs(sym_dict)))
        display(ReportText("Druga pochodna:"))
        display(SympyFormula(prdct_ddiff.subs(sym_dict)))
        display(
            ReportText(
                "Następnym krokiem jest podstawienie poszczególnych przewidywań do odpowiednich pochodnych w pierwotnym równaniu:"
            )
        )
        display(SympyFormula(sys_sub.subs(prdct_dict)))
        display(
            ReportText(
                "Można zauważyć powtarzające się elementy i wyłączyć je przed nawias:"
            )
        )
        display(SympyFormula(sys_sub.subs(prdct_dict).simplify()))
        display(
            ReportText(
                "Podstawowa zasada matematyki mówi, że jeśli iloczyn dwóch czynników jest równy zero, to przynajmniej jeden z tych czynników musi być równy zero. W metodzie przewidywań nie rozpatruje się jednak dwóch przypadków, ponieważ gdy stała całkowania, a zatem funkcja eksponenty, jest równa zeru to zadanie staje się trywialne i nie wymaga dalszej analizy. Dlatego skupiamy się na drugim czynniku, powszechnie nazywanym równaniem charakterystycznym, które dostarczy nam pierwiastków niezbędnych do przewidzenia ostatecznego równania. Równanie poddawane dalszej analizie ma następująca postać:"
            )
        )
        display(SympyFormula(r_eq))
        display(
            ReportText(
                "Po rozwiązaniu pewnej liczby zadań można zauważyć pewną zależność wynikającą z powyższego przewidywania, co pozwoli na pominięcie tego kroku i podstawiania odpowiednich wartości zgodnie z poniższymi podstawieniami:"
            )
        )
        display(SympyFormula(r0_eq))
        display(SympyFormula(r1_eq))
        display(SympyFormula(r2_eq))


class RootsAnalysisComponent(ReportComponent):

    title = "Rozwiązanie równania charakterystycznego i ogólna analiza jego rozwiązania w kontekście rozwiązywania równań różniczkowych drugiego rzędu"

    def append_elements(self):

        system = self.reported_object._hom_equation()
        sym_dict = system._get_dvars_symbols()
        y = system.dvars[0].subs(sym_dict)
        x = system.ivar
        char = system.fundamental_matrix()[0]
        char_eq = Eq(char, 0)
        roots = system._eigenvalues()
        C1 = Symbol("C_1", positive=True)
        C2 = Symbol("C_2", positive=True)
        r0 = Symbol("r_0")
        r1 = Symbol("r_1")
        a = Symbol("alpha")
        b = Symbol("beta")
        i = Symbol("i")

        r1_img = Eq(r0, a + i * b)
        r2_img = Eq(r1, a - i * b)

        pred_real = Eq(y, C1 * exp(r0 * x) + C2 * exp(r1 * x))
        pred_one = Eq(y, C1 * exp(r0 * x) + C2 * x * exp(r0 * x))
        pred_imagine = Eq(y, exp(a * x) * (C1 * cos(b * x) + C2 * sin(b * x)))

        sol_eq = Eq(y, system._ode_solution().rhs[0])

        display(
            ReportText(
                f"W równaniach różniczkowych drugiego rzędu równanie charakterystyczne przyjmuje postać równania kwadratowego. Takie równanie może w rozwiązaniu dać trzy typy pierwiastków i w zależności od danego typu dobiera się rozwiązanie ogólne tego równania. Gdy policzona delta jest większa od zera, otrzymujemy dwa pierwiastki rzeczywiste o postaci ${latex(r0)}$ i ${latex(r1)}$, a wtedy przewidywane rozwiązanie przyjmuje postać:"
            )
        )
        display(SympyFormula(pred_real))
        display(
            ReportText(
                f"Gdy wyliczona delta jest równa 0, otrzymujemy jeden pierwiastek o postaci ${latex(r0)}$, a wtedy przewidywane rozwiązanie przyjmuje postać:"
            )
        )
        display(SympyFormula(pred_one))
        display(
            ReportText(
                "Gdy wyliczona delta jest mniejsza od 0, otrzymujemy dwa pierwiastki zespolone o postaci:"
            )
        )
        display(SympyFormula(r1_img))
        display(SympyFormula(r2_img))
        display(
            ReportText(
                f"Gdzie ${latex(a)}$ to część rzeczywista pierwiastka, a ${latex(b)}$ to część urojona. Przewidywane rozwiązanie w takim wypadku przyjmuje następującą postać:"
            )
        )
        display(SympyFormula(pred_imagine))
        display(ReportText("Równanie charakterystyczne w danym przykładzie ma postać:"))
        display(SympyFormula(char_eq))
        display(ReportText("Pierwiastki tego równania przyjmują postać:"))
        i = 0
        lista = system._roots_list()
        for elem in lista:
            if elem != 0:
                display(SympyFormula(Eq(Symbol(f"r_{i}"), elem), marker=None))
                i += 1
        display(
            ReportText(
                "Po rozpoznaniu typu pierwiastków możemy zapisać rozwiązanie jako:"
            )
        )
        display(SympyFormula(sol_eq))
