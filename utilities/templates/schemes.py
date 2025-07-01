from ..report import TikZPicture


class TuningTikZDiagram(TikZPicture):

    def _scheme_desc(self):

        code = r"""
\coordinate (origo) at (0,0);
\node (B1) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, rounded corners=0.5cm] at (origo) {Begin tuning};
\node (B2) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B1.south) {Pump out the \\liquid from TMD to \\the skyscraper.};
\draw [thick, ->] ([yshift=0cm]B1.south) -- ([xshift=0cm,yshift=0cm]B2.north) node[midway,above] {};
\node (B3) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-2cm]B2.south) {Examine the \\natural frequency \\of the dynamic \\mass of the \\skyscraper.};
\draw [thick, ->] ([yshift=0cm]B2.south) -- ([xshift=0cm,yshift=0cm]B3.north) node[midway,above] {};
\node (B4) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-2.5cm]B3.south) {Calculate and \\transfer the \\amount of liquid \\that provides \\lower amplitudes \\of the damped \\structure.};
\draw [thick, ->] ([yshift=0cm]B3.south) -- ([xshift=0cm,yshift=0cm]B4.north) node[midway,above] {};
\node (B5) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, diamond, aspect=1.5] at ([xshift=0cm,yshift=-5.5cm]B4.south) {Compare the data with the model, \\is value of amplitude of linear acceleration \\of TMD equal or higher than in model?};
\draw [thick, ->] ([yshift=0cm]B4.south) -- ([xshift=0cm,yshift=0cm]B5.north) node[midway,above] {};
\node (B6) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4.2cm,yshift=0cm]B5.east) {TMD is tuned properly.};
\node (B7) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=-4.2cm,yshift=0cm]B5.west) {Retune the TMD.};
\draw [thick, ->] ([xshift=0cm]B5.west) -- ([xshift=0cm,yshift=0cm]B7.east) node[midway,above] {Yes};
\draw [thick, ->] ([xshift=0cm]B5.east) -- ([xshift=0cm,yshift=0cm]B6.west) node[midway,above] {No};
\node (B8) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, rounded corners=0.5cm] at ([xshift=0cm,yshift=5cm]B7.north) {End tuning};
\draw [thick, ->] ([xshift=0cm]B7.north) -- ([xshift=0cm,yshift=0cm]B8.south) node[midway,above] {};
\draw [thick, ->] ([xshift=0cm]B6.north) |- ([xshift=0cm,yshift=0cm]B3.east) node[midway,above] {};
"""

        return code


class DiagnosticsTikZDiagram(TikZPicture):

    def _scheme_desc(self):

        code = r"""
\coordinate (origo) at (0,0);
\node (B1) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, rounded corners=0.5cm] at (origo) {Begin diagnostics};
\node (B2) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B1.south) {Estimate the TMD's mass.};
\draw [thick, ->] ([yshift=0cm]B1.south) -- ([xshift=0cm,yshift=0cm]B2.north) node[midway,above] {};
\node (B3) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-3cm]B2.south) {Gain data from \\the linear \\acceleration sensor \\on the TMD};
\draw [thick, ->] ([yshift=0cm]B2.south) -- ([xshift=0cm,yshift=0cm]B3.north) node[midway,above] {};
\node (B4) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, diamond, aspect=1.5] at ([xshift=0cm,yshift=-3cm]B3.south) {$a_p {real} \geqslant  a_p {mod}$};
\draw [thick, ->] ([yshift=0cm]B3.south) -- ([xshift=0cm,yshift=0cm]B4.north) node[midway,above] {};
\node (B5) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4.2cm,yshift=0cm]B4.east) {Retune the TMD. \\Structure health \\is compromised.};
\node (B6) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=-4.2cm,yshift=0cm]B4.west) {TMD is tuned properly. \\Structure health \\is not compromised.};
\draw [thick, ->] ([xshift=0cm]B4.west) -- ([xshift=0cm,yshift=0cm]B6.east) node[midway,above] {Yes};
\draw [thick, ->] ([xshift=0cm]B4.east) -- ([xshift=0cm,yshift=0cm]B5.west) node[midway,above] {No};
\node (B7) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=5cm]B6.north) {Wait for another \\diagnostic interval.};
\draw [thick, ->] ([xshift=0cm]B6.north) -- ([xshift=0cm,yshift=0cm]B7.south) node[midway,above] {};
\node (B8) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, rounded corners=0.5cm] at ([xshift=0cm,yshift=5cm]B5.north) {End diagnostics};
\draw [thick, ->] ([xshift=0cm]B5.north) -- ([xshift=0cm,yshift=0cm]B8.south) node[midway,above] {};
\draw [thick, ->] ([xshift=0cm]B7.north) |- ([xshift=0cm,yshift=0cm]B1.west) node[midway,above] {};
"""

        return code


class TuningTikZDiagram1(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    
    

    \draw
    (0,0) to[american voltage source, l_=$U_{oc}$] (0,-3)
    (0,0) to[R,european,-*,l=$R_{0}$] (6,0) node[label={right:$+$}] {}
    (0,-3) to[short, -*] (6,-3) node[label={right:$-$}] {}
;

"""

        return code


class TheveninTikZDiagramVer1(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    
    \draw
    (0,0) to[american voltage source, l_=$U_{oc}$] (0,-2.5)
    (0,0) to[R,european, l=$R_{0}$] (3,0)
    
    (3,0) -- (3,1.5)
    to[C, l=$C_{1}$] (6,1.5) -- (6,0)
    (3,0) to[R,european, l=$R_{1}$] (6,0)
    
    (6,0) to[short, -*] (7,0) node[label={right:$+$}] {}
    
    (0,-2.5) to[short, -*] (7,-2.5) node[label={right:$-$}] {}
;

"""

        return code


class TheveninTikZDiagramVer2(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    \draw
    (0,0) to[american voltage source, l_=$U_{oc}$] (0,-3)
    (0,0) to[R,european, l=$R_{0}$] (3,0)
    
    (3,0) -- (3,1.5)
    to[C, l=$C_{1}$] (6,1.5) -- (6,0)
    (3,0) to[R,european, l=$R_{1}$] (6,0)
    
    (6,0) -- (7,0)
    
    (7,0) -- (7,1.5)
    to[C, l=$C_{2}$] (10,1.5) -- (10,0)
    (7,0) to[R,european, l=$R_{2}$] (10,0)
    
    
    
    (10,0) to[short, -*] (11,0) node[label={right:$+$}] {}
    
    (0,-3) to[short, -*] (11,-3) node[label={right:$-$}] {}
;

"""

        return code


class TheveninTikZDiagram(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    
    \draw
    (0,0) to[american voltage source, l_=$U_{oc}$] (0,-3)
    (0,0) to[R,european, l=$R_{0}$] (3,0)
    
    (3,0) to[short, i_=$I_{TH}$] (3,1.5)
    to[C, l=$C_{1}$] (6,1.5) -- (6,0)
    (3,0) to[R,european, l=$R_{1}$,a= + $ U_{TH} $ -] (6,0)
    
    (6,0) to[short, -*, i_=$I_{L}$] (7,0) node[label={right:$+$}] {}
    
    (0,-3) to[short, -*] (7,-3) node[label={right:$-$}] {}
;
"""

        return code


class TikZDiagram(TikZPicture):

    def _scheme_desc(self):

        code = r"""
\coordinate (origo) at (0,0);
\node (B1) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, rounded corners=0.5cm] at (origo) {Podaj funkcję limitu stanu, \\wartości parametrów, cenę, \\oczekiwany indeks niezawodności, \\funkcję kosztów zmiany parametrów};
\node (B2) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B1.south) {Wykonaj obliczenia przy użyciu \\metody odwrotnej niezawodności};
\draw [thick, ->] ([yshift=0cm]B1.south) -- ([xshift=0cm,yshift=0cm]B2.north) node[midway,above] {};
\node (B3) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B2.south) {Zbadaj wrażliwość parametrów};
\draw [thick, ->] ([yshift=0cm]B2.south) -- ([xshift=0cm,yshift=0cm]B3.north) node[midway,above] {};
\node (B4) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B3.south) {Na podstawie funkcji kosztów określ \\cenę dla aktualnej tolerancji wykonania};
\draw [thick, ->] ([yshift=0cm]B3.south) -- ([xshift=0cm,yshift=0cm]B4.north) node[midway,above] {};
\node (B5) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center, diamond, aspect=1.5] at ([xshift=0cm,yshift=-3cm]B4.south) {Porównaj aktualną cenę \\z ceną poprzednią};
\draw [thick, ->] ([yshift=0cm]B4.south) -- ([xshift=0cm,yshift=0cm]B5.north) node[midway,above] {};
\node (B6) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4.5cm,yshift=0cm]B5.east) {ALGORYTM \\OPTYMALIZACYJNY \\dobór nowych parametrów};
\node (B7) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=0cm,yshift=-1.5cm]B5.south) {Zapisz wartości};
\node (B8) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center,, rounded corners=0.5cm] at ([xshift=-4.2cm,yshift=0cm]B5.west) {Podaj wartości \\parametrów i zakończ};
\draw [thick, ->] ([xshift=0cm]B5.west) -- ([xshift=0cm,yshift=0cm]B8.east) node[midway,above] {Równa};
\draw [thick, ->] ([xshift=0cm]B5.south) -- ([xshift=0cm,yshift=0cm]B7.north) node[midway,right] {Niższa cena};
\draw [thick, ->] ([xshift=0cm]B5.east) -- ([xshift=0cm,yshift=0cm]B6.west) node[midway,above] {Wyższa cena};
\draw [thick, ->] ([xshift=0cm]B6.north) |- ([xshift=0cm,yshift=0cm]B2.east) node[midway,above] {};
\draw [thick, ->] ([xshift=0cm]B7.east) -| ([xshift=0cm,yshift=0cm]B6.south) node[midway,above] {};
"""

        return code
