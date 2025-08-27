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


class ForcedDampedSpringMassSystemScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

\coordinate (origo) at (0,0);

\tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.3cm,post length=0.3cm,segment length=0.3cm}]

\tikzstyle{damper}=[thick,decoration={markings,
mark connection node=dmp,
mark=at position 0.5 with
{
\node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
\draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
\draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
}
}, decorate]

\tikzstyle{ground}=[fill,pattern=north east lines,draw=none,minimum width=0.75cm,minimum height=0.3cm]



%draw axes
%\fill[black] (origo) circle (0.05);


\node (M) [draw,outer sep=0pt,thick,minimum width=3cm, minimum height=1.5cm,yshift=2cm] at (origo) {$m$};
%\node (m1_label) at ($(M)+(0:0.5cm)$) {$m$};

\draw [spring]    ([xshift=-1cm]M.south) --++(0,-3cm) node (k_tire_end){} node[midway,left] {$k$};

    \draw [damper] ([xshift=1cm]M.south) --++(0,-3cm) node (c_tire_end){}  node[midway,right=0.25cm] {$c$};

%\fill [black] (M.center) circle (0.05);

% \node (beam) [fill=gray,anchor=north,xshift=0cm,yshift=0cm,minimum width=2cm,minimum height=0.05pt] at (origo) {} node[right]{$m_{beam}=0$};

% \draw [ultra thick] (k_sup_spring_left.west) -- (k_sup_spring_right.east) node[midway] (beam_center) {};


%         \node (road_sur) [fill=black,anchor=north,xshift=0cm,yshift=-2cm,minimum width=2cm,minimum height=0.05cm] at (beam) {};

%beam_center.center

%    \draw[spring] ([xshift=-0.0cm]k_sup_spring_left.center) -- ++(0,-1.5cm) node (k_tire_end) {} node[midway,left] {$k_{tire}$};

%        \draw[damper] ([xshift=+0.5cm](k_sup_spring_right.center) -- ++(0,-2cm) node (c_tire_end) {} node[midway,right=0.25cm] {$c_{tire}$};

\draw[ultra thick] ([xshift=-0.0cm]k_tire_end.west) -- ([xshift=+0.0cm]c_tire_end.east) node (force_attachment_point) {};


\draw [thin] (force_attachment_point) -- +(1.5,0cm) coordinate (force_leader_edn);
\draw [-latex,ultra thick] (force_leader_edn) -- +(0,1.5cm) node[above right] (force_t)  {$F(t)=F \cos(\Omega t)$};

"""

        return code


class ForcedSpringMassSystemScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

\coordinate (origo) at (0,0);

\tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.3cm,post length=0.3cm,segment length=0.3cm}]

\tikzstyle{damper}=[thick,decoration={markings,
mark connection node=dmp,
mark=at position 0.5 with
{
\node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
\draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
\draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
}
}, decorate]

\tikzstyle{ground}=[fill,pattern=north east lines,draw=none,minimum width=0.75cm,minimum height=0.3cm]



%draw axes
%\fill[black] (origo) circle (0.05);


\node (M) [draw,outer sep=0pt,thick,minimum width=3cm, minimum height=1.5cm,yshift=2cm] at (origo) {$m$};
%\node (m1_label) at ($(M)+(0:0.5cm)$) {$m$};

\draw [spring]    ([xshift=0cm]M.south) --++(0,-3cm) node (k_tire_end){} node[midway,left] {$k$};

%    \draw [damper] ([xshift=1cm]M.south) --++(0,-3cm) node (c_tire_end){}  node[midway,right=0.25cm] {$c$};

%\fill [black] (M.center) circle (0.05);

% \node (beam) [fill=gray,anchor=north,xshift=0cm,yshift=0cm,minimum width=2cm,minimum height=0.05pt] at (origo) {} node[right]{$m_{beam}=0$};

% \draw [ultra thick] (k_sup_spring_left.west) -- (k_sup_spring_right.east) node[midway] (beam_center) {};


%         \node (road_sur) [fill=black,anchor=north,xshift=0cm,yshift=-2cm,minimum width=2cm,minimum height=0.05cm] at (beam) {};

%beam_center.center

%    \draw[spring] ([xshift=-0.0cm]k_sup_spring_left.center) -- ++(0,-1.5cm) node (k_tire_end) {} node[midway,left] {$k_{tire}$};

%        \draw[damper] ([xshift=+0.5cm](k_sup_spring_right.center) -- ++(0,-2cm) node (c_tire_end) {} node[midway,right=0.25cm] {$c_{tire}$};

\draw[ultra thick] ([xshift=-0.5cm]k_tire_end.west) -- ([xshift=+0.5cm]k_tire_end.east) node (force_attachment_point) {};


\draw [thin] (force_attachment_point) -- +(1.5,0cm) coordinate (force_leader_edn);
\draw [-latex,ultra thick] (force_leader_edn) -- +(0,1.5cm) node[above right] (force_t)  {$F(t)=F \cos(\Omega t)$};

"""

        return code


class GearModelOscillatorScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

\coordinate (origo) at (0,0);

\tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.3cm,post length=0.3cm,segment length=0.3cm}]

\tikzstyle{damper}=[thick,decoration={markings,
mark connection node=dmp,
mark=at position 0.5 with
{
\node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
\draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
\draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
}
}, decorate]

\tikzstyle{ground}=[fill,pattern=north east lines,draw=none,minimum width=0.75cm,minimum height=0.3cm]



%draw axes
%\fill[black] (origo) circle (0.05);


\node (M) [draw,outer sep=0pt,thick,minimum width=3cm, minimum height=1.5cm,yshift=2cm] at (origo) {$m_{red}$};
%\node (m1_label) at ($(M)+(0:0.5cm)$) {$m$};

\draw [spring]    ([xshift=-1cm]M.south) --++(0,-3cm) node (k_tire_end){} node[midway,left] {$k_{m}$};

    \draw [damper] ([xshift=1cm]M.south) --++(0,-3cm) node (c_tire_end){}  node[midway,right=0.25cm] {$c_{m}$};

%\fill [black] (M.center) circle (0.05);

% \node (beam) [fill=gray,anchor=north,xshift=0cm,yshift=0cm,minimum width=2cm,minimum height=0.05pt] at (origo) {} node[right]{$m_{beam}=0$};

% \draw [ultra thick] (k_sup_spring_left.west) -- (k_sup_spring_right.east) node[midway] (beam_center) {};


%         \node (road_sur) [fill=black,anchor=north,xshift=0cm,yshift=-2cm,minimum width=2cm,minimum height=0.05cm] at (beam) {};

%beam_center.center

%    \draw[spring] ([xshift=-0.0cm]k_sup_spring_left.center) -- ++(0,-1.5cm) node (k_tire_end) {} node[midway,left] {$k_{m}$};

%        \draw[damper] ([xshift=+0.5cm](k_sup_spring_right.center) -- ++(0,-2cm) node (c_tire_end) {} node[midway,right=0.25cm] {$c_{m}$};

\draw[ultra thick] ([xshift=-0.0cm]k_tire_end.west) -- ([xshift=+0.0cm]c_tire_end.east) node (force_attachment_point) {};


\draw [thin] (force_attachment_point) -- +(0,0) coordinate (force_leader_edn);
\draw [-latex,ultra thick] (force_leader_edn) -- +(0,-1.5cm) node[above right] (force_t)  {$F$};

\draw [thin] (0,-2.8) -- +(0,1) coordinate (force_leader_edn) {} node[midway,left] {$x$};
\draw[->]        (0,-2)   -- (0,-2.8);

"""

        return code


class HelicalGearModelScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    \tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.2cm,post length=0.2cm,segment length=0.2cm}]

    \tikzstyle{damper}=[thick,decoration={markings,
    mark connection node=dmp,
    mark=at position 0.5 with
    {
    \node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
    \draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
    \draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
    }
    }, decorate]

    % Draw the first circle with radius 1
    \draw[thick] (0,0) circle (1);
    \node at (0,0) {$O_{1}$};
    \node at (0, -1.3) {$J_{1}$};
    \node at (0, -1.6) {$r_{b1}$};

    % Draw the second circle with radius 2
    \draw [thick] (6,0) circle (2);
    \node at (6,0) {$O_{2}$};
    \node at (6, -2.3) {$J_{2}$};
    \node at (6, -2.6) {$r_{b2}$};

    % Draw the curved arrow for the first circle
    \draw[->, red, thick] (-1.3,0) arc[start angle=-10, end angle=45, radius=-1];
    \node[red] at (-1.2,-1.2) {$T_1$};

    % Draw the curved arrow for the second circle
    \draw[->, red, thick] (8.5,0) arc[start angle=0, end angle=45, radius=2];
    \node[red] at (8.5,1) {$T_2$};

    %Drawing the phi1 coordinates arrows
    \draw[<-, thick] (-1.3,0) arc[start angle=-10, end angle=-45, radius=-1];
    \node  at (-1.6,0.6) {$\phi_{1}(t)$};

    %Drawing the phi2 coordinate arrow
    \draw[->, thick] (8.5,0) arc[start angle=0, end angle=-45, radius=2];
    \node  at (8.9,-0.6) {$\phi_{2}(t)$};

    %Drawing the spring and damper connection
   \draw [thick] (0.28,-0.97) -- (2.3,0.17);
   \draw [thick] (4.8,1.6) -- (3.37,0.78);
   \draw [spring]  (2.02,0.571) -- (3.06,1.212);
   \draw [damper]  (2.705,-0.383) -- (3.745,0.293);
   \draw (2.02,0.571) -- (2.705,-0.383);
   \draw (3.06,1.212) -- (3.745,0.293);
   \node at (2.5,1.4) {$k_{m}(t)$};
   \node at (3.5,-0.6) {$c_{m}(t)$};


    """
        return code


class HelicalGearVerticalModelScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    \tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.2cm,post length=0.2cm,segment length=0.2cm}]

    \tikzstyle{damper}=[thick,decoration={markings,
    mark connection node=dmp,
    mark=at position 0.5 with
    {
    \node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
    \draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
    \draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
    }
    }, decorate]

    % Draw the first circle with radius 1
    \draw[thick] (0,0) circle (1);
    \node at (0,0) {$O_{1}$};
    \node at (0, -1.3) {$J_{1}$};
    \node at (0, -1.6) {$r_{b1}$};

    % Draw the second circle with radius 2
    \draw [thick] (0,6) circle (2);
    \node at (0,6) {$O_{2}$};
    \node at (0, -2.3) {$J_{2}$};
    \node at (0, -2.6) {$r_{b2}$};

    % Draw the curved arrow for the first circle
    \draw[->, red, thick] (-1.3,0) arc[start angle=-10, end angle=45, radius=-1];
    \node[red] at (-1.2,-1.2) {$T_1$};

    % Draw the curved arrow for the second circle
    \draw[->, red, thick] (0,8.5) arc[start angle=0, end angle=45, radius=2];
    \node[red] at (1,8.5) {$T_2$};

    %Drawing the phi1 coordinates arrows
    \draw[<-, thick] (-1.3,0) arc[start angle=-10, end angle=-45, radius=-1];
    \node  at (-1.6,0.6) {$\phi_{1}(t)$};

    %Drawing the phi2 coordinate arrow
    \draw[->, thick] (1,8.5) arc[start angle=0, end angle=-45, radius=2];
    \node  at (-0.6,8.5) {$\phi_{2}(t)$};

    %Drawing the spring and damper connection
   %\draw [thick] (0.28,-0.97) -- (2.3,0.17);
   %\draw [thick] (4.8,1.6) -- (3.37,0.78);
   %\draw [spring]  (2.02,0.571) -- (3.06,1.212);
   %\draw [damper]  (2.705,-0.383) -- (3.745,0.293);
   %\draw (2.02,0.571) -- (2.705,-0.383);
   %\draw (3.06,1.212) -- (3.745,0.293);
   %\node at (2.5,1.4) {$k_{m}(t)$};
   %\node at (3.5,-0.6) {$c_{m}(t)$};


    """
        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class CompoundPendulumScheme2(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    % Osie
    \draw[->] (-0.2,0) -- (4,0) node[right] {$x$};
    \draw[->] (0,0.2) -- (0,-4) node[below] {$y$};
    \node at (0.3,-1.7) {$\varphi$};


    % Wahadło obrócone
    \begin{scope}[rotate=25]
        \filldraw[fill=black!10, draw=black, thick] (0,0) -- (1.5,-3) arc[start angle=30, end angle=-190, radius=0.8] -- cycle;

    % Oś środka ciężkości
         \draw[dashed] (0,0) -- (1,-5);
         \filldraw[black] (0.5,-2.5) circle (1.5pt);

    % Oznaczenia
    \node at (0.8,-2.5) {$C$};
    \node at (2,-4.5) {$m, I, l$};
    %Moment wymuszający

    \end{scope}

    % Strzałka
    \draw[->] (0,-1.5) arc[start angle=-80,end angle=-50,radius=1.5];


    % Punkt zawieszenia bryły
    \fill (0,0) circle (2pt);



"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class ForcedCompoundPendulumScheme2(TikZPicture):

    def _scheme_desc(self):

        code = r"""

       % Osie
    \draw[->] (-0.2,0) -- (4,0) node[right] {$x$};
    \draw[->] (0,0.2) -- (0,-4) node[below] {$y$};
    \node at (0.3,-1.7) {$\varphi$};


    % Wahadło obrócone
    \begin{scope}[rotate=25]
        \filldraw[fill=black!10, draw=black, thick] (0,0) -- (1.5,-3) arc[start angle=30, end angle=-190, radius=0.8] -- cycle;

    % Oś środka ciężkości
         \draw[dashed] (0,0) -- (1,-5);
         \filldraw[black] (0.5,-2.5) circle (1.5pt);

    % Oznaczenia
    \node at (0.8,-2.5) {$C$};
    \node at (2,-4.5) {$m, I, l$};
    %Moment wymuszający
    \node at (2.5,-3.4) {$M_0 \cos(\Omega  t)$};
    \draw[->] (0.1,-3.5) arc[start angle=-130,end angle=-40,radius=1];
    \end{scope}

    % Strzałka
    \draw[->] (0,-1.5) arc[start angle=-80,end angle=-50,radius=1.5];


    % Punkt zawieszenia bryły
    \fill (0,0) circle (2pt);


"""

        return code


class CompoundPendulumScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

 % Osie
    \draw[->] (-0.2,0) -- (4,0) node[right] {$x$};
    \draw[->] (0,0.2) -- (0,-4) node[below] {$y$};
    \node at (0.2,-2) {$\varphi$};


    % Wahadło obrócone
    \begin{scope}[rotate=30]
        \filldraw[fill=black!10, draw=black, thick] (0,0) -- (1.5,-3) arc[start angle=-20, end angle=-160, radius=1.2] -- cycle;

         \draw[dashed] (0,0) -- (0,-4);
         \filldraw[black] (0,-2) circle (1.5pt);
    \end{scope}

    % Oznaczenia
    \node at (1.3,-1.7) {$C$};
    \node at (3.5,-3) {$m, I, l$};

    % Strzałka momentu zaczepiona na linii brzegowej bryły
    \draw[->] (0,-1.5) arc[start angle=-80,end angle=-50,radius=1.5];
   %\node at (3.8,-2.3) {$M_0 \cos(\Omega  t)$};
    %\draw[->] (1,-2.5) arc[start angle=-90,end angle=-20,radius=1];

    % Punkt zawieszenia bryły
    \fill (0,0) circle (2pt);


"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class ForcedCompoundPendulumScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    % Osie
    \draw[->] (-0.2,0) -- (4,0) node[right] {$x$};
    \draw[->] (0,0.2) -- (0,-4) node[below] {$y$};
    \node at (0.2,-2) {$\varphi$};


    % Wahadło obrócone
    \begin{scope}[rotate=30]
        \filldraw[fill=black!10, draw=black, thick] (0,0) -- (1.5,-3) arc[start angle=-20, end angle=-160, radius=1.2] -- cycle;

         \draw[dashed] (0,0) -- (0,-4);
         \filldraw[black] (0,-2) circle (1.5pt);
    \end{scope}

    % Oznaczenia
    \node at (1.3,-1.7) {$C$};
    \node at (3.5,-3) {$m, I, l$};

    % Strzałka momentu zaczepiona na linii brzegowej bryły
    \draw[->] (0,-1.5) arc[start angle=-80,end angle=-50,radius=1.5];
    \node at (3.8,-2.3) {$M_0 \cos(\Omega  t)$};
    \draw[->] (1,-2.5) arc[start angle=-90,end angle=-20,radius=1];

    % Punkt zawieszenia bryły
    \fill (0,0) circle (2pt);


"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class RollingHalfDiskScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    % Definicje
    \def\r{2}       % Promień półwalca
    \def\phi{-200}    % Kąt obrotu w stopniach
    % Półwalec obrócony o kąt \phi
    \begin{scope}[rotate around={\phi:(0,0)}]  % Rotacja wokół (0,0)
        \fill[gray!20] (-\r,0) arc[start angle=180, end angle=0, radius=\r] -- cycle;
        \draw[thick] (-\r,0) arc[start angle=180, end angle=0, radius=\r];
        \draw[thick] (-\r,0) -- (\r,0);
    \end{scope}

    % Środkiem masy
    \def\cx{0}  % Współrzędne C przed rotacją
    \def\cy{4*\r/3.1416} % Przybliżenie dla półwalca

    \begin{scope}[rotate around={\phi:(0,0)}]  % Rotacja punktu C
        \filldraw[black] (0,0.9) circle (1.5pt);

        \draw[dashed] (0,0) -- (\cx,2);

    % Kąt phi
    %\draw[->] (0,1.5) arc[start angle=0, end angle=60, radius=1];
    \draw[->] (0,1.5) arc[start angle=80,end angle=100,radius=1.5];
    \node at (-0.2,1.3) {$\varphi$};
    \end{scope}
    \node at (-0.55,-0.8) {$c$};

    % Oznaczenie masy i promienia
    \node at (-\r+0.3,1) {$m, r$};

    % Osie układu współrzędnych
    \draw[->] (-2.5,-2) -- (3,-2) node[right] {$x$};
    \draw[->] (-0.69,-2.5) -- (-0.69,1) node[above] {$y$};

    \draw[dashed] (0,0) -- (0,-2);


"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class ForcedRollingHalfDiskScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    % Definicje
    \def\r{2}       % Promień półwalca
    \def\phi{-200}    % Kąt obrotu w stopniach
    % Półwalec obrócony o kąt \phi
    \begin{scope}[rotate around={\phi:(0,0)}]  % Rotacja wokół (0,0)
        \fill[gray!20] (-\r,0) arc[start angle=180, end angle=0, radius=\r] -- cycle;
        \draw[thick] (-\r,0) arc[start angle=180, end angle=0, radius=\r];
        \draw[thick] (-\r,0) -- (\r,0);
        \draw[->] (1,-0.1) arc[start angle=300,end angle=250,radius=2];


    \end{scope}

    % Środkiem masy
    \def\cx{0}  % Współrzędne C przed rotacją
    \def\cy{4*\r/3.1416} % Przybliżenie dla półwalca

    \begin{scope}[rotate around={\phi:(0,0)}]  % Rotacja punktu C
        \filldraw[black] (0,0.9) circle (1.5pt);

        \draw[dashed] (0,0) -- (\cx,2);

    % Kąt phi
    \draw[->] (0,1.5) arc[start angle=80,end angle=100,radius=1.5];

    \node at (-0.2,1.3) {$\varphi$};
    \end{scope}
    \node at (-0.55,-0.8) {$c$};

    % Oznaczenie masy i promienia
    \node at (-\r+0.3,1) {$m, r$};


    % Osie układu współrzędnych
    \draw[->] (-2.5,-2) -- (3,-2) node[right] {$x$};
    \draw[->] (-0.69,-2.5) -- (-0.69,1) node[above] {$y$};

    \draw[dashed] (0,0) -- (0,-2);
    \node at (0.4,0.8) {$M_0 \cos(\Omega  t)$};

"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class RollingBarScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

% Semi-cylinder (half-circle)
    \fill[gray!20] (-2,0) arc[start angle=180,end angle=0,radius=2];
    \draw (-2,0) arc[start angle=180,end angle=0,radius=2];
    \draw (-2,0) -- (2,0);

    % Inclined rectangular board, rotated by phi = -20 degrees, extended length
    \begin{scope}[rotate around={-20:(0,0)}]
        \draw[fill=gray!30] (-3.5,2) rectangle (2.5,2.3);
        % Center of mass dot
        \filldraw (-0.4,2.15) circle (0.05);
    \end{scope}

    % Coordinate axes (moved to be on top)
    \draw[->] (0,-0.5) -- (0,3) node[above] {$y$};
    \draw[->] (-2.5,0) -- (2.5,0) node[right] {$x$};

    % Rotation angle phi
    \draw[thick,->] (0,0) -- (0.7,1.9);
    \node at (0.2,1) {$\varphi$};

    % Angle arc
    \draw[thick,->] (0,1.5) arc[start angle=90,end angle=70,radius=1.5];

    % Center of mass
    \node at (0.5,2.6) {$C$};

    % Radius
    \node at (-1.5,0.3) {$r$};

    % Board parameters
    \node at (2.8,2.0) {$m,l,h$};

"""

        return code


# Created and/or developed by Ania Fila (AnnFyla) & Michał Kacperek (MichalKacperek)
class ForcedRollingBarScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""
    % Semi-cylinder (half-circle)
    \fill[gray!20] (-2,0) arc[start angle=180,end angle=0,radius=2];
    \draw (-2,0) arc[start angle=180,end angle=0,radius=2];
    \draw (-2,0) -- (2,0);


    % Inclined rectangular board, rotated by phi
    \begin{scope}[rotate around={-20:(0,0)}]
        \draw[fill=gray!30] (-3.5,2) rectangle (2.5,2.3);
        % Center of mass dot
        \filldraw (-0.4,2.15) circle (0.05);
        %\draw[-] (-0.4,2.15) -- (1,2)
        %\draw[->] (1,2) arc[start angle=90,end angle=70,radius=4];
    \end{scope}

    % Coordinate axes (moved to be on top)
    \draw[->] (0,-0.5) -- (0,3) node[above] {$y$};
    \draw[->] (-2.5,0) -- (2.5,0) node[right] {$x$};

    % Rotation angle phi
    \draw[thick,->] (0,0) -- (0.7,1.9);
    \node at (0.2,1) {$\varphi$};

    % Angle arc
    \draw[thick,->] (0,1.5) arc[start angle=90,end angle=70,radius=1.5];

    % Center of mass
    \node at (0.5,2.6) {$C$};

    % Radius
    \node at (-1.5,0.3) {$r$};

    % Board parameters
    \node at (2.8,2.0) {$m,l,h$};
    \node at (1.5,3.5) {$M_0 \cos(\Omega  t)$};
    \draw[->] (0.5,3) arc[start angle=90,end angle=70,radius=4];

"""

        return code


class TrolleyWithElasticPendulumScheme(TikZPicture):
    def _scheme_desc(self):
        code = r"""
\tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.3cm,post length=0.3cm,segment length=0.3cm}]
\tikzstyle{damper}=[thick,decoration={markings,
    mark connection node=dmp,
    mark=at position 0.5 with
    {
    \node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
    \draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
    \draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
    }
    }, decorate]
\tikzstyle{ground}=[fill,pattern=north east lines,draw=none,minimum width=0.75cm,minimum height=0.3cm]
\node (M) [draw,outer sep=0pt,thick,minimum width=3cm, minimum height=1.5cm] {$m_1$};
    \node (l1_label) at ($(M.south)+(305:2cm)$)   {$l_{1},k_{e}$};
    \draw[thick,gray,dashed] (M.south) -- ++(-90:2.5);
    \draw[spring,black] (M.south) --  ++(-70:3.5) coordinate (bob21) node[near end, right]{$m_{1}$};
    \fill (bob21) circle (0.2);
\node (ground) [ground,anchor=north west,xshift=-2.4cm,yshift=-0.25cm,minimum width=6cm] at (M.south west) {};
\node (wall_1) [ground,anchor=east,xshift=-2cm,yshift=0cm,minimum height=2cm ,minimum width=0.4cm] at (M.west) {};
\draw (wall_1.south east) -- (wall_1.north east);
\draw (wall_1.south east) -- (ground.north east);
\draw [thick] (M.south west) ++ (0.4cm,-0.125cm) circle (0.125cm)  (M.south east) ++ (-0.4cm,-0.125cm) circle (0.125cm);
\draw [spring]  (M.180) ++ (0cm,0.5cm) -- ($ (wall_1.0 )+(0cm,0.5cm)  $);
\draw [damper]  (M.180) ++ (0cm,-0.5cm) -- ($ (wall_1.0 )+(0cm,-0.5cm)  $);
\draw [thin] (M.east) (0,1.5) -- (0,0.75);
\draw [-latex,ultra thick] (M.east) ++ (-1.5cm,1.5cm) -- +(1.5cm,0cm);
\node (k1) at (-2.5cm,1cm) {$k_1$};
\node (c1) at (-2.5cm,0cm) {$c_1$};
\node(ft) at (0.75cm,1.75cm) {$f(t)$};

"""
        return code


class DoublePendulumMeasureScheme(TikZPicture):
    def _scheme_desc(self):
        code = r"""

\node (B1) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center,] at (0,0) {Wymuszenie};
\node (B2) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4cm]B1.east) {Main member};
\node (UAcc) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4cm]B2.east) {Upper accelerometer};

\draw [thick, ->] (B1.east) -- (B2.west);
\draw [thick,double, ->] (B2.east) -- (UAcc.west);

\node (B3) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-4cm]B2.south) {Vibration absorber};

\node (LAcc) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-4cm]UAcc.south) {Lower accelerometer};

\draw [thick, <->] (B2.south) -- (B3.north);
\draw [thick,double, ->] (B3.east) -- (LAcc.west);

\node (DAQ) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-1.75cm]UAcc.south) {DAQ};

\draw [thick, ->] (LAcc.north) -- (DAQ);
\draw [thick, ->] (UAcc.south) -- (DAQ);
"""

        return code


# Created by Anna Mackojc
class ColouredDoublePendulumMeasureScheme(TikZPicture):
    def _scheme_desc(self):
        code = r"""

\definecolor{excitationColor}{RGB}{255, 102, 102} % Red
\definecolor{mainColor}{RGB}{102, 178, 255} % Blue
\definecolor{accColor}{RGB}{255, 178, 102} % Orange
\definecolor{ptmdColor}{RGB}{153, 204, 102} % Green
\definecolor{daqColor}{RGB}{204, 153, 255} % Purple

% Nodes
\node (B1) [draw, fill=excitationColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at (0,0) {\textbf{Excitation}};

\node (B2) [draw, fill=mainColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at ([xshift=4cm]B1.east) {\textbf{Main Member}};

\node (UAcc) [draw, fill=accColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at ([xshift=4cm]B2.east) {\textbf{Upper Accelerometer}};

% Arrows
\draw [thick, ->, >=latex, black] (B1.east) -- (B2.west);
\draw [thick, double, ->, >=latex, black] (B2.east) -- (UAcc.west);

% Vibration Absorber
\node (B3) [draw, fill=ptmdColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at ([yshift=-4cm]B2.south) {\textbf{Vibration Absorber (PTMD)}};

\node (LAcc) [draw, fill=accColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at ([yshift=-4cm]UAcc.south) {\textbf{Lower Accelerometer}};

% Arrows
\draw [thick, <->, >=latex, black] (B2.south) -- (B3.north);
\draw [thick, double, ->, >=latex, black] (B3.east) -- (LAcc.west);

% DAQ System
\node (DAQ) [draw, fill=daqColor!30, thick, minimum width=3.2cm, minimum height=1.2cm, align=center, rounded corners] at ([yshift=-1.75cm]UAcc.south) {\textbf{DAQ System}};

% Connections to DAQ
\draw [thick, ->, >=latex, black] (LAcc.north) -- (DAQ);
\draw [thick, ->, >=latex, black] (UAcc.south) -- (DAQ);
"""

        return code


# Created by Anna Mackojc
class SegmentBeamScheme(TikZPicture):
    def _scheme_desc(self):
        code = r"""

     \coordinate (origo) at (0,0);
     \coordinate (pivot) at (1,5);

    % draw axes
    \fill[black] (origo) circle (0.05);
    \draw[thick,gray,->] (origo) -- ++(15,0) node[black,right] {$x$};
    \draw[thick,gray,->] (origo) -- ++(0,4) node (mary) [black,right] {$y$};

    \filldraw[fill=blue!30, draw=black, thick, rotate around = {-30:(origo)}] (origo) rectangle ++ (3,0.25) coordinate (bob1) node[near end, below right]{};
    \filldraw[fill=blue!30, draw=black, thick, rotate around = {-90:(origo)}] (bob1) ++ (0,0) rectangle ++ (0.25,5) coordinate (bob2) node[anchor=north east]{};
    \filldraw[fill=blue!30, draw=black, thick, rotate around = {30:(origo)}] (bob2) ++ (0.15,-.05) rectangle ++ (3,0.25) coordinate (bob3) node[near end, right]{};

    \coordinate (new_origo) at (2,0);
    \filldraw[fill=blue!30, draw=black, thick, rotate around = {-30:(origo)}] (new_origo) rectangle ++ (3.5,0.25) coordinate (bob1) node[near end, below right]{};
    \filldraw[fill=blue!30, draw=black, thick, rotate around = {-90:(origo)}] (bob1) ++ (0,0) rectangle ++ (0.25,4) coordinate (bob2) node[anchor=north east]{};
    \filldraw[fill=blue!30, draw=black, thick, rotate around = {30:(origo)}] (bob2) ++ (0.15,-.05) rectangle ++ (3.5,0.25) coordinate (bob3) node[near end, right]{};


    \draw [shift={(bob1)},yshift=-5,domain=2*pi+310/180*pi:295/180*pi+7*pi,variable=\t,smooth,samples=75] plot ({\t r}: {0.0006*\t*\t});
    \draw [shift={(bob2)},xshift=1,yshift=4,domain=2*pi+315/180*pi:305/180*pi+7*pi,variable=\t,smooth,samples=75] plot ({\t r}: {0.0006*\t*\t});
"""

        return code


# Created by Krzysztof Twardoch
# Verified by BCh


class TrolleyWithPendulumMeasureScheme(TikZPicture):
    def _scheme_desc(self):
        code = r"""

\node (B1) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center,] at (0,0) {Kinematic excitation \\ (Vibrating table)};
\node (B2) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4cm]B1.east) {Main member};
\node (UAcc) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([xshift=4cm]B2.east) {Upper accelerometer};

\draw [thick, ->] (B1.east) -- (B2.west);
\draw [thick,double, ->] (B2.east) -- (UAcc.west);

\node (B3) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-4cm]B2.south) {Vibration absorber};

% \node (LAcc) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-4cm]UAcc.south) {Lower accelerometer};

\draw [thick, <->] (B2.south) -- (B3.north);
% \draw [thick,double, ->] (B3.east) -- (LAcc.west);

\node (DAQ) [draw,outer sep=0pt,thick,minimum width=3.0cm, minimum height=1.0cm, align=center] at ([yshift=-1.75cm]UAcc.south) {DAQ};

% \draw [thick, ->] (LAcc.north) -- (DAQ);
\draw [thick, ->] (UAcc.south) -- (DAQ);

"""
        return code


class RLScheme(TikZPicture):

    def _scheme_desc(self):

        code = r"""

    \draw
    (0,0) to[battery, l_=$V$] (0,-3)
    (0,0) to[resistor, l=$R$] (6,0)
    (6,-3) to[american inductor, l_=$L$] (6,0)
    (0,-3) -- (6,-3)
    [->] (0,0) -- (0.5,0) node[midway, above] {$\dot q$};
;

"""

        return code


class DynPyTree(TikZPicture):

    def _scheme_desc(self):

        code = r"""
% DynPy main block
\draw [line width=0.8pt, rounded corners=9.6] (5,13.75) rectangle (7.5,12.5);
\node [font=\LARGE] at (6.25,13.25) {DynPy};

% Horizontal lines to other modules
\draw [line width=0.8pt] (5,13.25) -- (3.75,13.25);
\draw [line width=0.8pt] (7.5,13.25) -- (8.75,13.25);

% Vertical line to models
\draw [line width=0.8pt] (6.25,12.5) -- (6.25,11.25);

% Boxes for other modules
\draw [line width=0.8pt, rounded corners=9.6] (8.75,13.75) rectangle (11.25,12.5);
\draw [line width=0.8pt, rounded corners=9.6] (5,11.25) rectangle (7.5,10);
\draw [line width=0.8pt, rounded corners=9.6] (1.25,12.5) rectangle (3.75,13.75);
\node at (6.25,10.75) {models};
\node at (10,13.25) {solvers};
\node at (2.5,13.25) {utilities};

% Mechanics and gears.py
\draw [line width=0.8pt, rounded corners=9.6] (5,7.5) rectangle (7.5,8.75);
\node at (6.25,8.25) {mechanics};
\draw [line width=0.8pt] (6.25,10) -- (6.25,8.75);

\draw [line width=0.8pt, rounded corners=9.6] (5,5) rectangle (7.5,6.25);
\node at (6.25,5.75) {gears.py};
\draw [line width=0.8pt] (6.25,7.5) -- (6.25,6.25);

% Bottom modules: nonlinear.py and report.py
\draw [line width=0.8pt, rounded corners=9.6] (8.75,10) rectangle (11.25,11.25);
\node [font=\large] at (10,10.75) {nonlinear.py};

\draw [line width=0.8pt, rounded corners=9.6] (1.25,10) rectangle (3.75,11.25);
\node at (2.5,10.75) {report.py};

% Arrows
\draw [line width=0.8pt, ->, >=Stealth] (2.5,12.5) -- (2.5,11.25);
\draw [line width=0.8pt, ->, >=Stealth] (10,12.5) -- (10,11.25);
\draw [line width=0.8pt, ->, >=Stealth] (6.25,7.5) -- (6.25,6.25);
"""
        return code


class TrolleyWithPendulumTestStandScheme(TikZPicture):

    def _scheme_desc(self):

        half_force_label = r"{\color{violet} $\frac{F}{2}$}"
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # PL
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # text1 = r'{\color{violet} Silnik elektryczny}'
        # text2 = r'{\color{violet} Przekładnia pasowa}'
        # text3 = r'{\color{violet} Building frame model}' # Obiekt badań: Model ramy budynku => Podatna rama
        # text4 = r'{\color{violet} APTMD}'
        # text5 = r'{\color{violet} Akcelerometry}'
        # text6 = r'{\color{violet} Stół wibracyjny}' # Stół ruchomy
        # text7 = r'{\color{violet} Korbowód}'
        # text8 = r'{\color{violet} Falownik}'
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # EN
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        text1 = r"{\color{violet} Electric Motor}"
        text2 = r"{\color{violet} Belt Drive}"
        text3 = r"{\color{violet} Flexible Frame}" # Test Object: Building Frame Model => Flexible Frame
        text4 = r"{\color{violet} APTMD}"
        text5 = r"{\color{violet} Accelerometers}"
        text6 = r"{\color{violet} Shake–Table}"
        text7 = r"{\color{violet} Connecting Rod}"
        text8 = r"{\color{violet} Inverter}"  # Controller
        # ––––––––––––––––––––––––––––––––––––––––––––––––

        pos1_x = 0
        pos1_y = 0

        del_row = -6.1
        del_2row = 0.4

        code = f"""

    %\draw[ultra thick,dashed,black] (3.9,-0.25) -- (3.9,-6.2);

    %\draw[ultra thick,dashed,black] (7.9,-1.5) -- (7.9,-6.2);
    %\draw[ultra thick,dashed,black] (11.9,-1.5) -- (11.9,-6.2);


    %\\node[inner sep=7pt,fill=white,align=center] (text) at (1.25,-2.25) {{label}};

    \\node[inner sep=5pt,fill=white] (pic1) at (7.25cm,-0.5cm) {{\includegraphics[width=12cm]{{../images/test_stand_twp.png}}}};
    %%%%%%  MID COL
    \\node[inner sep=2pt,align=center,draw] (text1) at ({pos1_x+4},{pos1_y+3.9}) {text1};
    \\node[inner sep=2pt,align=center,draw] (text2) at ({pos1_x+3.5},{pos1_y+2.7}) {text2};
    \\node[inner sep=2pt,align=center,draw] (text3) at ({pos1_x+6.6},{pos1_y+4.7}) {text3};
    \\node[inner sep=2pt,align=center,draw] (text4) at ({pos1_x+9.25},{pos1_y+2}) {text4};
    \\node[inner sep=2pt,fill=white,align=center,draw] (text5) at ({pos1_x+12.2},{pos1_y+2.2}) {text5};
    \\node[inner sep=2pt,align=center,draw] (text6) at ({pos1_x+10.8},{pos1_y-3.8}) {text6};
    \\node[inner sep=2pt,align=center,draw] (text7) at ({pos1_x+7.7},{pos1_y-3.5}) {text7};
    \\node[inner sep=2pt,align=center,draw] (text8) at ({pos1_x+7.15},{pos1_y-5}) {text8};

    \\draw[ultra thick,red,->] (text1.south) -- (6,2.8);
    \\draw[ultra thick,red,->] (text2.south) -- (4.3,1.7);
    \\draw[ultra thick,red,->] (text3.south) -- (8,3.8);
    \\draw[ultra thick,red,->] (text4.north) -- (9.83,2.6);
    \\draw[ultra thick,red,->] (text5.north) -- (12.11,4.3);  %A1
    \\draw[ultra thick,red,->] (text5.south) -- (11.33,-0.1); %A2
    \\draw[ultra thick,red,->] (text6.north) -- (10,-0.25);
    \\draw[ultra thick,red,->] (text7.north) -- (8.2,0.3);
    \\draw[ultra thick,red,->] (text8.north) -- (5.6,-3.9);
    """
        return code

class TrolleyWithPendulumTestStandSchemeNo(TikZPicture):

    def _scheme_desc(self):

        half_force_label = r"{\color{violet} $\frac{F}{2}$}"
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # PL
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # text1 = r'{\color{violet} 1}" # Silnik elektryczny
        # text2 = r'{\color{violet} 2}" # Przekładnia pasowa
        # text3 = r'{\color{violet} 3}" # Building frame model # Obiekt badań: Model ramy budynku => Podatna rama
        # text4 = r'{\color{violet} 4}" # APTMD
        # text5 = r'{\color{violet} 5}" # Stół wibracyjny # Stół ruchomy
        # text7 = r'{\color{violet} 6}" # Korbowód
        # text5 = r'{\color{violet} 7}" # Akcelerometry
        # text8 = r'{\color{violet} 8}" # Falownik
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        # EN
        # ––––––––––––––––––––––––––––––––––––––––––––––––
        text1 = r"{\color{violet} \LARGE 1}" # Electric Motor
        text2 = r"{\color{violet} \LARGE 2}" # Belt Drive
        text3 = r"{\color{violet} \LARGE 3}" # Flexible Frame # Test Object: Building Frame Model => Flexible Frame
        text4 = r"{\color{violet} \LARGE 4}" # APTMD
        text5 = r"{\color{violet} \LARGE 5}" # Shake–Table
        text6 = r"{\color{violet} \LARGE 6}" # Connecting Rod
        text7 = r"{\color{violet} \LARGE 7}" # Accelerometers
        text8 = r"{\color{violet} \LARGE 8}" # Inverter # Controller
        # ––––––––––––––––––––––––––––––––––––––––––––––––

        pos1_x = 0
        pos1_y = 0

        del_row = -6.1
        del_2row = 0.4

        code = f"""

    %\draw[ultra thick,dashed,black] (3.9,-0.25) -- (3.9,-6.2);

    %\draw[ultra thick,dashed,black] (7.9,-1.5) -- (7.9,-6.2);
    %\draw[ultra thick,dashed,black] (11.9,-1.5) -- (11.9,-6.2);


    %\\node[inner sep=7pt,fill=white,align=center] (text) at (1.25,-2.25) {{label}};

    \\node[inner sep=5pt,fill=white] (pic1) at (7.25cm,-0.5cm) {{\includegraphics[width=12cm]{{../images/test_stand_twp.png}}}};
    %%%%%%  MID COL
    \\node[inner sep=2pt,align=center,draw] (text1) at ({pos1_x+4.9},{pos1_y+3.9}) {text1};
    \\node[inner sep=2pt,align=center,draw] (text2) at ({pos1_x+3.9},{pos1_y+2.7}) {text2};
    \\node[inner sep=2pt,align=center,draw] (text3) at ({pos1_x+7.4},{pos1_y+4.6}) {text3};
    \\node[inner sep=2pt,align=center,draw] (text4) at ({pos1_x+9.1},{pos1_y+3.5}) {text4};
    \\node[inner sep=2pt,align=center,draw] (text5) at ({pos1_x+10.3},{pos1_y-3.8}) {text5};
    \\node[inner sep=2pt,align=center,draw] (text6) at ({pos1_x+7.7},{pos1_y-3.7}) {text6};
    \\node[inner sep=2pt,fill=white,align=center,draw] (text7) at ({pos1_x+12.3},{pos1_y+2.2}) {text7};
    \\node[inner sep=2pt,align=center,draw] (text8) at ({pos1_x+6.6},{pos1_y-5}) {text8};

    \\draw[ultra thick,red,->] (text1.south) -- (6,2.8);
    \\draw[ultra thick,red,->] (text2.south) -- (4.3,1.7);
    \\draw[ultra thick,red,->] (text3.south) -- (8,3.8);
    \\draw[ultra thick,red,->] (text4.south) -- (9.8,2.8);
    \\draw[ultra thick,red,->] (text5.north) -- (10,-0.25);
    \\draw[ultra thick,red,->] (text6.north) -- (8.2,0.3);
    \\draw[ultra thick,red,->] (text7.north) -- (12.11,4.3);  %A1
    \\draw[ultra thick,red,->] (text7.south) -- (11.33,-0.1); %A2
    \\draw[ultra thick,red,->] (text8.north) -- (5.65,-3.9);
    """
        return code

class NationalDAQSystemScheme(TikZPicture):

    def _scheme_desc(self):

        half_force_label = r"{\color{violet} $\frac{F}{2}$}"
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # PL
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # text1 = r'{\color{violet} NI cDAQ-9174}' # System Akwizycji Danych
        # text2 = r'{\color{violet} Karta pomiarowa (NI 9230)}'
        # text3 = r'{\color{violet} Akcelerometr}'
        # text4 = r'{\color{violet} Przewód pomiarowy}'
        # text5 = r'{\color{violet} Zasilacz}'
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        # EN
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
        text1 = r"{\color{violet} NI cDAQ-9174}"  # Data Acquisition System
        text2 = r"{\color{violet} Measurement Card NI 9230}"
        text3 = r"{\color{violet} Accelerometer}"
        text4 = r"{\color{violet} Measurement Cable}"
        text5 = r"{\color{violet} Power supply}"
        # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

        pos1_x = 0
        pos1_y = 0

        del_row = -6.1
        del_2row = 0.4

        code = f"""

    \\node[inner sep=5pt,fill=white] (pic1) at (7.25cm,-0.5cm) {{\includegraphics[width=12cm]{{../images/measurement_kit.png}}}};
    %%%%%%  MID COL
    \\node[inner sep=2pt,align=center,draw] (text1) at ({pos1_x+4.5},{pos1_y+3}) {text1};
    \\node[inner sep=2pt,align=center,draw] (text2) at ({pos1_x+4.5},{pos1_y+4.4}) {text2};
    \\node[inner sep=2pt,align=center,draw] (text3) at ({pos1_x+3.8},{pos1_y-5}) {text3};
    \\node[inner sep=2pt,align=center,draw] (text4) at ({pos1_x+11.4},{pos1_y-4.5}) {text4};
    \\node[inner sep=2pt,fill=white,align=center,draw] (text5) at ({pos1_x+3.6},{pos1_y+1.5}) {text5};

    \\draw[ultra thick,red,->] (text1.south) -- (6.4,2.24);
    \\draw[ultra thick,red,->] (text2.south) -- (7.17,3.5);
    \\draw[ultra thick,red,->] (text3.north) -- (6.28,-2.95);
    \\draw[ultra thick,red,->] (text3.south) -- (6.1,-5.8);
    \\draw[ultra thick,red,->] (text4.north) -- (12,-1.2); %C1
    \\draw[ultra thick,red,->] (text4.north) -- (8.3,-3.3); %C2
    \\draw[ultra thick,red,->] (text5.south) -- (3.15,0.1);
    """
        return code
class SuspensionModelOscillatorScheme(GearModelOscillatorScheme):

    def _scheme_desc(self):

        code = r"""

    \coordinate (origo) at (0,0);

    \tikzstyle{spring}=[thick,decorate,decoration={zigzag,pre length=0.3cm,post length=0.3cm,segment length=0.3cm}]

    \tikzstyle{damper}=[thick,decoration={markings,
    mark connection node=dmp,
    mark=at position 0.5 with
    {
    \node (dmp) [thick,inner sep=0pt,transform shape,rotate=-90,minimum width=15pt,minimum height=3pt,draw=none] {};
    \draw [thick] ($(dmp.north east)+(2pt,0)$) -- (dmp.south east) -- (dmp.south west) -- ($(dmp.north west)+(2pt,0)$);
    \draw [thick] ($(dmp.north)+(0,-5pt)$) -- ($(dmp.north)+(0,5pt)$);
    }
    }, decorate]

    \tikzstyle{ground}=[fill,pattern=north east lines,draw=none,minimum width=0.75cm,minimum height=0.3cm]



    %draw axes
    %\fill[black] (origo) circle (0.05);


    \node (M) [draw,outer sep=0pt,thick,minimum width=3cm, minimum height=1.5cm,yshift=2cm] at (origo) {$m_{red}$};
    %\node (m1_label) at ($(M)+(0:0.5cm)$) {$m$};

    \draw [spring]    ([xshift=-1cm]M.south) --++(0,-3cm) node (k_tire_end){} node[midway,left] {$k_{s}$};

        \draw [damper] ([xshift=1cm]M.south) --++(0,-3cm) node (c_tire_end){}  node[midway,right=0.25cm] {$c_{s}$};

    %\fill [black] (M.center) circle (0.05);

    % \node (beam) [fill=gray,anchor=north,xshift=0cm,yshift=0cm,minimum width=2cm,minimum height=0.05pt] at (origo) {} node[right]{$m_{beam}=0$};

    % \draw [ultra thick] (k_sup_spring_left.west) -- (k_sup_spring_right.east) node[midway] (beam_center) {};


    %         \node (road_sur) [fill=black,anchor=north,xshift=0cm,yshift=-2cm,minimum width=2cm,minimum height=0.05cm] at (beam) {};

    %beam_center.center

    %    \draw[spring] ([xshift=-0.0cm]k_sup_spring_left.center) -- ++(0,-1.5cm) node (k_tire_end) {} node[midway,left] {$k_{s}$};

    %        \draw[damper] ([xshift=+0.5cm](k_sup_spring_right.center) -- ++(0,-2cm) node (c_tire_end) {} node[midway,right=0.25cm] {$c_{s}$};

    \draw[ultra thick] ([xshift=-0.0cm]k_tire_end.west) -- ([xshift=+0.0cm]c_tire_end.east) node (force_attachment_point) {};


    \draw [thin] (force_attachment_point) -- +(0,0) coordinate (force_leader_edn);
    \draw [-latex,ultra thick] (force_leader_edn) -- +(0,-1.5cm) node[above right] (force_t)  {$F$};

    \draw [thin] (0,-2.8) -- +(0,1) coordinate (force_leader_edn) {} node[midway,left] {$u$};
    \draw[->]        (0,-2)   -- (0,-2.8);

    """

        return code