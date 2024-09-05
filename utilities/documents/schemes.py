from ..report import TikZPicture


class TuningTikZDiagram(TikZPicture):

    def _scheme_desc(self):

        code=r"""
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

        code=r"""
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

        
        code=r"""

    
    

    \draw
    (0,0) to[american voltage source, l_=$U_{oc}$] (0,-3)
    (0,0) to[R,european,-*,l=$R_{0}$] (6,0) node[label={right:$+$}] {}
    (0,-3) to[short, -*] (6,-3) node[label={right:$-$}] {}
;

"""

        return code


class TheveninTikZDiagramVer1(TikZPicture):

    def _scheme_desc(self):

        
        code=r"""

    
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

        
        code=r"""

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

        
        code=r"""

    
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

        
        code=r"""

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

        
        code=r"""

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

        
        code=r"""

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

        
        code=r"""
        
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

        
        code=r"""
        
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
    