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


\draw [thin] (force_attachment_point) -- +(0,-1.5) coordinate (force_leader_edn);
\draw [-latex,ultra thick] (force_leader_edn) -- +(-0,0cm) node[above right] (force_t)  {$F$};

\draw [thin] (force_attachment_point) -- +(2,0) coordinate (force_leader_edn);


"""

        return code