from sympy import *
import sympy as sym
import itertools as itools
from pylatex import (Document, Package, Command
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )
from pylatex.base_classes.containers import Container
#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)

from ..adaptable import *
from ..report import *
from ..report import display as ip_display


def display(obj):

    #ip_display(IPMarkdown('=='*25 +'\n \n \n ==== if u see it let me (@bogumilchilinski) know' ))
    return ip_display(obj)


class MultivariableTaylorSeries(Expr):
    
    def __new__(cls,expr, variables,*args, n=2, x0=None):
        
        obj=super().__new__(cls,expr,variables,*args)
        obj._vars = variables
        obj._order = n
        obj._op_point = x0
        
        obj._expr_symbol = None
        
        return obj

    def _set_default_op_point(self):
        '''
        It sets 0 as op_point if x0 (look __new__) is not provided. For lack of op_point, the result is the same as MacLaurina series 
        '''
        
        if self._op_point is None:
            self._op_point = {coord:0 for coord in self._vars}
        
        return self._op_point
    
    def _args_shifted(self,*args):
        
        self._set_default_op_point()
        
        args_shifted = {
            arg: arg - arg_shift
            for arg, arg_shift in self._op_point.items()
        }
        
        return args_shifted

    def _diff_orders_dict(self):
            
        order_max = self._order
        args = self._vars
        args_shifted = self._args_shifted()
        
            
        diff_orders_list = sum([
            list(itools.combinations_with_replacement(args, order))
            for order in range(1, order_max + 1, 1)
        ], [])
        
        
        diff_orders_dict = {
            comp: (sym.Mul(*comp).subs(args_shifted) / sym.Mul(*[
                sym.factorial(elem) for elem in sym.Poly(sym.Mul(*comp), *args).terms()[0][0]
            ])).doit()
            for comp in diff_orders_list
        }

        return diff_orders_dict

    
    def _diff_symbols_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        expr=self.args[0]
        op_point=self._op_point
        
        
        
        return {S.Zero:Subs(expr,list(op_point.keys()),list(op_point.values())),**{args_tmp:Subs(Derivative(expr,*args_tmp,evaluate=False),list(op_point.keys()),list(op_point.values())) 
            for args_tmp, poly in diff_orders_dict.items()}}

    def _diff_expr_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        expr=self.args[0]
        op_point=self._op_point
        
        return {S.Zero:expr.subs(op_point),**{args_tmp:expr.diff(*args_tmp).subs(op_point)
            for args_tmp, poly in diff_orders_dict.items()}}
    
    def _components_dict(self):
        
        diff_orders_dict=self._diff_orders_dict()
        derivatives_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return {
                Subs(expr,list(op_point.keys()),list(op_point.values())):expr.subs(op_point),
                **{derivatives_dict[args_tmp] : expr.diff(*args_tmp).subs(op_point) for args_tmp, poly in diff_orders_dict.items()}
        }
        
    def _series(self):
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return expr.subs(op_point).doit()+Add(*[expr.diff(*args_tmp).subs(op_point).doit() * poly for args_tmp, poly in diff_orders_dict.items()],evaluate=False)
    
    def _symbolic_sum(self):
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return Subs(expr,list(op_point.keys()),list(op_point.values()))+Add(*[Mul(diff_dict[args_tmp]  ,poly,evaluate=True) for args_tmp, poly in diff_orders_dict.items()],evaluate=False)
    
    def _latex(self,*args):
        
        diff_orders_dict=self._diff_orders_dict()
        diff_dict=self._diff_symbols_dict()
        
        expr=self.args[0]
        op_point=self._op_point
        
        return '+'.join([latex(Mul(diff_dict[args_tmp]  ,poly,evaluate=True)) for args_tmp, poly in diff_orders_dict.items()])

    
    def calculation_steps(self,expr_symbol=None,form=None):

        obj = self
        
        if expr_symbol is None:
            obj._expr_symbol = self.args[0]
            
        obj_sym = self.__class__(expr_symbol,self._vars, n=self._order, x0=self._op_point)
        
        expr_dict=(self._diff_expr_dict())
        diffs_dict=(obj_sym._diff_symbols_dict())
        

        
        
        return [Eq(diffs_dict[key],expr_dict[key].doit())   for  key  in diffs_dict.keys()]
    
    def __str__(self,*args):
        return (self.args[0]).__str__()
    
    def __repr__(self,*args):
        return (self.args[0]).__repr__()

    
    
    
    def doit(self,*args):
        return self._series()


class ReportComponent(Subsection):

    latex_name = 'subsection'
    packages=[
              Package('standalone'),
              Package('siunitx')
             ]
    
    title='Report generic component'

    def __init__(self, reported_object, title=None, numbering=False, *, label=True, **kwargs):
        """
        Args
        ----
        title: str
            The section title.
        numbering: bool
            Add a number before the section title.
        label: Label or bool or str
            Can set a label manually or use a boolean to set
            preference between automatic or no label
        """



        self.reported_object = reported_object #it's forced by pylatex library
        
        
        if title is None:
            title = self.dynamic_title()
        
        
        
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)
        CurrentContainer(self)
        
        if self.title is not None:
            ip_display(IPMarkdown(f'## {self.title}'))
        
        self.append_elements()
        
        
    def dynamic_title(self):
        return self.title
        
    def append_elements(self):
        pass

    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame

    @property
    def reported_object(self):

        from ...solvers.linear import ODESystem
        return self._reported_object       
#         if isinstance(self._reported_object, ODESystem):
#             return self._reported_object
#         else:
#             return self._reported_object._ode_system
            
    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj

    @property
    def _system(self):
        print('Kod do poprawienia #################### bo stosujesz starą zmienną self._system')
        print('zamień se self._system na self.reported_object')
        
        return self.reported_object
    
class ODESystemOperations(ReportComponent):
    
    def __init__(self, reported_object):
        from ...solvers.linear import ODESystem
        
        if isinstance(reported_object, ODESystem):
            ode_sys = reported_object
            ode_dummy = self.get_dummy()
            reported_object = (ode_sys, ode_dummy)

        super().__init__(reported_object)


    def get_dummy(self):
        from ...solvers.linear import ODESystem
        
        t = Symbol('t')
        m = Symbol('m',postive=True)
        g = Symbol('g',postive=True)
        c = Symbol('c',postive=True)
        b = Symbol('b',postive=True)
        x = Function('x')(t)
        
        return  ODESystem(odes=Matrix([m * x.diff(t,t)]), dvars=Matrix([x]), odes_rhs=Matrix([- c * x.diff(t) - b * x.diff(t,t) ], ode_order=2))
    
    def get_zero_ode_sys(self):
        from ...solvers.linear import ODESystem
        
        t = Symbol('t')
        x = Function('x')(t)
        return ODESystem(odes=0, dvars=x)

class TitlePageComponent(Environment):
    
    latex_name='titlepage'
    
    def __init__(self, system=None, options=None, arguments=None, start_arguments=None,
                 **kwargs):
        r"""
        Args
        ----
        options: str or list or  `~.Options`
            Options to be added to the ``\begin`` command
        arguments: str or list or `~.Arguments`
            Arguments to be added to the ``\begin`` command
        start_arguments: str or list or `~.Arguments`
            Arguments to be added before the options
        """

        self.system = system
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        
        
        super().__init__(options=options, arguments=arguments, start_arguments=start_arguments,**kwargs)
        
        if self.system is not None:

        
            system = self.system


            
            self.append(NoEscape('\centering'))

            self.append(NoEscape('\\Huge DRGANIA MECHANICZNE \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))

            
            
            if len(system.q)==1:
                dof_str = 'JEDNYM STOPNIU SWOBODY'
            else:
                dof_str = 'WIELU STOPNIACH SWOBODY'
                
            if system._dissipative_potential==0 or system._dissipative_potential is None:
                damping_str = 'NIETŁUMIONE'
            else:
                damping_str = 'TŁUMIONE'

           
            self.append(NoEscape(f'\\Large {damping_str} UKŁADY O {dof_str} \n \n'))    
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(NoEscape(f'{system._label} \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(Command('MyAuthor'))
            self.append(NoEscape(f'\\par'))
            self.append(Command('vspace',arguments='1cm'))
            
            #self.append(Command('vspace',arguments='1cm'))
            #self.append(NoEscape(f'\\protect\\par'))
            #self.append(NewLine())
            self.append(Command('MyDate'))

 


    
class ExemplaryPictureComponent(ReportComponent):
    
    title="Przykład rzeczywistego obiektu"
    packages=[Package('float')] 
                                
    def append_elements(self):
        
        system = self._system

        display(ReportText(f'''Ilustracja przedstawia rzeczywisty obiekt mechaniczny, będący przedmiotem modelowania i analizy dynamicznej.
                            '''))

        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example(),width='8cm')

        display(ReportText(f'''Model dynamiczny układu określa się na podstawie analizy rozważanego przypadku. 
        Należy pamiętać, że stopień odwzorowania (poziom abstrakcji) modelu zależy od tego do czego planuje się go używać.
                            '''))            

        
class SchemeComponent(ExemplaryPictureComponent):
    title="Schemat układu"


    def append_elements(self):
        
        system = self._system

        display(ReportText(f'''Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, wyznaczony na podstawie uprzedniej analizy rzeczywistego obiektu.
                            '''))
        
        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._scheme(),width='8cm')

        display(ReportText(f'''Analizując przedstawiony układ można stwierdzić, że jego liczba stopni swobody to {len(system.q)}.
                            '''))



class NumericalAnalysisComponent(ExemplaryPictureComponent):
    title="Symulacja numeryczna"

    _default_ics = None
    _default_parameters = {}
    _default_tspan = np.linspace(0,1,1000)
    
    @classmethod
    def set_default_ics(cls,default_ics=None):
        
        if default_ics is None:
            cls._default_ics = default_ics
        
        return cls

    @classmethod
    def set_default_parameters(cls,default_parameters=None):
        
        if default_parameters is None:
            cls._default_parameters = default_parameters
        
        return cls
    
    @classmethod
    def set_default_tspan(cls,default_tspan=None):
        
        if default_tspan is None:
            cls._default_tspan = default_tspan
        
        return cls
    

    def append_elements(self):
        
        system = self._system

        display(ReportText(f'''Dla Damiana :P
                            '''))
        if self._default_ics is None:
            ics = list(system.Y*0.0)
        else:
            ics = self._default_ics
            
        sym_res = system.subs(self._default_parameters).numerized().compute_solution(self._default_tspan,ics)

        LatexDataFrame.formatted(sym_res).plotted(preview=True)
        

        display(ReportText(f'''Dla Damiana :P
                            '''))

class KineticEnergyComponent(ReportComponent):
    
    title="Energia kinetyczna"
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

        


        display(ReportText(f'''
                           Energia kinetyczna układu wyrażona jest wzorem:
                           
                           '''))
        

        display(SympyFormula( Eq(Symbol('T'),
                     dyn_sys_lin._kinetic_energy) , marker=None))
               
        display(ReportText(f'''
                           Wyznaczona wielkość określa energię układu wynikającą z jego własności inercyjnych (energię zmagazynowaną w elementach bezwładnych).
                           
                           '''))
    
    
    
class PotentialEnergyComponent(ReportComponent):
    
    title="Energia potencjalna"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Energia potencjalna układu wyrażona jest wzorem:
                           '''))

        display(SympyFormula( Eq(Symbol('V'),
                     dyn_sys_lin._potential_energy ), marker=None))

        display(ReportText(f'''
                           Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt.
                           '''))        
        
class DissipationComponent(ReportComponent):  # Marcel
    
    title="Dyssypacyjna funkcja Rayleigh'a"
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys

        

        display(ReportText(f'''
                           Energia rozpraszana tłumieniem wyrażona jest wzorem:
                           
                           '''))

        display(SympyFormula( Eq(Symbol('D'),
             dyn_sys_lin._dissipative_potential) , marker=None ))

        display(ReportText(f'''
                           Podana zależność stanowi potencjał dysynpacyjny Rayleigh'a, 
                           który poddany różniczkowaniu względem wektora prędkości uogólnionych pozwala na określenie sił wiskotycznego tłumienia.
                           '''))
        
class SpringForce(ReportComponent):
    
    title="Siła od sprężyny"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Siła od sprężyny wyrażona jest wzorem:
                           '''))
        display(SympyFormula( Eq(Symbol('F'),-1*system._left_mount.stiffness*system.z,evaluate=False)))
        
        display(ReportText(f'''
                           zamiast x używam steady solution
                           '''))

        display(SympyFormula( Eq(Symbol('F'),-1*system._left_mount.stiffness*system.steady_solution(),evaluate=False)))


        display(ReportText(f'''
                           Siła od sprężyny, zwana również siłą naciągu, pojawia sie przy ściskaniu lub rozciaganiu. Siła, która działa jest przeciwnie skierowana do ruch i chce przywrócić do pierwotnego jej położenia. Zależy od sztywności sprężyny k oraz od tego o ile została rozciagnieta bądź skrócona x.
                           '''))
        
class DamperForce(ReportComponent):
    
    title="Siła tłumienia"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Siła tłumienia wyrażona jest wzorem:
                           '''))

        display(SympyFormula( Eq(Symbol('F'),-1*system.c*system.v,evaluate=False)))

        display(ReportText(f'''
                           zastępuje prędkość jako pochodna steady stolution po czasie:
                           '''))

        display(SympyFormula( Eq(Symbol('F'),-1*system.c*system.steady_solution().diff(system.ivar),evaluate=False)))

        display(ReportText(f'''
                           Siła tłumienia zmniejsza amplitude drgań, ma zwrot przeciwny do prędkości. Zależy od współczynnika tłumienia b oraz od prędkości v.
                           '''))                
class LogarithmicDecrement(ReportComponent):
    
    title="logarytmiczny dekrement tłumienia liczony z amplitud"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Logarytmiczny dekrement tłumienia wyrażona jest wzorem przy pomocy amplitud:
                           '''))

        display(SympyFormula( Eq(Symbol('delta'),log(system.A_n1/system.A_n2))))


        display(ReportText(f'''
                           Zaprezentowana zależność opisuje bezwymaiarową wielkość charakteryzującą intensywność tłumienia drgań swobodnych w przypadku podkrytycznym. Jest to wielkość stała dla rozpaywanego                                  układu drgającego i nie zależy od warrunków początkowych. Zależy natomiast od maksymalnych wychyleń w chwilach różniących się o okres drań tłumionych.
                           '''))
        display(ReportText(f'''
                           Logarytmiczny dekrement tłumienia u wyrażona jest wzorem przy pomocy okresu drgań tłumionych:
                           '''))
        display(SympyFormula( Eq(Symbol('delta'),2*pi*system.damping_coefficient()*(system.natural_frequencies())**-1,evaluate=False)))
        
        display(ReportText(f'''
                           Zaprezentowana zależność opisuje bezwymaiarową wielkość charakteryzującą intensywność tłumienia drgań swobodnych w przypadku podkrytycznym. Jest to wielkość stała dla rozpaywanego układu drgającego i nie zależy od warrunków początkowych.Zależy natomiast od okresu drań tłumionych i współczynnika h.
                           '''))


class LagrangianComponent(ReportComponent):
    
    title="Lagrangian (funkcja Lagrange'a)  układu"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        
        print(system._scheme())
        
        mrk_lagrangian_nonlin=Marker('lagrangianNL',prefix='eq')

        #display(ReportText(f'''The following model is considered. The system's Lagrangian is described by the formula ({Ref(mrk_lagrangian_nonlin).dumps()}):
        #                    '''))
        display(ReportText(f'''Lagrangian systemu wyrażony jest wzorem ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):
                            '''))

        display((SympyFormula(  Eq(Symbol('L'),dyn_sys.L.expand()[0])  , marker=mrk_lagrangian_nonlin )  ))
        
        q_sym =[ Symbol(f'{coord}'[0:-3]) for coord in dyn_sys.q]
        
        diffL_d=lambda coord: Symbol(f'\\frac{{ \\partial L}}{{  \\partial {vlatex(coord)}  }}')
        diffD_d=lambda coord: Symbol(f'\\frac{{ \\partial D}}{{  \\partial {vlatex(diff(coord))}  }}')
        d_dt_diffL_d=lambda coord: Symbol(f'\\frac{{ \\mathrm{{d}}  }}{{  \\mathrm{{d}} {vlatex(system.ivar)}  }} {vlatex(diffL_d(coord))} ')        

        display(ReportText(f'''Równania Eulera-Lagrange'a dla rozważanego przypadku są nastęujące: 
                            '''))
        
        for coord in dyn_sys.q:
            display((SympyFormula(  Eq(d_dt_diffL_d(coord.diff(system.ivar)) - diffL_d(coord) + diffD_d(coord),Symbol(f'Q_{{ {vlatex(coord)} }}^N'))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        
        
        display(ReportText(f'''Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: 
                            '''))
        
        for coord in dyn_sys.Y:
            display((SympyFormula(  Eq(diffL_d(coord),dyn_sys.L.expand()[0].diff(coord))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            


        for coord in dyn_sys.q.diff(system.ivar):
            display((SympyFormula(  Eq(d_dt_diffL_d(coord),dyn_sys.L.expand()[0].diff(coord).diff(system.ivar))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        for coord in dyn_sys.q:
            display((SympyFormula(  Eq(diffD_d(coord),(S.One/2 * diff(dyn_sys.q.transpose())* dyn_sys_lin.damping_matrix()* diff(dyn_sys.q))[0].diff(diff(coord)))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            #display(Markdown(f'\\begin{equation}    \\end{equation}').reported())
        #with doc_model.create(DMath()) as eq:
        #    eq.append(NoEscape(latex(Derivative(Symbol('L'),q_sym[0],evaluate=False))))
        #    eq.append(NoEscape('='))
        #    eq.append(NoEscape(vlatex(dyn_sys.L.expand()[0].diff(dyn_sys.q[0]))))

        mrk_gov_eq_nonlin=Marker('gov_eq_nonlin_sys',prefix='eq')

        #display(ReportText(f'''The governing equations of the system have a following form ({Ref(mrk_gov_eq_nonlin).dumps()}):
        #                    '''))

        display(ReportText(f'''Wyniki przedstawionych operacji wykorzystuje się wyznaczenia równań ruchu układu.
                            '''))


        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame
    
    
class GoverningEquationComponent(ReportComponent):
    
    title="Równania ruchu"
    
    def append_elements(self):
        
        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        

        
        if len(system.q)==1:
            markers_str=f"przedstawia zależność ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})"            
        else:
            markers_str=f"przedstawiają zależności ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))})"

        display(ReportText(f'''
                           Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru.
                           Równania ruchu układu (nieliniowe w ogólnym przypadku) {markers_str}:
                           '''))
        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=None))

        display(ReportText(f'''
                           Wyznaczone równania stanowią matematyczny opis dynamiczny właściwości układu.
                           Dalsza analiza pozwala na skuteczną analizę działania modelowanego obiektu i określenie jego parametrów mechanicznych.
                           '''))


class LinearizationComponent(ReportComponent):
    
    title="Linearyzacja równań ruchu"
    
    def append_elements(self):
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()

        coords=tuple(list(dyn_sys.Y) + list(dyn_sys.q.diff(t,t)))
        op_point = {coord: 0 for coord in coords}

        #display(self._op_points(hint=hint, subs=True))
        op_point.update(dyn_sys._op_points(subs=True)[0])


        mrk_lagrangian_nonlin = Marker('lagrangLin',prefix='eq')
        mrk_lagrangian_lin = Marker('lagrangLin',prefix='eq')

        display(ReportText(
                f'''Linearyzaja równań polega na znalezieniu ich rozwinięcia w szereg Taylora względem współrzędnych, prędkości i przyspieszeń uogólnionych w otoczeniu punktu równowagi.
                Celem uproszczenia wprowadzono następujące oznaczenia:'''))

        for coord in coords:
            display((SympyFormula(  Eq(Symbol(vlatex(coord)),Symbol(latex(coord))) , marker=mrk_lagrangian_lin  )  )) 



        display(ReportText(
                f'''Punkty równowagi rozważanego układu są następujące:
                            '''))          

        for eq_coord,val in op_point.items():
            display((SympyFormula(  Eq(eq_coord,val) , marker=mrk_lagrangian_lin  )  ))


        diffL_d=lambda coord: Symbol(latex(Derivative(Symbol('L'),Symbol(vlatex(coord))))  )



        for no,eom in enumerate(dyn_sys._eoms):




            eq_sym=Symbol(f'RR_{latex(dyn_sys.q[no])}')


            display(ReportText(f'''Równanie ruchu dla współrzędnej ${latex(dyn_sys.q[no])}$ można przestawić jako:
                                '''))

            display((SympyFormula(  Eq(eq_sym,eom,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))


            display(ReportText(
                f'''Formalnie należy obliczyć pochodne cząstkowe wielkości uogólnionych ze składników równań Lagrange'a:
                            '''))


            display((SympyFormula(  Eq(MultivariableTaylorSeries(eq_sym,coords,n=1,x0=op_point)._symbolic_sum(),0) , marker=None,backend=latex  )  ))

            diff_list=MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).calculation_steps(expr_symbol=eq_sym)

            display(ReportText(
                f'''Poszczególne pochodne mają następującą postać:
                            '''))

            for diff_eq in diff_list:

                display((SympyFormula(  diff_eq , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display(ReportText(f'''Po podstawieniu obliczonych pochodnych, otrzumuje się następujące zlinearyzowane równanie:
                                '''))
            display((SympyFormula(  Eq(MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).doit().expand().simplify().expand(),0,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))



        AutoBreak.latex_backend = latex_store


    
class FundamentalMatrixComponent(ReportComponent):
    
    title="Wyznaczanie macierzy fundamentalnej"
    
    def append_elements(self):
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(f'''Z równań ruchu wyznaczono macierz mas i sztywności układu:
                                '''))
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker='a' )  ))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker='a')  ))

        Delta = Symbol('\Delta')

        display(ReportText(f'''Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu ${latex(Delta)}$, przedstawiają się następująco:
                                '''))

        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker='a'  )  ))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker='a',backend=latex  )  ))

        display(ReportText(f'''Macierz fundamentalna pozwala określić rozwiązanie ustalone. Natomiast bazując na równaniu charakterystycznym określa się częstości własne układu.
                                '''))

        AutoBreak.latex_backend = latex_store

class GeneralSolutionComponent(ReportComponent):
    
    title="Rozwiązanie ogólne"
    
    def append_elements(self):

        from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(f'''Rozwiązanie ogólne przedstawia wyrażenie:
                                '''))
        display((SympyFormula(  Eq(Symbol('X'),HarmonicOscillator(dyn_sys_lin.linearized(
                                            )).general_solution().n(3),
                                            evaluate=False) , marker='a',backend=latex  )  ))

        display(ReportText(f'''Rozwiązanie ogólne opisuje ruch analizowanego układu (przedstawia przemieszczenie w funkcji czasu) i wynika z rozważań dotyczących drgań swobodnych układu.
                                '''))

        AutoBreak.latex_backend = latex_store

        
        
class FrequencyResponseFunctionComponent(ReportComponent):
    
    title="Charakterystyka Amplitudowo-Częstotliwościowa"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           funkcja odpowiedzi częstotliwościowej:
                           '''))

        display(SympyFormula( Eq(Symbol('frf'),
                     dyn_sys.frequency_response_function() ), marker=None))

        display(ReportText(f'''
                           jest to suma kwadratów amplitud pod pierwiastkiem
                           '''))  
        
FRFComponent = FrequencyResponseFunctionComponent

class FrequencyResponseFunctionComponentToSecond(ReportComponent):
    
    title="sens mocy układu"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           sens mocy układu:
                           '''))

        display(SympyFormula( Eq(Symbol('frf^2'),
                     dyn_sys.frequency_response_function().doit()**2 ), marker=None))

        display(ReportText(f'''
                           Charakterystyka Amplitudowo-Częstotliwościowa podniesiona do kwadratu
                           '''))  
        


        

class SteadySolutionComponent(ReportComponent):
    
    title="Rozwiązanie szczególne"
    _phi=False
    def append_elements(self,phi=_phi):

        from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator
        
        system = self._system
        ReportText.set_directory('./SDAresults')

        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex_store
        
        t=system.ivar
        

        dyn_sys=system
        dyn_sys_lin=dyn_sys.linearized()


        display(ReportText(f'''Rozwiązanie szczególne przedstawia wyrażenie:
                                '''))

        display((SympyFormula(  Eq(Symbol('X_s'),
                            HarmonicOscillator(dyn_sys_lin.linearized(
                            )).steady_solution().n(3),
                            evaluate=False) , marker='b',backend=latex  )  ))

        AutoBreak.latex_backend = latex_store

        display(ReportText(f'''Rozwiązanie szczególne związane jest obecnością wielkości wymuszających ruch (drgania) analizowanego układu.
                                '''))
        
class MaxStaticForce(ReportComponent):
    
    title="Maksymalna siła statyczna"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Wartość maksymalna siły statycznej działającej na pojedynczy element mocujący:
                           '''))

        display(SympyFormula( Eq(Symbol('F_s'),
                     dyn_sys.max_static_force().doit() ), marker=None))

        display(ReportText(f'''
                           Wartość maksymalna siły statycznej działającej na pojedynczy element mocujący.
                           '''))  
        
class MaxDynamicForce(ReportComponent):
    
    title="Maksymalna siła statyczna"

    def append_elements(self):

        system = self._system
        dyn_sys=system
        dyn_sys_lin = dyn_sys


        display(ReportText(f'''
                           Wartość maksymalna siły dynamicznej działającej na pojedynczy element mocujący:
                           '''))

        display(SympyFormula( Eq(Symbol('F_d'),
                     dyn_sys.max_dynamic_force().doit() ), marker=None))

        display(ReportText(f'''
                           Wartość maksymalna siły dynamicznej działającej na pojedynczy element mocujący:
                           ''')) 