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

from ..report import *
from ..adaptable import *

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



class LagrangianComponent(Section):
    
    def __init__(self, system,title='Analiza dynamiczna układu drgającego', numbering=False, *, label=True, **kwargs):
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

        self._system=system
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)




        self.packages.append(Package('booktabs'))
        self.packages.append(Package('float'))
        self.packages.append(Package('standalone'))
        self.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(self)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(self)
        

        #self.append(Section('Analiza dynamiczna układu drgającego',numbering=False))
            

#         if code:
#             display(ReportText(f'''Rozpatrywany układ został opisany (sformalizowany) przy pomocy odpowiedniej klasy dziedziczącej po typie nadrzędnym "ComposedModel". Kod zaprezentowanej struktury jest następujący:
#                                 '''))

#             doc_model.content_separator='\n'
#             with doc_model.create(PyVerbatim()) as verb:
#                 print(inspect.getsource(system.__class__))
#                 verb.append(NoEscape(inspect.getsource(system.__class__))  )
        
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
        
        diffL_d=lambda coord: Symbol(latex(Derivative(Symbol('L'),Symbol(vlatex(coord))))  )
        
        display(ReportText(f'''Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: 
                            '''))
        
        for coord in dyn_sys.Y:
            display((SympyFormula(  Eq(diffL_d(coord),dyn_sys.L.expand()[0].diff(coord))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            
        d_dt_diffL_d=lambda coord: Symbol(latex(Derivative(diffL_d(coord)  , system.ivar ))  )

        for coord in dyn_sys.q.diff(system.ivar):
            display((SympyFormula(  Eq(d_dt_diffL_d(coord),dyn_sys.L.expand()[0].diff(coord).diff(system.ivar))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        
        #with doc_model.create(DMath()) as eq:
        #    eq.append(NoEscape(latex(Derivative(Symbol('L'),q_sym[0],evaluate=False))))
        #    eq.append(NoEscape('='))
        #    eq.append(NoEscape(vlatex(dyn_sys.L.expand()[0].diff(dyn_sys.q[0]))))

        mrk_gov_eq_nonlin=Marker('gov_eq_nonlin_sys',prefix='eq')

        #display(ReportText(f'''The governing equations of the system have a following form ({Ref(mrk_gov_eq_nonlin).dumps()}):
        #                    '''))

        display(ReportText(f'''
                           Wykorzystując obliczone pochodne, wyznacza się równania ruchu na podstawie odpowiedniego wzoru.
                           Równania ruchu układu (nieliniowe w ogólnym przypadku) przedstawiają zależności ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))}):
                           '''))
        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=mrk_gov_eq_nonlin ))

        #mrk_lagrangian_lin=Marker('lagrangian_lin_sys',prefix='eq')
        #
        #display(ReportText(f'''Linearized model of the system can be useful as an initial step of the analysis. 
        #                        It enables to find a simplified solution in the neighborhood of the critical point.
        #                        Such an approach introduces some error but the solution has qualitative compability of the exact and linearized result. The simplified Lagrangian formula ({Ref(mrk_lagrangian_lin).dumps()}) is as follows:
        #                    '''))

        #display((SympyFormula(  Eq(Symbol('L'),dyn_sys_lin.L.expand()[0]) , marker=mrk_lagrangian_lin  )  ))


        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame
    
class LinearizationComponent(Section):
    
    def __init__(self, system,title='Linearyzacja równań ruchu', numbering=False, *, label=True, **kwargs):
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

        self._system=system
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)




        self.packages.append(Package('booktabs'))
        self.packages.append(Package('float'))
        self.packages.append(Package('standalone'))
        self.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(self)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(self)
        
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



#             display(ReportText(f'''Kolejne pochodne wynikające z zastosowania równań Eulera-Lagrange'a są nastęujące: 
#                                    ({Ref(mrk_lagrangian_nonlin).dumps()}):
#                                 '''))

        #op_point = {coord  for coord in self.Y + list(self.q.diff(self.ivar))}



        #display((SympyFormula(  Eq(Symbol('L'),dyn_sys_lin.L.expand()[0]) , marker=mrk_lagrangian_lin  )  ))

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

        display(ReportText(f'''Z równań ruchu wyznaczono macierz mas i sztywności układu:
                                '''))
        display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

        display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

        Delta = Symbol('\Delta')

        display(ReportText(f'''Macierz fundamentalna, na podstawie której wyznaczono równanie charakterystyczne rozważanego układu ${latex(Delta)}$, przedstawiają się następująco:
                                '''))

        display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))
        display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

#         display(ReportText(f'''Rozwiązanie równania charakterystycznego (dwukwadratowego) pozwala obliczyć częstości drgań własnych układu:
#                                 '''))
#         for no,omega in enumerate([omega for omega in (dyn_sys_lin).natural_frequencies().doit() if omega !=0]):
#             display((SympyFormula( Eq(Symbol(f'omega_0{no+1}'),omega) , marker=mrk_lagrangian_lin,backend=latex  )  ))

        AutoBreak.latex_backend = latex_store

    
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +=list(self.packages)
        frame+=(list(self))
        return frame
    
class SchemeComponent(Section):
    
    def __init__(self, system,title='Przykładowy obiekt', numbering=False, *, label=True, **kwargs):
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

        self._system=system
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)




        self.packages.append(Package('booktabs'))
        self.packages.append(Package('float'))
        self.packages.append(Package('standalone'))
        self.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(self)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(self)
        

        #self.append(Section('Analiza dynamiczna układu drgającego',numbering=False))
            
#         display(ReportText(f'''Ilustracja przedstawia schemat rzeczywisty obiekt mechaniczny, będący przedmiotem modelowania i analizy dynamicznej.
#                             '''))
        display(ReportText(f'''Ilustracja przedstawia schemat rzeczywistego obiektu mechanicznego, będący przedmiotem modelowania i analizy dynamicznej.
                            '''))        
        print(system._scheme())
#         with self.create(Figure(position='H')) as fig:
#             fig.add_image(system._real_example())
        
        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._scheme(),width='6cm')
            
            



        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +=list(self.packages)
        frame+=(list(self))
        return frame
    
class ExemplaryPictureComponent(Section):
    
    def __init__(self, system,title='Schemat modelu dynamicznego', numbering=False, *, label=True, **kwargs):
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

        self._system=system
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)




        self.packages.append(Package('booktabs'))
        self.packages.append(Package('float'))
        self.packages.append(Package('standalone'))
        self.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(self)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(self)
        

        #self.append(Section('Analiza dynamiczna układu drgającego',numbering=False))
            
        display(ReportText(f'''Ilustracja przedstawia rzeczywisty obiekt mechaniczny, będący przedmiotem modelowania i analizy dynamicznej.
                            '''))
      
        print(system._scheme())
        with self.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example(),width='6cm')
        



        
    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +=list(self.packages)
        frame+=(list(self))
        return frame