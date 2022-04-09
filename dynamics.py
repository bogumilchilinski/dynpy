from typing import Type
from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                    hessian, Function, flatten, Tuple, im, pi, latex,dsolve,solve,
                    fraction,factorial,Derivative, Integral,Expr,Subs, Mul, Add)

from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from .utilities.timeseries import TimeSeries, TimeDataFrame

from collections import ChainMap

from IPython.display import display, Image
import base64

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3

from .solvers.numerical import OdeComputationalCase

from .solvers.linear import LinearODESolution, FirstOrderODE,MultivariableTaylorSeries,FirstOrderODESystem

from .solvers.nonlinear import WeakNonlinearProblemSolution, MultiTimeScaleMethod


from pylatex import Document, Section, Subsection, Subsubsection, Itemize, Package, HorizontalSpace, Description, Marker, Ref, Marker, Figure, Command, NewPage, LargeText, HugeText, MediumText, Center
from pylatex.base_classes import Environment
from pylatex.section import Paragraph, Chapter
from pylatex.utils import italic, NoEscape

from .utilities.report import (SystemDynamicsAnalyzer,DMath,ReportText,SympyFormula, AutoBreak, PyVerbatim)
from .utilities.templates.document import *
from .utilities.components import mechanics as mech_comp


from .utilities.adaptable import AutoMarker
import inspect
import copy



    

def multivariable_taylor_series(expr, args, n=2, x0=None):
    '''
    Computes the multivariable Taylor series of expresion expr for the given order n in the neighbourhood of x0 point.

    Args:
        expr (:obj:Expr): Expresion describing the function under consideration.
        args (array_like): Set of the arguments.
        to be continued ...


    Examples:
        Taylor series expansion of the second order for the function $y+x**2$ is as follows:

        >>> multivariable_taylor_series(y+x**2,args=[x,y],n=2,x0=[0,0])
        y+x**2
    '''

    order_max = n
    op_point = x0

    args_shifted = {
        arg: arg - arg_shift
        for arg, arg_shift in op_point.items()
    }

    # def diff_orders(var_list, order_arg): return [
    #     args for tmp in range(order_arg)
    # ]

    #     term_tmp_list=sum([list(itools.combinations_with_replacement(var_list,ord_tmp)) for ord_tmp in range(1,4)],[])

    #     term_tmp_dict= {comp:Mul(*comp)/Mul(*[factorial(elem)  for elem  in Poly(Mul(*comp),*var_list).terms()[0][0]]).doit()  for comp in term_tmp_list}
    #     term_tmp_dict

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

    # print(diff_orders_list)

    return (sum([
        expr.diff(*args_tmp).subs(op_point).doit() * poly
        for args_tmp, poly in diff_orders_dict.items()
    ]) + expr.subs(op_point)).doit()



    


def scalar_fun_quadratic_form(expr, coordinates, op_point):
    '''Vector of deformation '''
    u = (Matrix(coordinates) -
         Matrix(coordinates).subs(op_point, simultaneous=True)).doit()

    constant_term = (expr.subs(op_point, simultaneous=True))
    linear_term = sum((Matrix([expr]).jacobian(coordinates).subs(
        op_point, simultaneous=True)) * u)
    quad_tem = sum(
        S.One / 2 * u.T *
        ((hessian(expr, coordinates).subs(op_point, simultaneous=True))) * u)

    return constant_term + linear_term + quad_tem


class DynamicSymbol(Function):
    
    def _repr_latex_(self):
        print(self.args)
        return super()._repr_latex_()[0:-2]
    



class DynamicDerivative(Derivative):
    
    def _repr_latex_(self):
        
        print(self.args)
        return super()._repr_latex_()[0:-2]
    
    

class LagrangesDynamicSystem(me.LagrangesMethod):
    '''Lagrange's method object

    The object generates equation of motion after passing Lagrangian and generalized coordinates. After initialization it can be performed several operations on the object to find desired answers.

    Arguments
    =========
    Lagrangian: Symbol object
        The lagrangian equation - subtraction of potential and kinematic energy

    qs=None (optional): dynamicsymbol object
        Generalized coordinates

    forcelist=None (optional): (tuples) (Point,Vector) object, (ReferenceFrame,Vector) object
        Forces acting on the dynamic system

    bodies=None (optional): Point object, RigitBody object
        Bodies represented as points or rigit bodies in reference frame

    frame=None (optional): ReferenceFrame object
        Reference frame of the dinamic system

    hol_coneqs=None (optional): array-like
        The holonomic constraint equations

    nonhol_coneqs=None (optional): array-like
        The nonholonomic constraint equations

    label=None (optional): string
        Label of the class instance. Default label: '{Class name} with {length of qs} DOF'

    ivar=None (optional): Symbol object
        Independent variable

    evaluate=True (optional):
        Evaluates the dinamic system

    Example
    =======

    Rod propped up with spring on the right side - constant k2 and spring and dumper on the left side - constants k and c.

    >>>t = symbols('t')
    >>>g,k,m1,m2,l,c = symbols('g k m1 m2 l c')
    >>>x1,x2 = dynamicsymbols('x_1 x_2')    # Generalized coordinates
    >>>val = {g:9.81,k:20,m1:2,m2:2,c:10}    # Assumed values for constants to performe numerical calculation
    >>>m = m1 + m2    # Mass of rod determination
    >>>h = (x1+x2)/2    # High determination for potential energy
    >>>phi = (x1-x2)/l    # Angle determination of rod rotation
    >>>I = S.One/12 *m*l**2    # Moment of inertia of rod determination
    >>>x = (x1+x2)/2    # Rod vertical displacement determination
    >>>T = S.Half*m*x.diff(t)**2 + S.Half*I*phi.diff(t)**2     # Kinetic Energy equation
    >>>V = S.Half*k*x1**2 + S.Half*2*k*x2**2 + m*g*h    # Potential Energy equation
    >>>L = T - V    # Lagrangian calculation
    >>>N=ReferenceFrame('N')    # Defining of reference frame for coordinate system
    >>>P1=Point('P_1')    # Defining point in space
    >>>P2=Point('P_2')
    >>>P1.set_vel(N,x1.diff(t)*N.x)    # Set velocity of point P1 in reference system N on axis x
    >>>P2.set_vel(N,x2.diff(t)*N.x)
    >>>forcelist = [(P1,-c*x1.diff(t)*N.x)]    # External forces determination and seting reference frame N and axis x
    >>>rod = dyn.LagrangesDynamicSystem(L,qs=[x1,x2],forcelist=forcelist,frame=N)    # Inicjalization of LagrangesDynamicSystem instance
    >>>rod_numerical = rod.numerized(parameter_values=val)

    Firstly there are defined symbols with SymPy Symbol class and generalized coordinates which are defined with dynamicsymbols method. Next the mass m, hight h, angle phi, and moment of inertia I are determined to make use of it in energys equations. 
    Subsequently potential energy V and kinetic enetgy T are evaluated and defined to calcualte Lagrangian L. When it has been done, one can start to create ReferenceFrame N. In the frame points P1 and P2 are created where there will be applied velocity on x axis in above case. Additonali external forces are assined to the variable forcelist - it is necessery to chose the point and axies. 
    Finaly one can define instance rod of class LagrangesDynamicSystem by making use of calculation and variables made previously. In addition there is possibility to create instant of class OdeComputationalCase (rod_numerical) that allows to perform numerical operations on the EOM

    '''
    _hint=[]
    _default_subs_method='direct'
    scheme_name = 'engine.png'
    real_name = 'engine_real.PNG'
    reportclass= ExampleTemplate

    @classmethod
    def _scheme(cls):

        path = __file__.replace('.py', '/images/') + cls.scheme_name
        
        path = './dynpy/models/images/' + cls.scheme_name

        return path

    @classmethod
    def _real_example(cls):

        path = __file__.replace('.py', '/images/') + cls.real_name

        path = './dynpy/models/images/' + cls.real_name
        
        return path

    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return Image(base64.b64decode(encoded_string))




    def __init__(self,
                 Lagrangian,
                 qs=None,
                 forcelist=None,
                 bodies=None,
                 frame=None,
                 hol_coneqs=None,
                 nonhol_coneqs=None,
                 label=None,
                 ivar=sym.Symbol('t'),
                 evaluate=True,
                 system = None):
        """
        Supply the following for the initialization of DynamicSystem in the same way as LagrangesMethod
        """

        self._kinetic_energy = 0
        self._potential_energy = 0
        self._dissipative_potential = 0

        if system:
            # print(system._kinetic_energy)
            Lagrangian=system
            system=None
            
            
            # self._kinetic_energy = Lagrangian._kinetic_energy
            # self._potential_energy = Lagrangian._potential_energy

        if isinstance(Lagrangian, LagrangesDynamicSystem):

            self._kinetic_energy = Lagrangian._kinetic_energy
            self._potential_energy = Lagrangian._potential_energy
            self._dissipative_potential = Lagrangian._dissipative_potential


        if isinstance(Lagrangian, me.LagrangesMethod):
#             print('standart init')
            
            bodies = Lagrangian._bodies
            frame = Lagrangian.inertial
            forcelist = Lagrangian.forcelist
            hol_coneqs = (Lagrangian._hol_coneqs)
            nonhol_coneqs = list(
                Lagrangian.coneqs)[len((Lagrangian._hol_coneqs)):]
            qs = Lagrangian.q
            Lagrangian = sum(Lagrangian._L)



        
        self.ivar = ivar
        #         self.forcelist=forcelist
        self.system = system
        self.frame = frame

        super().__init__(Lagrangian=Lagrangian,
                         qs=qs,
                         forcelist=forcelist,
                         bodies=bodies,
                         frame=frame,
                         hol_coneqs=hol_coneqs,
                         nonhol_coneqs=nonhol_coneqs)

            
        self.L = self._L

        if evaluate == True:
            self.__governing_equations = self.form_lagranges_equations()
        else:
            self.__governing_equations = None

        self.governing_equations = self.__governing_equations

        if label == None:
            label = self.__class__.__name__ + ' with ' + str(len(self.q)) + 'DOF'

        self._label = label
        self._given_data={}
        
        self._nonlinear_base_system=None

    @classmethod
    def from_system(cls, system):
        
        kwargs=system._kwargs()
        
        new_system=cls(**kwargs)
        
        new_system._kinetic_energy = system._kinetic_energy
        new_system._potential_energy = system._potential_energy
        new_system._dissipative_potential = system._dissipative_potential

        return new_system
        #LM=me.LagrangesMethod(Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs)
    
    def symbols_description(self):
        self.sym_desc_dict = {

            tuple(self.q): r'generalized coordinates of the system,',
            self.ivar: r'independent variable (time),',
        }

        return self.sym_desc_dict
    
    
    @property
    def Y(self):
        return Matrix(list(self.q) + list(self.q.diff(self.ivar)))

    
    
    
    def _kwargs(self):
        """
        Returns all keyword arguments that an instance has
        """
        return {
            'bodies': self.bodies,
            'frame': self.frame,
            'forcelist': self.forcelist,
            'hol_coneqs': self._hol_coneqs,
            'nonhol_coneqs': list(self.coneqs)[len((self._hol_coneqs)):],
            'qs': self.q,
            'Lagrangian': self.lagrangian(),
            'label': self._label,
            'ivar': self.ivar,
            'system': self.system,
        }

    def __add__(self, other):
        """
        Returns the sum of provided instances in form of one class instance
        """
        self_dict = self._kwargs()
        other_dict = other._kwargs()

        self_dict[
            'Lagrangian'] = self_dict['Lagrangian'] + other_dict['Lagrangian']

        self_dict['qs'] = list(self_dict['qs']) + list(
            coord
            for coord in other_dict['qs'] if not coord in self_dict['qs'])

#         print(self_dict['qs'])

        def list_build(x): return [x] if x else []

        def forcelist_build(x): return x if x else []

        self_dict['forcelist'] = forcelist_build(
            self_dict['forcelist']) + forcelist_build(other_dict['forcelist'])

        self_dict['bodies'] = list_build(self_dict['bodies']) + list_build(
            other_dict['bodies'])
        
        
        if not self_dict['frame']:
            self_dict['frame']=other_dict['frame']

        if not self_dict['frame']:
            self_dict['frame'] = other_dict['frame']
            
        systems_sum=LagrangesDynamicSystem(**self_dict)
        systems_sum._given_data={**other._given_data,**self._given_data}
        
        systems_sum._kinetic_energy = sum([energy for energy in [self._kinetic_energy,other._kinetic_energy] if energy is not None])
        systems_sum._potential_energy = sum([energy for energy in [self._potential_energy,other._potential_energy] if energy is not None])
        systems_sum._dissipative_potential = sum([energy for energy in [self._dissipative_potential,other._dissipative_potential] if energy is not None])
        
        # print(systems_sum._kinetic_energy)
        # print(systems_sum._potential_energy)

        return systems_sum

    def shranked(self, *args):
        """
        Returns class instance with reduced Degrees of Freedom
        """
        

        
        self_dict = self._kwargs()
        self_dict['qs'] = flatten(args)
        


        new_sys=LagrangesDynamicSystem(**self_dict)
        

        
        return type(self)(0,system=new_sys)

    def remove(self, *args):

        bounded_coordinates = sym.flatten(args)

        self_dict = self._kwargs()
        self_dict['qs'] = [
            coord for coord in self.q if coord not in bounded_coordinates
        ]

        return LagrangesDynamicSystem(**self_dict)

    def subs(self, *args, **kwargs):
        """
        Returns class instance with substituted numerical values
        """

        given_data=args[0]
        
        hol_coneqs = list(self._hol_coneqs)

        if 'method' in kwargs.keys():
            method = kwargs['method']
        else:
            method = self._default_subs_method

        if method == 'as_constrains':
            if len(args) == 2:
                args = ([(args[0], args[1])]),

            elif len(args) == 1:
                if isinstance(args[0], dict):
                    # print(args[0])
                    args = list(args[0].items()),

            constrains = [
                lhs - rhs for lhs, rhs in args[0] if any(
                    coord in (lhs - rhs).atoms(Function) for coord in self.q)
            ]
            args = ([(lhs, rhs) for lhs, rhs in args[0]
                     if not any(coord in (lhs - rhs).atoms(Function)
                                for coord in self.q)], )

            hol_coneqs = hol_coneqs + constrains



        old_points = [point for point, force in self.forcelist]
        new_forces = [
            force.subs(*args, **kwargs) for point, force in self.forcelist
        ]

        frame = self.frame

        lagrangian_subs = self.lagrangian().subs(*args, **kwargs)

        new_points = [me.Point(str(point)) for point in old_points]

        for new_point, old_point in zip(new_points, old_points):
            new_point.set_vel(self.frame,
                              old_point.vel(frame).subs(*args, **kwargs))


        forces_subs = list(zip(new_points, new_forces))

        nonhol_coneqs_subs = list(self.coneqs)[len(
            (self._hol_coneqs)):]  # ,ivar=self.ivar

        #print(forces_subs)
        #print(self.forcelist)

        # print(type(self))

        new_system=LagrangesDynamicSystem(Lagrangian=lagrangian_subs,
                          qs=self.q,
                          forcelist=forces_subs,
                          bodies=self._bodies,
                          frame=self.frame,
                          hol_coneqs=hol_coneqs,
                          nonhol_coneqs=nonhol_coneqs_subs,
                          label=f'{self._label} for {args} and {kwargs}'
                          )
        
        #display(new_system._eoms)
        
        #print('subs is ran for '+str(type(self)))
        
        new_sys = type(self)(0,system=new_system)#(f'{self._label} for {args} and {kwargs}')
        new_sys._label=f'{self._label} for {args} and {kwargs}'

        #print('self._kinetic_energy')
        #print(self._kinetic_energy)
        #print(self._potential_energy)
        #print(self._potential_energy)
        
        new_sys._kinetic_energy = (self._kinetic_energy*S.One).subs(*args, **kwargs)
        new_sys._potential_energy = (self._potential_energy*S.One).subs(*args, **kwargs)
#         if hasattr(self._dissipative_potential,'subs'):
        new_sys._dissipative_potential = (self._dissipative_potential*S.One).subs(*args, **kwargs) 
#         else:
#             new_sys._dissipative_potential = self._dissipative_potential
        

        new_sys._given_data=given_data
        new_sys._nonlinear_base_system = copy.copy(self._nonlinear_base_system)
        
        #print(new_sys)
        #display(new_system._eoms)
        #display(new_system.forcelist)
        return new_sys






    def __preview(self, expr, preview=None, preview_mode=None):
        '''
        Private method which processes preview printing. The method is or isn't used depending on user's bool 'preview' input.
        '''
        display(expr)
        pass



    def __call__(self, *args, label=None):
        """
        Returns the label of the object or class instance with reduced Degrees of Freedom.
        """

        
        if isinstance(args[0], str):
            if label:
                self._label=label
            else:
                self._label=args[0]
                
                
            return self
        
        else:
            #print(flatten(*args))
            return self.shranked(*args)
            


    def __str__(self):

        return self._label

    def __repr__(self):

        return self.__str__()

    def lagrangian(self):
        """
        Returns the system Lagrange function defined within LagrangeMethod instance. Defined as self._L = Matrix([sympify(Lagrangian)]) in Sympy documentation.
        """

        return sum(self.L)

    def rhs_eq(self):
        """
        Returns the right-hand side of the equations of motion in the form of equations matrix. Output is provided basing on LagrangesMethod object.
        """

        return Matrix([
            Eq(comp.diff(self.ivar), self.rhs()[no], evaluate=False)
            for no, comp in enumerate(list(self.Y))
        ])

    def equilibrium_equation(self, static_disp_dict=None):
        """
        Finds the equilibrium conditions of the considered problem based on the system governing equations stated in the class instance.
        """
        if static_disp_dict is None:
            static_disp_dict = {
                q_tmp:
                Symbol(
                    str(q_tmp).replace('(' + str(self.ivar) + ')', '{^s}')
                )
                for q_tmp in self.q
            }

        self.q_0 = static_disp_dict
        
        eq_eqns=self.governing_equations.subs(static_disp_dict)
        
        
        
        trig_comps=eq_eqns.atoms(sin,cos)
        #trig_comp={}
        

        return self.governing_equations.subs(static_disp_dict).subs({comp:0 for comp  in trig_comps if comp.has(self.ivar)})

    
    def static_load(self):
        """
        Finds the static load of the considered problem based on the system governing equations stated in the class instance.
        """
        eqns = self.equilibrium_equation()
        

        return eqns.subs({coord:0 for coord in self.q_0.values()})
    
    @property
    def report_components(self):
        
        comp_list=[
        mech_comp.TitlePageComponent,
        mech_comp.SchemeComponent,
        mech_comp.ExemplaryPictureComponent,
        mech_comp.KineticEnergyComponent,
        mech_comp.PotentialEnergyComponent,
        mech_comp.LagrangianComponent,
        mech_comp.GoverningEquationComponent,
        mech_comp.FundamentalMatrixComponent,
        mech_comp.GeneralSolutionComponent,
        mech_comp.SteadySolutionComponent,
            
            
        ]
        
        return comp_list
    
    @property
    def report(self):

      
        sys=self
        doc = ExampleTemplate()
        # doc.append(mech_comp.TitlePageComponent(sys))
        # doc.append(mech_comp.SchemeComponent(sys))
        # doc.append(mech_comp.ExemplaryPictureComponent(sys))
        # doc.append(mech_comp.KineticEnergyComponent(sys))
        # doc.append(mech_comp.PotentialEnergyComponent(sys))
        # doc.append(mech_comp.LagrangianComponent(sys))
        # doc.append(mech_comp.GoverningEquationComponent(sys))
        # doc.append(mech_comp.FundamentalMatrixComponent(sys))
        # doc.append(mech_comp.GeneralSolutionComponent(sys))
        # doc.append(mech_comp.SteadySolutionComponent(sys))
    
    
        return doc
    
    
    def calculations_steps(self,preview=True,system=None,code=False,documentclass=Document,lang='pl'):

        if lang=='pl':
            doc_model = self._calculations_steps_pl(preview=preview,system=system,code=code,documentclass=documentclass)
        else:
            doc_model = self._calculations_steps_en(preview=preview,system=system,code=code,documentclass=documentclass)
        
        
        return doc_model
    
    def _calculations_steps_pl(self,preview=True,system=None,code=False,documentclass=Document):

        
        
        doc_model = documentclass('model')

        doc_model.packages.append(Package('booktabs'))
        doc_model.packages.append(Package('float'))
        doc_model.packages.append(Package('standalone'))
        doc_model.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(doc_model)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(doc_model)
        
        if system is None:
            system = self


        doc_model.append(NewPage())
        doc_model.append(Section('Analiza dynamiczna układu drgającego',numbering=False))
            
        display(ReportText(f'''Ilustracja przedstawia rzeczywisty obiekt mechaniczny, będący przedmiotem modelowania i analizy dynamicznej.
                            '''))
        
        print(system._scheme())
        with doc_model.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example())
        
        with doc_model.create(Figure(position='H')) as fig:
            fig.add_image(system._scheme())
            
            
        if code:
            display(ReportText(f'''Rozpatrywany układ został opisany (sformalizowany) przy pomocy odpowiedniej klasy dziedziczącej po typie nadrzędnym "ComposedModel". Kod zaprezentowanej struktury jest następujący:
                                '''))

            doc_model.content_separator='\n'
            with doc_model.create(PyVerbatim()) as verb:
                print(inspect.getsource(system.__class__))
                verb.append(NoEscape(inspect.getsource(system.__class__))  )
        
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        
        print(self._scheme())
        
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
            
        d_dt_diffL_d=lambda coord: Symbol(latex(Derivative(diffL_d(coord)  , self.ivar ))  )

        for coord in dyn_sys.q.diff(self.ivar):
            display((SympyFormula(  Eq(d_dt_diffL_d(coord),dyn_sys.L.expand()[0].diff(coord).diff(self.ivar))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        
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

        #mrk_gov_eq_lin=Marker('gov_eq_lin_sys',prefix='eq')
        
        return doc_model
    
    def _calculations_steps_en(self,preview=True,system=None,code=False,documentclass=Document):

        
        
        doc_model = documentclass('model')

        doc_model.packages.append(Package('booktabs'))
        doc_model.packages.append(Package('float'))
        doc_model.packages.append(Package('standalone'))
        doc_model.packages.append(Package('siunitx'))

        
        

        ReportText.set_container(doc_model)
        ReportText.set_directory('./SDAresults')

        #LatexDataFrame.set_picture_mode(True)
        #LatexDataFrame.set_directory('./SDAresults')

        SympyFormula.set_container(doc_model)
        
        if system is None:
            system = self
            
        doc_model.append(Section('Analysis of dynamic system',numbering=False))
            
        display(ReportText(f'''The figure presents real object, which is subject of demonstrated analysis.
                            '''))
        
        print(system._scheme())
        with doc_model.create(Figure(position='H')) as fig:
            fig.add_image(system._real_example())
        
        with doc_model.create(Figure(position='H')) as fig:
            fig.add_image(system._scheme())
            
            
        if code:
            display(ReportText(f'''Considered object was described (formalized) with utilization of appropriate class. The class inherits from "ComposedModel" type. The code for this system is as follows:
                                '''))

            doc_model.content_separator='\n'
            with doc_model.create(PyVerbatim()) as verb:
                print(inspect.getsource(system.__class__))
                verb.append(NoEscape(inspect.getsource(system.__class__))  )
        
        dyn_sys=system
        dyn_sys_lin = dyn_sys.linearized()

        
        print(self._scheme())
        
        mrk_lagrangian_nonlin=Marker('lagrangianNL',prefix='eq')

        #display(ReportText(f'''The following model is considered. The system's Lagrangian is described by the formula ({Ref(mrk_lagrangian_nonlin).dumps()}):
        #                    '''))
        display(ReportText(f'''The Lagrangian of the system under consideration is described by formula ({AutoMarker(Eq(Symbol('L'),dyn_sys.L.expand()[0]))}):
                            '''))

        display((SympyFormula(  Eq(Symbol('L'),dyn_sys.L.expand()[0])  , marker=mrk_lagrangian_nonlin )  ))
        
        q_sym =[ Symbol(f'{coord}'[0:-3]) for coord in dyn_sys.q]
        
        diffL_d=lambda coord: Symbol(latex(Derivative(Symbol('L'),Symbol(vlatex(coord))))  )
        
        display(ReportText(f'''The next step is to calculate derivatives appering in the Euler-Lagranges equations. Subsequent derivatives are as follows: 
                            '''))
        
        for coord in dyn_sys.Y:
            display((SympyFormula(  Eq(diffL_d(coord),dyn_sys.L.expand()[0].diff(coord))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
            
        d_dt_diffL_d=lambda coord: Symbol(latex(Derivative(diffL_d(coord)  , self.ivar ))  )

        for coord in dyn_sys.q.diff(self.ivar):
            display((SympyFormula(  Eq(d_dt_diffL_d(coord),dyn_sys.L.expand()[0].diff(coord).diff(self.ivar))  , marker=mrk_lagrangian_nonlin,backend=vlatex )  ))
        
        #with doc_model.create(DMath()) as eq:
        #    eq.append(NoEscape(latex(Derivative(Symbol('L'),q_sym[0],evaluate=False))))
        #    eq.append(NoEscape('='))
        #    eq.append(NoEscape(vlatex(dyn_sys.L.expand()[0].diff(dyn_sys.q[0]))))

        mrk_gov_eq_nonlin=Marker('gov_eq_nonlin_sys',prefix='eq')

        #display(ReportText(f'''The governing equations of the system have a following form ({Ref(mrk_gov_eq_nonlin).dumps()}):
        #                    '''))

        display(ReportText(f'''
                           The computed derivatives were applied for calculations of equations of motion.
                           Obtained formulas are described by equations ({AutoMarker(Eq(dyn_sys._eoms[0].simplify().expand(),0))})-({AutoMarker(Eq(dyn_sys._eoms[-1].simplify().expand(),0))}):
                           '''))
        
        for eq in dyn_sys._eoms:
            display(SympyFormula( Eq(eq.simplify().expand(),0) , marker=mrk_gov_eq_nonlin ))

        #mrk_lagrangian_lin=Marker('lagrangian_lin_sys',prefix='eq')
        #
        #                        It enables to find a simplified solution in the neighborhood of the critical point.
        #                        Such an approach introduces some error but the solution has qualitative compability of the exact and linearized result. The simplified Lagrangian formula ({Ref(mrk_lagrangian_lin).dumps()}) is as follows:
        #                    '''))

        #display((SympyFormula(  Eq(Symbol('L'),dyn_sys_lin.L.expand()[0]) , marker=mrk_lagrangian_lin  )  ))

        #mrk_gov_eq_lin=Marker('gov_eq_lin_sys',prefix='eq')
        
        return doc_model
    
    
    def external_forces(self):
        """
        Returns Matrix with external forces
        """
        return self.governing_equations.subs(
            {gen_coord: 0
             for gen_coord in self.Y}).doit()

    def reactions(self):
        multipliers = self.solve_multipliers()

        return [self._eoms.jacobian([lam])*value for lam, value in multipliers.items()]

    def generalized_momentum(self, dict=True):

        if dict:
            momentum_dict = {
                q_tmp:
                Symbol(
                    'p_{'+str(q_tmp).replace('(' + str(self.ivar) + ')', '')+'}')
                for q_tmp in self.q
            }

            return {momentum_sym: self.lagrangian().diff(coord.diff(self.ivar)) for coord, momentum_sym in momentum_dict.items()}

    def qdot_from_p(self, dict=True):

        if dict:
            qdot = self.q.diff(self.ivar)

            mom_subs = solve(
                [lhs-rhs for lhs, rhs in self.generalized_momentum().items()], list(qdot), dict=True)

        return mom_subs

    def hamiltonian(self):

        if dict:
            momentum_dict = self.generalized_momentum()

        qdot = self.q.diff(self.ivar)

        mom_subs = self.qdot_from_p()

        ham_sum = sum(
            [coord*mom for coord, mom in zip(qdot, momentum_dict.keys())])

        return Eq(Symbol('H'), (ham_sum-self.lagrangian()).subs(mom_subs[0]))

    def default_ics(self,critical_point=False):
        return {coord:0 for coord in self.Y}
    
    def _op_points(self,
                   static_disp_dict=None,
                   dict=True,
                   subs=False,
                   hint=None, 
                   *args,
                   **kwargs):
        """
        Provides the interface for critical points evaluation (solves the equlibrium conditions and returns its roots).
        """
        
        if hint is None: hint = self._hint
            
        #display(hint)

        eqns_to_solve = self.equilibrium_equation(
            static_disp_dict=static_disp_dict).doit()
        #print('op point - new')
        #display(eqns_to_solve )
        
        hint= [eq.subs(self.q_0) for eq  in hint]
        
        combined_eqns=list(eqns_to_solve) + list(flatten([hint]))
        #display(combined_eqns)
                                                 
        
        roots = solve(combined_eqns,
                      list(self.q_0.values()),
                      dict=dict)
        
        #display(roots)

        if subs:
            roots = [{
                lhs: rhs
                for lhs, rhs in zip(self.q, root_dict.values())
            } for root_dict in roots]

        #         if type(roots) is dict:
        #             roots = list(roots.values())

        #         #print(roots)
        #         roots = list(Matrix(roots).doit())

        return roots

    def critical_points(self, hint=None, static_disp_dict=None, dict=True):
        '''
        Provides the interface for critical points evaluation (solves the equlibrium conditions and returns its roots).
        '''

        roots_list = self._op_points(hint=hint,
                                     static_disp_dict=static_disp_dict,
                                     dict=dict)

        return [[Eq(coord, solution) for coord, solution in root_dict.items()]
                for no, root_dict in enumerate(roots_list)]

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        return self.mass_matrix

    def system_parameters(self, parameter_values=None):
        '''
        Recognises system parameters as symbols which are not independent variable or its relations and returns the Tuple (Sympy object) containing these elements. It does not take into account undefined functions (instances of Function class) of independent variable.
        '''
#         print('system parameters @')
#         print(self)
        
#         print(self.lagrangian())
#         print(self.q)
        
        params = self.rhs().free_symbols
        params.remove(self.ivar)
        if parameter_values == None:
            return list(params)
        else:
            return {
                param: parameter_values[no]
                for no, param in enumerate(params)
            }

    def computational_case(self, parameter_values=None):
        ''''
        Returns the instance of ComputatutionalOdeCase related to described within this class system.
        '''
        params = self.system_parameters(parameter_values=parameter_values)
        return {
            'odes_system': self.rhs().doit(),
            'ivar': self.ivar,
            'dvars': self.Y,
            'params': params,
            'label': 'Numerical model of ' + self.__str__(),
        }

    
    
    
    def solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        pass

    def approximated(self, n=3, x0=None, op_point=False, hint=[], label=None):
        """
        Returns approximated N-th order function calculated with Taylor series method as an instance of the class
        """

        # print('x0',x0)
        if not x0:
            x0 = {coord: 0 for coord in self.Y}

        #display(self._op_points(hint=hint, subs=True))
        if op_point:
            x0.update(self._op_points(hint=hint, subs=True)[0])
            #print('current op')
            #display(self._op_points(hint=hint, subs=True))

        lagrangian_approx = multivariable_taylor_series(self.lagrangian(),
                                                        self.Y,
                                                        n=n + 1,
                                                        x0=x0)
        
        approx_sys =LagrangesDynamicSystem(lagrangian_approx,
                                      self.q,
                                      forcelist=self.forcelist,
                                      frame=self.frame,
                                      label=label,
                                      ivar=self.ivar)

        approx_sys._nonlinear_base_system = self
        
        return approx_sys


    def linearized(self, x0=None, op_point=False, hint=[], label=None):
        """
        Returns approximated first order function calculated with Taylor series method as an instance of the class. It enables to obtain linearized output.
        Arguments:
        =========
            System = Created system based on symbolical represent of mechanical parts of it
            
            op_point - boolean, which points out if the operating point will be evaluated
            
            x0 - setting operating point
            
            hint - (optional) Adds additional equation to equilibrium condition and calculate op_point as equilibrium system.
            
            label=None (optional): string
                Label of the class instance. Default label: '{Class name} with {length of qs} DOF'

        Example:
        =======
        Creating the examplary system. A mass oscillating up and down while being held up by a spring with a spring constant kinematicly 

        >>> t = symbols('t')
        >>> m, g, l = symbols('m, g, l')
        >>> qs = dynamicsymbols('varphi') 
        >>> Pendulum()

        Creating linerized system in symbolic pattern
        >>> System_linearized = Sytem.linearized()

        """
        linearized_sys = self.approximated(n=1,
                                           op_point=op_point,
                                           hint=hint,
                                           x0=x0)

        
        lin_sys = LinearDynamicSystem(linearized_sys.lagrangian(),
                                   self.q,
                                   forcelist=self.forcelist,
                                   frame=self.frame,
                                   label=label,
                                   ivar=self.ivar)
        
        lin_sys._nonlinear_base_system = self
        
        return lin_sys

    @property
    def _eoms(self):
        """
        Returns Equations of Motion of the system
        """
        return self.governing_equations
    
    def _to_acc(self,expand=True):
        """
        Returns Equations of Motion of the system
        """
        self_dict = self._kwargs()
        self_dict['Lagrangian'] = self.lagrangian()/self.inertia_matrix()[0]
        if expand:
            self_dict['Lagrangian']=        self_dict['Lagrangian'].expand()


        new_sys=LagrangesDynamicSystem(**self_dict)
        

        
        return type(self)(0,system=new_sys)
    

    def numerized(self, parameter_values={}, FFT = None,label=None,backend='fortran',**kwargs):
        '''
        Takes values of parameters, substitute it into the list of parameters and changes list it into a Tuple. Returns instance of class OdeComputationalCase.
        Arguments:
        =========
            System = Created system based on symbolical represent of mechanical parts of it

        Example:
        =======
        Creating the examplary system. A mass oscillating up and down while being held up by a spring with a spring constant k

        >>> t = symbols('t')
        >>> m, k = symbols('m, k')
        >>> qs = dynamicsymbols('z') 
        >>> System = SDoFHarmonicOscillator(m,k, qs=[z]) 

        Defining the list of values to substitute
        >>> val ={
                M: 1000,
                k1: 1000
            }
        >>> System numeric = Sytem.numerized(parameters_values = val)

        - In default returned numerized system is in the time domain but can be represented in the frequency domain if it is desired
        >>> System numeric = Sytem.numerized(parameters_values = val, FFT = True)

        - if necessary the created numerized system can be solved in order to represent displacement, velocity, or acceleration 
        '''
        if not FFT:
            data_Tuple = Tuple(*self.system_parameters()).subs(parameter_values)
            computed_case = self.computational_case(parameter_values=data_Tuple)
            if label:
                print('dynpy label',label)
                computed_case['label']=label
                print('computed case',computed_case)

            return OdeComputationalCase(**computed_case,backend=backend, evaluate=True)







class LinearDynamicSystem(LagrangesDynamicSystem):
    '''

    '''


    def calculations_steps(self,preview=True,system=None,code=False,documentclass=Document,lang='pl'):

        if lang=='pl':
            doc_model = self._calculations_steps_pl(preview=preview,system=system,code=code,documentclass=documentclass)
        else:
            doc_model = self._calculations_steps_en(preview=preview,system=system,code=code,documentclass=documentclass)
        
        
        return doc_model        
    
    
    def _calculations_steps_pl(self,preview=True,system=None,code=False,documentclass=Document):
        
        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex
        
        t=self.ivar
        
        if self._nonlinear_base_system is None:
            doc_model=super()._calculations_steps_pl(preview=preview,code=code,documentclass=documentclass)
            
        else:
            doc_model= super()._calculations_steps_pl(preview=preview,system=self._nonlinear_base_system,code=code,documentclass=documentclass)


            
        
        
            if system is None:
                system = self

            
                
            dyn_sys=self._nonlinear_base_system
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

            display(ReportText(f'''Rozwiązanie równania charakterystycznego (dwukwadratowego) pozwala obliczyć częstości drgań własnych układu:
                                    '''))
            for no,omega in enumerate([omega for omega in HarmonicOscillator(dyn_sys_lin).natural_frequencies().doit() if omega !=0]):
                display((SympyFormula( Eq(Symbol(f'omega_0{no+1}'),omega) , marker=mrk_lagrangian_lin,backend=latex  )  ))

        AutoBreak.latex_backend = latex_store
        return doc_model
    

    
    def _calculations_steps_en(self,preview=True,system=None,code=False,documentclass=Document):
        
        latex_store=AutoBreak.latex_backend
        AutoBreak.latex_backend = latex
        
        t=self.ivar
        
        if self._nonlinear_base_system is None:
            doc_model=super()._calculations_steps_en(preview=preview,code=code,documentclass=documentclass)
            
        else:
            doc_model= super()._calculations_steps_en(preview=preview,system=self._nonlinear_base_system,code=code,documentclass=documentclass)


            
        
        
            if system is None:
                system = self

            
                
            dyn_sys=self._nonlinear_base_system
            dyn_sys_lin=dyn_sys.linearized()
            
            coords=tuple(list(dyn_sys.Y) + list(dyn_sys.q.diff(t,t)))
            op_point = {coord: 0 for coord in coords}

            #display(self._op_points(hint=hint, subs=True))
            op_point.update(dyn_sys._op_points(subs=True)[0])


            mrk_lagrangian_nonlin = Marker('lagrangLin',prefix='eq')
            mrk_lagrangian_lin = Marker('lagrangLin',prefix='eq')
            
            display(ReportText(
                    f'''Linearization of governing equations is based on Taylor series expansion for generalized coordinates,velocities and accelerations in the neighborhood of the stationary point.
                    The following conditions were introduced, in order to simplify equations of motions:'''))
            
            for coord in coords:
                display((SympyFormula(  Eq(Symbol(vlatex(coord)),Symbol(latex(coord))) , marker=mrk_lagrangian_lin  )  )) 
            
            
            
            display(ReportText(
                    f'''Equilibrium points of the system are as follows:
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
                
                
                display(ReportText(f'''Equation of motion for the coordinate ${latex(dyn_sys.q[no])}$ has a following form:
                                    '''))
                
                display((SympyFormula(  Eq(eq_sym,eom,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

                
                display(ReportText(
                    f'''The partial derivatives of Lagranges equations components  have to be calculated, in order to obtained linearized equations:
                                '''))


                display((SympyFormula(  Eq(MultivariableTaylorSeries(eq_sym,coords,n=1,x0=op_point)._symbolic_sum(),0) , marker=None,backend=latex  )  ))
                
                diff_list=MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).calculation_steps(expr_symbol=eq_sym)
                
                display(ReportText(
                    f'''The derivatives have a following form:
                                '''))
                
                for diff_eq in diff_list:
                
                    display((SympyFormula(  diff_eq , marker=mrk_lagrangian_lin,backend=latex  )  ))
                    
                display(ReportText(f'''The linearized equation is obtained by the substition of calculated derivatives:
                                    '''))
                display((SympyFormula(  Eq(MultivariableTaylorSeries(eom,coords,n=1,x0=op_point).doit().expand().simplify().expand(),0,evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))
                
            display(ReportText(f'''The inertia matrix $M$ and stiffness matrix $K$ are as follows:
                                    '''))
            display((SympyFormula(  Eq(Symbol('M'),dyn_sys_lin.inertia_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display((SympyFormula(  Eq(Symbol('K'),dyn_sys_lin.stiffness_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))
            
            Delta = Symbol('\Delta')
            
            display(ReportText(f'''Fundamental matrix (needed to obtain characteristic equation ${latex(Delta)}$), has a following representation:
                                    '''))

            display((SympyFormula(  Eq(Symbol('A'),dyn_sys_lin.fundamental_matrix(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))
            display((SympyFormula(  Eq(Delta,dyn_sys_lin.fundamental_matrix().det().expand().simplify().simplify().expand(),evaluate=False) , marker=mrk_lagrangian_lin,backend=latex  )  ))

            display(ReportText(f'''The solution of characteristic equation (biquatradic polynomial) enables to determine natural frequencies of the system:
                                    '''))
            for no,omega in enumerate([omega for omega in HarmonicOscillator(dyn_sys_lin).natural_frequencies().doit() if omega !=0]):
                display((SympyFormula( Eq(Symbol(f'omega_0{no+1}'),omega.expand().simplify()) , marker=mrk_lagrangian_lin,backend=latex  )  ))

        AutoBreak.latex_backend = latex_store
        return doc_model
    
    
    
    def stiffness_matrix(self):
        '''
        Returns the system stiffness matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized coordinates vector.
        '''
        return self.governing_equations.jacobian(self.q)

    def fundamental_matrix(self, freq=Symbol('omega', positive=True)):
        '''
        Method returns a fundamental matrix of the system built from inertia and stiffness matrices. Takes one optional argument.

        Args:
            freq (optional, obj:Symbol): This argument lets user to choose a symbol for frequency representation. Set to 'omega' by default, 'positive=True' ensures its affiliation to real numbers domain.
        '''
        return -freq**2 * self.inertia_matrix() +freq *sym.I * self.damping_matrix()   + self.stiffness_matrix()


    def damping_matrix(self):
        '''
        Returns the system damping matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized velocities vector.
        '''
        return self.governing_equations.jacobian(self.u)

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        main_matrix = self.rhs().jacobian(self.Y)

        return (main_matrix).diagonalize()[1]

    def modes(self):
        '''
        Determines the system vibration modes matrix (eigenmodes) based on the system free vibration frequencies.
        '''
        #         return (((self.inertia_matrix().inv() *
        #                   self.stiffness_matrix()).diagonalize()[0]))

        main_matrix = self.rhs().jacobian(self.Y)

        return (main_matrix).diagonalize()[0][len(self.q):, :]


class HarmonicOscillator(LinearDynamicSystem):
    """
    This object allows for a determination of any dynamic system by providing methods that serve to generate the equations of motion, solution, 
    natural frequencies, eigenmodes, FRF of the mechanical system and many others that are discussed further within the documentation. 
    For the initialization of HarmonicOscillator at least one argument is neccessary.

    Arguments
    =========
    Lagrangian: Symbol object
        The lagrangian equation - subtraction of potential and kinematic energy

    qs=None (optional): dynamicsymbol object
        Generalized coordinates

    forcelist=None (optional): (tuples) (Point,Vector) object, (ReferenceFrame,Vector) object
        Forces acting on the dynamic system

    bodies=None (optional): Point object, RigitBody object
        Bodies represented as points or rigit bodies in reference frame

    frame=None (optional): ReferenceFrame object
        Reference frame of the dinamic system

    hol_coneqs=None (optional): array-like
        The holonomic constraint equations

    nonhol_coneqs=None (optional): array-like
        The nonholonomic constraint equations

    label=None (optional): string
        Label of the class instance. Default label: '{Class name} with {length of qs} DOF'

    ivar=None (optional): Symbol object
        Independent variable

    evaluate=True (optional):
        Evaluates the dinamic system

    Example
    =======
    A mass oscillating up and down while being held up by a spring with a spring constant k

    >>> t = symbols('t')
    >>> m, k = symbols('m, k')
    >>> qs = dynamicsymbols('z') # Generalized Coordinates 
    >>> T = S.Half*m*z.diff(t)**2 # Kinetic Energy 
    >>> V = S.Half*k*z**2 # Potential Energy 
    >>> L = T - V # Lagrangian Calculation
    >>> N = ReferenceFrame('N') # Defining of reference frame for coordinate system
    >>> P = Point('P') # Defining point in space
    >>> P.set_vel(N, z.diff(t)*N.y) # Set velocity of point P in reference system N on axis z
    >>> Forcelist = [(P,f*sin(omega*t)*N.y)] # external forces on the system 
    >>> mass = dyn.HarmonicOscillator(dyn.LagrangesDynamicSystem(L, qs=[z], frame=N)) # Initialization of LagrangesDynamicSystem instance

    -We define the symbols and dynamicsymbols
    -Kinetic energy T and potential energy v are evaluated to calculate the lagrangian L
    -Reference frame was created with point P defining the position and the velocity determined on the z axis
    -external forces assigned 
    -finally we determine the instance of the system using class LagrangeDynamicSystem
    
    -damped natural frequencies, im eigenvals, and natual frequencies determine the system damped natural frequencies matrix, output is obtained from inertia matrix and stiffness matrix.
    -__solve solves the problem in a symbolic way. 
    -general solution gives the symbolic general solution and returns the matrix of the solution. 
    -steady solution computes the steady solution amplitude for the system defined in the instance of this class.
    -frequency response function returns the FRF of the system for the given excitation amplitude
    -dynamic amplification factor returns the DAF of the system for the given excitation amplitude
    -critical frequencies determines the critical frequency of the system
    -cycles numer determines the cycles number fopr the given time
    -exitation amplitude computes the excitation amplitude causing steady solution of the system with the given level of vibration
    -spring force determines the force in the elastic connector with the utilization of the analytical solution (it has to be done first)
    -DampedHarmonicOscillator solves the natrual frequency and vibration modes matrix
    """

    def damped_natural_frequencies(self):
        '''
        Determines the system damped natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        damped_freqs = [sqrt(nat_freq**2 - (self.damping_coefficient()[0]/2)**2)
                        for nat_freq in self.natural_frequencies() if nat_freq]

        return diag(*damped_freqs)

    def im_eigenvals(self):
        '''
        Determines the system damped natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        natural_freqs = list({(im((eigen).doit()))
                             for eigen in self.eigenvalues() if not eigen == 0})

        return diag(*natural_freqs)

    def natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        return (((self.inertia_matrix().inv() *
                  self.stiffness_matrix()).diagonalize()[1])).applyfunc(sqrt)

    def damping_coefficient(self):

        return ((self.inertia_matrix().inv() *
                 self.damping_matrix()))

    def logarithmic_decrement(self):

        log_dec = [(2*pi*self.damping_coefficient()[0]/2) /
                   damp_freq for damp_freq in self.damped_natural_frequencies()]

        return log_dec

    def canonical_governing_equation(
        self,
        natural_freq=None,
        damping_const=None,
        damped_freq=None,
    ):
        '''
        Does nothing and prints 'Done'.
        '''
        if natural_freq == None:
            self.omega0 = Symbol('omega_0', positive=True)

        if damping_const == None:
            self.h = Symbol('h', positive=True)

        if damped_freq == None:
            self.omegah = Symbol('omega_h', positive=True)
        
        self._canonical_governing_equation = Matrix(self.u).diff(
            self.ivar) + 2 * self.h * Matrix(self.u) + (self.omegah**2 +
                                                   self.h**2) * Matrix(self.q)
        print('Done')
        return self._canonical_governing_equation

    def __solve(self, initial_conditions=None):
        '''
        Private method solving the problem in the symbolic way and sets the related attribute.
        '''
        if len(self.q) == 1:
            solution = dsolve(sum(self.__governing_equations),
                              sum(self.q),
                              ics=initial_conditions)
        else:
            solution = dsolve((self.__governing_equations), (self.q),
                              ics=initial_conditions)

        return solution

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and returns matrix of solution (in the form of equations (objects of Eq class)).
        '''

        eoms = self._eoms
        subs_dict = {}
        if len(self.q) == 1:
            eoms = eoms-self.stiffness_matrix()*Matrix(self.q) + ((Symbol('omega_h', positive=True)
                                                                   ** 2+(self.damping_coefficient()[0]/2)**2)*self.inertia_matrix()*Matrix(self.q))
            subs_dict = {Symbol('omega_h', positive=True): sqrt(
                self.natural_frequencies()[0]**2-(self.damping_coefficient()[0]/2)**2)}
#             print('len',len(self.q))
#             display(subs_dict)

        return LinearODESolution(eoms, ivar=self.ivar, dvars=self.q).general_solution(
            initial_conditions=initial_conditions).subs(subs_dict).doit()



    def steady_solution(self, initial_conditions=None):
        """

        """
        return LinearODESolution(self._eoms, ivar=self.ivar, dvars=self.q).steady_solution(
            initial_conditions=initial_conditions)

    def solution(self, initial_conditions=None):

        return self.general_solution(initial_conditions)+self.steady_solution(initial_conditions)

#         eoms = self._eoms
#         subs_dict={}
#         if len(self.q)==1:
#             eoms = eoms-self.stiffness_matrix()*Matrix(self.q)+ ((Symbol('omega_h',positive=True)**2+(self.damping_coefficient()[0]/2)**2)*self.inertia_matrix()*Matrix(self.q) )
#             subs_dict={Symbol('omega_h',positive=True):sqrt(self.natural_frequencies()[0]**2-(self.damping_coefficient()[0]/2)**2)}
#             print('len',len(self.q))
#             display(subs_dict)

#         return LinearODESolution(eoms,ivar=self.ivar,dvars=self.q).genral_solution(
#             initial_conditions=initial_conditions).subs(subs_dict).doit()

    def steady_solution_amp(
            self,
            cos_amp=None,
            sin_amp=None,
            excitation_freq=None,
    ):
        ''''
        Computes the steady solution amplitude for the system defined in the instance of this class.
        '''

        self.Omega = excitation_freq
        omg = excitation_freq

        fund_mat = -self.inertia_matrix(
        ) * omg**2 + sym.I * omg * self.damping_matrix(
        ) + self.stiffness_matrix()
        steady_sol_amp = (fund_mat.inv() * Matrix(cos_amp)), (fund_mat.inv() *
                                                              Matrix(sin_amp))

        return steady_sol_amp

    def frequency_response_function(self,
                                    excitation_freq=Symbol('Omega',
                                                           positive=True)):
        '''
        Returns the Frequency Response Function of the system for the given excitation amplitude working correctly for systems with defined stiffenes matrix.
        '''


        self.Omega = excitation_freq
        omg = self.Omega

        #solution = self.steady_solution()[0].expand()
        sin_fun=list(self.external_forces().atoms(sin))
        cos_fun=list(self.external_forces().atoms(cos))

        if len(sin_fun) == 1: 
        
            omg_sin=list(set((sin_fun[0]).args)-{self.ivar})[0]

            comp_sin = self.external_forces().applyfunc(lambda comp: comp.coeff(sin_fun[0]))
        else:
            comp_sin=(self.external_forces()*S.Zero).doit()
            omg_sin=S.Zero
            
            
        if len(cos_fun) == 1: 
            omg_cos=list(set(cos_fun[0].args)-{self.ivar})[0]

            comp_cos = self.external_forces().applyfunc(lambda comp: comp.coeff(cos_fun[0]))
        else:
            comp_cos=(self.external_forces()*S.Zero).doit()
            omg_cos=S.Zero

        if omg_sin == 0:
            omg = omg_cos
        else:
            omg = omg_sin
            
        if excitation_freq is not None:
            omg = excitation_freq
        else:
            omg = omg.coeff(self.ivar)
            
        fund_mat = -self.inertia_matrix(
        ) * omg**2 + sym.I * omg * self.damping_matrix(
        ) + self.stiffness_matrix()

        amp=(fund_mat.inv() * (comp_sin)).T*(fund_mat.inv() * (comp_sin))  +  (fund_mat.inv() * (comp_cos)).T*(fund_mat.inv() * (comp_cos))

        return sqrt(amp[0])

    def dynamic_amplification_factor(self,
                                     excitation_amp=None,
                                     excitation_freq=Symbol('Omega', positive=True)):
        '''
        Returns the Dynamic Amplification Factor of the system for the given excitation amplitude (working correctly for single degree of freedom systems).
        '''
        frf = self.frequency_response_function(excitation_freq=excitation_freq)

#         display(self.external_forces())

        sin_comp = self.external_forces()[0].coeff(
            sin(excitation_freq * self.ivar))
        cos_comp = self.external_forces()[0].coeff(
            cos(excitation_freq * self.ivar))

#         display(sin_comp)

        maximal_force = sqrt(sin_comp**2 + cos_comp**2).subs(
            excitation_freq, 1).simplify()

        static_def = maximal_force / sum(self.stiffness_matrix())

        _dummy = symbols('_dummy', positive=True)

        nat_freq = self.natural_frequencies()[0]

        daf = (frf.subs(excitation_freq, _dummy * nat_freq) /
               static_def).simplify().subs(_dummy, excitation_freq / nat_freq)

        return daf

    def critical_frequencies(self, excitation_freq=None):
        '''
        Determines the critical frequency of the system (working correctly for single degree of freedom systems).
        '''
        frf_max_cond = (self.frequency_response_function(
            excitation_freq=excitation_freq)**2).diff(excitation_freq)
        # display(frf_max_cond)

        return solve(frf_max_cond, excitation_freq)

    def cycles_number(self, operational_time=None, excitation_freq=None):
        '''
        Determines the cycles number fopr the given time.
        '''
        excitation_period = 2 * pi / excitation_freq

        return operational_time / excitation_period

    def excitation_amplitude(self,
                             excitation_freq=None,
                             excitation_amp=None,
                             steady_vib_amp=None):
        '''
        Computes the excitation amplitude causing steady solution of the system with the given level of vibration.
        '''
        if steady_vib_amp == None:
            steady_vib_amp = Symbol('A', positive=True)

        self.A = steady_vib_amp

        #general_solution = self.solution().rhs

        exct_amp = (solve(
            self.frequency_response_function(excitation_freq=excitation_freq) -
            steady_vib_amp, excitation_amp))

        return exct_amp

    def spring_force(self,
                     spring_stiffness=None,
                     spring_position=None,
                     initial_conditions=None):
        ''''
        Determines the force in the elastic connector with the utilization of the analytical solution (it has to be done first).
        '''
        solution = self.solution(initial_conditions=initial_conditions).rhs

        return spring_stiffness * spring_position * solution

    def small_parameter(self, order=3):

        return self._eoms[0].diff(self.q[0], order).subs(self.q[0], 0)/factorial(order)


class DampedHarmonicOscillator(HarmonicOscillator):
    def solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and returns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        pass


class UndampedHarmonicOscillator(HarmonicOscillator):
    def natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        return (((self.inertia_matrix().inv() *
                  self.stiffness_matrix()).diagonalize()[1])).applyfunc(sqrt)

    def vibration_modes(self):
        '''
        Determines the system vibration modes matrix (eigenmodes) based on the system free vibration frequencies.
        '''
        return (((self.inertia_matrix().inv() *
                  self.stiffness_matrix()).diagonalize()[0]))


class WeakNonlinearOscillator(HarmonicOscillator):
    def __init__(self,
                 Lagrangian,
                 qs=None,
                 forcelist=None,
                 bodies=None,
                 frame=None,
                 hol_coneqs=None,
                 nonhol_coneqs=None,
                 label=None,
                 order=4,
                 eps=sym.Symbol('varepsilon'),
                 ivar=sym.Symbol('t')):

        nonlinear_system = LagrangesDynamicSystem(Lagrangian=Lagrangian,
                                                  qs=qs,
                                                  forcelist=forcelist,
                                                  bodies=bodies,
                                                  frame=frame,
                                                  hol_coneqs=hol_coneqs,
                                                  nonhol_coneqs=nonhol_coneqs,
                                                  label=label,
                                                  ivar=ivar,
                                                  evaluate=False)

        # display(list(nonlinear_system.linearized().lagrangian().atoms(Function)))

        stationary_subs_dict = {
            func: 0
            for func in list(nonlinear_system.linearized().lagrangian().atoms(
                Function)) if func not in nonlinear_system.q
        }

        lin_lagrangian = nonlinear_system.linearized().lagrangian().subs(
            stationary_subs_dict).doit()

        nonlin_lagrangian = (
            (nonlinear_system.approximated(order).lagrangian() -
             lin_lagrangian) / eps).doit()

        # display(nonlin_lagrangian)

        self._eps = Symbol(latex(eps), positive=True)

        super().__init__(
            Lagrangian=lin_lagrangian - self._eps * nonlin_lagrangian,
            qs=nonlinear_system.q,
            forcelist=nonlinear_system.forcelist,
            bodies=nonlinear_system.bodies,
            frame=nonlinear_system.frame,
            hol_coneqs=nonlinear_system._hol_coneqs,
            nonhol_coneqs=list(
                nonlinear_system.coneqs)[len((nonlinear_system._hol_coneqs)):],
            ivar=ivar)

        self.order = order

    @property
    def eps(self):
        return self._eps

    def __str__(self):
        return str(self.order) + '-order approximated ' + super().__str__()

    def solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        pass

    def zeroth_approximation(self):

        subscripts_dict = {
            q_tmp: dynamicsymbols(str(q_tmp).replace('(t)', '') + str('_0'))
            for q_tmp in self.q
        }
        print(subscripts_dict)
        
        return HarmonicOscillator(
            Lagrangian=self.linearized().lagrangian().subs(
                self.eps, 0).subs(subscripts_dict).doit(),
            qs=self.q.subs(subscripts_dict),
            forcelist=self.forcelist,
            bodies=self._bodies,
            frame=self.frame,
            hol_coneqs=self._hol_coneqs,
            nonhol_coneqs=list(self.coneqs)[len((self._hol_coneqs)):],
            ivar=self.ivar)
