from dynpy import core
from dynpy import dynamic

from sympy import *
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
from sympy.simplify.fu import TR8, TR10, TR7, TR3

import numpy as np
import itertools as itools
import scipy.integrate as solver
from .timeseries import *
from collections import ChainMap
from IPython.display import display

class OdeComputationalCase:
    '''
    This object allows for a fully numerical investigation on the dynamic system - by supplying methods such as formation of numerical right-hand sides of ordinary differential equations, preparing the input for scipy 'solve_ivp' integration function returned as a dictionary of numerical odes, initial conditions and integration method used, the object provides a comprehansive tool that can be utilised to determine mechanical system's behaviour numerically. Other methods are discussed in details further in this document.
    
    Arguments
    =========
    
    odes_system: Symbol object
        Ordinary differential equation in symbolic form
    
    ivar=None (optional): Symbol object
        Independent variable
    
    dvars: Symbol object
        Derivative symbol
    
    t_span: TimeSeries object
        Time span
    
    params: Symbol object
        
    
    params_values: float
        
    
    ic_point:
        
    
    evaluate=False (optional): bool
        Evaluate the ODE equation, False as a default
    
    label=None (optional): string
        Labels the instance. The default label is: '{Class name} with {length of dvars} equations'
    
    '''
    def __init__(self,
                 odes_system=[],
                 ivar=None,
                 dvars=[],
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 evaluate=False,
                 label=None):


        #if label==None:

        self.odes_system = odes_system
        self.ivar = ivar
        self.dvars = dvars
        self.ic_point = ic_point

        if isinstance(params, dict):
            self.params_values = params
            self.params = list(self.params_values.keys())

        elif (all(isinstance(elem, Symbol) for elem in params)):
            self.params = params
            self.params_values = params_values

        elif (all(
                isinstance(elem, tuple) and len(elem) == 2
                for elem in params)):
            self.params_values = {var: value for var, value in params}
            self.params = list(self.params_values.keys())

        else:
            self.params = self.odes_system.free_symbols
            self.params.remove(self.ivar)
            self.params_values = params_values

        self.t_span = t_span

        if evaluate:
            self.form_numerical_rhs()
        else:
            self.__numerical_odes = None

        if label == None:
            label = self._label = self.__class__.__name__ + ' with ' + str(
                len(self.dvars)) + ' equations'
        self._label = label

    def __call__(self, label=None):

        #         if len(args)>0:
        #             if isinstance(args[0],str):
        #                 label=args[0]
        #             else:
        #                 q=args

        #         self.qs=
        self._label = label
        return self

    def __str__(self):
        #         if self._label==None:
        #             self._label = self.__class__.__name__ + ' with ' + str(len(self.dvars)) + ' equations'
        #self._label = self.__class__.__name__ + ' with ' + str(len(self.dvars)) + ' equations'

        return self._label

    def __repr__(self):

        return self.__str__()

    def __fortran_odes_rhs(self):
        '''
        Generates the bininary code related to symbolical expresions declered in the __init__ method. The function object is returned where the all arguments create a tuple.
        '''
        subs_dict = {
            var: Symbol('temp_sym_' + str(i))
            for i, var in enumerate(self.dvars)
        }

        args_list = [self.ivar] + list(subs_dict.values()) + self.params

        return autowrap((self.odes_system.subs(subs_dict, simultaneous=True)),
                        args=args_list)
        #return autowrap(  (msubs(self.odes_system,subs_dict)),args=args_list)

    def form_numerical_rhs(self):
        '''
        Generates and returns the bininary code related to symbolical expresions declered in the __init__ method. Ready-to-use function object in the compact form f(t,y,params).
        '''
        odes_rhs = self.__fortran_odes_rhs()

        self.__numerical_odes = lambda t, y, *args, **kwargs: np.asarray(
            (odes_rhs(t, *y, *args, **kwargs))).reshape(y.shape)
        return self.__numerical_odes

    def solve_ivp_input(self,
                        t_span=None,
                        ic_list=None,
                        t_eval=None,
                        params_values=None,
                        method='RK45'):
        '''
        Returns the dictionary containing the necessary argument of solve_ivp integrator from scipy.integrate module.
        '''
        if type(ic_list) == type(None):
            ic_list = list(Matrix(self.dvars).subs(self.ic_point))

        if type(t_span) == type(None):
            t_span = self.t_span

        if type(params_values) == type(None):
            params_values = tuple(Matrix(self.params).subs(self.params_values))

        case_odes = self.__numerical_odes

        return {
            'fun': case_odes,
            't_span': [t_span[0], t_span[-1]],
            'y0': ic_list,
            't_eval': t_eval,
            'method': method,
            'args': params_values
        }

    def compute_solution(self,
                         t_span=None,
                         ic_list=None,
                         t_eval=None,
                         params_values=None,
                         method='RK45'):
        '''
        Returns the result of the computations of solve_ivp integrator from scipy.integrate module.
        '''
        solution = solver.solve_ivp(
            **self.solve_ivp_input(t_span=t_span,
                                   ic_list=ic_list,
                                   t_eval=t_eval,
                                   params_values=params_values,
                                   method=method))

        solution_tdf = TimeDataFrame(
            data={key: solution.y[no, :]
                  for no, key in enumerate(self.dvars)},
            index=t_span)

        solution_tdf.index.name = 't'

        return solution_tdf

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

        #display(list(nonlinear_system.linearized().lagrangian().atoms(Function)))

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

        #display(nonlin_lagrangian)

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
        print(self.args)
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

class DampedHarmonicOscillator(HarmonicOscillator):
    def solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
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
