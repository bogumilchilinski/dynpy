from sympy import Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq#, Tuple
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
#from sympy import *
#import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from sympy.utilities.lambdify import lambdify
from ..utilities.timeseries import TimeSeries, TimeDataFrame
from scipy.misc import derivative
from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3
from timer import timer



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
                 ic_point=None,
                 evaluate=False,
                 label=None,
                 backend='fortran',
                 ):


        #if label==None:
        self._backend = backend

        self.odes_system = odes_system
        self.ivar = ivar
        self.dvars = dvars
        self.ic_point = ic_point
        if type(ic_point) == type(None):
            self.ic_point = [0]*len(self.dvars)
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
    
    @property
    def ics_dvars(self):



        return self.dvars
    

    def __fortran_odes_rhs(self):
        '''
        Generates the bininary code related to symbolical expresions declered in the __init__ method. The function object is returned where the all arguments create a tuple.
        '''
        subs_dict = {
            var: Symbol('temp_sym_' + str(i))
            for i, var in enumerate(self.dvars)
        }

        args_list = [self.ivar] + list(subs_dict.values()) + self.params

        return autowrap(((self.odes_system).subs(subs_dict, simultaneous=True)),
                        args=args_list)


    def __numpy_odes_rhs(self):
        '''
        Generates the numpy code related to symbolical expresions declered in the __init__ method. The function object is returned where the all arguments create a tuple.
        '''
        subs_dict = {
            var: Symbol('temp_sym_' + str(i))
            for i, var in enumerate(self.dvars)
        }

        args_list = [self.ivar] + list(subs_dict.values()) + self.params

        return lambdify( args_list ,
                         ((self.odes_system).subs(subs_dict, simultaneous=True)).n(),
                         [{'sin':np.sin,'cos':np.cos},'numpy']
                        )
    
        # return autowrap(  (msubs(self.odes_system,subs_dict)),args=args_list)

    def form_numerical_rhs(self):
        '''
        Generates and returns the bininary code related to symbolical expresions declered in the __init__ method. Ready-to-use function object in the compact form f(t,y,params).
        '''
        if self._backend == 'numpy':
            odes_rhs = self.__numpy_odes_rhs()        
        else:
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
#         print('input parames')
#         display(params_values)
        if type(ic_list) == type(None):
            ic_list = list(Matrix(self.dvars).subs(self.ic_point))
        if len(ic_list) != len(self.dvars):
            raise IndexError('Number of initial conditions is not correct.')
        if type(t_span) == type(None):
            t_span = self.t_span

        if type(t_eval) == type(None) and type(t_span) != type(None) :
            t_eval = t_span
            
        if type(params_values) == type(None):
            
            self.params = list(self.odes_system.free_symbols)
            self.params.remove(self.ivar)
            
            print(self.params)
            print(self.params_values)
            
            params_values = tuple(Matrix(self.params).subs(self.params_values))
            
        if isinstance(params_values,dict):
            self.params_values = params_values
            self.params = list(self.params_values.keys())
            
            self.params = self.odes_system.free_symbols
            self.params.remove(self.ivar)
            
            self.params=list(self.params)
            
            
#             print('check')
#             display(Matrix(self.params))
#             print(*self.params)
#             print(*map(type,(self.params)))
#             display(params_values)
#             display(Tuple(*self.params).subs(params_values))
            
#             print([ par in  params_values  for par in self.params ])
            
#             print(*params_values.keys())
            
            params_values = tuple(Matrix(self.params).subs(params_values))
            

        case_odes = self.__numerical_odes
        
        print('params case')
        display(self.params)
        display(params_values)

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
        
        with timer() as t:
            solution = solver.solve_ivp(
                **self.solve_ivp_input(t_span=t_span,
                                       ic_list=ic_list,
                                       t_eval=t_eval,
                                       params_values=params_values,
                                       method=method))

            solution_tdf = TimeDataFrame(
                data={key: solution.y[no, :]
                      for no, key in enumerate(self.dvars)}, index=t_span)



            velocities = self.dvars[int(len(self.dvars)/2) :]
            for vel in velocities:
                solution_tdf[vel].to_numpy()
                gradient = np.gradient(solution_tdf[vel].to_numpy(),t_span)
            solution_tdf[vel.diff(self.ivar)] = gradient
            print('_'*100,t.elapse)
            comp_time=t.elapse



        solution_tdf._set_comp_time(comp_time)
        solution_tdf.index.name = self.ivar
        return solution_tdf


    
    