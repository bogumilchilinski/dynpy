from sympy import Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Dict, ImmutableMatrix, latex#, Tuple
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
#from sympy import *
#import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from sympy.utilities.lambdify import lambdify
from ..utilities.adaptable import TimeSeries, TimeDataFrame,NumericalAnalysisDataFrame
from scipy.misc import derivative
from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3
import time
import pandas as pd


class OdeComputationalCase:
    '''
    This object allows for a fully numerical investigation on the dynamic system - by supplying methods such as formation of numerical right-hand sides of ordinary differential equations, preparing the input for scipy 'solve_ivp' integration function returned as a dictionary of numerical odes, initial conditions and integration method used, the object provides a comprehansive tool that can be utilised to determine mechanical system's behaviour numerically. Other methods are discussed in details further in this document.

    Arguments
    =========
    odes_system: Symbol object
        System of first order ordinary differential equations in symbolic form

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
    
    _cached_odes = {}
    
    
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
        self._default_ics = None
        self.ivar = ivar
        self.dvars = dvars
        self.ic_point = ic_point
        if type(ic_point) == type(None):
            self.ic_point = [0]*len(self.dvars)
        if isinstance(params, (dict, Dict)):
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
            self._evaluated=True
        else:
            self.__numerical_odes = None
            self._evaluated=False

        if label == None:
            label = self._label = self.__class__.__name__ + ' with ' + str(
                len(self.dvars)) + ' equations'
        self._label = label

    @property
    def parameters(self):
        return self.odes_system.free_symbols-{self.ivar}
        
    def default_ics(self,critical_point=False):
        
        if isinstance(self._default_ics,dict):
            ics_instance={coord:self._default_ics[coord] for coord in self.dvars if coord in self._default_ics}
            
            return {**{coord:0 for coord in self.dvars},**ics_instance}
        else:
            return {coord:0 for coord in self.dvars}
        
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

        if self._default_ics is not None:
            return self._label + f' with ics {self._default_ics}'
        else:
            return self._label + ' without ics'

    def __repr__(self):

        return self.__str__()

    def _repr_latex_(self):
        if self._default_ics is None:
            return f'${latex(self.odes_system)} for {latex(self.dvars)}$'
        else:
            return f'${latex(self.odes_system)} for {latex(self.dvars)} with {latex(self.default_ics())}$'

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
        
        ivar_temp = Symbol('temp_ivar')
        #subs_dict[self.ivar] = ivar_temp

        #Zmiana Franek
        args_list = [ivar_temp] + list(subs_dict.values()) + self.params
        #args_list = [ivar_temp] + list(subs_dict.values()) + list(self.params)

#         display(self.odes_system.subs(subs_dict, simultaneous=True))

        #Zmiana Franek
        return autowrap(((self.odes_system).subs({self.ivar:ivar_temp,**subs_dict}, simultaneous=True)),
                        args=args_list)
#         return autowrap(((self.odes_system).subs({self.ivar[0]:ivar_temp,**subs_dict}, simultaneous=True)),
#                         args=args_list)
        

    def __numpy_odes_rhs(self):
        '''
        Generates the numpy code related to symbolical expresions declered in the __init__ method. The function object is returned where the all arguments create a tuple.
        '''
        subs_dict = {
            var: Symbol('temp_sym_' + str(i))
            for i, var in enumerate(self.dvars)
        }
        #Zmiana Franek
        args_list = [self.ivar] + list(subs_dict.values()) + self.params
        #args_list = list(self.ivar) + list(subs_dict.values()) + list(self.params)

        return lambdify( args_list ,
                         ((self.odes_system).subs(subs_dict, simultaneous=True)).doit().n(),
                         [{'sin':np.sin,'cos':np.cos},'numpy']
                        )
    
        # return autowrap(  (msubs(self.odes_system,subs_dict)),args=args_list)

    def form_numerical_rhs(self):
        '''
        Generates and returns the bininary code related to symbolical expresions declered in the __init__ method. Ready-to-use function object in the compact form f(t,y,params).
        '''
        
        odes_key = (ImmutableMatrix(self.odes_system),tuple(self.dvars))
        if odes_key in self.__class__._cached_odes:
            odes_rhs = self.__class__._cached_odes[odes_key]
        else:
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

        #Zmiana Franek
        if ic_list is None:
        #if ic_list is None:
            #Zmiana Franek
            ic_list = list(Matrix(self.dvars).subs(self.ic_point))
            #ic_list = self.ic_point
        if len(ic_list) != len(self.dvars):
            raise IndexError('Number of initial conditions is not correct.')
        
        if type(t_span) == type(None):
            t_span = self.t_span

        if type(t_eval) == type(None) and type(t_span) != type(None) :
            t_eval = t_span
            
        if type(params_values) == type(None):
            
            #Zmiana Franek
            self.params = list(self.odes_system.free_symbols - {self.ivar,Symbol('t')})
            # self.params = list(self.odes_system.free_symbols - {self.ivar[0],Symbol('t')})
            
            print(self.params)
            print(self.params_values)
            
            params_values = tuple(Matrix(self.params).subs(self.params_values))
            
        if isinstance(params_values,dict):
            self.params_values = params_values
            self.params = list(self.params_values.keys())
            
            self.params = self.odes_system.free_symbols - {self.ivar,Symbol('t')}
            
            
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

        
    def with_ics(self, ics=None, ivar0=0, sol0=None):
        self._default_ics = ics
       
        return self
        
    def compute_solution(self,
                         t_span=None,
                         ic_list=None,
                         t_eval=None,
                         params_values=None,
                         method='RK45',
                         derivatives=False):
        '''
        Returns the result of the computations of solve_ivp integrator from scipy.integrate module.
        '''
        if ic_list is None:
            ic_list=self._default_ics

        if not self._evaluated:
            self.form_numerical_rhs()
            self._evaluated=True


        t_0 = time.time()

        
        velocities = self.dvars[int(len(self.dvars)/2) :]
        if len(t_span)>1:
            solution = solver.solve_ivp(
                **self.solve_ivp_input(t_span=t_span,
                                    ic_list=ic_list,
                                    t_eval=t_eval,
                                    params_values=params_values,
                                    method=method))

            solution_tdf = TimeDataFrame(
                data={key: solution.y[no, :]
                    for no, key in enumerate(self.dvars)}, index=t_span)



            
            for vel in velocities:
                solution_tdf[vel].to_numpy()
                gradient = np.gradient(solution_tdf[vel].to_numpy(),t_span)
                solution_tdf[vel.diff(self.ivar)] = gradient
        else:
            solution_tdf = TimeDataFrame(
                data={key: ic_list[no]
                    for no, key in enumerate(self.dvars)}, index=t_span)
            for vel in velocities:

                solution_tdf[vel.diff(self.ivar)] = 0.0
            
        t_e = time.time()
        t_d = t_e-t_0
        print('_'*100,t_d)
        comp_time=t_d

        solution_tdf._set_comp_time(comp_time)
        solution_tdf.index.name = self.ivar
        return solution_tdf

    def numerized(self,*args,**kwargs):

        return self
        
    def _as_na_df(self):
        
        coords = list(self.dvars)
        num_cases = pd.MultiIndex.from_product([[self],[Eq(Symbol('a'),10)], coords])

        if len(self._default_ics)==len(coords):
            ics = self._default_ics
        else:
            ics = None

        return NumericalAnalysisDataFrame(data=None,
                    index=pd.Index([0.0],name=self.ivar),
                    columns=num_cases,
                    )

   
