from sympy import *
from sympy.physics.mechanics import *
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from .timeseries import *

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8,TR10,TR7

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

    diff_orders = lambda var_list, order_arg: [
        args for tmp in range(order_arg)
    ]

    #     term_tmp_list=sum([list(itools.combinations_with_replacement(var_list,ord_tmp)) for ord_tmp in range(1,4)],[])

    #     term_tmp_dict= {comp:Mul(*comp)/Mul(*[factorial(elem)  for elem  in Poly(Mul(*comp),*var_list).terms()[0][0]]).doit()  for comp in term_tmp_list}
    #     term_tmp_dict

    diff_orders_list = sum([
        list(itools.combinations_with_replacement(args, order))
        for order in range(1, order_max + 1, 1)
    ], [])
    diff_orders_dict = {
        comp: (Mul(*comp).subs(args_shifted) / Mul(*[
            factorial(elem) for elem in Poly(Mul(*comp), *args).terms()[0][0]
        ])).doit()
        for comp in diff_orders_list
    }

    #print(diff_orders_list)

    return (sum([
        expr.diff(*args_tmp).subs(op_point).expand().doit() * poly
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


class OdeComputationalCase:
    '''
    This object allows for a fully numerical investigation on the dynamic system - by supplying methods such as formation of numerical right-hand sides of ordinary differential equations, preparing the input for scipy 'solve_ivp' integration function returned as a dictionary of numerical odes, initial conditions and integration method used, the object provides a comprehansive tool that can be utilised to determine mechanical system's behaviour numerically. Other methods are discussed in details further in this document.
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
        '''
        Supply the following arguments for the initialization of OdeComputationalCase:
        
        Args:
        '''

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


class LinearODESolution:
    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 equation_type=None):
        '''
        Supply the following arguments for the initialization of OdeComputationalCase:
        
        Args:
        '''
        if isinstance(odes_system, LinearDynamicSystem):
            odes_system = odes_system._eoms

        self.odes_system = odes_system

        self.governing_equations = self.odes_system
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
        self.__numerical_odes = None

        self.eq_type = equation_type

    def stiffness_matrix(self):
        '''
        Returns the system stiffness matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized coordinates vector.
        '''
        return self.governing_equations.jacobian(self.dvars)

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_ddot = list(sym.Matrix(self.dvars).diff(self.ivar, 2))

        return self.governing_equations.jacobian(dvars_ddot)

    def damping_matrix(self):
        '''
        Returns the system damping matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_dot = list(sym.Matrix(self.dvars).diff(self.ivar, 1))

        return self.governing_equations.jacobian(dvars_dot)

    def external_forces(self):
        return self.odes_system.subs(
            {gen_coord: 0
             for gen_coord in self.dvars}).doit()

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        C = numbered_symbols('C', start=1)

        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        Y_mat = Matrix(self.dvars)

        diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        t_sol = self.ivar

        solution = [
            next(C) * modes[:, i] * sin(sym.sqrt(eigs[i, i]) * t_sol) +
            next(C) * modes[:, i] * cos(sym.sqrt(eigs[i, i]) * t_sol)
            for i, coord in enumerate(self.dvars)
        ]

        return sum(solution, Matrix([0] * len(Y_mat)))

    def steady_solution(self, initial_conditions=None):

        ext_forces = self.external_forces()

        #         sin_components=ext_forces.atoms(sin)
        #         cos_components=ext_forces.atoms(cos)
        components = ext_forces.atoms(sin, cos)

        display(components)

        steady_sol = Matrix([0 for gen_coord in self.dvars])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
            display(omg)
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

            #display(amp_vector)

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (fund_mat.inv() * amp_vector) * comp

            ext_forces -= amp_vector * comp
        #print(ext_forces.doit().expand())
        if ext_forces.doit().expand() != sym.Matrix(
            [0 for gen_coord in self.dvars]):
            steady_sol += sym.dsolve(
                (self.governing_equations - self.external_forces() +
                 ext_forces).expand().doit(), self.dvars)

        return steady_sol

    def solution(self, initial_conditions=None):
        return self.general_solution(
            initial_conditions=initial_conditions) + self.steady_solution(
                initial_conditions=initial_conditions)


class WeakNonlinearProblemSolution(LinearODESolution):
    def __init__(self,
                 odes_system,
                 ivar=Symbol('t'),
                 dvars=[],
                 eps=Symbol('varepsilon'),
                 t_span=[],
                 params=[],
                 params_values={},
                 ic_point={},
                 equation_type=None):

        super().__init__(odes_system=odes_system,
                         ivar=ivar,
                         dvars=dvars,
                         t_span=t_span,
                         params=params,
                         params_values=params_values,
                         ic_point=ic_point,
                         equation_type=equation_type)
        self.eps = eps

    def approximation_function(self, order):

        dvars_no = len(self.dvars)

        return Matrix([
            me.dynamicsymbols('Y_' + str(dvar_no) + str(order))
            for dvar_no, dvar in enumerate(self.dvars)
        ])

    def predicted_solution(self, components_no=3,dict=False,equation=False):

        dvars_no = len(self.dvars)

        solution=sum((self.approximation_function(comp_ord) * self.eps**comp_ord
                    for comp_ord in range(components_no)),
                   sym.zeros(dvars_no, 1))

        if equation:
            solution=Matrix([Eq(lhs,rhs) for lhs,rhs in zip(self.dvars,list(solution))])
            
        if dict:
            solution={lhs:rhs for lhs,rhs in zip(self.dvars,list(solution))}
            
        return solution
        
        return 
    
    def zeroth_approximation(self,dict=False,equation=False):
        
        eoms=Matrix(self.odes_system).subs(self.eps,0)
        
        solution=LinearODESolution(eoms,ivar=self.ivar,dvars=self.dvars).solution()
        
        if equation:
            solution=Matrix([Eq(lhs,rhs) for lhs,rhs in zip(self.approximation_function(order=0),list(solution))])
            
        if dict:
            solution={lhs:rhs for lhs,rhs in zip(self.approximation_function(order=0),list(solution))}
            
        return solution


    def first_approximation(self,dict=False,equation=False):
        
        eoms=Matrix(self.odes_system).subs(self.predicted_solution(2,dict=True)).diff(self.eps).subs(self.eps,0).expand()
        
        eoms=(  eoms.subs(self.zeroth_approximation(dict=True)).expand() ).applyfunc(lambda eqn:  TR10(TR8(TR10( eqn ).expand()).expand()).expand())
        
        display(eoms)
        
        solution=LinearODESolution(eoms,ivar=self.ivar,dvars=self.approximation_function(order=1)).solution()
        
        if equation:
            solution=Matrix([Eq(lhs,rhs) for lhs,rhs in zip(self.approximation_function(order=0),list(solution))])
            
        if dict:
            solution={lhs:rhs for lhs,rhs in zip(self.approximation_function(order=0),list(solution))}
            
        return solution
    
class LagrangesDynamicSystem(me.LagrangesMethod):
    '''
    Arguments
    =========
    Lagrangian: object
        
    
    qs (optional): object: Symbols
        
    
    forcelist (optional): list(object)
        
    
    bodies (optional):
    
    
    frame (optional):
    
    
    hol_coneqs (optional):
    
    
    nonhol_coneqs (optional):
    
    
    label (optional):
    
    
    ivar (optional):
    
    
    evaluate (optional):
    
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
                 evaluate=True):
        """
        Supply the following for the initialization of DynamicSystem in the same way as LagrangesMethod
        """

        if isinstance(Lagrangian, me.LagrangesMethod):
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
        self.Y = list(self.q) + list(self.u)

        if label == None:
            label = self.__class__.__name__ + ' with ' + str(len(
                self.q)) + 'DOF'

        self._label = label

        #LM=me.LagrangesMethod(Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs)
    def _kwargs(self):
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
        }

    def __add__(self, other):

        self_dict = self._kwargs()
        other_dict = other._kwargs()

        self_dict[
            'Lagrangian'] = self_dict['Lagrangian'] + other_dict['Lagrangian']

        self_dict['qs'] = list(self_dict['qs']) + list(
            coord
            for coord in other_dict['qs'] if not coord in self_dict['qs'])

        print(self_dict['qs'])

        list_build = lambda x: flatten([x]) if x else []

        self_dict['forcelist'] = list_build(
            self_dict['forcelist']) + list_build(other_dict['forcelist'])
        self_dict['bodies'] = list_build(self_dict['bodies']) + list_build(
            other_dict['bodies'])

        return LagrangesDynamicSystem(**self_dict)

    def shranked(self, *args):
        self_dict = self._kwargs()
        self_dict['qs'] = flatten(args)

        return LagrangesDynamicSystem(**self_dict)

    def remove(self, *args):

        bounded_coordinates = flatten(args)

        self_dict = self._kwargs()
        self_dict['qs'] = [
            coord for coord in self.q if coord not in bounded_coordinates
        ]

        return LagrangesDynamicSystem(**self_dict)

    def subs(self, *args, **kwargs):

        if 'method' in kwargs.keys():
            method = kwargs['method']
        else:
            method = 'as_constrains'

        if method == 'as_constrains':
            if len(args) == 2:
                args = ([(args[0], args[1])]),

            elif len(args) == 1:
                if isinstance(args[0], dict):
                    #print(args[0])
                    args = list(args[0].items()),

            constrains = [
                lhs - rhs for lhs, rhs in args[0] if any(
                    coord in (lhs - rhs).atoms(Function) for coord in self.q)
            ]
            args = ([(lhs, rhs) for lhs, rhs in args[0]
                     if not any(coord in (lhs - rhs).atoms(Function)
                                for coord in self.q)], )

#         print(constrains)
#         print(args)

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

        #         print(forces_subs)
        #         print(self.forcelist)

        return type(self)(Lagrangian=lagrangian_subs,
                          qs=self.q,
                          forcelist=forces_subs,
                          bodies=self._bodies,
                          frame=self.frame,
                          hol_coneqs=self._hol_coneqs,
                          nonhol_coneqs=list(
                              self.coneqs)[len((self._hol_coneqs)):],
                          ivar=self.ivar)


#     def rhs(self):
#         return self.rhs

    def __preview(self, expr, preview=None, preview_mode=None):
        display(expr)
        pass
        '''
        Private method which processes preview printing. The method is or isn't used depending on user's bool 'preview' input.
        '''

    def __call__(self, label=None):
        self._label = label

        return self

    def __str__(self):

        return self._label

    def __repr__(self):

        return self.__str__()

    def lagrangian(self):
        '''
        Returns the system Lagrange function defined within LagrangeMethod instance. Defined as self._L = Matrix([sympify(Lagrangian)]) in Sympy documentation.
        '''

        return sum(self.L)

    def rhs_eq(self):
        '''
        Returns the right-hand side of the equations of motion in the form of equations matrix. Output is provided basing on LagrangesMethod object.
        '''

        return Matrix([
            Eq(comp, self.rhs()[no], evaluate=False) for no, comp in enumerate(
                list(self.q.diff(self.ivar)) + list(self.q.diff(self.ivar, 2)))
        ])

    def equilibrium_equation(self, static_disp_dict=None):
        '''
        Finds the equilibrium conditions of the considered problem based on the system governing equations stated in the class instance.
        '''
        if static_disp_dict is None:
            static_disp_dict = {
                q_tmp:
                Symbol(str(q_tmp).replace('(' + str(self.ivar) + ')', '_0'))
                for q_tmp in self.q
            }

        self.q_0 = static_disp_dict

        return self.governing_equations.subs(static_disp_dict)

    def external_forces(self):
        return self.governing_equations.subs(
            {gen_coord: 0
             for gen_coord in self.Y}).doit()

    def _op_points(self,
                   hint=[],
                   static_disp_dict=None,
                   dict=True,
                   *args,
                   **kwargs):
        '''
        Provides the interface for critical points evaluation (solves the equlibrium conditions and returns its roots).
        '''

        eqns_to_solve = self.equilibrium_equation(
            static_disp_dict=static_disp_dict).doit()
        #         display(args,kwargs)
        roots = solve(list(eqns_to_solve) + list(flatten([hint])),
                      list(self.q_0.values()),
                      dict=dict)

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

    def approximated(self, n=3, x0=None, label=None):

        lagrangian_approx = multivariable_taylor_series(
            self.lagrangian(), self.Y, n=n, x0={coord: 0
                                                for coord in self.Y})

        return LagrangesDynamicSystem(lagrangian_approx,
                                      self.q,
                                      forcelist=self.forcelist,
                                      frame=self.frame,
                                      label=label,
                                      ivar=self.ivar)

    def linearized(self, x0=None, label=None):

        linearized_sys = self.approximated(n=2, x0=x0)

        return LinearDynamicSystem(linearized_sys.lagrangian(),
                                   self.q,
                                   forcelist=self.forcelist,
                                   frame=self.frame,
                                   label=label,
                                   ivar=self.ivar)

    @property
    def _eoms(self):
        return self.governing_equations

    def numerized(self, parameter_values=None):
        '''
        Takes values of parameters, substitute it into the list of parameters and changes list it into a Tuple. Returns instance of class OdeComputationalCase.
        '''
        data_Tuple = Tuple(*self.system_parameters()).subs(parameter_values)
        computed_case = self.computational_case(parameter_values=data_Tuple)

        return OdeComputationalCase(**computed_case, evaluate=True)


class LinearDynamicSystem(LagrangesDynamicSystem):
    '''

    '''
    def stiffness_matrix(self):
        '''
        Returns the system stiffness matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized coordinates vector.
        '''
        return self.governing_equations.jacobian(self.q)

    def fundamental_matrix(self, freq=Symbol('omega', positive=True)):
        return -freq**2 * self.inertia_matrix() + self.stiffness_matrix()
        '''
        Method returns a fundamental matrix of the system built from inertia and stiffness matrices. Takes one optional argument.
        
        Args:
            freq (optional, obj:Symbol): This argument lets user to choose a symbol for frequency representation. Set to 'omega' by default, 'positive=True' ensures its affiliation to real numbers domain.
        '''

    def damping_matrix(self):
        '''
        Returns the system damping matrix, which is based on the equations of motion of the Lagrange's system. Matrix is obtained from jacobian which is called with system's generalized velocities vector.
        '''
        return self.governing_equations.jacobian(self.u)

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        return ((self.inertia_matrix().inv() *
                 self.stiffness_matrix()).diagonalize()[1])

    def modes(self):
        '''
        Determines the system vibration modes matrix (eigenmodes) based on the system free vibration frequencies.
        '''
        return (((self.inertia_matrix().inv() *
                  self.stiffness_matrix()).diagonalize()[0]))


class HarmonicOscillator(LinearDynamicSystem):
    '''
    This object allows for a determination of any dynamic system by providing methods that serve to generate the equations of motion, solution, 
    natural frequencies, eigenmodes, FRF of the mechanical system and many others that are discussed further within the documentation. 
    For the initialization of HarmonicOscillator at least one argument is neccessary.
    '''
    '''

    Arguments
    ==========
    LagrangesSystem : obj 
        Formerly prepared LagrangesMethod object as explained in Sympy's 'Lagrange’s Method in Physics/Mechanics' module.

    preview (optional) : bool
        False by default. Allows user to display the provided formula, e.g. for checking its correctness.

    preview_mode (optional) : str
        Functional with 'preview' set as True. Lets user to choose the method of rendering expressions. Set to LaTeX by default.

    ivar (optional) : obj:Symbol
        This argument lets user to change the default 't' symbol for time derivative to other, for compatibility with LagrangesSystem expressions.

    For first step, LagrangesMethod object has to be formulated. For its proper definition check the Sympy LagrangesMethod documentation.

    Example

    ==========

    >>>import sympy as sym
    >>>from sympy import *
    >>>from sympy.physics.vector import dynamicsymbols
    >>>from sympy.physics.mechanics import *
    >>>t, m1, m2, k1, k2, c1, c2, F= symbols('t m1 m2 k1 k2 c1 c2 F(s)')
    >>>x1, x2 =dynamicsymbols('x1 x2')
    >>>dx1,dx2=dynamicsymbols('x1 x2', 1)
    >>>N1=ReferenceFrame('N1')
    >>>P1=Point('P1')
    >>>P2=Point('P2')
    >>>P1.set_vel(N1, dx1 * N1.x)
    >>>P2.set_vel(N1, dx2 * N1.x)
    >>>R1=1/2*c2*(dx1-dx2)**2
    >>>R2=1/2*c1*dx2**2
    >>>R=R1+R2
    >>>dR1=R.diff(dx1)
    >>>dR2=R.diff(dx2)
    >>>Pa1 = Particle('Pa1', P1, m1)
    >>>Pa2 = Particle('Pa2', P2, m2)
    >>>Pa2.potential_energy=1/2*k1*x2**2+1/2*k2*(x1-x2)**2
    >>>Pa1.potential_energy=1/2*k2*(x1-x2)**2
    >>>V=Pa2.potential_energy+Pa1.potential_energy
    >>>L1 = Lagrangian(N1,Pa1)
    >>>L2 = Lagrangian(N1,Pa2)
    >>>L=L1+L2
    >>>FL = [(P2, -dR2*N1.x),(P1, -dR1*N1.x)]
    >>>lall=LagrangesMethod(L,[x1,x2],forcelist=FL,frame=N1)

    At this point we will work only using the defined class. Below an instance called instance_one is created. From now on the instance name can and should be used to call any desired method.

    >>>instance_one=HarmonicOscillator(lall)

    '''
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

    def natural_frequency(self):
        '''
        Determines the expression representing the system natural frequency (working correctly for single degree of freedom systems).
        '''
        return sqrt(sum(self.stiffness_matrix()) / sum(self.inertia_matrix()))

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

        self._canonical_governing_equation = Matrix(u).diff(
            self.ivar) + 2 * self.h * Matrix(u) + (self.omegah**2 +
                                                   self.h**2) * Matrix(q)
        print('Done')
        return self.__governing_equation

    def __solve(self, initial_conditions=None):
        '''
        Private method solving the problem in the symbolic way and sets the related attribute.
        '''
        if len(self.q) == 1:
            solution = dsolve(sum(self.__governing_equations),
                              sum(q),
                              ics=initial_conditions)
        else:
            solution = dsolve((self.__governing_equations), (q),
                              ics=initial_conditions)

        return solution

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''
        C = numbered_symbols('C', start=1)

        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        Y_mat = Matrix(self.q)

        diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        t_sol = self.ivar

        solution = [
            next(C) * modes[:, i] * sin(sym.sqrt(eigs[i, i]) * t_sol) +
            next(C) * modes[:, i] * cos(sym.sqrt(eigs[i, i]) * t_sol)
            for i, coord in enumerate(self.q)
        ]

        return sum(solution, Matrix([0] * len(Y_mat)))

    def steady_solution(self, initial_conditions=None):
        """
        
        """

        ext_forces = self.external_forces()

        #         sin_components=ext_forces.atoms(sin)
        #         cos_components=ext_forces.atoms(cos)
        components = ext_forces.atoms(sin, cos)

        display(components)

        steady_sol = Matrix([0 for gen_coord in self.q])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
            display(omg)
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

            #display(amp_vector)

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (fund_mat.inv() * amp_vector) * comp

            ext_forces -= amp_vector * comp
        #print(ext_forces.doit().expand())
        if ext_forces.doit().expand() != sym.Matrix(
            [0 for gen_coord in self.q]):
            steady_sol += sym.dsolve(
                (self.governing_equations - self.external_forces() +
                 ext_forces).expand().doit(), self.q)

        return steady_sol

    def solution(self, initial_conditions=None):
        return self.general_solution(initial_conditions=initial_conditions)


#     def solution(self, initial_conditions=None):
#         '''
#         Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
#         '''

#         if initial_conditions == None:
#             if self.__general_solution == None:
#                 self.__general_solution = self.__solve()

#             solution = self.__general_solution

#         else:
#             if self.__particular_solution == None:
#                 self.__particular_solution = self.__solve(
#                     initial_conditions=initial_conditions)

#             solution = self.__particular_solution

#         return solution

#     def integration_constants(self,natural_freq=None,initial_conditions=None):
#         if natural_freq==None:
#             self.omega0=Symbol('omega_0',positive=True)

#         solution=self.solution(initial_conditions=initial_conditions).rhs
#         =initial_conditions
#         C1_comp=solution.coeff(exp**(imag*omega0*t))
#         C2_comp=solution.coeff(exp**(-imag*omega0*t))

    def steady_solution_amp(
            self,
            cos_amp=None,
            sin_amp=None,
            excitation_freq=Symbol('Omega'),
    ):
        ''''
        Computes the steady solution amplitude for the system defined in the instance of this class.
        '''
        if excitation_freq == None:
            excitation_freq = Symbol('Omega', positive=True)

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
        Returns the Frequency Response Function of the system for the given excitation amplitude (working correctly for single degree of freedom systems).
        '''
        if excitation_freq == None:
            excitation_freq = Symbol('Omega', positive=True)

        self.Omega = excitation_freq

        general_solution = self.solution().rhs

        n, d = fraction(general_solution)
        comp_sin = n.expand().coeff(sin(excitation_freq * self.ivar)) / d
        comp_cos = n.expand().coeff(cos(excitation_freq * self.ivar)) / d

        frf_expr = (sqrt((comp_sin**2 + comp_cos**2).simplify()))

        return frf_expr

    def dynamic_amplification_factor(self,
                                     excitation_amp=None,
                                     excitation_freq=None):
        '''
        Returns the Dynamic Amplification Factor of the system for the given excitation amplitude (working correctly for single degree of freedom systems).
        '''
        frf = self.frequency_response_function(excitation_freq=excitation_freq)
        eom = sum(self.governing_equations)

        sin_comp = eom.expand().coeff(sin(excitation_freq * self.ivar))
        cos_comp = eom.expand().coeff(cos(excitation_freq * self.ivar))

        maximal_force = sqrt(sin_comp**2 + cos_comp**2).subs(
            excitation_freq, 1).simplify()

        static_def = maximal_force / sum(self.stiffness_matrix())

        _dummy = symbols('_dummy', positive=True)

        nat_freq = self.natural_frequency()

        daf = (frf.subs(excitation_freq, _dummy * nat_freq) /
               static_def).simplify().subs(_dummy, excitation_freq / nat_freq)

        return daf

    def critical_frequencies(self, excitation_freq=None):
        '''
        Determines the critical frequency of the system (working correctly for single degree of freedom systems).
        '''
        frf_max_cond = (self.frequency_response_function(
            excitation_freq=excitation_freq)**2).diff(excitation_freq)
        #display(frf_max_cond)

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

        general_solution = self.solution().rhs

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
    pass
