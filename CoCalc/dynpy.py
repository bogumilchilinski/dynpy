from sympy import *
from sympy.physics.mechanics import *
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
from IPython.display import display
import sympy.physics.mechanics as me


def multivariable_taylor_series(expr, args, n=2, x0=None):
    order_max = n
    op_point = x0

    args_shifted = {
        arg: arg - arg_shift
        for arg, arg_shift in op_point.items()
    }
    diff_orders = lambda var_list, order_arg: [
        args for tmp in range(order_arg)
    ]

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

    return (sum([
        expr.diff(*args_tmp).subs(op_point).expand().doit() * poly
        for args_tmp, poly in diff_orders_dict.items()]) + expr.subs(op_point)).doit()

def scalar_fun_quadratic_form(expr, coordinates, op_point):
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
        self._label = label
        return self

    def __str__(self):
        return self._label

    def __repr__(self):
        return self.__str__()

    def __fortran_odes_rhs(self):
        subs_dict = {
            var: Symbol('temp_sym_' + str(i))
            for i, var in enumerate(self.dvars)
        }

        args_list = [self.ivar] + list(subs_dict.values()) + self.params

        return autowrap((self.odes_system.subs(subs_dict, simultaneous=True)),
                        args=args_list)

    def form_numerical_rhs(self):
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

# -------------------------------------------------------------------------- Linear ODE Solution
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
        return self.governing_equations.jacobian(self.dvars).subs({coord:0 for coord  in self.dvars}).doit()

    def inertia_matrix(self):
        dvars_ddot = list(sym.Matrix(self.dvars).diff(self.ivar, 2))
        return self.governing_equations.jacobian(dvars_ddot).subs({coord:0 for coord  in self.dvars}).doit()

    def damping_matrix(self):
        dvars_dot = list(sym.Matrix(self.dvars).diff(self.ivar, 1))
        return self.governing_equations.jacobian(dvars_dot).subs({coord:0 for coord  in self.dvars}).doit()

    def external_forces(self):
        return self.odes_system.subs(
            {gen_coord: 0
             for gen_coord in self.dvars}).doit()

    def general_solution(self, initial_conditions=None):
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

        ext_forces = self.external_forces().expand().applyfunc(lambda row: (TR8(row).expand())  )
        components = ext_forces.atoms(sin, cos)
        steady_sol = Matrix([0 for gen_coord in self.dvars])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (fund_mat.inv() * amp_vector) * comp

            ext_forces -= (amp_vector * comp).expand()
        
        const_mat=Matrix([sum((expr for expr in comp.expand().args if not expr.has(self.ivar)),0)  for  comp in ext_forces.doit()])
                
        steady_sol += (self.stiffness_matrix().inv() * const_mat).doit()
        ext_forces -=const_mat
        
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




# ---------------------------------------------------------------------- Lagranes Dynamic System
class LagrangesDynamicSystem(me.LagrangesMethod):

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
                    args = list(args[0].items()),

            constrains = [
                lhs - rhs for lhs, rhs in args[0] if any(
                    coord in (lhs - rhs).atoms(Function) for coord in self.q)
            ]
            args = ([(lhs, rhs) for lhs, rhs in args[0]
                     if not any(coord in (lhs - rhs).atoms(Function)
                                for coord in self.q)], )

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

        return type(self)(Lagrangian=lagrangian_subs,
                          qs=self.q,
                          forcelist=forces_subs,
                          bodies=self._bodies,
                          frame=self.frame,
                          hol_coneqs=self._hol_coneqs,
                          nonhol_coneqs=list(
                              self.coneqs)[len((self._hol_coneqs)):],
                          ivar=self.ivar)

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

        roots = solve(list(eqns_to_solve) + list(flatten([hint])),
                      list(self.q_0.values()),
                      dict=dict)

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
            self.lagrangian(), self.Y, n=n+1, x0={coord: 0
                                                for coord in self.Y})

        return LagrangesDynamicSystem(lagrangian_approx,
                                      self.q,
                                      forcelist=self.forcelist,
                                      frame=self.frame,
                                      label=label,
                                      ivar=self.ivar)

    def linearized(self, x0=None, label=None):

        linearized_sys = self.approximated(n=1, x0=x0)

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