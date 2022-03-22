from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq,
                   hessian, Function, flatten, Tuple, im, re, pi, latex,
                   dsolve, solve, fraction, factorial, Add, Mul, exp,
                   numbered_symbols, integrate, ImmutableMatrix)

from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from ..utilities.timeseries import TimeSeries, TimeDataFrame

from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3


class FirstOrderODE:
    """
    This class represents the system of first order defferential equations. Allows to get the solution of given system.
    """

    _const_list = set()

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
        Supply the following arguments for the initialization of FirstOrderODE
        
        Args:
        '''
        #         if isinstance(odes_system, LinearDynamicSystem):
        #             self.dvars = odes_system.q
        #             odes_system = odes_system._eoms

        self.odes_system = odes_system

        #         print('solver init')

        #         display(self.odes_system)

        self.governing_equations = self.odes_system
        self.ivar = ivar

        if len(dvars) > 0:
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

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        return odes.jacobian(reg_vars).subs({coord: 0
                                             for coord in reg_vars}).doit()

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        dvars_ddot = list(sym.Matrix(reg_vars).diff(self.ivar, 2))

        result = odes.jacobian(dvars_ddot).subs(
            {coord: 0
             for coord in reg_vars}).doit()

        #         print('------------------- linear mat ------------------')
        #         display(self.odes_system,self.dvars,dvars_ddot,result)

        #         print('------------------- linear mat ------------------')

        return result

    def damping_matrix(self):
        '''
        Returns the system damping matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        dvars_dot = list(sym.Matrix(reg_vars).diff(self.ivar, 1))

        result = odes.jacobian(dvars_dot).subs(
            {coord: 0
             for coord in reg_vars}).doit()

        return result

    def external_forces(self):
        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        return odes.subs({gen_coord: 0 for gen_coord in self.dvars}).doit()

    def diagonalize(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        ode_sys = self.odes_system

        #display(ode_sys)

        main_matrix = ode_sys.jacobian(self.dvars).subs(
            {dvar: 0
             for dvar in self.dvars})

        inertia_matrix = ode_sys.jacobian(self.dvars).subs(
            {dvar.diff(self.ivar): 0
             for dvar in self.dvars}) #.diff(self.ivar)

        #         display(inertia_matrix,main_matrix)

        linear_odes = inertia_matrix.inv() * (
            main_matrix * sym.Matrix(self.dvars) +
            ode_sys.subs({dvar: 0
                          for dvar in self.dvars}).doit())

        linear_odes_dict = {
            coord: eq
            for coord, eq in zip(
                self.dvars,
                linear_odes,
            )
        }

        const_odes_list = {}
        regular_odes_list = {}

        for coord, ode in (linear_odes_dict.items()):
            if ode == 0:
                const_odes_list[coord] = ode
            else:
                regular_odes_list[coord] = ode


#         print('regular odes')
#         display(*list(regular_odes_list.values()))

        if isinstance(self.ivar, Function):
            regular_vars = Matrix(list(
                regular_odes_list.values())).atoms(Function) - {self.ivar}
        else:
            regular_vars = Matrix(list(
                regular_odes_list.values())).atoms(Function)

        regular_vars = regular_odes_list.keys()

        #         display(regular_vars)

        #display( row  if row  in Matrix(regular_odes_list).rows  )
        #         display( list(regular_odes_list.values()) )
        regular_main_matrix = Matrix(list(
            regular_odes_list.values())).jacobian(list(regular_vars))
        #singular_odes=[no  for no in  const_odes_list]
        #         print('diagonalize')
        #         print([Matrix(list(regular_odes_list.values())) ,list(regular_vars) ,const_odes_list])
        #         print('diagonalize')

        return [
            Matrix(list(regular_odes_list.values())),
            list(regular_vars), const_odes_list
        ]

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        odes, vars = self.diagonalize()[0:2]

        return odes.jacobian(vars).diagonalize()[1]

    def eigenmodes(self):
        '''
        Determines the system eigenmodes matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        odes, vars = self.diagonalize()[0:2]
        #         display(odes,vars)

        return odes.jacobian(vars).diagonalize()[0]

    def damped_natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        natural_freqs = list({(im((eigen))).doit()
                              for eigen in self.eigenvalues()
                              if not eigen == 0})

        return diag(*natural_freqs)

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''

        #         print('------------------- linear gen sol ----------------')
        #         display(self.odes_system,self.dvars)

        #         display(self.inertia_matrix(),self.stiffness_matrix())
        #         print('------------------- linear gen sol ----------------')

        C = numbered_symbols('C', start=1)

        C_list = []

        for i in range(len(self.dvars) * 2):

            C_list += [next(C)]

        args_list = self.dvars[0].args

        if len(self.dvars[0].args) > 1:

            params = {*args_list} - {self.ivar}

            C_list = [Function(str(C_tmp)) for C_tmp in C_list]
            C_list = [(C_tmp)(*params) for C_tmp in C_list]

        #         print('o tu')
        self.__class__._const_list |= set(C_list)

        modes, eigs = self.eigenmodes(), self.eigenvalues()
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        #         display( (const_part))

        #Y_mat = Matrix(self.dvars)

        #diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        # t_sol = self.ivar

#         display(modes,eigs)
        solution = [
            C_list[i] * modes[:, i] * exp(eigv * self.ivar)
            for i, eigv in enumerate([eigv for eigv in eigs if not eigv == 0])
        ]
        #         print('sol')
        #         display(solution)
        #         print('sol')

        const_dict = ({
            coord: C_list[no + len(solution)]
            for no, coord in enumerate(const_part.keys())
        })

        sol_dict = {
            coord: eq
            for coord, eq in zip(
                reg_vars, list(sum(solution, Matrix([0] * len(solution[0])))))
        }

        return {**sol_dict, **const_dict}

    def steady_solution(self, initial_conditions=None):

        odes = self.diagonalize()[0]
        const_part = self.diagonalize()[2]
        reg_vars = self.diagonalize()[1]

        ext_forces = self.external_forces().expand().applyfunc(
            lambda row: (TR8(row).expand()))

        #         sin_components=ext_forces.atoms(sin)
        #         cos_components=ext_forces.atoms(cos)

        #        display('sin cos reco',)
        components = ext_forces.atoms(sin, cos)

#         display('ext_forces', ext_forces)
#         display(odes)
#         display(reg_vars)

        steady_sol = Matrix([0 for gen_coord in self.dvars])

        for comp in components:

            omg = (comp.args[0].diff(self.ivar)).doit()
#             print('omg')
#             display(omg)
            amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

#             display(self.inertia_matrix())
#             display(self.damping_matrix())

            fund_mat = -self.inertia_matrix(
            ) * omg**2 + sym.I * omg * self.damping_matrix(
            ) + self.stiffness_matrix()

            steady_sol += (
                sym.re(fund_mat.inv()).doit() * amp_vector) * comp + (
                    sym.im(fund_mat.inv()).doit() * amp_vector) * TR3(
                        comp.subs(omg * self.ivar, omg * self.ivar - pi / 2))

            #steady_sol += ((fund_mat.inv()) * amp_vector) * comp

            ext_forces -= (amp_vector * comp).expand()
        #print(ext_forces.doit().expand())

        # const_elems = lambda expr: [
        #     expr for expr in comp.expand().args if not expr.has(self.ivar)
        #     if isinstance(comp, Add)
        # ]

        const_mat = Matrix([
            sum((expr for expr in comp.expand().args if not expr.has(self.ivar)
                 if isinstance(comp, Add)), 0) for comp in ext_forces.doit()
        ])

        #        display('const',const_mat)

        steady_sol += (self.stiffness_matrix().inv() * const_mat).doit()

        ext_forces -= const_mat

        #        display('res',ext_forces)

        if ext_forces.doit().expand() != sym.Matrix(
            [0 for gen_coord in self.dvars]):
            #             print('o tu')
            #             display((self.governing_equations - self.external_forces() +
            #                      ext_forces).expand().doit())
            eqns_res = (self.governing_equations - self.external_forces() +
                        ext_forces).expand().doit()

            if len(self.dvars) == 1:
                steady_sol += Matrix([
                    sym.dsolve(eqns_res[0], self.dvars[0]).rhs.subs({
                        Symbol('C1'):
                        0,
                        Symbol('C2'):
                        0
                    })
                ])
            else:
                steady_sol += sym.dsolve(eqns_res, self.dvars)

        return {coord: -eq for coord, eq in zip(reg_vars, steady_sol)}

    def solution(self, initial_conditions=None):
        return {
            key:
            self.general_solution(initial_conditions=initial_conditions)[key] +
            self.steady_solution(initial_conditions=initial_conditions)[key]
            for key in self.general_solution(
                initial_conditions=initial_conditions).keys()
        }


class LinearODESolution:

    _const_list = set()
    _cache={}
    

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
        #         if isinstance(odes_system, LinearDynamicSystem):
        #             self.dvars = odes_system.q
        #             odes_system = odes_system._eoms

        self.odes_system = odes_system

        #         print('solver init')

        #         display(self.odes_system)

        self.governing_equations = self.odes_system
        self.ivar = ivar

        if len(dvars) > 0:
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
        return self.governing_equations.jacobian(self.dvars).subs(
            {coord: 0
             for coord in self.dvars}).doit()

    def inertia_matrix(self):
        '''
        Returns the system inertia matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_ddot = list(sym.Matrix(self.dvars).diff(self.ivar, 2))

        result = self.governing_equations.jacobian(dvars_ddot).subs(
            {coord: 0
             for coord in self.dvars}).doit()

        #         print('------------------- linear mat ------------------')
        #         display(self.odes_system,self.dvars,dvars_ddot,result)

        #         print('------------------- linear mat ------------------')

        return result
    
    def fundamental_matrix(self, freq=Symbol('omega', positive=True)):
        '''
        Method returns a fundamental matrix of the system built from inertia and stiffness matrices. Takes one optional argument.

        Args:
            freq (optional, obj:Symbol): This argument lets user to choose a symbol for frequency representation. Set to 'omega' by default, 'positive=True' ensures its affiliation to real numbers domain.
        '''
        return -freq**2 * self.inertia_matrix() + self.stiffness_matrix()

    def damping_matrix(self):
        '''
        Returns the system damping matrix which is based on the equations of motion of the Lagrange's system. mass_matrix is an argument of sympy.physics.mechanics.lagrange module.
        '''
        dvars_dot = list(sym.Matrix(self.dvars).diff(self.ivar, 1))

        return self.governing_equations.jacobian(dvars_dot).subs(
            {coord: 0
             for coord in self.dvars}).doit()

    def external_forces(self):
        return self.odes_system.subs(
            {gen_coord: 0
             for gen_coord in self.dvars}).doit()

    def eigenvalues(self):
        '''
        Determines the system eigenvalues matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        #         display(self.governing_equations)

        q_dot = (Matrix(self.dvars).diff(self.ivar))

        ode_sys = Matrix([
            q_dot,
            self.inertia_matrix().inv() *
            (-self.stiffness_matrix() * Matrix(self.dvars) -
             self.damping_matrix() * q_dot)
        ])

        #         display(ode_sys)

        main_matrix = ode_sys.jacobian(list(self.dvars) + list(q_dot))

        #         display(main_matrix)

        return (main_matrix).diagonalize()[1]

    def damped_natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''

        natural_freqs = list({(im((eigen))).doit()
                              for eigen in self.eigenvalues()
                              if not eigen == 0})

        return diag(*natural_freqs)
    
    def natural_frequencies(self):
        '''
        Determines the system natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        
        natural_freqs=list({(sqrt(abs(eigen**2))) for eigen  in self.eigenvalues() if not eigen==0 })
        
        return diag(*natural_freqs)

    def general_solution(self, initial_conditions=None):
        '''
        Solves the problem in the symbolic way and rteurns matrix of solution (in the form of equations (objects of Eq class)).
        '''

        #         print('------------------- linear gen sol ----------------')
        #         display(self.odes_system,self.dvars)

        #         display(self.inertia_matrix(),self.stiffness_matrix())
        #         print('------------------- linear gen sol ----------------')

        C = numbered_symbols('C', start=1)

        C_list = []

        for i in range(len(self.dvars) * 2):

            C_list += [next(C)]

        args_list = self.dvars[0].args

        if len(self.dvars[0].args) > 1:

            params = {*args_list} - {self.ivar}

            C_list = [Function(str(C_tmp)) for C_tmp in C_list]
            C_list = [(C_tmp)(*params) for C_tmp in C_list]

        #         print('o tu')
        #         display(self.odes_system)

        self.__class__._const_list |= set(C_list)

        modes, eigs = ((self.inertia_matrix().inv() *
                        self.stiffness_matrix()).diagonalize())

        eigs = self.eigenvalues()

        Y_mat = Matrix(self.dvars)

        #diff_eqs = Y_mat.diff(self.ivar, 2) + eigs * Y_mat

        t_sol = self.ivar

        solution = [(C_list[2 * i] * modes[:, i] *
                     sin(im(eigs[2 * i + 1, 2 * i + 1]).doit() * t_sol) +
                     C_list[2 * i + 1] * modes[:, i] *
                     cos(im(eigs[2 * i + 1, 2 * i + 1]).doit() * t_sol)) *
                    exp(re(eigs[i, i]).doit() * t_sol)
                    for i, coord in enumerate(self.dvars)]

        return sum(solution, Matrix([0] * len(Y_mat)))

    def steady_solution(self, initial_conditions=None):
        if (tuple(self.dvars),ImmutableMatrix(self.governing_equations)) in self.__class__._cache:


            print('cache for Linear ODE')
            steady_sol=self.__class__._cache[(tuple(self.dvars),ImmutableMatrix(self.governing_equations))]
            
        else:
            print('new solution')
            ext_forces = self.external_forces().expand().applyfunc(
                lambda row: (TR8(row).expand()))

            #         sin_components=ext_forces.atoms(sin)
            #         cos_components=ext_forces.atoms(cos)

            #        display('sin cos reco',)
            components = [
                comp for comp in ext_forces.atoms(sin, cos) if comp.has(self.ivar)
            ]

            #        display('ext_forces',ext_forces)

            steady_sol = Matrix([0 for gen_coord in self.dvars])

            for comp in components:
                #display(comp)

                omg = (comp.args[0].diff(self.ivar)).doit()
                #            display(omg)
                amp_vector = Matrix([row.coeff(comp) for row in ext_forces])

                #display(amp_vector)

                fund_mat = -self.inertia_matrix(
                ) * omg**2 + sym.I * omg * self.damping_matrix(
                ) + self.stiffness_matrix()

                #print('++linea.py++'*10)
                #display(fund_mat)
                #print('++linea.py++'*10)

                steady_sol += (
                    sym.re(fund_mat.inv()).doit() * amp_vector) * comp + (
                        sym.im(fund_mat.inv()).doit() * amp_vector) * TR3(
                            comp.subs(omg * self.ivar, omg * self.ivar - pi / 2))

                #steady_sol += ((fund_mat.inv()) * amp_vector) * comp

                ext_forces -= (amp_vector * comp).expand()
            #print(ext_forces.doit().expand())

            # const_elems = lambda expr: [
            #     expr for expr in comp.expand().args if not expr.has(self.ivar)
            #     if isinstance(comp, Add)
            # ]

            const_mat = Matrix([
                sum((expr for expr in comp.expand().args if not expr.has(self.ivar)
                     if isinstance(comp, Add)), 0) for comp in ext_forces.doit()
            ])

            #        display('const',const_mat)

            steady_sol += (self.stiffness_matrix().inv() * const_mat).doit()

            ext_forces -= const_mat

            #        display('res',ext_forces)

            if ext_forces.doit().expand() != sym.Matrix(
                [0 for gen_coord in self.dvars]):
                #             print('o tu')
                #             display((self.governing_equations - self.external_forces() +
                #                      ext_forces).expand().doit())
                eqns_res = (self.governing_equations - self.external_forces() +
                            ext_forces).expand().doit()

                if len(self.dvars) == 1:
                    steady_sol += Matrix([
                        sym.dsolve(eqns_res[0], self.dvars[0]).rhs.subs({
                            Symbol('C1'):
                            0,
                            Symbol('C2'):
                            0
                        })
                    ])
                else:
                    steady_sol += sym.dsolve(eqns_res, self.dvars)
                    
            self.__class__._cache[(tuple(self.dvars),ImmutableMatrix(self.governing_equations))] = steady_sol

                
        return -steady_sol

    def solution(self, initial_conditions=None):
        return self.general_solution(
            initial_conditions=initial_conditions) + self.steady_solution(
                initial_conditions=initial_conditions)
