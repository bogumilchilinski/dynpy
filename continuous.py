from sympy import (flatten, SeqFormula, Function, Symbol, symbols, Eq, Matrix,
                   S, oo, dsolve, solve, Number, pi, cos, sin, Tuple, Derivative)

from sympy.physics.vector.printing import vpprint, vlatex


class ContinuousSystem:
    def __init__(self,
                 L,
                 q,
                 bc_dict=None,
                 t_var=Symbol('t'),
                 spatial_var=Symbol('x'),
                 derivative_order=2,
                 label=None,
                 system=None,
                 **kwargs):

        if system:
            t_var = system.t
            L = system.L
            q = system.q
            spatial_var = system.r
            derivative_order = system.diff_ord
            bc_dict = system.bc_dict

        self.t = t_var
        self.L = L
        self.q = q
        self.r = flatten((spatial_var, ))
        self.diff_ord = derivative_order
        self.bc_dict = bc_dict
        self._bc = Symbol('BC')

        self._sep_expr = Symbol('k', positive=True)

        if label == None:
            label = self.__class__.__name__ + ' on ' + str(self.q)

        self._label = label

    @property
    def BC(self):
        return self._bc

    def __call__(self, *args, label=None):
        """
        Returns the label of the object or class instance with reduced Degrees of Freedom.
        """

        if isinstance(args[0], str):
            if label:
                self._label = label
            else:
                self._label = args[0]

            return self

    def __str__(self):

        return self._label

    def __repr__(self):

        return self.__str__()

    def subs(self, *args, **kwargs):
        L_new = self.L.subs(*args)
        q_new = self.q.subs(*args)

        bc_trap = self.BC.subs(*args)

        if bc_trap == self.BC:
            bc_trap = self.bc_dict

        new_system = ContinuousSystem(L=L_new,
                                      q=q_new,
                                      bc_dict=bc_trap,
                                      t_var=self.t,
                                      spatial_var=self.r,
                                      derivative_order=self.diff_ord,
                                      label=self._label)

        return type(self)(0, q_new, system=new_system)

    @property
    def _eoms(self):
        return self.governing_equation

    def inertia_force(self):
        q = self.q
        t = self.t
        L = self.L

        return L.diff(q.diff(t)).diff(t)

    def restoring_force(self):
        q = self.q
        t = self.t
        L = self.L

        return sum([(-1)**(order + 2) *
                    L.diff(q.diff(r_var, order + 1)).diff(r_var, order + 1)
                    for r_var in self.r
                    for order in range(self.diff_ord)]) - L.diff(q)

    def governing_equation(self):

        return self.inertia_force() + self.restoring_force()

    def eom_coeff(self, expr):

        return self.governing_equation().coeff(expr)

    def apply_separation(self,
                         time_comp=Function('T')(Symbol('t')),
                         spatial_comp=Function('X')(Symbol('x'))):

        return self.governing_equation().subs(self.q,
                                              time_comp * spatial_comp).doit()

    def separated_vars_eqn(self,
                           time_comp=Function('T')(Symbol('t')),
                           spatial_comp=Function('X')(Symbol('x'))):

        eqn_with_subs = self.apply_separation(time_comp, spatial_comp)

        return Eq((eqn_with_subs.coeff(spatial_comp) / (time_comp)),
                  -eqn_with_subs.coeff(time_comp) / spatial_comp)

    def spatial_eqn(self,
                    sep_expr=None,
                    spatial_comp=Function('X')(Symbol('x'))):

        if not sep_expr:
            sep_expr = self._sep_expr

        separated_eqn_rhs = self.separated_vars_eqn(
            spatial_comp=spatial_comp).rhs

        return Eq(separated_eqn_rhs, sep_expr)

    def time_eqn(self, sep_expr=None, time_comp=Function('T')(Symbol('t'))):

        if not sep_expr:
            sep_expr = self._sep_expr

        separated_eqn_lhs = self.separated_vars_eqn(time_comp=time_comp).lhs

        return Eq(separated_eqn_lhs, sep_expr)

    def spatial_general_solution(self,
                                 sep_expr=None,
                                 spatial_comp=Function('X')(Symbol('x'))):

        if not sep_expr:
            sep_expr = self._sep_expr

        spatial_ode = self.spatial_eqn(sep_expr, spatial_comp)

        return dsolve(spatial_ode,
                      spatial_comp)  #.rewrite(cos).expand().simplify()

    def fundamental_matrix(self,
                           bc_dict=None,
                           sep_expr=None,
                           spatial_comp=Function('X')(Symbol('x'))):

        if not sep_expr:
            sep_expr = self._sep_expr

        if bc_dict:
            self.bc_dict = bc_dict
        else:
            bc_dict = self.bc_dict

        spatial_sol = self.spatial_general_solution(sep_expr=sep_expr,
                                                    spatial_comp=spatial_comp)

        fun_eqns = Matrix([
            spatial_sol.rhs.subs(spatial_comp.args[0],
                                 flatten([key.args[-1]])[0]) - val
            for key, val in bc_dict.items()
        ])

        matrix_comps_list = []

        for key, val in bc_dict.items():

            #display(list(key.atoms(Symbol,Number) - spatial_comp.atoms(Symbol)))
            free_sym = list(
                key.atoms(Symbol, Number) - spatial_comp.atoms(Symbol))[-1]
            #print({spatial_sol.lhs:spatial_sol.rhs, spatial_sol.lhs.subs(spatial_comp.args[0],free_sym):spatial_sol.rhs, })
            matrix_comps_list += [(key.subs({
                spatial_sol.lhs:
                spatial_sol.rhs,
                spatial_sol.lhs.subs(spatial_comp.args[0], free_sym):
                spatial_sol.rhs.subs(spatial_comp.args[0], free_sym),
            }).doit() - val)]

        fun_eqns = Matrix(matrix_comps_list)

        return fun_eqns.jacobian(symbols('C1:' + str(len(bc_dict) + 1)))

    def char_poly(self,
                  bc_dict=None,
                  sep_expr=None,
                  spatial_comp=Function('X')(Symbol('x'))):

        if not sep_expr:
            sep_expr = self._sep_expr

        if bc_dict:
            self.bc_dict = bc_dict
        else:
            bc_dict = self.bc_dict

        return self.fundamental_matrix(bc_dict, sep_expr,
                                       spatial_comp).det().simplify()


#     def eigenvalues(self, bc_dict=None, sep_expr=Symbol('k',positive=True), arg=Symbol('k',positive=True), spatial_comp=Function('X')(Symbol('x')), index=Symbol('n', integer=True, positive=True)):

#         if bc_dict:
#             self.bc_dict=bc_dict
#         else:
#             bc_dict=self.bc_dict

#         root = solve(self.char_poly(bc_dict, sep_expr, spatial_comp), arg)[0]

#         spatial_span = list(root.free_symbols)[0]

#         return SeqFormula(root+(index-1)/spatial_span*pi, (index, 0, oo))

    def eigenvalues(self,
                    bc_dict=None,
                    sep_expr=None,
                    arg=Symbol('k', positive=True),
                    spatial_comp=Function('X')(Symbol('x')),
                    index=Symbol('n', integer=True, positive=True)):

        if not sep_expr:
            sep_expr = self._sep_expr

        if bc_dict:
            self.bc_dict = bc_dict
        else:
            bc_dict = self.bc_dict

        roots = solve(self.char_poly(bc_dict, sep_expr, spatial_comp), arg)
        print(roots)

        if len(roots) == 1:
            spatial_span = roots[0]
        else:
            spatial_span = roots[1] - roots[0]

        return SeqFormula(roots[0] + (index - 1) * spatial_span,
                          (index, 0, oo))

    def eigenmodes(self,
                   mode_no,
                   bc_dict=None,
                   sep_expr=None,
                   arg=Symbol('k', positive=True),
                   spatial_comp=Function('X')(Symbol('x')),
                   index=Symbol('n', integer=True, positive=True)):

        if not sep_expr:
            sep_expr = self._sep_expr

        if bc_dict:
            self.bc_dict = bc_dict
        else:
            bc_dict = self.bc_dict

        C_list = list(symbols('C1:' + str(len(bc_dict) + 1)))

        eig_value = self.eigenvalues(bc_dict, sep_expr, arg, spatial_comp,
                                     index).formula.expand().simplify()

        mode_eqn = self.fundamental_matrix(bc_dict, sep_expr,
                                           spatial_comp) * Matrix(C_list)

        eig_aid = ({
            comp: 0
            for comp in (mode_eqn.atoms(sin, cos))
            if comp.subs(arg, eig_value).subs(index, mode_no).n() < 0.001
        })

        display(mode_eqn.subs(eig_aid))
        display(mode_eqn[1:])

        display(eig_aid)

        mode_subs = solve(
            mode_eqn.applyfunc(lambda x: x.subs(eig_aid))[:-1], C_list)
        #         mode_subs = solve(mode_eqn[1:], C_list)

        #         display(mode_subs)

        return self.spatial_general_solution(sep_expr=sep_expr, spatial_comp=spatial_comp).rhs.subs(arg, eig_value).subs(mode_subs).subs(arg, eig_value).subs({c_var: 1 for c_var in C_list}).subs(index, mode_no).n(2)\


class PlaneStressProblem:
    def __init__(self,
                 disp_func=[Function('\\mathit{u}')(Symbol('r')), 0],
                 stress_tensor=Matrix(2, 2, [
                     Function('\\sigma_r')(Symbol('r')), 0, 0,
                     Function('\\sigma_\\varphi')(Symbol('r'))
                 ]),
                 bc_dict=None,
                 coords=[Symbol('r'), Symbol('\\varphi')],
                 E=Symbol('E', positive=True),
                 nu=Symbol('\\nu', positive=True),
                 label=None,
                 system=None,
                 **kwargs):

        if system:
            disp_func = system.u
            stress_tensor = system.stress
            coords = system.coords
            E = system.E_module
            nu = system.poisson
            bc_dict = system.bc_dict

        self.u = Matrix(disp_func)
        self.stress = stress_tensor
        self.E_module = E
        self.poisson = nu
        self.coords = coords
        self.dim = max(self.u.shape)
        self._equilibrium_eqn = None
        self._integrantion_cs = None
        self._bc = Symbol('BC')

        if label == None:
            label = self.__class__.__name__ + ' on ' + str(self.coords)

        self._label = label

    @property
    def BC(self):
        return self._bc

    def __call__(self, *args, label=None):
        """
        Returns the label of the object or class instance with reduced Degrees of Freedom.
        """

        if isinstance(args[0], str):
            if label:
                self._label = label
            else:
                self._label = args[0]

            return self

    def __str__(self):

        return self._label

    def __repr__(self):

        return self.__str__()

    def subs(self, *args, **kwargs):
        E_new, nu_new, u_new, stress_new = Tuple(self.E_module, self.poisson,
                                                 self.u,
                                                 self.stress).subs(*args)

        bc_trap = self.BC.subs(*args)

        if bc_trap == self.BC:
            bc_trap = self.bc_dict

        new_system = PlaneStressProblem(disp_func=u_new,
                                        stress_tensor=stress_new,
                                        bc_dict=bc_trap,
                                        coords=self.coords,
                                        E=E_new,
                                        nu=nu_new,
                                        label=self._label,
                                        system=None,
                                        **kwargs)

        return type(self)(system=new_system)

    @property
    def _eoms(self):
        return self.governing_equation

    def strain_field(self, disp_func=None):
        if disp_func is not None:
            self.u = disp_func
            self.dim = max(self.u.shape)

        self.strain = Matrix(self.dim, self.dim, [
            Derivative(self.u[0], self.coords[0]), 0, 0,
            self.u[0] / self.coords[0]
        ])
        return self.strain

    def stress_field(self, stress_tensor=None):
        if stress_tensor is not None:
            self.stress = stress_tensor

        return self.stress

    def constitutive_equation(self, E=None, nu=None):
        if E is not None:
            self.E_module = E

        if nu is not None:
            self.poisson = nu

        return Matrix(self.dim, self.dim, [
            1 / self.E_module, -self.poisson / self.E_module,
            -self.poisson / self.E_module, 1 / self.E_module
        ])

    def strain_load_equation(self):

        return self.constitutive_equation() * Matrix(
            [self.stress[0, 0], self.stress[-1, -1]])

    def load_strain_equation(self, dict=None):

        load_mat = self.constitutive_equation().inv() * Matrix(
            [self.strain_field()[0, 0],
             self.strain_field()[-1, -1]])

        if dict is True:
            return {
                self.stress[0, 0]: load_mat[0],
                self.stress[-1, -1]: load_mat[-1]
            }
        else:
            return load_mat

    def load_equilibrium(self, volumetric_load):

        return Derivative(self.stress[0] * self.coords[0],
                          self.coords[0]) - self.stress[-1] + volumetric_load * self.coords[0]

    def equilibrium_eqn(self, volumetric_load, subs_dict=None):

        stress_r, stress_phi = list(self.load_strain_equation())

        self._equilibrium_eqn = self.load_equilibrium(
            volumetric_load=volumetric_load).doit().subs(
                self.load_strain_equation(
                    dict=True)).subs(subs_dict).expand().simplify()

        return self._equilibrium_eqn

    def solution(self, volumetric_load, dvar, ics=None, ivar=Symbol('r')):

        if ics is None:
            return dsolve(self._equilibrium_eqn, dvar)
        else:
            return dsolve(self._equilibrium_eqn, dvar, ics=ics)

    def integration_constants(self,
                              volumetric_load,
                              dvar,
                              ics,
                              ivar=Symbol('r')):

        sol = self.solution(volumetric_load=volumetric_load, dvar=dvar)

        conds_eqns = [
            Eq((key.subs(sol.lhs, sol.rhs).doit()), cond_value).doit()
            for key, cond_value in ics.items()
        ]

        self._integrantion_cs = solve(conds_eqns,
                                      symbols('C1:' + str(self.dim + 1)))

        return self._integrantion_cs

    def particular_solution(self, volumetric_load, dvar, ics,
                            ivar=Symbol('r')):

        sol = self.solution(volumetric_load=volumetric_load, dvar=dvar)

        return sol.subs(self._integrantion_cs)

    def loads(self, volumetric_load, dvar, ics, ivar=Symbol('r')):

        particular_sol = self.particular_solution(volumetric_load,
                                                  dvar,
                                                  ics,
                                                  ivar=Symbol('r'))

        return self.load_strain_equation().subs(particular_sol.lhs,
                                                particular_sol.rhs)
