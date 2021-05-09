from sympy import (flatten, SeqFormula, Function,
                   Symbol, symbols, Eq, Matrix, S, oo, dsolve)

from sympy.physics.vector.printing import vpprint, vlatex


class ContinuousSystem:

    def __init__(self, L, q,bc_dict=None, t_var=Symbol('t'), spatial_var=Symbol('x'), derivative_order=2,label=None,**kwargs):
        self.t = t_var
        self.L = L
        self.q = q
        self.r = flatten((spatial_var,))
        self.diff_ord = derivative_order
        self.bc_dict=bc_dict
        self._bc=Symbol('BC')


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
                self._label=label
            else:
                self._label=args[0]
                
                
            return self


    def __str__(self):

        return self._label

    def __repr__(self):

        return self.__str__()

    def subs(self,*args):
        L_new=self.L.subs(*args)
        q_new=self.q.subs(*args)

        bc_trap=self.BC.subs(*args)

        if bc_trap==self.BC:
            bc_trap=self.bc_dict

        return type(self)( L=L_new, q=q_new,bc_dict=bc_trap, t_var=self.t, spatial_var=self.r, derivative_order=self.diff_ord,label=self._label)


    def inertia_force(self):
        q = self.q
        t = self.t
        L = self.L

        return L.diff(q.diff(t)).diff(t)

    def restoring_force(self):
        q = self.q
        t = self.t
        L = self.L

        return sum([(-1)**(order+2)*L.diff(q.diff(r_var, order+1)).diff(r_var, order+1) for r_var in self.r for order in range(self.diff_ord)])-L.diff(q)

    def governing_equation(self):

        return self.inertia_force()+self.restoring_force()

    def eom_coeff(self, expr):

        return self.governing_equation().coeff(expr)

    def apply_separation(self, time_comp=Function('T')(Symbol('t')), spatial_comp=Function('X')(Symbol('x'))):

        return self.governing_equation().subs(self.q, time_comp*spatial_comp).doit()

    def separated_vars_eqn(self, time_comp=Function('T')(Symbol('t')), spatial_comp=Function('X')(Symbol('x'))):

        eqn_with_subs = self.apply_separation(time_comp, spatial_comp)

        return Eq((eqn_with_subs.coeff(spatial_comp)/(time_comp)), -eqn_with_subs.coeff(time_comp)/spatial_comp)

    def spatial_eqn(self, sep_expr, spatial_comp=Function('X')(Symbol('x'))):

        separated_eqn_rhs = self.separated_vars_eqn(
            spatial_comp=spatial_comp).rhs

        return Eq(separated_eqn_rhs, sep_expr)

    def time_eqn(self, sep_expr, time_comp=Function('T')(Symbol('t'))):

        separated_eqn_lhs = self.separated_vars_eqn(time_comp=time_comp).lhs

        return Eq(separated_eqn_lhs, sep_expr)

    def spatial_general_solution(self, sep_expr, spatial_comp=Function('X')(Symbol('x'))):

        spatial_ode = self.spatial_eqn(sep_expr, spatial_comp)

        return dsolve(spatial_ode, spatial_comp)

    def fundamental_matrix(self, bc_dict=None, sep_expr=Symbol('k'), spatial_comp=Function('X')(Symbol('x'))):

        if bc_dict:
            self.bc_dict=bc_dict
        else:
            bc_dict=self.bc_dict

        spatial_sol = self.spatial_general_solution(
            sep_expr=sep_expr, spatial_comp=spatial_comp)

        fun_eqns = Matrix([spatial_sol.rhs.subs(spatial_comp.args[0], flatten(
            [key.args[-1]])[0])-val for key, val in bc_dict.items()])

        matrix_comps_list = []

        for key, val in bc_dict.items():

            #display(list(key.atoms(Symbol,Number) - spatial_comp.atoms(Symbol)))
            free_sym = list(key.atoms(Symbol, Number) -
                            spatial_comp.atoms(Symbol))[-1]
            #print({spatial_sol.lhs:spatial_sol.rhs, spatial_sol.lhs.subs(spatial_comp.args[0],free_sym):spatial_sol.rhs, })
            matrix_comps_list += [(key.subs({spatial_sol.lhs: spatial_sol.rhs, spatial_sol.lhs.subs(
                spatial_comp.args[0], free_sym):spatial_sol.rhs.subs(spatial_comp.args[0], free_sym), }).doit()-val)]

        fun_eqns = Matrix(matrix_comps_list)

        return fun_eqns.jacobian(symbols('C1:'+str(len(bc_dict)+1)))


    def char_poly(self, bc_dict=None, sep_expr=Symbol('k'), spatial_comp=Function('X')(Symbol('x'))):

        if bc_dict:
            self.bc_dict=bc_dict
        else:
            bc_dict=self.bc_dict

        return self.fundamental_matrix(bc_dict, sep_expr, spatial_comp).det().simplify()

    def eigenvalues(self, bc_dict=None, sep_exp=Symbol('k'), arg=Symbol('k'), spatial_comp=Function('X')(Symbol('x')), index=Symbol('n', integer=True, positive=True)):
        
        if bc_dict:
            self.bc_dict=bc_dict
        else:
            bc_dict=self.bc_dict

        root = solve(self.char_poly(bc_dict, sep_expr, spatial_comp), arg)[0]

        spatial_span = list(root.free_symbols)[0]

        return SeqFormula(root+(index-1)/spatial_span*pi, (index, 0, oo))

    def eigenmodes(self, mode_no, bc_dict=None, sep_expr=Symbol('k'), arg=Symbol('k'), spatial_comp=Function('X')(Symbol('x')), index=Symbol('n', integer=True, positive=True)):

        if bc_dict:
            self.bc_dict=bc_dict
        else:
            bc_dict=self.bc_dict

        C_list = list(symbols('C1:'+str(len(bc_dict)+1)))

        eig_value = self.eigenvalues(
            bc_dict, sep_expr, arg, spatial_comp, index).formula.expand().simplify()

        mode_eqn = self.fundamental_matrix(
            bc_dict, sep_expr, spatial_comp).subs(arg, eig_value)*Matrix(C_list)

        mode_subs = solve(mode_eqn[:-1], C_list)

        return self.spatial_general_solution(sep_expr=sep_expr, spatial_comp=spatial_comp).rhs.subs(mode_subs).subs(k, eig_value).subs({c_var: 1 for c_var in C_list}).subs(index, mode_no)

