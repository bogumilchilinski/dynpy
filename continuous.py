from numpy import positive
from sympy import (flatten, SeqFormula, Function, Symbol, symbols, Eq, Matrix,
                   S, oo, dsolve, solve, Number, pi, cos, sin, Tuple,
                   Derivative,sqrt)

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
                 time_comp=Function('T')(Symbol('t')),
                 spatial_comp=Function('X')(Symbol('x')),
                 mode_symbol=Symbol('k',positive=True),
                 **kwargs):

        if system:
            t_var = system.t
            L = system.L
            q = system.q
            spatial_var = system.r
            derivative_order = system.diff_ord
            bc_dict = system.bc_dict
            time_comp = system.time_comp
            spatial_comp = system.spatial_comp
            mode_symbol = system._mode_symbol



        self.t = t_var
        self.L = L
        self.q = q
        self.r = flatten((spatial_var, ))
        self.diff_ord = derivative_order
        self.bc_dict = bc_dict
        self._bc = Symbol('BC')
        self._mode_symbol = mode_symbol

        self.time_comp=time_comp
        self.spatial_comp=spatial_comp

        self._sep_expr = Symbol('k', positive=True)

        if label == None:
            label = self.__class__.__name__ + ' on ' + str(self.q)

        self._label = label
        self._given_data={}
        
        self._nonlinear_base_system=None
        

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
        
        given_data=args[0]
        
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
        
        new_sys = type(self)(L_new, q_new, system=new_system)
        new_sys._given_data=given_data
        new_sys._nonlinear_base_system = self._nonlinear_base_system

        return new_sys

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
                    arg=None,
                    spatial_comp=Function('X')(Symbol('x')),
                    index=Symbol('n', integer=True, positive=True)):

        if not sep_expr:
            sep_expr = self._sep_expr

        if bc_dict:
            self.bc_dict = bc_dict
        else:
            bc_dict = self.bc_dict

        if arg is None:
            arg = self._mode_symbol


        roots = solve(self.char_poly(bc_dict, sep_expr, spatial_comp), arg)
#         print('+++++++++++++ roots from eigs')
#         print(roots)

        if len(roots) == 1:
            spatial_span = roots[0]
            first_root=roots[0]
        else:
            spatial_span = roots[1] - roots[0]
            first_root=roots[0]
            
        seq=SeqFormula(first_root + (index -1) * spatial_span,
                          (index, 1, oo))

#         display(seq)
        
        return seq

    def eigenmodes(self,
                   mode_no,
                   bc_dict=None,
                   sep_expr=None,
                   arg=None,
                   spatial_comp=Function('X')(Symbol('x')),
                   index=Symbol('n', integer=True, positive=True)):

        if arg is None:
            arg = self._mode_symbol

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
        #display( eig_value)
        #display( eig_value)

        eig_aid = ({
            comp: 0
            for comp in (mode_eqn.atoms(sin, cos))
            if comp.subs(arg, eig_value).subs(index, mode_no).n() < 0.01
        })

#         print('mode_eqn')
#         display(mode_eqn.subs(eig_aid))
#         display(mode_eqn)
#         print('aid')
#         display(eig_aid)

        mode_subs = solve(
            mode_eqn.applyfunc(lambda x: x.subs(eig_aid)), C_list)
        #         mode_subs = solve(mode_eqn[1:], C_list)

#         print('solve for Cs')
#         display(mode_subs)

        

        gen_sol=self.spatial_general_solution(sep_expr=sep_expr, spatial_comp=spatial_comp).rhs
#         display(gen_sol)
        
        gen_sol=gen_sol.subs(arg, eig_value).subs(eig_aid)
#         display(gen_sol)
        
        return gen_sol.subs(mode_subs).subs(arg, eig_value).subs({c_var: 1 for c_var in C_list}).subs(index, mode_no).n(2)

    def phase_velocity(self):


        X_fun=self.spatial_comp
        T_fun=self.time_comp

        nat_freq=Symbol('omega',positive=True)

        #display(self.time_eqn())

        nat_freq_eq=self.time_eqn().subs({T_fun.diff(self.t,2) :-T_fun*nat_freq **2 })
        
        nat_freq_expr = solve( nat_freq_eq.lhs - nat_freq_eq.rhs, nat_freq **2 )

        return sqrt(max(nat_freq_expr))/self._mode_symbol


class PlaneStressProblem:
    def __init__(self,
                 #disp_func=[Function('\\mathit{u}')(Symbol('r')), 0],
                 disp_func,
                 stress_tensor=Matrix(2, 2, [
                     Function('\\sigma_r')(Symbol('r')), 0, 0,
                     Function('\\sigma_\\varphi')(Symbol('r'))
                 ]),
                 bc_dict=None,
                 coords=[Symbol('r'), Symbol('\\varphi')],
                 E=Symbol('E', positive=True),
                 nu=Symbol('\\nu', positive=True),
                 D=Symbol('D', positive=True),                 
                 label=None,
                 system=None,
                 volumetric_load=None,
                 **kwargs):

        if system:
            disp_func = system.u
            stress_tensor = system.stress
            coords = system.coords
            E = system.E_module
            nu = system.poisson
            D=system.D
            bc_dict = system.bc_dict
            volumetric_load=system.volumetric_load
            
#         print('__init')
#         print(type(self))
#         display(disp_func)

        self.u = Matrix(disp_func)
        
#         print('__init')
#         print(type(self))
#         display(self.u)
        
        
        self.stress = stress_tensor
        self.E_module = E
        self.poisson = nu
        self.D = D
        self.coords = coords
        self.dim = max(self.u.shape)
        self._equilibrium_eqn = None
        self._integrantion_cs = None
        self._bc = Symbol('BC')
        self.bc_dict = bc_dict
        self.volumetric_load = volumetric_load
        #print('__init__',self.bc_dict)
        if label == None:
            label = self.__class__.__name__ + ' on ' + str(self.coords)

        self._label = label
        self._subs_dict = {}
        self._given_data={}

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
        
#         print('u_old')
#         display(self.u)
        given_data=args[0]
    
        E_new, nu_new, u_new, stress_new, vol_load_new, D_new = Tuple(
            self.E_module, self.poisson, self.u, self.stress,
            self.volumetric_load,self.D).subs(*args)


        bc_trap = self.BC.subs(*args)

        if bc_trap == self.BC:
            bc_trap = self.bc_dict


        new_system = type(self)(disp_func=u_new,
                                        stress_tensor=stress_new,
                                        bc_dict=bc_trap,
                                        volumetric_load=vol_load_new,
                                        coords=self.coords,
                                        E=E_new,
                                        nu=nu_new,
                                        D=D_new,
                                        label=self._label,
                                        system=None,
                                        **kwargs)

        
        new_sys = type(self)(system=new_system)
        new_sys._given_data=given_data
        
        return new_sys

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

    def load_equilibrium(self, volumetric_load=None):

#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        return Derivative(
            self.stress[0] * self.coords[0], self.coords[0]
        ) - self.stress[-1] + self.volumetric_load * self.coords[0]

    def equilibrium_eqn(self, volumetric_load=None, subs_dict=None):

        if subs_dict:
            self._subs_dict = subs_dict

#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        stress_r, stress_phi = list(self.load_strain_equation())

        
        
        self._equilibrium_eqn = self.load_equilibrium(
            volumetric_load=self.volumetric_load).doit().subs(
                self.load_strain_equation(dict=True)).subs(
                    self._subs_dict).expand().simplify()

        return self._equilibrium_eqn

    def solution(self, volumetric_load, dvar, ics=None, ivar=Symbol('r')):


        #         display(ics)
#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        
            
        if not self._equilibrium_eqn:
            self._equilibrium_eqn=self.equilibrium_eqn( volumetric_load=self.volumetric_load, subs_dict=self._subs_dict)

        C_list=symbols('C1:' + str(self.dim + 1))
            
        if ics is None:
            #display(self._equilibrium_eqn)
            
            return (Eq(dvar,-(ivar*(self._equilibrium_eqn.subs(dvar,0)/ivar).doit().integrate(ivar)).integrate(ivar)/self.D/ivar  + C_list[0]/ivar + C_list[1]*ivar  )  )
            
            #return dsolve(self._equilibrium_eqn, dvar)
        else:
            return dsolve(self._equilibrium_eqn, dvar, ics=ics)

    def integration_constants(self,
                              volumetric_load,
                              dvar,
                              ics,
                              ivar=Symbol('r')):
        
#         display('vol',self.volumetric_load)

#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        sol = self.solution(volumetric_load=self.volumetric_load, dvar=dvar)

        conds_eqns = [
            ((key.subs(sol.lhs, sol.rhs).doit()) - cond_value).doit()
            for key, cond_value in ics.items()
        ]

#         display(conds_eqns)
        
        self._integrantion_cs = solve(conds_eqns,
                                      symbols('C1:' + str(self.dim + 1)))

        return self._integrantion_cs

    def particular_solution(self, volumetric_load, dvar, ics,
                            ivar=Symbol('r')):

#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        sol = self.solution(volumetric_load=self.volumetric_load, dvar=dvar)

        if not self._integrantion_cs:
            self._integrantion_cs=self.integration_constants(
                              volumetric_load=self.volumetric_load,
                              dvar=dvar,
                              ics=ics,
                              ivar=ivar)
        
        return sol.subs(self._integrantion_cs)

    def loads(self, volumetric_load, dvar, ics, ivar=Symbol('r')):

#         if volumetric_load:
#             self.volumetric_load = volumetric_load

        particular_sol = self.particular_solution(self.volumetric_load,
                                                  dvar,
                                                  ics,
                                                  ivar=Symbol('r'))

        return self.load_strain_equation().subs(particular_sol.lhs,
                                                particular_sol.rhs)
