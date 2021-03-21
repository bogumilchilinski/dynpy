from sympy import *
from sympy.physics.mechanics import *
from sympy.physics.vector.printing import vpprint, vlatex
import sympy as sym
from sympy.utilities.autowrap import autowrap, ufuncify
import numpy as np
import itertools as itools
import scipy.integrate as solver
from .utilities.timeseries import *

from collections import ChainMap

from IPython.display import display

import sympy.physics.mechanics as me

from sympy.simplify.fu import TR8, TR10, TR7, TR3




from .solvers.numerical import OdeComputationalCase

from .solvers.linear import LinearODESolution

from .solvers.nonlinear import WeakNonlinearProblemSolution, MultiTimeScaleMethod



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

        if label == None:
            label = self.__class__.__name__ + ' with ' + str(len(
                self.q)) + 'DOF'

        self._label = label

        #LM=me.LagrangesMethod(Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs)
    @property
    def Y(self):
        return Matrix(list(self.q) + list(self.q.diff(self.ivar)))

    def _kwargs(self):
        """
        Returns all key words arguments that an instance has
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

        print(self_dict['qs'])

        list_build = lambda x: [x] if x else []

        self_dict['forcelist'] = list_build(
            self_dict['forcelist']) + list_build(other_dict['forcelist'])
        self_dict['bodies'] = list_build(self_dict['bodies']) + list_build(
            other_dict['bodies'])
        
        if not self_dict['frame']:
            self_dict['frame']=other_dict['frame']

        return LagrangesDynamicSystem(**self_dict)

    def shranked(self, *args):
        """
        Returns class instance with reduced Degrees of Freedom
        """
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
        """
        Returns class instance with substituted numerical values
        """

        hol_coneqs=list(self._hol_coneqs)
        
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
            
            hol_coneqs=hol_coneqs+constrains
            
            

#         print(hol_coneqs)
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

        #print(type(self))
        
        
        
        return type(self)(Lagrangian=lagrangian_subs,
                          qs=self.q,
                          forcelist=forces_subs,
                          bodies=self._bodies,
                          frame=self.frame,
                          hol_coneqs=hol_coneqs,
                          nonhol_coneqs=list(
                              self.coneqs)[len((self._hol_coneqs)):],
                          ivar=self.ivar)


#     def rhs(self):
#         return self.rhs

    def __preview(self, expr, preview=None, preview_mode=None):
        '''
        Private method which processes preview printing. The method is or isn't used depending on user's bool 'preview' input.
        '''
        display(expr)
        pass
        

    def __call__(self, label=None):
        """
        Returns a label of the object
        """
        self._label = label

        return self

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
                Symbol(str(q_tmp).replace('(' + str(self.ivar) + ')', '_0'))
                for q_tmp in self.q
            }

        self.q_0 = static_disp_dict

        return self.governing_equations.subs(static_disp_dict)

    def external_forces(self):
        """
        Returns Matrix with external forces
        """
        return self.governing_equations.subs(
            {gen_coord: 0
             for gen_coord in self.Y}).doit()
    
    def reactions(self):
        multipliers=self.solve_multipliers()
        
        return [self._eoms.jacobian([lam])*value for lam,value  in multipliers.items()]
    
    def generalized_momentum(self,dict=True):
        


        if dict:
            momentum_dict = {
                q_tmp:
                Symbol('p_{'+str(q_tmp).replace('(' + str(self.ivar) + ')', '')+'}')
                for q_tmp in self.q
                }

            return {momentum_sym:self.lagrangian().diff(coord.diff(self.ivar)) for coord,momentum_sym  in momentum_dict.items()}
    
    def hamiltonian(self,dict=True):
        
        if dict:
            momentum_dict = self.generalized_momentum()
        print(momentum_dict)
        ham_dict = {Symbol('H_'+str(momentum_sym)):coord.diff(self.ivar) * self.generalized_momentum()[coord] - self.lagrangian() for coord,momentum_sym in momentum_dict.items()}
        print(ham_dict)
        return {Symbol('H'):ham_dict.subs(solve([Eq(momentum_sym,self.lagrangian().diff(coord.diff(self.ivar)))],coord.diff(self.ivar))) for coord,momentum_sym in momentum_dict.items()}

#         ham = self.q.diff(self.ivar) * self.lagrangian().diff(self.q.diff(self.ivar)) - self.lagrangian()
        
#         return {Symbol('H'):ham.subs(solve([Eq(momentum_sym,self.lagrangian().diff(coord.diff(self.ivar)))],coord.diff(self.ivar))) for coord,momentum_sym in momentum_dict.items()}
        
    def _op_points(self,
                   static_disp_dict=None,
                   dict=True,
                   subs=False,
                   hint=[],
                   *args,
                   **kwargs):
        """
        Provides the interface for critical points evaluation (solves the equlibrium conditions and returns its roots).
        """

        eqns_to_solve = self.equilibrium_equation(
            static_disp_dict=static_disp_dict).doit()
        #         display(args,kwargs)
        roots = solve(list(eqns_to_solve) + list(flatten([hint])),
                      list(self.q_0.values()),
                      dict=dict)

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

        #print('x0',x0)
        if not x0:
            x0 = {coord: 0 for coord in self.Y}

        if op_point:
            x0.update(self._op_points(hint=hint, subs=True)[0])

        lagrangian_approx = multivariable_taylor_series(self.lagrangian(),
                                                        self.Y,
                                                        n=n + 1,
                                                        x0=x0)

        return LagrangesDynamicSystem(lagrangian_approx,
                                      self.q,
                                      forcelist=self.forcelist,
                                      frame=self.frame,
                                      label=label,
                                      ivar=self.ivar)

    def linearized(self, x0=None, op_point=False, hint=[], label=None):
        """
        Returns the same result as def approximated() but only for first order functions
        """

        linearized_sys = self.approximated(n=1,
                                           op_point=op_point,
                                           hint=hint,
                                           x0=x0)

        return LinearDynamicSystem(linearized_sys.lagrangian(),
                                   self.q,
                                   forcelist=self.forcelist,
                                   frame=self.frame,
                                   label=label,
                                   ivar=self.ivar)

    @property
    def _eoms(self):
        """
        Returns Equations of Motion of the system
        """
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

        main_matrix = self.rhs().jacobian(self.Y)

        return (main_matrix).diagonalize()[1]

    def modes(self):
        '''
        Determines the system vibration modes matrix (eigenmodes) based on the system free vibration frequencies.
        '''
        #         return (((self.inertia_matrix().inv() *
        #                   self.stiffness_matrix()).diagonalize()[0]))

        main_matrix = self.rhs().jacobian(self.Y)

        return (main_matrix).diagonalize()[0][len(self.q):,:]


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
    def damped_natural_frequencies(self):
        '''
        Determines the system damped natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        damped_freqs = [sqrt(nat_freq**2 - (self.damping_coefficient()[0]/2)**2) for nat_freq in self.natural_frequencies()]

        return diag(*damped_freqs)
    
    def im_eigenvals(self):
        '''
        Determines the system damped natural frequencies matrix (in the diagonal form). Output is obtained from inertia matrix and stiffness matrix.
        '''
        
        natural_freqs=list({(im((eigen).doit())) for eigen  in self.eigenvalues() if not eigen==0 })
        
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
        
        log_dec = [(2*pi*self.damping_coefficient()[0]/2) / damp_freq for damp_freq in self.damped_natural_frequencies()]
        
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

        eoms = self._eoms
        subs_dict={}
        if len(self.q)==1:
            eoms = eoms-self.stiffness_matrix()*Matrix(self.q)+ ((Symbol('omega_h',positive=True)**2+(self.damping_coefficient()[0]/2)**2)*self.inertia_matrix()*Matrix(self.q) )
            subs_dict={Symbol('omega_h',positive=True):sqrt(self.natural_frequencies()[0]**2-(self.damping_coefficient()[0]/2)**2)}
#             print('len',len(self.q))
#             display(subs_dict)
                    
        return LinearODESolution(eoms,ivar=self.ivar,dvars=self.q).general_solution(
            initial_conditions=initial_conditions).subs(subs_dict).doit()

    def steady_solution(self, initial_conditions=None):
        """
        
        """
        return LinearODESolution(self._eoms,ivar=self.ivar,dvars=self.q).steady_solution(
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

        general_solution = self.steady_solution()[0]
        
        comp_sin = general_solution.coeff(sin(excitation_freq * self.ivar)) 
        comp_cos = general_solution.coeff(cos(excitation_freq * self.ivar))        

        n_sin, d = fraction(comp_sin)
        n_cos, d = fraction(comp_cos)

#         print(n_sin)
#         print(n_cos)
#         print(d)

        if len(self.q) == 1:  # if single degree of freedom
            frf_expr = ((sqrt((n_sin**2 + n_cos**2).simplify())) /d.doit()).simplify()   # sDoF
        else: # DoF > 1
            frf_expr = inv( stiffness_matrix() - excitation_freq**2 * inertia_matrix() ) # mDoF

        return frf_expr

    def dynamic_amplification_factor(self,
                                     excitation_amp=None,
                                     excitation_freq=Symbol('Omega',positive=True)):
        '''
        Returns the Dynamic Amplification Factor of the system for the given excitation amplitude (working correctly for single degree of freedom systems).
        '''
        frf = self.frequency_response_function(excitation_freq=excitation_freq)
        
#         display(self.external_forces())
        
        sin_comp = self.external_forces()[0].coeff(sin(excitation_freq * self.ivar))
        cos_comp = self.external_forces()[0].coeff(cos(excitation_freq * self.ivar))

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



