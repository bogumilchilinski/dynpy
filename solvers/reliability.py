from sympy import Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, diag, Eq, Dict, ImmutableMatrix, latex, solve#, Tuple



# Amadeusz Radomski

class Reliability():
    '''
    Class bases on First Order Reliability Method. By asumming mean values and standard deviations of the parameters, the class calculates new design point parameters until reliabiity index convergence.
    '''

    def __init__(self,
                limit_state_function,
                mean_values_dict,
                standard_deviations_dict,
                parameters_list,
                dependent_variable=None,
                new_design_point_coords=None,
                ):

        self.limit_state_function = limit_state_function
        self.mean_values_dict = mean_values_dict
        self.parameters_list = parameters_list
        self.standard_deviations_dict = standard_deviations_dict

        if dependent_variable is not None: self.dependent_variable = dependent_variable
        else: self.dependent_variable = self.parameters_list[-1]

        if new_design_point_coords is not None:

            dictionary={}
            for i, param in enumerate(self.parameters_list):
                dictionary[param] = new_design_point_coords[i]

            self.new_design_point_coords = dictionary

        else:
            self.new_design_point_coords = mean_values_dict



    def dependent_variable_dict(self, iteration=None):
        '''
        Method returns dict containing calulated dependent parameter expression in original space.
        '''

        mean_val_dict = self.new_design_point_coords.copy()

        func = Eq(0, self.limit_state_function.rhs, evaluate=False)

        func_sol = (solve(func, self.dependent_variable)[0])#subs(mean_val_dict)

#         func_sol_eq = Eq(self.dependent_variable, func_sol, evaluate=False)

        mean_val_dict[self.dependent_variable] = func_sol

        return mean_val_dict

    def design_point_in_original_space(self, iteration=None):
        '''
        Method defines coordinates of the design point in the original space.
        '''

        x_star=Symbol(f'x__*',positive=True)

        mean_val_dict = self.dependent_variable_dict()

        design_point_mat = Matrix([self.parameters_list]).subs(mean_val_dict)

        eq1 = Eq(x_star, design_point_mat, evaluate=False)

        eq2 = Eq(x_star, eq1.rhs.subs(mean_val_dict).n(), evaluate=False)

        return eq2


    def design_point_in_reduced_space(self, iteration=None):
        '''
        Method defines coordinates of the design point in the reduced space.
        '''

        x_star_prim=Symbol('''x__'*''',positive=True)

        dp_in_og_space = self.design_point_in_original_space().rhs

        coords_list=[]

        for i,coord in enumerate(dp_in_og_space):

            reduced_coord = ((coord-self.parameters_list[i].subs(self.mean_values_dict))/self.parameters_list[i].subs(self.standard_deviations_dict)).n()
            coords_list.append(reduced_coord)

        design_point_mat = Matrix([coords_list])

        eq = Eq(x_star_prim, design_point_mat.T, evaluate=False)

        return eq

    def design_point_derivatives_vector(self, iteration=None):
        '''
        Method defines vector containing the derivatives of the performance function evaluated at the design point.
        '''
        A=Symbol('A',positive=True)

        mean_val_dict = self.dependent_variable_dict()

        diff_list=[]

        for i,parameter in enumerate(self.parameters_list):

            calulated_diff = (self.limit_state_function.rhs.diff(parameter).subs(mean_val_dict).subs(mean_val_dict)*parameter.subs(self.standard_deviations_dict)).n()
            diff_list.append(calulated_diff)

        diff_design_point_vector = Matrix([diff_list])

        eq = Eq(A, diff_design_point_vector.T, evaluate=False)

        return eq

    def estimate_reliability_index(self, iteration=None):
        '''
        Method estimates the reliability index.
        '''

        beta_hl=Symbol('beta_HL')

        A=self.design_point_derivatives_vector().rhs
        x_prim_star = self.design_point_in_reduced_space().rhs

        eq = Eq(beta_hl, ((-(A.T*x_prim_star)[0])/sqrt((A.T*A)[0,0])).n(), evaluate=False)

        return eq

    def vector_of_directional_cosines(self, iteration=None):
        '''
        Method estimates the vector of directional cosines.
        '''

        A=self.design_point_derivatives_vector().rhs

        alpha=Symbol('alpha',positive=True)

        directional_cosines_vector = (A/(sqrt(A.T * A))[0,0]).n()

        eq = Eq(alpha, directional_cosines_vector, evaluate=False)

        return eq

    def new_design_point_in_reduced_space(self, iteration=None):
        '''
        Method estimates the new design point in reduced space for n-1 variables.
        '''

        alpha= self.vector_of_directional_cosines().rhs
        beta_hl = self.estimate_reliability_index().rhs

        x_star_prim=Symbol('''x__'*''',positive=True)


        new_design_point_list=[]

        for i,parameter in enumerate(self.parameters_list):

            new_design_point = (-alpha[i,0]*beta_hl).n()
            new_design_point_list.append(new_design_point)

        new_design_point_vector = Matrix([new_design_point_list])

        new_design_point_vector.col_del(i)

        eq = Eq(x_star_prim, new_design_point_vector.T, evaluate=False)

        return eq

    def new_design_point_in_original_space(self, iteration=None):
        '''
        Method estimates the new design point in original space.
        '''

        new_design_point_in_reduced_space = self.new_design_point_in_reduced_space().rhs

        x_star=Symbol('''x__*''',positive=True)

        new_design_point_list=[]
        design_point_dict={}

        for i,parameter in enumerate(self.parameters_list):

            if i < len(new_design_point_in_reduced_space):
                new_design_point = (parameter.subs(self.mean_values_dict)+new_design_point_in_reduced_space[i,0]*parameter.subs(self.standard_deviations_dict)).n()
                new_design_point_list.append(new_design_point)
                design_point_dict[parameter] = new_design_point
            else:
                break


        new_design_point_vector = Matrix([new_design_point_list])

        sol = (self.dependent_variable.subs(self.dependent_variable_dict()).subs(design_point_dict)).n()

        new_design_point_vector = new_design_point_vector.col_insert(i, Matrix([[sol]]))


        eq = Eq(x_star, new_design_point_vector.T, evaluate=False)

        return eq

    def next_step(self):

        lsf=self.limit_state_function

        return Reliability(limit_state_function=lsf,
                mean_values_dict=self.mean_values_dict,
                standard_deviations_dict=self.standard_deviations_dict,
                parameters_list=self.parameters_list,
                dependent_variable=self.dependent_variable,
                new_design_point_coords=self.new_design_point_in_original_space().rhs,)


    def search_for_stable_index(self, n=None, treshold=None):

        i=0
        if n is None and treshold is None: n=5
        elif treshold is not None: n = 100

        display(f'iteracja {i+1}')

        list_for_coords = []
        list_for_param = []
        list_for_new_mdp = []

        og_space = self.design_point_in_original_space()
        display(og_space)
#         display(self.design_point_in_reduced_space())
#         display(self.design_point_derivatives_vector())
        rel_index_initial=self.estimate_reliability_index()
        display(rel_index_initial)
#         display(self.vector_of_directional_cosines())
#         display(self.new_design_point_in_reduced_space())

        new_og_space = self.new_design_point_in_original_space()
        display(new_og_space)

        list_for_coords.append(og_space.rhs)
        list_for_param.append(rel_index_initial.rhs)
        list_for_new_mdp.append(new_og_space.rhs)

        new_reliability_instance = self.next_step()

        for i in range(n-1):
            display(f'iteracja {i+2}')

            og_space = new_reliability_instance.design_point_in_original_space()
            display(og_space)
            rel_index = new_reliability_instance.estimate_reliability_index()
            display(rel_index)
            new_og_space = new_reliability_instance.new_design_point_in_original_space()
            display(new_og_space)

            list_for_coords.append(og_space.rhs)
            list_for_param.append(rel_index.rhs)
            list_for_new_mdp.append(new_og_space.rhs)

            if abs(rel_index.rhs - rel_index_initial.rhs) <= treshold:
                display(f'Różnica dwóch ostatnich indeksów niezawodności jest mniejsza bądź równa progowi równemu {treshold}')
                break

            new_reliability_instance = new_reliability_instance.next_step()
            rel_index_initial = rel_index

            i=i+1


        return {'Original space':list_for_coords, 'Estimated reliability index':list_for_param, 'New original space':list_for_new_mdp}


    def mean_values_sensivity_analysis(self, mpp, pf, beta):

        list_for_sens=[]

        for i,parameter in enumerate(self.parameters_list):

            sens_symbol = Symbol(f's_\mu {parameter}',positive=True)

            sensitivity = pf*mpp.T[0,i]/(beta*parameter.subs(self.standard_deviations_dict))

            eq=Eq(sens_symbol, sensitivity, evaluate=False)

            display(eq)

            list_for_sens.append([eq])

        return list_for_sens

    def standard_deviations_sensivity_analysis(self, mpp, pf, beta):

        list_for_sens=[]

        for i,parameter in enumerate(self.parameters_list):

            sens_symbol = Symbol(f's_\sigma {parameter}',positive=True)

            sensitivity = pf*(mpp.T[0,i])**2/(beta*parameter.subs(self.standard_deviations_dict))

            eq=Eq(sens_symbol, sensitivity, evaluate=False)

            display(eq)

            list_for_sens.append([eq])

        return list_for_sens
    
class InverseReliability(Reliability):
    '''
    Class bases on First Order Reliability Method. By asumming expected reliability index and input, the class calculates dependent variable value.
    '''

    def __init__(self,
                limit_state_function,
                mean_values_dict,
                standard_deviations_dict,
                parameters_list,
                assumed_reliability_index,
                dependent_variable=None,
                new_design_point_coords=None,
                ):

        self.limit_state_function = limit_state_function
        self.mean_values_dict = mean_values_dict
        self.parameters_list = parameters_list
        self.standard_deviations_dict = standard_deviations_dict
        self.assumed_reliability_index = assumed_reliability_index

        if dependent_variable is not None: self.dependent_variable = dependent_variable
        else: self.dependent_variable = self.parameters_list[-1]

        if new_design_point_coords is not None:

            dictionary={}
            for i, param in enumerate(self.parameters_list):
                dictionary[param] = new_design_point_coords[i]

            self.new_design_point_coords = dictionary

        else:
            self.new_design_point_coords = mean_values_dict

    def estimate_value(self, deviation=False):
        '''
        Method estimates the parameter value.
        '''

        beta_hl=Symbol('beta_HL')

        if deviation is False:
#             symbol = self.dependent_variable.subs(self.new_design_point_coords)
            symbol = self.dependent_variable.subs(self.mean_values_dict)
        else:
            symbol = self.dependent_variable.subs(self.standard_deviations_dict)

        eq = self.estimate_reliability_index()
#         display(eq)
#         display(self.new_design_point_coords)
#         display(symbol)
        sol = solve(eq, symbol)[0].subs({beta_hl: self.assumed_reliability_index})
    
#         if sol < 0:
#             sol = 0

        return Eq(symbol, sol, evaluate=False)

    def next_iteration(self):

        coords = self.new_design_point_in_original_space().rhs.subs(self.estimate_value().lhs,self.estimate_value().rhs).T

#         display(coords)
#         display(type(coords))

        new_design_point_list = []
        for i,parameter in enumerate(coords):

            if i < len(coords):
                new_design_point = coords[0,i]
                new_design_point_list.append(new_design_point)
            else:
                break

        new_design_point_vector = Matrix([new_design_point_list])

        variable_to_append = self.dependent_variable.subs(self.new_design_point_coords)

        new_design_point_vector = new_design_point_vector.col_insert(i, Matrix([[variable_to_append]]))

#         coords[-1] = self.dependent_variable.subs(self.new_design_point_coords)

        next_iteration_step = InverseReliability(limit_state_function=self.limit_state_function,
                                         mean_values_dict=self.mean_values_dict,
                                         standard_deviations_dict=self.standard_deviations_dict,
                                         parameters_list=self.parameters_list,
                                         dependent_variable=self.dependent_variable,
                                         assumed_reliability_index=self.assumed_reliability_index,
                                         new_design_point_coords=new_design_point_vector)

        return next_iteration_step


    def search_for_stable_parameter(self, n=None, treshold=None):

        i=0
        if n is None and treshold is None: n=5
        elif treshold is not None: n = 100

        display(f'iteracja {i+1}')

        list_for_coords = []
        list_for_param = []
        list_for_new_mdp = []

        og_space = self.design_point_in_original_space()

#         display(self.design_point_in_reduced_space())
#         display(self.design_point_derivatives_vector())
        param_value_initial=self.estimate_value()
        display(param_value_initial)

#         display(self.vector_of_directional_cosines())
#         display(self.new_design_point_in_reduced_space())

        new_og_space = self.new_design_point_in_original_space().rhs.subs(self.estimate_value().lhs,self.estimate_value().rhs)
        new_og_space = Eq(self.new_design_point_in_original_space().lhs, new_og_space, evaluate=False)

        list_for_coords.append(og_space.rhs)
        list_for_param.append(param_value_initial.rhs)
        list_for_new_mdp.append(new_og_space.rhs)

        display(new_og_space)

        new_inverse_reliability_instance = self.next_iteration()

        for i in range(n-1):
            display(f'iteracja {i+2}')

            og_space = new_inverse_reliability_instance.design_point_in_original_space()#.rhs.subs(inv_form_test_y.estimate_value().lhs,inv_form_test_y.estimate_value().rhs)
#             og_space = Eq(new_reliability_instance.new_design_point_in_original_space().lhs, og_space, evaluate=False)
            display(og_space)
            param_value=new_inverse_reliability_instance.estimate_value()
            display(param_value)
            new_og_space = new_inverse_reliability_instance.new_design_point_in_original_space().rhs.subs(new_inverse_reliability_instance.estimate_value().lhs,new_inverse_reliability_instance.estimate_value().rhs)
            new_og_space = Eq(new_inverse_reliability_instance.new_design_point_in_original_space().lhs, new_og_space, evaluate=False)
            display(new_og_space)

            list_for_coords.append(og_space.rhs)
            list_for_param.append(param_value.rhs)
            list_for_new_mdp.append(new_og_space.rhs)

            if abs(param_value.rhs-param_value_initial.rhs) <= treshold:
                display(f'Różnica dwóch ostatnich wartości parametrów jest mniejsza bądź równa progowi równemu {treshold}')
                break

            new_inverse_reliability_instance = new_inverse_reliability_instance.next_iteration()
            param_value_initial = param_value

            i=i+1

        return {'Original space':list_for_coords, 'Estimated parameter':list_for_param, 'New original space':list_for_new_mdp}

    def validate(self, n=None, treshold=None):
        '''
        Method validates the result for obtaining assumed reliability index.
        '''

        mean_values_dict = {}
        param_list = self.parameters_list

        for i,parameter in enumerate(param_list):

            if i < len(param_list):
                mean_values_dict[parameter] = param_list[i].subs(self.new_design_point_coords)
            else:
                break


        mean_values_dict[self.dependent_variable] = self.estimate_value().rhs



        validate_algorithm = Reliability(limit_state_function=self.limit_state_function,
                                         mean_values_dict=mean_values_dict,
                                         standard_deviations_dict=self.standard_deviations_dict,
                                         parameters_list=self.parameters_list,
                                         dependent_variable=self.dependent_variable)

        return validate_algorithm.search_for_stable_index(n=n, treshold=treshold)

class InverseReliabilityDeviation(InverseReliability):
    '''
    Class bases on First Order Reliability Method. By asumming expected reliability index and input, the class calculates dependent variable value.
    '''

    def __init__(self,
                limit_state_function,
                mean_values_dict,
                standard_deviations_dict,
                parameters_list,
                assumed_reliability_index,
                dependent_variable=None,
                new_design_point_coords=None,
                ):

        self.limit_state_function = limit_state_function
        self.mean_values_dict = mean_values_dict
        self.parameters_list = parameters_list
        self.standard_deviations_dict = standard_deviations_dict
        self.assumed_reliability_index = assumed_reliability_index

        if dependent_variable is not None: self.dependent_variable = dependent_variable
        else: self.dependent_variable = self.parameters_list[-1]

        if new_design_point_coords is not None:

            dictionary={}
            for i, param in enumerate(self.parameters_list):
                dictionary[param] = new_design_point_coords[i]

            self.new_design_point_coords = dictionary

        else:
            self.new_design_point_coords = mean_values_dict

    def estimate_value(self, deviation=True):
        '''
        Method estimates the parameter value.
        '''

        beta_hl=Symbol('beta_HL')

        if deviation is False:
#             symbol = self.dependent_variable.subs(self.new_design_point_coords)
            symbol = self.dependent_variable.subs(self.mean_values_dict)
        else:
            symbol = self.dependent_variable.subs(self.standard_deviations_dict)

        eq = self.estimate_reliability_index()
#         display(eq)
#         display(self.new_design_point_coords)
#         display(symbol)
        sol = solve(eq, symbol)[1].subs({beta_hl: self.assumed_reliability_index})
    
#         if sol < 0:
#             sol = 0

        return Eq(symbol, sol, evaluate=False)

    def next_iteration(self):

        coords = self.new_design_point_in_original_space().rhs.subs(self.estimate_value().lhs,self.estimate_value().rhs).T

#         display(coords)
#         display(type(coords))

        new_design_point_list = []
        for i,parameter in enumerate(coords):

            if i < len(coords):
                new_design_point = coords[0,i]
                new_design_point_list.append(new_design_point)
            else:
                break

        new_design_point_vector = Matrix([new_design_point_list])

        variable_to_append = self.dependent_variable.subs(self.new_design_point_coords)

        new_design_point_vector = new_design_point_vector.col_insert(i, Matrix([[variable_to_append]]))

#         coords[-1] = self.dependent_variable.subs(self.new_design_point_coords)

        next_iteration_step = InverseReliabilityDeviation(limit_state_function=self.limit_state_function,
                                         mean_values_dict=self.mean_values_dict,
                                         standard_deviations_dict=self.standard_deviations_dict,
                                         parameters_list=self.parameters_list,
                                         dependent_variable=self.dependent_variable,
                                         assumed_reliability_index=self.assumed_reliability_index,
                                         new_design_point_coords=new_design_point_vector)

        return next_iteration_step


    def search_for_stable_parameter(self, n=None, treshold=None):

        i=0
        if n is None and treshold is None: n=5
        elif treshold is not None: n = 100

        display(f'iteracja {i+1}')

        list_for_coords = []
        list_for_param = []
        list_for_new_mdp = []

        og_space = self.design_point_in_original_space()
        display(og_space)
        
        list_for_coords.append(og_space.rhs)
#         display(self.design_point_in_reduced_space())
#         display(self.design_point_derivatives_vector())
        param_value_initial=self.estimate_value()
        display(param_value_initial)
        
        list_for_param.append(param_value_initial.rhs)
#         display(self.vector_of_directional_cosines())
#         display(self.new_design_point_in_reduced_space())

        new_og_space = self.new_design_point_in_original_space().rhs.subs(self.estimate_value().lhs,self.estimate_value().rhs)
        new_og_space = Eq(self.new_design_point_in_original_space().lhs, new_og_space, evaluate=False)
        display(new_og_space)

        list_for_new_mdp.append(new_og_space.rhs)
        
        new_inverse_reliability_instance = self.next_iteration()
        
        


        for i in range(n-1):
            display(f'iteracja {i+2}')

            og_space = new_inverse_reliability_instance.design_point_in_original_space()#.rhs.subs(inv_form_test_y.estimate_value().lhs,inv_form_test_y.estimate_value().rhs)
#             og_space = Eq(new_reliability_instance.new_design_point_in_original_space().lhs, og_space, evaluate=False)
            list_for_coords.append(og_space.rhs)
            display(og_space)
            param_value=new_inverse_reliability_instance.estimate_value()
            display(param_value)
            list_for_param.append(param_value.rhs)
            new_og_space = new_inverse_reliability_instance.new_design_point_in_original_space().rhs.subs(new_inverse_reliability_instance.estimate_value().lhs,new_inverse_reliability_instance.estimate_value().rhs)
            new_og_space = Eq(new_inverse_reliability_instance.new_design_point_in_original_space().lhs, new_og_space, evaluate=False)
            display(new_og_space)
            list_for_new_mdp.append(new_og_space.rhs)

            if abs(param_value.rhs-param_value_initial.rhs) <= treshold:
                display(f'Różnica dwóch ostatnich wartości parametrów jest mniejsza od progu równego {treshold}')
                break

            new_inverse_reliability_instance = new_inverse_reliability_instance.next_iteration()
            param_value_initial = param_value

            i=i+1

        return {'Original space':list_for_coords, 'Estimated parameter':list_for_param, 'New original space':list_for_new_mdp}

    def validate(self, n=None, treshold=None):
        '''
        Method validates the result for obtaining assumed reliability index.
        '''

        mean_values_dict = {}
        param_list = self.parameters_list

        for i,parameter in enumerate(param_list):

            if i < len(param_list):
                mean_values_dict[parameter] = param_list[i].subs(self.new_design_point_coords)
            else:
                break


        mean_values_dict[self.dependent_variable] = self.estimate_value().rhs



        validate_algorithm = Reliability(limit_state_function=self.limit_state_function,
                                         mean_values_dict=mean_values_dict,
                                         standard_deviations_dict=self.standard_deviations_dict,
                                         parameters_list=self.parameters_list,
                                         dependent_variable=self.dependent_variable)

        return validate_algorithm.search_for_stable_index(n=n, treshold=treshold)

class ReliabilityAlgorithm(InverseReliability):
    '''
    Class bases on First Order Reliability Method. By asumming expected reliability index and input, the class calculates new design point parameters.
    '''

    def __init__(self,
                limit_state_function,
                mean_values_dict,
                standard_deviations_dict,
                parameters_list,
                assumed_reliability_index,
                dependent_variable=None,
                new_design_point_coords=None,
                ):

        self.limit_state_function = limit_state_function
        self.mean_values_dict = mean_values_dict
        self.parameters_list = parameters_list
        self.standard_deviations_dict = standard_deviations_dict
        self.assumed_reliability_index = assumed_reliability_index

        if dependent_variable is not None: self.dependent_variable = dependent_variable
        else: self.dependent_variable = self.parameters_list[-1]

        if new_design_point_coords is not None:

            dictionary={}
            for i, param in enumerate(self.parameters_list):
                dictionary[param] = new_design_point_coords[i]

            self.new_design_point_coords = dictionary

        else:
            self.new_design_point_coords = mean_values_dict


    @classmethod
    def from_inverse_reliability_method(cls, InverseReliability):
        return cls(InverseReliability.limit_state_function, InverseReliability.mean_values_dict, InverseReliability.standard_deviations_dict, InverseReliability.parameters_list, InverseReliability.assumed_reliability_index, InverseReliability.dependent_variable, None)


    def toFORM(self):

        self.validate(n=1)
