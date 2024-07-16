from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset, atan)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin
from ...continuous import ContinuousSystem, PlaneStressProblem

import base64
import random
import IPython as IP
import numpy as np
import inspect
import pandas as pd

ix = pd.IndexSlice

from .pendulum import Pendulum, PendulumKinematicExct
from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST
from .tmd import TunedMassDamperRelativeMotion
from ...utilities.components.mech import en as mech_comp
from .trolley import VariableMassTrolleyWithPendulumRayleighDamping

from pylatex.utils import italic, NoEscape

from pint import UnitRegistry
ureg = UnitRegistry()

from functools import cached_property, lru_cache

class AdaptableTunedFluidDamper(VariableMassTrolleyWithPendulumRayleighDamping):
    
    trolley_max_value = Symbol('tv_max')
    trolley_min_value = Symbol('tv_min')
    trolley_max_spec_value = Symbol('tsv_max')
    pendulum_max_value = Symbol('pv_max')
    pendulum_min_value = Symbol('pv_min')
    pendulum_max_spec_value = Symbol('psv_max')
    mainplot_tikz_data_graph = Symbol('mx_dg')
    mainplot_tikz_value_graph = Symbol('mx_vg')
    spectrum_tikz_data_graph = Symbol('sx_dg')
    spectrum_tikz_graph = Symbol('sx_vg')
    mainplot_phi_tikz_data_graph = Symbol('mp_dg')
    mainplot_phi_tikz_graph = Symbol('mp_vg')
    spectrum_phi_tikz_data_graph = Symbol('sp_dg')
    spectrum_phi_tikz_graph = Symbol('sp_vg')

    def simulation_performance(self, params_dict, parameter, param_span, t_span, ic_list):
        
        na_df = self._ode_system.subs(params_dict).numerical_analysis(parameter=parameter, param_span=param_span)
                
        simulation_result = na_df.compute_solution(t_span = t_span, ic_list = ic_list)
        
        return simulation_result.droplevel(0, axis=1)
    
    def simulation_results_presentation(self, simulation_result, dynamic_variable, caption = None):
        
        trolley_max_value = Symbol('tv_max')
        trolley_min_value = Symbol('tv_min')
        trolley_max_spec_value = Symbol('tsv_max')
        pendulum_max_value = Symbol('pv_max')
        pendulum_min_value = Symbol('pv_min')
        pendulum_max_spec_value = Symbol('psv_max')
        mainplot_tikz_data_graph = Symbol('mx_dg')
        mainplot_tikz_value_graph = Symbol('mx_vg')
        spectrum_tikz_data_graph = Symbol('sx_dg')
        spectrum_tikz_graph = Symbol('sx_vg')
        mainplot_phi_tikz_data_graph = Symbol('mp_dg')
        mainplot_phi_tikz_graph = Symbol('mp_vg')
        spectrum_phi_tikz_data_graph = Symbol('sp_dg')
        spectrum_phi_tikz_graph = Symbol('sp_vg')

        display(dynamic_variable[0])
        
        trolley_max = simulation_result.loc[:, ix[:, dynamic_variable[0]]].max().values[0].round(2)
        trolley_min = simulation_result.loc[:, ix[:, dynamic_variable[0]]].min().values[0].round(2)
        trolley_max_spec = simulation_result.loc[:, ix[:, dynamic_variable[0]]].to_frequency_domain().single_sided_rms().truncate(0, 0.6).max().values[0].round(2)

        pendulum_max = simulation_result.loc[:, ix[:, dynamic_variable[1]]].max().values[0].round(2)
        pendulum_min = simulation_result.loc[:, ix[:, dynamic_variable[1]]].min().values[0].round(2)
        pendulum_max_spec = simulation_result.loc[:, ix[:, dynamic_variable[1]]].to_frequency_domain().single_sided_rms().truncate(0, 0.6).max().values[0].round(2)

        mainplot_tikz_data = simulation_result.loc[:,ix[:,dynamic_variable[0]]][list((simulation_result.loc[:,ix[:,dynamic_variable[0]]].columns))]
        
        mainplot_tikz = mainplot_tikz_data.to_pylatex_tikz(height=NoEscape(r'8cm'),width=NoEscape(r'16cm')).in_figure(caption = f"Trolley displacement for various amount of {caption}")

        spectrum_tikz_data = (simulation_result.loc[:,ix[:,dynamic_variable[0]]].to_frequency_domain().single_sided_rms( ).truncate(0, 1))[list(((simulation_result.loc[:,ix[:,dynamic_variable[0]]].to_frequency_domain().single_sided_rms( ).truncate(0, 1)).columns))]
        
        spectrum_tikz = spectrum_tikz_data.to_pylatex_tikz(height=NoEscape(r'8cm'),width=NoEscape(r'16cm')).in_figure(caption = f"Frequency spectra of trolley displacement for various amount of {caption}")

        mainplot_phi_tikz_data = simulation_result.loc[:,ix[:,dynamic_variable[1]]][list((simulation_result.loc[:,ix[:,dynamic_variable[1]]].columns))]
        
        mainplot_phi_tikz = mainplot_phi_tikz_data.to_pylatex_tikz(height=NoEscape(r'8cm'),width=NoEscape(r'16cm')).in_figure(caption = f"Pendulum angular displacement for various amount of {caption}")

        spectrum_phi_tikz_data = (simulation_result.loc[:,ix[:,dynamic_variable[1]]].to_frequency_domain().single_sided_rms( ).truncate(0, 1))[list(((simulation_result.loc[:,ix[:,dynamic_variable[1]]].to_frequency_domain().single_sided_rms( ).truncate(0, 1)).columns))]
        
        spectrum_phi_tikz = spectrum_phi_tikz_data.to_pylatex_tikz(height=NoEscape(r'8cm'),width=NoEscape(r'16cm')).in_figure(caption = f"Frequency spectra of pendulum angular displacement for various amount of {caption}")
        
        results_presentation = {
            trolley_max_value: trolley_max,
            trolley_min_value: trolley_min,
            trolley_max_spec_value: trolley_max_spec,
            pendulum_max_value: pendulum_max,
            pendulum_min_value: pendulum_min,
            pendulum_max_spec_value: pendulum_max_spec,
            mainplot_tikz_data_graph: mainplot_tikz_data,
            mainplot_tikz_value_graph: mainplot_tikz,
            spectrum_tikz_data_graph: spectrum_tikz_data,
            spectrum_tikz_graph: spectrum_tikz,
            mainplot_phi_tikz_data_graph: mainplot_phi_tikz_data,
            mainplot_phi_tikz_graph: mainplot_phi_tikz,
            spectrum_phi_tikz_data_graph: spectrum_phi_tikz_data,
            spectrum_phi_tikz_graph: spectrum_phi_tikz,
        }
        
        return results_presentation

    def investigation_of_initial_frequency_influence(self, ic_list=None, t_span=None, params_list=None, params_dict=None):

        if ic_list is None: ic_list =[0.0, 0.0, 0.0, 0.0]
        if t_span is None: t_span=np.linspace(0,100,2001)
        if params_list is None: params_list = [0.494,0.594,0.694,0.794,0.894]
        if params_dict is None:
            params_dict= {
                self.m_p: 10,
                self.m_f: 70,
                self.b: 2,
                self.F: 50,
                self.flow_coeff: 10,
                self.t_0: 30,
                self.g: 9.81,
                self.l: 6,
                self.k: 80,
                self.m_t: 80,
                self.m_tf:70,
                self.m_pf:0
            }

        parameter=self.Omega
        
        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)

    def adaptable_tuned_mass_damper_validation(self):

        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.b: 2,
            self.F: 50,
            self.flow_coeff: 10,
            self.Omega: 0.794,
            self.t_0: 30,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.t_0
        params_dict[self.t_0] = self.t_0
        
        params_dict_2 = params_dict
        params_dict_2[self.t_0] = self.t_0
        
        params_dict_2.update({self.m_t:  8.0*params_dict_2[self.m_p],
                            self.m_tf: params_dict_2[self.m_f] + 0.000000001*params_dict_2[self.t_0],
                            self.m_pf: 0.000000001*params_dict_2[self.t_0]})

        params_dict.update({self.m_t:  8.0*params_dict[self.m_p],
                            self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })
        
        param_min=110
        param_max=120
        cases_no = 1

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)
        
        ic_list = [0.0,0.0,0.0,0.0]
        
        simple_model_na_df = self.simulation_performance(params_dict = params_dict_2, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)
        
        advanced_model_na_df = self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)
        
        return pd.concat([simple_model_na_df, advanced_model_na_df], axis=1)

    
    def investigation_of_damping_coefficient_influence(self):
        
        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.b: 15,
            self.F: 50,
            self.flow_coeff: 10,
            self.t_0: 30,
            self.Omega: 0.694,#3.14,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.b

        params_dict.update({self.m_t:  8.0*params_dict[self.m_p],
                            self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })
        
        params_dict[parameter] = parameter

        param_min=4
        param_max=20
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)
        
        ic_list = [0.0,0.0,0.0,0.0]
        
        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)

    def investigation_of_damping_coefficient_damping_capabilities(self):
        
        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.b: 10,
            self.F: 50,
            self.flow_coeff: 10,
            self.t_0: 110,
            self.Omega: 0.694,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.b

        params_dict.update({self.m_t:  8.0*params_dict[self.m_p],
                            self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })
        
        params_dict[parameter] = parameter

        param_min=4
        param_max=20
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)
        
        ic_list = [0.0,0.0,0.0,0.0]
        
        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)
    
    def investigation_of_activation_time_influence(self):
        
        
        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.b: 10,
            self.F: 50,
            self.flow_coeff: 10,
            self.Omega: 0.694,
            self.t_0: 110,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.t_0
        
        params_dict[self.t_0] = self.t_0

        params_dict.update({self.m_t:  8.0*params_dict[self.m_p],
                            self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })
        
        param_min=-30
        param_max=90
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)
        
        ic_list = [0.0,0.0,0.0,0.0]
        
        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)

    
    def investigation_of_trolley_mass_influence(self):

        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.b: 10,
            self.F: 50,
            self.flow_coeff: 10,
            self.Omega: 0.694,
            self.t_0: 30,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.m_t

        params_dict[parameter] = parameter

        params_dict.update({self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })

        param_min=60
        param_max=100
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)

        ic_list = [0.0,0.0,0.0,0.0]

        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)

    def investigation_of_pendulum_mass_influence(self):

        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.m_t: 80,
            self.b: 10,
            self.F: 50,
            self.flow_coeff: 10,
            self.Omega: 0.694,
            self.t_0: 30,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.m_p

        params_dict[parameter] = parameter

        params_dict.update({self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })

        param_min=5
        param_max=25
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)

        ic_list = [0.0,0.0,0.0,0.0]

        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)
    
    def investigation_of_fluid_mass_influence(self):

        params_dict= {
            self.m_p: 10,
            self.m_f: 70,
            self.m_t: 80,
            self.b: 10,
            self.F: 50,
            self.flow_coeff: 10,
            self.Omega: 0.694,
            self.t_0: 30,
            self.g: 9.81,
            self.l: 6,
            self.k: 80,
        }

        t_span=np.linspace(0,100,2001)

        parameter=self.m_f

        params_dict[parameter] = parameter

        params_dict.update({self.m_tf: params_dict[self.m_f]*((S.One/2-atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi)),
                            self.m_pf: params_dict[self.m_f]*((S.One/2+atan(params_dict[self.flow_coeff]*(self.ivar-params_dict[self.t_0]))/pi))
                           })

        param_min=50
        param_max=90
        cases_no = 5

        params_list = []
        for value in np.linspace(param_min, param_max, cases_no).round(3):
            params_list.append(value)

        ic_list = [0.0,0.0,0.0,0.0]

        return self.simulation_performance(params_dict = params_dict, parameter = parameter, param_span = params_list, t_span = t_span, ic_list = ic_list)