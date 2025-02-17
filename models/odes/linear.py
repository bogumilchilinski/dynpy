from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.perturbational import MultiTimeScaleSolution
from dynpy.models.mechanics import ForcedSpringMassSystem
import numpy as np
import matplotlib.pyplot as plt

class LinearFirstOrder(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z = Function('z')(t)
        a=Symbol('a',positive=True)
        ode_eq=z.diff(t)+a*z-sin(10*t)

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=1)

        return odes
    
class LinearFirstOrderNonHomo(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z = Function('z')(t)
        a = Symbol('a',positive=True)
        ode_eq = z.diff(t)+a*z-cos(t)*Heaviside(t-10)

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=1)

        return odes

class LinearSecondOrder(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_static',positive=True)
        fd=Symbol('f_dynamic',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z+fd*cos(omg*t)+fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes

class LinearSecondOrderExpSin(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        a=Symbol('a',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z-fd*exp(-a*t)*sin(omg*t)-fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    
class LinearSecondOrderExpCos(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        a=Symbol('a',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z-fd*exp(-a*t)*cos(omg*t)-fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    
class LinearSecondOrderHeavisideComp(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+omg0**2*z-fd*Heaviside(sin(omg*t))-fs#+2*h*z.diff(t)

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    
class LinearSecondOrderSinComp(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_static',positive=True)
        fd=Symbol('f_dynamic',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z+fd*sin(omg*t)+fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z-fd*exp(-h*t)*sin(omg*t)-fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    
class LinearSecondOrderExpCos(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z-fd*exp(-h*t)*cos(omg*t)-fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    
class LinearSecondOrderHeavisideComp(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        fs=Symbol('f_s',positive=True)
        fd=Symbol('f_d',positive=True)
        omg=Symbol('Omega',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z-fd*Heaviside(omg*t)-fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes

class LinearHorizontalSystem(ODESystem):

    
    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes



#gotowe!!!
class SpringMassSystem(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        c=Symbol('c',positive=True)
        k=Symbol('k',positive=True)
        m=Symbol('m',positive=True)
        F=Symbol('F',positive=True)
        Omega=Symbol('Omega', positive=True)
        ode_eq=Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,F*sin(Omega*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SpringMassSmallParameter(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        eps=Symbol('\epsilon')
        omega=Symbol('\omega',positive=True)
        ode_eq=Eq(x.diff(t,2) + eps*diff(x) + omega**2*x,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SpringMassDamperMSM(MultiTimeScaleSolution):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        eps=Symbol('\epsilon')
        omega=Symbol('\omega',positive=True)
        ode_eq=Eq(x.diff(t,2) + eps*diff(x) + omega**2*x,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x, eps=eps)

        return odes

class SpringMassDamperExcitedMSM(MultiTimeScaleSolution):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        eps=Symbol('\epsilon')
        omega=Symbol('\omega',positive=True)
        omega_0 = Symbol('omega_0',positive = True)
        Omega = Symbol('Omega', positive = True)
        A = Symbol('A', positive = True)
        delta = Symbol('delta', positive = True)
        k=Symbol('k',positive=True)
        m=Symbol('m',positive=True)
        mu=Symbol('mu',positive=True)
        
        MSM_eq=Eq(x.diff(t,2) + eps*diff(x) + omega_0**2*x,eps*mu*sin(Omega * t))
        MSM = cls(MSM_eq.lhs-MSM_eq.rhs,dvars=x, eps=eps)
#         exc_dict = {omega_0:Omega**2+eps*delta}
#         MSM_sol_sub = MSM.solution.subs(exc_dict)
        
        return MSM
    
class SpringMassExcitedMSM(MultiTimeScaleSolution):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        eps=Symbol('\epsilon')
        omega=Symbol('\omega',positive=True)
        omega_0 = Symbol('omega_0',positive = True)
        Omega = Symbol('Omega', positive = True)
        A = Symbol('A', positive = True)
        delta = Symbol('delta', positive = True)
        k=Symbol('k',positive=True)
        m=Symbol('m',positive=True)
        
        MSM_eq=Eq(x.diff(t,2) + omega_0**2*x,eps*sin(Omega * t))
        MSM = cls(MSM_eq.lhs-MSM_eq.rhs,dvars=x, eps=eps)
#         exc_dict = {omega_0:Omega**2+eps*delta}
#         MSM_sol_sub = MSM.solution.subs(exc_dict)
        
        return MSM


class HarmonicOscillator(ODESystem):


    @classmethod
    def from_reference_data(cls):
        system=dyn_sys=ForcedSpringMassSystem()
        sim_time=250
        t_span = np.linspace(1,sim_time,sim_time+1)
        ic_list = [0.0,0.0]
        sym_num=18
        sym_list=np.linspace(0,int((sym_num*2)-2),int(sym_num))
        k_list=np.linspace(0.1,1.5,sym_num)
        omg=Symbol("omega")
        param_dict={system.F:1*cos(omg*system.ivar),system.k:omg,system.g:0,system.m:1}
        na_df = dyn_sys.subs(param_dict).eoms.numerical_analysis(parameter=omg,param_span=k_list, t_span = t_span)
        na_df_sol = na_df.compute_solution(t_span = t_span, ic_list = [0.1,0.1]).iloc[:,sym_list]
        plot=na_df_sol.plot()
        plt.show()
        max_plot=na_df_sol.set_axis(k_list,axis=1).max().plot()
        plt.show()

        return plot, max_plot

class ParametricOscillator(ODESystem):


    @classmethod
    def from_reference_data(cls):
        system=dyn_sys=ForcedSpringMassSystem()
        sim_time=250
        t_span = np.linspace(1,sim_time,sim_time+1)
        ic_list = [0.0,0.0]
        sym_num=18
        sym_list=np.linspace(0,int((sym_num*2)-2),int(sym_num))
        k_list=np.linspace(0.1,1.5,sym_num)
        omg=Symbol("omega")
        param_dict_k={system.k:1+0.01*cos(omg*system.ivar),system.F:1*cos(omg*system.ivar),system.g:0,system.m:1}
        na_df_k = dyn_sys.subs(param_dict_k).eoms.numerical_analysis(parameter=omg,param_span=k_list, t_span = t_span)
        na_df_sol_k = na_df_k.compute_solution(t_span = t_span, ic_list = [0.1,0.1]).iloc[:,sym_list]
        plot=na_df_sol_k.plot()
        plt.show()
        max_plot=na_df_sol_k.set_axis(k_list,axis=1).max().plot()
        plt.show()

        return plot, max_plot


