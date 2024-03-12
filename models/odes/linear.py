from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.nonlinear import MultiTimeScaleSolution

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

        odes = cls(ode_eq,z,t,ode_order=2)

        return odes

class LinearSecondOrderMSM(ODESystem):
    
    @classmethod
    def from_reference_data(cls):
        F=Symbol('F')
        c=Symbol('c',positive=True)
        k=Symbol('k',positive=True)
        eps=Symbol('varepsilon')
        t=Symbol('t')
        z=Function('z')(t)
        a0=Symbol('a_0',positive=True)
        a1=Symbol('a_1',positive=True)
        omg=Symbol('omega',positive=True)
        m=Symbol('m',positive=True)
        eom=-F+c*z.diff(t)+k*(eps*(a0*sin(omg*t)+a1*sin(2*omg*t)+1)*z+m*z.diff(t,2))


        odes = cls(eom,dvars=z,ivar=t, order=2)._as_msm()
        
        return odes



class LinearEzSystem(ODESystem):

    
    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h=Symbol('h',positive=True)
        omg0=Symbol('omega_0',positive=True)
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes
    

#klasa rozwijana 

class SpringMassSystem(ODESystem):

    
    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        c=Symbol('c',positive=True)
        k=Symbol('k',positive=True)
        m=Symbol('m',positive=True)
        ode_eq=Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes
