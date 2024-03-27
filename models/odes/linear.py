from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.nonlinear import MultiTimeScaleSolution

class LinearFirstOrder(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z = Function('z')(t)
        a=Symbol('a',positive=True)
        ode_eq=z.diff(t)+a*z-sin(10*t)

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
        ode_eq=Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes
