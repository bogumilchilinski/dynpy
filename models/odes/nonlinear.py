from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.nonlinear import MultiTimeScaleSolution


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