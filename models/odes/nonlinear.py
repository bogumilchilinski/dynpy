from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.nonlinear import MultiTimeScaleSolution
from dynpy.models.mechanics.trolleys import MDoFTMD1

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
    
class MDoFTMDMSM(MultiTimeScaleSolution):
    
    @classmethod
    def from_reference_data(cls):

        t=Symbol('t')
        xb=Function("x_b")(t)
        xe=Function("x_e")(t)
        eps= Symbol('varepsilon',positive=True)
        delta =  Symbol('delta',positive= True)
        delta1 =  Symbol('delta_1',positive= True)
        delta2 =  Symbol('delta_2',positive= True)
        omg1 =  Symbol('omega_1',positive= True)
        omg2 =  Symbol('omega_2',positive= True)

        omg =  Symbol('Omega',positive= True)

        K_mat = Matrix([[omg**2,0],[0,omg**2]])

        stiff_part = K_mat * Matrix([xb,xe])

        odes = Matrix([xb.diff(t,2),xe.diff(t,2)]) + stiff_part +eps*Matrix([xb*(omg1**2-omg**2)+xb,xe*(omg2**2-omg**2)])+eps*Matrix([1*xe*cos(omg/10*t)+cos(omg*t),0*xe])


        odesy = cls(odes,Matrix([xb,xe]),order=2,eps=eps)#._as_msm()
        
        return odesy
class WeakNonLinearOscillatorMSM(MultiTimeScaleSolution):
    
    @classmethod
    def from_reference_data(cls):
        tau = Symbol('t')
        omega,p = symbols('omega p',positive=True)
        delta,eps, A = symbols('delta varepsilon A', positive=True)
        x = Function('x')(tau)

        p=1

        ode = ODESystem(odes=Matrix([x.diff(tau,2) + ((S.One/p)**2*(1)+delta*eps)*x +A*eps*cos(tau)+(eps)*x.diff(tau)]), dvars=Matrix([x]), ivar=tau, ode_order=2)


        odes = cls(ode.odes[0], ode.dvars[0], ivar=ode.ivar, omega=S.One, order=2,eps=eps)
        
        return odes
class CubicWeakOscillatorMSM(MultiTimeScaleSolution):
    
    @classmethod
    def from_reference_data(cls):
        tau = Symbol('t')
        omega,p = symbols('omega p',positive=True)
        delta,eps, A = symbols('delta varepsilon A', positive=True)
        x = Function('x')(tau)

        p=1

        ode = ODESystem(odes=Matrix([x.diff(tau,2) + ((S.One/p)**2*(1)+delta*eps)*x +A*eps*cos(tau)+(eps)*x.diff(tau)+x*eps**3]), dvars=Matrix([x]), ivar=tau, ode_order=2)


        odes = cls(ode.odes[0], ode.dvars[0], ivar=ode.ivar, omega=S.One, order=2,eps=eps)
        
        return odes