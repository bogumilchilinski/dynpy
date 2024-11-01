from sympy import *
from dynpy.solvers.linear import ODESystem
from dynpy.solvers.nonlinear import MultiTimeScaleSolution

class LinearFirstOrder(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z = Function('z')(t)
        a = 10
        ode_eq=z.diff(t)+a*z-sin(10*t)

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=1)

        return odes

class LinearSecondOrder(ODESystem):

    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h= 3^2
        omg0= 14
        fs= 7
        fd= 5
        omg= 2
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z+fd*cos(omg*t)+fs

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes



class LinearHorizontalSystem(ODESystem):

    
    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        z= Function('z')(t)
        h= 3^2
        omg0= 14
        ode_eq=z.diff(t,2)+2*h*z.diff(t)+omg0**2*z

        odes = cls(ode_eq,dvars=z,ivar=t,ode_order=2)

        return odes



#gotowe!!!
class SpringMassSystem(ODESystem):

    
    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        c= 3
        k= 2
        m= 1
        ode_eq=Eq(m*x.diff(t,2)+c*x.diff(t)+k*x,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes
#przyklady FODE
class BasicFODE(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-2*x+1,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEexpExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,exp(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEexpResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,exp(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEsinExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,sin(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEsinResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,sin(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEcosExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x-1,cos(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEcosResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,cos(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEexpSumExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,exp(2*t)+exp(3*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes

class FODEsinSumExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,sin(t)+sin(3*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=1)

        return odes
#PRZYKLADY 2ND ORDER ODE
class BasicSODE(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)+x.diff(t)+x+1,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class BasicSODEWithoutXDot(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)+x+1,0)

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes
    
    
class SODEexpExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,exp(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEexpResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,exp(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEsinExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,sin(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEsinResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,sin(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEcosExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,cos(2*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEcosResonance(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,cos(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEexpSumExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,exp(2*t)+exp(3*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEsinSumExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,sin(t)+sin(3*t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class FODEsincosExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t)-x+1,sin(t)*cos(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes

class SODEsincosExcitation(ODESystem):


    @classmethod
    def from_reference_data(cls):
        t = Symbol('t')
        x= Function('x')(t)
        ode_eq=Eq(x.diff(t,2)-2*x.diff(t)+x+1,sin(t)*cos(t))

        odes = cls(ode_eq.lhs-ode_eq.rhs,dvars=x,ivar=t,ode_order=2)

        return odes