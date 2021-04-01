from .systems import ComposedSystem
from sympy.physics.mechanics import dynamicsymbols
from sympy import (Symbol, symbols, Matrix,sin,cos)


class DDoFVessel(ComposedSystem):
    def __init__(self,
                 m_vessel, 
                 I_5,
                 qs=dynamicsymbols('H, Phi'),
                 wave_level=dynamicsymbols('W'),
                 wave_slope=dynamicsymbols('S'),
                 rho=Symbol('rho'),
                 g=Symbol('g'),
                 A_wl=Symbol('A_wl'),
                 V=Symbol('V'),
                 GM_L=Symbol('GM_L'),
                 CoB=Symbol('CoB'),
                 CoF=Symbol('CoF')
                 ):

        # vessel mass and stiffness matrix
        M_matrix = Matrix([[m_vessel, 0], [0, I_5]])

        K_matrix = Matrix([[rho*g*A_wl, -rho*g*A_wl*(CoF - CoB)],
                          [-rho*g*A_wl*(CoF - CoB), rho*g*V*GM_L]])
        
        #generalized displacements and velocities
        q=Matrix(qs)+Matrix([wave_level,wave_slope])
        dq=q.diff(self.ivar)
        
        # lagrangian components definition
        T = 1/2 * sum(dq.T * M_matrix * dq)
        V = 1/2 * sum(Matrix(qs).T * K_matrix *  Matrix(qs))


        super().__init__(Lagrangian=T-V,qs=qs)


class TDoFCompensatedPayload(ComposedSystem):
    pass

