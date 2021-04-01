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
        self.T = 1/2 * sum(dq.T * M_matrix * dq)
        self.V = 1/2 * sum(Matrix(qs).T * K_matrix *  Matrix(qs))


        super().__init__(Lagrangian=self.T-self.V,qs=qs)


class TDoFCompensatedPayload(ComposedSystem):

    def __init__(self,qs):
        self.T = 1/2*m_p*v**2 + 1/2*m_c*v_c**2
        self.V = 1/2*k_w*(h_c+h_c_eq)**2 + 1/2*k_c*(h+h_eq-(h_c+h_c_eq))**2 - m_p*g*z - m_c*g*z_c 
        super().__init__(self.T-self.V,qs=qs)

