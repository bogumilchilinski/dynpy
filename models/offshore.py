from .systems import ComposedSystem
from sympy.physics.mechanics import dynamicsymbols
from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S)


class DDoFVessel(ComposedSystem):
    def __init__(self,
                 m_vessel=Symbol('M_vessel'),
                 I_5=Symbol('I_5'),
                 qs=dynamicsymbols('H, Phi'),
                 wave_level=dynamicsymbols('W'),
                 wave_slope=dynamicsymbols('S'),
                 rho=Symbol('rho'),
                 g=Symbol('g'),
                 A_wl=Symbol('A_wl'),
                 V=Symbol('V'),
                 GM_L=Symbol('GM_L'),
                 CoB=Symbol('CoB'),
                 CoF=Symbol('CoF'),
                 ivar=Symbol('t')
                 ):

        # vessel mass and stiffness matrix
        M_matrix = Matrix([[m_vessel, 0], [0, I_5]])

        K_matrix = Matrix([[rho*g*A_wl, -rho*g*A_wl*(CoF - CoB)],
                          [-rho*g*A_wl*(CoF - CoB), rho*g*V*GM_L]])

        # generalized displacements and velocities
        q = Matrix(qs)+Matrix([wave_level, wave_slope])
        dq = q.diff(self.ivar)

        # lagrangian components definition
        self.T = 1/2 * sum(dq.T * M_matrix * dq)
        self.V = 1/2 * sum(Matrix(qs).T * K_matrix * Matrix(qs))

        super().__init__(Lagrangian=self.T-self.V, qs=qs)


class TDoFCompensatedPayload(ComposedSystem):

    def __init__(self,
                 m_p=Symbol('m_p'),
                 k_w=Symbol('k_w'),
                 l_0=Symbol('l_0'),
                 qs=dynamicsymbols('varphi h h_c'),
                 y_e=dynamicsymbols('y_e'),
                 z_e=dynamicsymbols('z_e'),
                 m_c=Symbol('m_c'),
                 k_c=Symbol('k_c'),
                 l_c=Symbol('l_c'),
                 g=Symbol('g', positive=True),
                 h_eq=Symbol('h_eq'),
                 h_c_eq=Symbol('h_ceq'),
                 ivar=Symbol('t')
                 ):

        phi, h, h_c = qs

        y = (h+h_eq+l_0+l_c)*sin(phi)+y_e
        z = (h+h_eq+l_0+l_c)*cos(phi)+z_e
        v = sqrt(diff(y, ivar)**2+diff(z, ivar)**2)

        y_c = (h_c+h_c_eq+l_0)*sin(phi)+y_e
        z_c = (h_c+h_c_eq+l_0)*cos(phi)+z_e
        v_c = sqrt(diff(y_c, ivar)**2+diff(z_c, ivar)**2)

        self.T = S.One/2*m_p*v**2 + S.One/2*m_c*v_c**2

        self.V = (S.One/2*k_w*(h_c+h_c_eq)**2 +
                  1/2*k_c * (h+h_eq-(h_c+h_c_eq))**2
                  - m_p*g*z - m_c*g*z_c)

        super().__init__(self.T-self.V, qs=qs)
