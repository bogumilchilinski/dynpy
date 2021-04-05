from .systems import ComposedSystem
from sympy.physics.mechanics import dynamicsymbols
from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S)


class DDoFVessel(ComposedSystem):
    def __init__(self,
                 m_vessel=Symbol('M_vessel', positive=True),
                 I_5=Symbol('I_5', positive=True),
                 qs=dynamicsymbols('H, Phi'),
                 wave_level=dynamicsymbols('W'),
                 wave_slope=dynamicsymbols('S'),
                 rho=Symbol('rho', positive=True),
                 g=Symbol('g', positive=True),
                 A_wl=Symbol('A_wl'), positive=True,
                 V=Symbol('V', positive=True),
                 GM_L=Symbol('GM_L', positive=True),
                 CoB=Symbol('CoB', positive=True),
                 CoF=Symbol('CoF', positive=True),
                 ivar=Symbol('t')
                 ):

        # vessel mass and stiffness matrix
        M_matrix = Matrix([[m_vessel, 0], [0, I_5]])

        K_matrix = Matrix([[rho*g*A_wl, -rho*g*A_wl*(CoF - CoB)],
                          [-rho*g*A_wl*(CoF - CoB), rho*g*V*GM_L]])

        # generalized displacements and velocities
        q = Matrix(qs)+Matrix([wave_level, wave_slope])
        dq = q.diff(ivar)

        # lagrangian components definition
        self.T = 1/2 * sum(dq.T * M_matrix * dq)
        self.V = 1/2 * sum(Matrix(qs).T * K_matrix * Matrix(qs))
        
        self.sym_desc_dict={                 
                    m_vessel:'mass of vessel \si{[\kilogram]}',
                 I_5:'moment of inertia of \num{5}-th degree (with respect to \(y\) axis, determined by the radius of gyration) \si{[\kilo\gram\metre\squared]}',
                 qs:'generalized coordinates',
                 wave_level:'???',
                 wave_slope:'???',
                 rho:'fluid density \si{[\kilo\gram/\cubic\metre]}',
                 g:'acceleration of gravity \si{[\metre/\second\squared]}',
                 A_wl:'wetted area \si{[\metre\squared]}',
                 V:'submerged volume of the vessel \si{[\cubic\metre]}',
                 GM_L:'longitudinal metacentric height \si{[\metre]}',
                 CoB:'centre of buoyancy \si{[\metre]}',
                 CoF:'centre of floatation \si{[\metre]}',
                 ivar:'independent variable',
              }

        super().__init__(Lagrangian=self.T-self.V, qs=qs,ivar=ivar)


class TDoFCompensatedPayload(ComposedSystem):

    def __init__(self,
                 m_p=Symbol('m_p', positive=True),
                 k_w=Symbol('k_w', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 qs=dynamicsymbols('varphi h h_c'),
                 y_e=dynamicsymbols('y_e'),
                 z_e=dynamicsymbols('z_e'),
                 m_c=Symbol('m_c', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 l_c=Symbol('l_c', positive=True),
                 g=Symbol('g', positive=True),
                 h_eq=Symbol('h_eq', positive=True),
                 h_c_eq=Symbol('h_ceq', positive=True),
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

        super().__init__(self.T-self.V, qs=qs,ivar=ivar)

        
