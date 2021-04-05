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

        self.m_vessel = m_vessel
        self.I_5 = I_5
        self.wave_level = wave_level
        self.wave_slope = wave_slope
        self.rho = rho
        self.g = g
        self.A_wl = A_wl
        self.V = V
        self.GM_L = GM_L
        self.CoB = CoB
        self.CoF = CoF

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

        super().__init__(Lagrangian=self.T-self.V, qs=qs, ivar=ivar)

    def symbols_description(self):
        self.sym_desc_dict = {self.m_vessel: r'mass of vessel \si{[\kilogram]}',
                              self.I_5: r'moment of inertia of \num{5}-th degree (with respect to \(y\) axis, determined by the radius of gyration) \si{[\kilo\gram\metre\squared]}',
                              self.q: r'generalized coordinates',
                              self.wave_level: r'???',
                              self.wave_slope: r'???',
                              self.rho: r'fluid density \si{[\kilo\gram/\cubic\metre]}',
                              self.g: r'acceleration of gravity \si{[\metre/\second\squared]}',
                              self.A_wl: r'wetted area \si{[\metre\squared]}',
                              self.V: r'submerged volume of the vessel \si{[\cubic\metre]}',
                              self.GM_L: r'longitudinal metacentric height \si{[\metre]}',
                              self.CoB: r'centre of buoyancy \si{[\metre]}',
                              self.CoF: r'centre of floatation \si{[\metre]}',
                              self.ivar: r'independent time variable',
                              }

        return self.sym_desc_dict


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

        self.m_p = m_p
        self.k_w = k_w
        self.l_0 = l_0
        self.qs = qs
        self.y_e = y_e
        self.z_e = z_e
        self.m_c = m_c
        self.k_c = k_c
        self.l_c = l_c
        self.g = g
        self.h_eq = h_eq
        self.h_c_eq = h_c_eq
        self.ivar = ivar

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

        super().__init__(self.T-self.V, qs=qs, ivar=ivar)

    def symbols_description(self):
        self.sym_desc_dict = {self.m_p: '???',
                              self.k_w: '???',
                              self.l_0: '???',
                              self.q:   'generalized coordinates',
                              self.y_e: '???',
                              self.z_e: '???',
                              self.m_c: '???',
                              self.k_c: '???',
                              self.l_c: '???',
                              self.g:   '???',
                              self.h_eq:   '???',
                              self.h_c_eq: '???',
                              self.ivar:   'independent time variable',
                              }

        return self.sym_desc_dict
