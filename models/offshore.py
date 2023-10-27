from .mechanics.principles import ComposedSystem
from sympy.physics.mechanics import dynamicsymbols
from sympy import (Symbol, symbols, Matrix, sin, cos, diff, sqrt, S, Eq, pi)
from ..dynamics import HarmonicOscillator

import pint
ureg = pint.UnitRegistry()


class DDoFVessel(ComposedSystem):

    scheme_name = 'vessel.jpg'

    def __init__(self,
                 m_vessel=Symbol('M_vessel', positive=True),
                 I_5=Symbol('I_5', positive=True),
                 qs=dynamicsymbols('H, Phi'),
                 wave_level=dynamicsymbols('W'),
                 wave_slope=dynamicsymbols('S'),
                 rho=Symbol('rho', positive=True),
                 g=Symbol('g', positive=True),
                 A_wl=Symbol('A_wl',positive=True),
                 V=Symbol('V', positive=True),
                 GM_L=Symbol('GM_L', positive=True),
                 CoB=Symbol('CoB', positive=True),
                 CoF=Symbol('CoF', positive=True),
                 A_h=Symbol('A_h',positive=True),
                 Phi_h=Symbol('Phi_h',positive=True),
                 omega=Symbol('omega',positive=True),
                 ivar=Symbol('t'),
                 **kwargs
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
        self.A_h = A_h
        self.Phi_h = Phi_h
        self.omega = omega
        

        # vessel mass and stiffness matrix
        M_matrix = Matrix([[m_vessel, 0], [0, I_5]])

        K_matrix = Matrix([[rho * g * A_wl, -rho * g * A_wl * (CoF - CoB)],
                           [-rho * g * A_wl * (CoF - CoB),
                            rho * g * V * GM_L]])

        # generalized displacements and velocities
        
#         wave_level = A_h*cos(omega*self.ivar + Phi_h)
#         wave_slope = diff(wave_level,self.ivar)*(omega/g)
        
        #dh_wave_slope = diff(wave_level,t,t)*(omega/g)   # TU JEST TEN MAGICZNY FORMUŁ, KTÓRY CHYBA JEDNAK POWINNIŚMY WPROWADZIĆ BO dq go tak nie policzy XD #
        
        q = Matrix(qs) + Matrix([wave_level, wave_slope])
        dq = q.diff(ivar)

        # lagrangian components definition
        self.kinetic_energy = S.Half * sum(dq.T * M_matrix * dq)
        self.potential_energy = S.Half * sum(
            Matrix(qs).T * K_matrix * Matrix(qs))

        super().__init__(Lagrangian=self.kinetic_energy -
                         self.potential_energy,
                         qs=qs,
                         ivar=ivar,**kwargs)

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m_vessel: r'mass of vessel \si{[\kilogram]},',
#             self.I_5:
#             r'moment of inertia of \num{5}-th degree (with respect to \(y\) axis, determined by the radius of gyration) \si{[\kilo\gram\metre\squared]},',
#             tuple(self.q): r'generalized coordinates,',
#             self.wave_level: r'???,',
#             self.wave_slope: r'???,',
#             self.rho: r'fluid density \si{[\kilo\gram/\cubic\metre]},',
#             self.g: r'acceleration of gravity \si{[\metre/\second\squared]},',
#             self.A_wl: r'wetted area \si{[\metre\squared]},',
#             self.V: r'submerged volume of the vessel \si{[\cubic\metre]},',
#             self.GM_L: r'longitudinal metacentric height \si{[\metre]},',
#             self.CoB: r'centre of buoyancy \si{[\metre]},',
#             self.CoF: r'centre of floatation \si{[\metre]},',
#             self.ivar: r'independent time variable.',
#         }

#         return self.sym_desc_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m_vessel: r'mass of vessel,',
            self.I_5:
            r'moment of inertia of 5-th degree (with respect to y axis, determined by the radius of gyration),',
            tuple(self.q): r'vessel generalized coordinates,',
            self.wave_level: r'wave level,',
            self.wave_slope: r'wave slope,',
            self.rho: r'fluid density,',
            self.g: r'acceleration of gravity,',
            self.A_wl: r'wetted area,',
            self.V: r'submerged volume of the vessel,',
            self.GM_L: r'longitudinal metacentric height,',
            self.CoB: r'centre of buoyancy,',
            self.CoF: r'centre of floatation.',

        }

        return self.sym_desc_dict
    
    def numerical_data(self):
        numbers_dict={
            self.I_5:21*self.m_vessel,
            self.m_vessel: 1e7,
#             self.wave_level: 0.5 * cos(2/7 * pi * self.ivar),
#             self.wave_slope: 0.5 * cos(2/7 * pi * self.ivar),
            self.rho:1025,
            self.g: 9.81,
            self.A_wl:4025,
            self.V:25000,
            self.GM_L:290,
            self.CoB:3.6,
            self.CoF:13,
            self.A_h:0.45,
            self.Phi_h:-2.95,
        }
        
        return numbers_dict    
    
    def units(self):
        units_dict={

            self.q[1]:ureg.radian,
            self.q[0]:ureg.meter,
            self.q[1].diff(self.ivar):ureg.radian/ureg.second,
            self.q[0].diff(self.ivar):ureg.meter/ureg.second,

           }
        
        return units_dict    
    
class TDoFCompensatedPayload(ComposedSystem):

    scheme_name = '3dofs_new.PNG'

    def __init__(self,
                 phi=dynamicsymbols('varphi'),
                 h=dynamicsymbols('h'),
                 h_c=dynamicsymbols('h_c'),
                 m_p=Symbol('m_p', positive=True),
                 k_w=Symbol('k_w', positive=True),
                 l_0=Symbol('l_0', positive=True),
                 qs=dynamicsymbols('varphi h h_c'),
                 y_e=dynamicsymbols('y_e'),
                 z_e=dynamicsymbols('z_e'),
                 m_c=Symbol('m_c', positive=True),
                 k_c=Symbol('k_c', positive=True),
                 l_c=Symbol('l_c', positive=True),
                 l_w=Symbol('l_w', positive=True),
                 g=Symbol('g', positive=True),
                 h_eq=Symbol('h_eq', positive=True),
                 h_ceq=Symbol('h_ceq', positive=True),
                 ivar=Symbol('t'),
                 **kwargs):

        self.m_p = m_p
        self.k_w = k_w
        self.l_0 = l_0
        self.y_e = y_e
        self.z_e = z_e
        self.m_c = m_c
        self.k_c = k_c
        self.l_c = l_c
        self.l_w = l_w
        self.g = g
        self.h_eq = h_eq
        self.h_ceq = h_ceq
        self.ivar = ivar
        self.phi = phi
        self.h = h
        self.h_c = h_c
        self.qs = qs

        phi, h, h_c = qs

        y = (h + h_eq + l_0 + l_c) * sin(phi) + y_e
        z = (h + h_eq + l_0 + l_c) * cos(phi) + z_e

        y_c = (h_c + h_ceq + l_0) * sin(phi) + y_e
        z_c = (h_c + h_ceq + l_0) * cos(phi) + z_e

        y, z, y_c, z_c = dynamicsymbols('y_{p},z_{p},y_c,z_c')

        v_c = sqrt(diff(y_c, ivar)**2 + diff(z_c, ivar)**2)
        v = sqrt(diff(y, ivar)**2 + diff(z, ivar)**2)

        positions_dict = {
            y: (h + h_eq + l_0 + l_c) * sin(phi) + y_e,
            z: (h + h_eq + l_0 + l_c) * cos(phi) + z_e,
            y_c: (h_c + h_ceq + l_0) * sin(phi) + y_e,
            z_c: (h_c + h_ceq + l_0) * cos(phi) + z_e,
        }

        self.kinetic_energy = S.Half * m_p * v**2 + S.Half * m_c * v_c**2

        self.potential_energy = (S.Half * k_w * (h_c + h_ceq)**2 +
                                 S.Half* k_c * (h + h_eq - (h_c + h_ceq))**2 -
                                 m_p * g * z - m_c * g * z_c)

        super(HarmonicOscillator,self).__init__((self.kinetic_energy - self.potential_energy).subs(positions_dict).doit(),
                         qs=qs,
                         ivar=ivar,**kwargs)


    def single_dof_transformation(self,p=None):
        
        if p is None:
            p = S.One
        
        
        omega, delta, eps = symbols('omega, delta, varepsilon',positive=True)
        A_y, A_z, Phi_y, Phi_z = symbols('A_y, A_z, Phi_y, Phi_z,',positive=True)
        
        sdof_payload = self(self.phi)
        
        nonlin_mathieu_model = sdof_payload.subs({self.h:0,self.h_c:0,self.m_c:0}).subs({self.h_eq:self.l_c,self.h_ceq:self.l_c,self.l_0:self.l_w-self.l_c})

        mathieu_nonlin_approx = nonlin_mathieu_model.approximated(n=3)

        mathieu_acc_form = mathieu_nonlin_approx._to_acc()

        mathieu_basis = mathieu_acc_form.subs({self.l_c:0,
                                               self.g:omega**2*self.l_w})

        mathieu_instability_eq = mathieu_basis.subs(S.Half*omega**2*self.phi**2,(S.Half)**3*omega**2*self.phi**2*(1+eps*delta),simultenouos=True)

        mathieu_transformed_eq = mathieu_instability_eq.subs({self.l_w:1/eps,
                                                              self.phi.diff(self.ivar)*self.phi**2*self.y_e.diff(self.ivar):0,
                                                              self.z_e:A_z*cos(omega/p*self.ivar+Phi_z),
                                                              self.y_e:A_y*cos(omega/p*self.ivar+Phi_y),
                                                              A_y*self.phi**2*eps:0
                                                             })

        return mathieu_transformed_eq._eoms
        
        

#     def symbols_description(self):
#         self.sym_desc_dict = {
#             self.m_p: r'mass of payload \si{[\kilogram]}',
#             self.k_w: r'wire stiffness \si{[\newton\per\meter]}',
#             self.l_0: r'length of the lifting cable \si{[\metre]}',
#             tuple(self.q): r'generalized coordinates',
#             self.y_e:
#             r'lateral displacement at crane tip obtained from RAOs (a regular wave excitation) \si{[\metre]}',
#             self.z_e:
#             r'vertical displacement at crane tip obtained from RAOs (a regular wave excitation) \si{[\metre]}',
#             self.m_c: r'mass of compensator \si{[\kilogram]}',
#             self.k_c:
#             r'stiffness of heave compensator \si{[\newton\per\meter]}',
#             self.l_c:
#             r'length of the attached compensating element \si{[\metre]}',
#             self.g: r'acceleration of gravity \si{[\metre/\second\squared]}',
#             self.h_eq: r'equilibrium point of payload \si{[\metre]}',
#             self.h_ceq: r'equilibrium point of compensator \si{[\metre]}',
#             self.ivar: r'independent time variable',
#         }

#         return self.sym_desc_dict

    def symbols_description(self):

        parent_symbols_dict=super().symbols_description()

        self.sym_desc_dict = {** parent_symbols_dict , **{
            self.m_p: r'mass of payload,',
            self.k_w: r'wire stiffness,',
            self.l_0: r'length of the lifting cable,',
            self.l_w: r'length of the lifting cable-payload system,',
            tuple(self.q): r'payload generalized coordinates,',
            self.y_e:
            r'lateral displacement at crane tip obtained from RAOs (a regular wave excitation),',
            self.z_e:
            r'vertical displacement at crane tip obtained from RAOs (a regular wave excitation),',
            self.m_c: r'mass of compensator,',
            self.k_c:
            r'stiffness of heave compensator,',
            self.l_c:
            r'length of the attached compensating element,',
            self.g: r'acceleration of gravity,',
            self.h_eq: r'equilibrium point of payload,',
            self.h_ceq: r'equilibrium point of compensator,'
        }}

        return self.sym_desc_dict

    def numerical_data(self):
        numbers_dict={
            self.h_eq: ((self.m_c+self.m_p)/self.k_w+(self.m_p)/self.k_c)*self.g,
            self.h_ceq: self.g*(self.m_c+self.m_p)/self.k_w,
            self.m_p: 1e5,
            self.k_w: 4.09e6,
            self.l_0: self.l_w - self.l_c,
            self.y_e: 0.5 * cos(2/7 * pi * self.ivar -2.396),
            self.z_e: 1.5 * cos(2/7 * pi * self.ivar + 0.315),
            self.m_c: 1e4,
            self.k_c: 2.5e6,
            self.l_c:5,
            self.l_w:48,
            self.g: 9.81,
        }
        
        return numbers_dict    
    

    
    def units(self):
        units_dict={
            self.q[2]:ureg.meter,
            self.q[1]:ureg.meter,
            self.q[0]:ureg.radian,
            self.q[2].diff(self.ivar):ureg.meter/ureg.second,            
            self.q[1].diff(self.ivar):ureg.meter/ureg.second,
            self.q[0].diff(self.ivar):ureg.radian/ureg.second,

           }
        
        return units_dict 
# class PayloadVesselSystem(ComposedSystem):
    
#     def __init__(self,
#                  y_e=dynamicsymbols('y_e'),
#                  z_e=dynamicsymbols('z_e'),
#                  wave_level=dynamicsymbols('W'),
#                  wave_slope=dynamicsymbols('S'),
# #                  payload=TDoFCompensatedPayload(),
# #                  vessel=DDoFVessel(),
#                  system=None
#                 ):
        
#         self.payload = TDoFCompensatedPayload(y_e,z_e)
#         self.vessel = DDofVessel(wave_level,wave_slope)
        
#         system = self.payload + self.vessel
        
#         super().__init__(system)


