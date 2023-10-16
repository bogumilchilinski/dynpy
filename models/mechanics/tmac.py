from sympy import (Symbol, symbols, Matrix, sin, cos, asin, diff, sqrt, S,
                   diag, Eq, hessian, Function, flatten, Tuple, im, pi, latex,
                   dsolve, solve, fraction, factorial, Subs, Number, oo, Abs,
                   N, solveset)

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point
from sympy.physics.vector import vpprint, vlatex
from ...dynamics import LagrangesDynamicSystem, HarmonicOscillator, mech_comp

from ..elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame, base_origin, EngineMount
from ..continuous import ContinuousSystem, PlaneStressProblem
from dynpy.models.mechanics.tmd import TunedMassDamper

import base64
import random
import IPython as IP
import numpy as np
import inspect

from .principles import ComposedSystem, NonlinearComposedSystem, base_frame, base_origin, REPORT_COMPONENTS_LIST


class CrankSystem(ComposedSystem):

    scheme_name = 'crank_mechanismlow.jpg'
    real_name = 'crank_slider_real.jpg'

    
    _default_folder_path = "./dynpy/models/images/"



    def _detail_scheme(self):

        self.preview()

        return self._path
    
    
    def preview(self, example=False):
        

        if self._path:
            path = self._path
            
        else:

            plt.figure(figsize=(10,10))


            plt.xlim(-0.5,1)
            plt.ylim(-0.25,1.25)
            plt.grid(True)
            
            if (self._given_data)=={}:
            
                path = self.__class__._default_folder_path + self.__class__.scheme_name
            else:
                

                
                h_num=float(self.h.subs(self._given_data))
                r_num=float(self.r.subs(self._given_data))
                phi_num=float(self.phi.subs(self._given_data))
                a_num=float(self.a.subs(self._given_data))
                b_num=float(self.b.subs(self._given_data))                
                
                pOx = 0
                pOy = 0
                plt.text(-0.05+(pOx),(-0.05+pOy),'O')
                
                pAx=0
                pAy=h_num
                plt.text(-0.05+(pAx),(+0.05+pAy),'A')
                
                plt.plot([pOx,pAx],[pOy,pAy],'b',linewidth=2)
                plt.text(-0.05+(pOx+pAx)/2,(pOy+pAy)/2,'h')
                
                pBx = pAx + r_num*np.sin(phi_num)
                pBy = pAy + r_num*np.cos(phi_num)
                plt.text(+0.05+(pBx),(+0.00+pBy),'B')
                
                plt.plot([pBx,pAx],[pBy,pAy],'r',linewidth=2)
                plt.text((pBx+pAx)/2,0.05+(pBy+pAy)/2,'r')
                
                lBO = ((pBx-pOx)**2 + (pBy-pOy)**2)**0.5
                
                p1x = 1.1*(h_num+r_num)*pBx/lBO
                p1y = 1.1*(h_num+r_num)*pBy/lBO
                
                plt.plot([p1x,pOx],[p1y,pOy],'g',linewidth=2)

                pCx = a_num*pBx/lBO
                pCy = a_num*pBy/lBO
                plt.text(+0.05+(pCx),(+0.05+pCy),'C')
                
                plt.plot([pOx,pCx],[pOy,pCy],'m',linewidth=2)
                plt.text(0.05+(pOx+pCx)/2,(pOy+pCy)/2,'a')
                
#                 pDx = pCx+(b_num**2 - pCy**2)**0.5        # poprawna wartość położenia D
                pDx = (pCx+(b_num**2 - pCy**2)**0.5)*1.15   # celowo wprowadzony błąd
                pDy = 0
                plt.text(+0.05+(pDx),(+0.00+pDy),'D')
                
                plt.plot([pCx,pDx],[pCy,pDy],'k',linewidth=2)
                plt.text((pCx+pDx)/2,0.05+(pCy+pDy)/2,'b')
            
                path = self.__class__._default_folder_path + 'previews/' + self.__class__.__name__ + str(
                                            next(self.__class__._case_no)) + '.png'

                plt.savefig(path)
                self._path=path

                plt.close()

        

#        print('check' * 100)
        print(self._path)
#        print('check' * 100)
        plt.close()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))
    
    
    
    def __init__(self,
                 I=Symbol('I', positive=True),
                 r=Symbol('r', positive=True),
                 h=Symbol('h', positive=True),
                 a=Symbol('a', positive=True),
                 b=Symbol('b', positive=True),
                 phi=dynamicsymbols('\\varphi'),
                 beta=dynamicsymbols('beta'),
#                  alpha=dynamicsymbols('alpha'),
                 **kwargs):

        self.I = I
        self.h=h
        self.r=r
        self.phi = phi
        self.beta = beta
#         self.alpha = alp
        self.a = a
        self.b = b

        self.crank = MaterialPoint(I, phi, qs=[phi])
        composed_system = (self.crank)

        super().__init__(composed_system,**kwargs)

        
    def _crank_horizontal_pos(self):
        return None
        
    @property
    def _dbeta(self):
        beta=atan(self.r/self.l*self.phi) #it's probably wrong - has to be checked
        return beta.diff(self.ivar)
    
    @property
    def _displacement_d(self):

        return -(-self.b*sqrt((self.h**2 - self.h**2*self.a**2/self.b**2 + 2*self.h*self.r*cos(self.phi) - 2*self.h*self.r*self.a**2*cos(self.phi)/self.b**2 + self.r**2 - self.r**2*self.a**2*cos(self.phi)**2/self.b**2)/(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2)) - self.r*self.a*sin(self.phi)/sqrt(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2))
    
    @property
    def _velocity_b2(self):
        return (self.phi.diff(self.ivar)*self.r)
    @property
    def _velocity_b2b3(self):
        return (sqrt(self.h**2 + 2*self.h*self.r*cos(self.phi) + self.r**2)).diff(self.ivar)
    @property
    def _velocity_b3(self):
        gamma=asin(self._velocity_b2b3/self._velocity_b2)
        return self._velocity_b2*cos(gamma)
    @property
    def _velocity_d(self):
        return self._displacement_d.diff(self.ivar)
    @property
    def _velocity_c(self):
        omega3=self._velocity_b3/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi))
        return omega3*self.a
    @property
    def _acceleration_b2(self):
        return self.phi.diff(self.ivar)**2*self.r
    @property
    def _acceleration_b3n(self):
        return (self._velocity_b3**2/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))
    @property
    def _acceleration_cn(self):
        return (self._velocity_b3**2*(self.a/(sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))**2))
    @property
    def _acceleration_d(self):
        return self._velocity_d.diff(self.ivar)

    
    @property
    def _omega_3(self):
        return (self._velocity_b3/sqrt(self.r**2 + self.h**2 - 2*self.r*self.h*cos(pi-self.phi)))
    @property
    def linkage_ang_velocity(self):
        return self._dbeta

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I, l, m0, k0 = symbols('E I0 l_beam m_0 k_0', positive=True)
        
        default_data_dict = {
            self.I: [20 * I, 30 * I,60 * I,50 * I,40 * I,],
            self.h:[1,2,3],
            self.r:[0.1,0.2,0.3],
            self.a:[10],
            self.b:[20],
            self.phi : [10,20,30],
        }

        return default_data_dict


class SDOFWinchSystem(ComposedSystem):

    scheme_name = 'winch_complete.png'
    real_name = 'winch_cross_sec.png'

    I_s=Symbol('I_s', positive=True)
    I_k=Symbol('I_k', positive=True)
    I_1=Symbol('I_1', positive=True)
    i_1=Symbol('i_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    i_2=Symbol('i_2', positive=True)
    I_3=Symbol('I_3', positive=True)
    I_4=Symbol('I_4', positive=True)
    I_b=Symbol('I_b', positive=True)
    A=Symbol('A', positive=True)
    B=Symbol('B', positive=True)
    ivar=Symbol('t')
    phi=dynamicsymbols('\\varphi')
    dphi=dynamicsymbols('\\varphi',1)
    D=Symbol('D', positive=True)
    mu=Symbol('\\mu', positive=True)
    M_s=Symbol('M_s', positive=True)
    M_T=Symbol('M_T', positive=True)
    G=Symbol('G', positive=True)
    g=Symbol('g',positive=True)
    alpha=Symbol('\\alpha')
    pos1=phi
    qs=[phi]
    phi_1=phi
    phi_2=phi_1*i_1
    phi_3=phi_2*i_2
    def __init__(self,
                 I_s=None,
                 I_k=None,
                 I_1=None,
                 i_1=None,
                 I_2=None,
                 i_2=None,
                 I_3=None,
                 I_4=None,
                 I_b=None,
                 A=None,
                 B=None,
                 ivar=Symbol('t'),
                 phi=None,
                 dphi=None,
                 D=None,
                 mu=None,
                 M_s=None,
                 M_T=None,
                 G=None,
                 g=None,
                 alpha=None,
                 phi_1=None,
                 phi_2=None,
                 phi_3=None,
                 qs=None,
                 **kwargs):
        

        
        if I_s is not None: self.I_s = I_s
        if I_k is not None: self.I_k = I_k
        if I_1 is not None: self.I_1 = I_1
        if i_1 is not None: self.i_1 = i_1
        if I_2 is not None: self.I_2 = I_2
        if I_3 is not None: self.I_3 = I_3
        if i_2 is not None: self.i_2 = i_2
        if I_4 is not None: self.I_4 = I_4
        if I_b is not None: self.I_b = I_b
        if D is not None: self.D = D
        if G is not None: self.G = G
        if g is not None: self.g = g
        if mu is not None: self.mu = mu
        if M_T is not None: self.M_T = M_T
        if M_s is not None: self.M_s = M_s
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if phi_3 is not None: self.phi_3 = phi_3
        if alpha is not None: self.alpha = alpha
        if A is not None: self.A=A
        if B is not None: self.B=B
        if dphi is not None: self.dphi = dphi
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if phi_3 is not None: self.phi_3 = phi_3
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._engine = Disk(self.I_s, self.phi_1, qs=self.qs)
#         self.disk_k = Disk(I_k, phi_1, qs=qs)
        self._disk_1 = Disk(self.I_1, self.phi_1, qs=self.qs)
        self._disk_2 = Disk(self.I_2, self.phi_2, qs=self.qs)
        self._disk_3 = Disk(self.I_3, self.phi_2, qs=self.qs)
        self._disk_4 = Disk(self.I_4, self.phi_3, qs=self.qs)
        self._disk_B = Disk(self.I_b, self.phi_3, qs=self.qs)
        self._load = MaterialPoint(self.G/self.g, pos1=self.phi_3* self.D/2, qs=self.qs)
        self._friction_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * cos(self.alpha)*self.mu , qs=self.qs)
        self._force = Force(-self.M_T,pos1=self.phi_3,qs=self.qs)
        self._drive = Force(self.M_s,pos1=self.phi_1,qs=self.qs)
        self._gravity_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * sin(self.alpha) , qs=self.qs)


        components['engine'] = self._engine
        components['disk_1'] = self._disk_1
        components['disk_2'] = self._disk_2
        components['disk_3'] = self._disk_3
        components['disk_4'] = self._disk_4
        components['disk_B'] = self._disk_B
        components['load'] = self._load
        components['friction_comp'] = self._friction_comp
        components['force'] = self._force
        components['drive'] = self._drive
        components['gravity_comp'] = self._gravity_comp
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
            self.I_s: [20,30,40],
            self.I_k: [20,30,40],
            self.I_1: [20,30,40],
            self.I_2: [20,30,40],
            self.I_3: [20,30,40],
            self.I_4: [20,30,40],
            self.I_b: [20,30,40],
            self.i_1: [20,30,40],
            self.i_2: [20,30,40],
            self.D: [10,20,30],
            self.G: [10,20,30],
            self.g: [10,20,30],
            self.mu : [10,20,30],
            self.M_T : [10,20,30],
            self.M_s : [10,20,30],
#             self.phi_1: [10,20,30],
#             self.phi_2: [10,20,30],
#             self.phi_3  : [10,20,30],
            self.alpha  : [10,20,30],
        }

        return default_data_dict
        
    def steady_angular_velocity(self):
        obj=self
        eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)
        return omega_steady_state
    
    def reduced_torque(self):
        obj=self
        
        red_tor = -(obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())

        return delta_1

    def reduced_inertia(self):
        obj = self
        
        ans = obj.inertia_matrix()[0]

        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.02)-1)*I_r).n(5))
        return ik
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_3
        return st_val
    

    
class SDOFWinchSystemTest(ComposedSystem):


    scheme_name = 'winch_complete.png'
    real_name = 'winch_cross_sec.png'
    
    I_s=Symbol('I_s', positive=True)
    I_k=Symbol('I_k', positive=True)
    I_1=Symbol('I_1', positive=True)
    i_1=Symbol('i_1', positive=True)
    I_2=Symbol('I_2', positive=True)
    I_b=Symbol('I_b', positive=True)
    A=Symbol('A', positive=True)
    B=Symbol('B', positive=True)
    ivar=Symbol('t')
    phi=dynamicsymbols('\\varphi')
    dphi=dynamicsymbols('\\varphi',1)
    D=Symbol('D', positive=True)
    mu=Symbol('\\mu', positive=True)
    M_s=Symbol('M_s', positive=True)
    M_T=Symbol('M_T', positive=True)
    G=Symbol('G', positive=True)
    g=Symbol('g',positive=True)
    alpha=Symbol('\\alpha')
    pos1=phi
    qs=[phi]
    phi_1=phi
    phi_2=phi_1*i_1

    def __init__(self,
                 I_s=None,
                 I_k=None,
                 I_1=None,
                 i_1=None,
                 I_2=None,
                 I_b=None,
                 A=None,
                 B=None,
                 ivar=Symbol('t'),
                 phi=None,
                 dphi=None,
                 D=None,
                 mu=None,
                 M_s=None,
                 M_T=None,
                 G=None,
                 g=None,
                 alpha=None,
                 phi_1=None,
                 phi_2=None,
                 qs=None,
                 **kwargs):
        

        
        if I_s is not None: self.I_s = I_s
        if I_k is not None: self.I_k = I_k
        if I_1 is not None: self.I_1 = I_1
        if i_1 is not None: self.i_1 = i_1
        if I_2 is not None: self.I_2 = I_2
        if I_b is not None: self.I_b = I_b
        if D is not None: self.D = D
        if G is not None: self.G = G
        if g is not None: self.g = g
        if mu is not None: self.mu = mu
        if M_T is not None: self.M_T = M_T
        if M_s is not None: self.M_s = M_s
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if alpha is not None: self.alpha = alpha
        if A is not None: self.A=A
        if B is not None: self.B=B
        if dphi is not None: self.dphi = dphi
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        self._init_from_components(**kwargs)

    @property
    def components(self):

        components = {}

        self._engine = Disk(self.I_s, self.phi_1, qs=self.qs)
#         self.disk_k = Disk(I_k, phi_1, qs=qs)
        self._disk_1 = Disk(self.I_1, self.phi_1, qs=self.qs)
        self._disk_2 = Disk(self.I_2, self.phi_2, qs=self.qs)
#         self._disk_3 = Disk(self.I_3, self.phi_2, qs=self.qs)
#         self._disk_4 = Disk(self.I_4, self.phi_3, qs=self.qs)
        self._disk_B = Disk(self.I_b, self.phi_2, qs=self.qs)
        self._load = MaterialPoint(self.G/self.g, pos1=self.phi_2* self.D/2, qs=self.qs)
#         self._friction_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_3* self.D/2 * cos(self.alpha)*self.mu , qs=self.qs)
        self._force = Force(-self.M_T,pos1=self.phi_2,qs=self.qs)
        self._drive = Force(self.M_s,pos1=self.phi_1,qs=self.qs)
        self._gravity_comp = GravitationalForce(self.G/self.g, self.g, pos1=self.phi_2* self.D/2 * sin(self.alpha) , qs=self.qs)


        components['engine'] = self._engine
        components['disk_1'] = self._disk_1
        components['disk_2'] = self._disk_2
#         components['disk_3'] = self._disk_3
#         components['disk_4'] = self._disk_4
        components['disk_B'] = self._disk_B
        components['load'] = self._load
#         components['friction_comp'] = self._friction_comp
        components['force'] = self._force
        components['drive'] = self._drive
        components['gravity_comp'] = self._gravity_comp
        
        return components

    def symbols_description(self):
        self.sym_desc_dict = {
            self.I: r'crank moment of inertia',
            # self.k_beam: r'Beam stiffness',
            # self.g: r'gravitational field acceleration'
        }

        return self.sym_desc_dict

    def get_default_data(self):
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
            self.I_s: [20,30,40],
            self.I_k: [20,30,40],
            self.I_1: [20,30,40],
            self.I_2: [20,30,40],
            self.I_b: [20,30,40],
            self.i_1: [20,30,40],
            self.D: [10,20,30],
            self.G: [10,20,30],
            self.g: [10,20,30],
            self.mu : [10,20,30],
            self.M_T : [10,20,30],
            self.M_s : [10,20,30],
            self.phi_1: [10,20,30],
            self.phi_2: [10,20,30],
            self.alpha  : [10,20,30],
        }

        return default_data_dict
        
    def steady_angular_velocity(self):
        obj=self
        eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)
        return omega_steady_state
    
    def reduced_torque(self):
        obj=self
        
        red_tor = -(obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())

        return delta_1

    def reduced_inertia(self):
        obj = self
        
        ans = obj.inertia_matrix()[0]

        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.02)-1)*I_r).n(5))
        return ik
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_2
        return st_val

#DONE - Kuba
class SDOFDrivetrainVehicleSystem(ComposedSystem):


    #scheme_name = 'vehicle_drivetrain.png'
    #real_name = 'vehicle_drivetrain_real.PNG'

    R=Symbol('R')
    M_m=Symbol('M_m', positive=True)
    omega_m=Symbol('\\omega_m', positive=True)
    I_m=Symbol('I_m', positive=True)
    I=Symbol('I', positive=True)
    I_k=Symbol('I_k', positive=True)
    omega_1=Symbol('\\omega_1', positive=True)
    I_a=Symbol('I_a', positive=True)
    z_1=Symbol('z_1', positive=True)
    I_1=Symbol('I_1', positive=True)
    z_2=Symbol('z_2', positive=True)
    I_2=Symbol('I_2', positive=True)      
    I_b=Symbol('I_b', positive=True)
    z_3=Symbol('z_3', positive=True)
    I_3=Symbol('I_3', positive=True)
    z_4=Symbol('z_4', positive=True)
    I_4=Symbol('I_4', positive=True)
    omega_3=Symbol('omega_3', positive=True)
    I_c=Symbol('I_c',positive=True)
    I_d=Symbol('I_d', positive=True)
    I_w=Symbol('I_w', positive=True)
    r=Symbol('r', positive=True)
    m=Symbol('m', positive=True)
    i_1=Symbol('i_1',positive=True)
    i_2=Symbol('i_2',positive=True)
    A=Symbol('A', positive=True)
    B=Symbol('B', positive=True)
    phi=dynamicsymbols('\\varphi')
    phi_1=dynamicsymbols('\\varphi_1')
    phi_2=dynamicsymbols('\\varphi_2')
    phi_3=dynamicsymbols('\\varphi_3')
    dphi=dynamicsymbols('\\varphi',1)
    alpha=Symbol('alpha')
    c=Symbol('c')
    g=Symbol('g')    
    
    def __init__(self,
                 R=None,
                 M_m=None,
                 omega_m=None,
                 I_m=None,
                 I=None,
                 I_k=None,
                 omega_1=None,
                 I_a=None,
                 z_1=None,
                 I_1=None,
                 z_2=None,
                 I_2=None,
                 I_b=None,
                 z_3=None,
                 I_3=None,
                 z_4=None,
                 I_4=None,
                 omega_3=None,
                 I_c=None,
                 I_d=None,
                 I_w=None,
                 r=None,
                 m=None,
                 i_1=None,
                 i_2=None,
                 A=None,
                 B=None,
                 phi=None,
                 phi_1=None,
                 phi_2=None,
                 phi_3=None,
                 dphi=None,
                 alpha=None,
                 c=None,
                 g=None,
                 ivar=Symbol('t'),
                 **kwargs):
        


        if M_m is not None: self.M_m = M_m
        if omega_m is not None: self.omega_m = omega_m
        if I_1 is not None: self.I_1 = I_1
        if i_1 is not None: self.i_1 = i_1
        if I_2 is not None: self.I_2 = I_2
        if I_3 is not None: self.I_3 = I_3
        if i_2 is not None: self.i_2 = i_2
        if I_4 is not None: self.I_4 = I_4
        if I_a is not None: self.I_a = I_a
        if I_b is not None: self.I_b = I_b
        if I_c is not None: self.I_c = I_c
        if I_d is not None: self.I_d = I_d
        if I_k is not None: self.I_k = I_k
        if I_m is not None: self.I_m = I_m
        if I_w is not None: self.I_w = I_w
        if phi_1 is not None: self.phi_1 = phi_1
        if phi_2 is not None: self.phi_2 = phi_2
        if phi_3 is not None: self.phi_3 = phi_3
        if alpha is not None: self.alpha = alpha
        if A is not None: self.A = A
        if B is not None: self.B = B
        if dphi is not None: self.dphi = dphi
        if r is not None: self.r = r
#        if D is not None: self.D = 2*self.r
        if m is not None: self.m = m
        if c is not None: self.c = c
        if g is not None: self.g = g
        self.ivar=ivar
        self.qs = [self.phi] 

        self._init_from_components(**kwargs)
    #@cached_property  //not defined
    @property
    def components(self):
        components = {}
        
#        pos1=phi
        #qs=[self.phi]
        self.phi_1=self.phi
##         i_1=z_2/z_1
##         i_2=z_4/z_3
        self.phi_2=self.phi_1/self.i_1
        self.phi_3=self.phi_2/self.i_2
  
        
        self.engine_force = Force(self.M_m, pos1=self.phi_1,qs=[self.qs])
        self.engine_inertia = Disk(self.I_m, self.phi_1, qs=[self.qs])
#         self.flywheel = Disk(I, phi_1, qs=qs)
        self.clutch = Disk(self.I_c, self.phi_1, qs=[self.qs])
        self.driveshaft_1 = Disk(self.I_a, self.phi_1, qs=[self.qs])
        self.transmission_i = Disk(self.I_1, self.phi_1, qs=[self.qs])
        self.transmission_o = Disk(self.I_2, self.phi_2, qs=[self.qs])
        self.driveshaft_2 = Disk(self.I_b, self.phi_2, qs=[self.qs])
        self.diff_i = Disk(self.I_3, self.phi_2, qs=[self.qs])
        self.diff_o = Disk(self.I_4, self.phi_3, qs=[self.qs])
        self.axleshaft = Disk(self.I_d, self.phi_3, qs=[self.qs])
        self.wheels = Disk(4*self.I_w, self.phi_3, qs=[self.qs])
        self.mass = MaterialPoint(self.m, pos1=self.phi_3*self.r , qs=[self.qs])
        self.gravity = Force(-self.m*self.g*sin(self.alpha), pos1=self.phi_3*self.r,qs=[self.qs])#GravitationalForce(m, g, pos1=phi_3* r * sin(alpha) , qs=qs)
        self.air_force = Damper(self.c,pos1=self.phi_3*self.r,qs=[self.phi_1]) #(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=base_frame)

        components['_engine_force'] = self.engine_force
        components['_engine_inertia'] = self.engine_inertia
        components['_clutch'] = self.clutch
        components['_driveshaft_1'] = self.driveshaft_1
        components['_transmission_i'] = self.transmission_i
        components['_transmission_o'] = self.transmission_o
        components['_driveshaft_2'] = self.driveshaft_2
        components['_diff_i'] = self.diff_i
        components['_diff_o'] = self.diff_o
        components['_axleshaft'] = self.axleshaft
        components['_wheels'] = self.wheels
        components['_mass'] = self.mass
        components['_gravity'] = self.gravity
        components['_air_force'] = self.air_force


        
        return components                


    def get_default_data(self):
        
        E, I_s,I_k,I_b,I_1,I_2,I_3,I_4,i_1, i_2,D,G,g,mu,M_T,M_s= symbols('E I_s I_k I_b I_1 I_2 I_3 I_4 i_1 i_2 D G g \mu M_T M_s', positive=True)
        
        default_data_dict = {
#             self.I_s: [20,30,40],
#             self.I_k: [20,30,40],
#             self.I_1: [20,30,40],
#             self.I_2: [20,30,40],
#             self.I_3: [20,30,40],
#             self.I_4: [20,30,40],
#             self.I_b: [20,30,40],
#             self.i_1: [20,30,40],
#             self.i_2: [20,30,40],
#             self.D: [10,20,30],
#             self.G: [10,20,30],
#             self.g: [10,20,30],
#             self.mu : [10,20,30],
#             self.M_T : [10,20,30],
#             self.M_s : [10,20,30],
#             self.phi_1: [10,20,30],
#             self.phi_2: [10,20,30],
#             self.phi_3  : [10,20,30],
        }

        return default_data_dict
    
    def linear_vehicle_velocity(self):
        obj=self
        lin_veh_vel=omega_steady_state*vel_ratio
        return lin_veh_vel
    
    def velocity_ratio(self):
        obj=self
        vel_ratio=self.phi_3*self.r/self.phi
        return vel_ratio
    

    def steady_angular_velocity(self):
        obj=self


        eom_sol = obj._eom_solution()
        #omega_steady_state =  (eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo))

        #display(obj._eoms[0])
        #eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
        #display(eoms_num)
        #eom_sol = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
        omega_steady_state =  N(eom_sol.rhs.diff(obj.ivar).subs(obj.ivar,oo),3)

        return omega_steady_state
    
    def steady_angular_velocity_rev_per_min(self):
        obj=self
        ang_vel = obj.steady_angular_velocity()
        return ang_vel* 60/360
    
    def reduced_torque(self):
        obj=self
        red_tor = (obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)]))
        return red_tor
    
    def delta_1(self):
        obj=self
        M_Z = obj._eoms[0].doit().subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data).subs([(obj.q[0].diff(obj.ivar,obj.ivar),0),(obj.q[0].diff(obj.ivar),0)])
        delta_1 = 0.4*(obj.A-Abs(M_Z))/(self.steady_angular_velocity()**2*self.reduced_inertia())
#         0.4*(num_data.loc[case_no,'A']-Abs(M_Z))/(omega_ust**2*I_r)
        return delta_1

    def reduced_inertia(self):
        obj = self
        ans = obj.inertia_matrix()[0]
        return ans

    def flywheel_moment_of_inertia(self):
        obj=self
        delta_1=self.delta_1()
        I_r=self.reduced_inertia()

        ik=((((delta_1/0.004)-1)*I_r).n(5))
        return ik
    
    def startup_mass_velocity(self):
        obj=self
        st_val=obj.D/2*obj.phi_3
        return st_val
    
    def reducted_mass(self):
        obj=self
        
        return obj.flywheel_moment_of_inertia()/(obj.startup_mass_velocity()*obj.phi1)**2
    
    def reducted_force(self):
        obj=self
        
        return obj.reduced_torque()/obj.startup_mass_velocity()


    def _eom_solution(self):
        obj=self

        if obj.__cached_solution is None:
            eoms_num = N(obj._eoms[0].subs(obj.M_s,obj.A-obj.B*obj.dphi).subs(obj._given_data),3)
            EOM_dsolve = dsolve(eoms_num, obj.q[0], ics={obj.q[0].subs(obj.ivar, 0): 0, obj.q[0].diff(obj.ivar).subs(obj.ivar, 0): 0})
            obj.__cached_solution = EOM_dsolve
        else:
            EOM_dsolve = obj.__cached_solution

        return EOM_dsolve        
    
