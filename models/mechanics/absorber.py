from sympy import * # provides mathematical interface for symbolic calculations

from dynpy import * # enables mechanical models for mathematical modelling
from dynpy.dynamics import LagrangesDynamicSystem, HarmonicOscillator
from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point





from dynpy.models.mechanics.pendulum import Pendulum
from dynpy.models.elements import MaterialPoint, Spring, GravitationalForce, Disk, RigidBody2D, Damper, PID, Excitation, Force, base_frame,base_origin



class ComposedSystem(HarmonicOscillator):
    """Base class for all systems

    """
    scheme_name = 'damped_car_new.PNG'
    real_name = 'car_real.jpg'
    _default_args = ()

    @classmethod
    def _scheme(cls):

        path = __file__.replace('systems.py', 'images/') + cls.scheme_name

        return path

    @classmethod
    def _real_example(cls):

        path = __file__.replace('systems.py', 'images/') + cls.real_name

        return path

    @classmethod
    def preview(cls, example=False):
        if example:
            path = cls._real_example()

        else:
            path = cls._scheme()

        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))

    def calculations_steps(self,preview=True,system=None,code=False):
        
#         latex_store=AutoBreak.latex_backend
#         AutoBreak.latex_backend = latex
        
        print('zlo')
        print(inspect.getsource(self.__class__))

        
        doc_model=super().calculations_steps(preview=True,code=code)
        
        
#         AutoBreak.latex_backend = latex_store
        return doc_model
            

    
    
    def get_default_data(self):
        return None

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()


        
        if default_data_dict:
            parameters_dict = {
                key: random.choice(items_list)
                for key, items_list in default_data_dict.items()
            }
        else:
            parameters_dict=None

        return parameters_dict
    
    def linearized(self, x0=None, op_point=False, hint=[], label=None):

        return type(self).from_system(super().linearized(x0=x0,op_point=op_point,hint=hint,label=label))

t,f = symbols('t, f')
g, m1, m2, lam, l1, l2=symbols('g, m_alpha, m_beta, lambda, l_alpha, l_beta')

    
class VariableMassDoublePendulum(ComposedSystem):
    
    t0=Symbol('t_0')
    flow_coeff = Symbol('lambda')
    
    

    scheme_name = 'MDOFTriplePendulum.PNG'
    real_name = 'TriplePendulum_real.jpg'
    
    
    

    def __init__(self,
                 m=Symbol('m', positive=True),
                 m1=Symbol('m_alpha_0', positive=True),
                 m_f1=dynamicsymbols('m_alpha'),
                 m2=Symbol('m_beta_0', positive=True),
                 m_f2=dynamicsymbols('m_beta'),
                 m_0 = Symbol('m_0', positive=True),
                 l_1=Symbol('l_alpha', positive=True),
                 l_2=Symbol('l_beta', positive=True),
                 l_3=Symbol('l_3', positive=True),
                 w=Symbol('w', postive=True),
                 k=Symbol('k', postive=True),
                 ls=Symbol('l_s', postive=True),
                 u=Symbol('u', postive=True),
                 x0=Symbol('x_0', postive=True),
                 l0=Symbol('l_0', postive=True),
                 u0=Symbol('u_0', postive=True),
                 x_k=Symbol('x_k', postive=True),
                 b=Symbol('b', postive=True),
                 omega=Symbol('omega',positive=True),
                 f=Symbol('f',positive=True),
                 rho_f=Symbol('rho_f',positive=True),
                 h_fa=Symbol('h_f_alpha',positive=True),
                 h_fb=Symbol('h_f_beta',positive=True),
                 A_pipe1=Symbol('A_alpha',positive=True),
                 A_pipe2=Symbol('A_beta',positive=True),
                 flow_coeff = flow_coeff,
                 g=Symbol('g', positive=True),
                 t0=t0,
                 phi_1=dynamicsymbols('alpha'),
                 phi_2=dynamicsymbols('beta'),
                 #phi_3=dynamicsymbols('varphi_3'),
                 #phi_u=dynamicsymbols('varphi_u'),
                 #phi_l=dynamicsymbols('varphi_l'),
                 qs=dynamicsymbols('alpha, beta'),
                 ivar=Symbol('t'),
                 m_t1=Symbol('m_alpha_t',positive=True),
                 m_t2=Symbol('m_beta_t',positive=True),
                 m_t=Symbol('m_t',positive=True),
                 **kwargs):

        self.m = m
        self.m1 = m1
        self.m2 = m2
        self.m_f1 = m_f1
        self.m_f2 = m_f2
        self.m_0 = m_0
        self.l_1 = l_1
        self.l_2 = l_2
        self.l_3 = l_3
        self.w = w
        self.k = k
        self.ls = ls
        #self.u = u
        self.x0 = x0
        self.l0 = l0
        self.u0 = u0
        self.x_k = x_k
        self.b = b
        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.omega = omega
        self.f = f
        self.rho_f = rho_f
        self.h_fa = h_fa
        self.A_pipe1 = A_pipe1
        self.A_pipe2 = A_pipe2
        self.h_fb = h_fb
        #self.phi_3 = phi_3
        #self.phi_u = phi_u
        #self.phi_l = phi_l
        self.g = g
        self.t0=t0
        self.flow_coeff=flow_coeff
        self.m_t1=m_t1
        self.m_t2=m_t2
        self.m_t=m_t
        
        t = ivar
        self.dphi_1 = phi_1.diff(t)
        self.dphi_2 = phi_2.diff(t)
        

        self.trans_expr=((S.One/2-atan(flow_coeff*(t-t0))/pi))
        
        #g, m1, m2, lam, l1, l2=symbols('g, m_alpha, m_beta, lambda, l_alpha, l_beta')

        # water masses as symbols
        m1_0 , m2_0 = m1 , m2
        m1_pipe , m2_pipe = m_f1 , m_f2
        self.m1_pipe , self.m2_pipe = m1_pipe , m2_pipe

        A_pipe1=3.14*(0.14)**2 # m2
        A_pipe2=3.14*(0.14)**2 #  m2

        rho_f=1000 # kg/m3

        h_f1=0
        h_f2=0

        phi1, phi2=qs

        xc1=l1*sin(phi_1)
        yc1=l1*cos(phi_1)

        xfc1=S.Zero;
        yfc1=S.Zero;

        x_1=l1*sin(phi_1)+l2*sin(phi_2)
        y_2=l1*cos(phi_1)+l2*cos(phi_2)

        xfc2=S.Zero;
        yfc2=S.Zero;

        ### velocities cal

#         w, k, ls, u = symbols('w, k, ls, u')
#         x0, l0 = symbols('x0, l0',positive=True)

        w=self.omega
        x_k=u0*sin(w*t)

        xA=sin(phi1)*l_1/2**2
        yA=cos(phi1)*l_1/2**2
        xB=x_k
        yB=l_1

        ls=sqrt((xA-xB)**2-(yA-yB)**2)
        u=ls-l0

        Vs=k/2*(ls-l0)**2

        #linearized form of L

        N = ReferenceFrame('N')
        P = Point('P')
        P1 = Point('P1')

        P.set_vel(N, phi1.diff(t) * N.x)
        P1.set_vel(N, phi2.diff(t) * N.x)
        fl = [(P, -b * phi1.diff(t) * N.x), (P1, -b * phi2.diff(t) * N.x)]

        x_2 = sin(phi_1)*l_1 + sin(phi_2)*l_2
        y_2 = cos(phi_1)*l_1 + cos(phi_2)*l_2


        self.Pendulum1 = Pendulum(m1, g, l_1, angle=phi_1, qs=[phi_1]) + Pendulum(m1_pipe, g, l_1, angle=phi_1, qs=[phi_1])

        #second arm
        self.material_point_11 = MaterialPoint(m2+m2_pipe, x_2, qs=[phi_1, phi_2])
        self.material_point_21 = MaterialPoint(m2+m2_pipe, y_2, qs=[phi_1, phi_2])
        self.gravity_1 = GravitationalForce(m2, g, pos1=-y_2, qs=[phi_2]) + GravitationalForce(m2_pipe, g, pos1=-y_2, qs=[phi_2])
        self.spring = Spring(k*0,ls-l0,qs=[phi_1]).linearized()
        self.damper_1 = Damper(b,phi_1,qs=[phi1])
        self.damper_2 = Damper(b,phi_2,qs=[phi2])
        self.force = Force(u0*sin(w*t)*k*l_1,phi_1,qs=[phi1])
#         self.material_point_12 = MaterialPoint(m3, x_3, qs=[phi_1, phi_2, phi_3])
#         self.material_point_22 = MaterialPoint(m3, y_3, qs=[phi_1, phi_2, phi_3])
#         self.gravity_2 = GravitationalForce(m3, g, pos1=-y_3, qs=[phi_3])

        system = self.Pendulum1 + self.material_point_11 + self.material_point_21 + self.gravity_1 + self.spring + self.damper_1 + self.damper_2+self.force #+ self.gravity_2
        super().__init__(system(qs),**kwargs)
        
        self._T=self.material_point_11.lagrangian()+self.material_point_21.lagrangian() # sum([comp.lagrangian()  for  comp   in [self.material_point_21,self.material_point_11]])
        self._V=-(self.spring.lagrangian()+self.gravity_1.lagrangian())
        

    def get_default_data(self):


        m0, l0 = symbols('m_0 l_0', positive=True)

        default_data_dict = {
            self.m1: [S.Half * m0, 1 * m0, 2 * m0, 1 * m0, S.Half * m0],
            self.m2: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],
            self.m3: [1 * m0, 2 * m0, S.Half * m0, 1 * m0, 2 * m0],

            self.l_1: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_2: [1 * l0, 2 * l0, S.Half * l0, 2 * l0, S.Half * l0],
            self.l_3: [2 * l0, 4 * l0, S.Half * l0, 2 * l0, S.Half * l0],

            self.phi_1:[self.phi_u,0],
            self.phi_2:[self.phi_u,self.phi_l],
            self.phi_3:[self.phi_l],
        }

        return default_data_dict

    def get_random_parameters(self):

        default_data_dict = self.get_default_data()

        parameters_dict = {
            key: random.choice(items_list)
            for key, items_list in default_data_dict.items()
        }

        if parameters_dict[self.phi_2] == parameters_dict[self.phi_3]:

            parameters_dict[self.phi_2] = self.phi_u


        return parameters_dict

    def symbols_description(self):
        self.sym_desc_dict = {
            self.m: r'mass of the system',
            self.m1: r'initial mass of the upper pendulum member',
            self.m2: r'initial mass of the lower pendulum member',
            self.m_f1: r'mass of the upper pendulum',
            self.m_f2: r'mass of the lower pendulum',
            self.m_f1.diff(t): r'mass fluid flow of the upper pendulum',
            self.m_f2.diff(t): r'mass fluid flow of the lower pendulum',
            self.m_t: r'total mass of the system',
            self.m_t2: r'total mass of the lower pendulum',
            self.m_0: r'transferred damping mass',
            self.l_1: r'length of the upper pendulum',
            self.l_2: r'length of the lower pendulum',
            self.l_3: r'length of the damper pendulum',
            self.w: r'angular velocity',
            self.k: r'spring stiffness',
            self.ls: r'length of the spring',
            #self.u: r'difference of spring length and initial spring length',
            self.x0: r'initial position of the system',
            self.l0: r'initial spring length',
            self.u0: r'initial difference of spring length and starting spring length',
            self.x_k: r'forcing force',
            self.b: r'damping coefficient value',
            self.omega: r'excitation frequency',
            #self.f: r'????',
            self.g: r'acceleration of gravity',
            self.rho_f: r'fluid density',
            self.phi_1: r'angle of the upper pendulum swing',
            self.phi_2: r'angle of the lower pendulum swing',
            self.phi_1.diff(t): r'generalized velocity of the upper pendulum member',
            self.phi_2.diff(t): r'generalized velocity of the lower pendulum member',
            self.phi_1.diff(t,t): r'generalized acceleration of the upper pendulum member',
            self.phi_2.diff(t,t): r'generalized acceleration of the lower pendulum member',
            self.t0: r'activation time',
            self.flow_coeff: r'fluid flow rate',
            self.ivar: r'time',
            self.h_fa:r'centre of gravity of the higher member',
            self.h_fb:r'centre of gravity of the lower member',
            self.A_pipe1:r'cross-section area of the upper member',
            self.A_pipe2:r'cross-section area of the lower member',
        }
        return self.sym_desc_dict
    
    def flow_equation(self):
        return_eq=Eq(self.m2,(self.m_0-self.m_0*self.trans_expr).expand())
        
        return return_eq
    
    def upper_mass_equation(self):
        return_eq=Eq(self.m1,(self.m_0*self.trans_expr).expand())
        
        return return_eq