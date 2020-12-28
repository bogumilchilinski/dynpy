import importlib
import dynpy as dyn
importlib.reload(dyn)
from sympy.physics.vector import dynamicsymbols
from sympy import *
from sympy.physics.mechanics import *
import sympy as sym

mechanics_printing()



t=Symbol('t') #independent variable - time

### system parametes
m_fr, m_rear,m_3, k_r, k_rt, k_f, k_ft = symbols('m_fr, m_r, M, k_r, k_rt, k_f, k_ft',positive=True)
m,k,g,F_1,F_2,Omega,F,R,v0,u0= symbols('m,k,g,F_1,F_2,Omega, F_ R, v_0,u_0',positive=True)
I_ch, I_w , z_c3,l_fr,l_rear= symbols('I_chair, I_wheel, z_c3, l_fr, l_r',positive=True)
m_RC, I_RC, l_RC, k_RC, phi0 = symbols('m_RC, I_RC, l_RC, k_RC, varphi_0',positive=True)
m_w, I_wrc, r_w, k_w, k_fix,k_tire = symbols('m_w, I_wRC, R_w, k_w, k_fix, k_t',positive=True)
c,c_mu,c_lam= symbols('c,c_mu, c_lambda',positive=True)

pm=symbols('PM')

var_static = {m_fr:'masa koło przednie wózka',
              m_rear:'masa koło tylne wózka',
              m_3:'???',
              k_r:'sztywność koła tylnego wózka',
              k_rt:'sztywność ogumienia koła tylnego wózka',
              k_f:'sztywność koła przedniego wózka',
              k_ft:'sztywność ogumienia koła przedniego wózka',
              m:'masa wózka',
              k:'sztywność wózka',
              g:'przyspieszenie ziemskie',
              F_1:'siła ???',
              F_2:'siła ???',
              Omega:'???',
              F:'siła ???',
              R:'promień tylnego koła wózka',
              v0:'prędkość początkowa wózka',
              u0:'???',
              I_ch:'moment bezwładności wózka',
              I_w:'moment bezwładności koła',
              z_c3:'pionowa odległość środka ciężkości od środka obrotu wózka',
              l_fr:'długość początkowa sprężyny przedniego zawieszenia',
              l_rear:'długość początkowa sprężyny tylnego zawieszenia',
              m_RC:'masa wahacza napędu RapidChair',
              I_RC:'moment bezwładności napędu RapidChair',
              l_RC:'długość wahacza napędu RapidChair',
              k_RC:'sztywność napędu RapidChair',
              phi0:'początkowy kąt położenia osi wahacza względem płaszczyzny poziomej',
              m_w:'masa koła napędu RapidChair',
              I_wrc:'moment bezwładności koła napędu RapidChair',
              r_w:'promień koła napędu RapidChir', #do sprawdzenia, nie jestem pewien
              k_w:'sztywność koła napędu RapidChair',
              k_fix:'sztywność mocowania',
              k_tire:'sztywność ogumienia koła napędu RapidChair',
              c:'???',
              c_mu:'???',
              c_lam:'???',
              pm:'???'}

t_l, delta_t,l_ramp,l_bumps = symbols('t_l,dt,l_ramp, l_bumps',positive=True)
t0 = symbols('t_0')

# generalized coordinates
x,z_fr,z_rear,z, phi,phi_rc,theta,z_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc')
dx,dz_fr,dz_rear,dz,dphi,dphi_rc,dtheta,dz_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, dz_wrc', 1)

var_dynamic = {x:'przemieszczenie poziome wózka',
               z_fr:'przemieszczenie pionowe koła przedniego wózka',
               z_rear:'przemieszczenie pionowe koła tylnego wózka',
               z:'przemieszczenie pionowe punktu mocowania koła wózka',
               phi:'przemieszczenie kątowe wahacza wózka',
               phi_rc:'przemieszczenie kątowe napędu RapidChair',
               theta:'przemieszczenie kątowe koła napędu RapidChair',
               dx:'prędkość liniowa wózka',
               dz_fr:'prędkość liniowa koła przedniego wózka',
               dz_rear:'prędkość liniowa koła tylnego wózka',
               dz:'prędkość liniowa wózka',
               dphi:'prędkość kątowa wahacza wózka',
               dphi_rc:'prędkość kątowa napędu RapidChair',
               dtheta:'prędkość kątowa koła napędu RapidChair'}

q=[x,phi,z,z_fr,z_rear,phi_rc]
u=[dx,dphi,dz,dz_fr,dz_rear,dphi_rc]
Y = q + u
nDOF=len(q)


# u_rr=u0/2*(1+(Heaviside(x-t0-2*R,0.5) - Heaviside(x-t0-delta_t-2*R,0.5) ))#road profile given in time domain
# u_rf=u0/2*(1+(Heaviside(x-t0,0.5) - Heaviside(x-t0-delta_t,0.5) ))#road profile given in time domain

u_rf_bump=u0/2*(1+2/pi*atan(5*(x-t0)))-u0/2*(1+2/pi*atan(5*(x-t0-t_l)))
u_rr_bump=u_rf_bump.subs(x,x-2*R) #road profile given in time domain
 #road profile given in time domain

u_rf_mul_bumps=sum(u_rf_bump.subs(x,x-ii*l_bumps) for ii in range(5))
u_rr_mul_bumps=u_rf_mul_bumps.subs(x,x-2*R) #road profile given in time domain
 #road profile given in time domain
    
    
u_rf_ramp=u0/l_ramp*x
u_rr_ramp=u_rf_ramp.subs(x,x-2*R)

u_road=Function('u_r')
u_rf=u_road(x)
u_rr=u_road(x-2*R)

u_rr=u0*(sin(x-2*R))  #road profile given in time domain
u_rf=u_rr.subs(x,x+2*R) #road profile given in time domain
#












#u_fr=Function('u_fr')(phi,z_fr)

# u_fr=u_fr_expr

# u_rear_expr=sqrt((z+R*sin(phi)+l_rear-z_rear)**2 + (R-R*cos(phi))**2)-l_rear
# #u_rear=Function('u_rear')(phi,z_fr)

# u_rear=u_rear_expr



# D=S.One/2*c*2*(dx**2) + c/2*sum(vel**2 for vel in u)







x_m3=x+z_c3*phi
x_RC=x+l_RC/2*cos(phi0+phi_rc)
z_RC=z+R*phi+l_RC/2*sin(phi0+phi_rc)


# Dynamic description of the chair
T_body = S.One/2*m_3*x_m3.diff(t)**2 + S.One/2 * m_3* dz**2 # Ek ramy wózka w ruchu postępowym
T_rot =  S.One/2 * I_ch* dphi**2  # Ek ramy wózka w ruchu obrotowym
T_rear = S.One/2 * m_rear* dz_rear**2 + S.One/2 * m_rear * dx**2    # Ek tylniego koła Ir (mr)
T_fr = S.One/2 * m_fr* dz_fr**2 + S.One/2 * m_fr * dx**2    # Ek przedniego koła If (mf)
T_wheel = S.One/2 * I_w * (dx/R)**2

# Dynamic description of the drive attachment
# T_RC = S.One/2 * m_RC* (x_RC.diff(t)**2+z_RC.diff(t)**2) #RC drive kinetic energy
T_RC_arm=S.One/2 * m_RC* (x_RC.diff(t)**2+z_RC.diff(t)**2)
T_RC_wheel_z=S.One/2 * m_w * dz_wrc**2
T_RC_wheel_rev=S.One/2 * I_wrc * dtheta**2

T_RC_3dof=T_RC_arm+T_RC_wheel_z+T_RC_wheel_rev

V_chair_g=m_fr*g*z_fr + m_rear*g*z_rear + m_3 * g*(z + z_c3*(cos(phi))) # EPG wózka względem punktu cięzkości wózka

# nonlinear springs deflection
u_rear=sqrt((z+R*sin(phi)+l_rear-z_rear)**2 + (R-R*cos(phi))**2)-l_rear#
u_fr=sqrt((z-R*sin(phi)+l_fr-z_fr)**2 + (R-R*cos(phi))**2)-l_fr

# linear springs deflection
u_rear=((z+R*phi+l_rear -z_rear)**2 )-l_rear#
u_fr=  ((z-R*phi+l_fr   -z_fr)**2  )-l_fr

# spring potentials
V_rear = S.One/2*k_rt*(z_rear-u_rr)**2  + S.One/2*k_r*(u_rear)**2 # Ep tylnich prętów ramy względem ich sprężystości
V_fr = S.One/2*k_ft*(z_fr-u_rf)**2  + S.One/2*k_f*(u_fr)**2       # Ep przednich prętów ramy względem ich sprężystości



V_rc_g=-m_RC*g*z_RC
V_rc_fix =  S.One/2 * k_fix * phi_rc**2 + S.One/2 *k_tire * (l_RC*phi_rc)**2 
V_rc_w = m_w*g*z_wrc

V_RC_3dof=V_rc_fix+V_rc_g+V_rc_w

T_chair5dof = T_body + T_rot + T_rear + T_fr + T_wheel
V_chair5dof = V_chair_g + V_rear + V_fr


T_chair5dof_rc3dof = T_chair5dof + T_RC_3dof
V_chair5dof_rc3dof = V_chair5dof + V_RC_3dof


dof2 = {z_rear:z,z_fr:z,phi:0} # 2DoF chair
dof3_norev = {z_rear:z,z_fr:z} # 3 DoF w/o rev chair
dof3_rev={z_fr:z,phi:0} #3 DoF w/rev chair

rcdof1={z_wrc:0,theta:0} #RC 2dof
rcdof2={z_wrc:0}


D=c_mu*T_5dof
# D+= c_lam* + S.One/2*k_r*(((z+R*(phi)-z_rear)**2 ).diff(t))**2

# D+= c_lam* S.One/2*k_f*(((z-R*(phi)-z_fr)**2 ).diff(t))**2

D=D.doit().expand()


N = ReferenceFrame('N')
Pl = Point('P_l')
Pr = Point('P_r')
Pl.set_vel(N, (u[0]+u[1]*R*sin(Omega*t).diff(t)) * N.x+(u[2]*R*cos(Omega*t).diff(t)) * N.y)
Pr.set_vel(N, u[0] * N.x)


points_list=[Point('P_'+str(no))  for no,vel in enumerate(u)]
[points_list[no].set_vel(N,vel*N.x)  for no,vel in enumerate(u)]

FL = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*F*(cos(Omega*t)) * N.x+sin(Omega*t) * N.y)]+[(points_list[no],-D.diff(vel)*N.x)  for no,vel in enumerate(u)]   #,(Pr, R*F*cos(Omega*t) * N.x)]
FL = [(Pl, 2*F*N.x)]+[(points_list[no],-D.diff(vel)*N.x)  for no,vel in enumerate(u)]


L_chair5dof = T_chair5dof - V_chair5dof
# L_5dof=

# L_chair5dof_rc3dof = T_chair5dof_rc3dof-V_chair5dof_rc3dof


# L_rc1dof=

L_default_5dof= L_chair5dof

L_rc_3dof=T_RC_3dof-V_RC_3dof

L_chair5dof_rc3dof=T_chair5dof_rc3dof-V_chair5dof_rc3dof



qs_5dof = x, z_rear, z_fr, z, phi
qs_rc_3dof=z_wrc,theta,phi_rc
qs_chair5dof_rc_3dof=qs_5dof+qs_rc_3dof
# qs_3dof = x, z, phi
# qs_2dof = x, z




class Chair5DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_default_5dof, qs=qs_5dof, forcelist=FL, bodies=None, frame=N,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)

class RapidChair3DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_rc_3dof, qs=qs_rc_3dof, forcelist=FL, bodies=None, frame=N,
                   hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):

        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
             hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)

class Chair5DOFwithRC3DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_rc_3dof, qs=qs_rc_3dof, forcelist=FL, bodies=None, frame=N,
                   hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):

        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
             hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)


chair_5dof = Chair5DOF()('Chair 5DOF model')
chair_3dof = chair_5dof.subs(dof3_rev,method='direct').remove([phi,z_fr])('Chair 3DOF model with revolution')
chair_3dof = chair_5dof.subs(dof3_norev,method='direct').remove([z_rear,z_fr])('Chair 3DOF model without revolution')       
chair_2dof = chair_5dof.subs(dof2,method='direct').remove([z_rear,z_fr,phi])('Chair 2DOF model')

rapidchair_3dof = RapidChair3DOF()('RapidChair 3DOF model')
rapidchair_2dof=rapidchair_3dof.subs(rcdof2,method='direct').remove([z_wrc])('RapidChair 2DOF model')
rapidchair_1dof=rapidchair_3dof.subs(rcdof1,method='direct').remove([z_wrc,theta])('RapidChair 1DOF model')


chair5dof_rc3dof=Chair5DOFwithRC3DOF()
# class Chair3dof(dyn.LagrangesDynamicSystem):
#     def __init__(self, Lagrangian=L_default_3dof, qs=qs_3dof, forcelist=FL, bodies=None, frame=N,
#                        hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
#         super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
#                  hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)

        
# class Chair2dof(dyn.LagrangesDynamicSystem):
#     def __init__(self, Lagrangian=L_default_2dof, qs=qs_2dof, forcelist=FL, bodies=None, frame=N,
#                        hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
#         super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
#                  hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)
        


