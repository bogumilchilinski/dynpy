import importlib
from . import dynamics as dyn
importlib.reload(dyn)
from . import dynamics as dyn

from sympy.physics.vector import dynamicsymbols
from sympy import *
from sympy.physics.mechanics import *
import sympy as sym
import numpy as np
import pylatex
from pylatex import Command, NoEscape
import pint
ureg = pint.UnitRegistry()

mechanics_printing()



t=Symbol('t') #independent variable - time
a_sim=Symbol('a_sim',positive=True)
ivar=Symbol('t')

### system parameters
m_fr, m_rear,m_3, k_r, k_rt, k_f, k_ft, k_rot,amplitude,length,speed,axw = symbols('m_fr, m_r, M, k_r, k_rt, k_f, k_ft k_rot, A, L, v, W_wb ',positive=True)
m,k,g,F_1,F_2,Omega,F,R,v0,u0,l_l,l_r= symbols('m,k,g,F_1,F_2,Omega, F_0, R, v_0,u_0,l_l,l_r',positive=True)
I_ch, I_w , z_c3,l_fr,l_rear= symbols('I_chair, I_wheel, z_c3, l_fr, l_r',positive=True)
m_RC, I_RC, l_RC, k_RC, phi0 = symbols('m_RC, I_RC, l_RC, k_RC, varphi_0',positive=True)
m_w, I_wrc, r_w, k_w, k_fix,k_tire = symbols('m_w, I_wRC, R_w, k_w, k_fix, k_t',positive=True)
c,c_mu,c_lam= symbols('c,c_mu, c_lambda',positive=True)
a_ox, a_oz, a_rz ,a_rcz=symbols('a_ox, a_oz, a_rz ,a_rcz')

pm=symbols('PM')

t_l, delta_t,l_ramp,l_bumps = symbols('t_l,dt,l_ramp, l_bumps',positive=True)
t0 = symbols('t_0')




var_static = {
              m_fr:'masa przedniego koła wózka',
              m_rear:'masa tylnego koła wózka',
              m_3:'masa konstrukcji wózka wraz z niepełnosprawnym',
              k_r:'sztywność tylnego koła wózka',
              k_rt:'sztywność ogumienia tylnego koła wózka',
              k_f:'sztywność przedniego koła wózka',
              k_ft:'sztywność ogumienia przedniego koła wózka',
              m:'masa wózka',
              k:'sztywność wózka',
              g:'przyspieszenie ziemskie',
              F_1:'siła wymuszająca',
              F_2:'siła wymuszająca',
              Omega:'częstość wymuszenia',
              F:'siła przykładana do wózka przez niepełnosprawnego',
              R:'promień tylnego koła wózka',
              v0:'prędkość początkowa wózka',
              u0:'ammplituda profilu drogi',
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
              pm:'czas działania siły na wózek',
             }








# generalized coordinates
x,z_fr,z_rear,z, phi,phi_rc,theta,z_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc')
dx,dz_fr,dz_rear,dz,dphi,dphi_rc,dtheta,dz_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc', 1)




var_dynamic = {x:'przemieszczenie poziome wózka',
               z_fr:'przemieszczenie pionowe koła przedniego wózka',
               z_rear:'przemieszczenie pionowe koła tylnego wózka',
               z:'przemieszczenie pionowe punktu mocowania koła wózka',
               z_wrc:'przemieszczenie pionowe koła napędu RapidChair',
               phi:'przemieszczenie kątowe wahacza wózka',
               phi_rc:'przemieszczenie kątowe napędu RapidChair',
               theta:'przemieszczenie kątowe koła napędu RapidChair',
               dx:'prędkość liniowa wózka',
               dz_fr:'prędkość liniowa koła przedniego wózka',
               dz_rear:'prędkość liniowa koła tylnego wózka',
               dz:'prędkość liniowa wózka',
               dz_wrc:'prędkość liniowa koła napędu RapidChair',
               dphi:'prędkość kątowa wahacza wózka',
               dphi_rc:'prędkość kątowa napędu RapidChair',
               dtheta:'prędkość kątowa koła napędu RapidChair',
              }


symbols_description={**var_static,**var_dynamic}


q=[x,phi,z,z_fr,z_rear,phi_rc]
u=[dx,dphi,dz,dz_fr,dz_rear,dphi_rc]
Y = q + u
nDOF=len(q)


# u_rr=u0/2*(1+(Heaviside(x-t0-2*R,0.5) - Heaviside(x-t0-delta_t-2*R,0.5) ))#road profile given in time domain
# u_rf=u0/2*(1+(Heaviside(x-t0,0.5) - Heaviside(x-t0-delta_t,0.5) ))#road profile given in time domain

u_rf_bump_on=u0/2*(1+2/pi*atan(5*(t-t0)))-u0/2*(1+2/pi*atan(5*(t-t0-t_l)))
u_rf_bump=S.Zero
u_rr_bump_on=u_rf_bump.subs(t,t-1) 
u_rr_bump=S.Zero#road profile given in time domain
 #road profile given in time domain

u_rf_mul_bumps_on=sum(u_rf_bump_on.subs(x,x-ii*l_bumps) for ii in range(5))
u_rf_mul_bumps=0
u_rr_mul_bumps_on=u_rf_mul_bumps_on.subs(x,x-2*R)
u_rr_mul_bumps=0#road profile given in time domain
 #road profile given in time domain
    
    
# u_rf_ramp=u0/l_ramp*x
u_rf_ramp=S.Zero
# u_rr_ramp=u_rf_ramp.subs(x,x-2*R)
u_rr_ramp=S.Zero

u_road=Function('u_r')
# u_rf=u_road(x)
u_rf=S.Zero
# u_rr=u_road(x-2*R)
u_rr=S.Zero

# u_rr=u0*(sin(x-2*R))  #road profile given in time domain
u_rr=S.Zero
# u_rf=u_rr.subs(x,x+2*R) #road profile given in time domain
u_rf=S.Zero

#
#u_fr=Function('u_fr')(phi,z_fr)

# u_fr=u_fr_expr

# u_rear_expr=sqrt((z+R*sin(phi)+l_rear-z_rear)**2 + (R-R*cos(phi))**2)-l_rear
# #u_rear=Function('u_rear')(phi,z_fr)

# u_rear=u_rear_expr
# D=S.One/2*c*2*(dx**2) + c/2*sum(vel**2 for vel in u)




x_m3=x+z_c3*phi
x_RC=x+l_RC/2*cos(phi0+phi_rc)
z_RC=z+0.1*R*phi+l_RC/2*sin(phi0+phi_rc)


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
# u_rear=sqrt((z+R*sin(phi)+l_rear-z_rear)**2 + (R-R*cos(phi))**2)-l_rear#
u_rear=0
# u_fr=sqrt((z-R*sin(phi)+l_fr-z_fr)**2 + (R-R*cos(phi))**2)-l_fr
u_fr=0


# linear springs deflection
u_rear=((z+R*phi+l_rear -z_rear)**2 )-l_rear#
#u_rear=0
u_fr=  ((z-R*phi+l_fr   -z_fr)**2  )-l_fr
#u_fr=0



# spring potentials
V_rear = S.One/2*k_rt*(z_rear-u_rr)**2  + S.One/2*k_r*(u_rear)**2 # Ep tylnich prętów ramy względem ich sprężystości
V_fr = S.One/2*k_ft*(z_fr-u_rf)**2  + S.One/2*k_f*(u_fr)**2       # Ep przednich prętów ramy względem ich sprężystości



V_rc_g=-m_RC*g*z_RC
V_rc_fix =  S.One/2 * k_fix * phi_rc**2 + S.One/2 *k_tire * (z_wrc)**2 
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


# D_chair5dof=(c_mu*T_chair5dof).doit().expand()
D_chair5dof=(c_mu*T_chair5dof).doit().expand()

D_rapidchair3dof=(c_mu*T_RC_3dof).doit().expand()


D_chair5dof_rc3dof=(c_mu*T_chair5dof_rc3dof).doit().expand()

# D+= c_lam* S.One/2*k_r*(((z+R*(phi)-z_rear)**2 ).diff(t))**2

# D+= c_lam* S.One/2*k_f*(((z-R*(phi)-z_fr)**2 ).diff(t))**2


qs_5dof = x, z_rear, z_fr, z, phi
qs_rc_3dof=z_wrc,theta,phi_rc
qs_chair5dof_rc_3dof=qs_5dof+qs_rc_3dof




N = ReferenceFrame('N')
Pl = Point('P_l')
Pr = Point('P_r')
Pl.set_vel(N, (u[0]+u[1]*R*sin(Omega*t).diff(t)) * N.x+(u[2]*R*cos(Omega*t).diff(t)) * N.y)
Pr.set_vel(N, u[0] * N.x)




points_list=[Point('P_'+str(no))  for no,vel in enumerate(Matrix(qs_chair5dof_rc_3dof).diff(t))]
[points_list[no].set_vel(N,vel*N.x)  for no,vel in enumerate(Matrix(qs_chair5dof_rc_3dof).diff(t))]



# FL = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*F*(cos(Omega*t)) * N.x+sin(Omega*t) * N.y)]+[(points_list[no],-D_chair5dof.diff(vel)*N.x)  for no,vel in enumerate(u)]   #,(Pr, R*F*cos(Omega*t) * N.x)]
FL_chair5dof = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*2*F*N.x)]+[(points_list[no],-D_chair5dof.diff(vel)*N.x)  for no,vel in enumerate(Matrix(qs_5dof).diff(t))]


FL_rapidchair3dof = [(points_list[no],-D_rapidchair3dof.diff(vel)*N.x)  for no,vel in enumerate(Matrix(qs_rc_3dof).diff(t))]


FL_rapidchair3dof = [(points_list[no],-D_rapidchair3dof.diff(vel)*N.x)  for no,vel in enumerate(Matrix(qs_rc_3dof).diff(t))]


FL_chair5dof_rc3dof = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*2*F*N.x)]+[(points_list[no],-D_chair5dof_rc3dof.diff(vel)*N.x)  for no,vel in enumerate(Matrix(qs_chair5dof_rc_3dof).diff(t))]

FL_chair5dof_rc3dof = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*2*F*N.x)]+[(points_list[no],-D_chair5dof_rc3dof.diff(vel)*N.x)  for no,vel in enumerate(Matrix(qs_chair5dof_rc_3dof).diff(t))]




L_chair5dof = T_chair5dof - V_chair5dof
# L_5dof=

# L_chair5dof_rc3dof = T_chair5dof_rc3dof-V_chair5dof_rc3dof


# L_rc1dof=

L_default_5dof= L_chair5dof

L_rc_3dof=T_RC_3dof-V_RC_3dof

L_chair5dof_rc3dof=T_chair5dof_rc3dof-V_chair5dof_rc3dof





V_rear_road = S.One/2*k_rt*(z_rear-u_rr_bump_on)**2  + S.One/2*k_r*(u_rear)**2 # Ep tylnich prętów ramy względem ich sprężystości
V_fr_road = S.One/2*k_ft*(z_fr-u_rf_bump_on)**2  + S.One/2*k_f*(u_fr)**2       # Ep przednich prętów ramy względem ich sprężystości
T_chair5dof = T_body + T_rot + T_rear + T_fr + T_wheel
V_chair5dof_road = V_chair_g + V_rear_road + V_fr_road
L_chair5dof_withroad = T_chair5dof - V_chair5dof_road





class Chair5DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_default_5dof, qs=qs_5dof, forcelist=FL_chair5dof, bodies=None, frame=N,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t'),**kwargs):
        
        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar,**kwargs)
        

        
        
        
class RapidChair3DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_rc_3dof, qs=qs_rc_3dof, forcelist=FL_rapidchair3dof, bodies=None, frame=N,
                   hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t'),**kwargs):

        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
             hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar,**kwargs)

        
        
        
class Chair5DOFwithRC3DOF(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_chair5dof_rc3dof, qs=qs_chair5dof_rc_3dof, forcelist=FL_chair5dof_rc3dof, bodies=None, frame=N,
                   hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t'),**kwargs):

        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
             hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar,**kwargs)


chair_5dof = Chair5DOF()('Chair 5DOF model')
chair_5dof_lin=chair_5dof.linearized()('Chair linearized 5DOF model')
chair_3dof_rev = chair_5dof.subs(dof3_rev,method='direct').shranked(x, z_rear, z)('Chair 3DOF model with revolution')
chair_3dof_norev = chair_5dof.subs(dof3_norev,method='direct').shranked(x, z, phi)('Chair 3DOF model without revolution')       
chair_2dof = chair_5dof.subs(dof2,method='direct').shranked(x, z)('Chair 2DOF model')




rapidchair_3dof = RapidChair3DOF()('RapidChair 3DOF model')
rapidchair_2dof=rapidchair_3dof.subs(rcdof2,method='direct').shranked(theta,phi_rc)('RapidChair 2DOF model')
rapidchair_1dof=rapidchair_3dof.subs(rcdof1,method='direct').shranked([phi_rc])('RapidChair 1DOF model')





chair5dof_rc3dof=Chair5DOFwithRC3DOF()

### system parametes
m_fr, m_rear,m_3, k_r, k_rt, k_f, k_ft, k_rot = symbols('m_fr, m_r, M, k_r, k_rt, k_f, k_ft k_rot',positive=True)
m,k,g,F_1,F_2,Omega,F,R,v0,u0,l_l,l_r= symbols('m,k,g,F_1,F_2,Omega, F_0 R, v_0,u_0,l_l,l_r',positive=True)
I_ch, I_w , z_c3,l_fr,l_rear= symbols('I_chair, I_wheel, z_c3, l_fr, l_r',positive=True)
m_RC, I_RC, l_RC, k_RC, phi0 = symbols('m_RC, I_RC, l_RC, k_RC, varphi_0',positive=True)
m_w, I_wrc, r_w, k_w, k_fix,k_tire = symbols('m_w, I_wRC, R_w, k_w, k_fix, k_t',positive=True)
c,c_mu,c_lam= symbols('c,c_mu, c_lambda',positive=True)
a_ox, a_oz, a_rz ,a_rcz=symbols('a_ox, a_oz, a_rz ,a_rcz')

pm=symbols('PM')

t_l, delta_t,l_ramp,l_bumps = symbols('t_l,dt,l_ramp, l_bumps',positive=True)
t0 = symbols('t_0')

# generalized coordinates
x,z_fr,z_rear,z, phi,phi_rc,theta,z_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc')
dx,dz_fr,dz_rear,dz,dphi,dphi_rc,dtheta,dz_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc', 1)
# ===========================SIMPLE MODEL=========================================================================================================
T_body = S.One/2*m_3*x_m3.diff(t)**2 + S.One/2 * m_3* dz**2 # Ek ramy wózka w ruchu postępowym
T_rot =  S.One/2 * I_ch* dphi**2  # Ek ramy wózka w ruchu obrotowym
T_rear = S.One/2 * m_rear* dz_rear**2 + S.One/2 * m_rear * dx**2    # Ek tylniego koła Ir (mr)
T_fr = S.One/2 * m_fr* dz_fr**2 + S.One/2 * m_fr * dx**2    # Ek przedniego koła If (mf)
T_wheel = S.One/2 * I_w * (dx/R)**2

V_rear = S.One/2*k_rt*(z_rear-amplitude*cos(2*pi/(length/speed)*t))**2# Ep tylnich prętów ramy względem ich sprężystości
V_fr = S.One/2*k_ft*(z_fr-amplitude*cos(2*pi/(length/speed)*(t-axw/speed)))**2     # Ep przednich prętów ramy względem ich sprężystości


#+++++++++++++++++++++
V_pion= S.One/2 * k_r*(z-z_rear)**2 + S.One/2 * k_f*(z-z_fr)**2  ##### te
#+++++++++++++++++++++

##                k_rear                                    k_fr
V_obr = S.One/2 * k_r * (z+l_l*phi-z_rear)**2 + S.One/2 * k_f * (z-l_r*phi-z_fr)**2 - S.One/2 * m_3 * z_c3 * g *phi**2

Ekinetic=T_body+T_rot+T_rear+T_fr+T_wheel
Epotential=V_rear+V_fr+V_pion+V_obr

L_simple=Ekinetic-Epotential


Dlam_chair5dof=((chair_5dof_lin.q.diff(t)).T*Matrix([Epotential]).jacobian(chair_5dof_lin.q).jacobian(chair_5dof_lin.q)*chair_5dof_lin.q.diff(t)*c_lam)[0]
Dmu_chair5dof=((chair_5dof_lin.q.diff(t)).T*Matrix([Ekinetic]).jacobian(chair_5dof_lin.q.diff(t)).jacobian(chair_5dof_lin.q.diff(t))*chair_5dof_lin.q.diff(t)*c_mu)[0]


FL_SIMPLECHAIR = [(Pl, 0.5*(1*sign(cos(Omega*t)-pm  )+1)*2*F*N.x)]+[(points_list[nom],(-Dmu_chair5dof.diff(velm)-Dlam_chair5dof.diff(velm))*N.x)  for nom,velm in enumerate(Matrix(qs_5dof).diff(t))]
# +[(points_list[non],-Dlam_chair5dof.diff(veln)*N.x)  for non,veln in enumerate(Matrix(qs_5dof).diff(t))]
class SimpleChair5DOF(dyn.HarmonicOscillator):
    
    m_fr, m_rear,m_3, k_r, k_rt, k_f, k_ft, k_rot = symbols('m_fr, m_r, M, k_r, k_rt, k_f, k_ft k_rot',positive=True)
    m,k,g,F_1,F_2,Omega,F,R,v0,u0,l_l,l_r= symbols('m,k,g,F_1,F_2,Omega, F_0 R, v_0,u_0,l_l,l_r',positive=True)
    I_ch, I_w , z_c3,l_fr,l_rear= symbols('I_chair, I_wheel, z_c3, l_fr, l_r',positive=True)
    m_RC, I_RC, l_RC, k_RC, phi0 = symbols('m_RC, I_RC, l_RC, k_RC, varphi_0',positive=True)
    m_w, I_wrc, r_w, k_w, k_fix,k_tire = symbols('m_w, I_wRC, R_w, k_w, k_fix, k_t',positive=True)
    c,c_mu,c_lam= symbols('c,c_mu, c_lambda',positive=True)
    a_ox, a_oz, a_rz ,a_rcz=symbols('a_ox, a_oz, a_rz ,a_rcz')

    pm=symbols('PM')

    t_l, delta_t,l_ramp,l_bumps = symbols('t_l,dt,l_ramp, l_bumps',positive=True)
    t0 = symbols('t_0')
    
    # generalized coordinates
    x,z_fr,z_rear,z, phi,phi_rc,theta,z_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc')
    dx,dz_fr,dz_rear,dz,dphi,dphi_rc,dtheta,dz_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc', 1)
    
    def __init__(self,
                 m_fr=m_fr,
                 m_rear=m_rear,
                 m_3=m_3,
                 k_r=k_r,
                 k_rt=k_rt,
                 k_f=k_f,
                 k_ft=k_ft,
                 k_rot=k_rot,
                 qs=[x,z_fr,z_rear,z, phi],
                 forcelist=FL_SIMPLECHAIR,
                 bodies=None, frame=N,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t'),**kwargs):
        

   
        
        
        if m_3 is None:
            m_3 = self.__class__.m_3

        self.m_3 = m_3
        
        m,k,g,F_1,F_2,Omega,F,R,v0,u0,l_l,l_r= symbols('m,k,g,F_1,F_2,Omega, F_0 R, v_0,u_0,l_l,l_r',positive=True)
        I_ch, I_w , z_c3,l_fr,l_rear= symbols('I_chair, I_wheel, z_c3, l_fr, l_r',positive=True)
        m_RC, I_RC, l_RC, k_RC, phi0 = symbols('m_RC, I_RC, l_RC, k_RC, varphi_0',positive=True)
        m_w, I_wrc, r_w, k_w, k_fix,k_tire = symbols('m_w, I_wRC, R_w, k_w, k_fix, k_t',positive=True)
        c,c_mu,c_lam= symbols('c,c_mu, c_lambda',positive=True)
        a_ox, a_oz, a_rz ,a_rcz=symbols('a_ox, a_oz, a_rz ,a_rcz')

        pm=symbols('PM')

        t_l, delta_t,l_ramp,l_bumps = symbols('t_l,dt,l_ramp, l_bumps',positive=True)
        t0 = symbols('t_0')
        
        
        # generalized coordinates
        x,z_fr,z_rear,z, phi,phi_rc,theta,z_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc')
        dx,dz_fr,dz_rear,dz,dphi,dphi_rc,dtheta,dz_wrc = dynamicsymbols('x,z_fr,z_r,z, varphi,varphi_RC, theta, z_wrc', 1)
        
        # ===========================SIMPLE MODEL=========================================================================================================
        T_body = S.One/2*m_3*x_m3.diff(t)**2 + S.One/2 * m_3* dz**2 # Ek ramy wózka w ruchu postępowym
        T_rot =  S.One/2 * I_ch* dphi**2  # Ek ramy wózka w ruchu obrotowym
        T_rear = S.One/2 * m_rear* dz_rear**2 + S.One/2 * m_rear * dx**2    # Ek tylniego koła Ir (mr)
        T_fr = S.One/2 * m_fr* dz_fr**2 + S.One/2 * m_fr * dx**2    # Ek przedniego koła If (mf)
        T_wheel = S.One/2 * I_w * (dx/R)**2

        V_rear = S.One/2*k_rt*(z_rear-amplitude*cos(2*pi/(length/speed)*t))**2# Ep tylnich prętów ramy względem ich sprężystości
        V_fr = S.One/2*k_ft*(z_fr-amplitude*cos(2*pi/(length/speed)*(t-axw/speed)))**2     # Ep przednich prętów ramy względem ich sprężystości


        #+++++++++++++++++++++
        V_pion= S.One/2 * k_r*(z-z_rear)**2 + S.One/2 * k_f*(z-z_fr)**2  ##### te
        #+++++++++++++++++++++

        ##                k_rear                                    k_fr
        V_obr = S.One/2 * k_r * (z+l_l*phi-z_rear)**2 + S.One/2 * k_f * (z-l_r*phi-z_fr)**2 - S.One/2 * m_3 * z_c3 * g *phi**2

        Ekinetic=T_body+T_rot+T_rear+T_fr+T_wheel
        Epotential=V_rear+V_fr+V_pion+V_obr

        L_simple=Ekinetic-Epotential
        
        
        super().__init__( Lagrangian=L_simple, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar,**kwargs)

    def get_param_values(self):
        default_data_dict={F:150,
                   
                   c_mu:0.0001,
                   c_lam:0.0001,
                   l_l:0.2,
                   l_r:0.4,
                   
                   k_f:607500,
                   k_ft:475000,
                   k_r:580000,
                   k_rt:400000,
                   m_3:75,
                   I_ch:9.479342+m_3*0.39*0.39,
                   m_rear:1.5,
                   m_fr:0.6,
                   pm:0.1,
                   Omega:0.055*np.pi*2,
                   R:0.3,
                   z_c3:0.4,
                   g:9.81,
                   I_w:m_rear*R**2,
                   l_fr:0.2,
                   l_rear:0.01,
                   u0:0.005,

                   l_bumps:0.15,
                   amplitude:0.0165,
                   length:0.19,
                   speed:1.7,
                   axw:0.47}
        
        return default_data_dict       
    
    def get_table_values(self):
        table_data_dict={F:150,
                   
                   c_mu:0.0001,
                   c_lam:0.0001,
                   l_l:0.2,
                   l_r:0.4,
                   
                   k_f:607500,
                   k_ft:475000,
                   k_r:580000,
                   k_rt:400000,
                   m_3:75,
                   I_ch:20.8868,
                   m_rear:1.5,
                   m_fr:0.6,
                   pm:0.1,
                   Omega:0.3454,
                   R:0.3,
                   z_c3:0.4,
                   g:9.81,
                   I_w:0.135,
                   l_fr:0.2,
                   l_rear:0.01,
                   u0:0.005,

                   l_bumps:0.15,
                   amplitude:0.0165,
                   length:0.19,
                   speed:1.7,
                   axw:0.47}
        return table_data_dict       
    def numerical_model(self):
        pass

        
# ===========================SIMPLE MODEL=========================================================================================================
chair_5dof_simple = SimpleChair5DOF()('Simple chair 5DOF model')






ureg.define('dimles = 0 = - ')
ureg.dimles

# units_dict={
#             c:ureg.kilogram/ureg.second,
#             c_mu:S.One/ureg.second,
#             c_lam:S.One/ureg.second,
#             m_3:ureg.kilogram,    
#             m_w:ureg.kilogram,
#             m_fr:ureg.kilogram,
#             m_rear:ureg.kilogram,
#             m:ureg.kilogram,
#             m_RC:ureg.kilogram,
#             k_r:ureg.newton/ureg.meter,
#             k_rt:ureg.newton/ureg.meter,
#             k_f:ureg.newton/ureg.meter,
#             k_w:ureg.newton/ureg.meter,
#             k_fix:ureg.newton/ureg.meter,
#             k_tire:ureg.newton/ureg.meter,
#             k_ft:ureg.newton/ureg.meter,
#             k:ureg.newton/ureg.meter,
#             k_RC:ureg.newton/ureg.meter,
#             g:ureg.meter/ureg.second**2,
#             F_1:ureg.newton,
#             F_2:ureg.newton,
#             F:ureg.newton,
#             u0:ureg.meter,
#             Omega:ureg.radian/ureg.second,
#             R:ureg.meter,
#             R_curve:ureg.meter,
#             r_w:ureg.meter,
#             v0:ureg.meter/ureg.second,
#             I_ch:ureg.kilogram*ureg.meter**2,
#             I_w:ureg.kilogram*ureg.meter**2,
#             I_RC:ureg.kilogram*ureg.meter**2,
#             I_wrc:ureg.kilogram*ureg.meter**2,
#             dphi:ureg.radian/ureg.second,
#             dphi_rc:ureg.radian/ureg.second,
#             z_c3:ureg.meter,
#             z_fr:ureg.meter,
#             l_fr:ureg.meter,
#             l_rear:ureg.meter,
#             l_RC:ureg.meter,
#             phi0:ureg.radian,
#             x:ureg.meter,
#             y:ureg.meter,
#             z:ureg.meter,
#             z_rear:ureg.meter,
#             z_fr:ureg.meter,
#             z_wrc:ureg.meter,
#             phi:ureg.radian,
#             alpha:ureg.radian,
#             phi_rc:ureg.radian,
#             theta:ureg.radian,
#             dx:ureg.meter/ureg.second,
#             dy:ureg.meter/ureg.second,
#             dz_fr:ureg.meter/ureg.second,
#             dz_rear:ureg.meter/ureg.second,
#             dz:ureg.meter/ureg.second,
#             dtheta:ureg.radian/ureg.second,
#             dalpha:ureg.radian/ureg.second,
#             dz_wrc:ureg.meter/ureg.second,
#             pm:ureg.dimles,
#             t_l:ureg.second,
#             delta_t:ureg.second,
#             l_ramp:ureg.meter,
#             l_bumps:ureg.meter,
#             t0:ureg.second,
#             'dimensionless':ureg.meter/ureg.meter,
# #             ao_x_max:ureg.gram,
#             a_ox:ureg.gram,
# #             ao_x_min:ureg.gram,
# #             ao_x_idmax:ureg.second,
# #             ao_x_idmin:ureg.second,
# #             ao_z_max:ureg.gram,
#             a_oz:ureg.gram,
# #             ao_z_min:ureg.gram,
# #             ao_z_idmax:ureg.second,
# #             ao_z_idmin:ureg.second,
# #             ao_rx_max:ureg.gram,
# #             ao_rx_min:ureg.gram,
# #             ao_rx_idmax:ureg.second,
# #             ao_rx_idmin:ureg.second,
# #             ao_rz_max:ureg.gram,
# #             ao_rz_min:ureg.gram,
# #             ao_rz_idmax:ureg.second,
# #             ao_rz_idmin:ureg.second,
#             a_rz:ureg.gram,
#             a_rcz:ureg.gram,
#             a_sim:ureg.meter/ureg.second/ureg.second,
#             dx.diff(t):ureg.meter/ureg.second/ureg.second,
#             dy.diff(t):ureg.meter/ureg.second/ureg.second,
#             dz_fr.diff(t):ureg.meter/ureg.second/ureg.second,
#             dz_rear.diff(t):ureg.meter/ureg.second/ureg.second,
#             dz.diff(t):ureg.meter/ureg.second/ureg.second,
#             dtheta.diff(t):ureg.radian/ureg.second/ureg.second,
#             dphi.diff(t):ureg.radian/ureg.second/ureg.second,
#             dz_wrc.diff(t):ureg.meter/ureg.second/ureg.second,
#             l_l:ureg.meter,
#             l_r:ureg.meter,
#             amplitude:ureg.meter,
#             length:ureg.meter,
#             speed:ureg.meter/ureg.second,
#             axw:ureg.meter,
#             ivar:ureg.second
#            }



# model_names_dict={
#     chair_2dof:'chair_model_2dof',
#     chair_3dof_norev:'chair_model_3dof_norev',
#     chair_3dof_rev:'chair_model_3dof_rev',
#     chair_5dof:'chair_model_5dof',
#     chair_5dof_lin:'chair_model_5dof',
#     rapidchair_1dof:'rapidchair_model_1dof',
#     rapidchair_2dof:'rapidchair_model_2dof',
#     rapidchair_3dof:'rapidchair_model_3dof',
#     chair5dof_rc3dof:'chair_and_rapidchair_model'
#     }





# dof_names_dict={chair_2dof:'o dwóch stopniach swobody',
#                 chair_3dof_norev:'o trzech stopniach swobody bez obrotu',
#                 chair_3dof_rev:'o trzech stopniach swobody z obrotem',
#                 chair_5dof:'o pięciu stopniach swobody (nieliniowy)',
#                 chair_5dof_lin:'o pięciu stopniach swobody (zlinearyzowany)',
#                 rapidchair_1dof:'o jednym stopniu swobody',
#                 rapidchair_2dof:'o dwóch stopniach swobody',
#                 rapidchair_3dof:'o trzech stopniach swobody',
#                 chair5dof_rc3dof:'o ośmiu stopniach swobody'
#                }


