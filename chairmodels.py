import importlib
import tools as dyn
importlib.reload(dyn)
from sympy.physics.vector import dynamicsymbols
from sympy.physics.mechanics import *
import sympy as sym
from sympy import *
mechanics_printing()

coords='x, z_r, z_f, z, phi'
phi, x, z_r, z_f, z = dynamicsymbols('phi, x, z_r, z_f, z')
m_3,m_r,m_f,k_f,k_r,R,I_r,I_f,g,r,t,k_r_t, k_f_t, z_C= symbols('m_3, m_r, m_f, k_f, k_r, R, I_r, I_f, g, r, t, k_r_t, k_f_t, z_C', positive=True)


u_front=(x-phi*R)**2
v_front=u_front.diff(t)


T1 = S.One/2*m_3*diff(x)**2    # Ek wózka po x
T2 = S.One/2*I_r*diff(phi)**2  # Ek tylniego koła Ir (mr)
T3 = S.One/2*I_f*diff(phi)**2  # Ek przedniego koła If (mf)


x_r = (x - ((sin(phi)+cos(phi))**2)*R) # geometria wzgledem przemieszczenia - tył
x_f = (x + ((sin(phi)+cos(phi))**2)*r) # geometria względem przemieszczeni - przód


v_r = (x - ((sin(phi)+cos(phi))**2)*R).diff(x) + (x - ((sin(phi)+cos(phi))**2)*R).diff(phi) # geometria wzgledem predkosci - tył
v_f = (x + ((sin(phi)+cos(phi))**2)*r).diff(x) + (x + ((sin(phi)+cos(phi))**2)*r).diff(phi) # geometria względem predkosci - przód


T4 = S.One/2*m_r*(v_r**2) # Ek bezwładności tyłu wózka wzlędem całkowitego przemieszczenia wózka
T5 = S.One/2*m_f*(v_f**2) # Ek bezwładności przodu wózka wzlędem całkowitego przemieszczenia wózka


V1 = m_3*g*z*sin(phi) # EPG wózka względem punktu cięzkości wózka
V2 = (S.One/2)*k_r*(x_r-z_r)**2+(S.One/2)*k_r_t*z_r**2 # Ep tylnich prętów ramy względem ich sprężystości
V3 = (S.One/2)*k_f*(x_f-z_f)**2+(S.One/2)*k_f_t*z_f**2 # Ep przednich prętów ramy względem ich sprężystości


T_5dof = T1 + T2 + T3 + T4 + T5
V_5dof = V1 + V2 + V3

dof3 = {z_r:z,z_f:z} # 3 DoF
dof2 = {z_r:z,z_f:z,x:phi*z_C} # 2 DoF 

T_3dof = T_5dof
V_3dof = V1 + V2.subs(dof3) + V3.subs(dof3)

T_2dof = T1.subs(dof2) + T2 + T3 + T4.subs(dof2) + T5.subs(dof2)
V_2dof = V1 + V2.subs(dof2) + V3.subs(dof2)


L_5dof = T_5dof - V_5dof
L_3dof = T_3dof - V_3dof
L_2dof = T_2dof - V_2dof
L_5dof, L_3dof, L_2dof

L_default_5dof= L_5dof
L_default_3dof= L_3dof
L_default_2dof= L_2dof

qs_5dof = phi, x, z_r, z_f, z
qs_3dof = phi, x, z
qs_2dof = phi, z


class Chair5dof(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_default_5dof, qs=qs_5dof, forcelist=None, bodies=None, frame=None,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)

class Chair3dof(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_default_3dof, qs=qs_3dof, forcelist=None, bodies=None, frame=None,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)

        
class Chair2dof(dyn.LagrangesDynamicSystem):
    def __init__(self, Lagrangian=L_default_2dof, qs=qs_2dof, forcelist=None, bodies=None, frame=None,
                       hol_coneqs=None, nonhol_coneqs=None,label=None,ivar=sym.Symbol('t')):
        
        super().__init__( Lagrangian=Lagrangian, qs=qs, forcelist=forcelist, bodies=bodies, frame=frame,
                 hol_coneqs=hol_coneqs, nonhol_coneqs=nonhol_coneqs,label=label,ivar=ivar)
        


