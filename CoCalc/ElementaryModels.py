from sympy import S, Symbol, diff
from dynpy import LagrangesDynamicSystem
from sympy.physics.mechanics import *

class MaterialPoint(LagrangesDynamicSystem):
    def __init__(self, m, pos1, pos2=0, qs=None, ivar=Symbol('t')):
        
        if qs == None:
            qs = [pos1]
        else:
            qs = [pos1,pos2]

            if pos2.diff(ivar) is False:
                print(pos2==0)
            else:
                print(pos2.diff(ivar))
        
        Lagrangian = S.Half * m * (pos1-pos2).diff(ivar)**2
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Spring(LagrangesDynamicSystem):
    def __init__(self, k, pos1, pos2=0, qs=None, ivar=Symbol('t')):
        
        if qs == qs:
            qs = [pos1]
        else:
            qs = [pos1,pos2] # array list type
        
        Lagrangian = -((S.Half * (k * (pos1 - pos2)**2)))
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Disk(LagrangesDynamicSystem):
    def __init__(self, I, pos1, ivar=Symbol('t')):
        
        qs=[pos1]
        
        Lagrangian = (S.Half * (I * (pos1.diff(ivar))**2))
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class RigidBody2D(LagrangesDynamicSystem):
    def __init__(self, m, I, pos_lin, pos_rot, ivar=Symbol('t')):
        
        qs = [pos_lin, pos_rot] # array list type
        
        Lagrangian = S.Half * m * pos_lin.diff(ivar)**2 + S.One / 2 * I * pos_rot.diff(ivar)**2      
        super().__init__(Lagrangian=Lagrangian, qs=qs, ivar=ivar)

        
class Damper(LagrangesDynamicSystem):
    def __init__(self, c, pos1, pos2=0, qs=None, ivar=Symbol('t'), frame=ReferenceFrame('N')):
        
        if qs == qs:
            qs = [pos1]
        else:
            qs = [pos1, pos2] # array list type
        
        dpos1 = diff(pos1, ivar)
        dpos2 = diff(pos2, ivar)
        
        P = Point('P')
        P.set_vel(frame, (dpos1 - dpos2) * frame.x)
        
        D = (((S.Half) * c * (dpos1 - dpos2)**2).diff(dpos1))
        
        forcelist = [(P, -D*frame.x)]     
        super().__init__(0, qs=qs, forcelist=forcelist, frame=frame, ivar=ivar)