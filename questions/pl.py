# import sympy as sym
# from sympy import *
# import sympy.physics.mechanics as mech
# from .dynamics import *
from .en import *

class LagrangianMCA(LagrangianMCA):
    question = 'Lagrangian dla rozważanego układu można wyrazić następującym równaniem:'
    
class SDoFGoverningEquationMCA(SDoFGoverningEquationMCA):
    question = 'Równanie ruchu dla układu przedstawionego na rysunku wyraża następujący wzór:'
    
class SDoFApproximatedGoverningEquationMCA(SDoFApproximatedGoverningEquationMCA):
    question = 'Nieliniowe przybliżone równanie ruchu opisuje następujący wzór:'
    
class SmallParameterMCA(SmallParameterMCA):
    question = 'Mały parametr układu wyraża się wzorem:'
    
class ResonanceCurveMCA(ResonanceCurveMCA):
    question = 'Zależność pomiędzy amplitudą a częstością drgań własnych dla drgań swobodnych rozważanego układu wyraża wzór:'