from . import en

class LagrangianMCA(en.LagrangianMCA):
    question = 'Lagrangian dla rozważanego układu można wyrazić następującym równaniem:'
    
class SDoFGoverningEquationMCA(en.SDoFGoverningEquationMCA):
    question = 'Równanie ruchu dla układu przedstawionego na rysunku wyraża następujący wzór:'
    
class SDoFApproximatedGoverningEquationMCA(en.SDoFApproximatedGoverningEquationMCA):
    question = 'Nieliniowe przybliżone równanie ruchu opisuje następujący wzór:'
    
class SmallParameterMCA(en.SmallParameterMCA):
    question = 'Mały parametr układu wyraża się wzorem:'
    
class ResonanceCurveMCA(en.ResonanceCurveMCA):
    question = 'Zależność pomiędzy amplitudą a częstością drgań własnych dla drgań swobodnych rozważanego układu wyraża wzór:'