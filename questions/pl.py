from . import en

class LagrangianMCA(en.LagrangianMCA):
    question = 'Lagrangian dla rozważanego układu można wyrazić następującym równaniem:'
    
class SDoFGoverningEquationMCA(en.SDoFGoverningEquationMCA):
    question = 'Równanie ruchu dla układu przedstawionego na rysunku wyraża następujący wzór:'
    
class GoverningEquationMCA(en.GoverningEquationMCA):
    question = 'Równania dynamiki dla układu przedstawionego na rysunku wyrazić można następującym układem równań:'
    
class FundamentalMatrixMCA(en.FundamentalMatrixMCA):
    question = 'Wyznacz macierz fundamentalną układu:'
    
class OmegaMCA(en.OmegaMCA):
    question = 'Wyznacz częstości drgań własnych układu:'
    
class SDoFApproximatedGoverningEquationMCA(en.SDoFApproximatedGoverningEquationMCA):
    question = 'Nieliniowe przybliżone równanie ruchu opisuje następujący wzór:'
    
class SmallParameterMCA(en.SmallParameterMCA):
    question = 'Mały parametr układu wyraża się wzorem:'
    
class ResonanceCurveMCA(en.ResonanceCurveMCA):
    question = 'Zależność pomiędzy amplitudą a częstością drgań własnych dla drgań swobodnych rozważanego układu wyraża wzór:'

class SDoFFreeGalerkinGoverningEquationIntegralMCA(en.SDoFFreeGalerkinGoverningEquationIntegralMCA):
    question = 'Wyrażenie przedstawiające szereg pozwalajacy na określenie częstości własnych badanego układu to:'
    
class SDoFFreeGalerkinGoverningEquationIntegralTR8MCA(en.SDoFFreeGalerkinGoverningEquationIntegralTR8MCA):
    question = 'Odpowiednie tożsamości trygonometryczne zostały zastosowane w równaniu:'

class SmallParameterNMSA(en.SmallParameterNMSA):
    question = 'Określ wartość małego parametru dla następujących danych \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'
    
class FrequencyForAmplitudeNMSA(en.FrequencyForAmplitudeNMSA):
    question = 'Wyznacz częstotliwość nielinowych drgań dla amplitudy \(a=0.1 m\), jeśli \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

class NonLinLowerAmplitudeValueNMSA(en.NonLinLowerAmplitudeValueNMSA):
    question = 'Wyznacz mniejszą amplitudę nielinowych drgań dla częstości \(\omega=3.5 \\frac{rad}{s}\), jeśli \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'
    
class NonLinHigherAmplitudeValueNMSA(en.NonLinHigherAmplitudeValueNMSA):
    question = 'Wyznacz większą amplitudę nielinowych drgań dla częstości \(\omega=3.5 \\frac{rad}{s}\) , jeśli \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

class SystemResponseAmplitudeC2MCA(en.SystemResponseAmplitudeC2MCA):
    question = 'Wyznacz amplitudę C2 drgań ustalonych rozważanego układu:'
    
class LowerNaturalFreqValueNMSA(en.LowerNaturalFreqValueNMSA):
    question = 'Wyznacz wartość niższej częstości drgań własnych układu dla zadanych parametrów: \(m=150 kg\), \(k=900 \\frac{N}{m}\) (wartość zaokrąglić do drugiego miejsca po przecinku):'

class HigherNaturalFreqValueNMSA(en.HigherNaturalFreqValueNMSA):
    question = 'Wyznacz wartość wyższej częstości drgań własnych układu dla zadanych parametrów: \(m=150 kg\), \(k=900 \\frac{N}{m}\) (wartość zaokrąglić do drugiego miejsca po przecinku):'

class SteadyAmplitudesRatioNMSA(en.SteadyAmplitudesRatioNMSA):
    question = 'Ile razy amplituda drgań ustalonych pierwszego wózka jest większa od amplitudy drgań ustalonych drugiego wózka przy zadanych parametrach: \(m=150 kg\), \(k=900 \\frac{N}{m}\), \(\\omega=15 \\frac{rad}{s}\):'

class Test(en.Test):
    question = 'Test:'
    
class SecondModeFirstComponentValueNMSA(en.SecondModeFirstComponentValueNMSA):
    question = 'Wyznacz wartość pierwszej składowej drugiej postaci drgań dla zadanych parametrów \(m=150 kg\), \(k=900 \\frac{N}{m}\):'
    