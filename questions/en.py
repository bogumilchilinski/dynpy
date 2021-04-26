import sympy as sym
from sympy import *
import sympy.physics.mechanics as mech
from sympy.simplify.fu import TR8

from dynpy import HarmonicOscillator

from ..moodle import *
import base64

t = Symbol('t')

########################### MDOF ########################


class GoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj:
                 (Eq(obj._eoms.doit().expand().applyfunc(simplify).expand(),
                     Matrix([0] * len(obj.q)),
                     evaluate=False)),
                 **kwargs):

        self.title = 'Podaj równiania ruchu układu:'
        self.title = 'Determine the system of equations of motion for the analysed system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFGoverningEquationMCA(MechanicalSystemAnswer):
    question = 'Determine equation of motion of the system:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj:
                 (Eq(obj._eoms.doit().expand().applyfunc(simplify).expand()[0],
                     0,
                     evaluate=False)),
                 **kwargs):

        #         self.title = 'Równanie ruchu dla układu przedstawionego na rysunku wyraża następujący wzór:'
        #         self.title = 'Determine equation of motion of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SDoFLinearizedGoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (Eq(obj.linearized()._eoms.doit(
                 ).expand().applyfunc(simplify).expand()[0],
                                                  0,
                                                  evaluate=False)),
                 **kwargs):

        self.title = 'Małe drgania układu opisuje równanie:'
        #         self.title = 'Determine equation of motion of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFApproximatedGoverningEquationMCA(MechanicalSystemAnswer):
    question = 'Determine equation of motion of the system:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(HarmonicOscillator(
                     obj.approximated())._eoms[0].doit().expand(),
                                                 0,
                                                 evaluate=False),
                 **kwargs):

        #         self.title = 'Nieliniowe przybliżone równanie ruchu opisuje nastęujący wzór:'
        #self.title = 'Determine equation of motion of the system:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SDoFFreeGalerkinGoverningEquationIntegralMCA(MechanicalSystemAnswer):
    question = 'Determine equation of motion of the system:'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(
                Integral((
                    (HarmonicOscillator(obj.approximated())._eoms[0].subs(
                        obj.q[0],
                        Symbol('a', positive=True) * sin(Symbol('omega') * t)).
                     doit().expand())).ratsimp() * sin(Symbol('omega') * t),
                         (t, 0, Symbol('T'))), 0),
            **kwargs):

        #         self.title = 'Nieliniowe przybliżone równanie ruchu opisuje nastęujący wzór:'
        #self.title = 'Determine equation of motion of the system:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SDoFFreeGalerkinGoverningEquationIntegralTR8MCA(MechanicalSystemAnswer):
    question = 'Determine equation of motion of the system:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Integral(((TR8(
                         (HarmonicOscillator(obj.approximated())._eoms[0].subs(
                             obj.q[0],
                             Symbol('a', positive=True) * sin(
                                 Symbol('omega') * t)).doit().expand())
                     )).ratsimp().doit().expand()) * sin(Symbol('omega') * t),
                              (t, 0, Symbol('T'))), 0),
                 **kwargs):

        #         self.title = 'Nieliniowe przybliżone równanie ruchu opisuje nastęujący wzór:'
        #self.title = 'Determine equation of motion of the system:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class CriticalPointsMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Matrix(obj.critical_points()[0:2]),  #
            **kwargs):

        self.title = 'Określ punkty równowagi rozważanego układu:'
        self.title = 'Determine the equilibrium points of the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class LinearizedGoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj:
                 (Eq(HarmonicOscillator(obj.linearized())._eoms.doit(),
                     Matrix([0] * len(obj.q)),
                     evaluate=False)),
                 **kwargs):

        self.title = 'Liniowe równania ruchu dla układu przedstawionego na rysunku można wyrazić następującym układem równań:'
        self.title = 'Linear equations of motion for the system presented might be expressed as the following system of equations:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class LagrangianMCA(MechanicalSystemAnswer):
    question = 'Choose the relation which describes the Lagrangian for the considered system:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('L'),
                                                 obj.lagrangian().doit()),
                 **kwargs):

        #         self.title = 'Lagrangian dla rozważanego układu można wyrazić następującym równaniem:'
        #         self.title = 'Choose corect dependance which determines Lagrangian of the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class ExternalForcesMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('F'), obj.external_forces(), evaluate=False),
                 **kwargs):

        self.title = 'Określ wektor sił w układzie:'
        self.title = 'Specify a forcing vector of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class LinearizedLagrangianMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(
                Symbol('L'),
                HarmonicOscillator(obj.linearized()).lagrangian()
            ),  # op_point to remove | True or False - doesn't matter
            **kwargs):

        self.title = 'Lagrangian dla małych drgań układu można wyrazić następującym równaniem:'
        self.title = 'Lagrangian for small vibrations of the system can be expressed as:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class OmegaMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [
                     Eq(Symbol('omega_0'), eig_val)
                     for eig_val in HarmonicOscillator(obj.linearized()).
                     natural_frequencies().doit().expand().doit()
                     if eig_val != 0
                 ],
                 **kwargs):

        self.title = 'Określ częstości drgań swobodnych występujących w układzie:'
        self.title = 'Determine the frequencies of free vibrations occuring in the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFOmegaMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [(
                     eig_val) for eig_val in HarmonicOscillator(obj.linearized(
                     )).natural_frequencies() if eig_val != 0][0],
                 **kwargs):

        self.title = 'Określ częstości drgań swobodnych występujących w układzie:'
        self.title = 'Determine the frequency of the free vibrations occuring in the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class FirstModeMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: HarmonicOscillator(
                     obj.linearized()).modes()[:, 0].n(3),
                 **kwargs):

        self.title = 'Określ pierwszą postać drgań układu:'
        self.title = 'Provide a first mode of a system'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SecondModeMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: HarmonicOscillator(
                     obj.linearized()).modes()[:, 1].n(3),
                 **kwargs):

        self.title = 'Określ drugą postać drgań układu:'
        self.title = 'Provide a second mode of a system'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class GeneralSolutionMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('X'),
                                            HarmonicOscillator(obj.linearized(
                                            )).general_solution().n(3),
                                            evaluate=False),
            **kwargs):

        self.title = 'Wyznacz rozwiązanie ogólne dla rozważanego układu:'
        self.title = 'Determine a general solution of ODEs for the investigated system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFGeneralSolutionMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(obj.q[0],
                                            HarmonicOscillator(obj.linearized(
                                            )).general_solution().n(3)[0],
                                            evaluate=False),
            **kwargs):

        self.title = 'Wyznacz rozwiązanie ogólne dla rozważanego układu:'
        self.title = 'Determine a general solution of ODEs for the investigated system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SteadySolutionMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('X_s'),
                                            HarmonicOscillator(obj.linearized(
                                            )).steady_solution().n(3),
                                            evaluate=False),
            **kwargs):

        self.title = 'Wyznacz rozwiązanie szczególne dla rozważanego układu:'
        self.title = 'Determine a particular solution of ODEs for the system under investigation:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFSteadySolutionMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(obj.q[0],
                                            HarmonicOscillator(obj.linearized(
                                            )).steady_solution().n(3)[0],
                                            evaluate=False),
            **kwargs):

        self.title = 'Wyznacz rozwiązanie szczególne dla rozważanego układu:'
        self.title = 'Determine a particular solution of ODEs for the system under investigation:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class InertiaMatrixMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('M'), (
                     HarmonicOscillator(obj.linearized()).inertia_matrix()),
                                                 evaluate=False),
                 **kwargs):
        self.title = 'Określ macierz bezwładności układu:'
        self.title = 'For the system under consideration determine the inertia matrix:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class StiffnessMatrixMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('K'), (
                     HarmonicOscillator(obj.linearized()).stiffness_matrix()),
                                                 evaluate=False),
                 **kwargs):
        self.title = 'Określ macierz sztywności układu:'
        self.title = 'For the system under consideration determine the stiffness matrix:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class DampingMatrixMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('C'), (
                     HarmonicOscillator(obj.linearized()).damping_matrix()),
                                                 evaluate=False),
                 **kwargs):
        self.title = 'Określ macierz tłumienia układu:'
        self.title = 'For the system under consideration determine the damping matrix:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class LeftSpringForceMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj:
        Eq(Symbol('F_sl'),
           ((HarmonicOscillator(obj.linearized()).stiffness_matrix()[0] +
             HarmonicOscillator(obj.linearized()).stiffness_matrix()[1])) *
           HarmonicOscillator(obj.linearized()).general_solution().n(3)[0],
           evaluate=False),
            **kwargs):
        self.title = 'Podaj wartość siły dynamicznej w lewej sprężynie:'
        self.title = 'Specify a value of dynamic force in the left spring:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class RightSpringForceMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj:
        Eq(Symbol('F_sr'),
           (-HarmonicOscillator(obj.linearized()).stiffness_matrix()[1] *
            (HarmonicOscillator(obj.linearized()).general_solution().n(3)[0] -
             HarmonicOscillator(obj.linearized()).general_solution().n(3)[1])),
           evaluate=False),
            **kwargs):
        self.title = 'Podaj wartość siły dynamicznej w prawej sprężynie:'
        self.title = 'Specify a value of dynamic force in the right spring:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class PeriodMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [
                     2 * pi / ((freq_val)) for freq_val in HarmonicOscillator(
                         obj.linearized()).natural_frequencies()
                     if freq_val != 0
                 ],
                 **kwargs):
        self.title = 'Podaj wartość okresu:'
        self.title = 'Specify the value of natural periods of vibrations:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class FrequencyResponseFunctionMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('FRF'),
                                            HarmonicOscillator(obj.linearized(
                                            )).frequency_response_function(),
                                            evaluate=False),
            **kwargs):
        self.title = 'Podaj wartość FRF:'
        self.title = 'Calculate the FRF:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


FRFMCA = FrequencyResponseFunctionMCA


class FRFforOmega0(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj:
        Eq(Symbol('FRF(\omega_0)'),
           HarmonicOscillator(obj.linearized()).frequency_response_function()
           .subs(Omega,
                 sqrt(obj.stiffness_matrix()[0] / obj.inertia_matrix()[0])),
           evaluate=False),
            **kwargs):
        self.title = 'Wyznaczać wartość FRF dla częstości drgań swobodnych nietłumionych:'
        self.title = 'Determine the FRF value for undamped natural frequency:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class DynamicAmplificationFactorMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            #                  answer_generator=lambda obj: Eq(Symbol(
            #                      'DAF'),HarmonicOscillator(obj.linearized()).stiffness_matrix()[0]/(HarmonicOscillator(obj.linearized()).stiffness_matrix()[0]-Omega**2*HarmonicOscillator(obj.linearized()).inertia_matrix()[0]),
            #                                                  evaluate=False),
            answer_generator=lambda obj:
        Eq(Symbol('DAF'), obj.dynamic_amplification_factor(), evaluate=False),
            **kwargs):
        self.title = 'Podaj wartość DAF:'
        self.title = 'Calculate the DAF:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


DAFMCA = DynamicAmplificationFactorMCA


class PotentialEnergyMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('V'),
                     -obj.lagrangian().subs(
                         {coord: 0
                          for coord in Matrix(obj.q).diff(t)}),
                     evaluate=False),
                 **kwargs):
        self.title = 'Podaj wartość energi potencjalnej:'
        self.title = 'Estimate value of the potential energy:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class KineticEnergyMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('T'),
                     obj.lagrangian() - obj.lagrangian().subs(
                         {coord: 0
                          for coord in Matrix(obj.q).diff(t)}),
                     evaluate=False),
                 **kwargs):
        self.title = 'Podaj wartość energi kinetycznej:'
        self.title = 'Estimate value of the kinetic energy:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class DampingFactorMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('c'), obj.damping_matrix()[0], evaluate=False),
                 **kwargs):
        self.title = 'Podaj wartość współczynnika tłumienia układa:'
        self.title = 'What is the value of a damping factor?'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFDampedOmegaMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('omega_h'),
                     list(
                         HarmonicOscillator(obj.linearized()).
                         damped_natural_frequencies())[0]),
                 **kwargs):

        self.title = 'Określ częstość tłumionych drgań swobodnych występujących w układzie:'
        self.title = 'Determine the system natural frequency of damped vibration:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class DampedOmegaMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [
                     Eq(Symbol('omega_h'), om_eq) for om_eq in list(
                         HarmonicOscillator(obj.linearized()).
                         damped_natural_frequencies()) if om_eq
                 ],
                 **kwargs):

        self.title = 'Określ częstości tłumionych drgań swobodnych występujących w układzie:'
        self.title = 'Determine the system natural frequencies of damped vibration:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SmallParameterMCA(MechanicalSystemAnswer):
    question = 'Small parameter for the cosidered system is expressed as:'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('varepsilon'), (
                obj.small_parameter() / obj.inertia_matrix()[0]).simplify()),
            **kwargs):

        #         self.title = 'Mały parametr układu wyraża się wzorem:'
        #self.title = '''Find the formula that represents the system's Lagrangian:'''
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class ResonanceCurveMCA(MechanicalSystemAnswer):
    question = 'A relation between the amplitude and natural frequency of free vibrations for the considered systen is given as:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('omega')**2,
                     (HarmonicOscillator(obj.linearized()).natural_frequencies(
                     )[0]**2 + S('3') / 4 * obj.small_parameter() / obj.
                      inertia_matrix()[0] * Symbol('a')**2).expand()),
                 **kwargs):

        #         self.title = 'Zależność pomiędzy amplitudą a częstością drgań własnych dla drgań swobodnych rozważanego układu wyraża wzór:'
        #self.title = '''Find the formula that represents the system's Lagrangian:'''
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class FundamentalMatrixMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('A'),
                                            (HarmonicOscillator(obj.linearized(
                                            )).fundamental_matrix()),
                                            evaluate=False),
            **kwargs):
        self.title = 'A fundamental matrix for the system under consideration is given by:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class FundamentalMatrixDeterminantMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('\Delta'),
                                            (HarmonicOscillator(obj.linearized(
                                            )).fundamental_matrix().det()),
                                            evaluate=False),
            **kwargs):
        self.title = 'Determine the characteristic polynomial for the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class GeneralizedMomentumsMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('p'),
                                                 obj.generalized_momentum()),
                 **kwargs):

        self.title = 'What are generalised momentums of the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class EquilibriumEquationsMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(obj.equilibrium_equation().doit(),
                                            Matrix([0] * len(obj.q))),
            **kwargs):

        self.title = 'What are equilibrium equations for the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SDoFEquilibriumEquationsMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     obj.equilibrium_equation().doit()[0], 0),
                 **kwargs):

        self.title = 'What are equilibrium equations for the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class InputAmplitudeMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('F_0'),
                                            Symbol('z_0') /
                                            (obj.stiffness_matrix()[0] /
                                             (obj.stiffness_matrix()[0] - Omega
                                              **2 * obj.inertia_matrix()[0])),
                                            evaluate=False),
            **kwargs):
        self.title = 'Jaka musi być amplituda wejścia, aby układ drgał z amplitudą \(z_{0} \)?'
        self.title = 'Determine the value of an input amplitude ensuring the system vibrations of the \(\z_{0}\) amplitude:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SpringForceMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('F'), (
                (obj.stiffness_matrix()[0]) * obj.steady_solution().n(3)[0]),
                                            evaluate=False),
            **kwargs):
        self.title = 'Podaj wartość siły dynamicznej w jednej ze sprężyn:'
        self.title = 'Determine the value of a dynamic force occurring in the spring:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class EigenvalsMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [
                     Eq(Symbol('lambda'), eig_val) for eig_val in
                     (HarmonicOscillator(obj.linearized())).eigenvalues()
                     if eig_val != 0
                 ],
                 **kwargs):

        self.title = 'Określ częstość drgań swobodnych występujących w układzie:'
        self.title = 'Determine the eigenvalues for the considered system:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SmallParameterNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Specify a small parameter value for the following data \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj:
                 ((obj.small_parameter() / obj.inertia_matrix()[0]).simplify()
                  ).subs({
                      Symbol('k_0', positive=True): 300,
                      Symbol('m_0', positive=True): 500,
                      Symbol('d', positive=True): 0.5,
                      Symbol('l_0', positive=True): 0.3
                  }).n(3),
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class FrequencyForAmplitudeNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find the nonlinear vibration frequency for the given amplitude \(a=0.1 m\), if \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: (solve(
                (HarmonicOscillator(obj.linearized()).natural_frequencies()[0]
                 **2 + S('3') / 4 * obj.small_parameter() / obj.inertia_matrix(
                 )[0] * Symbol('a')**2).subs({
                     Symbol('k_0', positive=True): 300,
                     Symbol('m_0', positive=True): 500,
                     Symbol('d', positive=True): 0.5,
                     Symbol('l_0', positive=True): 0.3,
                     Symbol('a'): 0.1
                 }) - Symbol('omega')**2, Symbol('omega')))[1].n(3),
            **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class NonLinLowerAmplitudeValueNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find a lower value of the amplitude of nonlinear vibrations for frequency \(\omega=3.5 \\frac{rad}{s}\), if \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: (solve(
                (HarmonicOscillator(obj.linearized()).natural_frequencies()[0]
                 **2 + S('3') / 4 * obj.small_parameter() / obj.inertia_matrix(
                 )[0] * Symbol('a')**2) - Symbol('omega', positive=True)**2,
                Symbol('a')))[0].subs({
                    Symbol('k_0', positive=True): 300,
                    Symbol('m_0', positive=True): 500,
                    Symbol('d', positive=True): 0.5,
                    Symbol('l_0', positive=True): 0.3,
                    Symbol('omega', positive=True): 3.5
                }).n(3),
            **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class NonLinHigherAmplitudeValueNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find a higher value of the amplitude of nonlinear vibrations for frequency \(\omega=3.5 \\frac{rad}{s}\), if \(m=500\) \(kg\), \(k=300 \\frac{N}{m}\), \(d=0.5\) \(m\), \(l_0=0.3\) \(m\):'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: (solve(
                (HarmonicOscillator(obj.linearized()).natural_frequencies()[0]
                 **2 + S('3') / 4 * obj.small_parameter() / obj.inertia_matrix(
                 )[0] * Symbol('a')**2) - Symbol('omega', positive=True)**2,
                Symbol('a')))[1].subs({
                    Symbol('k_0', positive=True): 300,
                    Symbol('m_0', positive=True): 500,
                    Symbol('d', positive=True): 0.5,
                    Symbol('l_0', positive=True): 0.3,
                    Symbol('omega', positive=True): 3.5
                }).n(3),
            **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SystemResponseAmplitudeC2MCA(MechanicalSystemAnswer):
    question = 'Determine the \(C2\) amplitude of a steady solution:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [Eq(Symbol('C2'),(list(
                     obj.steady_solution_amp(obj.external_forces().subs(
                         t, 0), Matrix([0, 0]))[0])[1]).simplify())],
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class LowerNaturalFreqValueNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find the value of the lower natural frequency of the system for the given parameters \(m=150 kg\), \(k=900 \\frac{N}{m}\:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: min([
                     eig for eig in HarmonicOscillator(obj.linearized()).
                     natural_frequencies().subs({
                         Symbol('k', positive=True): 900,
                         Symbol('m', positive=True): 500
                     }).applyfunc(lambda comp: round(comp.n(), 2)) if eig != 0
                 ]),
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class HigherNaturalFreqValueNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find the value of the higher natural frequency of the system for the given parameters \(m=150 kg\), \(k=900 \\frac{N}{m}\:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: max([
                     eig for eig in HarmonicOscillator(obj.linearized()).
                     natural_frequencies().subs({
                         Symbol('k', positive=True): 900,
                         Symbol('m', positive=True): 500
                     }).applyfunc(lambda comp: round(comp.n(), 2)) if eig != 0
                 ]),
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SteadyAmplitudesRatioNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'How many times the steady amplitude of the first trolley is greater than the amplitude of steady vibrations of the second trolley if parameters \(m=150 kg\), \(k=900 \\frac{N}{m}\), \(\\Omega=15 \\frac{rad}{s}\):'

    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: (abs((list(
                obj.steady_solution_amp(obj.external_forces(
                ).subs(t, 0), Matrix([0, 0]))[0])[0]).simplify() / (list(
                    obj.steady_solution_amp(obj.external_forces().subs(t, 0),
                                            Matrix([0, 0]))[0]
                )[1]).simplify())).subs({
                    Symbol('k', positive=True): 900,
                    Symbol('m', positive=True): 500,
                    Symbol('Omega', positive=True): 15
                }).n(3),
            **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class SecondModeFirstComponentValueNMSA(TitledNumericalAnswerForMechanicalSystem):
    question = 'Find the value of the first component of the second mode of vibration for the given parameters \(m=150 kg\), \(k=900 \\frac{N}{m}\):'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (HarmonicOscillator(
                     obj.linearized()).modes()[:,1]).subs({Symbol('k', positive=True): 900,
                    Symbol('m', positive=True): 500,}).applyfunc(lambda comp: round(comp.n(),3))[0],
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)

class Test(MechanicalSystemAnswer):
    question = 'Determine the \(C2\) amplitude of a steady solution:'

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (abs((list(
                obj.steady_solution_amp(obj.external_forces(
                ).subs(t, 0), Matrix([0, 0]))[0])[0]).simplify() / (list(
                    obj.steady_solution_amp(obj.external_forces().subs(t, 0),
                                            Matrix([0, 0]))[0]
                )[1]).simplify())).subs({
                    Symbol('k', positive=True): 900,
                    Symbol('m', positive=True): 500,
                    Symbol('Omega', positive=True): 15
                }).n(3),
                 **kwargs):

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)
