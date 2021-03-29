import sympy as sym
from sympy import *
import sympy.physics.mechanics as mech
from .dynamics import *

import base64
t=Symbol('t')

class EmbeddedAnswer:
    """ Class EmeddedAnswer allows to create a space for numerical answer.
    
    Class EmbeddedAnswer contains answer's value, error, possible inaccuracy, score and precision of decimal place.
    
    Parameters:
    ============
    value: float
        The user's answer.
    
    error: float
        value of allowed inaccuracy.
    
    relative_error: bool
        Gives oportunity to make an uncertainty in negative and positive way.
    
    precision: int
        Value of precision of decimal place.
    
    score: int
        Value of question score.
    """
    def __init__(self,value,error=0.1,relative_error=True,precision=4,score=1,question_str='NUMERICAL'):
        
        self.value=value
        self.error=error
        self.relative_error=relative_error
        self.precision=precision
        self.score=score
        self.question_str=question_str
        
    def to_string(self):
        
        answer_value=(round(self.value,self.precision))
        abs_error=round(answer_value*self.error,self.precision)
        
        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':='+ str(answer_value) +':'+str(abs_error) +'}' 
    
        return answer_string
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'Moodle'+ self.__class__.__name__ +self.to_string()
    def to_string_alt(self):

        answer_value=(round(self.value,self.precision))
        abs_error=round(answer_value*self.error,self.precision)

        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':%100%'+ str(answer_value) +':'+str(abs_error) +'}' 

        return answer_string

class EmbeddedShortAnswer(EmbeddedAnswer):    
    """ Class EmeddedShortAnswer allows to create a space for short answer. It inherits from class EmbeddedAnswer and creates blank place for string type variable.
    
    Class EmbeddedShortAnswer contains answer's string and score.
    
    Parameters:
    ============
    value: str
        The user's answer.
    
    score: int
        Value of question score.
    """
    def __init__(self,answer,score=1):
        super().__init__(value=answer,error=0.1,relative_error=True,precision=4,score=1,question_str='SHORTANSWER')

        
    def to_string(self):
        
        answer_value=self.value
        
        
        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':='+ str(answer_value) +'}' 
    
        return answer_string

    
class EmbeddedMultichoiceAnswer(EmbeddedAnswer):    
    """ Class EmeddedMultichoiceAnswer allows to create an answear with more than one correct value. It inherits from class EmbeddedAnswer and creates space for multichoice answer.
    
    Class EmbeddedMultichoiceAnswer contains correct answer's value,wrong answer's value and error.
    
    Parameters:
    ============
    correct_answers: str
        contains correct answers.
    
    wrong_answers: str
        contains wrong answers.
    
    score: int
        Value of question score.
    """
    def __init__(self,correct_answers,wrong_answers=None,score=1,**kwargs):
        super().__init__(value=correct_answers,error=0.1,relative_error=True,precision=4,score=1,question_str='MULTICHOICE_VS')
        self.correct_answers=correct_answers
        self.wrong_answers=wrong_answers
        
    def to_string(self):
        
        answer_value=self.value
        
        
        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':='+ '~='.join(self.correct_answers) +'~' + '~'.join(self.wrong_answers) +  '}' 
    
        return answer_string
    
class EmbeddedGraphics:
    """
    Class EmbeddedGraphics allows input graphics inside the question.
    
    Parameters:
    ============
    path_to_file: str
        Path to graphics
    
    """
    def __init__(self,path_to_file):
        self.filepath=path_to_file
        self.filename=path_to_file.split('/')[-1]
        
    def to_string(self):
        
        figure_str='''
                    <img src="@@PLUGINFILE@@/{figure_name}" />
                    '''.format(figure_name=self.filename)
    
        return figure_str
    
    def to_base64(self):
        with open(self.filepath, "rb") as image_file:
            encoded_pic= base64.b64encode(image_file.read()).decode('utf-8')
            
        return encoded_pic

    
    def to_base64_entry(self):
        
        figure_str='''
        <file name="{filename}" path="/" encoding="base64">
        {graphics}
        </file>
        '''.format(filename=self.filename,graphics=self.to_base64())

        return figure_str
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'Moodle'+ self.__class__.__name__ +self.to_string()
    
class EmbeddedMultichoiceMathAnswer(EmbeddedMultichoiceAnswer):    
    """
    Class EmeddedMultichoiceMathAnswer allows to create an answear with more than one correct value. It inherits from class EmbeddedMultichoiceAnswer and creates space for multichoice answer.
    
    Parameters:
    ============
    
    correct_answers: str
        Correct answer
   
    wrong_answer: str
       Wrong answer
   
    score: int
       Score for correct answer
   
    backend: str
        The way of displaying equations
    """

    def __init__(self,correct_answers,wrong_answers=None,score=1,backend=mech.vlatex,**kwargs):
        super().__init__(correct_answers=correct_answers,wrong_answers=wrong_answers,score=score,**kwargs)
        
        self.backend=backend
    
    def to_string(self):
        
        answer_value=self.value
        
        correct_answers_string='~='.join(['\('+self.backend(ans)+'\)' for  ans in self.correct_answers])
        wrong_answers_string='~'.join(['\('+self.backend(ans)+'\)' for  ans in self.wrong_answers])
        
        
        answer_string = '{'+ str(self.score)+ ':' +(self.question_str+':='+ correct_answers_string +'~' +  wrong_answers_string).replace('}','\}') +  '}' 
    
        return answer_string    
    
    
#class 
        
class Question:
    """
    Instance of this class creates complete question with answers within it. Additionaly it wraps question with needed xml commends.
    
    Parameters:
    ============
    
    question: str
        Question content
   
    id=1: int
       Number of question
   
    title='Cloze_question': str
       Title of question
   
    question_type='cloze': str
       Defines type of question
    """
    def __init__(self,entries_list,id=1,figure='',title='Cloze_question',question_type='cloze'):
        
        if not isinstance(entries_list,list):
            entries_list=[entries_list]
        
        self.entries_list=entries_list
        self.question_type=question_type
        self.title=title
        self.number=id
        self.figure=figure
        
    def to_string(self):
        question_name=self.title+"_"+str(self.number)
        
        xml_string_pre="""<!-- question: {q_no}  -->
          
          <question type="cloze">
            <name>
              <text>{q_name}</text>
            </name>
            <questiontext format="html">
            <text>
            <![CDATA[
            """.format(q_no=self.number,q_name=question_name,graphics=self.figure)

        question_string='\n'.join([str(entry)  for entry in self.entries_list])

        xml_string_post="""]]>
        
        </text>
        {graphics}
        </questiontext>
        
            <generalfeedback format="html">
              <text></text>
            </generalfeedback>
            <penalty>0.3333333</penalty>
            <hidden>0</hidden>
            <idnumber></idnumber>
          </question>
        """.format(graphics=self.figure)
        
        return xml_string_pre+question_string+xml_string_post
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'MoodleQuestion('+self.entries+')'

class Paragraph:
    """
    Object of this class wraps up question string with HTML paragraphs.
    
    Parameters
    ==========
    
    element : str
        Text to wrap with paragraphs
    
    """
    def __init__(self,element):
        self.elem=element
        
    def to_string(self):
        return "<p> {obj}  </p> ".format(obj=str(self.elem))
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'MoodleParagraph('+str(self.elem)+')'
    
class Category:
    """
    Object of this class creates and represents document file of type xml which contains question list of specific category.
    
    Parameters:
    ============
    
    name: str
        Name of category of questions
   
    questions: str
       Question list
   
    id: int
       Number id of first question
    """
    def __init__(self,name,questions,id=0):
        self.id=id
        self.name=name
        self.questions=questions
        
        
    def to_string(self):
        question_name=self.name+"_"+str(self.id)
        
        begin_string="""<?xml version="1.0" encoding="UTF-8"?>
                        <quiz>"""


        cat_string="""
        <!-- question: 0  -->
          <question type="category">
            <category>
              <text>$course$/top/Domyślna dla: 386477#129825/{cat}</text>
            </category>
            <info format="html">
              <text></text>
            </info>
            <idnumber></idnumber>
          </question>

        """.format(cat=self.name)

        end_string="</quiz>"

        
        
        return (begin_string+ cat_string +'\n'.join(self.questions)+end_string)
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'MoodleCategory('+self.question_str+')'
    
    
import base64

class MechanicalSystemAnswer(EmbeddedMultichoiceMathAnswer):
    
    question=None
    
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=None,
                 title=None,
                 **kwargs):
        
        if title:
            self.title = title
        else:
            self.title=type(self).question
            
        if not answer_generator:
            answer_generator= self.answer_entry
        
        self._correct_answers = [answer_generator(system) for system in sym.flatten([correct_system])]
        self._other_answers = [answer_generator(system) for system in other_systems]
        
        super().__init__(self._correct_answers, self._other_answers, **kwargs)
   
    def answer_entry(self,system):
        return system
        
        
    
    def to_string(self):
        return self.title + '\n' + super().to_string()

    def preview(self, backend=None):
        print(self.title)
        print('=' * 100)
        display(*self._correct_answers)
        print('=' * 100)
        print('x' * 100)
        display(*self._other_answers)
        print('x' * 100)

        
        
########################### MDOF ########################


class GoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (Eq(obj._eoms.doit().expand().applyfunc(simplify).expand(),Matrix([0]*len(obj.q)),evaluate=False)),
                 **kwargs):

        self.title = 'Podaj równiania ruchu układu:'
        self.title = 'Determine equation of motion of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)

        

class SDoFGoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (Eq(obj._eoms.doit().expand().applyfunc(simplify).expand()[0],0,evaluate=False)),
                 **kwargs):

        self.title = 'Podaj równiania ruchu układu:'
        self.title = 'Determine equation of motion of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)       

class SDoFLinearizedGoverningEquationMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: (Eq(obj.linearized()._eoms.doit().expand().applyfunc(simplify).expand()[0],0,evaluate=False)),
                 **kwargs):

        self.title = 'Podaj równiania ruchu układu:'
        self.title = 'Determine equation of motion of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)      
        
        
        
class CriticalPointsMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Matrix(obj.critical_points()[0:2]),#
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
                 answer_generator=lambda obj:(Eq( HarmonicOscillator(obj.linearized())._eoms.doit(),Matrix([0]*len(obj.q)),evaluate=False)),
                 **kwargs):

        self.title = 'Liniowe równania ruchu dla układu przedstawionego na rysunku można wyrazić następującym układem równań:'
        self.title = 'Linear equations of motion for the system presented on a picture might be expressed as following system of equations:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)
        

class LagrangianMCA(MechanicalSystemAnswer):
    question= 'Choose corect dependance which determines Lagrangian of the considered system:'
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('L'), obj.lagrangian()),
            **kwargs):

        self.title = 'Wskaż zależność określającą Lagrangian rozpatrywanego układu:'
        self.title = 'Choose corect dependance which determines Lagrangian of the considered system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=None,
                         **kwargs)


class LagrangianMCAPL(LagrangianMCA):
    question= 'Wskaż zależność określającą Lagrangian rozpatrywanego układu:'        
        
class ExternalForcesMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('F'), obj.external_forces(),evaluate=False),
            **kwargs):

        self.title = 'Określ wektor sił w układzie:'
        self.title = 'Specify a forcing vector of the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)

        
class LinearizedLagrangianMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('L'),
                     HarmonicOscillator(obj.linearized()).lagrangian()), # op_point to remove | True or False - doesn't matter
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
                     (eig_val) for eig_val in HarmonicOscillator(
                         obj.linearized()).natural_frequencies() if eig_val != 0
                 ],
                 **kwargs):

        self.title = 'Określ częstości drgań swobodnych występujących w układzie:'
        self.title = 'Determine the frequency of the free vibrations occuring in the system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)

class SDoFOmegaMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: [
                     (eig_val) for eig_val in HarmonicOscillator(
                         obj.linearized()).natural_frequencies() if eig_val != 0
                 ][0],
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


class SolutionMCA(MechanicalSystemAnswer):
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
        self.title = 'Determine a general solution of ODE for provided system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class SteadySolutionMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(
                     Symbol('X_s'),
                     HarmonicOscillator(obj.linearized()).steady_solution().n(3),
                     evaluate=False),
                 **kwargs):

        self.title = 'Wyznacz rozwiązanie szczególne dla rozważanego układu:'
        self.title = 'Determine a particular solution of ODE for provided system:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class InertiaMatrixMCA(MechanicalSystemAnswer):
    def __init__(
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('M'), (HarmonicOscillator(
                obj.linearized()).inertia_matrix()),
                                            evaluate=False),
            **kwargs):
        self.title = 'Określ macierz bezwładności układu:'
        self.title = 'Determine the inertia matrix:'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class StiffnessMatrixMCA(MechanicalSystemAnswer):
    def __init__(
        
            self,
            correct_system,
            other_systems,
            answer_generator=lambda obj: Eq(Symbol('K'), (HarmonicOscillator(
                obj.linearized()).stiffness_matrix()),
                                            evaluate=False),
            **kwargs):
        self.title = 'Określ macierz sztywności układu:'
        self.title = 'Determine the stiffness matrix of our system'

        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class LeftSpringForceMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda
                 obj: Eq(Symbol('F_sl'),
                         ((HarmonicOscillator(obj.linearized()).stiffness_matrix()[0] + HarmonicOscillator(obj.linearized()).stiffness_matrix()[
                             1])) * HarmonicOscillator(obj.linearized()).general_solution().n(3)[0],
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
            answer_generator=lambda obj: Eq(Symbol('F_sr'),
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
                     2*pi / ((freq_val))
                     for freq_val in HarmonicOscillator(obj.linearized()).natural_frequencies() if freq_val != 0
                 ][0],
                 **kwargs):
        self.title = 'Podaj wartość okresu:'
        self.title = 'Specify the value of a period:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)


class FrequencyResponseFunctionMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol(
                     'FRF'), HarmonicOscillator(obj.linearized()).frequency_response_function(),
                                            evaluate=False),
                 **kwargs):
        self.title = 'Podaj wartość FRF:'
        self.title = 'Calculate the FRF:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)

FRFMCA=FrequencyResponseFunctionMCA


class FRFforOmega0(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol(
                     'FRF(\omega_0)'), HarmonicOscillator(obj.linearized()).frequency_response_function().subs(Omega,sqrt(obj.stiffness_matrix()[0]/obj.inertia_matrix()[0])),
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
    def __init__(self,
                 correct_system,
                 other_systems,
#                  answer_generator=lambda obj: Eq(Symbol(
#                      'DAF'),HarmonicOscillator(obj.linearized()).stiffness_matrix()[0]/(HarmonicOscillator(obj.linearized()).stiffness_matrix()[0]-Omega**2*HarmonicOscillator(obj.linearized()).inertia_matrix()[0]),
#                                                  evaluate=False),
                 answer_generator=lambda obj: Eq(Symbol(
                      'DAF'),obj.dynamic_amplification_factor(),evaluate=False),
                 **kwargs):
        self.title = 'Podaj wartość DAF:'
        self.title = 'Calculate the DAF:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs)        
DAFMCA=DynamicAmplificationFactorMCA        

class PotentialEnergyMCA(MechanicalSystemAnswer):
    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=lambda obj: Eq(Symbol('V'), -obj.lagrangian().subs({ coord:0 for coord in Matrix(obj.q).diff(t)}),
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
                 answer_generator=lambda obj: Eq(Symbol('T'),obj.lagrangian() -obj.lagrangian().subs({ coord:0 for coord in Matrix(obj.q).diff(t)}),
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
                 answer_generator=lambda obj: Eq(Symbol('c'), obj.damping_matrix()[0]   ,
                                                 evaluate=False),
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
                 answer_generator=lambda obj: Eq(Symbol('omega_h'),list(HarmonicOscillator(obj.linearized()).damped_natural_frequencies())[0]),
                 **kwargs):

        self.title = 'Określ częstość tłumionych drgań swobodnych występujących w układzie:'
        self.title = 'Determine the system natural frequency of damped vibration:'
        super().__init__(correct_system,
                         other_systems,
                         answer_generator=answer_generator,
                         title=self.title,
                         **kwargs) 