import sympy as sym
import sympy.physics.mechanics as mech
import base64

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