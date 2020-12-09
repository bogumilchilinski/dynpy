import sympy as sym
import sympy.physics.mechanics as mech

class EmbeddedAnswer:
    
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
    def __init__(self,answer,score=1):
        super().__init__(value=answer,error=0.1,relative_error=True,precision=4,score=1,question_str='SHORTANSWER')

        
    def to_string(self):
        
        answer_value=self.value
        
        
        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':='+ str(answer_value) +'}' 
    
        return answer_string

    
class EmbeddedMultichoiceAnswer(EmbeddedAnswer):    
    def __init__(self,correct_answers,wrong_answers=None,score=1):
        super().__init__(value=correct_answers,error=0.1,relative_error=True,precision=4,score=1,question_str='MULTICHOICE_VS')
        self.correct_answers=correct_answers
        self.wrong_answers=wrong_answers
        
    def to_string(self):
        
        answer_value=self.value
        
        
        answer_string = '{'+ str(self.score)+ ':' +self.question_str+':='+ '~='.join(self.correct_answers) +'~' + '~'.join(self.wrong_answers) +  '}' 
    
        return answer_string
    
    
    
class EmbeddedMultichoiceMathAnswer(EmbeddedMultichoiceAnswer):    

    def __init__(self,correct_answers,wrong_answers=None,score=1,backend=mech.vlatex):
        super().__init__(correct_answers=correct_answers,wrong_answers=wrong_answers,score=score)
        
        self.backend=backend
    
    def to_string(self):
        
        answer_value=self.value
        
        correct_answers_string='~='.join(['\('+self.backend(ans)+'\)' for  ans in self.correct_answers])
        wrong_answers_string='~'.join(['\('+self.backend(ans)+'\)' for  ans in self.wrong_answers])
        
        
        answer_string = '{'+ str(self.score)+ ':' +(self.question_str+':='+ correct_answers_string +'~' +  wrong_answers_string).replace('}','\}') +  '}' 
    
        return answer_string    
    
        
        
class Question:
    
    def __init__(self,question,id=1,title='Cloze_question',question_type='cloze'):
        
        self.question_str=question
        self.question_type=question_type
        self.title=title
        self.number=id
        
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
            """.format(q_no=self.number,q_name=question_name)

        question_string=self.question_str

        xml_string_post="""]]>
        </text>
        </questiontext>
            <generalfeedback format="html">
              <text></text>
            </generalfeedback>
            <penalty>0.3333333</penalty>
            <hidden>0</hidden>
            <idnumber></idnumber>
          </question>
        """
        return xml_string_pre+question_string+xml_string_post
    
    def __str__(self):
        return self.to_string()
    
    def __repr__(self):
        return 'MoodleQuestion('+self.question_str+')'

class Category:
    def __init__(self,name,questions,id=0):
        
    def to_string(self):
        question_name=self.title+"_"+str(self.number)
        
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

        """.format(cat=name)

        end_string="</quiz>"

        return (begin_string+ cat_string +''.join(question_set)+end_string)