import sympy as sym
from sympy import *
import sympy.physics.mechanics as mech
from sympy.physics.vector.printing import vpprint, vlatex
import dynpy
from dynpy import HarmonicOscillator
import dynpy.models.systems as sys
import random as rand
from random import *
import IPython as IP
import base64
t = Symbol('t')


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

    def __init__(self, value, error=0.1, relative_error=True, precision=4, score=1, question_str='NUMERICAL'):

        self.value = value
        self.error = error
        self.relative_error = relative_error
        self.precision = precision
        self.score = score
        self.question_str = question_str

    def to_string(self):

        answer_value = (round(self.value, self.precision))
        abs_error = round(answer_value*self.error, self.precision)

        answer_string = '{' + str(self.score) + ':' + self.question_str + \
            ':=' + str(answer_value) + ':'+str(abs_error) + '}'

        return answer_string

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return 'Moodle' + self.__class__.__name__ + self.to_string()

    def to_string_alt(self):

        answer_value = (round(self.value, self.precision))
        abs_error = round(answer_value*self.error, self.precision)

        answer_string = '{' + str(self.score) + ':' + self.question_str + \
            ':%100%' + str(answer_value) + ':'+str(abs_error) + '}'

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

    def __init__(self, answer, score=1):
        super().__init__(value=answer, error=0.1, relative_error=True,
                         precision=4, score=1, question_str='SHORTANSWER')

    def to_string(self):

        answer_value = self.value

        answer_string = '{' + str(self.score) + ':' + \
            self.question_str+':=' + str(answer_value) + '}'

        return answer_string

class EmbeddedNumericalAnswer(EmbeddedAnswer):
    """ Class EmeddedNumericalAnswer allows to create a space for numerical answer. It inherits from class EmbeddedAnswer and creates blank place for string type variable.

    Class EmbeddedShortAnswer contains answer's string and score.

    Parameters:
    ============
    value: str
        The user's answer.

    score: int
        Value of question score.
    """

    def __init__(self, answer, score=1):
        super().__init__(value=answer, error=0.1, relative_error=True,
                         precision=4, score=1, question_str='NUMERICAL')


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

    def __init__(self, correct_answers, wrong_answers=None, score=1, **kwargs):
        super().__init__(value=correct_answers, error=0.1, relative_error=True,
                         precision=4, score=1, question_str='MULTICHOICE_VS')
        self.correct_answers = correct_answers
        self.wrong_answers = wrong_answers

    def to_string(self):

        answer_value = self.value

        answer_string = '{' + str(self.score) + ':' + self.question_str+':=' + '~='.join(
            self.correct_answers) + '~' + '~'.join(self.wrong_answers) + '}'

        return answer_string


class EmbeddedGraphics:
    """
    Class EmbeddedGraphics allows input graphics inside the question.

    Parameters:
    ============
    path_to_file: str
        Path to graphics

    """

    def __init__(self, path_to_file):
        self.filepath = path_to_file
        self.filename = path_to_file.split('/')[-1]

    def to_string(self):

        figure_str = '''
                    <img src="@@PLUGINFILE@@/{figure_name}" />
                    '''.format(figure_name=self.filename)

        return figure_str

    def to_base64(self):
        with open(self.filepath, "rb") as image_file:
            encoded_pic = base64.b64encode(image_file.read()).decode('utf-8')

        return encoded_pic

    def to_base64_entry(self):

        figure_str = '''
        <file name="{filename}" path="/" encoding="base64">
        {graphics}
        </file>
        '''.format(filename=self.filename, graphics=self.to_base64())

        return figure_str

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return 'Moodle' + self.__class__.__name__ + self.to_string()


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

    def __init__(self, correct_answers, wrong_answers=None, score=1, backend=mech.vlatex, **kwargs):
        super().__init__(correct_answers=correct_answers,
                         wrong_answers=wrong_answers, score=score, **kwargs)

        self.backend = backend

    def to_string(self):

        answer_value = self.value

        correct_answers_string = '~='.join(
            ['\('+self.backend(ans)+'\)' for ans in self.correct_answers])
        wrong_answers_string = '~'.join(
            ['\('+self.backend(ans)+'\)' for ans in self.wrong_answers])

        answer_string = '{' + str(self.score) + ':' + (self.question_str+':=' +
                                                       correct_answers_string + '~' + wrong_answers_string).replace('}', '\}') + '}'

        return answer_string


# class

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

    def __init__(self, entries_list, id=1, figure='', title='Cloze_question', question_type='cloze'):

        if not isinstance(entries_list, list):
            entries_list = [entries_list]

        self.entries_list = entries_list
        self.question_type = question_type
        self.title = title
        self.number = id
        self.figure = figure

    def to_string(self):
        question_name = self.title+"_"+str(self.number)

        xml_string_pre = """<!-- question: {q_no}  -->
          
          <question type="cloze">
            <name>
              <text>{q_name}</text>
            </name>
            <questiontext format="html">
            <text>
            <![CDATA[
            """.format(q_no=self.number, q_name=question_name, graphics=self.figure)

        question_string = '\n'.join([str(entry)
                                    for entry in self.entries_list])

        xml_string_post = """]]>
        
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

    def __init__(self, element):
        self.elem = element

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

    def __init__(self, name, questions, id=0):
        self.id = id
        self.name = name
        self.questions = questions

    def to_string(self):
        question_name = self.name+"_"+str(self.id)

        begin_string = """<?xml version="1.0" encoding="UTF-8"?>
                        <quiz>"""

        cat_string = """
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

        end_string = "</quiz>"

        return (begin_string + cat_string + '\n'.join(self.questions)+end_string)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return 'MoodleCategory('+self.question_str+')'


class MechanicalSystemAnswer(EmbeddedMultichoiceMathAnswer):

    question = None

    def __init__(self,
                 correct_system,
                 other_systems,
                 answer_generator=None,
                 title=None,
                 **kwargs):

        

        if title:
            self.title = title
        else:
            self.title = type(self).question

        if not answer_generator:
            answer_generator = self.answer_entry

        self.answer_generator=answer_generator
        self._correct_answers = [answer_generator(
            system) for system in sym.flatten([correct_system])]
        self._other_answers = [answer_generator(
            system) for system in other_systems]

        super().__init__(self._correct_answers, self._other_answers, **kwargs)

    def answer_entry(self, system):
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

class TitledNumericalAnswerForMechanicalSystem(EmbeddedNumericalAnswer):

    question = None

    def __init__(self,
                 correct_system,
                 answer_generator=None,
                 title=None,
                 **kwargs):
#         print('hello')

        if title:
            self.title = title
        else:
            self.title = type(self).question

        if not answer_generator:
            answer_generator = self.answer_entry

        self._correct_answers = answer_generator(correct_system)

        super().__init__(self._correct_answers, **kwargs)

    def to_string(self):
        return self.title + '\n' + super().to_string()

    def preview(self, backend=None):
        print(self.title)
        print('=' * 100)
        display(self._correct_answers)
        print('=' * 100)

        

class QuizOn(sys.ComposedSystem):
    
    def __init__(self,*args,**kwargs):
       
        self.default_data_dict=args[0].get_default_data()
        print(self.default_data_dict)
        self.system=HarmonicOscillator(args[0])
        
        
        self._scheme_path=args[0]._scheme()
        self._real_example_path=args[0]._real_example()
        super().__init__(*args,**kwargs)
        

    def generate_dict(self,cases_no=5,param_range=None):
        sys_par=self.system_parameters()
        self.param_range=param_range
        sym_list=[]
        
        if self.param_range==None:
            self.param_range=self.default_data_dict


        if self.param_range==None:

            for caseno in range(cases_no):
                print('jestes tu :(')
                sym_dict={}
                temp_dict={}
                for num,sym in enumerate(sys_par):
                    split_name_default=str(sys_par[num]).split('_',1)

                    sym_dict[sym]=Symbol(str(rand.randrange(1,20,1))+split_name_default[0]+'_0',positive=True)

#                 for key,val in param_range.items(): 
#                     if isinstance(val,list)==True:
#                         for elem in val:
#                             if isinstance(elem,Expr):
#                                 sym_dict[key]=rand.choice(val)
#                             else:
#                                 split_name_dict=str(key).split('_',1)
#                                 sym_dict[key]=Symbol(str(rand.choice(self.param_range[key])))
#                     else:
#                         split_name_dict=str(key).split('_',1)
#                         sym_dict[key]=Symbol(str(rand.randrange(val[0],val[1],val[2]))+split_name_dict[0]+'_0')

                sym_list.append(sym_dict)   
        else: 
            for caseno in range(cases_no):
                sym_dict={}
                temp_dict={}
                for key,val in param_range.items(): 
                        if isinstance(val,list)==True:
                            for elem in val:
                                if isinstance(elem,Expr):
                                    sym_dict[key]=rand.choice(val)
                                else:
                                    split_name_dict=str(key).split('_',1)
                                    sym_dict[key]=Symbol(str(rand.choice(self.param_range[key])),positive=True)
                        else:
                            split_name_dict=str(key).split('_',1)
                            sym_dict[key]=Symbol(str(rand.randrange(val[0],val[1],val[2]))+split_name_dict[0]+'_0',positive=True)

                sym_list.append(sym_dict)   
            
        


        
        return sym_list
#         for num,sym in enumerate(symbols_list):
#             sym_dict={symbols_list[sym]:'a'}

    def generate_cases_data(self,cases_no=5,param_range=None,question_list=None):

        if question_list:
            unique_sym_list=[],
            data_dict_list=[]

            while len(unique_sym_list)<=cases_no:
                data_dict=self.generate_dict(cases_no=1, param_range=param_range )[0]
                case_sys=self.system.subs(data_dict)
                for qs in question_list:
                
                    

                    answer_formula=qs.answer_genrator(case_sys)
                    if not answer_formula in unique_sym_list:
                        unique_sym_list.append(answer_formula)
                        data_dict_list.append(data_dict)
            
            return data_dict_list

        else:
            return self.generate_dict(cases_no=cases_no, param_range=param_range )


    def preview(self,example=False):
        if example:
            path=self._real_example_path
             
        else:
            path=self._scheme_path
            
        with open(f"{path}", "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read())
        image_file.close()

        return IP.display.Image(base64.b64decode(encoded_string))

        
    def prepare_quiz(self,question_list,subs_dict,preview=True):
        self.question_list=question_list
        self.title=str(self)
        self._preview=preview
        self.ds_list=[HarmonicOscillator(self.system.subs(subs)) for subs in subs_dict]
        ds_list=self.ds_list
        self.subs_dict=subs_dict
        question_cat=self.title
        preview_on=self.preview

        subs_dicts_list=self.subs_dict

        cases_no


        question_set=[]
        for case_no in range(len(ds_list)):
            case = 'case'+str(case_no)

            fig=EmbeddedGraphics(self._scheme_path)
            pic=EmbeddedGraphics(self._real_example_path)



            question_string = '\n \n'
            question_string += '<p>'+self.title+'</p>  \n \n'

            question_string+=str(pic)

            question_string +='<p> Based on the model presented, study the dynamics of the considered system. Perform the analysis for the following data.  \(' + '\),\('.join([vlatex(Eq(lhs,rhs,evaluate=False) ) for lhs,rhs in                 self.subs_dict[case_no].items()]) + '\). Compute:</p> \n \n'
            
            question_string+=str(fig)

            for no,question in enumerate(self.question_list):

                counter=no
                question_name=question_cat+"_zestaw_"+case#+str(counter+1001)



                question_string+=str(Paragraph("{q_no}.".format(q_no=str(counter+1))))
                question_string+="<p>"+ (str(question(ds_list[case_no],ds_list[0:case_no]+ds_list[case_no+1:],score=2))) +  " </p> \n"

                question_string+='\n \n'



            question_set.append(str(Question([question_string],id=counter+1000,figure=fig.to_base64_entry()+pic.to_base64_entry(),title=question_name)))


        dump=(Category(name=question_cat,questions=question_set)).to_string()

        with open(question_cat+'_questions.xml', 'w+') as q_file:
            q_file.write(dump)

        if self._preview:
            for no,question in enumerate(self.question_list):

                qs = question(ds_list[0:1], ds_list[1:])

                qs.preview()

                print(str(qs))