from pylatex import (Document, Package, Command, NewPage, Tabularx
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )
#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)


class CaseTemplate(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=['B. Chiliński, A. Mackojć'],options=['R']),
                  Command('fancyhead', arguments=['Mechanical Vibration, 2021'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
                  #Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
                  #Command('newcommand{\dyplom}', arguments=['Inżynierska']),
                  #Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
                  #Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
                  #Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
                  #Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
                  #Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
                  #Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
                  #Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
                  #Command('newcommand{\\title}', arguments=['Wpisać tytuł pracy dyplomowej po angielsku']),
                  #Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
                  #Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
                  #Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
                  #Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
    ]


# \renewcommand{\headrulewidth}{1pt}
# \renewcommand{\footrulewidth}{1pt}
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )

class ExampleTemplate(Document):
    
    latex_name = 'document'
    packages = [
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  Command('pagestyle', arguments=['fancy']),
                  Command('fancyhf', arguments=['']),
                  Command('fancyhead',  arguments=['B. Chiliński, A. Mackojć, D. Sierociński'],options=['R']),
                  Command('fancyhead', arguments=['Drgania mechaniczne, 2022'],options=['L']),
                  Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
                  
        
#                   Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
#                   Command('newcommand{\dyplom}', arguments=['Inżynierska']),
#                   Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
#                   Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
#                   Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
#                   Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
#                   Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
#                   Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
#                   Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
        
                  
#                   Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
#                   Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
#                   Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
#                   Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
    ]


# \renewcommand{\headrulewidth}{1pt}
# \renewcommand{\footrulewidth}{1pt}
    
    
    def __init__(self,
                 title='Basic title',
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
        self.title=title
        #self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
        self.packages.append(Command('author', arguments=['Bogumił Chiliński, Anna Mackojć, Damian Sierociński']))
        self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
        #self.append(Command('maketitle'))
        self.append(NewPage())
        # tu implementować co tam potrzeba
        
class SzymonTemplate(Document):
    
    latex_name = 'document'
    packages = [
                  Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('microtype'),
                  Package('authoraftertitle'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  Package('titlesec'),
                  Package('fancyhdr'),
                  
        
                  Command('newcommand{\praca}', arguments=['Praca dyplomowa']),
                  Command('newcommand{\dyplom}', arguments=['Inżynierska']),
                  Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']),
                  Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']),
                  Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']),
                  Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']),
                  Command('newcommand{\promotor}', arguments=['Wpisać promotora']),
                  Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']),
                  Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']),
                  Command('newcommand{\\album}', arguments=['303596']),
                  Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']),
                  Command('newcommand{\\rok}', arguments=['Rok składania pracy']),
                  Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']),
                  Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']),
    ]


# \renewcommand{\headrulewidth}{1pt}
# \renewcommand{\footrulewidth}{1pt}
    
    
    def __init__(self,
                 title='Basic title',
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
        self.title=title
        #self.packages.append(Command('title', arguments=[NoEscape(self.title)]))
        self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
        self.packages.append(Command('newcommand{\praca}', arguments=['Praca dyplomowa']))
        self.packages.append(Command('newcommand{\dyplom}', arguments=['Inżynierska']))
        self.packages.append(Command('newcommand{\kierunek}', arguments=['Wpisać kierunek']))
        self.packages.append(Command('newcommand{\specjalnosc}', arguments=['Wpisać specjalność']))
        self.packages.append(Command('newcommand{\\autor}', arguments=['Imię i nazwisko autora']))
        self.packages.append(Command('newcommand{\opiekun}', arguments=['Wpisać opiekuna']))
        self.packages.append(Command('newcommand{\promotor}', arguments=['Wpisać promotora']))
        self.packages.append(Command('newcommand{\konsultant}', arguments=['Wpisać konsultanta']))
        self.packages.append(Command('newcommand{\\tytul}', arguments=['Wpisać tytuł pracy dyplomowej po polsku']))
        self.packages.append(Command('newcommand{\\album}', arguments=['303596']))
        self.packages.append(Command('newcommand{\supervisor}', arguments=['dr inż. Bogumił Chiliński']))
        self.packages.append(Command('newcommand{\\rok}', arguments=['Rok składania pracy']))
        self.packages.append(Command('newcommand{\kluczowe}', arguments=['Słowa kluczowe: Wpisać słowa kluczowe po polsku']))
        self.packages.append(Command('newcommand{\keywords}', arguments=['Keywords: Wpisać słowa kluczowe po angielsku']))
        #self.append(Command('maketitle'))
        self.append(NewPage())
        # tu implementować co tam potrzeba
        
class QuotationTemplate(Document):
    packages = [
                  Package('geometry',options=['left=1in','top=1in','right=1in','bottom=1in']),
                  Package('inputenc',options=['utf8']),
                  Package('hyperref',options=['colorlinks']),
                  Package('graphicx'),
                  Package('tabularx'),
                  Package('multirow'),
                  Package('ragged2e'),
                  Package('hhline'),
                  Package('array'),
                  Package('booktabs'),
                  Package('polski',options=['MeX']),]

    def __init__(self,
                 title='Wycena zlecenia',
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
#         self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
# #         self.append(Command('maketitle'))
#         self.append(NewPage())
        # tu implementować co tam potrzeba
        
class QuotationTemplate(Document):
    packages = [
                  Package('geometry',options=['left=1cm','top=1cm','right=1cm','bottom=1cm']),
                  Package('inputenc',options=['utf8']),
                  Package('hyperref',options=['colorlinks']),
                  Package('graphicx'),
                  Package('tabularx'),
                  Package('multirow'),
                  Package('ragged2e'),
                  Package('hhline'),
                  Package('array'),
                  Package('booktabs'),
                  Package('polski',options=['MeX']),]

    def __init__(self,
                 title='Wycena zlecenia',
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
#         self.packages.append(Command('date', arguments=[NoEscape('\\today')]))
# #         self.append(Command('maketitle'))
#         self.append(NewPage())
        # tu implementować co tam potrzeba
        
class ScheduleTemplate(QuotationTemplate):
    pass
    
class LandscapeScheduleTemplate(ScheduleTemplate):
    def __init__(self,
                 title='Wycena zlecenia',
                 default_filepath='default_filepath',
                 *,
                 documentclass='report',
                 document_options='landscape',
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,#['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
                 data=None):

        super().__init__(
#             title=title,
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )
#         label=self.label
#         self.title=title
        self.preamble.append(Command('hypersetup', 'urlcolor=blue'))
        self.preamble.append(Command('newcolumntype', arguments='R',options='1',extra_arguments=NoEscape('>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}')))
        
class BeamerTemplate(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  #Package('titlesec'),
                  #Package('fancyhdr'),
                  #Command('pagestyle', arguments=['fancy']),
                  Command('author', arguments=['Szymon Kozłowski & Bogumił Chiliński']),
                  #Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
                  #Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
                  #Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
                  #Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']), 
                  Command('usetheme', arguments=['Madrid']),
                  Command('graphicspath', arguments=[NoEscape('{../}')])

            ]
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 title=None,
                 documentclass='beamer',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data, 
            )
        self.title = title
        
        if self.title is not None:
            self.packages.append(Command('title', arguments=[self.title]))
        self.append(Command('frame', arguments=[NoEscape(r'\titlepage')]))

class PosterTemplate(Document):

    latex_name = 'document'
    packages = [
                  Package('microtype'),
                  Package('polski',options=['MeX']),
                  #Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
                  Package('listings'),
                  #Package('titlesec'),
                  #Package('fancyhdr'),
                  #Command('pagestyle', arguments=['fancy']),
                  Command('author', arguments=['Szymon Kozłowski & Bogumił Chiliński']),
                  #Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
                  #Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
                  #Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
                  #Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']), 
                  Command('usetheme', arguments=['Simple']),
                  Command('graphicspath', arguments=[NoEscape('{../}')])

            ]
    
    
    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 title=None,
                 documentclass='tikzposter',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data, 
            )
        self.title = title
        title_with_box = NoEscape(r'\parbox{0.9\linewidth}{\centering ' +self.title+ ' }')
        
        if self.title is not None:
            self.packages.append(Command('title', arguments=[title_with_box]))

        #self.append(Command('frame', arguments=[NoEscape(r'\titlepage')]))
        self.append(Command('maketitle'))