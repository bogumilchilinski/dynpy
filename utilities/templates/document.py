from pylatex import (Document, Package, Command
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
                 documentclass='article',
                 document_options=None,
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=True,
                 page_numbers=True,
                 indent=None,
                 geometry_options=['lmargin=25mm', 'rmargin=25mm',  'top=20mm', 'bmargin=25mm', 'headheight=50mm'],
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

