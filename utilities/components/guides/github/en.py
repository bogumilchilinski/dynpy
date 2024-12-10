from  ...mechanics import *
from  ...mechanics import ReportComponent as BaseReportComponent


from . import pl

import pandas as pd
import numpy as np
from sympy import lambdify

from sympy import *
from pandas import *
from sympy.physics.mechanics import dynamicsymbols

from .....solvers.linear import *
from .....dynamics import *

months_list = ['January', 'February', 'March','April','May','June','July','August','September','October','November','December']

from  ...mechanics import display
import datetime


miesiace_list = ['styczeń', 'luty', 'marzec','kwiecień','maj','czerwiec','lipiec','sierpień','wrzesień','październik','listopad','grudzień']

average_temp_list = [-1.9,-0.8,3.2,9.3,14.6,18,20.1,19.5,14.7,9.3,4.8,0.5]

Eg_daily_list_Wh_per_meter2 =[600,1000,3000,3800,4800,5400,5300,4900,3300,1700,700,500]

Eg_daily_list_kWh_per_meter2 = [0.6,1,3,3.8,4.8,5.3,4.9,3.3,1.7,0.7,0.5]

day_length_in_months_hours = [8.3,10.0,11.8,13.9,15.7,16.7,16.3,14.7,12.7,10.7,8.8,7.8]

data_atmospheric_conditions = {'Day length in a month [h]':day_length_in_months_hours,'Daily energy intensity [${kWh/m^2}$]':Eg_daily_list_Wh_per_meter2,'Average temperature [$^{\circ}$C]':average_temp_list}

df = pd.DataFrame(index = months_list,data = data_atmospheric_conditions)

class ReportComponent(BaseReportComponent):



    @property
    def reported_object(self):

        from .....solvers.linear import ODESystem
        
        if isinstance(self._reported_object, ODESystem):
            return self._reported_object

        else:
            return self._reported_object
            
    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj

    @property
    def _system(self):
        print('Kod do poprawienia #################### bo stosujesz starą zmienną self._system')
        print('zamień se self._system na self.reported_object')
        
        return self.reported_object

#git        
        
class CocalcLoginComponent(ReportComponent):
    
    title="Logowanie do CoCalc"


    def append_elements(self):

        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage


        display(ReportText('''Link do CoCalca: [LINK](https://cocalc.com/) '''))

        pic1 = Picture('./dynpy/utilities/components/guides/images/ColacSignIn.png', caption = "Przed Zalogowaniem",width='9cm')

        display(pic1)

        display(ReportText('Teraz należy zalogować się'))

        pic2 = Picture('./dynpy/utilities/components/guides/images/MetLogow.png', caption = "Metody logowania",width='9cm')

        display(pic2)

        pic3 = Picture('./dynpy/utilities/components/guides/images/Zalogowany.png', caption = "Zalogowany",width='9cm')

        display(pic3)

        display(ReportText('''Poprzez tego linka można dołączyć do projektu Ongoing: [ONGOING](https://cocalc.com/app?project-invite=dS62jsbRJvcfj2Mu) '''))
        


#zostawic to
class CocalcFolderComponent(ReportComponent):
    
    title="Tworzenie Folderu i Jupytera"


    def append_elements(self):
        
        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage

        display(ReportText('Tworzenie Folderu i Jupytera'))
        
        pic4 = Picture('./dynpy/utilities/components/guides/images/OtwarcieOngo.png', width='9cm')
        
        display(pic4)
        
        display(ReportText('W wyszukiwarce szukamy "Ongoing" i otwieramy projekt  '))
        
        pic5 = Picture('./dynpy/utilities/components/guides/images/PoWejsciuWOngo.png', width='9cm')
        
        display(pic5)
        
        display(ReportText('Otworzenie folderu UCZESTNICY  '))
        
        
        display(ReportText('Nadajemy nazwę folderu od naszego imieniem i nazwiskiem jak na zalączonym obrazku i klikamy create folder:'))
        
        pic6 = Picture('./dynpy/utilities/components/guides/images/NazwanieFolderu.png', width='9cm')
        
        display(pic6)
        
        
        display(ReportText('Po utworzeniu folderu klikamy przycisk NEW i utworzymy pierwszego Jupytera. W tym celu po nazwaniu Jupytera wybieramy zakladke Jupyter Notebook.'))
        
        
        pic7 = Picture('./dynpy/utilities/components/guides/images/WyborKernela.png', width='9cm')
        
        display(pic7)
        
        display(ReportText('Wybieramy Python 3 (system - wide)  '))
        display(ReportText('Po wybraniu Kernela utworzony zostanie nasz Jupyter.'))
        
        pic8 = Picture('./dynpy/utilities/components/guides/images/JupekUtworzony.png', width='9cm')
        
        display(pic8)


#podstawy raportowania
        



        
code_dynsys_list_str = '''
#komendy importujące biblioteki
import sympy 
from sympy import Symbol

from dynpy.models.mechanics.pendulum import Pendulum

#komenda wywołująca preview
Pendulum().preview()

#komenda wywołująca pomocne informacje dotyczące wahadła
help(Pendulum)

'''


class CocalcDynSysListComponent(ReportComponent):
    
    title="Lista systemów dynamicznych"


    def append_elements(self):

        system = self.reported_object # it's useless in the case of permanent content - it's commented for future usage
        
        display(ReportText(' Możemy również wyświetlić podgląd takiego wahadła używając poniższej funkcji:  '))
        display(Picture('./dynpy/utilities/components/guides/images/preview.jpg'))
        display(ReportText(' Co więcej, używając poniższej komendy help, możemy otrzymać wszelkie informacje na temat konkretnego modelu:  '))
        display(Picture('./dynpy/utilities/components/guides/images/help.jpg'))

        display(ReportText('''Kod wywołania pomocy jest następujący:'''))
        
        display(GuideCode(   code_dynsys_list_str   ))


        display(ReportText('Poniżej znajduje się lista systemów mechanicznych, których można użyć:  '))
        from dynpy.models import mechanics


        moduly=inspect.getmembers(mechanics, inspect.ismodule)


        for nazwa,modul in moduly:

            mods_tmp = [name for name,cls in inspect.getmembers(modul, inspect.isclass) if issubclass(cls,mechanics.principles.ComposedSystem)]
            classes_str=',\n - '.join(mods_tmp)
            if mods_tmp != []:
                display(ReportText( '\\par' + f'-{nazwa}'+ '\\par'))
                display(ReportText(f'- {classes_str}'))


    
#git
    
class IssuePreparationComponent(pl.IssuePreparationComponent):
    
    title="Issue preparation component"


##########################
##SYNCHRONIZACJA GITHUB### 
##########################

#git

class GitSynchroPanelAccessComponent(ReportComponent):

    title="Dostep do panelu synchronizacji Github"

    @property
    def reported_object(self):


        default_data = {'classname':'ReportingModuleIntroComponent',
                       'module':'guide.en.py',
                       'field':'guide or report',
                       'target':'`ODESystem` class',
                       'issue_no':359,
                       }


        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data

        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj
    
    
    def append_elements(self):
#       variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        display(ReportText(f''' Repozytoria Github'a zmieniają się dosyć dynamicznie, a w przypadku, gdy korzysta z niego więcej niż jedna osoba dobrze by było aby każdy miał dostęp do najnowszej wersji. W tym celu przeprowadza się regularnie synchronizację danego rezpozytorium. '''))

        display(ReportText('Jak się dostać do panelu synchronizacji?'))

        display(ReportText(f''' Wybierz projekt, który chcesz zsynchronizować.'''))
        display(Picture('./dynpy/utilities/components/guides/images/1.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('11cm'),caption='Wybór projektu z listy'))

        display(ReportText(f'''Przejdź do foledru, który chcesz zsynchronizować. W naszym przypadku zazwyczaj aktualizujemy moduły.'''))
        display(Picture('./dynpy/utilities/components/guides/images/2.png', position = 'H', height= NoEscape('18cm'), width = NoEscape('18cm'),caption='Wybór folderu do synchronizacji'))

        display(ReportText(f'''Wybierz ikonę VS Code - zostaniesz przeniesiony do nowego okna.'''))
        display(Picture('./dynpy/utilities/components/guides/images/3.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('11cm'),caption='Przycisk przenoszący do VS Code'))

        display(ReportText(f'''Będąc już w panelu VS Code przjedź do zakładki 'Source control' (trzecia od góry), to tu będą przeprowadzane wszystkie operacje.'''))
        display(Picture('./dynpy/utilities/components/guides/images/4.png', position = 'H', height= NoEscape('18cm'), width = NoEscape('18cm'),caption='Source control w panelu VS Code'))

#git

class GitSynchroIntroComponent(GitSynchroPanelAccessComponent):

    title="Wprowadzenie do synchronizacji repozytorium Github"

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here

        display(ReportText('Synchronizacja repozytorium.'))

        display(ReportText('Przypadek 1. : "Commit" po lewej stronie NIE działa'))

        display(ReportText(f'''Jeśli "Commit" nie działa to kliknij od razu "Sync changes" lub kółeczka synchronizacji na górze'''))
        display(Picture('./dynpy/utilities/components/guides/images/10.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='Kółeczko synchronizacji'))
        display(Picture('./dynpy/utilities/components/guides/images/9.png', position = 'H', height= NoEscape('8cm'), width = NoEscape('6cm'),caption='Przycisk, na którym pojawia się "Commit" albo "Sync changes"'))

        display(ReportText('Przypadek 2. : "Commit" po lewej stronie wyświetla się normalnie'))

        display(ReportText(f'''Kliknij "Commit".'''))
        display(Picture('./dynpy/utilities/components/guides/images/5.png', position = 'H', height= NoEscape('8cm'), width = NoEscape('6cm'),caption='Przycisk "Commit"'))

        display(ReportText(f'''Otworzy się opis Commit'a - odkomentuj linijki z modyfikacjami, czyli wszystko od około 4. linii. '''))
        display(Picture('./dynpy/utilities/components/guides/images/6.png', position = 'H', height= NoEscape('13cm'), width = NoEscape('13cm'),caption='''Opis Commit'a '''))

        #display(ObjectCode(example_commit_str))

        display(ReportText(f'''Po odkomentowaniu dopisz tytuł synchronizacji oraz krótkie opisy poszczególnych zmian - jest to ważne w kwestii komunikacji między osobami pracującymi nad projektem. '''))
        display(Picture('./dynpy/utilities/components/guides/images/7.png', position = 'H', height= NoEscape('13cm'), width = NoEscape('11cm'),caption=''' Przykladowy opis Commit'a '''))
        display(ReportText('Przykładowy opisany Commit:'))
        #display(ObjectCode(example_commented_commit_str))
        display(ReportText('Spis przykładowych komenatrzy do Committów:'))
        #display(ObjectCode(example_commit_comments_str))


        display(ReportText(f''' Dla pewności poprawności wprowadzonych zmian warto potwierdzić Commita na chatcie do synchronizacji na Slacku'''))
        display(Picture('./dynpy/utilities/components/guides/images/8_5.png', position = 'H', height= NoEscape('6cm'), width = NoEscape('6cm'),caption='''Kanał do synców '''))

        display(ReportText(f'''Zatwierdź Commit'a klikjąc tick w prawym górnym rogu. Jeśli przycisku nie ma, przełącz karty i spróbuj ponownie.'''))
        display(Picture('./dynpy/utilities/components/guides/images/8.png', position = 'H', height= NoEscape('20cm'), width = NoEscape('20cm'),caption='Zatwierdzanie zmian'))

        display(ReportText(f'''Po zatwierdzeniu, ze strony Github'a takie zmiany wyglądają następująco '''))
        display(Picture('./dynpy/utilities/components/guides/images/7_5.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='''Zatwierdzone zmiany z perspektywy Github'a'''))
        display(ReportText('Wystąpienie komunikatu merge'))

        display(ReportText(f'''Taki komunikat będzie oznacza błąd na wyższym szczeblu - nie ruszaj nic więcej i odezwij się do Pana Chilińskiego lub Pana Sierocińskiego. '''))
        display(Picture('./dynpy/utilities/components/guides/images/11.png', position = 'H', height= NoEscape('11cm'), width = NoEscape('8cm'),caption='Error'))
        
       #git
        
class UsageOfGitHubInterfacesComponent(GitSynchroPanelAccessComponent):

    title="Wprowadzenie do GitHubInterface"

    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']
    
        #implement reporting activieties here
        from dynpy.utilities.creators import GitHubInterface
        display(ReportText('Obsluga rytynowych czynności w dynpy jest wykonanywana przez klasę GitHubInterface'))
        display(ObjectCode(GitHubInterface))
        display(ReportText('Więcej informacji możesz uzyskać w helpie'))
        display(ObjectCode("help(GitHubInterface)"))
        display(ObjectCode(help(GitHubInterface)))


#github

class IssueFeedbackComponent(ReportComponent):

    title="Komentarz zamykajacy issue z komponentami"

    @property
    def reported_object(self):
        from .....dynamics import LagrangesDynamicSystem

        default_data = {'classname':'Component',
                       'module':'dynpy.utilities.components.guides.en',
                       'field':'None',
                       'target':'`ReportComponent` class',
                       'issue_no':359,
                       }

        if isinstance(self._reported_object, dict):
            return {**default_data,**self._reported_object}

        if isinstance(self._reported_object, ReportComponent) or isinstance(self._reported_object, LagrangesDynamicSystem):
            return {**default_data,'classname':self._reported_object.__class__.__name__,'module':self._reported_object.__class__.__module__}
        
        elif isinstance(self._reported_object, str):
            return {**default_data,'classname':self._reported_object}

        elif self._reported_object is None:
            return default_data
        
        else:
            return self._reported_object


    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj


    def append_elements(self):
        #variables provided by `reported_object` arg
        classname = self.reported_object['classname']
        class_module = self.reported_object['module']
        class_field = self.reported_object['field']
        target = self.reported_object['target']

        display(ReportText('The code was checked with the following call: '))
        #display(ObjectCode(git_com_str.format(classname=classname)))
        #display(ObjectCode(git_com_str.format(classname=classname,module=class_module,class_field=class_field)))

class GithubIssueReportComponent(ReportComponent):

    title="Details of GitHub issue"

    def dynamic_title(self):
        
        
        issue = self.reported_object
        if issue==None:
            from github.Issue import Issue
            issue = Issue(requester = 'lsikor', headers = {'title' : 'TEST'}, attributes = {'body':'Przykladowy opis - nie został podany żaden argument', 'title':'TEST', 'number': 999}, completed = False)
        
        return "Details of " + issue.title + "- No: " + str(issue.number)
    
    def append_elements(self):
        
        issue = self.reported_object
        if issue==None:
            from github.Issue import Issue
            issue = Issue(requester = 'lsikor', headers = {'title' : 'TEST'}, attributes = {'body':'Przykladowy opis - nie został podany żaden argument', 'title':'TEST', 'number': 999}, completed = False)
        display(Markdown((f'Issue title: {issue.title} Issue number: {issue.number}')))
        display(ReportText('\\newline'))
        if issue.body is None:
            display(Markdown("No issue description"))
        else:
#             display(Markdown("Issue description: " + issue.body))
            display(IPMarkdown("Issue description: " + issue.body))
        display(ReportText('\\newline'))
        display(ReportText('-'*100))
        display(ReportText('\\newline'))

class UsageOfMeetingCreatorComponent(ReportComponent):
    
    
    title="Github meeting creator class"

    def append_elements(self):
        from .....utilities.creators import MeetingIssueCreator
        display(ReportText('This class is used to create meeting issues through API. The class consists of the following code: \\newline'))
        system = self.reported_object
        display(ObjectCode(MeetingIssueCreator))
        display(ReportText(' \\newline To learn more about this class use help. \\newline'))

git_synchro_intro_call_str=("""
from dynpy.utilities.documents.guides import GithubSynchroGuide
from dynpy.utilities.components.guides.en import *

GithubSynchroGuide();

""")


git_synchro_intro_str =('''Guide consist of bellow components which can be called separatelly:
*GithubSynchroIntroComponent
*GitSynchroPanelAccessComponent
*GitSynchroIntroComponent
*UsageOfGitHubInterfacesComponent
*UsageOfMeetingCreatorComponent
''')

class GithubSynchroIntroComponent(ReportComponent):

    title = "Introduction to basics of GithubSynchronization panel"

    def append_elements(self):

        from dynpy import mechanics


        display(ReportText('This guide concers basics of GithubSynchronization panel.'))
        #display(ReportText('Basic call of it is as follows and runs default dynamic system which is `ForcedSpringMassSystem'))
        display(ObjectCode(git_synchro_intro_call_str))
        
        display(ReportText(git_synchro_intro_str))


cocalc_setup_str=("""
from dynpy.utilities.documents.guides import IntroToCocalcGuide
from dynpy.utilities.components.guides.en import *

IntroToCocalcGuide(None);

""")

cocalcus_str =('''Guide consist of bellow components which can be called separatelly:
*CocalcUsageComponent
*JupyterSetUpComponent
*CocalcFolderComponent
*CocalcDynSysListComponent
*ReportingBasicsComponent
''')

class CocalcUsageComponent(ReportComponent):

    title = "Introduction to usage of Cocalc"

    def append_elements(self):

        from dynpy.models import mechanics


        display(ReportText('This guide concers basics of Cocalc enviroment capabilities that are crucial for explorating `DynPy` library.'))
#         display(ReportText('Basic call of it is as follows and runs default dynamic system which is `ForcedSpringMassSystem'))
        display(ObjectCode(cocalc_setup_str))
        
        display(ReportText(cocalcus_str))
