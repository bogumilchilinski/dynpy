from  ..mechanics import *



class ReportComponent(Subsection):

    latex_name = 'subsection'
    packages=[
              Package('standalone'),
              Package('siunitx')
             ]
    
    title='Report generic component'

    def __init__(self, reported_object, title=None, numbering=False, *, label=True, **kwargs):
        """
        Args
        ----
        title: str
            The section title.
        numbering: bool
            Add a number before the section title.
        label: Label or bool or str
            Can set a label manually or use a boolean to set
            preference between automatic or no label
        """



        self.reported_object = reported_object #it's forced by pylatex library
        
        
        if title is None:
            title = self.title
        
        super().__init__(title=title, numbering=numbering, label=label, **kwargs)
        CurrentContainer(self)
        

        self.append_elements()
        
    def append_elements(self):
        pass

    def as_frame(self):
        frame=Frame(title=self.title,options=['allowframebreaks'])
        #frame.packages +(self.packages)
        frame+=(list(self))
        return frame

    @property
    def reported_object(self):

        from ....solvers.linear import ODESystem
        
        if isinstance(self._reported_object, ODESystem):
            return self._reported_object

        else:
            return None
            
    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object=obj

    @property
    def _system(self):
        print('Kod do poprawienia #################### bo stosujesz starą zmienną self._system')
        print('zamień se self._system na self.reported_object')
        
        return self.reported_object


class TitlePageComponent(Environment):
    
    latex_name='titlepage'
    
    def __init__(self, system=None, options=None, arguments=None, start_arguments=None,
                 **kwargs):
        r"""
        Args
        ----
        options: str or list or  `~.Options`
            Options to be added to the ``\begin`` command
        arguments: str or list or `~.Arguments`
            Arguments to be added to the ``\begin`` command
        start_arguments: str or list or `~.Arguments`
            Arguments to be added before the options
        """

        self.system = system
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        
        
        super().__init__(options=options, arguments=arguments, start_arguments=start_arguments,**kwargs)
        
        if self.system is not None:

        
            system = self.system


            
            self.append(NoEscape('\centering'))

            self.append(NoEscape('\\Huge DRGANIA MECHANICZNE \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))

            
            
            if len(system.dvars)==1:
                dof_str = 'JEDNYM STOPNIU SWOBODY'
            else:
                dof_str = 'WIELU STOPNIACH SWOBODY'
                
            if system._dissipative_potential==0 or system._dissipative_potential is None:
                damping_str = 'NIETŁUMIONE'
            else:
                damping_str = 'TŁUMIONE'

           
            self.append(NoEscape(f'\\Large {damping_str} UKŁADY O {dof_str} \n \n'))    
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(NoEscape(f'{system._label} \n \n'))
            
            self.append(Command('vspace',arguments='1cm'))
            
            self.append(Command('MyAuthor'))
            self.append(NoEscape(f'\\par'))
            self.append(Command('vspace',arguments='1cm'))
            
            #self.append(Command('vspace',arguments='1cm'))
            #self.append(NoEscape(f'\\protect\\par'))
            #self.append(NewLine())
            self.append(Command('MyDate'))



