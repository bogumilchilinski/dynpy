from  ..mechanics import *



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

            
            
            if len(system.q)==1:
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


class ODESystemComponent(ReportComponent):
    
    title="Differential quations"
    @property
    def header_text(self):
        #"Energia potencjalna układu wyrażona jest wzorem:"
        return "The investigated differential equations are as follows"

        
    @property
    def footer_text(self):
        #"Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."
        return "Eigenvalues and modal method is used to solve this type equations"

    def append_elements(self):

        system = self._system
        dyn_sys=system



        display(ReportText(  self.header_text   ))

        display(SympyFormula(  system))

        display(ReportText(  self.footer_text   ))


class VariablesComponent(ReportComponent):
    
    title="Differential quations"
    @property
    def header_text(self):
        #"Energia potencjalna układu wyrażona jest wzorem:"
        return "The investigated differential equations are as follows"

        
    @property
    def footer_text(self):
        #"Zaprezentowana zależność opisuje oddziaływanie potencjalnych pól sił w których znajduje się obiekt."
        return "Eigenvalues and modal method is used to solve this type equations"

    def append_elements(self):

        system = self._system
        dyn_sys=system



        display(ReportText(  self.header_text   ))

        display(SympyFormula(  system.ivar))
        display(SympyFormula(  system.dvar))

        display(ReportText(  self.footer_text   ))
        
#########################olc components from perioua module
    
     
