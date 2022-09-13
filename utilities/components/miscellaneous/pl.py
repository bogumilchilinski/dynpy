from ..mechanics import *
from pylatex import Document, Package, Command, NewPage, Tabularx, VerticalSpace


class CompanyDataComponent(ReportComponent):
    
    title=" "
    packages=[Package('float')] 
                                
    def append_elements(self):
        
        system = self._system

        

        
        with self.create(Tabularx("l X l",width_argument=NoEscape(r"\textwidth"))) as company_name:


        
            company_data = r"""
\hspace{-8pt} \multirow{5}{*}{\includegraphics[height=2.1cm]{yellowaylogo.jpg}} & \textbf{Yelloway Sp. z o. o.} & \hskip12pt\multirow{5}{*}{\begin{tabular}{r}\footnotesize\bf WYCENA \\[-0.8ex] \footnotesize 01/2022 \\[-0.4ex] \footnotesize\bf DATA \\[-0.8ex] \footnotesize \MakeUppercase{\today} \end{tabular}}\hspace{-6pt} \\
   & & \\
   & Paweł Kamiński & \\
   & +48 667-799-660 & \\
   & pawel.kaminski.yelloway@gmail.com & \\
"""
            self.append(NoEscape(company_data))
            self.append(VerticalSpace(size='3cm'))
