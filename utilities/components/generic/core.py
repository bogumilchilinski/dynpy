import itertools as itools

import sympy as sym
from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    Package,
)
from pylatex.base_classes.containers import Container

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape  # italic,
from sympy import *

from ...adaptable import *
from ...report import *
from ...report import display as ip_display


def display(obj):

    # ip_display(IPMarkdown('=='*25 +'\n \n \n ==== if u see it let me (@bogumilchilinski) know' ))
    return ip_display(obj)




class Component(Subsection):

    latex_name = "subsection"
    packages = [Package("standalone"), Package("siunitx")]

    title = "Report generic component"

    def __init__(
        self, reported_object, title=None, numbering=False, *, label=True, **kwargs
    ):
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

        self.reported_object = reported_object  # it's forced by pylatex library

        if title is None:
            title = self.dynamic_title()

        super().__init__(title=title, numbering=numbering, label=label, **kwargs)
        CurrentContainer(self)

        if self.title is not None:
             if self.elements() is None:
                 ip_display(IPMarkdown(f"## {self.title}"))

        self.append_elements()

    def dynamic_title(self):
        return self.title

    def append_elements(self):

        elems = self.elements()

        if elems is not None:

            for element in self.elements().values():
            
                if isinstance(element,(ReportText, SympyFormula)):
                    element._repr_markdown_()
                elif isinstance(element,AdaptableDataFrame):
                    element.reported()
                else:    
                    self.append(element)
                    
            return True
        else:
            return None

    def _repr_markdown_(self):


        md_list =[]

        repr_methods = {
            "_repr_markdown_": lambda obj: obj._repr_markdown_(),
            "_repr_html_": lambda obj: obj._repr_html_(),
        }

        
        for element in self.elements().values():
            for method_name, fn in repr_methods.items():
                if hasattr(element, method_name):
                    elem = fn(element)
                    break
            else:
                elem = element.__repr__()
            md_list.append(elem)

        md_str = '# ' + self.dynamic_title() + "\n\n\n" 
        md_str +=  '\n'.join(md_list)
        
        
        return md_str


    def elements(self):
        return None


    def as_frame(self):
        frame = Frame(title=self.title, options=["allowframebreaks"])
        # frame.packages +(self.packages)
        frame += list(self)
        return frame

    @property
    def reported_object(self):

        return self._reported_object

    @reported_object.setter
    def reported_object(self, obj):
        self._reported_object = obj

    @property
    def _system(self):
        print(
            "Kod do poprawienia #################### bo stosujesz starą zmienną self._system"
        )
        print("zamień se self._system na self.reported_object")

        return self.reported_object

