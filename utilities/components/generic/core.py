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




class Component(Section):
    """
    Base class for report components rendered as document sections.

    This class represents a high-level report component that:
    - acts as a LaTeX `Section`,
    - optionally renders itself in Jupyter (Markdown/HTML),
    - aggregates sub-elements returned by the `elements()` method.

    Subclasses should override the `elements()` method to define
    the content of the component.

    Notes
    -----
    * The `elements()` method is expected to return either:
        - a dictionary of named report elements, or
        - None if the component has no body.
    * Keys of the dictionary are informational only.
    * Values may include objects such as:
        - ReportText
        - SympyFormula
        - other pylatex-compatible elements
    * Elements providing `_repr_markdown_` or `_repr_html_` will be rendered directly in Jupyter.

    The component title may be defined statically via the `title`
    class attribute or dynamically by overriding `dynamic_title()`.

    Parameters
    ----------
    reported_object : object
        Domain object that the component reports on.
        Required by the pylatex integration.
    title : str, optional
        Section title. If None, `dynamic_title()` is used.
    numbering : bool, optional
        Whether the section should be numbered.
    label : bool or str, optional
        Section label:
        - True  -> automatic label
        - False -> no label
        - str   -> custom label

    Examples
    --------
    Minimal custom report component:

    >>> from sympy import Symbol, Eq
    >>> from dynpy.utilities.components.generic.core import Component
    >>>
    >>> class MyReportComponent(Component):
    ...     title = "Report generic component"
    ...
    ...     def elements(self):
    ...         dict = {}
    ...   
    ...     comp_dict = {}
    ...
    ...     comp_dict['text1'] = ReportText('header '*10 + '\n\n')
    ...     comp_dict['text2'] = ReportText('body'*50+ '\n\n')
    ...     
    ...     comp_dict['eq'] = SympyFormula(Eq(Symbol('E'),Symbol('mc^2')))
    ...     
    ...     comp_dict['text3'] = ReportText('foot'*10+ '\n\n')
    ...     comp_dict['other comp'] = ReportText('foot'*10+ '\n\n')
    ...      
    ...     return comp_dict

    Extending behavior
    ------------------
    Recommended extension points:
    - `elements()`       – define component content
    - `dynamic_title()` – compute title dynamically
    - `as_frame()`       – convert the component to a Beamer frame
    
    """
    
    
    
    
    numbering = True

    latex_name = "section"
    packages = [Package("standalone"), Package("siunitx")]

    title = "Report generic component"

    def __init__(
        self, reported_object, title=None, numbering=None, *, label=True, **kwargs
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

        self._elems_added = False
        self.append_elements()
        

    def dynamic_title(self):
        return self.title

    def append_elements(self):
        
        # TODO: method should be removed to avoid calling elements() during initialization, but for now it is needed to display title before elements

        elems = self.elements()

        if self._elems_added is False:
            
            self._elems_added = True
            if elems is not None:

                for element in self.elements().values():
                
                    if isinstance(element,(ReportText, SympyFormula)):
                        
                        CurrentContainer
                        element._repr_markdown_()
                    else:    
                        self.append(element)
                        
                
                return True
            else:
                return None
        else:
            
            return None

    def _repr_markdown_(self):


        md_list =[]

        repr_methods = {
            "_repr_markdown_": lambda obj: obj._repr_markdown_(),
            "_repr_html_": lambda obj: obj._repr_html_(),
        }


        # prevents componentent appending on each display call
        empty_sec = Section('Dummy sec')
        CurrentContainer(empty_sec)
        
        for element in self:
            # for method_name, fn in repr_methods.items():
            #     # if hasattr(element, method_name):
            #     #     elem = fn(element)

            #     # else:
            if hasattr(element,"_repr_markdown_"):
                elem = element._repr_markdown_()
                
            elif hasattr(element,"_repr_html_"):
                elem = element._repr_html_()            
            
            else:                
                elem = (str(element))
                
            md_list.append(elem)

        md_str = '# ' + self.dynamic_title() + "\n\n\n" 
        md_str +=  '\n\n'.join(md_list)
        
        CurrentContainer(self)
        
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

