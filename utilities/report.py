import copy
import datetime as dtime
import inspect
import os
import random
from collections.abc import Iterable

import pandas as pd

# import numpy as np
# import pandas as pd
# import pint
import sympy.physics.mechanics as me
from IPython.display import Latex
from IPython.display import Markdown as IPMarkdown
from IPython.display import display
from pylatex import (
    Alignat,
    Axis,
    Center,
    Command,
    Description,
    Document,
    Eqref,
    Figure,
    FlushLeft,
    FlushRight,
    Label,
    LongTable,
    Marker,
    Math,
    NewLine,
    NewPage,
    Package,
    Plot,
    Quantity,
    Ref,
    Section,
    Subsection,
    Table,
    Tabular,
    TextColor,
    TikZ,
)
from pylatex.base_classes import Environment, Options
from pylatex.base_classes.command import CommandBase
from pylatex.package import Package
from pylatex.section import Chapter
from pylatex.utils import NoEscape, italic
from sympy import (
    Add,
    Derivative,
    Eq,
    Equality,
    Expr,
    Function,
    GreaterThan,
    Heaviside,
    ImmutableMatrix,
    LessThan,
    Matrix,
    Mul,
    Number,
    StrictGreaterThan,
    StrictLessThan,
    Subs,
    Symbol,
    atan,
    cos,
    exp,
    latex,
    sign,
    sin,
    symbols,
    tan,
)
from sympy.core.relational import Relational
from sympy.physics.vector.printing import vlatex, vpprint

from .adaptable import *
from .timeseries import TimeDataFrame, TimeSeries

# import matplotlib.pyplot as plt
# import matplotlib.ticker as plticker


pol_lang_dict = {
    "ą": "a",
    "Ą": "A",
    "ć": "c",
    "Ć": "C",
    "ę": "e",
    "Ę": "E",
    "ł": "l",
    "Ł": "L",
    "ó": "o",
    "Ó": "O",
    "ń": "n",
    "Ń": "N",
    "ś": "s",
    "Ś": "S",
    "ż": "z",
    "Ż": "Z",
    "ź": "z",
    "Ź": "Z",
}


def fix_code(str_to_fix):

    for old_char, new_char in pol_lang_dict.items():
        if old_char in str_to_fix:
            # print(f'old code for {old_char}')
            # print(str_to_fix)

            str_to_fix = str_to_fix.replace(old_char, new_char)

            # print('new code')
            # print(str_to_fix.replace(old_char, new_char))
    return str_to_fix


def plots_no():
    num = 0
    while True:
        yield num
        num += 1


plots_no_gen = plots_no()


class RandomDescription:
    def __init__(self, *args, random_or_order=None, **kwargs):

        self.pool_list = args
        self.random_or_order = random_or_order

    def to_string(self):

        if hasattr(self.random_or_order, "__contains__"):
            self.order = list(random_or_order)
        else:
            self.selected_sentences = [random.choice(pool) for pool in self.pool_list]

        return " ".join(self.selected_sentences)

    def __str__(self):

        return self.to_string()

    def __repr__(self):
        return self.__class__.__name__ + ": " + self.to_string()


class ReportModule:
    r"""
    Basic class for maintaining global options of a report module. It provides methods for setting options common with every class inheriting from ReportModule instance.

    Arguments
    =========
        container: obj
            Python's build-in list or Pylatex container object such as Document() or Section().
        path: str
            Path for saving plots.
    Methods
    =======
        set_container(cls,container=None):
            Classmethod; sets Pylatex container object such as Document() or Section(). None by default.
        set_caption(cls,caption=''):
            Classmethod; sets caption for images plotted under current instance. Empty string by default.
        set_directory
            Classmethod; sets directory for saving generated plots. Directory must already exist. Default path: './SDA_results'.
        set_units_dict
            Classmethod; sets dictionary which matches symbols to their units. Empty dict by default.
    Example
    =======
        >>>from modules.utilities.report import ReportModule
        >>>from pylatex import Document
        >>>from sympy import Symbol
        >>>import pint

        >>>ureg = pint.UnitRegistry()
        >>>m=Symbol('m')
        >>>doc=Document
        >>>unit_dict={m:ureg.kilogram}

        >>>RM=ReportModule()

        >>>RM.set_container(doc)

        >>>RM.set_caption('This is caption.')

        >>>RM.set_directory('./my_directory')

        >>>RM.set_units_dict(unit_dict)
    """

    _cls_container = []  # Section('Dummy')
    cls_container = []  # Section('Dummy')
    _container = []  # Section('Dummy')
    cls_path = "."
    _caption = "Figure describes the numerical data"
    _label = "fig:markerIsMissing"
    _units = {}
    _autoreport = False
    # _frame = TimeDataFrame()
    _frame = None  # pd.DataFrame()
    _list = []
    _subplot = False
    _hold = False
    _out_formatter = None  # BaseFrameFormatter  # lambda data: data
    _height = NoEscape(r"6cm")

    @classmethod
    def set_output_formatter(cls, formatter=None):
        # def set_output_formatter(cls, formatter=BaseFrameFormatter):
        cls._out_formatter = formatter
        return cls

    @classmethod
    def set_container(cls, container=None):
        cls.cls_container = container
        return cls

    @classmethod
    def set_caption(cls, caption=""):
        cls._caption = caption
        return cls

    @classmethod
    def set_reporting_mode(cls, mode=True):
        cls._autoreporting = mode

        return cls

    @classmethod
    def set_ploting_mode(cls, subplots=False):
        cls._subplot = subplots

        return cls

    @classmethod
    def set_plot_height(cls, height=None):
        cls._height = height

        return cls

    def _get_str_key(self):
        return self.dumps()

    @classmethod
    def _reset_storage(cls, *args, **kwargs):

        cls._storage = {}
        cls._dict = {}
        cls._frame = TimeDataFrame()
        cls._list = []
        cls._subplot = False

        return cls

    @classmethod
    def set_directory(cls, path="./SDA_results"):

        cls.cls_path = path
        return cls

    @classmethod
    def set_units_dict(cls, units={}):

        cls._units = units
        return cls

    #     @property
    #     def _cls_container(self):
    #         return self.__class__._cls_container

    #     @property
    #     def _container(self):
    #         return self.__class__._cls_container

    def __init__(
        self, container=None, path=None, autoreporting=False, output_formatter=None
    ):
        if container:
            self._container = container
        else:
            self._container = type(self).cls_container

        if path:
            self._path = path
        else:
            self._path = type(self).cls_path

        self._autoreporting = autoreporting

        self._storage = None
        self._story_point = []
        self._frame = pd.DataFrame()

        self._last_result = None
        self._list = []

        if output_formatter:

            self._out_format = output_formatter
        else:
            self._out_format = self.__class__._out_formatter

        # print(f'Report module init - formatter is {self._out_format}')

    def _apply_formatter(self, data):

        # print(type(self._out_format))
        if (
            self._out_format is BaseFrameFormatter,
            PivotFrameSummary,
            FFTFrameFormatter,
            PivotPlotFrameSummary,
        ):
            # print('Base frmatter is working')

            # print('data.index', data.index)
            # print('data.index.name', data.index.name)

            result = self._out_format(data)()

            #             #print('#'*100)
            #             display(result)
            #             #print(result.index.name)
            if not result.index.name:
                result.index.name = ""

            return result
        else:
            # print('callable is working')
            return self._out_format(data)

    @property
    def frame(self):

        if (self._frame) is not None:

            time_frame = TimeDataFrame(self._frame)

        else:

            time_frame = TimeDataFrame(self.__class__._frame)

        if time_frame.columns != [] and time_frame.columns != pd.Index([]):

            time_frame.columns = pd.MultiIndex.from_tuples(time_frame.columns)
            time_frame.columns.set_names(
                ["model", "parameter", "coordinate"], inplace=True
            )

        return time_frame

    def set_frame(self, key, frame):
        # print('self frame is modified')
        self._frame[key] = frame

        return frame

    def set_class_frame(self, key, frame):
        # print('class frame is modified')

        self.__class__._frame[key] = frame

        return frame

    def clear_frame(self, obj=True):

        if not self.__class__._hold:

            if obj and (self._frame is not None):
                # print('self cleaner')
                self._frame = type(self._frame)()
                # self._frame = None
            if (not obj) and (self.__class__._frame is not None):
                # print('class cleaner')
                self.__class__._frame = type(self.__class__._frame)()
            # self.__class__._frame =  None

        return None

    def __str__(self):
        return self._container.__str__()

    def __repr__(self):
        return self._container.__repr__()

    def reported(self):
        self.cls_container.append(self)

        return copy.deepcopy(self)

    def _report(self):
        self.cls_container.append(self)


class CurrentContainer:
    """
    A container class that initializes and configures various components for handling LaTeX reports,
    visualizations, and symbolic mathematics. This class serves as a centralized configuration manager
    for different utilities within a reporting or symbolic computation framework.

    Attributes:
        mdpi (bool): A class-level attribute to define whether MDPI-specific configurations are applied.
                     Defaults to False.

    Methods:
        __init__(container): Initializes the container and sets up the configurations for various components.
        set_filename_prefix(prefix): Sets a prefix for file names generated by TikZPlot and returns a copy of the instance.
        set_descriptions_registry(symbols_description_dict): Configures a registry to store and manage symbol descriptions.
    """

    mdpi = False

    def __init__(self, container: object):
        """
        Initializes the CurrentContainer instance with the given container.

        Args:
            container (object): The main container object to be used for managing various components. This object
                                is expected to provide functionality for report generation, visualization, and symbolic
                                computation.

        Sets up the following components:
            - ReportText: For generating textual content in reports.
            - Picture: For managing visualizations.
            - SympyFormula: For handling symbolic formulas.
            - Markdown: For generating Markdown content.
            - LatexDataFrame: For handling dataframes in LaTeX format with picture mode enabled.
            - ObjectCode: For managing code representations.
            - Block, AlertBlock, ExampleBlock: For creating different types of blocks in LaTeX documents.
            - SymbolsDescription: For managing symbolic descriptions.

        Example:
            >>> container = MyContainer()  # Replace with an actual container class instance
            >>> current_container = CurrentContainer(container)
        """

        from .adaptable import LatexDataFrame

        self._container = container

        ReportText.set_container(self._container)
        Picture.set_container(self._container)
        SympyFormula.set_container(self._container)

        Markdown.set_container(self._container)

        LatexDataFrame.set_container(self._container)
        LatexDataFrame.set_picture_mode(True)

        ObjectCode.set_container(self._container)

        Picture.set_container(self._container)

        Block.set_container(self._container)
        AlertBlock.set_container(self._container)
        ExampleBlock.set_container(self._container)

        SymbolsDescription.set_container(self._container)

    def set_filename_prefix(self, prefix: str) -> "CurrentContainer":
        """
        Sets a prefix for filenames generated by TikZPlot and returns a copy of the CurrentContainer instance.

        Args:
            prefix (str): The prefix to be used for filenames in TikZPlot.

        Returns:
            CurrentContainer: A copy of the instance with the updated filename prefix.

        Example:
            >>> current_container = CurrentContainer(container)
            >>> updated_container = current_container.set_filename_prefix("report_")
        """

        TikZPlot._prefix = prefix

        return copy.copy(self)

    def set_descriptions_registry(self, symbols_description_dict: dict) -> None:
        """
        Configures the DescriptionsRegistry with a given dictionary of symbol descriptions.

        Args:
            symbols_description_dict (dict): A dictionary where keys are symbols (e.g., strings or objects) and values
                                             are their respective descriptions (e.g., strings).

        Functionality:
            - Resets the existing descriptions registry.
            - Sets new descriptions using the provided dictionary.

        Example:
            >>> descriptions = {
            ...     "x": "Position in meters",
            ...     "v": "Velocity in m/s"
            ... }
            >>> current_container.set_descriptions_registry(descriptions)
        """

        DescriptionsRegistry().reset_registry()
        DescriptionsRegistry().set_descriptions(symbols_description_dict)


class DateCommand(CommandBase):
    """
    A class representing a custom date LaTeX command.

    This class represents a custom LaTeX command named
    ``DateCommand``.
    """

    _latex_name = "date"


# class AbstractFrameFormatter(AdaptableDataFrame):
#     _applying_func=lambda x: (x*100)

#     @classmethod
#     def _apply_func(cls,func=None,**kwargs):

#         if func:
#             ops_to_apply=func
#         elif cls._applying_func:
#             ops_to_apply=cls._applying_func
#         else:
#             ops_to_apply=lambda obj: obj

#         return ops_to_apply

#     def __new__(cls,data=None, index=None, columns=None, dtype=None, copy=None,**kwargs):

#         ops_to_apply=cls._apply_func(**kwargs)

#         return ops_to_apply(BasicFormattedFrame(data=data, index=index, columns=columns, dtype=dtype, copy=copy))

# class AbstractSeriesFormatted(AdaptableSeries):
#     _applying_func=lambda x: (x*100)

#     @classmethod
#     def _apply_func(cls,func=None,**kwargs):

#         if func:
#             ops_to_apply=func
#         elif cls._applying_func:
#             ops_to_apply=cls._applying_func
#         else:
#             ops_to_apply=lambda obj: obj

#         return ops_to_apply

#     # def __new__(cls,data=None, index=None, columns=None, dtype=None, copy=None,**kwargs):

#     #     ops_to_apply=cls._apply_func(**kwargs)

#     #     return ops_to_apply(BasicFormattedSeries(data=data, index=index, columns=columns, dtype=dtype, copy=copy))


#         if self.__class__._subplot==False:
#             for coord, data in data_dict.items():

#                 data.plot()
#                 plt.ylabel(coord)

#                 filepath=f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                 ########### for tikz
#                 #ndp=DataPlot('wykres_nowy',position='H',preview=False)

#                 #it should be replaced with data.rename
#     #             #print(data)

#                 data.columns=[type(self)._formatter(label) for label in data.columns ]
#                 #print(type(self)._units)
#                 y_unit_str=f'{(type(self)._units[coord]):Lx}'.replace('[]','')

#                 ndp=data.to_standalone_figure(filepath,colors_list=colors_list,height=NoEscape(r'7cm'),width=NoEscape(r'12cm'),y_axis_description=NoEscape(f',ylabel=${vlatex(coord)}$,y unit={y_unit_str} ,x unit=\si{{\second}}'),legend_pos=legend_pos+','+f'legend columns= {legend_columns}' )
#                 #ndp.add_data_plot(filename=f'{self._path}/{self.__class__.__name__}_data_{next(plots_no_gen)}.png',width='11cm')

#                 ########### for tikz
#                 #ndp.append(data.to_pylatex_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))

#                 ndp.add_caption(NoEscape(f'{type(self)._caption}'))

#                 plt.show()

#                 #print('marker - plot')
#                 #print(self.last_marker)
#                 #print(type(self)._last_marker)
#             #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#                 ndp.append(Label(type(self)._last_marker))

#                 if analysis:
#                     analysis._container.append(ndp)
#                 else:
#                     filepath=f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                     #latex_code=(TimeDataFrame(data).to_tikz_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))
#                     self._container.append(ndp)

#         else:
#             for coord, data in data_dict.items():
#                 data.plot(subplots=self.__class__._subplot,ylabel=coord)
#                 filepath=f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                 ########### for tikz
#                 #ndp=DataPlot('wykres_nowy',position='H',preview=False)

#                 #it should be replaced with data.rename
#     #             #print(data)

#                 data.columns=[type(self)._formatter(label) for label in data.columns ]
#                 #print(type(self)._units)
#                 y_unit_str=f'{(type(self)._units[coord]):Lx}'.replace('[]','')

#                 ndp=data.to_standalone_figure(filepath,subplots=self.__class__._subplot,colors_list=colors_list,height=NoEscape(r'6cm'),width=NoEscape(r'0.9\textwidth'),y_axis_description=NoEscape(f',ylabel=${vlatex(coord)}$,y unit={y_unit_str} ,x unit=\si{{\second}}'),legend_pos=legend_pos+','+f'legend columns= {legend_columns}' )
#                 #ndp.add_data_plot(filename=f'{self._path}/{self.__class__.__name__}_data_{next(plots_no_gen)}.png',width='11cm')

#                 ########### for tikz
#                 #ndp.append(data.to_pylatex_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))

#                 ndp.add_caption(NoEscape(f'{type(self)._caption}'))

#                 plt.show()

#                 #print('marker - plot')
#                 #print(self.last_marker)
#                 #print(type(self)._last_marker)
#             #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#                 ndp.append(Label(type(self)._last_marker))

#                 if analysis:
#                     analysis._container.append(ndp)
#                 else:
#                     filepath=f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                     #latex_code=(TimeDataFrame(data).to_tikz_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))
#                     self._container.append(ndp)
#             return ndp


class ReportEntry:
    r"""
    This class creates a report section with a title provided by a user.

    Arguments
    =========
    block_title: str
        String to create a block title.

    Methods
    =======

    Example
    =======
    """

    def __init__(self, block_title, **kwargs):
        self._block_title = block_title

    def __call__(self, analysis, **kwargs):
        sec = Section(self._block_title)

        analysis._container.append(sec)

        return sec


class ReportText(ReportModule):
    r"""
    This class appends a user defined text to the existing document container.

    Arguments
    =========
    text: str
        String that will be appended to the defined container.
    key_dict: dict
        Dictionary containing entries for string format method.

    Methods
    =======

    Example
    =======
    """

    _conjuction_dict = None
    _color = None
    _alignment = None

    @classmethod
    def set_text_color(cls, color=None):

        cls._color = color

        return cls

    @classmethod
    def set_text_alignment(cls, align=None):

        if align == "right":

            cls._alignment = FlushRight

        elif align == "center":

            cls._alignment = Center

        elif align == "left":

            cls._alignment = FlushLeft

        else:
            cls._alignment = align

        return cls

    def __init__(self, text=None, alignment=None, key_dict={}, **kwargs):

        self._text = "Figures {first_marker}-{last_marker}"

        if text:
            self._text = text

        if alignment:
            self._alignment = alignment

        try:
            self._text = self._text.format(**DataStorage._dict)

        except:
            pass
            # print('.w')
        finally:
            self._text = self._text

        super().__init__()

        self._text_formatted = self._text

    def __call__(self, analysis):

        # print(self._text)

        analysis._container.append(NoEscape(self._text))

        return self._text

    def __str__(self):

        return self._text

    def __repr__(self):

        # return (self._text)
        return self._text

    @property
    def _alignment_env(self):

        align = self._alignment

        if align == "right":

            _alignment = FlushRight

        elif align == "center":

            _alignment = Center

        elif align == "left":

            _alignment = FlushLeft

        else:
            _alignment = align

        return _alignment

    def _repr_markdown_(self):
        if self.__class__._color:
            print(f"Color text {self.__class__._color}")
        return self.reported()._text

    def reported(self):

        formatted_text = self._text_formatted

        if self._conjuction_dict is not None and isinstance(
            self._conjuction_dict, dict
        ):
            for dict_key in self._conjuction_dict.keys():
                formatted_text = formatted_text.replace(
                    dict_key, self._conjuction_dict[dict_key]
                )

        if self.__class__._color:

            text_to_add = TextColor(self.__class__._color, NoEscape(formatted_text))
        else:
            text_to_add = NoEscape(formatted_text)

        #         text_to_add =NoEscape(self._text_formatted)

        if self._alignment_env:
            align_env = self._alignment_env()
            align_env.append(text_to_add)

            text_to_add = align_env

        self.cls_container.append(text_to_add)

        return copy.copy(self)

    def with_ai(text, api_key):
        """
        Method which helps with creating the summary of given text for preparing a report using the generative artificial intelligence model.

         -----
         Example:

         >>>display(AiInterfaceOutput(ReportText.summary_with_ai(“text”, “api_key”

         where:
         text - text given for the example
         api_key - specific API key
        """

        from dynpy.utilities.tools import AiInterface

        ai_chat = AiInterface(api_key)

        prompt_for_summary = "Summarize this text: " + text

        response = ai_chat.chat(prompt_for_summary)

        return response

    def call_code(self):
        # Declaration of command vector used in call_code method
        self.command_vector = []

        if self._text != None:
            self.command_vector.append("text = ''' {} '''".format(self._text))
        if self._alignment != None:
            self.command_vector.append("alignment = '''{}  '''".format(self._alignment))

        vector = "ReportText("

        for element in self.command_vector:
            vector += str(element) + ", "

        return vector[:-2] + ")"


import re
from typing import List, Tuple


class LatexSanitizer:
    # ---------- 1. twarde podmiany Unicode --------------------------- #
    _UNICODE = {
        "–": "-",
        "—": "---",
        "“": "``",
        "”": "''",
        "‘": "`",
        "’": "'",
        "…": r"\ldots{}",
        "→": r"$\rightarrow$",
        "✓": r"$\checkmark$",
        "•": r"$\bullet$",
        "€": r"\euro{}",
    }

    # ---------- 2. uciekanie znaków w TEKŚCIE ------------------------- #
    _ESC_TEXT = str.maketrans(
        {
            "#": r"\#",
            "$": r"\$",
            "%": r"\%",
            "&": r"\&",
            "_": r"\_",
            "{": r"\{",
            "}": r"\}",
            "~": r"\textasciitilde{}",
            "^": r"\^{}",
        }
    )

    # ---------- 3. regex na „fragment kodu” --------------------------- #
    _RE_CODE = re.compile(
        r"""(?<!\\)(               # nie poprzedzony ukośnikiem
            __\w+__\.py |          #  __init__.py, __main__.py…
            [A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*\(\)? |   # func(), pkg.mod()
            [A-Za-z_]\w*(?:\.[A-Za-z_]\w*)+          # pkg.mod.attr
        )""",
        re.VERBOSE,
    )

    # ---------- 4. uciekanie MINIMALNE wewnątrz \texttt{…} ------------ #
    _ESC_CODE = str.maketrans(
        {
            "_": r"\_",
            "#": r"\#",
            "%": r"\%",
            "&": r"\&",
        }
    )

    # ================================================================= #
    @classmethod
    def sanitize(cls, text: str) -> str:
        segments: List[Tuple[bool, str]] = []  # (is_code, fragment)

        # --- A. rozdziel kod / tekst --------------------------------- #
        last = 0
        for m in cls._RE_CODE.finditer(text):
            if m.start() > last:
                segments.append((False, text[last : m.start()]))
            segments.append((True, m.group(0)))
            last = m.end()
        segments.append((False, text[last:]))

        # --- B. przetwarzaj osobno ----------------------------------- #
        out_parts: List[str] = []
        for is_code, frag in segments:
            if is_code:
                out_parts.append(frag.translate(cls._ESC_CODE))
            else:
                # 1. Unicode
                for orig, repl in cls._UNICODE.items():
                    frag = frag.replace(orig, repl)
                # 2. zwykłe znaki LaTeX
                frag = frag.translate(cls._ESC_TEXT)
                out_parts.append(frag)

        return "".join(out_parts)


class Aspect(Section, ReportModule):
    packages = [  # Package('float'),
        # Package('fvextra'),
        # Command('floatname', arguments=['listing'],extra_arguments=['Listing']),
    ]

    latex_name = "subparagraph"
    content_separator = "\n"


class Verbatim(Environment):
    pass


class Minted(Environment):
    packages = [Package("minted")]
    content_separator = "\n"


class LstListing(Environment):
    packages = [
        Package("listings"),
        Package("float"),
        Package("fvextra"),
        Package("xcolor"),
        Command("floatname", arguments=["listing"], extra_arguments=["Listing"]),
        Command(
            "definecolor{codegreen}", arguments=["rgb"], extra_arguments=["0,0.6,0"]
        ),
        Command(
            "definecolor{codegray}", arguments=["rgb"], extra_arguments=["0.5,0.5,0.5"]
        ),
        Command(
            "definecolor{codepurple}",
            arguments=["rgb"],
            extra_arguments=["0.58,0,0.82"],
        ),
        Command(
            "definecolor{backcolour}", arguments=["rgb"], extra_arguments=["0,0,0"]
        ),
        Command(
            "lstdefinestyle{dynPyStyle}",
            arguments=[
                NoEscape(
                    r"commentstyle=\color{codegreen},keywordstyle=\color{magenta},numberstyle=\tiny\color{codegray},stringstyle=\color{codepurple}, basicstyle=\ttfamily\footnotesize,breakatwhitespace=false,breaklines=true,captionpos=b,keepspaces=true,numbers=left,numbersep=5pt,showspaces=false,showstringspaces=false,showtabs=false,tabsize=2,language=Python"
                )
            ],
        ),
        # Command('lstset', arguments = [NoEscape('style=dynPyStyle,language=Python')]),
    ]

class BibliographyManager(Environment):

    r'''
    Description how to use it:

    from dynpy.utilities.report import *
    biblio_mngr=BiblographyManager(arguments='biblio.bib')

    biblio_entries=NoEscape(
    r"""
    @misc{DynPi,
    author={GitHub},
    title="bogumilchilinski/dynpy",
    url={https://github.com/bogumilchilinski/dynpy},
    note="Accessed:2024-05-04"
    }
    """)

    biblio_mngr.append(biblio_entries)

    '''

    latex_name = 'filecontents'
    content_separator = "\n"

class PyVerbatim(Environment):
    packages = [
        Package("pythontex"),
        Package("fvextra"),
        Command(
            "fvset",
            arguments=[NoEscape(r"breaklines, commandchars=\\\{\},numbers=left")],
        ),
    ]
    content_separator = "\n"
    # inspect.getsource(system.__class__)


class Listing(Figure, ReportModule):
    packages = [
        Package("float"),
        Package("fvextra"),
        Command("floatname", arguments=["listing"], extra_arguments=["Listing"]),
    ]

    content_separator = "\n"

    _default_width = NoEscape("0.8\\textwidth")
    _position = "H"

    def __init__(
        self,
        verbatim=None,
        position=None,
        caption=None,
        width=None,
        height=None,
        marker=None,
        **kwargs,
    ):
        """
        Args
        ----
        position: str
            Define the positioning of a floating environment, for instance
            ``'h'``. See the references for more information.
        width: str
            Documentation entry

        References
        ----------
            * https://www.sharelatex.com/learn/Positioning_of_Figures
        """

        self.verbatim = verbatim
        self.caption = caption

        if width is not None:
            self.width = width
        else:
            self.width = self._default_width

        self.height = height

        self.marker = marker

        if position is not None:
            self._position = position

        super().__init__(position=self._position, **kwargs)

        if self.verbatim is not None:
            self.append(self.verbatim)

        if self.caption is not None:
            self.add_caption(NoEscape(self.caption))

        self.marker = marker


#         if self.marker is not None:
#             self.append(Label(self.marker))
#         else:
#             auto_mrk = AutoMarker(self).marker
#             marker = auto_mrk

#             self.marker = marker

#             self.append(Label(marker))


class BeamerTemplate(Document):

    latex_name = "document"
    packages = [
        Package("microtype"),
        Package("polski", options=["MeX"]),
        # Package('geometry',options=['lmargin=25mm', 'rmargin=25mm',  'top=30mm', 'bmargin=25mm', 'headheight=50mm']),
        Package("listings"),
        # Package('titlesec'),
        Package("fancyhdr"),
        Command("pagestyle", arguments=["fancy"]),
        Command("fancyhf", arguments=[""]),
        # Command('fancyhead', arguments=[NoEscape('\includegraphics[height=1.5cm]{./images/logoPOWER.jpg}')],options=['C']),
        # Command('fancyfoot', arguments=['BCh&KT'],options=['R']),
        # Command('fancyfoot', arguments=['Practical Python, 2022'],options=['L']),
        # Command('fancyfoot', arguments=[NoEscape('\\thepage')],options=['C']),
        Command("usetheme", arguments=["Madrid"]),
    ]

    def __init__(
        self,
        default_filepath="default_filepath",
        *,
        title=None,
        documentclass="beamer",
        document_options=None,
        fontenc="T1",
        inputenc="utf8",
        font_size="normalsize",
        lmodern=False,
        textcomp=True,
        microtype=True,
        page_numbers=True,
        indent=None,
        geometry_options=None,
        data=None,
    ):

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
            self.packages.append(Command("title", arguments=[self.title]))
        self.append(Command("frame", arguments=[NoEscape(r"\titlepage")]))


class Frame(Environment, ReportModule):
    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self, title=None, options=None, arguments=None, start_arguments=None, **kwargs
    ):
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

        self.title = title
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.title is not None:
            self.append(Command("frametitle", arguments=[self.title]))


class TitlePage(Environment, ReportModule):
    pass


class Frame(Environment, ReportModule):
    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self, title=None, options=None, arguments=None, start_arguments=None, **kwargs
    ):
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

        self.title = title
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.title is not None:
            self.append(Command("frametitle", arguments=[self.title]))


class CodeFrame(Environment, ReportModule):
    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self, title=None, options=None, arguments=None, start_arguments=None, **kwargs
    ):
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

        _latex_name = "frame"

        self.title = title
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments
        self._latex_name = _latex_name

        if options is not None:
            options = ["containsverbatim"] + options

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.title is not None:
            self.append(Command("frametitle", arguments=[self.title]))


class Markdown(Environment, ReportModule):
    packages = [
        # Package('markdown'),
        # Package('hyperref'),
    ]

    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self,
        markdown=None,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
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

        self.markdown = markdown
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.markdown is not None:
            self.append(NoEscape(markdown))

    def dumps(self):

        import pypandoc as ppd

        return NoEscape(ppd.convert_text(self.markdown, to="latex", format="md"))

    def _repr_markdown_(self):
        return self.reported().markdown

    def reported(self):

        import pypandoc as ppd

        latex_code = ppd.convert_text(self.markdown, to="latex", format="md")

        if "\\begin{verbatim}" in latex_code:

            env_name = ObjectCode._latex_name

            latex_code = latex_code.replace(
                "\\begin{verbatim}", f"\\begin{{{env_name}}} \n"
            ).replace("\\end{verbatim}", f"\n \\end{{{env_name}}}")

            latex_code = fix_code(latex_code)

            self.cls_container.packages |= ObjectCode.packages

        self.cls_container.packages |= self.packages

        self.cls_container.append(NoEscape(latex_code))

        return copy.copy(self)

    def call_code(self):
        # Declaration of command vector used in call_code method
        self.command_vector = []

        if self.markdown != None:
            self.command_vector.append("markdown = {}".format(self.markdown))
        if self.options != None:
            self.command_vector.append("options = {}".format(self.options))
        if self.arguments != None:
            self.command_vector.append("arguments = {}".format(self.arguments))
        if self.start_arguments != None:
            self.command_vector.append(
                "start_arguments = {}".format(self.start_arguments)
            )

        vector = "Markdown("

        for element in self.command_vector:
            vector += str(element) + ", "

        print(vector[:-2] + ")")

    @classmethod
    def set_mdpi(cls):
        cls.packages = [
            #                 Package('markdown'),
        ]


class Block(Environment, ReportModule):

    packages = [
        Package("markdown"),
    ]
    content_separator = "\n"

    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self,
        header=None,
        data=None,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
        r"""
        Args
        ----
        arguments: str or `list`
            The arguments for the container command
        options: str, list or `~.Options`
            The options for the preamble command
        data: str or `~.LatexObject`
            The data to place inside the preamble command
        """

        self.data = data
        self.options = options
        self.arguments = arguments
        self.header = header

        if type(self.header) is not list:
            header_list = [self.header]
        else:
            header_list = self.header

        super().__init__(
            options=options,
            arguments=header_list,
            start_arguments=start_arguments,
            **kwargs,
        )

        if self.data is not None:
            self.append(data)

    def _repr_markdown_(self):

        self.cls_container.append(self)

        return self.dumps()

    def reported(self):

        # latex_code=ppd.convert_text(self.markdown,to='latex',format='md')

        self.cls_container.append(self)

        return copy.deepcopy(self)


class AlertBlock(Block):
    pass


class ExampleBlock(Block):
    pass


class Picture(Figure, ReportModule):
    """
    A class that represents a figure environment in a LaTeX document.

    This class is used to insert and customize an image into the document, with options for caption, size, and position.

    Attributes:
        separate_paragraph (bool): Whether the float is positioned in a separate paragraph.
        _default_width (str): Default width of the image.
        packages (List[Union[Package, Command]]): Required LaTeX packages and commands for the Picture class.
        _position (str): Default position for the floating environment.
        _preview_default_size (str): Default size for previewing the image.

    Exemplary Usage:
        >>> Picture('./assets/exemplary-image.png', caption='Exemplary caption.')

    Costomised exemplary usage:
        >>> Picture(
        >>>     image='./assets/exemplary-image.png',
        >>>     caption='Exemplary caption.',
        >>>     position='h',
        >>>     width='0.5\\textwidth',
        >>>     height='10cm',
        >>>     marker='fig:example'
        >>> )
    """

    separate_paragraph = True
    _default_width = NoEscape("0.8\\textwidth")
    packages = [
        Package("float"),
        Package("graphicx"),
        Command("graphicspath{{../}}"),
    ]

    _position = None
    _preview_default_size = 600 # in pixels

    @classmethod
    def _settle_dynpy(cls):
        from ..dynamics import LagrangesDynamicSystem

        if os.path.exists("./dynpy"):

            return False

        elif not os.path.exists("./._dynpy_env/dynpy"):
            LagrangesDynamicSystem._settle_dynpy()
            return True
        else:
            return False

    #     @classmethod
    #     def _as_picture(cls, position=None, caption=None,width=None,height=None,marker=None, **kwargs):

    #         pic_path = "./._dynpy_env" + cls._default_folder_path[1:] + cls.scheme_name
    #         cls._settle_dynpy()
    #         return Picture(pic_path, position=position, caption=caption, width=width, height=height, marker=marker, **kwargs)

    @classmethod
    def set_default_position(cls, position: str = "") -> None:
        cls._position = position

    @classmethod
    def set_preview_default_size(cls, size: str) -> "Picture":
        cls._preview_default_size = size
        return cls

    _latex_name = "figure"

    def __init__(
        self,
        image=None,
        position="H",
        caption=None,
        width=None,
        height=None,
        marker=None,
        **kwargs,
    ):
        """
        Initialize the Picture class.

        Args:
            image (Optional[str]): Picture or path to the image file.
            position (Optional[str]): LaTeX position specifier for the figure.
            caption (Optional[str]): Caption for the image.
            width (Optional[str]): Width of the image (e.g., '0.5\\textwidth').
            height (Optional[str]): Height of the image.
            marker (Optional[str]): Marker for referencing the figure in LaTeX.

        References:
            https://www.sharelatex.com/learn/Positioning_of_Figures
        """

        if isinstance(image, Picture):
            self.image = image.image
        elif isinstance(image, str):
            self.image = image
        else:
            self.image = None

        self.caption = caption
        self.preview_size = None

        if width is not None:
            self.width = width
        else:
            self.width = self._default_width

        self.height = height

        self.marker = marker

        if position is not None:
            self._position = position

        super().__init__(position=self._position, **kwargs)

        if self.image is not None:
            if not os.path.exists("./dynpy") and "./dynpy/" in self.image:

                type(self)._settle_dynpy()
                self.image = self.image.replace("./dynpy", "./._dynpy_env/dynpy")

            self.add_image(NoEscape(self.image), width=self.width)

        if self.caption is not None:
            self.add_caption(NoEscape(self.caption))

        if self.marker is not None:
            self.append(Label(self.marker))
        else:
            auto_mrk = AutoMarker(self).marker
            marker = auto_mrk

            self.marker = marker

            self.append(Label(marker))

    @property
    def preview_size(self):
        if self._preview_size is not None:

            return self._preview_size
        else:
            return self.__class__._preview_default_size

    @preview_size.setter
    def preview_size(self, size):
        if isinstance(size, float) or isinstance(size, int):
            size = NoEscape(f"{size}cm")
        if isinstance(size, str):
            size = NoEscape(size)

        self._preview_size = size

    def _get_str_key(self) -> str:
        """
        Generate a string key for the image based on its attributes.

        Returns:
            str: A string representation of the image key.
        """

        return self.image + "_caption:" + str(self.caption)

    def __repr__(self) -> str:
        """
        Represent the Picture object as a string.

        Returns:
            str: A string representation of the Picture object.
        """

        if self.image is not None:
            path = self.image
        else:
            path = "nothing to preview :("

        repr_string = f"""
        \n++++++++++ IMAGE +++++++
        \n++                    ++
        \n++  path:{path}       ++
        \n++                    ++
        \n++++++++++ IMAGE +++++++

        """

        return repr_string

    def _repr_markdown_(self) -> str:
        """
        Generate a Markdown representation of the Picture object.

        Returns:
            str: A Markdown representation for preview.
        """

        self.cls_container.append(self)

        if self.caption is None:
            caption = ""
        else:
            caption = self.caption

        if self.image is not None:
            path = self.image
            

            if "pdf" in path:
                # from wand.image import Image as WImage
                

                # img = WImage(filename=path, resolution=144)

                # hsize, vsize = img.size

                # img.resize(hsize, vsize)
                # display(img)
                # display(f"Fig. X: {caption}")

                import fitz
                
                from IPython.display import Image, display

                size = self.preview_size
                doc = fitz.open(path)
                page = doc.load_page(0)
                pixmap = page.get_pixmap(dpi=144)
                #pixmap.save(path.replace('pdf','png'))

                display(Image(data= pixmap.tobytes(output='png'),width=size))

                display(f"Fig. X: {caption}")
                return ""
            else:
                

                from IPython.display import Image,display
                size = self.preview_size


                display(Image(filename=path,
                              width=size,
                              ))
                display(f"Fig. X: {caption}")

                return ""
        else:
            return f"Nothing to plot \n \n Fig. X: {caption}"


    # def _repr_html_(self) -> str:
    #     """
    #     Generate a Markdown representation of the Picture object.

    #     Returns:
    #         str: A Markdown representation for preview.
    #     """

    #     self.cls_container.append(self)

    #     if self.caption is None:
    #         caption = ""
    #     else:
    #         caption = self.caption

    #     if self.image is not None:
    #         path = self.image

    #         if "pdf" in path:
    #             from wand.image import Image as WImage

    #             img = WImage(filename=path, resolution=144)

    #             hsize, vsize = img.size

    #             img.resize(hsize, vsize)
    #             display(img)
    #             display(f"Fig. X: {caption}")
    #             return ""
    #         else:
    #             size = self.preview_size

    #             from IPython.display import Image
    #             size = self.preview_size


    #             display(Image(filename=path,
    #                           #width=size,
    #                           ))
    #             display(f"Fig. X: {caption}")

    #             return ""
    #             #return f'<img src="{path}" alt="{caption}" style="width: 700px;" width="{size}"/>    \n \n Fig. X: {caption}'
    #             # return f'![image preview]({path}) \n \n Fig. X: {caption}'
    #     else:
    #         return f"Nothing to plot \n \n Fig. X: {caption}"


    def reported(self):
        self.cls_container.append(self)

        return copy.deepcopy(self)

    def call_code(self):

        # Declaration of command vector used in call_code method
        self.command_vector = []

        if self.image != None:
            self.command_vector.append('image = "{}"'.format(self.image))
        if self._position != None:
            self.command_vector.append('position = "{}"'.format(self._position))
        if self.caption != None:
            self.command_vector.append('caption = "{}"'.format(self.caption))
        if self.width != None:
            self.command_vector.append('width = NoEscape(r"{}")'.format(self.width))
        if self.height != None:
            self.command_vector.append("height = {}".format(self.height))
        if self.marker != None:
            self.command_vector.append('marker = "{}"'.format(self.marker))

        vector = "Picture("

        for element in self.command_vector:
            vector += str(element) + ", "

        print(vector[:-2] + ")")


class PltPlot(Picture, ReportModule):
    """
    A class that represents a figure environment in a LaTeX document.

    This class is used to insert and customize an image into the document, with options for caption, size, and position.

    Attributes:
        separate_paragraph (bool): Whether the float is positioned in a separate paragraph.
        _default_width (str): Default width of the image.
        packages (List[Union[Package, Command]]): Required LaTeX packages and commands for the Picture class.
        _position (str): Default position for the floating environment.
        _preview_default_size (str): Default size for previewing the image.

    Exemplary Usage:
        >>> Picture('./assets/exemplary-image.png', caption='Exemplary caption.')

    Costomised exemplary usage:
        >>> Picture(
        >>>     image='./assets/exemplary-image.png',
        >>>     caption='Exemplary caption.',
        >>>     position='h',
        >>>     width='0.5\\textwidth',
        >>>     height='10cm',
        >>>     marker='fig:example'
        >>> )
    """


    _latex_name = "figure"
    _default_path = "./tikzplots"
    _floats_no_gen = plots_no()
    _default_dpi = 300

    def __init__(
        self,
        caption=None,
        position="H",
        width=None,
        height=None,
        image=None,
        marker=None,
        **kwargs,
    ):
        """
        Initialize the Picture class.

        Args:
            image (Optional[str]): Picture or path to the image file.
            position (Optional[str]): LaTeX position specifier for the figure.
            caption (Optional[str]): Caption for the image.
            width (Optional[str]): Width of the image (e.g., '0.5\\textwidth').
            height (Optional[str]): Height of the image.
            marker (Optional[str]): Marker for referencing the figure in LaTeX.

        References:
            https://www.sharelatex.com/learn/Positioning_of_Figures
        """

        if isinstance(image, Picture):
            self.image = image.image
        elif isinstance(image, str):
            self.image = image
        else:
            self.image = self._save_plt_plot()

        self.caption = caption
        self.preview_size = None

        if width is not None:
            self.width = width
        else:
            self.width = self._default_width

        self.height = height

        self.marker = marker

        if position is not None:
            self._position = position

        super(Figure,self).__init__(position=self._position, **kwargs)

        if self.image is not None:
            if not os.path.exists("./dynpy") and "./dynpy/" in self.image:

                type(self)._settle_dynpy()
                self.image = self.image.replace("./dynpy", "./._dynpy_env/dynpy")

            self.add_image(NoEscape(self.image), width=self.width)

        if self.caption is not None:
            self.add_caption(NoEscape(self.caption))

        if self.marker is not None:
            self.append(Label(self.marker))
        else:
            auto_mrk = AutoMarker(self).marker
            marker = auto_mrk

            self.marker = marker

            self.append(Label(marker))


    
    def _save_plt_plot(self):
            
        import matplotlib.pyplot as plt

        fig_no = next(self._floats_no_gen)

        pic_filename = f"{self._default_path}/{self.__class__.__name__}{fig_no}.png"

        plt.savefig(pic_filename, bbox_inches='tight',dpi=self._default_dpi)
        plt.close()

        #self.append(Command('includegraphics', arguments=[pic_filename]))

        return pic_filename


class StarredPicture(Picture):
    _latex_name = "figure*"


class TikZFigure(Picture):
    """
    A class that represents a figure environment in a TikZPoster document.

    This class is used to insert and customize an image into the document, with options for caption, size, and position.

    Attributes:
        separate_paragraph (bool): Whether the float is positioned in a separate paragraph.
        _default_width (str): Default width of the image.
        packages (List[Union[Package, Command]]): Required LaTeX packages and commands for the Picture class.
        _position (str): Default position for the floating environment.
        _preview_default_size (str): Default size for previewing the image.

    Exemplary Usage:
        >>> TikZFigure('./assets/exemplary-image.png', caption='Exemplary caption.')

    Costomised exemplary usage:
        >>> TikZFigure(
        >>>     image='./assets/exemplary-image.png',
        >>>     caption='Exemplary caption.',
        >>>     position='h',
        >>>     width='0.5\\textwidth',
        >>>     height='10cm',
        >>>     marker='fig:example'
        >>> )
    """

    _latex_name = "tikzfigure"

    def __init__(
        self,
        image=None,
        position=None,
        caption=None,
        width=None,
        height=None,
        marker=None,
        **kwargs,
    ):
        """
        Initialize the Picture class.

        Args:
            image (Optional[str]): Path to the image file.
            position (Optional[str]): LaTeX position specifier for the figure.
            caption (Optional[str]): Caption for the image.
            width (Optional[str]): Width of the image (e.g., '0.5\\textwidth').
            height (Optional[str]): Height of the image.
            marker (Optional[str]): Marker for referencing the figure in LaTeX.

        References:
            https://www.sharelatex.com/learn/Positioning_of_Figures
        """

        super().__init__(
            image=image,
            position=caption,
            caption=None,
            width=width,
            height=height,
            marker=marker,
            **kwargs,
        )


class ObjectCode(LstListing, ReportModule):
    _latex_name = "lstlisting"
    default_title = "Listing"

    _default_header = None
    _style = "dynPyStyle"

    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self,
        source=None,
        caption=None,
        style=None,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
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

        self.source = source

        if style is not None:
            self._style = style
        self.arguments = arguments
        self.start_arguments = start_arguments
        self._caption = caption

        if self._caption is not None:
            caption_str = f"{{{ self._caption }}}"
            options_c = [NoEscape(f"caption={caption_str}")]
        else:
            options_c = []

        if self._style is not None:
            options_s = [NoEscape(f"style={self._style}")]
        else:
            options_s = []

        options_added = options_c + options_s

        if len(options_added) != 0:
            if options is not None:
                options = options + options_added
            else:
                options = options_added

        self.options = options

        super().__init__(
            options=options, arguments=None, start_arguments=None, **kwargs
        )

        self.code_type

    def __repr__(self):

        # self.reported()

        code = self.code_type

        repr_string = f"\n CODE\n {code} \n CODE"

        return repr_string

    @property
    def code_type(self):

        if self.source is None:
            code = "nothing to print :("

        elif type(self.source) == str:
            code = self.source
            # self.append(NoEscape('   \n   '+code+'   \n   '))

        else:

            code = inspect.getsource(self.source)
            # self.append(NoEscape('   \n   '+code+'   \n   '))

        if self._default_header is not None:

            return fix_code(self._default_header + "\n" + code)
        else:
            return fix_code(code)

    def reported(self):

        code = self.code_type

        self.append(NoEscape("   \n   " + code + "   \n   "))

        lst_env = Aspect(self.default_title)

        # lst_env.append(self)
        lst_env = self
        self.cls_container.append(lst_env)

        return copy.copy(self)

    def _repr_markdown_(self):

        self.reported()
        # print(self.code_type)

        return "```python \n" + self.code_type + "\n``` "

    def _as_string(self):

        return self.code_type

    def call_code(self):
        # Declaration of command vector used in call_code method
        self.command_vector = []

        if self.source != None:
            self.command_vector.append("source = {}".format(self.source))
        if self.options != None:
            self.command_vector.append("options = {}".format(self.options))
        if self.arguments != None:
            self.command_vector.append(
                "start_arguments = {}".format(self.start_arguments)
            )
        if self.start_arguments != None:
            self.command_vector.append("width = {}".format(self.width))
        if self._caption != None:
            self.command_vector.append("caption = {}".format(self._caption))

        vector = "ObjectCode("

        for element in self.command_vector:
            vector += str(element) + ", "

        print(vector[:-2] + ")")


guide_code_header = """
#!!! REMEMBER TO PASTE AND RUN ALL PREVIOUS CELL
#!!! REMEMBER TO PASTE AND RUN ALL PREVIOUS CELL
"""


class GuideCode(ObjectCode):

    _default_header = guide_code_header


class AiInterfaceOutput(ObjectCode):

    packages = [
        # Command('lstset', arguments = [NoEscape('style=nonumbers,language=Python')]),
        Command(
            NoEscape(
                r"lstdefinestyle{nonumbers}{numbers=none, breakatwhitespace=false, breaklines=true,captionpos=b, keepspaces=true, numbersep=5pt, showspaces=false,showstringspaces=false,showtabs=false,tabsize=2, basicstyle=\ttfamily\small}"
            )
        ),
        # Command('lstset', arguments = [NoEscape('style=nonumbers')]),
    ]

    _style = "nonumbers"


class TikZPicture(Environment, ReportModule):

    _latex_name = "tikzpicture"

    packages = [
        Package("tikz"),
        Package("polski", options=["MeX"]),
        Package("fontenc", options=["T1"]),
        #         Package('babel',options=['polish']),
        Package("amssymb"),
        Package("rotating"),
        NoEscape(
            r"\usetikzlibrary{snakes,calc,positioning,patterns,angles,quotes,decorations.pathmorphing,decorations.markings,through}"
        ),
        NoEscape(r"\usetikzlibrary{arrows,decorations,decorations.text}"),
        NoEscape(r"\usetikzlibrary{patterns.meta,math}"),
        NoEscape(r"\usetikzlibrary{shapes.geometric}"),
        NoEscape(r"\usepackage{circuitikz}"),
        #         Package('fontenc')
        # Package('hyperref'),
    ]

    _height = 5
    _width = 0.9
    _ylim = None

    _subplots_gap = 0.0
    _subplots_horizontal_gap = 0.5
    _default_grid = None
    _grid = None

    _in_figure = False
    _figure_gen = lambda: Figure(position="H")
    _image_parameters = {"width": None}

    _floats_no_gen = plots_no()
    _default_path = "./tikzplots"
    _filename = None
    _prefix = None

    _picture = True
    _caption = "Default caption"

    r"""A base class for LaTeX environments.
    This class implements the basics of a LaTeX environment. A LaTeX
    environment looks like this:
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False

    def __init__(
        self,
        tikzcode=None,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
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

        self._tikzcode = tikzcode
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        self.append(self.tikzcode)

    @property
    def tikzcode(self):

        if self._tikzcode is None:
            code = self._scheme_desc()
        else:
            code = self._tikzcode

        return NoEscape(code)

    def _scheme_desc(self):

        code = """
        Lack of code xD
        """

        return code

    def _repr_markdown_(self):

        return "\t" + self.reported().tikzcode.replace("\n", "\n \t")

    def reported(self):

        TikZPicture.cls_container.append(self)

        return copy.copy(self)

    def in_figure(
        self,
        caption=None,
        
        position=None,
        
        width=None,
        height=None,
        marker=None,
        filename=None,
        **kwargs,
    ):

        from .documents import tikz

        ReportCache.update_existing_files(self.__class__._default_path)

        obj = copy.copy(self)
        obj._in_figure = True

        if caption is not None:
            obj._caption = caption


        filename = self.filename

        standalone_plot = tikz.TikzStandalone(default_filepath=filename)
        standalone_plot.append(self)

        

        fig = self.__class__._figure_gen()
        fig.packages.append(Package("float"))
        #         display(fig.packages)
        #         fig.packages.remove(Package('fontenc',options=['T1']))
        #         fig.packages.remove(Package('inputenc',options=['utf8']))
        #         fig.packages.remove(Package('lmodern'))
        #         fig.packages.append(Command(NoEscape(r'''\if\ifxetex T\else\ifluatex T\else F\fi\fi T%
        #                                                   \usepackage{fontspec}
        #                                                 \else
        #                                                   \usepackage[T1]{fontenc}
        #                                                   \usepackage[utf8]{inputenc}
        #                                                   \usepackage{lmodern}
        #                                                 \fi''')))

        img_params = self.__class__._image_parameters

        if width is None:
            width = self.__class__._image_parameters["width"]

        if self._picture:
            from .report import Picture

            key = standalone_plot.dumps()

            if key in ReportCache._file_names:
                filename = ReportCache._file_names[key]

            else:

                from .creators import PdfLatexGenerator

                # standalone_plot.generate_pdf(
                #     filename, clean_tex=False
                # )  # ,compiler_args=['--lualatex']) #ZMIANA
                
                
                PdfLatexGenerator(standalone_plot).generate_file()
                
                ReportCache._file_names[key] = filename

            fig = Picture(
                filename + ".pdf",
                caption=caption,
                width=width,
                height=height,
                marker=marker,
                **kwargs,
            )
        else:
            standalone_plot.generate_tex(filename)

            fig.append(Command(command="input", arguments=filename))

        return fig

    @property
    def filename(self):

        if self._filename is None and self._prefix is None:
            filename = f"{self.__class__._default_path}/plot{self.__class__.__name__}{next(self.__class__._floats_no_gen)}"

        elif self._filename is None and self._prefix is not None:
            filename = f"{self.__class__._default_path}/{self._prefix}_plot{self.__class__.__name__}{next(self.__class__._floats_no_gen)}"

        else:
            filename = self._filename

        if filename in list(ReportCache._file_names.values()):
            filename = self.filename

        return filename


class SympyFormula(ReportModule):
    r"""
    This class appends a sympy expression to the existing document container.

    Arguments
    =========
    expr: str
        String that will be appended to the defined container.
    key_dict: dict
        Dictionary containing entries for string format method.
    marker: pylatex.labelref.Marker object
        User-defined marker for labeling and referencing math expressions.
    backend: func
        Type of backend used for processing provided Sympy expressions to LaTeX formula.
    **kwargs

    Methods
    =======

    Example
    =======
    """

    _color = None
    _break_mode = "autobreak"

    @classmethod
    def set_text_color(cls, color=None):

        cls._color = color

        return cls

    def __init__(self, expr, key_dict={}, marker=None, backend=vlatex, **kwargs):

        self._text = "Figures {first_marker}-{last_marker}"
        self._backend = backend
        self.latex_backend = backend

        self._marker = marker

        if not expr == None:
            self._expr = expr

        super().__init__(**kwargs)

        if self._break_mode == "autobreak":
            if isinstance(expr, (Matrix, ImmutableMatrix)):

                self._eq = Equation()
                self._eq.append(NoEscape(self._backend(self._expr)))

            elif isinstance(
                expr,
                (
                    Eq,
                    Relational,
                    StrictGreaterThan,
                    LessThan,
                    StrictLessThan,
                    GreaterThan,
                ),
            ):

                if isinstance(expr.lhs, (Matrix, ImmutableMatrix)) or isinstance(
                    expr.rhs, (Matrix, ImmutableMatrix)
                ):

                    self._eq = Equation()
                    self._eq.append(NoEscape(self._backend(self._expr)))

                else:
                    self._eq = Align()
                    with self._eq.create(AutoBreak()) as eq:
                        eq.latex_backend = self.latex_backend
                        eq.append_formula(expr)

            else:

                self._eq = Align()
                with self._eq.create(AutoBreak()) as eq:
                    eq.latex_backend = self.latex_backend
                    eq.append_formula(expr)

            if self._marker:
                marker = self._marker
            else:
                auto_mrk = AutoMarker(self._expr).marker
                marker = auto_mrk

            self._eq.append(Label(marker))
            self._marker = marker

        elif self._break_mode == "eq":
            self._eq = Equation()
            self._eq.append(NoEscape(self._backend(self._expr)))

            if self._marker:
                marker = self._marker
            else:
                auto_mrk = AutoMarker(self._expr).marker
                marker = auto_mrk

            self._eq.append(Label(marker))

    @property
    def latex_backend(self):
        return self._latex_backend

    @latex_backend.setter
    def latex_backend(self, backend):
        self._latex_backend = backend

    def __call__(self, analysis):

        display(self._expr)

        analysis._container.append(self._eq)

        return self._text

    def __str__(self):

        return self._backend(self._expr)

    #     def __repr__(self):

    #         display(self._expr)

    #         return ''

    def _repr_markdown_(self):
        # return f'\\begin{{equation}}\n  {self.reported()._backend(self._expr)} \n \\end{{equation}}'
        return f"$${self.reported()._backend(self._expr)}$$"

    def reported(self):

        if self.__class__._color:

            self._container.append(TextColor(self.__class__._color, self._eq))
        else:
            self._container.append(self._eq)

        return copy.copy(self)
        # return self

    def _latex(self):

        return self._backend(self._expr)


# class PlotTestResult:
#     r'''
#     The class creates a plot from data provided by DataStorage and appends it to an existing document container instance.

#     Arguments
#     =========
#     *args

#     keys_map: dict
#         Dictionary containing keys for mapping.
#     **kwargs


#     Methods
#     =======

#     Example
#     =======
#     '''
#     def __init__(self, *args, keys_map=None, **kwargs):

#         data_to_plot = DataStorage._storage

#         self._keys_map = {key: key for key in data_to_plot.keys()}

#         if keys_map:
#             self._keys_map = keys_map

#         self._data_to_plot = data_to_plot

#     def __call__(self, analysis, *args, **kwargs):
#         step_val = analysis._current_value

#         new_key = self._keys_map[step_val]

#         if not 'first_mrk' in DataStorage._dict.keys():
#             DataStorage._dict['first_mrk'] = Marker('first_mrk', 'fig')

#         DataStorage._dict['last_mrk'] = Marker('last_mrk', 'fig')

#         self._data_to_plot[new_key].plot()

#         ndp = DataPlot('wykres_nowy', position='H', preview=False)
#         ndp.add_data_plot(filename=f'Wykres_alpha_{next(plots_no_gen)}.png',
#                           width='11cm')
#         ndp.add_caption(
#             NoEscape(
#                 f'''Summary plot: simulation results for parameter - pomiary'''
#             ))
#         plt.show()

#         analysis._container.append(ndp)

#         return self._data_to_plot[new_key]


# class SystemDynamicsAnalyzer(ReportModule):
#     r'''
#     This is a computational block that runs a simulations on provided dynamic system. The class has methods responsible for data preparation, system analysis and reporting.

#     Arguments
#     =========
#     dynamic_system: Lagrange's method object
#         Dynamic model prepared basing on Sympy's Lagrange's method object.
#     reference_data: dict
#         Dictionary containing default values of systems's parameters.
#     report_init: list
#         List containing objects called at an initial part of report.
#     report_step: list
#         List containing objects called after the report_init.
#     report_end: list
#         List containing objects called after the report_step.

#     Methods
#     =======

#     Example
#     =======
#     '''
#     def __init__(self,
#                  dynamic_system,
#                  reference_data={},
#                  report_init=[ReportEntry('Report Beginning')],
#                  report_step=[
#                      SimulationalBlock(np.linspace(0, 300, 1000)).do_simulation
#                  ],
#                  report_end=[ReportEntry('Report End')],
#                  **kwargs):

#         self._dynamic_system = dynamic_system
#         self._reference_data = reference_data

#         self._init_steps = report_init
#         self._loop_steps = report_step
#         self._end_steps = report_end

#         self._fig_no = plots_no()

#         self._container = []
#         super().__init__()

#         #print(self._last_result)

#     def prepare_data(self, parameter, parameter_range=None):

#         self._parameter = parameter
#         self._parameter_range = parameter_range

#         #print('prepare data')

#         if isinstance(self._parameter, dict):
#             analysis_span_list = []
#             for key, value in parameter.items():

#                 if isinstance(value, list):

#                     for num, val in enumerate(value):

#                         analysis_span = {**self._reference_data, **{key: val}}
#                         analysis_span_list.append(analysis_span)
#                         #print(analysis_span_list)

#                 else:
#                     raise TypeError('Each dictionary value should be a list.')
#             self._analysis_span = analysis_span_list
#             self.value = value
#         else:
#             analysis_span = [{
#                 **self._reference_data,
#                 **{
#                     self._parameter: param_value
#                 }
#             } for param_value in parameter_range]
#             ##print(analysis_span)
#             self._analysis_span = analysis_span

#         return analysis_span

#     def analyze_system(self, t_span, container=None):

#         if container:
#             self._container = container
#         if self._dynamic_system:
#             solution_list = []

#             self.init_report()

#             for num, case_data in enumerate(self._analysis_span):
#                 self.num = num
#                 data_for_plot = self.analysis_step(case_data=case_data,
#                                                    t_span=t_span,
#                                                    ics_list=None)

#                 solution_list += [(case_data, data_for_plot)]

#             self.report_end()

#             self.solution_list = solution_list
#             return solution_list
#         else:
#             self.init_report()
#             ##print(self._analysis_span)

#             return (self._analysis_span)

#     def analysis_step(self, case_data, t_span, ics_list=None):

#         self._current_value = case_data[self._parameter]
#         self._current_data = case_data
#         ##print(self._current_data)
#         for action in self._loop_steps:
#             self._current_result = action(analysis=self)

#         self.report_step(self._current_result)

#         return self._current_result

#     def init_report(self, result_to_report=None):

#         for action in self._init_steps:
#             self._current_result = action(analysis=self)

#         return self._current_result

#     def report_step(self, result_to_report, container_type=None):

#         return self._current_result

#     def report_end(self, result_to_report=None, container_type=None):

#         for action in self._end_steps:
#             self._current_result = action(analysis=self)

#         return self._current_result


class InlineMath(Math):
    """A class representing a inline math environment."""

    def __init__(self, formula, escape=False, backend=vlatex):
        r"""
        Args
        ----
        data: list
            Content of the math container.
        inline: bool
            If the math should be displayed inline or not.
        escape : bool
            if True, will escape strings
        """

        self.escape = escape
        self.formula = vlatex(formula)
        self.backend = backend

        super().__init__(inline=True, data=backend(formula), escape=escape)


class SymbolsList(NoEscape):
    def __new__(cls, symbols_list, backend=vlatex):
        r"""
        Args
        ----
        symbols_list: list
            List of the Symbol objects to convert and append.
        backend: function
            Callable which is used to convert Symbol from list to its latex representation.
        escape : bool
            if True, will escape strings
        """

        list_str = f", ".join([f"\\( {backend(sym)} \\)" for sym in symbols_list])

        # return super(SymbolsList,cls).__new__(cls,list_str)
        return list_str


class NumbersList(NoEscape):
    def __new__(cls, numbers_list, backend=vlatex):
        r"""
        Args
        ----
        symbols_list: list
            List of the Symbol objects to convert and append.
        backend: function
            Callable which is used to convert Symbol from list to its latex representation.
        escape : bool
            if True, will escape strings
        """

        list_str = f", ".join([f"\\( {sym} \\)" for sym in numbers_list])

        return list_str


class DescriptionsRegistry:
    """
    A class representing LaTeX description registry of Symbols.

    classmethods:
    - set_descriprtions() - argument is the dictionary containing the descriptions of the symbols
    - set_description_mode() - argument is string containing the name of the chosen mode
    available modes:
        - single or one_occurence - each symbol makes only one occurence in the document (default)
        - all - for each equations all of the symbols are shown - and repeated

    Example of usage:
    from sympy import *
    from dynpy.utilities.report import SympyFormula, DescriptionsRegistry, SymbolsDescription

    a=Symbol("a", positive=True)
    b=Symbol("b", positive=True)
    c=Symbol("c", positive=True)

    eq=Eq(a**2+b**2,c**2)
    symbol_desc={
        a: 'Length of the side a',
        b: 'Length of the side b',
        c: 'Length of the side c'
    }

    DescriptionsRegistry.set_descriptions(symbol_desc)
    """

    _descriptions = {}
    _described_elements = {}
    _description_mode = "one occurrence"

    @classmethod
    def set_descriptions(cls, description={}):
        cls._descriptions = description
        return cls

    @classmethod
    def set_description_mode(cls, description_mode="single"):
        cls._description_mode = description_mode
        return cls

    @classmethod
    def reset_registry(cls):
        cls._described_elements = {}
        return cls

    def __init__(self, description_dict=None, method="add"):
        if description_dict is not None:
            self._descriptions = description_dict

    def _get_description(self, items, description_mode=None):
        symbols_list = self._descriptions.keys()

        missing_symbols_desc = {sym: "???" for sym in items if sym not in symbols_list}

        self._descriptions = {**self._descriptions, **missing_symbols_desc}
        self.__class__._descriptions = {
            **self.__class__._descriptions,
            **missing_symbols_desc,
        }

        # syms_to_desc={sym for sym in items if sym not in self._described_elements.keys()}
        if (
            self._description_mode == "one occurrence"
            or self._description_mode == "single"
        ):
            syms_to_desc = {
                sym: self._descriptions[sym]
                for sym in items
                if sym not in self._described_elements.keys()
            }
        elif self._description_mode == "all":
            syms_to_desc = {sym: self._descriptions[sym] for sym in items}

        else:
            syms_to_desc = {sym: self._descriptions[sym] for sym in items}

        self._described_elements = {**self._described_elements, **syms_to_desc}
        self.__class__._described_elements = {
            **self.__class__._described_elements,
            **syms_to_desc,
        }

        return syms_to_desc


class SymbolsDescription(Description, ReportModule):
    """A class representing LaTeX description environment of Symbols explained in description_dict.
    Example of usage:
    from sympy import *
    from dynpy.utilities.report import SympyFormula, DescriptionsRegistry, SymbolsDescription

    a=Symbol("a", positive=True)
    b=Symbol("b", positive=True)
    c=Symbol("c", positive=True)

    eq=Eq(a**2+b**2,c**2)
    symbol_desc={
        a: 'Length of the side a',
        b: 'Length of the side b',
        c: 'Length of the side c'
    }

    DescriptionsRegistry.set_descriptions(symbol_desc)
    display(SympyFormula(eq))
    display(SymbolsDescription(expr=eq))
    """

    _latex_name = "description"
    cls_container = []
    _description_head = "where:"

    @classmethod
    def set_container(cls, container=[]):
        cls.cls_container = container
        return cls

    @classmethod
    def set_default_header(cls, header="where:"):
        cls._description_head = header
        return cls

    def __init__(
        self,
        description_dict=None,
        expr=None,
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
        self.description_dict = description_dict
        self.expr = expr
        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        self._added_symbols = self._symbols_to_add_dict()

        self.add_items(self._added_symbols)

    def _symbols_to_add_dict(self):

        description_dict = self.description_dict
        expr = self.expr

        if description_dict is not None and expr is not None:

            symbols_set = expr.atoms(Symbol, Function, Derivative) - expr.atoms(
                sin, cos, tan, exp, Heaviside
            )

            symbols_to_add = {
                sym: desc
                for sym, desc in description_dict.items()
                if sym in symbols_set
            }

            return symbols_to_add

        elif (description_dict is None) and (expr is not None):

            symbols_set = set()
            if isinstance(expr, Iterable):
                for elem in expr:
                    symbols_set |= elem.atoms(
                        Symbol, Function, Derivative
                    ) - elem.atoms(sin, cos, tan, exp, atan, sign, Heaviside)
            # print(symbols_set)
            else:
                symbols_set |= expr.atoms(Symbol, Function, Derivative) - expr.atoms(
                    sin, cos, tan, exp, atan, sign, Heaviside
                )

            description_dict = DescriptionsRegistry()._get_description(symbols_set)

            # print(description_dict)

            symbols_to_add = {sym: desc for sym, desc in description_dict.items()}

            ##print('symbols')
            # print(symbols_to_add)

            return symbols_to_add

        elif description_dict is not None and expr is None:

            symbols_to_add = {sym: desc for sym, desc in description_dict.items()}

            return symbols_to_add

        else:
            return {}

    def reported(self, container=None):
        if container:
            self._container = container
        else:
            self._container = type(self).cls_container

        #         self._container.append(self)
        entries = [
            f"${vlatex(key)}$ - {value}" for key, value in self._added_symbols.items()
        ]

        # display(entries)
        # print(entries)

        end_sign = "."
        if len(entries) == 0:
            end_sign = ""

        if len(entries) != 0:

            self._container.append(NoEscape(self._description_head) + "\n")
            self._container.append(self)
        # return (self._text)

        # return copy.copy(self)

    def add_items(self, description_dict):

        end_symbol = "."

        if len(description_dict.keys()) > 0:
            last_key = list(description_dict.keys())[-1]
            end_symbol = ";"

        for label, entry in description_dict.items():

            if label == last_key:
                end_symbol = "."

            #             self.add_item(NoEscape(InlineMath(vlatex(label)).dumps()),
            #                           NoEscape(f'- {vlatex(entry)}{end_symbol}'))
            self.add_item(
                NoEscape(InlineMath((label)).dumps()),
                NoEscape(f"- {(entry)}{end_symbol}"),
            )

    def _repr_markdown_(self):
        self._container = type(self).cls_container

        entries = [
            f"${vlatex(key)}$ - {value}" for key, value in self._added_symbols.items()
        ]

        end_sign = "."
        head = self._description_head + "  \n"

        if len(entries) == 0:
            end_sign = ""
            head = ""

        text = head + ",  \n".join(entries) + end_sign

        if len(entries) != 0:

            self._container.append(NoEscape(self._description_head + "\n"))
            self._container.append(copy.copy(self))

        return text

    def __repr__(self):

        # self.reported()

        entries = [
            f"${vlatex(key)}$ - {value}" for key, value in self._added_symbols.items()
        ]

        end_sign = "."
        head = self._description_head + "  \n"

        if len(entries) == 0:
            end_sign = ""
            head = ""

        text = head + ",  \n".join(entries) + end_sign

        self._container = type(self).cls_container
        #         if len(entries) != 0:

        #            self._container.append(NoEscape(self._description_head)+ 'try' + '\n')
        #             self._container.append(copy.deepcopy(self))

        return text


class Align(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package("mathtools")]
    escape = False
    content_separator = "\n"


class Equation(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package("mathtools")]
    escape = False
    content_separator = "\n"


class CompositeMatrix(Equation):
    r"""A class representing a composite matrix surrounded by Equation environment
    .. code-block:: latex
        \begin{environment_name}
            Some content that is in the environment
        \end{environment_name}
    The text that is used in the place of environment_name is by default the
    name of the class in lowercase.
    However, this default can be overridden in 2 ways:
    1. setting the _latex_name class variable when declaring the class
    2. setting the _latex_name attribute when initialising object
    """

    #: Set to true if this full container should be equivalent to an empty
    #: string if it has no content.
    omit_if_empty = False
    _latex_name = "equation"

    def __init__(
        self,
        matrix=None,
        repr_sym="A",
        options=None,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
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

        self.matrix = matrix
        self.repr_sym = "A"
        self.options = options
        self.arguments = arguments
        self.start_arguments = start_arguments

        super().__init__(options=None, arguments=None, start_arguments=None, **kwargs)

        if self.matrix is not None:

            display(self._sym_matrix)
            self.append(NoEscape(latex(self._sym_matrix)))

    @property
    def _sym_matrix(self):
        subs_list = [
            (elem, Symbol(f"{self.repr_sym}_{no}"))
            for no, elem in enumerate(self.matrix)
            if elem != 0
        ]

        return self.matrix.subs(subs_list)

    def __repr__(self):

        if self.inspected_obj is not None:
            code = inspect.getsource(self.inspected_obj)
        else:
            code = "nothing to print :("

        repr_string = f"""
        \n++++++++++ CODE +++++++
        \n{code}
        \n++++++++++ CODE +++++++

        """

        return repr_string

    def reported(self):
        self.cls_container.append(self)

        return copy.deepcopy(self)


class DMath(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package("flexisym"), Package("breqn")]
    escape = False
    content_separator = "\n"


class AutoBreak(Environment):
    """A class to wrap LaTeX's alltt environment."""

    packages = [Package("mathtools"), Package("autobreak")]
    escape = False
    content_separator = " "
    _latex_backend = vlatex

    @property
    def latex_backend(self):
        return self._latex_backend

    @latex_backend.setter
    def latex_backend(self, backend):
        self._latex_backend = backend

    def _split_expr(self, expr):

        if isinstance(expr, Equality):

            elems = [expr.lhs, Symbol("="), expr.rhs]

        elif isinstance(expr, StrictGreaterThan):

            elems = [expr.lhs, Symbol(">"), expr.rhs]

        elif isinstance(expr, LessThan):

            elems = [expr.lhs, Symbol(r"\leq"), expr.rhs]

        elif isinstance(expr, GreaterThan):

            elems = [expr.lhs, Symbol(r"\geq"), expr.rhs]

        elif isinstance(expr, StrictLessThan):

            elems = [expr.lhs, Symbol("<"), expr.rhs]

        elif isinstance(expr, Add):

            elems = list(expr.args)
            # elems = sum(([obj,Symbol('+')] for obj in elems),[]  )[0:-1]
            # elems=[]

            # for obj in list(expr.args):
            #    display(obj.args)
            #    elems += [obj,Symbol('+')]

        else:
            elems = [expr]

        # print('elems')
        # display(elems)

        if len(elems) > 1:
            # print('elems')
            # display(elems)
            elems = sum([self._split_expr(obj) for obj in elems], [])
        else:
            elems = [expr]

        return elems

    def append_formula(self, expr):

        terms = self._split_expr(expr)

        new_terms = []

        for no, obj in enumerate(terms):

            if (
                terms[no - 1] == Symbol("=")
                or terms[no - 1] == Symbol(">")
                or terms[no - 1] == Symbol(r"\leq")
                or terms[no - 1] == Symbol(r"\geq")
                or terms[no - 1] == Symbol("<")
            ):
                new_terms += [Symbol("\n"), obj, Symbol("\n")]

            elif isinstance(obj, Mul) and (
                (any([elem.is_negative for elem in obj.args]))
            ):
                # print(' any negative')
                # display(obj)
                new_terms += [Symbol("\n -"), -obj]

            elif isinstance(obj, Mul) and not (
                (any([elem.is_negative for elem in obj.args]))
            ):
                # print(' not any negative')
                # display(obj)
                new_terms += [Symbol("\n +"), obj]

            elif (
                obj == Symbol("=")
                or obj == Symbol(">")
                or obj == Symbol("\leq")
                or obj == Symbol("<")
                or obj == Symbol("\geq")
            ):
                new_terms += [obj]

            elif isinstance(obj, (Symbol, Function, Number, Derivative, Subs, Expr)):
                new_terms += [Symbol("\n +"), obj]

            else:
                new_terms += [obj]

        prev_term = ""

        if new_terms[0] == Symbol("\n -"):
            new_terms_loop = new_terms
        else:
            self.append("\n")
            new_terms_loop = new_terms[1:]

        for term in new_terms_loop:
            #             print('+++',self.__class__.latex_backend(term))
            #             print('++ prev +++',prev_term)
            #             self.append(self.__class__.latex_backend(term)) ######### TA LINIA ZAKOMENTOWANA /AMADI
            self.append(self.latex_backend(term))
            # prev_term=self.__class__.latex_backend(term)
        self.append("\n")


############# TEST


# class EqRef(Environment):
#         """A class to wrap LaTeX's alltt environment."""

#     packages = [Package('breqn'), Package('flexisym')]
#     escape = False
#     content_separator = " "


# class ReportModelAnalysis:

#     chapter_name = 'Model analityczny układ o swobody'

#     def __init__(self, model=None):
#         geometry_options = {
#             "margin": "1cm",
#         }
#         doc = Document(documentclass='report',
#                        geometry_options=geometry_options)
#         #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

#         doc.packages.append(Package('fontenc'))
#         #doc.packages.append(Package('fontspec'))
#         doc.packages.append(Package('standalone'))
#         doc.packages.append(Package('siunitx'))
#         doc.packages.append(Package('amsmath'))
#         doc.packages.append(Package('float'))
#         doc.packages.append(Package('tikz'))
#         doc.packages.append(Package('pgfplots'))
#         doc.packages.append(Package('preamble'))
#         doc.packages.append(Package('polski', options=["MeX"]))
#         doc.packages.append(Command('usepgfplotslibrary', arguments='units'))
#         #         doc.packages.append(Package('preamble'))
#         self.doc = doc

#         self.model = model

#         self.chapter_name = type(self).chapter_name

#     def lagranges(self):

#         doc = self.doc

#         model = self.model

#         #musisz też dopisać o jednym stopniu - to będzie ten indeks
#         with doc.create(Chapter(self.chapter_name)):

#             with doc.create(DMath()):
#                 doc.append('T=')

#                 doc.append(
#                     vlatex(model.lagrangian() +
#                            model.lagrangian().subs({vel: 0
#                                                     for vel in model.u})))
#             with doc.create(DMath()):
#                 doc.append('V=')

#                 doc.append(
#                     vlatex(-self.model.lagrangian().subs(
#                         {vel: 0
#                          for vel in model.u})))

#             with doc.create(DMath()):
#                 doc.append('L=')

#                 doc.append(vlatex(self.model.lagrangian()))
#         return doc

#     def stat_and_dyn(self):
#         doc = self.doc

#         model = self.model

#         with doc.create(Chapter(self.chapter_name)):

#             for i in range(len(model.q)):
#                 with doc.create(DMath()):
#                     doc.append(
#                         vlatex(self.model.equilibrium_equation().doit()[i]))
#             for lhs, rhs in zip(
#                     symbols('k_1:{ndof}_1:{ndof}'.format(ndof=len(model.q) +
#                                                          1)),
#                 (list(model.linearized().stiffness_matrix()))):
#                 with doc.create(DMath()):
#                     doc.append(vlatex([Eq(lhs, rhs)]))
#             for i in range(len(model.q)):
#                 with doc.create(DMath()):
#                     doc.append(vlatex(self.model._eoms[i]))
#             for lhs, rhs in zip(
#                     symbols('m_1:{ndof}_1:{ndof}'.format(ndof=len(model.q) +
#                                                          1)),
#                 (list(model.inertia_matrix()))):
#                 with doc.create(DMath()):
#                     doc.append(vlatex([Eq(lhs, rhs)]))
#         return doc

#     def doit(self):

#         doc = self.doc

#         model = self.model

#         with doc.create(Chapter(self.chapter_name)):
#             doc.append(
#                 '''Model przyjęty do weryfikacji z danymi doświadaczalnymi stanowił złożenie modelu o trzech stopniach oraz układu napędowego o jednym stopniu swobody. Funkcja Lagrangea tego układu miała następującą postać:'''
#             )

#             with doc.create(DMath()):
#                 doc.append('L=')

#                 doc.append(vlatex(model.lagrangian()))

#             doc.append(
#                 '''Bazując na wyznaczonej funkcji Lagrange'a określono równania równowagi (specjalny przypadek równań Lagrange'a drugiego rodzaju")'''
#             )

#             for no, balance_eq in enumerate(model.equilibrium_equation()):
#                 with doc.create(DMath()):

#                     #doc.append(vlatex(chair_dict['linearized_ramp'].q[no])+ ':~~~~'  +  vlatex(Eq(balance_eq.doit(),0) )   )
#                     doc.append(vlatex(Eq(balance_eq.doit(), 0)))

#             doc.append(
#                 '''Bazując na wyznaczonej funkcji Lagrange'a określono równania równowagi (specjalny przypadek równań Lagrange'a drugiego rodzaju)'''
#             )

#             for no, eom in enumerate(model._eoms):
#                 with doc.create(DMath()):

#                     #doc.append(vlatex(chair_dict['linearized_ramp'].q[no])+ ':~~~~'  +  vlatex(Eq(balance_eq.doit(),0) )   )
#                     doc.append(vlatex(Eq(eom.doit(), 0)))

#         return doc


# class RandomDescription:
#     def __init__(self, *args, random_or_order=None, **kwargs):

#         self.pool_list = args
#         self.random_or_order = random_or_order

#     def to_string(self):

#         if hasattr(self.random_or_order, '__contains__'):
#             self.order = list(random_or_order)
#         else:
#             self.selected_sentences = [
#                 random.choice(pool) for pool in self.pool_list
#             ]

#         return ' '.join(self.selected_sentences)

#     def __str__(self):

#         return self.to_string()

#     def __repr__(self):
#         return self.__class__.__name__ + ': ' + self.to_string()


################################################## IDK WHAT'S IT


# class DataStorage:
#     r'''
#     This class represents data collector that stores the results able to obtain within other computational blocks.
#     It ensures easy data transferring in order to do efficient analysis and data processing.

#     Arguments
#     =========
#         data_set: dict
#             Dictionary of data to be deposited into storage.
#     Methods
#     =======
#         reset_storage(cls,*args,**kwargs):
#             Classmethod; Cleans storage.
#     Example
#     =======
#         >>>elements=['first','second','third','fourth','fifth']
#         >>>new_dict={'abs':312}

#         >>>DS1=DataStorage(new_dict)

#         >>>DataStorage._storage=new_dict

#         >>>DataStorage._dict['x']='123'

#         >>>DataStorage._list=['a','b','c']

#         >>>DataStorage._plot_markers_dict={elem:Marker('plot1'+elem ,'fig')   for elem in elements}

#         >>>DataStorage._subplot_markers_dict={elem:Marker('subplot1'+elem ,'fig')   for elem in elements}

#         >>>DataStorage.reset_storage()
#     '''

#     _storage = {}
#     _dict = {}
#     _plot_markers_dict = {}
#     _subplot_markers_dict = {}
#     _list = []
#     first_marker = None
#     last_marker = None

#     #_frame = TimeDataFrame()
#     last_result = None
#     _last_result = None

#     @property
#     def frame(self):

#         if self._frame:
#             time_frame = TimeDataFrame(self._frame)

#             time_frame.columns = pd.MultiIndex.from_tuples(
#                 list(time_frame.columns))
#             time_frame.columns.set_names(['model', 'parameter', 'coordinate'],
#                                          inplace=True)

#         else:
#             time_frame = TimeDataFrame(self.__class__._frame)

#             time_frame.columns = pd.MultiIndex.from_tuples(
#                 list(self.__class__._frame.columns))
#             time_frame.columns.set_names(['model', 'parameter', 'coordinate'],
#                                          inplace=True)

#         return time_frame

#     def set_frame(self, key, frame):
#         self._frame[key] = frame

#         return frame

#     def set_class_frame(self, key, frame):

#         self.__class__._frame[key] = frame

#         return frame

#     def __init__(self, data_set=None):

#         self._data_set = data_set
#         self._marker_dict = {}
#         self._marker_dict_sub = {}

#         self._storage = self._data_set
#         self._story_point = []
#         self._frame = data_set
#         self._last_result = None
#         self._list

#         _frame = TimeDataFrame()
#         last_result = None
#         _last_result = None

#     @classmethod
#     def reset_storage(cls, *args, **kwargs):

#         cls._storage = {}
#         cls._dict = {}
#         cls._list = []

#         return cls


# class NumericalDataSet:
#     def __init__(self,
#                  numerical_data,
#                  key=None,
#                  *args,
#                  keys_map=None,
#                  label='Experimental data',
#                  **kwargs):

#         data_to_plot = numerical_data

#         self._key = key

#         self._keys_map = {key: key for key in data_to_plot.keys()}

#         if keys_map:
#             self._keys_map = keys_map

#         self._data_to_plot = data_to_plot
#         self._label = label

#         self.dvars = list(list(data_to_plot.values())[0].columns)

#     def __call__(self, analysis, *args, **kwargs):
#         step_val = analysis._current_value

#     def __str__(self):
#         return self._label

#     def numerized(self, params_values={}, **kwargs):

#         if params_values == {}:
#             return copy.copy(self)
#         else:

#             return copy.copy(self)

#     def compute_solution(self,
#                          t_span=None,
#                          ic_list=None,
#                          t_eval=None,
#                          params_values=None,
#                          method='RK45'):

#         #         #print('compute_solution')
#         #         display(params_values)
#         #         display(ic_list)

#         if ic_list:
#             #print('class ics has been taken')
#             self.ics = ic_list

#         display(self._data_to_plot)

#         return self._data_to_plot[params_values[self._key]]


# class SimulationalBlock(ReportModule):
#     r'''
#     It is computational module which enables to perform numerical simulations of the dynamic system.
#     Class provides several methods devoted for data processing, ploting and reporting.

#     Arguments
#     =========
#     t_span: iterable
#         Time span.
#     ics_list: iterable
#         List containing values of initial conditions.
#     dynamic_system: Lagrange's method object
#     reference_data: dict
#         Dictionary containing default values of systems's parameters.
#     **kwargs

#     Methods
#     =======
#     frame(self):
#         Property; Sets time frame.
#     reset_storage(cls):
#         Classmethod; cleans storage including TimeDataFrame.
#     set_t_span(cls,t_span):
#         Classmethod; sets time span.
#     label_formatter(self,analysis=None,label_generator=None):
#         Provides label generated basing on current simulation or defined dictionary.
#     show_eoms(self,analysis,**kwargs):

#     do_simulation(self,analysis,**kwargs):

#     simulation_result(self,analysis):

#     plot_result(cls,analysis):


#     Example
#     =======
#     '''

#     #_frame.columns=pd.MultiIndex()

#     _list = []
#     _dict = {}
#     general_t_span = []

#     _model = None
#     _reference_data = None
#     _hold = False
#     last_result = []

#     _ref_data = None

#     def holded(self, hold=True):
#         self.__class__._hold = hold

#         return self.copy()

#     def label_formatter(self, analysis=None, label_generator=None):

#         if analysis:
#             var = analysis._parameter
#             value = analysis._current_value

#             if self._dynamic_system:
#                 system = self._dynamic_system
#             else:
#                 system = analysis._dynamic_system

#             #print('system label', str(system))
#             label = Eq(var, value, evaluate=False), str(system)

#         else:

#             label = ((var, value) for var, value in self._ref_data.items())

#         return label

#     @classmethod
#     def reset_storage(cls):
#         cls._frame = TimeDataFrame()

#         cls._list = []
#         cls._dict = {}
#         return cls

#     @classmethod
#     def set_t_span(cls, t_span):

#         cls.general_t_span = t_span

#         return cls

#     def __init__(self,
#                  t_span=None,
#                  ics_list=None,
#                  dynamic_system=None,
#                  reference_data=None,
#                  label=None,
#                  **kwargs):

#         self._t_span = t_span
#         self._ics_list = ics_list

#         if t_span is not None:
#             self._t_span = t_span
#         else:
#             self._t_span = type(self).general_t_span

#         self._numerical_system = None

#         self._dynamic_system = dynamic_system
#         self._ref_data = reference_data
#         #print(self.__class__, label)
#         self._model_label = label

#         super().__init__()

#     def show_eoms(self, analysis, **kwargs):

#         if self._ref_data:
#             case_data = self._ref_data
#             var = analysis._parameter
#             value = analysis._current_value

#             case_data[var] = value

#         else:
#             case_data = analysis._current_data

#         if self._dynamic_system:
#             dynamic_system = self._dynamic_system
#         else:
#             dynamic_system = analysis._dynamic_system

#         display(dynamic_system._eoms)
#         return dynamic_system

#     def do_simulation(self, analysis, **kwargs):

#         if self._ref_data:
#             case_data = self._ref_data
#             var = analysis._parameter
#             value = analysis._current_value

#             case_data[var] = value

#         else:
#             case_data = analysis._current_data

#         print('$' * 100)
#         display(case_data)

#         if self._dynamic_system:
#             dynamic_system = self._dynamic_system
#         else:
#             dynamic_system = analysis._dynamic_system

#         if self._model_label:
#             dynamic_system._label = self._model_label

#         #print('dynamic model name', dynamic_system._label)

#         if not self._numerical_system:

#             #             display(analysis._dynamic_system.system_parameters())
#             #             display(analysis._dynamic_system._eoms)

#             self._numerical_system = dynamic_system.numerized(
#                 parameter_values=case_data)

#         numerical_system = self._numerical_system
#         no_dof = len((numerical_system.dvars))

#         if not self._ics_list:
#             ics_list = [0] * no_dof
#         else:
#             ics_list = self._ics_list

#         #print('numer', numerical_system)

#         #print('ics_self', self._ics_list)
#         #print('ics', ics_list)

#         simulation_result = numerical_system.compute_solution(
#             t_span=self._t_span,
#             ic_list=ics_list,
#             t_eval=self._t_span,
#             params_values=case_data)

#         self._simulation_result = simulation_result

#         var = analysis._parameter
#         value = analysis._current_value

#         label = self.label_formatter(analysis)

#         DataStorage._storage[label] = simulation_result

#         self.__class__._frame[[
#             label + (coord, ) for coord in simulation_result.columns
#         ]] = simulation_result
#         self.__class__._list += [simulation_result]
#         self.__class__._last_result = type(self.__class__._frame)()
#         self.__class__._last_result[[
#             (coord, ) + label for coord in simulation_result.columns
#         ]] = simulation_result

#         self._frame[[label + (coord, ) for coord in simulation_result.columns
#                      ]] = simulation_result
#         self._list += [simulation_result]
#         self._last_result = type(self._frame)()
#         self._last_result[[
#             (coord, ) + label for coord in simulation_result.columns
#         ]] = simulation_result

#         ##print(DataStorage._list)

#         DataStorage._list += [simulation_result]

#         if (analysis) is not None:
#             analysis._frame[[
#                 label + (coord, ) for coord in simulation_result.columns
#             ]] = simulation_result
#             analysis._list += [simulation_result]
#             analysis._last_result = type(analysis._frame)()
#             analysis._last_result[[
#                 (coord, ) + label for coord in simulation_result.columns
#             ]] = simulation_result

#         return simulation_result

#     def simulation_result(self, analysis):
#         return self._simulation_result

#     @classmethod
#     def plot_result(cls, analysis):

#         last_result = DataStorage._list[-1]

#         plot_of_result = last_result.plot()

#         ndp = DataPlot('wykres_nowy', position='H', preview=False)

#         ndp.add_data_plot(
#             filename=
#             f'{cls.cls_path}/block_simulation_{next(plots_no_gen)}.png',
#             width='11cm')

#         ndp.add_caption(
#             NoEscape(
#                fr'''Summary plot: simulation results for parameter \({latex(analysis._parameter)}\)'''

#             ))
#         #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#         #ndp.append(Label(self.marker_dict[coord]))

#         analysis._container.append(ndp)

#         plt.show()

#         return plot_of_result


# class Summary(ReportModule):
#     def __init__(self,
#                  block=None,
#                  coordinate=None,
#                  caption=None,
#                  label=None,
#                  subplots=False,
#                  height=None,
#                  extra_commands=None):

#         if subplots:
#             self._subplot = subplots
#         else:
#             self._subplot = self.__class__._subplot

#         if coordinate:
#             #print('hello')
#             self._coord = coordinate

#         else:
#             self._coord = slice(None, None)

#         #print(f'the coord is {self._coord}')

#         super().__init__()

#         #         #print(self._frame)
#         #         #print(block,'abc')
#         #         #print(type(block))

#         self._block = block
#         if block:

#             self._frame = block._frame

#             self._last_result = block._last_result
#         if caption:
#             self._caption = caption
#         else:
#             self._caption = self.__class__._caption

#         if label:
#             self._label = label
#         else:
#             self._label = self.__class__._label

#         if height:
#             self._height = height
#         else:
#             self._height = self.__class__._height

#         if extra_commands is not None:
#             self._extra_commands = extra_commands
#         else:
#             self._extra_commands = None

#     def holded(self, hold=True):
#         self.__class__._hold = hold

#         return copy.copy(self)

#     def plot(self, *args, analysis=None, **kwargs):

#         #         #print('analysis')
#         #         #print(analysis)
#         if analysis:
#             #print('data is pushed to store - analysis')
#             self.set_frame(analysis._last_result.columns,
#                            analysis._last_result)
#             self.set_class_frame(analysis._last_result.columns,
#                                  analysis._last_result)

#             self.__class__._frame.index.name = analysis._last_result.index.name
#             self._frame.index.name = analysis._last_result.index.name

#             print('_' * 100, analysis._last_result._get_comp_time())

#         #print('summary plot - call')
#         #print((self._block), type((self._block)))

#         result_of_plot = None

#         if (self._block)._last_result is not None:
#             ##print()

#             result_to_add = type(self._block)._last_result
#             print('_' * 100, result_to_add._get_comp_time())
#             columns_to_add = result_to_add.columns

#             #print('plot index', result_to_add.index.name)

#             if self._frame is not None:
#                 #print('data is pushed to store - self -block')
#                 self.set_frame(columns_to_add, result_to_add)
#                 self._frame.index.name = result_to_add.index.name
#                 display(self._frame)

#                 self.set_class_frame(columns_to_add, result_to_add)
#                 self.__class__._frame.index.name = result_to_add.index.name
#                 display(self.__class__._frame)

#                 result_of_plot = self._frame.plot(*args, **kwargs)

#             else:

#                 #print('data is pushed to store -  class -block')

#                 self.set_class_frame(columns_to_add, result_to_add)
#                 result_of_plot = self.__class__._frame.plot(*args, **kwargs)
#                 self.__class__._frame.index.name = result_to_add.index.name

#         plt.clf()

#         return result_of_plot

#     def show(
#         self,
#         analysis=None,
#         legend_pos='north east',
#         legend_columns=1,
#         colors_list=[
#             'blue', 'red', 'green', 'orange', 'violet', 'magenta', 'cyan'
#         ],
#         extra_commands=None,
#         **kwargs,
#     ):

#         #print('show')
#         ##print(type(self)._frame)

#         result = type(self)()

#         if (type(self)._frame) is not None:

#             data = (type(self)._frame)
#             if not data.columns.empty:
#                 data.columns = pd.MultiIndex.from_tuples(data.columns)

#                 result = data[self._coord]

#         if ((self._frame)) is not None:

#             data = ((self)._frame)
#             if not data.columns.empty:
#                 data.columns = pd.MultiIndex.from_tuples(data.columns)

#                 result = data[self._coord]

#         #print('story in Sumary', (type(self)._frame).empty)
#         #print('out format', self._out_format)
#         #         #print(result)

#         formatter = self._out_format()

#         if self._coord != slice(None, None):
#             ylabel = f'$ {formatter._match_unit(self._coord)}$'
#         else:
#             ylabel = 'coords'

#         if not (self)._frame.empty:  #result.empty:

#             data = ((self)._frame)
#             data.columns = pd.MultiIndex.from_tuples(data.columns)

#             self._apply_formatter(data[self._coord]).plot(
#                 ylabel=ylabel, subplots=self._subplot)

#             #print('o tu - from self')
#             ##print(pd.MultiIndex.from_tuples(list(result.columns)))
#             plot_of_result = plt.show()
#             self.clear_frame()

#         elif not type(self)._frame.empty:  #result.empty:

#             data = ((self).__class__._frame)
#             data.columns = pd.MultiIndex.from_tuples(data.columns)

#             self._apply_formatter(data[self._coord]).plot(
#                 ylabel=ylabel, subplots=self._subplot)

#             #print('o tu - from cls')
#             ##print(pd.MultiIndex.from_tuples(list(result.columns)))
#             plot_of_result = plt.show()
#             self.clear_frame(obj=False)

#         if self._autoreporting:

#             filepath = f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#             units = self.__class__._units
#             #print('units', units)
#             if self._coord in units:
#                 y_unit_str = f'y unit = {units[self._coord]:~Lx}'.replace(
#                     '[]', '')

#             else:
#                 y_unit_str = ''

#             new_data = self._apply_formatter(data[self._coord])

#             if new_data.__class__._domain:
#                 ivar = new_data.__class__._domain
#             else:
#                 ivar = data[self._coord].index.name

#             if ivar in units:
#                 x_unit_str = f'x unit = {units[ivar]:~Lx}'.replace('[]', '')

#             else:
#                 x_unit_str = ''

#             #print('y_unit_str', y_unit_str)

#             if extra_commands is None:
#                 extra_commands = self._extra_commands

#             fig = new_data.to_standalone_figure(
#                 filepath,
#                 colors_list=colors_list,
#                 subplots=self._subplot,
#                 height=self._height,
#                 width=NoEscape(r'0.9\textwidth'),
#                 x_axis_description=
#                 f',xlabel=${NoEscape(vlatex(ivar))}$, {x_unit_str},'.replace(
#                     '$$', '$'),
#                 y_axis_description=
#                 f'ylabel=${NoEscape(vlatex(self._coord))}$, {y_unit_str},',
#                 legend_pos=legend_pos,
#                 extra_commands=extra_commands,
#             )
#             fig.add_caption(NoEscape(self._caption))
#             fig.append(Label(self._label))

#             self._container.append(fig)

#         return result

#     def summary(self, analysis=None, **kwargs):

#         #print('summary')

#         result = type(self)()

#         if (type(self)._frame) is not None:

#             data = (type(self)._frame)
#             if not data.columns.empty:
#                 data.columns = pd.MultiIndex.from_tuples(data.columns)

#                 result = data[self._coord]

#         if ((self._frame)) is not None:

#             data = ((self)._frame)
#             if not data.columns.empty:
#                 data.columns = pd.MultiIndex.from_tuples(data.columns)

#                 result = data[self._coord]

#         #print('story in Sumary', (type(self)._frame).empty)

#         #         #print(result)
#         if not (self)._frame.empty:  #result.empty:

#             data = ((self)._frame)
#             data.columns = pd.MultiIndex.from_tuples(data.columns)

#             display(self._apply_formatter(data[self._coord]))

#             #print('o tu - from self')
#             ##print(pd.MultiIndex.from_tuples(list(result.columns)))

#             self.clear_frame()

#         elif not type(self)._frame.empty:  #result.empty:

#             data = ((self).__class__._frame)
#             data.columns = pd.MultiIndex.from_tuples(data.columns)

#             display(self._apply_formatter(data[self._coord]))

#             #print('o tu - from cls')
#             ##print(pd.MultiIndex.from_tuples(list(result.columns)))

#             self.clear_frame(obj=False)


# #         #print(self._apply_formatter(data[self._coord]).to_latex(escape=False).replace('\\toprule','\\toprule \n \\midrule').replace('\\bottomrule','\\midrule \n \\bottomrule')   )
# #         #print(self._apply_formatter(data[self._coord]).to_latex(escape=False).replace('\\toprule','\\toprule \n \\midrule').__repr__() )
# #         display(self._apply_formatter(data[self._coord]).to_latex(escape=False).replace('\\toprule','\\toprule \n \\midrule') )

#         if self._autoreporting:
#             tab = Table(position='!htb')

#             #                         f'xyz_{next(plots_no_gen)}',
#             #                         subplots=self.__class__._subplot,

#             #                         height=NoEscape(r'6cm'),
#             #                         width=NoEscape(r'0.9\textwidth'),
#             #                         y_axis_description='',
#             #                         legend_pos='north west')))

#             tab.add_caption(NoEscape(self._caption))
#             tab.append(
#                 NoEscape(
#                     self._apply_formatter(data[self._coord]).style.to_latex(
#                         escape=False,
#                         caption='').replace('\\toprule',
#                                             '\\toprule \n \\midrule').replace(
#                                                 '\\bottomrule',
#                                                 '\\midrule \n \\bottomrule')))
#             tab.append(Label(self._label))

#             self._container.append(tab)

#         return self._apply_formatter(data[self._coord])

#     def prepare_markers(self,
#                         analysis=None,
#                         coordinate=None,
#                         xlim=None,
#                         **kwargs):

#         DataStorage._plot_markers_dict = {
#             elem: Marker(f'plot{self.__class__.__name__}', 'fig')
#             for elem in self._frame
#         }
#         DataStorage._subplot_markers_dict = {
#             elem: Marker(f'subplot{self.__class__.__name__}', 'fig')
#             for elem in self._frame
#         }
#         #         DataStorage.first_marker = list(
#         #             DataStorage._plot_markers_dict.values())[0]
#         #         DataStorage.last_marker = list(
#         #             DataStorage._plot_markers_dict.values())[-1]
#         #         self.last_marker = list(DataStorage._plot_markers_dict.values())[-1]
#         #         type(self)._last_marker = list(
#         #             DataStorage._plot_markers_dict.values())[-1]
#         #         #print('marker - def')
#         #         #print(self.last_marker)

#         return analysis


# class SimulationFFT:
#     r'''It is a class that provides Fast Fourier Transform techniques for formerly performed numerical simulations in time domain. Class supplies a method for plotting a double sided RMS.

#     Arguments
#     =========
#     *args

#     Methods
#     =======

#     Example
#     =======
#     '''
#     def __init__(self, *args):
#         self._args = args

#     @classmethod
#     def plot_fft(cls, analysis):

#         last_result = DataStorage._list[-1]

#         fft_result = last_result.to_frequency_domain().double_sided_rms()

#         fft_result.plot(xlim=(0, None), subplots=True, logy=False)
#         #         plt.yscale('log')
#         plt.show()
#         return fft_result

#     @classmethod
#     def fft(cls, analysis=None):

#         fft_of_storage = {
#             key: result.to_frequency_domain().double_sided_rms()
#             for key, result in DataStorage._storage.items()
#         }

#         DataStorage._storage = fft_of_storage

#         return fft_of_storage


# class AccelerationComparison(ReportModule):
#     r'''
#     It is computational block that prepares the comparison of particular coordinates regarding to changes of selected parameter.
#     Class provides several methods devoted for data processing, ploting and reporting.

#     Arguments
#     =========
#     t_span: iterable
#         Time span.
#     ics_list: iterable
#         List containing values of initial conditions.
#     data: dict
#         Dictionary consisting of data for comparison.
#     label: str
#         User-defined LaTeX label for generated plots.

#     Methods
#     =======

#     Example
#     =======
#     '''
#     _subplot = False
#     _height = r'7cm'

#     _story_point = None

#     general_t_span = None

#     _data_storage = {}

#     _last_marker = None

#     _formatter = lambda entry: f'${latex(entry[0].lhs)} = {round(entry[0].rhs/1000)} \\si {{\\tonne}} ({ (entry[0].rhs/10000000*100).n(2,chop=True)  } \\% m_v ) $'

#     @classmethod
#     def set_t_span(cls, t_span):

#         cls.general_t_span = t_span

#         return cls

#     @classmethod
#     def set_label_formatter(cls, formatter):

#         cls._formatter = formatter

#         return cls

#     @classmethod
#     def reset_storage(cls):

#         cls._story_point = {}
#         cls._data_storage = {}

#         return cls

#     def __init__(self, t_span=None, data=None, ics_list=None, label=None):

#         self.last_marker = None
#         self._t_span = t_span
#         self._ics_list = ics_list

#         if t_span is not None:
#             self._t_span = t_span
#         else:
#             self._t_span = type(self).general_t_span

#         if data:
#             self._data = data
#         else:
#             self._data = DataStorage._storage

#         if label:
#             self._label = label
#         else:
#             self._label = None

#         super().__init__(None)

#     def _prepare_data(self, coordinate=None, xlim=None):

#         if xlim:
#             data = {
#                 key: result.truncate(xlim[0], xlim[-1])
#                 for key, result in self._data.items()
#             }
#         else:
#             data = self._data


# #         #print('_______________test of plot_____________')
# #         #print(data)
# #         #print('_______________test of plot_____________')
# #         #print(data)

#         elements = list((data.values()))[0].columns
#         #         #print('frametype')
#         #         #print(type(list((data.values()))[0])())
#         summaries_dict = {
#             dynsym: type(list((data.values()))[0])()
#             for dynsym in elements
#         }

#         for key, result in data.items():
#             for coord in elements:
#                 summaries_dict[coord][key] = result[coord]
#         type(self)._story_point = summaries_dict

#         if coordinate:
#             return summaries_dict[coordinate]
#         else:
#             return summaries_dict

#     def prepare_summary(self, analysis=None, coordinate=None, xlim=None):

#         if analysis:
#             self._analysis = analysis

#         result = self._prepare_data(xlim=xlim)

#         elements = result.keys()

#         if self._label:
#             DataStorage._plot_markers_dict = {
#                 elem: Marker(f'plot{self.__class__.__name__}{self._label}',
#                              'fig')
#                 for elem in elements
#             }
#             DataStorage._subplot_markers_dict = {
#                 elem: Marker(f'subplot{self.__class__.__name__}{self._label}',
#                              'fig')
#                 for elem in elements
#             }
#         else:

#             DataStorage._plot_markers_dict = {
#                 elem:
#                 Marker(f'plot{self.__class__.__name__}{next(plots_no_gen)}',
#                        'fig')
#                 for elem in elements
#             }
#             DataStorage._subplot_markers_dict = {
#                 elem:
#                 Marker(f'subplot{self.__class__.__name__}{next(plots_no_gen)}',
#                        'fig')
#                 for elem in elements
#             }

#         DataStorage.first_marker = list(
#             DataStorage._plot_markers_dict.values())[0]
#         DataStorage.last_marker = list(
#             DataStorage._plot_markers_dict.values())[-1]
#         self.last_marker = list(DataStorage._plot_markers_dict.values())[-1]
#         type(self)._last_marker = list(
#             DataStorage._plot_markers_dict.values())[-1]
#         #print('marker - def')
#         #print(self.last_marker)

#         return result

#     def plot_summary(self,
#                      analysis=None,
#                      coordinate=None,
#                      xlim=None,
#                      subplots=_subplot,
#                      legend_pos='north east',
#                      legend_columns=1,
#                      colors_list=[
#                          'blue', 'red', 'green', 'orange', 'violet', 'magenta',
#                          'cyan'
#                      ],
#                      extra_commands=None,
#                      options=None):

#         self.subplots = subplots
#         if analysis:
#             self._analysis = analysis
#             self._parameter = analysis._parameter
#         else:
#             self._parameter = 'which name is missing.'

#         if coordinate:
#             if not isinstance(coordinate, list):
#                 coordinate = [coordinate]

#             data_dict = {
#                 coord: self._prepare_data(xlim=xlim)[coord]
#                 for coord in coordinate
#             }

#         else:
#             data_dict = self._prepare_data(xlim=xlim)

#         if self.__class__._subplot == False:
#             for coord, data in data_dict.items():

#                 data.plot()
#                 plt.ylabel(coord)

#                 filepath = f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                 ########### for tikz
#                 #ndp=DataPlot('wykres_nowy',position='H',preview=False)

#                 #it should be replaced with data.rename
#                 #             #print(data)

#                 data.columns = [
#                     type(self)._formatter(label) for label in data.columns
#                 ]
#                 #print(type(self)._units)
#                 y_unit_str = f'{(type(self)._units[coord]):Lx}'.replace(
#                     '[]', '')

#                 ndp = data.to_standalone_figure(
#                     filepath,
#                     colors_list=colors_list,
#                     height=NoEscape(r'7cm'),
#                     width=NoEscape(r'12cm'),
#                     y_axis_description=NoEscape(
#                         fr',ylabel=${vlatex(coord)}$,y unit={y_unit_str}, x unit=\si{{\second}}'

#                     ),
#                     legend_pos=legend_pos + ',' +
#                     f'legend columns= {legend_columns}',
#                     extra_commands=extra_commands,
#                     options=options)
#                 #ndp.add_data_plot(filename=f'{self._path}/{self.__class__.__name__}_data_{next(plots_no_gen)}.png',width='11cm')

#                 ########### for tikz
#                 #ndp.append(data.to_pylatex_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))

#                 ndp.add_caption(NoEscape(f'{type(self)._caption}'))

#                 plt.show()

#                 #print('marker - plot')
#                 #print(self.last_marker)
#                 #print(type(self)._last_marker)
#                 #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#                 ndp.append(Label(type(self)._last_marker))

#                 if analysis:
#                     analysis._container.append(ndp)
#                 else:
#                     filepath = f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                     #latex_code=(TimeDataFrame(data).to_tikz_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))
#                     self._container.append(ndp)

#         else:
#             for coord, data in data_dict.items():
#                 data.plot(subplots=self.__class__._subplot, ylabel=coord)
#                 filepath = f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                 ########### for tikz
#                 #ndp=DataPlot('wykres_nowy',position='H',preview=False)

#                 #it should be replaced with data.rename
#                 #             #print(data)

#                 data.columns = [
#                     type(self)._formatter(label) for label in data.columns
#                 ]
#                 #print(type(self)._units)
#                 y_unit_str = f'{(type(self)._units[coord]):Lx}'.replace(
#                     '[]', '')

#                 ndp = data.to_standalone_figure(
#                     filepath,
#                     subplots=self.__class__._subplot,
#                     colors_list=colors_list,
#                     height=NoEscape(r'6cm'),
#                     width=NoEscape(r'0.9\textwidth'),
#                     y_axis_description=NoEscape(
#                         fr',ylabel=${vlatex(coord)}$,y unit={y_unit_str} ,x unit=\si{{\second}}'
#                     ),
#                     legend_pos=legend_pos + ',' +
#                     f'legend columns= {legend_columns}',
#                     extra_commands=extra_commands,
#                     options=options)
#                 #ndp.add_data_plot(filename=f'{self._path}/{self.__class__.__name__}_data_{next(plots_no_gen)}.png',width='11cm')

#                 ########### for tikz
#                 #ndp.append(data.to_pylatex_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))

#                 ndp.add_caption(NoEscape(f'{type(self)._caption}'))

#                 plt.show()

#                 #print('marker - plot')
#                 #print(self.last_marker)
#                 #print(type(self)._last_marker)
#                 #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#                 ndp.append(Label(type(self)._last_marker))

#                 if analysis:
#                     analysis._container.append(ndp)
#                 else:
#                     filepath = f'{self._path}/{self.__class__.__name__}_tikz_{next(plots_no_gen)}'

#                     #latex_code=(TimeDataFrame(data).to_tikz_plot(filepath,colors_list=['blue','red','green','orange','violet','magenta','cyan'],height=NoEscape(r'5.5cm'),width=NoEscape(r'0.5\textwidth')))
#                     self._container.append(ndp)
#             return ndp

#     def plot_max_summary(self, analysis):

#         if not type(self)._story_point:
#             type(self)._story_point = self._prepare_data()

#         data_dict = type(self)._story_point

#         new_data_dict = {
#             key: data.abs().max()
#             for key, data in data_dict.items()
#         }

#         df = pd.DataFrame(new_data_dict)
#         df.plot()
#         plt.show()
#         df.plot(subplots=True)
#         plt.show()

#         ndp = DataPlot('wykres_nowy1', position='H', preview=False)
#         ndp.add_data_plot(
#             filename=
#             f'{self._path}/{self.__class__.__name__}_max_data_{next(plots_no_gen)}.png',
#             width='11cm')
#         ndp.add_caption(
#             NoEscape(
#                 fr'''Summary plot: simulation results for parameter \({latex(analysis._parameter)}\)'''
#             ))
#         #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#         #ndp.append(Label(self.marker_dict[coord]))

#         if analysis:
#             analysis._container.append(ndp)
#         else:
#             self._container.append(ndp)

#         return None

#     def plot_mean_summary(self, analysis):

#         if not type(self)._story_point:
#             type(self)._story_point = self._prepare_data()

#         data_dict = type(self)._story_point

#         new_data_dict = {
#             key: data.abs().mean()
#             for key, data in data_dict.items()
#         }

#         df = pd.DataFrame(new_data_dict)
#         df.plot()
#         plt.show()
#         df.plot(subplots=True)
#         plt.show()

#         ndp = DataPlot('wykres_nowy2', position='H', preview=False)

#         ndp.add_data_plot(filename=f'Wykres_mean_{next(plots_no_gen)}.png',
#                           width='11cm')
#         ndp.add_caption(
#             NoEscape(
#                 fr'''Summary plot: simulation results for parameter \({latex(analysis._parameter)}\)'''
#             ))
#         #ndp.add_caption(NoEscape(f'''Summary plot: simulation results for \({coord}\) coodinate and parameter \({latex(analysis._parameter)}\) values: {prams_vals_str} {units_dict[par]:~Lx}'''))
#         #ndp.append(Label(self.marker_dict[coord]))

#         if analysis:
#             analysis._container.append(ndp)
#         else:
#             self._container.append(ndp)

#         return None

#     def simulation_result(self, analysis):
#         return type(self)._story_point

#     @classmethod
#     def get_from_storage(cls, storage=DataStorage):

#         type(self)._story_point = storage._storage

#         return cls

#     @property
#     def data_storage(self):
#         return type(self)._data_storage

#     def plot_result(self, analysis, **kwargs):

#         self._simulation_result = type(self)._story_point

#         self._simulation_result.plot()
#         plt.show()

#         self._simulation_result.plot(subplots=True)
#         plt.show()

#         return self._simulation_result


# class FFTComparison(AccelerationComparison):
#     r'''
#     It is computational block that prepares the comparison of FFT signals of particular coordinates regarding to changes of selected parameter.
#     Class inherits from AccelerationComparison and provides several methods devoted for data processing, ploting and reporting.

#     Arguments
#     =========
#     t_span: iterable
#         Time span.
#     ics_list: iterable
#         List containing values of initial conditions.
#     data: dict
#         Dictionary consisting data for comparison.
#     label: str
#         User-defined LaTeX label for generated plots.

#     Methods
#     =======

#     Example
#     =======
#     '''
#     def _prepare_data(self, coordinate=None, xlim=None):

#         data = self._data

#         self._data = {
#             key: value.to_frequency_domain().double_sided_rms()
#             for key, value in data.items()
#         }

#         return super()._prepare_data(coordinate=None, xlim=xlim)


# class SummaryTable(ReportModule):
#     r'''
#     It is computational block that prepares the summary table of particular coordinates regarding to changes of selected parameter.
#     Class provides several methods devoted for data processing and reporting.

#     Arguments
#     =========
#     t_span: iterable
#         Time span.
#     ics_list: iterable
#         List containing values of initial conditions.
#     data: dict
#         Dictionary consisting data for comparison.
#     label: str
#         User-defined LaTeX label for generated plots.

#     Methods
#     =======

#     Example
#     =======
#     '''

#     _story_point = None

#     general_t_span = None

#     _data_storage = {}

#     _last_marker = None

#     _formatter = lambda entry: f'${latex(entry.lhs)} = {round(entry.rhs/1000)} \\si {{\\tonne}} ({ (entry.rhs/10000000*100).n(2,chop=True)  } \\% m_v ) $'

#     @classmethod
#     def set_t_span(cls, t_span):

#         cls.general_t_span = t_span

#         return cls

#     @classmethod
#     def set_label_foramatter(cls, formatter):

#         cls._formatter = formatter

#         return cls

#     @classmethod
#     def reset_storage(cls):

#         cls._story_point = {}
#         cls._data_storage = {}

#         return cls

#     def __init__(self,
#                  t_span=None,
#                  data=None,
#                  ics_list=None,
#                  label=None,
#                  **kwargs):

#         self.last_marker = None
#         self._t_span = t_span
#         self._ics_list = ics_list

#         if t_span is not None:
#             self._t_span = t_span
#         else:
#             self._t_span = type(self).general_t_span

#         if data:
#             self._data = data
#         else:
#             ##print(SimulationalBlock().frame)
#             self._data = SimulationalBlock().frame

#         if label:
#             self._label = label
#         else:
#             self._label = ''

#         super().__init__(None)

#     def _prepare_data(self, coordinate=None, level=None, xlim=None):

#         if xlim:
#             data = {
#                 key: result.truncate(xlim[0], xlim[-1])
#                 for key, result in self._data.items()
#             }
#         else:
#             data = self._data

# #         #print('_______________test of plot_____________')
# #         #print(data)
# #         #print('_______________test of plot_____________')
# #         #print(data)
#         elements = data.columns
#         #         #print('frametype')
#         #         #print(type(list((data.values()))[0])())
#         summaries_dict = data.swaplevel(0, -1, axis=1)

#         if coordinate:
#             return summaries_dict[coordinate]
#         else:
#             return summaries_dict

#     def prepare_summary(self,
#                         analysis=None,
#                         coordinate=None,
#                         xlim=None,
#                         **kwargs):

#         if analysis:
#             self._analysis = analysis

#         data_table = self._prepare_data(coordinate=coordinate,
#                                         xlim=xlim).abs().max().to_frame().T

#         models = data_table.columns.get_level_values(0).unique()
#         ##print(models)

#         #display([data_table[model].T for model in models])

#         result = type(data_table)()

#         for model in models:
#             result[model] = data_table[model].T

#         elements = result.keys()

#         DataStorage._plot_markers_dict = {
#             elem: Marker(f'plot{self.__class__.__name__}{self._label}', 'fig')
#             for elem in elements
#         }
#         DataStorage._subplot_markers_dict = {
#             elem: Marker(f'subplot{self.__class__.__name__}{self._label}',
#                          'fig')
#             for elem in elements
#         }
#         DataStorage.first_marker = list(
#             DataStorage._plot_markers_dict.values())[0]
#         DataStorage.last_marker = list(
#             DataStorage._plot_markers_dict.values())[-1]
#         self.last_marker = list(DataStorage._plot_markers_dict.values())[-1]
#         type(self)._last_marker = list(
#             DataStorage._plot_markers_dict.values())[-1]
#         #print('marker - def')
#         #print(self.last_marker)

#         return result

#     def show_summary(self,
#                      analysis=None,
#                      coordinate=None,
#                      xlim=None,
#                      legend_pos='north east',
#                      legend_columns=1,
#                      colors_list=[
#                          'blue', 'red', 'green', 'orange', 'violet', 'magenta',
#                          'cyan'
#                      ]):

#         #self.subplots=subplots
#         if analysis:
#             self._analysis = analysis
#             self._parameter = analysis._parameter
#         else:
#             self._parameter = 'which name is missing.'

#         data_table = self.prepare_summary(coordinate=coordinate, xlim=xlim)

#         display(data_table)

#         latex_table = NoEscape(data_table.style.to_latex())

#         if analysis:
#             analysis._container.append(NoEscape(latex_table))
#         else:

#             self._container.append(NoEscape(latex_table))

#         return data_table


# class BaseSeriesFormatter(pd.Series):
#     """
#         A base class for formatting data plots. This class extends TimeSeries and provides methods for
#         customizing the appearance and structure of data, including column names, units, labels, and legends.

#         Attributes:
#             _cols_name (str or None): Name of the columns to be used for formatting.
#             _domain (str or None): The domain for the data (if applicable).
#             _default_sep (str): The default separator for columns. Defaults to ', '.

#         Methods:
#             set_column_separator(sep): Sets the default column separator.
#             set_data_filtter(filtter): Sets a custom data filter function.
#             set_units_dict(units): Sets a dictionary of units for the data.
#             columns_name(name): Sets the column name(s) for the data.
#             set_multiindex(): Converts the index to a MultiIndex.
#             cut_coord(coord): Filters the data by a specific coordinate from a MultiIndex.
#             format_labels(): Formats the index labels, especially for MultiIndex.
#             set_ylabel(label): Sets the y-axis label for a plot.
#             set_xlabel(label): Sets the x-axis label for a plot.
#             set_legend(legend): Sets the legend for the plot.

#         Example:
#             >>> import pandas as pd
#             >>> from dynpy.utilities.report import BaseSeriesFormatter

#             >>> data = pd.Series([1, 2, 3], index=[("Temperature", "Day1"), ("Pressure", "Day2"), ("Humidity", "Day3")])
#             >>> formatter = BaseSeriesFormatter(data)

#             >>> formatter.set_column_separator(" | ")

#             >>> formatter.set_data_filtter(lambda df: df * 2)

#             >>> formatter.set_units_dict({"temperature": "C", "pressure": "Pa"})

#             >>> labels = formatter.format_labels()

#             >>> xlabel = formatter.set_xlabel()
#             >>> ylabel = formatter.set_ylabel()

#             >>> filtered_data = formatter.cut_coord(("Temperature", "Day1"))

#             >>> print(labels, xlabel, ylabel)
#     """


#     _cols_name = None
#     _domain = None
#     _default_sep = ', '

#     @classmethod
#     def set_column_separator(cls, sep: str = ', ') -> 'BaseSeriesFormatter':
#         """
#             Sets the default separator for columns.

#             Args:
#                 sep (str): The column separator to be used. Defaults to ', '.

#             Returns:
#                 BaseSeriesFormatter: The class with the updated column separator.
#         """

#         cls._default_sep = sep

#         return cls

#     @classmethod
#     def set_data_filtter(cls, filtter=lambda frame: frame.copy()) -> 'BaseSeriesFormatter':
#         """
#             Sets a custom data filter function to process the data.

#             Args:
#                 filtter (callable): A function that takes a DataFrame and returns a processed DataFrame.
#                                     Common examples include:
#                                     - Filtering out rows based on a condition.
#                                     - Selecting specific columns for processing.
#                                     - Scaling or normalizing data.
#                                     - Filling missing values.

#                                     Defaults to a function that returns a copy of the input frame.

#             Example:
#                 >>> import pandas as pd
#                 >>> from dynpy.utilities.report import BaseSeriesFormatter

#                 >>> data = pd.Series([1, 2, 3], index=[1, 2, 3])
#                 >>> formatter = BaseSeriesFormatter(data)

#                 >>> formatter.set_data_filtter(lambda s, *args: s * 2)

#                 >>> filtered_data = formatter._data_filtter(formatter)

#                 >>> print(filtered_data)

#             Returns:
#                 BaseSeriesFormatter: The class with the updated data filter.
#         """

#         cls._data_filtter = filtter

#         return cls

#     @classmethod
#     def set_units_dict(cls, units: dict = {}) -> 'BaseSeriesFormatter':
#         """
#             Sets a dictionary of units for the data.

#             Args:
#                 units (dict): A dictionary where keys are column names and values are unit strings.
#                               Defaults to an empty dictionary.

#             Returns:
#                 BaseSeriesFormatter: The class with the updated units dictionary.
#         """

#         cls._units = units

#         return cls

#     @classmethod
#     def columns_name(cls, name=None) -> 'BaseSeriesFormatter':
#         """
#             Sets the column name(s) for the data.

#             Args:
#                 name (str or list or None): The column name(s) to set. Defaults to None.

#             Returns:
#                 BaseSeriesFormatter: The class with the updated column name(s).
#         """

#         cls._cols_name = name

#         return cls

#     @property
#     def _constructor(self):
#         """
#             Returns the constructor for creating new BaseSeriesFormatter instances.

#             Returns:
#                 type: The BaseSeriesFormatter class.
#         """

#         return BaseSeriesFormatter

#     @property
#     def _constructor_expanddim(self):
#         """
#             Returns the constructor for creating a DataFrame-like object from the series.

#             Returns:
#                 type: The BaseFrameFormatter class (or equivalent).
#         """

#         return BaseFrameFormatter

#     @property
#     def _constructor_sliced(self):
#         """
#             Returns the constructor for creating a sliced series object.

#             Returns:
#                 type: The BaseSeriesFormatter class.
#         """

#         return BaseSeriesFormatter

#     def set_multiindex(self) -> 'BaseSeriesFormatter':
#         """
#             Converts the index to a MultiIndex.

#             Returns:
#                 BaseSeriesFormatter: A new instance with a MultiIndex applied to the index.

#             Example:
#                 >>> import pandas as pd
#                 >>> from dynpy.utilities.report import BaseSeriesFormatter

#                 >>> data = pd.Series([1, 2, 3], index=["A", "B", "C"])
#                 >>> formatter = BaseSeriesFormatter(data)

#                 >>> formatter.index = [("Group1", "A"), ("Group1", "B"), ("Group2", "C")]

#                 >>> multiindex_formatter = formatter.set_multiindex()

#                 >>> print(multiindex_formatter.index)
#         """

#         midx = pd.MultiIndex.from_tuples(self.index)
#         new_obj = self.__class__(self.copy())
#         new_obj.index = midx

#         return new_obj

#     def cut_coord(self, coord):
#         """
#             Filters the data by a specific coordinate from a MultiIndex.

#             Args:
#                 coord: The coordinate to filter by.

#             Returns:
#                 DataFrame or Series: A filtered view of the data corresponding to the specified coordinate.
#         """

#         new_frame = self.set_multiindex()

#         return new_frame[coord]

#     def format_labels(self):
#         """
#             Formats the index labels, especially for MultiIndex.

#             Returns:
#                 list: A list of formatted index labels if the index is a MultiIndex.
#         """

#         if isinstance(self.index, pd.MultiIndex):
#             return self.index.tolist()
#         else:
#             return self.index.tolist()

#     def set_ylabel(self, label: str = None) -> str:
#         """
#             Sets the y-axis label for a plot.

#             Args:
#                 label (str or None): The y-axis label to set. If None, generates labels from the index using LaTeX formatting.

#             Returns:
#                 str: The y-axis label. If no custom label is provided, it is derived from the first element of each index.

#             Example:
#                 >>> import pandas as pd

#                 >>> data = pd.Series([1, 2, 3], index=[("Temperature", "Zone1"), ("Pressure", "Zone2"), ("Humidity", "Zone3")])
#                 >>> formatter = BaseSeriesFormatter(data)

#                 >>> ylabel = formatter.set_ylabel()

#                 >>> print(ylabel)  # Outputs: Temperature
#         """

#         if label is None:
#             labels = self.format_labels()
#             if labels:
#                 label = f'$ {latex(labels[0])} )$'

#         return label

#     def set_xlabel(self, label: str = None) -> str:
#         """
#             Sets the x-axis label for a plot.

#             Args:
#                 label (str or None): The x-axis label to set. If None, the method attempts to generate a label
#                                      based on the index name or structure.

#             Returns:
#                 str: The x-axis label. If no custom label is provided, it defaults to the index name or a generic label.

#             Example:
#                 >>> import pandas as pd
#                 >>> from dynpy.utilities.report import BaseSeriesFormatter

#                 >>> data = pd.Series([1, 2, 3], index=pd.MultiIndex.from_tuples([("Time", "Start"), ("Time", "Middle"), ("Time", "End")], names=["Category", "Subcategory"],))
#                 >>> formatter = BaseSeriesFormatter(data)

#                 >>> xlabel = formatter.set_xlabel()

#                 >>> print(xlabel)  # Outputs: Category
#         """

#         if label is None:
#             if isinstance(self.index, pd.MultiIndex):
#                 label = self.index.names[0] if self.index.names else "Index"
#             else:
#                 label = self.index.name if self.index.name else "Index"

#         return label

#     def set_legend(self, legend: list = None) -> list:
#         """
#             Sets the legend for the plot.

#             Args:
#                 legend (list or None): A list of legend labels to use. If None, legends are automatically
#                                        generated based on column names or indices.

#             Returns:
#                 list: The final list of legend labels to be used.

#             Example:
#                 >>> import pandas as pd
#                 >>> from dynpy.utilities.report import BaseSeriesFormatter

#                 >>> # Create a Series with a name
#                 >>> data = pd.Series([1, 2, 3], name="Sample Series")
#                 >>> formatter = BaseSeriesFormatter(data)

#                 >>> # Case 1: Provide a custom legend
#                 >>> custom_legend = formatter.set_legend(["Custom Label"])
#                 >>> print(custom_legend)  # Outputs: ["Custom Label"]

#                 >>> # Case 2: Use automatic legend generation
#                 >>> auto_legend = formatter.set_legend()
#                 >>> print(auto_legend)  # Outputs: ["Sample Series"]

#                 >>> # Case 3: Mismatched legend length (not applicable here since Series expects only one label)
#                 >>> try:
#                 >>>     formatter.set_legend(["Too Many", "Labels"])
#                 >>> except ValueError as e:
#                 >>>     print(e)  # Outputs: "Legend for a Series must be a single label."

#         """

#         if legend is None:
#             if isinstance(self, pd.DataFrame):
#                 legend = self.columns.tolist()
#             elif isinstance(self, pd.Series):
#                 legend = [self.name] if self.name else ["Series"]

#         if isinstance(self, pd.DataFrame) and len(legend) != len(self.columns):
#             raise ValueError("Length of legend must match the number of columns in the data.")

#         if isinstance(self, pd.Series) and len(legend) != 1:
#             raise ValueError("Legend for a Series must be a single label.")

#         return legend


# class BaseFrameFormatter(pd.DataFrame):

#     _units = {}
#     _ylabel = None
#     _label_formatter = None
#     _data_filtter = lambda frame: frame.copy()
#     _cols_name = None
#     _domain = None
#     r'''
#     Basic class for formatting data plots. It provides methods for setting options

#     Arguments
#     =========

#     Methods
#     =======

#     Example
#     =======

#     '''

#     _default_sep = ', '

#     @classmethod
#     def set_column_separator(cls, sep=', '):
#         cls._default_sep = sep

#         return cls

#     @classmethod
#     def set_columns_name(cls, name=None):
#         cls._cols_name = name

#         return cls

#     @classmethod
#     def set_data_filtter(cls, filtter=lambda frame: frame.copy()):

#         cls._data_filtter = filtter
#         #print(cls._data_filtter)
#         return cls

#     @classmethod
#     def set_units_dict(cls, units={}):

#         cls._units = units
#         return cls

#     def _match_unit(self, sym, latex_printer=vlatex):

#         units = self.__class__._units

#         if isinstance(sym, Eq):
#             sym_check = sym.lhs
#         else:
#             sym_check = sym
#         #print(f'variable is {sym}')
#         #print(f'units are {units}')
#         if sym_check in units:
#             #print(f'{sym} is in units_dict')
#             #print('matched unit', f'{units[sym_check]:~L}')
#             return f'{latex_printer(sym)}~[{units[sym_check]:~L}]'

#             #return f'{latex(sym)}'
#         else:
#             return f'{latex_printer(sym)}'

#     def _format_label(self, obj):
#         if self.__class__._label_formatter:
#             return self.__class__._label_formatter(obj)
#         else:
#             if isinstance(obj, Iterable):
#                 obj = obj
#             else:
#                 obj = obj,

#             return self.__class__._default_sep.join(
#                 f'${self._match_unit(elem)}$' if isinstance(elem, (
#                     Expr, Eq)) else f'{elem}' for elem in obj if elem != 0)

#     def format_labels(self):
#         if isinstance(self.columns, pd.MultiIndex):

#             idx = self.columns.tolist()

#         else:
#             idx = self.columns
#         #print('idx', idx)
#         new_idx = []
#         for entry in idx:
#             #print('entry', entry, type(entry))
#             if entry in self.__class__._units and entry != 0:

#                 units = self.__class__._units

#                 #new_idx+=[ self._format_label(entry) + f'{units[entry]}']
#                 new_idx += [self._format_label(entry)]
#             elif entry != 0:
#                 new_idx += [self._format_label(entry)]

#         #print('new idx', new_idx)

#         self.__class__._ylabel = (list(idx)[0][0])

#         new_obj = self.copy()

#         new_obj.columns = new_idx

#         cols_name = self.__class__._cols_name

#         if cols_name:
#             idx_name = self._format_label(cols_name)
#             new_obj.columns.name = idx_name

#         return new_obj

#     @property
#     def _constructor(self):
#         return BaseFrameFormatter

#     @property
#     def _constructor_expanddim(self):
#         return BaseFrameFormatter

#     @property
#     def _constructor_sliced(self):
#         return BaseSeriesFormatter

#     def set_multiindex(self, names=None):

#         midx = pd.MultiIndex.from_tuples(self.columns, names=names)

#         new_obj = self.__class__(self)
#         new_obj.columns = midx

#         return new_obj

#     def cut_coord(self, coord):

#         new_frame = self.set_multiindex()

#         return new_frame[coord]

#     def format_index(self):
#         if isinstance(self.index, pd.MultiIndex):

#             idx = self.index.tolist()

#         else:
#             idx = self.index
#         #print('idx', idx)
#         new_idx = idx.copy()
#         #print('new idx', new_idx)

#         new_obj = self.copy()
#         new_obj.index = new_idx

#         #print('new_obj.index', new_obj.index.name)

#         if not new_obj.index.name:
#             new_obj.index.name = Symbol('f')

#         #print('new_obj.index', new_obj.index.name)

#         new_obj.index.name = f'${self._match_unit(new_obj.index.name)}$'
#         #new_obj.index.name = 'cos'

#         ##print('new_obj.index.name',new_obj.index.name)

#         return new_obj

#     def set_ylabel(self, label=None):
#         if isinstance(self.index, pd.MultiIndex):

#             if label == None:
#                 new_obj = self.copy()

#                 for label in new_obj.columns.tolist():
#                     if label in self.__class__._units:

#                         y_unit_str = f'[${type(self)._units[label]}$]'

#                     else:
#                         y_unit_str = ''

#                 label = f'$ {latex(label[0])} $, {y_unit_str}'

#                 new_obj.columns = label

#         return new_obj

#     def set_xlabel(self, label=None):
#         if isinstance(self.index, pd.MultiIndex):

#             if label == None:
#                 new_obj = self.copy()

#                 for label in new_obj.index.tolist():
#                     if label in self.__class__._units:

#                         x_unit_str = f'[${type(self)._units[label]}$]'

#                     else:
#                         x_unit_str = ''

#                 label = latex(self.index.name) + f'{x_unit_str}'

#                 new_obj.index.name = label

#         return new_obj


# #     def format_axis_label(self):

# #         new_obj = self.copy()

# #         new_obj.columns = self.set_ylabel
# #         new_obj.index = self.set_xlabel

# #         return new_obj

#     def set_legend(self, legend=[]):
#         if legend == None:
#             pass

#     def _filtter(self, filtter=None):
#         filtter = self.__class__._data_filtter

#         return filtter

#     def __call__(self):

#         filtter = self._filtter()
#         #print('on __call__', filtter)

#         #         display('old data',self)
#         #         display('new data',filtter(self))

#         return self.__class__(filtter(
#             self)).format_labels().format_index()  #.set_ylabel().set_xlabel()

#     def plot(self, *args, **kwargs):

#         if not 'ylabel' in kwargs:
#             kwargs['ylabel'] = self.__class__._ylabel

#         return super().plot(*args, **kwargs)


# class FFTSeriesFormatter(BaseSeriesFormatter):
#     _domain = Symbol('f')
#     r'''
#     Basic class for formatting data plots. It provides methods for setting options

#     Arguments
#     =========

#     Methods
#     =======

#     Example
#     =======

#     '''
#     _data_filtter = lambda frame: frame.copy()

#     @property
#     def _constructor(self):
#         return FFTSeriesFormatter

#     @property
#     def _constructor_expanddim(self):
#         return FFTFrameFormatter

#     @property
#     def _constructor_sliced(self):
#         return FFTSeriesFormatter

#     def format_index(self, domain=Symbol('f')):
#         if isinstance(self.index, pd.MultiIndex):

#             idx = self.index.tolist()

#         else:
#             idx = self.index
#         #print('idx', idx)
#         new_idx = idx.copy()
#         #print('new idx', new_idx)

#         new_obj = self.copy()
#         new_obj.index = new_idx

#         #print('new_obj.index', new_obj.index.name)

#         new_obj.index.name = Symbol('f')

#         #print('new_obj.index', new_obj.index.name)

#         new_obj.index.name = f'${self._match_unit(new_obj.index.name)}$'
#         #new_obj.index.name = 'cos'

#         ##print('new_obj.index.name',new_obj.index.name)

#         return new_obj


# class FFTFrameFormatter(BaseFrameFormatter):
#     _domain = Symbol('f')

#     _data_filtter = lambda obj: obj.to_frequency_domain().double_sided_rms(
#     ).truncate(0, 0.5)
#     r'''
#     Basic class for formatting data plots. It provides methods for setting options

#     Arguments
#     =========

#     Methods
#     =======

#     Example
#     =======

#     '''
#     @property
#     def _constructor(self):
#         return FFTFrameFormatter

#     @property
#     def _constructor_expanddim(self):
#         return FFTFrameFormatter

#     @property
#     def _constructor_sliced(self):
#         return FFTSeriesFormatter

#     def format_index(self, domain=Symbol('f')):
#         if isinstance(self.index, pd.MultiIndex):

#             idx = self.index.tolist()

#         else:
#             idx = self.index
#         #print('idx', idx)
#         new_idx = idx.copy()
#         #print('new idx', new_idx)

#         new_obj = self.copy()
#         new_obj.index = new_idx

#         #print('new_obj.index', new_obj.index.name)

#         new_obj.index.name = Symbol('f')

#         #print('new_obj.index', new_obj.index.name)

#         new_obj.index.name = f'${self._match_unit(new_obj.index.name)}$'
#         #new_obj.index.name = 'cos'

#         ##print('new_obj.index.name',new_obj.index.name)

#         return new_obj


# class PivotSeriesSummary(BaseSeriesFormatter):
#     _default_sep = ' \n '

#     @property
#     def _constructor(self):
#         return PivotSeriesSummary

#     @property
#     def _constructor_expanddim(self):
#         return PivotFrameSummary

#     @property
#     def _constructor_sliced(self):
#         return PivotSeriesSummary


# class PivotFrameSummary(BaseFrameFormatter):
#     #    _data_filtter = lambda frame: frame.abs().max().reset_index(level=1).pivot(columns=['level_1'])
#     #    _label_formatter = lambda entry: f'${latex(entry)}$'
#     _default_sep = ' \n '

#     @property
#     def _constructor(self):
#         return PivotFrameSummary

#     @property
#     def _constructor_expanddim(self):
#         return PivotFrameSummary

#     @property
#     def _constructor_sliced(self):
#         return PivotSeriesSummary

#     def _format_label(self, obj):
#         if self.__class__._label_formatter:
#             return self.__class__._label_formatter(obj)
#         else:
#             if isinstance(obj, Iterable):
#                 obj = obj
#             else:
#                 obj = obj,

#             return tuple([
#                 f'${self._match_unit(elem)}$'
#                 if isinstance(elem, (Expr, Eq)) else f'{elem}' for elem in obj
#                 if elem != 0
#             ])

#     def format_labels(self):
#         if isinstance(self.columns, pd.MultiIndex):

#             idx = self.columns.tolist()

#         else:
#             idx = self.columns
#         #print('idx', idx)
#         new_idx = []
#         for entry in idx:
#             #print('entry', entry, type(entry))
#             if entry in self.__class__._units and entry != 0:

#                 units = self.__class__._units

#                 #new_idx+=[ self._format_label(entry) + f'{units[entry]}']
#                 new_idx += [self._format_label(entry)]
#             elif entry != 0:
#                 new_idx += [self._format_label(entry)]

#         #print('new idx', new_idx)

#         self.__class__._ylabel = (list(idx)[0][0])

#         new_obj = self.copy()

#         new_obj.columns = pd.MultiIndex.from_tuples(new_idx)

#         cols_name = self.__class__._cols_name

#         if cols_name:
#             idx_name = self._format_label(cols_name)
#             new_obj.columns.name = idx_name

#         return new_obj

#     def format_index(self):
#         if isinstance(self.index, pd.MultiIndex):

#             idx = self.index.tolist()

#         else:
#             idx = self.index
#         #print('idx', idx)
#         new_idx = [entry.rhs.n(4) for entry in idx]

#         #print('new idx', new_idx)

#         #self.__class__._ylabel=f'${latex(list(idx)[0].lhs)}$'

#         new_obj = self.copy()
#         new_obj.index = new_idx

#         new_obj.index.name = f'${self._match_unit(list(idx)[0].lhs)}$'

#         return new_obj

#     def _filtter(self, filtter=None):
#         if self.columns.nlevels == 2:
#             filtter = lambda frame: frame.abs().max().reset_index(
#                 level=1).pivot(columns=['level_1'])
#         else:
#             filtter = lambda frame: frame.abs().max().reset_index().pivot(
#                 columns=['level_0', 'level_2'], index=['level_1'])[0]

#         return filtter


# class PivotPlotSeriesSummary(BaseSeriesFormatter):
#     _default_sep = ', '

#     @property
#     def _constructor(self):
#         return PivotPlotSeriesSummary

#     @property
#     def _constructor_expanddim(self):
#         return PivotPlotFrameSummary

#     @property
#     def _constructor_sliced(self):
#         return PivotPlotSeriesSummary


# class PivotPlotFrameSummary(BaseFrameFormatter):
#     #    _data_filtter = lambda frame: frame.abs().max().reset_index(level=1).pivot(columns=['level_1'])
#     #    _label_formatter = lambda entry: f'${latex(entry)}$'
#     _default_sep = ', '

#     @property
#     def _constructor(self):
#         return PivotPlotFrameSummary

#     @property
#     def _constructor_expanddim(self):
#         return PivotPlotFrameSummary

#     @property
#     def _constructor_sliced(self):
#         return PivotPlotSeriesSummary

#     def format_index(self):
#         if isinstance(self.index, pd.MultiIndex):

#             idx = self.index.tolist()

#         else:
#             idx = self.index
#         #print('idx', idx)
#         new_idx = [entry.rhs.n(4) for entry in idx]

#         #print('new idx', new_idx)

#         #self.__class__._ylabel=f'${latex(list(idx)[0].lhs)}$'

#         new_obj = self.copy()
#         new_obj.index = new_idx

#         new_obj.index.name = f'${self._match_unit(list(idx)[0].lhs)}$'

#         self.__class__._domain = list(idx)[0].lhs

#         return new_obj

#     def _filtter(self, filtter=None):
#         if self.columns.nlevels == 2:
#             filtter = lambda frame: frame.abs().max().reset_index(
#                 level=1).pivot(columns=['level_1'])
#         else:
#             filtter = lambda frame: frame.abs().max().reset_index().pivot(
#                 columns=['level_0', 'level_2'], index=['level_1'])[0]

#         return filtter


##################################################### COMMENTS FOR CHECK


# time_domain_summary_1 = [
#     'Przedstawione wyniki pokazują wpływ zmian badanego parametru \({param_name}\) na dynamikę wózka. ',
#     'Opracowane wykresy porównawcze pozwoliły na oszacowanie wpływu zmiennej \({param_name}\) na dynamikę rozpatrywanego układu. ',
#     'Zaprezentowane wyniki posłużyły ocenie wrażliwości modelowanego obiektu na badaną zmienną \({param_name}\). ',
#     'Przygotowane przebiegi służą zaprezentowaniu wyników symulacji numerycznych badanej zmiennej w funkcji wybranego parametru wpływowego \({param_name}\). ',
#     'Analizowany model wykazuje wrażliwość na zmianę rozważanego parametru \({param_name}\). ',
# ]

# time_domain_summary_2 = [
#     'Zaobserwowano wpływ badanego parametru \({x(t)}\) na drgania wózka, których maksymalna wartość amplitudy wyniosła {x(t)_max}. ',
#     'Analizowany przebieg parametru \({x(t)}\) osiąga wartość maksymalną równą {x(t)_max}. '
#     'Minimalna wartość dla przebiegu \({x(t)}\) jest równa {x(t)_min}. '
#     'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({x(t)}\) nie występują amplitudy o wartości mniejszej niż {x(t)_min} oraz większej niż {x(t)_max}. '
#     'Dokonując analizy wyników symulacji dla parametru \({x(t)}\) stwierdza się, że amplitudy drgań nie przekraczają wartości {x(t)_max}, a wartość mininalna wynosi {x(t)_min}. ',
# ]

# time_domain_summary_3 = [
#     'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi(t)}\) wynosi {varphi(t)_max}. ',
#     'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi(t)_max} dla współrzędnej \({varphi(t)}\). '
#     'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {varphi(t)_min}. ',
#     'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
# ]

# time_domain_summary_4 = [
#     'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({z(t)}\) wynosi {z(t)_max}. ',
#     'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {z(t)_max} dla współrzędnej \({z(t)}\). '
#     'Wartość minimalna przebiegu współrzędnej \({z(t)}\) jest równa {z(t)_min}. ',
#     'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({z(t)}\) nie przekracza {z(t)_max} i jest ograniczona z dołu przez {z(t)_min}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({z(t)}\) nie przekraczaa wartości {z(t)_max}, a wartość mininalna wynosi {z(t)_min}. ',
# ]

# time_domain_summary_5 = [
#     'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionej \({varphi_RC(t)}\) wynosi {varphi_RC(t)_max}. ',
#     'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się wartością maksymalną równą {varphi_RC(t)_max} dla współrzędnej \({varphi_RC(t)}\).'
#     'Wartość minimalna przebiegu współrzędnej \({varphi_RC(t)}\) jest równa {varphi_RC(t)_min}. ',
#     'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi_RC(t)}\) nie przekracza {varphi_RC(t)_max} i jest ograniczona z dołu przez {varphi_RC(t)_min}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi_RC(t)}\) nie przekraczaa wartości {varphi_RC(t)_max}, a wartość mininalna wynosi {varphi_RC(t)_min}. ',
# ]

# time_domain_summary_7 = [
#     'Zaobserwowane odczyty wynikające z wpływu \({param_name}\) występują w konketrnym czasie ruchu {t_val}. ',
#     'Odczytana maksymalna wartość amplitudy {max_val} wynika z działania \({param_name}\) wystąpiła konkretnie w czasie {t_val}. '
#     'Minimalny zaobserwowany wynik amplitudy {min_val} miał miejsce w {t_val} czasu trwania pomiaru. ',
# ]

# time_domain_summary_8 = [
#     'Wyniki symulacji pokazują, że wartość maksymalna dla współrzędnej uogólnionych są następujące: {max_val_list}. ',
#     'Analizowane sygnały reprezentujące odpowiedż czasową układu charakteryzują się następującymi wartościami maksymalnymi: {max_val_list}. '
#     'Wartość minimalna przebiegu współrzędnej \({varphi(t)}\) jest równa {min_val}. ',
#     'Przeprowadzone badania numeryczne pozwalają stwierdzić, że odpowiedź czasowa \({varphi(t)}\) nie przekracza {varphi(t)_max} i jest ograniczona z dołu przez {varphi(t)_min}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({varphi(t)}\) nie przekraczaa wartości {varphi(t)_max}, a wartość mininalna wynosi {varphi(t)_min}. ',
# ]

# freq_domain_summary = [
#     'Ponadto sprawdzono własności spektralne otrzymanych wyników symulacyjnych.',
#     'Według przeprowadzonych badań modelowych struktura częstotliwościowa nie zawiera komonentów większych niż {max_val_spec}.',
#     'Amplitudy drgań badanego parametru \({q_name}\) przyjmują maksymalną wartość równą {max_val_spec}.',
#     'Na podstawie analizy widmowej minimalna wartość częstotliwości rozpatrywanego parametru wynosi {min_val_spec}.',
# ]

# summary_bank_composed_model_gen = RandomDescription(
#     time_domain_summary_1,
#     time_domain_summary_2,
#     time_domain_summary_3,
#     time_domain_summary_4,
#     time_domain_summary_5,
# )

# summary_bank_chair_model_gen = RandomDescription(
#     time_domain_summary_1,
#     time_domain_summary_2,
#     time_domain_summary_3,
#     time_domain_summary_4,
#     #    time_domain_summary_5,
# )

# summary_bank_chair_model_gen = RandomDescription(
#     time_domain_summary_1,
#     time_domain_summary_2,
#     time_domain_summary_3,
#     time_domain_summary_4,
#     #    time_domain_summary_5,
# )

# summary_bank_drive_model_gen = RandomDescription(
#     time_domain_summary_1,
#     #     time_domain_summary_2,
#     #     time_domain_summary_3,
#     #     time_domain_summary_4,
#     time_domain_summary_5,
# )

# summary_bank_chair_move_model_gen = RandomDescription(
#     time_domain_summary_1,
#     time_domain_summary_2,
#     time_domain_summary_3,
#     #     time_domain_summary_4,
#     #   time_domain_summary_5,
# )

# summary_bank_drive_model = [
#     str(
#         RandomDescription(
#             time_domain_summary_1,
#             #             time_domain_summary_2,
#             #             time_domain_summary_3,
#             #             time_domain_summary_4,
#             time_domain_summary_5,
#         )) for obj in range(30)
# ]

# measurement_summary_1 = [
#     'Wykresy prezentowane są w formie pełnej, to znaczy od momentu włączenia do wyłączenia aparatury pomiarowej.',
#     'Początek wykresu odpowiada momentowi wyzwolenia pomiaru, natomiast jego koniec - wyłączeniu aparatury. ',
#     'Przedstawione wyniki odpowiadają całemu cyklowi pomiarowemu - od momentu włączenia, aż do wyłączenia aparatury. ',
#     'Prezentowane wykresy mają formę nieuciętą, to znaczy początek i koniec odpowiada włączeniu i wyłączeniu pomiaru. ',
#     'Na wykresach zaprezentowano przebiegi od momentu włączenia, aż do wyłączenia aparatury pomiarowej. ',
# ]

# measurement_summary_a_ox = [
#     'Dla czujnika \({a_ox}\), maksymalna wartość amplitudy przyśpieszenia wyniosła {a_ox_max}, natomiast minimalna - {a_ox_min}. ',
#     'W analizowanym przebiegu parametr \({a_ox}\) osiągnał wartość maksymalną równą {a_ox_max} przy wartości minimalnej {a_ox_min}. ',
#     'Minimalna wartość dla przebiegu \({a_ox}\) wyniosła {a_ox_min}, natomiast wartość maksymalna odpowiednio {a_ox_min}. ',
#     'Na podstawie przeprowadzonych badań stwierdzono, że w analizowanym przebiegu przyspieszeń \({a_ox}\) nie występują amplitudy o wartości mniejszej niż {a_ox_min} oraz większej niż {a_ox_max}. ',
#     'Dokonując analizy wyników przejazdu stwierdza się, że amplitudy przyśpieszeń \({a_ox}\) nie przekraczają wartości {a_ox_max}, a wartość mininalna wynosi {a_ox_min}.',
#     'Przeprowadzone badanie pomiarowe pozwalaja stwierdzić, że odpowiedź czasowa \({a_ox}\) nie przekracza {a_ox_max} i jest ograniczona z dołu przez {a_ox_min}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że amplituda przyspieszenia \({a_ox}\) nie przekraczaa wartości {a_ox_max}, a wartość mininalna wynosi {a_ox_min}. ',
# ]

# measurement_summary_a_oz = [
#     'Minimalna wartość amplitudy przyspieszenia z akcelerometru \({a_oz}\) wynosi {a_oz_min}, przy wartości maksymalnej {a_oz_max}. ',
#     'Przebieg parametru \({a_oz}\) osiąga wartość minimalną {a_oz_min} oraz maksymalną równą {a_oz_max}. ',
#     'Dla sygnału z akcelerometru \({a_oz}\) zarejestrowano wartość minimalną {a_oz_min}, natomiast maksimum odnotowano na poziomie {a_ox_max}. ',
#     'Z przeprowadzonej próby wynika, że przebieg sygnału \({a_oz}\) nie przekraacza przedziału wartości od {a_oz_min} do {a_oz_max}. ',
#     'Analizując przebieg sygnału z akcelerometru \({a_oz}\) stwierdza się, że amplituda maksymalna przyśpieszenia nie przekracza wartości {a_oz_max}, a wartość mininalna wynosi {a_oz_min}. ',
#     'Odpowiedź czasowa \({a_oz}\) nie przekracza wartości górnej {a_oz_max} oraz dolnej {a_oz_min}. ',
#     'W analizowanym przypadku, wartości graniczne sygnału z czujnika przyspieszeń \({a_oz}\) wyniosły {a_oz_min} (minimalna) oraz {a_oz_max} (maksymalna). ',
# ]

# measurement_summary_a_rz = [
#     'Sygnał \({a_rz}\) przyjął wartość minimialną {a_rz_min}, nie przekraczając górnej granicy {a_rz_max}. ',
#     'Analiza przebiegu sygnału \({a_rz}\) pozwala na określenie jego wartości minimalnej na poziomie {a_rz_min} oraz maksymalnej {a_rz_max}. ',
#     'Dla opisywanego przypadku, amplituda \({a_rz}\) nie przekracza granicy dolnej {a_rz_min} oraz granicy górnej {a_rz_max}. ',
#     'Akcelerometr \({a_rz}\) zarejestrował wartości przyspieszeń nie mniejsze niż {a_rz_min} i nie większe od {a_rz_max}. ',
#     'Z kolei z przebiegu sygnału \({a_rz}\) odczytano wartości minimalne i maksymalne odpowiednio {a_rz_min} oraz {a_rz_max}. ',
#     'Czujnik \({a_rz}\) zarejestrował sygnał w granicach od {a_rz_min} do {a_rz_max}. ',
#     'Dla czujnika przyspieszeń pionowych \({a_rz}\) odnotowano skrajne wartości amplitudy: minimalną {a_rz_min} oraz maksymalną {a_rz_max}. ',
# ]

# measurement_summary_a_rcz = [
#     'Akcelerometr umieszczony na wahaczu RapidChair (\({a_rcz}\)) zarejestrował wartość minimialną {a_rcz_min}, przy maksymalnej {a_rcz_max}. ',
#     'Czujnik rejestrujący sygnał przyspieszenia \({a_rcz}\) odnotował wartości w przedziale od {a_rcz_min} do {a_rz_max}. ',
#     'Analizowany przejazd wywołał na czujniku (\({a_rcz}\)) sygnał, dla którego amplituda maksymalna przyspieszenia nie przekracza wartości {a_rcz_max}, a wartość minimalna wynosi {a_rcz_min}. ',
#     'Przebieg wkresu przyspieszenia \({a_rcz}\) przyjumuje wartości nie mniejsze niż {a_rcz_min} i nie większe od {a_rcz_max}. ',
#     'Akcelerometr (\({a_rcz}\)) umieszczony na wahaczu napędu RapidChair zarejestrował przyspieszenia w zakresie od {a_rcz_min} do {a_rcz_max}. ',
#     'Odczytanie danych dla czujnika umieszczonego na wahaczu napędu RapidChair umożliwia określenie skrajnych wartości amplitudy przespieszeń pomiędzy {a_rcz_min} a {a_rcz_max}. ',
#     'Sygnał \({a_rcz}\) przyjmuje wartości nie mniejsze niż {a_rcz_min} oraz nie większe niż {a_rcz_max}.'
# ]

# meas_summary_bank = RandomDescription(measurement_summary_1,
#                                       measurement_summary_a_ox,
#                                       measurement_summary_a_oz,
#                                       measurement_summary_a_rz,
#                                       measurement_summary_a_rcz)

# measurement_summary_4 = [
#     'Zaobserwowano wpływ badanego parametru \({a_rx}\) na drgania wózka, których maksymalna wartość amplitudy porecntu (%) przyśpieszenia ziemskiego wyniosła {a_rx_max}. ',
#     'Analizowany przebieg parametru \({a_rx}\) osiąga wartość maksymalną równą {a_rx_max}. ',
#     'Minimalna wartość dla przebiegu \({a_rx}\) jest równa {a_rx_min}. ',
#     'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_rx}\) nie występują amplitudy o wartości mniejszej niż {a_rx_min} oraz większej niż {a_rx_max}. ',
#     'Dokonując analizy wyników symulacji dla parametru \({a_rx}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_rx_max}, a wartość mininalna wynosi {a_rx_min}.',
#     'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_rx}\) nie przekracza {a_rx_idmax} i jest ograniczona z dołu przez {a_rx_idmin}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_rx}\) nie przekraczaa wartości {a_rx_idmax}, a wartość mininalna wynosi {a_rx_idmin}. ',
# ]

# measurement_summary_5 = [
#     'Zaobserwowano wpływ badanego parametru \({a_rz}\) na drgania wózka, których maksymalna wartość przyśpieszenia wyniosła {a_rz_max}. ',
#     'Analizowany przebieg parametru \({a_rz}\) osiąga wartość maksymalną równą {a_rz_max}. ',
#     'Minimalna wartość dla przebiegu \({a_rz}\) jest równa {a_ox_min}. ',
#     'Na podstawie przeprowadzonych badań numerycznych można stwierdzić, że w analizowanym przebiegu parametru \({a_rz}\) nie występują amplitudy o wartości mniejszej niż {a_rz_min} oraz większej niż {a_rz_max}. ',
#     'Dokonując analizy wyników symulacji dla parametru \({a_rz}\) stwierdza się, że amplitudy przyśpieszeń nie przekraczają wartości {a_rz_max}, a wartość mininalna wynosi {a_rz_min}.',
#     'Przeprowadzone badania pomiarowe pozwalają stwierdzić, że odpowiedź czasowa \({a_rz}\) nie przekracza {a_rz_idmax} i jest ograniczona z dołu przez {a_rz_idmin}. ',
#     'Bazując na analizie przeprowadzonych badań zaobserwowano, że współrzędna \({a_rz}\) nie przekraczaa wartości {a_rz_idmax}, a wartość mininalna wynosi {a_rz_idmin}. ',
# ]

# # summary_bank_composed_model_measurements = [
# #     str(
# #         RandomDescription(
# #             measurement_summary_1,
# #             measurement_summary_2,
# #             measurement_summary_3,
# #             measurement_summary_4,
# #             measurement_summary_5,
# #         )) for obj in range(30)
# # ]

# sec_introduction_bank_0 = [
#     'Zaprezentowane symulacje wykonano dla następujących danych: {given_data}.',
#     'Do wykonania symulacji użyto danych jak następuje: {given_data}',
#     'Zakres danych użytych do procesu symulacji jest następujący: {given_data}.',
#     'Do opracowania wyników symulacyjnych posłużono się niniejszymi danymi: {given_data}.'
# ]
# sec_introduction_bank_1 = [
#     'Wykorzystując przygotowany model dynamiczny, przeprowadzono serię symulacji mających na celu ocenę wrażliwości modelu na zmiany wybranych parametrów. ',
#     'Przeprowadzenie analizy numerycznej wymaga znajomości modelu dynamicznego oraz zdeterminowania parametrów oddziałujących na zachowanie układu wskutek ich zmian. ',
#     'Wykorzystując przygotowane środowisko obliczeniowe, wykonano szereg symulacji, opracowano dane numeryczne i przedstawiono je w postaci następujących wykresów. ',
#     'Celem zbadania i oceny dynamiki rozpatrywanego systemu przeprowadzono symulacje numeryczne przedstawijące przebiegi zmiennych układu w funkcji czasu. '
#     'W oparciu o wybrany model dynamiczny obiektu, wykonano analizy numeryczne zaprezentowane na kolejnych wykresach. '
# ]

# sec_introduction_bank_2 = [
#     'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
#     'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
#     'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
#     'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
# ]

# sec_introduction_bank_3 = [
#     'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu możliwe będzie wprowadzenie zmian, które wpłyną na poprawę działania modelu. ',
#     'Charakter zmian przebiegów poszczególnch współrzędnych może posłużyć do opracowania wniosków, na podstawie których możliwe będzie ulepszenie funkcjonowania rozważanego układu. ',
#     'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu maja na celu ułatwienie analizy zachowania modelu dynamicznego oraz wprowadzenie ewentualnych poprawek na podstawie dokonanych obserwacji i spostrzeżeń. ',
#     'Dostrzegając wzajemne zależności między dynamicznymi parametrami wózka możliwe będzie wprowadzenie takich zmian, które w korzystny sposób wpłyną na odpowiedź układu. '
# ]

# sec_intro_composed_model_gen = RandomDescription(
#     sec_introduction_bank_1,
#     sec_introduction_bank_2,
#     sec_introduction_bank_3,
#     #sec_introduction_bank_0
# )

# sec_intro_composed_model = [
#     str(
#         RandomDescription(sec_introduction_bank_1, sec_introduction_bank_2,
#                           sec_introduction_bank_3)) for obj in range(30)
# ]

# introduction_bank_meas_2 = [  #poprawione
#     'W Tabeli {nr_rys} przedstawiano zakres zmienności parametrów, dla których przeprowadzono pomiary drgań struktury wózka.',
#     'Zestwienie parametrów, co do których założono potencjalny wpływ na drgania wózka, wraz z zakresami zmiennosci przedstawiono w Tabeli {nr_rys}.',
#     'W Tabeli {nr_rys} pokazano parametry układu wózek-pasażer, których zmienność założono jako potencjalnie wpływającą na charakter drgań.',
#     'Tabela {nr_rys} stanowi zestawienie zakresu wartości parametrów przyjętych jako potencjalnie wpływających na drgania wózka.',
# ]

# introduction_bank_meas_1 = [  #poprawione
#     'W niniejszym podrozdziale opisano próbę {string}. ',
#     'Opisane w dalszej kolejności testy dotyczą {string}. ',
#     'Następujące próby wykonano dla przypadku {string}. ',
#     'Przedmiotem próby opisywanej w tym podrozdziale jest przypadek {string}. ',
# ]

# introduction_bank_meas_3 = [  #poprawione
#     'Pomiary wykonano z uwzględnieniem napędu RapidChair w formie wahacza wleczonego. ',
#     'Do wózka inwalidzkiego dołączono napęd RapidChair w formie wahacza wleczonego. ',
#     'Badany układ stanowił wózek inwalidzki z dołączonym wahaczem wleczonym RapidChair. ',
#     'Pomiary wykonywano na układzie wózka inwalidzkiego z dołączonym wahaczem wleczonym RapidChair. ',
# ]

# introduction_bank_meas_4 = [  #poprawione (nie trzeba zmieniać zdań)
#     'Łącznie zaprezentowano wyniki {entries_no} przejazdów, które uznano za miarodajne. ',
#     'W niniejszej częsci raportu zamieszczono zapis z {entries_no} jazd testowych, co do których nie wyrażono wątpliwości względem otrzymanych wnyików. ',
#     'W raporcie zamieszczono wyniki jedynie z tych z przejazdów próbnych, których wyniki uznano za poprawne (łącznie {entries_no}). ',
#     'W dalszej części przedstawiono jedynie przebiegi drgań z {entries_no} pomiarów uznanych za miarodajne. '
# ]

# meas_intro_composed_model = RandomDescription(introduction_bank_meas_1,
#                                               introduction_bank_meas_2,
#                                               introduction_bank_meas_3,
#                                               introduction_bank_meas_4)

# ending_bank_meas_1 = [  #poprawione
#     'Dla każdego przypadku masy układu przeprowadzono co najmniej trzy przejazdy, dzięki czemu zmniejszono ryzyko uwzględnienia losowych błędów. ',
#     'Zgodnie z założeniami metodyki badawczej, aby zminimalizować losowe błędy, dla każdej masy wykonano co najmniej trzy próby. ',
#     'W celu zminimalizowania wpływu losowych błędów, wykonano co najmniej trzy przjeazdy testowe dla każdego przypadku masy układu. ',
#     'Pomiary wykonywano co najmniej trzykrotnie dla każdej przewidzianej masy, aby zminimalizować wpływ błędów losowych podczas analizy wyników. ',
# ]

# ending_bank_meas_2 = [  #poprawione
#     'Rozróżnienie wyników zapewniono poprzez uwzględnienie w nazewnictwie indeksu z numerem próby.',
#     'Wyodrębnienie poszczególnych przypadków umożliwiono stosując w nazwach indeksy z numerami prób.',
#     'Stąd indeks liczbowy przy nazwie próby oznacza kolejny przejazd z daną masą.',
#     'Kolejne przejazdy z daną masą rozróżnia się dzięki zastosowaniu indeksu liczbowego.',
# ]

# ending_bank_meas_3 = [  #poprawione (nie trzeba zmieniać zdań)
#     'Przebiegi czasowe zestawiono w kolejnej sekcji dokumentu.',
#     'Zestawienie wykresów czasowych zaprezentowano w następnym podrozdziale.',
#     'Wyniki przeprowadzonych testów w postaci wykresów czasowych znajdują się w kolejnym podrozdziale.',
#     'Zestawienie przebiegów czasowych zamieszczono w następnej sekcji dokumentu.'
# ]

# meas_ending_composed_model = RandomDescription(ending_bank_meas_1,
#                                                ending_bank_meas_3)

# intro_bank_meas_1 = [  #poprawione
#     'Na wykresie {nr_rys} przedstawiano zmiany wartości przyspieszeń drgań charakteryzujących ruch wózka.',
#     'Na rysunku {nr_rys} przedstawiono charakter przebiegu amplitud przyspieszeń drgan wózka. ',
#     'Wykres {nr_rys} pokazuje zmienność wartości przyspieszeń drgań wózka w trakcie przejazdu pomiarowego.',
#     'Rysunek {nr_rys} prezentuje wykres opisujący zmienność amplitudy przyspieszeń drgań wózka w trakcie jazdy pomiarowej.',
# ]

# intro_bank_meas_2 = [  #poprawione
#     'Pomiar został przeprowadzony dla następujących danych: {given_data}.',
#     'Zaprezentowane wyniki pomiarowe otrzymano dla danych równych: {given_data}.',
#     'Badanie pomiarowe zostało przeprowadzone dla następujących danych: {given_data}.',
#     'Wyniki pomiaru otrzymano dla: {given_data}.',
# ]

# intro_bank_meas_3 = [  #poprawione (nie trzeba zmieniać zdań)
#     'Zarejestrowane sygnały pochodzą z czujników umieszczonych w osi wózka (gdzie \({a_ox}\) - przyspieszenia wzdłużne, \({a_oz}\) - pionowe), na podnóżku (kierunek pionowy, \({a_rz}\)) oraz w wahaczu napędu RapidChair \({a_rcz}\).',
#     'Przebiegi uzyskano z akcelereometrów zamontowanych kolejno w osi wózka (gdzie \({a_ox}\) - przyspieszenia wzdłużne, \({a_oz}\)) - pionowe), na podnóżku (kierunek pionowy, \({a_rz}\))) oraz w wahaczu napędu RapidChair \({a_rcz}\)).',
#     'Oznaczenia \({a_ox}\), \({a_oz}\), \({a_rz}\) oraz \({a_rcz}\) odnoszą się odpowiednio do sygnałów z czujników w osi wózka (kolejno wzdłużne i pionowe), na podnóżku oraz w wahaczu RapidChar.',
#     'Przyjęto następujące oznaczenia: \({a_ox}\) - przyspieszenia wzdłużne czujnika umieszcoznego w osi wózka, \({a_oz}\) - przyspieszenia pionowe czujnika umieszcoznego w osi wózka, \({a_rz}\) - przyspieszenia pionowe czujnika na podnóżku oraz \({a_rcz}\) - przyspieszenia pionowe czujnika w wahaczu RapidChair.',
# ]

# intro_bank_meas_composed_model = RandomDescription(intro_bank_meas_1,
#                                                    intro_bank_meas_2,
#                                                    intro_bank_meas_3)
# introduction_bank_1 = [
#     'Na wykresie {nr_rys} przedstawiano zmiany wielkości dynamicznych charakteryzujących ruch obiektu. Po przeanalizowaniu można zanotować wzajemną zależność poszczególnych wielkości dynamicznych.',
#     'Charakter przebiegu wartości dynamicznych układu został przedstawiony na rysunku {nr_rys}.',
#     'Wykres {nr_rys} pokazuje zmienność parametrów dynamicznych wózka w trakcie symulownego przejazdu.',
#     'Rysunek {nr_rys} prezentuje wykres opisujący zmienność parametrów dynamicznych obiektu w trakcie jazdy po zamodelowanym torze.',
# ]

# introduction_bank_2 = [
#     'Symulacja została przeprowadzona dla następujących danych: {given_data}.',
#     'Zaprezentowane wyniki numeryczne otrzymano dla danych równych:{given_data}.',
#     'Eksperyment numeryczny został przeprowadzona dla następujących danych:{given_data}.',
#     'Wyniki stmulacji numerycznych otrzymano dla:{given_data}.',
# ]

# introduction_bank_3 = [
#     'Na podstawie spostrzeżeń wynikających z obserwacji odpowiedzi czasowych analizowanego układu stwierdza się oscylacyjny charakter odpowiedzi układu. Dynamika systemu ma stabilny charakter',
#     'Charakter zmian przebiegów poszczególnch współrzędnych jest stabilny. Otrzymane sygnały mieszczą się w zakresie dopuszczalnych uwarunkowaniami fizycznymi.',
#     'Zaprezentowane wykresy przedstawiające odpowiedzi czasowe układu mają drganiowy charakter i są zgodne z oczekiwaniami.',
#     'Bazując na otrzymanych wynikach można stwierdzić wzajemną zależność poszczególnych wielkości dynamicznych. Dodatkowo wielkości modelu stabilizują się w czasie. '
# ]

# intro_bank_composed_model_gen = RandomDescription(introduction_bank_1,
#                                                   introduction_bank_2,
#                                                   introduction_bank_3)

# intro_bank_composed_model = [
#     str(
#         RandomDescription(introduction_bank_1, introduction_bank_2,
#                           introduction_bank_3)) for obj in range(30)
# ]

# meas_comparation_bank_1 = [
#     'Na Rysunku {nr_rys} przedstawiono zestawienie wartości maksymalnych przyspieszeń drgań dla poszczególnych czujników w funkcji masy pasażera.',
#     'Wykres {nr_rys} reprezentuje zmienność wartości maksymalnych przyspieszeń drgań dla każdego z czujników, w odniesieniu do masy pasażera.',
#     'Rysunek {nr_rys} reprezentuje zestawienie zmienności maksymalnych amplitud przyspieszeń w zależności od masy testującego.'
# ]
# meas_comparation_bank_2 = [
#     'Między innymi na jego podstawie, uwzględniając informacje o poziomach minimalnych i średnich, dokonano oceny wpływu masy układu na poziom drgań dla danej próby. ',
#     'Posłużył on, wraz z informacjami o wartościach minimalnych i średnich, do określenia, w jakim stopniu masa układu wpływa na charakter drgań w zakresie całej przeprowadzonej próby. ',
#     'Wspólnie z danymi o wartościach średnich i minimalnych stanowił on podstawę do określenia wpływu masy układu na ogólny poziom amplitud przyspieszeń drgań na rozpatrywanej nawierzchni.'
# ]
# meas_comparation_bank_3 = [
#     'Opisane w dalszej kolejności poziomy wpływu masy na wielkość drgań odnoszą się do amplitud w punktach mocowania poszczególnych akcelerometrów.',
#     'Ogólny poziom drgań oceniano poprzez poziom wpływu masy na przyspieszenia drgań w każdym z punktów mocowania czujnika z osobna. ',
#     'Dalszej oceny wpływu masy pasaera {param_name} na poziom drgań dokonywano w każdym z punktów mocowania czujników z osobna.'
# ]
# meas_comparation_composed_model = RandomDescription(meas_comparation_bank_1,
#                                                     meas_comparation_bank_2,
#                                                     meas_comparation_bank_3)

# conclusion_bank_x = [
#     'Zauważa się {x(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({x(t)}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {x(t)_max} oraz {Derivative(x(t), t)_max} dla {x(t)_idmax} i {Derivative(x(t), t)_idmax}. ',
#     'Zaobserwowano {x(t)_inf} parametru \({param_name}\) na maksymalną wartość pokonanej drogi oraz osiąganą wartość prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {x(t)_max} oraz {Derivative(x(t), t)_max} dla maksymalnej wartości badanego parametru równej {x(t)_idmax}. ',
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({x(t)}\), gdzie największa wartość pokonanej drogi to {x(t)_max}. W konsekwencji zaobserwowano {x(t)_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(x(t), t)}\), dla której minimalna wartość wynosi {Derivative(x(t), t)_min}, a największą osiąganą wartością jest {Derivative(x(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(x(t), t)_idmin} oraz {Derivative(x(t), t)_idmax}. ',
# ]

# measurement_conclusion_bank_a_ox = [
#     'Zauważa się {a_ox_inf} parametru \({param_name}\) na poziom przyspieszeń \({a_ox}\). ',
#     'Zaobserwowano {a_ox_inf} parametru masy pasażera na maksymalną wartość przyspieszeń drgań z akcelerometru \({a_ox}\). ',
#     'Na podstawie wykresu oceniono, że parametr \({param_name}\) ma {a_ox_inf} na poziom drgań rejestrowanych przez czujnik \({a_ox}\). '
# ]
# measurement_conclusion_bank_a_oz = [
#     'W przypadku sygnału \({a_oz}\) oceniono, że masa pasażera {param_name} ma {a_oz_inf} na ogólny poziom drgań. ',
#     'Dla sygnału z akcelerometru \({a_oz}\) odnotowano {a_oz_inf} masy pasażera na amplitudę drgań. ',
#     'Dla przypadku sygnału z czujnika \({a_oz}\) oceniono {a_oz_inf} masy pasażera na poziom drgań. '
#     #'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rx}\), gdzie największa wartość pokonanej drogi to {a_rx_max}. W konsekwencji zaobserwowano {a_rx_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(a_rx, t)}\), dla której minimalna wartość wynosi {Derivative(a_rx, t)_min}, a największą osiąganą wartością jest {Derivative(a_rx, t)_max} odpowiednio dla wartości parmametru: {Derivative(a_rx, t)_idmin} oraz {Derivative(a_rx, t)_idmax}. ',
# ]
# measurement_conclusion_bank_a_rz = [
#     'Wykres reprezentuje {a_rz_inf} parametru na drgania pionowe przy podnóżku \({a_rz}\). ',
#     'Analizując wykres dla drgań pionowych podnózka \({a_rz}\) odnotowano {a_rz_inf} masy {param_name} na poziom amplitudy przyspieszeń. ',
#     'Reprezentacja porównawcza sygnału \({a_rz}\) wskazuje na {a_rz_inf} parametru {param_name} na przebieg drgań w jego funkcji. '
#     #'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rx}\), gdzie największa wartość pokonanej drogi to {a_rx_max}. W konsekwencji zaobserwowano {a_rx_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(a_rx, t)}\), dla której minimalna wartość wynosi {Derivative(a_rx, t)_min}, a największą osiąganą wartością jest {Derivative(a_rx, t)_max} odpowiednio dla wartości parmametru: {Derivative(a_rx, t)_idmin} oraz {Derivative(a_rx, t)_idmax}. ',
# ]
# measurement_conclusion_bank_a_rcz = [
#     'Ostatni przypadek - czujnika na wahaczu napędu RapidChair (\({a_rcz}\)) wskazuje na {a_rcz_inf} {param_name} na rejestrowane przez ten akcelerometr sygnały.',
#     'Rozpoznanie przeprowadzone dla akcelerometru (\({a_rcz}\)) wykazuje {a_rcz_inf} rozpatrywanego parametru na poziom drgań w tym punkcie. ',
#     'Ostatnie rozpoznanie, przeprowadzone dla czujnika na wahaczu napędu RapidCHair (sygnał \({a_rcz}\)) ukazuje {a_rcz_inf} masy pasażera na drgania struktury. '
#     #'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rx}\), gdzie największa wartość pokonanej drogi to {a_rx_max}. W konsekwencji zaobserwowano {a_rx_inf} analizowanej zmiennej na wartość prędkości liniowej \({Derivative(a_rx, t)}\), dla której minimalna wartość wynosi {Derivative(a_rx, t)_min}, a największą osiąganą wartością jest {Derivative(a_rx, t)_max} odpowiednio dla wartości parmametru: {Derivative(a_rx, t)_idmin} oraz {Derivative(a_rx, t)_idmax}. ',
# ]
# measurement_conclusion_bank_composed_model = RandomDescription(
#     measurement_conclusion_bank_a_ox, measurement_conclusion_bank_a_oz,
#     measurement_conclusion_bank_a_rz, measurement_conclusion_bank_a_rcz)

# conclusion_bank_varphi_rc = [
#     'Zaobserwowano {varphi_RC(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań i prędkość kątową napędu RC. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi_RC(t)_min} oraz większej niż {varphi_RC(t)_max} odpowiednio dla wartości parametru: {varphi_RC(t)_idmin} oraz {varphi_RC(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_max}.',
#     'Zmiana \({param_name}\) ma {varphi_RC(t)_inf} na wartość drgań i prędkość kątową napędu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi_RC(t)_min} oraz większej niż {varphi_RC(t)_max} odpowiednio dla wartości parametru: {varphi_RC(t)_idmin} oraz {varphi_RC(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_min}. ',
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({varphi_RC(t)}\), gdzie największa wartość amplitudy to {varphi_RC(t)_max}, a najmniejsza {varphi_RC(t)_min}. W konsekwencji zaobserwowano {varphi_RC(t)_inf} analizowanej zmiennej na wartość prędkości kątową napędu RC. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi_RC(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi_RC(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(varphi_RC(t), t)_idmin} oraz {Derivative(varphi_RC(t), t)_idmax}. ',
# ]

# conclusion_bank_varphi_zwrc = [
#     'Zaobserwowano {z_wrc(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań pionowych i prędkość napędu RC. Przemieszczenia pionowe nie przyjmują wartości mniejszej niż {z_wrc(t)_min} oraz większej niż {z_wrc(t)_max} odpowiednio dla wartości parametru: {z_wrc(t)_idmin} oraz {z_wrc(t)_idmax}. Dla prędkości drgań pionowych napędu minimalna wartość amplitudy to {Derivative(z_wrc(t), t)_min}, a największą osiąganą wartością jest {Derivative(z_wrc(t), t)_max}. ',
#     'Zmiana \({param_name}\) ma {z_wrc(t)_inf} na wartość drgań i prędkość kątową analizowanego układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {z_wrc(t)_min} oraz większej niż {z_wrc(t)_max} odpowiednio dla wartości parametru: {z_wrc(z_wrct)_idmin} oraz {z_wrc(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(z_wrc(t), t)_min}, a największą osiąganą wartością jest {Derivative(z_wrc(t), t)_min}. ',
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({varphi(t)}\), gdzie największa wartość amplitudy to {varphi_RC(t)_max}, a najmniejsza {varphi(t)_min}. W konsekwencji zaobserwowano {varphi(t)_inf} analizowanej zmiennej na wartość prędkości drgań kątowych badanego układu, dla którego minimalna wartość amplitudy wynosi {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(varphi(t), t)_idmin} oraz {Derivative(varphi(t), t)_idmax}. '
# ]

# conclusion_bank_phi = [
#     'Zmiana \({param_name}\) ma {varphi(t)_inf} na wartość drgań i prędkość kątową układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi(t)_min} oraz większej niż {varphi(t)_max} odpowiednio dla wartości parametru: {varphi(t)_idmin} oraz {varphi(t)_idmax}. Dla prędkości kątowej wózka minimalna wartość amplitudy to {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_min}. ',
#     'Zaobserwowano {varphi(t)_inf} rozpatrywanego parametru - \({param_name}\) na wartość drgań i prędkość kątową układu. Przemieszczenia kątowe nie przyjmują wartości mniejszej niż {varphi(t)_min} oraz większej niż {varphi(t)_max} odpowiednio dla wartości parametru: {varphi(t)_idmin} oraz {varphi(t)_idmax}. Dla prędkości kątowej napędu minimalna wartość amplitudy to {Derivative(varphi(t), t)_min}, a największą osiąganą wartością jest {Derivative(varphi(t), t)_max}. ',
#     'Zauważa się {varphi(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({varphi(t)}\) oraz odpowiadającej temu przemieszczniu kątowemu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {varphi(t)_max} oraz {Derivative(varphi(t), t)_max} dla {varphi(t)_idmax} i {Derivative(varphi(t), t)_idmax}. ',
# ]

# conclusion_bank_z = [
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({z(t)}\), gdzie największa wartość amplitudy to {z(t)_max}, a najmniejsza {z(t)_min}. W konsekwencji zaobserwowano {z(t)_inf} analizowanej zmiennej na wartość drgań prędkości drgań pionowych, dla których minimalna wartość amplitudy to {Derivative(z(t), t)_min}, a największą osiąganą wartością jest {Derivative(z(t), t)_max} odpowiednio dla wartości parmametru: {Derivative(z(t), t)_idmin} oraz {Derivative(z(t), t)_idmax}. ',
#     'Zauważa się {z(t)_inf} zmienności parametru \({param_name}\) dla współrzędnej \({z(t)}\) oraz odpowiadającej temu przemieszczniu prędkości. Stwierdzono wpływ badanego parametru, gdzie maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {z(t)_max} oraz {Derivative(z(t), t)_max} dla {z(t)_idmax} i {Derivative(z(t), t)_idmax}. ',
#     'Zaobserwowano {z(t)_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości. Maksymalne wartości dla wymienionych współrzędnych przyjmują odpowiednio {z(t)_max} oraz {Derivative(z(t), t)_max} dla maksymalnej wartości badanego parametru równej {z(t)_idmax}. ',
# ]

# measurement_conclusion_bank_z = [
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_oz}\), gdzie największa wartość amplitudy to {a_oz_max}, a najmniejsza {a_oz_min}.',
#     'Zauważa się {a_oz_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_oz}\) oraz odpowiadającej temu przemieszczniu prędkości.',
#     'Zaobserwowano {a_oz_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości.'
# ]

# measurement_conclusion_bank_z_rc = [
#     'Zmianę dynamiki pod wpływem zmienności parametru \({param_name}\) obserwuje się dla \({a_rz}\), gdzie największa wartość amplitudy to {a_rz_max}, a najmniejsza {a_rz_min}.',
#     'Zauważa się {a_rz_inf} zmienności parametru \({param_name}\) dla współrzędnej \({a_rz}\) oraz odpowiadającej temu przemieszczniu prędkości.',
#     'Zaobserwowano {a_rz_inf} parametru \({param_name}\) na maksymalną wartość drgań pionowych wózka oraz osiąganą wartość ich prędkości.',
# ]

# conclusion_bank_no_impact = [
#     'Nie odmotowano natomiast istotnego wpływu tego parametru na pozostałe stopnie swobody układu. ',
#     'Nie zauważono kluczowych zmian pozostałych zmiennych układu pod wpływem badanych wartości rozpatrywanego parametru. ',
#     'Pozostałe zmienne układu nie wykazują wpływu oddziałującego w sposób istoty na dynamikę rozpatrywanego systemu. ',
#     'Badany system dynamiczny nie wykazuje szczególnej wrażlowości na zmiany innych analizowanych parametrów. ',
# ]

# measurement_conclusion_bank_no_impact = [
#     'Nie odmotowano natomiast istotnego wpływu tego parametru na pozostałe stopnie swobody wózka. ',
#     'Nie zauważono kluczowych zmian pozostałych zmiennych wózka pod wpływem badanych wartości rozpatrywanego parametru. ',
#     'Pozostałe zmienne układu nie wykazują wpływu oddziałującego w sposób istoty na dynamikę rozpatrywanego wózka. ',
#     'Badany wózek nie wykazuje szczególnej wrażlowości na zmiany innych analizowanych parametrów. ',
# ]

# measurement_conclusion_bank_summary = [
#     'Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty pomiarów. Dla uzyskanych danych liczbowych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego wózka. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).'
# ]

# conclusion_bank_composed_model_gen = RandomDescription(
#     conclusion_bank_varphi_rc, conclusion_bank_x, conclusion_bank_z,
#     conclusion_bank_phi, conclusion_bank_no_impact)

# conclusion_bank_composed_model = [
#     str(
#         RandomDescription(conclusion_bank_varphi_rc, conclusion_bank_x,
#                           conclusion_bank_z, conclusion_bank_phi,
#                           conclusion_bank_no_impact)) for obj in range(30)
# ]

# conclusion_bank_chair_model_gen = RandomDescription(conclusion_bank_x,
#                                                     conclusion_bank_z,
#                                                     conclusion_bank_phi,
#                                                     conclusion_bank_no_impact)

# conclusion_bank_chair_move_model_gen = RandomDescription(
#     conclusion_bank_x, conclusion_bank_phi, conclusion_bank_no_impact)

# conclusion_bank_chair_model = [
#     str(
#         RandomDescription(conclusion_bank_x, conclusion_bank_z,
#                           conclusion_bank_phi, conclusion_bank_no_impact))
#     for obj in range(30)
# ]

# # measurement_conclusion_bank_chair_model = [
# #     str(
# #         RandomDescription(measurement_conclusion_bank_x,
# #                           measurement_conclusion_bank_x_rc,
# #                           measurement_conclusion_bank_z,
# #                           measurement_conclusion_bank_z_rc,
# #                           measurement_conclusion_bank_no_impact))
# #     for obj in range(30)
# # ]

# conclusion_bank_drive_model_gen = RandomDescription(conclusion_bank_varphi_rc,
#                                                     conclusion_bank_no_impact)

# conclusion_bank_drive_model = [
#     str(RandomDescription(conclusion_bank_varphi_rc,
#                           conclusion_bank_no_impact)) for obj in range(30)
# ]

# # measurement_conclusion_bank_drive_model = [
# #     str(
# #         RandomDescription(measurement_conclusion_bank_x_rc,
# #                           measurement_conclusion_bank_no_impact))
# #     for obj in range(30)
# # ]

# analysis_intro_ending = '''Tabela {nr_rys} przedstawia zakres parametrów przyjętych do wykonania symulacji numerycznych. Oceniono, że przyjęty zakres badnia odpowiada możliwym do uzyskania w praktyce wartością i przeprowadzone symulacje będę dobrze reprezentować dynamikę układu.
# '''

# simulations_summary_str = ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty symulacji numerycznych. Dla uzyskanych danych symulacyjnych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego układu. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''


# class PlottedData(Figure):
#     _latex_name = 'figure'

#     def __init__(self,
#                  numerical_data,
#                  fig_name='Name',
#                  *args,
#                  units_dict=None,
#                  preview=False,
#                  position=None,
#                  **kwargs):
#         super().__init__(position=position, **kwargs)

#         self._numerical_data = numerical_data
#         self.fig_name = str(fig_name)
#         self._latex_name = 'figure'
#         super()._latex_name
#         self.preview = preview
#         self._units_dict = units_dict

#     def add_data_plot(self,
#                       numerical_data=None,
#                       xlabel=None,
#                       ylabel=None,
#                       grid=True,
#                       subplots=True,
#                       num_yticks=None,
#                       fontsize=None,
#                       figsize=(10, 4)):
#         numerical_data = self._numerical_data

#         ax = numerical_data.plot(subplots=subplots,
#                                  figsize=figsize,
#                                  ylabel=ylabel,
#                                  xlabel=xlabel,
#                                  grid=grid,
#                                  fontsize=fontsize)

#         #plt.xlim(numerical_data.index[0],numerical_data.index[-1])

#         if num_yticks != None:
#             ticks_no = num_yticks
#             for axl in ax:

#                 ylimit = axl.set_ylim(bottom=round(np.floor(
#                     axl.get_ylim()[0])),
#                                       top=round(np.ceil(axl.get_ylim()[1])))

#                 axl.set_yticks(np.linspace(ylimit[0], ylimit[1], ticks_no))

#                 axl.plot()
#         #ax=solution_tmp.plot(subplots=True)

#         if self._units_dict:
#             label_formatter = lambda sym: '$' + vlatex(
#                 sym) + '$' + '[${val:~L}$]'.format(val=self._units_dict[sym])
#         else:
#             label_formatter = lambda sym: '$' + vlatex(sym) + '$'

#         label_formatter_without_SI = lambda sym: '$' + vlatex(sym) + '$'

#         if subplots:
#             ([
#                 ax_tmp.legend([label_formatter(sym)], loc='lower right')
#                 for ax_tmp, sym in zip(ax, numerical_data.columns)
#             ])

#             ([
#                 ax_tmp.set_ylabel(
#                     label_formatter_without_SI(
#                         sym)  #.replace(r '\' , ' ')#.replace( '\\' ,' ' )
#                 ) for ax_tmp, sym in zip(ax, numerical_data.columns)
#             ])

#         else:
#             ax.legend([[label_formatter(sym)]
#                        for sym in numerical_data.columns],
#                       loc='lower right')

#         #plt.legend(loc='lower right')
#         plt.savefig(self.fig_name + '.png')
#         self.add_image(self.fig_name, width=NoEscape('15cm'))

#         if self.preview == True:
#             plt.show()

#         plt.close()


# class DataPlot(Figure):
#     _latex_name = 'figure'

#     def __init__(self,
#                  fig_name='Name',
#                  *args,
#                  preview=False,
#                  position=None,
#                  **kwargs):
#         super().__init__(position=position, **kwargs)

#         self.fig_name = str(fig_name)
#         #self._latex_name='figure' #super()._latex_name
#         self.preview = preview

#     def add_data_plot(self, *args, filename=None, width='15cm', **kwargs):

#         import matplotlib.pyplot as plt
#         if not filename:
#             current_time = dtime.datetime.now().timestamp()
#             filename = f'autoadded_figure_{current_time}.png'

#         plt.savefig(filename, *args, **kwargs)
#         self.add_image(filename, width=NoEscape(width))

#         if self.preview == True:
#             plt.show()


# class ReportSection(Section):

#     _latex_name = 'section'

#     def __init__(self,
#                  title,
#                  results_frame=None,
#                  analysis_key='',
#                  analysis_name='dynamics',
#                  results_col_name='simulations',
#                  preview=True):

#         self.analysis_key = analysis_key
#         self.analysis_name = analysis_name
#         self.sims_name = results_col_name
#         self.preview = preview
#         super().__init__(title)

#     def influeance_level(self, data):

#         sensitivity_level = {
#             0.02: 'znikomy wpływ',
#             0.1: 'mały wpływ',
#             0.5: 'średni wpływ',
#             0.7: 'znaczny wpływ',
#             0.9: 'istotny wpływ'
#         }

#         indicator = round(
#             2 * float(
#                 abs((data.abs().max() - data.abs().min()) /
#                     (0.001 + data.abs().max() + data.abs().min()))), 1) / 2

#         select = sensitivity_level

#         a = np.asarray(list(select.keys()))

#         #print('indicator: ', indicator)
#         #print(select[a.flat[np.abs(a - indicator).argmin()]])
#         #print(select[a.flat[np.abs(a - indicator).argmin()]])

#         return select[a.flat[np.abs(a - indicator).argmin()]]

#     def match_unit_to_data(self, parameter, value, units_dict=None):

#         if parameter in units_dict:
#             return value * units_dict[parameter]
#         else:
#             return value * pint.Quantity(1)

#     def get_feature_dict(self,
#                          numerical_data=None,
#                          given_data_dict=None,
#                          units_dict={},
#                          marker=None,
#                          string=None):

#         feature_dict = {
#             'param_name': vlatex(self.analysis_key),
#             'string': string,
#             'entries_no': 'nnn'
#         }
#         if marker:
#             feature_dict.update({'nr_rys': Ref(marker).dumps()})

#         if type(numerical_data) != type(None):

#             column_names = numerical_data.columns
#             feature_dict.update(
#                 {str(name): vlatex(name)
#                  for name in column_names})
#             feature_dict.update({
#                 str(name) + '_max':
#                 '{val:~Lx}'.format(val=self.match_unit_to_data(
#                     name,
#                     numerical_data[name].round(1).max(),
#                     units_dict=units_dict))
#                 for name in column_names
#             })
#             feature_dict.update({
#                 str(name) + '_min':
#                 '{val:~Lx}'.format(val=self.match_unit_to_data(
#                     name,
#                     numerical_data[name].round(1).min(),
#                     units_dict=units_dict))
#                 for name in column_names
#             })
#             feature_dict.update({
#                 str(name) + '_idmax':
#                 '{val:~Lx}'.format(val=self.match_unit_to_data(
#                     name,
#                     numerical_data[name].round(1).idxmax(),
#                     units_dict=units_dict))
#                 for name in column_names
#             })
#             feature_dict.update({
#                 str(name) + '_idmin':
#                 '{val:~Lx}'.format(val=self.match_unit_to_data(
#                     name,
#                     numerical_data[name].round(1).idxmin(),
#                     units_dict=units_dict))
#                 for name in column_names
#             })

#         if type(given_data_dict) != type(None):
#             feature_dict.update({
#                 'given_data':
#                 '\(' + '\), \('.join([
#                     vlatex(lhs) +
#                     '={rhs:~Lx.2f}'.format(rhs=self.match_unit_to_data(
#                         lhs, rhs, units_dict=units_dict))
#                     for lhs, rhs in given_data_dict.items()
#                 ]) + '\)'
#             })

#         return feature_dict

#     def add_introduction(
#         self,
#         title='Zakres prowadzonej analizy',
#         numerical_data=None,
#         units_dict={},
#         initial_description='Initial description',  #random.choice(sec_intro_composed_model)
#         ending_summary='Ending summary',
#         string=None,
#         # random.choice(sec_intro_composed_model)
#     ):

#         with self.create(Subsection(NoEscape(title))) as subsec:

#             tab_mrk = Marker('given_data_' + str(self.analysis_name) + '_' +
#                              str(self.analysis_key),
#                              prefix='tab')

#             symbols_df = pd.DataFrame(data=[
#                 {
#                     'Zmienna': '$ {eq}  $'.format(eq=vlatex(key)),
#                     'Liczba przypadków': len(column.unique()),
#                     'Wartość minimalna': column.round(2).min(),
#                     'Wartość maksymalna': column.round(2).max(),
#                     'Jednostka': '$ {unit:~Lx}  $'.format(
#                         unit=units_dict[key]),
#                 } for key, column in numerical_data.items()
#             ])

#             format_dict = (self.get_feature_dict(
#                 numerical_data=None,
#                 given_data_dict=None,
#                 units_dict=units_dict,
#                 marker=tab_mrk,
#                 string=string,
#             ))

#             subsec.append(
#                 NoEscape(str(initial_description).format(**format_dict)))

#             with subsec.create(DataTable(symbols_df, position='H')) as table:

#                 table.add_caption(
#                     NoEscape(
#                         'Zestawienie parametrów modelu'.format(**format_dict)))
#                 table.add_table(numerical_data.T)
#                 table.append(Label(tab_mrk))

#             subsec.append(NoEscape(str(ending_summary).format(**format_dict)))

#             #subsec.append(NewPage())
#             #subsec.append(NoEscape('\\'))
#             subsec.append(NoEscape('\par'))

#     def add_simulational_results(
#             self,
#             title='Wyniki symulacji numerycznych',
#             numerical_data=None,
#             units_dict={},
#             initial_description='Initial description',  #random.choice(intro_bank_composed_model)
#             ending_summary='Ending summary',  #random.choice(summary_bank_composed_model)
#             ylabel='',
#             xlabel='',
#             subplots=True,
#             grid=True,
#             num_yticks=10,
#             fontsize=None,
#             figsize=(10, 4),
#             caption='Przebiegi czasowe modelu dla rozważanego zakresu danych',
#             plots_no=1):

#         with self.create(Subsection(title)) as subsec:

#             simulation_results_frame = numerical_data

#             for key, row in simulation_results_frame.iterrows():

#                 data_with_units = {
#                     parameter: value
#                     for parameter, value in row.items()
#                     if isinstance(parameter, Symbol)
#                 }

#                 current_time = dtime.datetime.now().timestamp()
#                 current_fig_mrk = Marker(
#                     ('data_plot_' + str(self.analysis_key) + '_' +
#                      str(1 + next(plots_no_gen))),
#                     prefix='fig')

#                 format_dict = {
#                     **(self.get_feature_dict(numerical_data=row['simulations'],
#                                              given_data_dict=data_with_units,
#                                              units_dict=units_dict,
#                                              marker=current_fig_mrk)),
#                     #**{str(name):vlatex(name) for name in row['simulations'].columns}
#                 }

#                 subsec.append(
#                     NoEscape(str(initial_description).format(**format_dict)))

#                 #                 print(
#                 #                     np.array_split(range(len(row['simulations'].columns)),
#                 #                                    plots_no))

#                 for no, control_list in enumerate(
#                         np.array_split(range(len(row['simulations'].columns)),
#                                        plots_no)):

#                     ##print(row['simulations'].iloc[:,int(control_list[0]):int(control_list[-1]+1)])
#                     #print('control list', control_list)

#                     current_time = dtime.datetime.now().timestamp()
#                     current_fig_mrk = Marker(
#                         ('data_plot_' + str(self.analysis_key) + '_' +
#                          str(next(plots_no_gen))),
#                         prefix='fig')

#                     format_dict = {
#                         **(self.get_feature_dict(numerical_data=row['simulations'],
#                                                  given_data_dict=data_with_units,
#                                                  units_dict=units_dict,
#                                                  marker=current_fig_mrk)),
#                         #**{str(name):vlatex(name) for name in row['simulations'].columns}
#                     }

#                     with subsec.create(
#                             PlottedData(row['simulations'].
#                                         iloc[:,
#                                              int(control_list[0]
#                                                  ):int(control_list[-1] + 1)],
#                                         './plots/fig_' + str(current_time) +
#                                         '_' + str(no),
#                                         position='H',
#                                         units_dict=units_dict,
#                                         preview=self.preview)) as fig:
#                         fig.add_data_plot(None,
#                                           ylabel=ylabel,
#                                           subplots=subplots,
#                                           xlabel=xlabel,
#                                           grid=grid,
#                                           figsize=figsize,
#                                           num_yticks=num_yticks,
#                                           fontsize=fontsize)
#                         fig.add_caption(
#                             NoEscape(
#                                 #'Przebiegi czasowe modelu dla danych: {given_data}.'.format(**format_dict)
#                                 caption.format(**format_dict)))
#                         fig.append(Label(current_fig_mrk))

#                 subsec.append(
#                     NoEscape(str(ending_summary).format(**format_dict)))

#                 #subsec.append(NewPage())
#                 #subsec.append(NoEscape('\\'))
#                 subsec.append(NoEscape('\par'))

#     def prepare_summary_data(self, simulation_results_frame):
#         summary_frame = pd.DataFrame(data=[
#             a.max() for a in simulation_results_frame['simulations'].values
#         ],
#                                      index=[
#                                          str(elem)
#                                          for elem in simulation_results_frame[
#                                              self.analysis_key].round(2)
#                                      ])

#         return summary_frame

#     def add_summary(
#             self,
#             title='Analiza otrzymanych wyników',
#             numerical_data=None,
#             units_dict={},
#             xlabel=' ',
#             figsize=(10, 4),
#             initial_description='Initial description',
#             ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
#     ):
#         ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty symulacji numerycznych. Dla uzyskanych danych symulacyjnych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego układu. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''

#         summary_frame = numerical_data

#         with self.create(Subsection(title)) as subsec:

#             current_time = dtime.datetime.now().timestamp()
#             summary_mrk = Marker('summary_' + str(self.analysis_name) + '_' +
#                                  str(self.analysis_key),
#                                  prefix='fig')

#             #             summary_frame = pd.DataFrame(
#             #                 data=[
#             #                     a.max()
#             #                     for a in simulation_results_frame['simulations'].values
#             #                 ],
#             #                 index=[
#             #                     str(elem) for elem in simulation_results_frame[
#             #                         self.analysis_key].round(2)
#             #                 ])

#             format_dict = (self.get_feature_dict(numerical_data=summary_frame,
#                                                  given_data_dict=None,
#                                                  units_dict=units_dict,
#                                                  marker=summary_mrk))

#             #print(format_dict)

#             subsec.append(
#                 NoEscape(str(initial_description).format(**format_dict)))

#             step_key_influence_dict = {
#                 str(name) + '_inf': self.influeance_level(summary_frame[name])
#                 for name in summary_frame.columns
#             }

#             #print(step_key_influence_dict)

#             with subsec.create(
#                     PlottedData(summary_frame,
#                                 './plots/fig_summary_' + str(current_time),
#                                 position='H',
#                                 units_dict=units_dict,
#                                 preview=True)) as fig:
#                 fig.add_data_plot(summary_frame,
#                                   xlabel=xlabel,
#                                   figsize=figsize)
#                 fig.add_caption(
#                     NoEscape(
#                         'Zestawienie wyników przeprowadzonej analizy.'.format(
#                             **format_dict)))
#                 fig.append(Label(summary_mrk))

#             subsec.append(
#                 NoEscape(
#                     str(ending_summary).format(
#                         **{
#                             **format_dict,
#                             **step_key_influence_dict,
#                             #**{str(name):vlatex(name) for name in row['simulations'].columns}
#                         })))

#             #subsec.append(NewPage())
#             #subsec.append(NoEscape('\\'))
#             subsec.append(NoEscape('\par'))

#     def add_optimization_result(
#             self,
#             title='Wyniki optymalizacji jednokryterialnej',
#             numerical_data=None,
#             units_dict={},
#             xlabel=' ',
#             figsize=(10, 4),
#             initial_description='Initial description',
#             ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
#     ):

#         self.add_summary(
#             title=title,
#             numerical_data=numerical_data,
#             units_dict=units_dict,
#             initial_description=initial_description,
#             xlabel=xlabel,
#             figsize=figsize,
#             ending_summary=
#             'Ending summary',  # random.choice(conclusion_bank_composed_model)
#         )

#     def add_measurement_summary(
#             self,
#             title='Analiza otrzymanych wyników pomiarów',
#             numerical_data=None,
#             units_dict={},
#             initial_description='Initial description',
#             ending_summary='Ending summary',  # random.choice(conclusion_bank_composed_model)
#     ):
#         ''' Dla rozważanego modelu dynamicznego wózka inwalidzkiego wraz z napędem RC przedstawiono efekty pomiarów. Dla uzyskanych danych liczbowych, przygotowano wykresy przedstawiające maksymalne wartości osiąganych amplitud w funkcji analizowanego parametru dla współrzędnych uogólnionych modelu oraz ich pierwszych pochodnych (przemieszczeń i prędkości). Opracowane wykresy porównawcze pozwoliły na określenie wpływu badanych parametrów na dynamikę rozpatrywanego wózka. Bazując na wynikach przerprowadzonych symulacji przygotowano zestawienie dla parametru \({param_name}\).  '''

#         simulation_results_frame = numerical_data
#         with self.create(Subsection(title)) as subsec:

#             current_time = dtime.datetime.now().timestamp()
#             summary_mrk = Marker('measurment_summary_' +
#                                  str(self.analysis_name) +
#                                  str(self.analysis_key),
#                                  prefix='fig')

#             summary_frame = pd.DataFrame(
#                 data=[
#                     a.max()
#                     for a in simulation_results_frame['simulations'].values
#                 ],
#                 index=[
#                     str(elem) for elem in simulation_results_frame[
#                         self.analysis_key].round(2)
#                 ])

#             format_dict = (self.get_feature_dict(numerical_data=None,
#                                                  given_data_dict=None,
#                                                  units_dict=units_dict,
#                                                  marker=summary_mrk))

#             subsec.append(
#                 NoEscape(str(initial_description).format(**format_dict)))

#             format_dict = (self.get_feature_dict(numerical_data=summary_frame,
#                                                  given_data_dict=None,
#                                                  units_dict=units_dict,
#                                                  marker=summary_mrk))

#             step_key_influence_dict = {
#                 str(name) + '_inf': self.influeance_level(summary_frame[name])
#                 for name in summary_frame.columns
#             }

#             with subsec.create(
#                     PlottedData(summary_frame,
#                                 './plots/fig_summary_' + str(current_time),
#                                 position='H',
#                                 preview=True)) as fig:
#                 fig.add_data_plot(summary_frame)
#                 fig.add_caption(
#                     NoEscape(
#                         'Wpływ zmiany rozważanego parametru na dynamikę układu.'
#                         .format(**format_dict)))
#                 fig.append(Label(summary_mrk))

#             subsec.append(
#                 NoEscape(
#                     str(ending_summary).format(
#                         **{
#                             **format_dict,
#                             **step_key_influence_dict,
#                             #**{str(name):vlatex(name) for name in row['simulations'].columns}
#                         })))

#             subsec.append(NewPage())
