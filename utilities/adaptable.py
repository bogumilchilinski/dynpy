import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Latex
from IPython.display import Markdown as IPMarkdown
from IPython.display import display
from numpy import fft
from pylatex import (
    Alignat,
    Axis,
    Command,
    Document,
    Figure,
    Label,
    Marker,
    Math,
    Package,
    Plot,
    Quantity,
    Ref,
    Section,
    Subsection,
    Table,
    Tabular,
    TikZ,
    TikZCoordinate,
)
from pylatex.base_classes import Environment, Options
from pylatex.utils import NoEscape, italic
from sympy import (
    Eq,
    Expr,
    Float,
    Function,
    ImmutableMatrix,
    Matrix,
    Number,
    Symbol,
    lambdify,
    latex,
    symbols,
)
from sympy.core.relational import Relational
from sympy.physics import units
from sympy.physics.mechanics import vlatex

# from .documents import tikz


# from sympy.physics import units
# ureg = units

ureg = units

import os


class ReportCache:

    _file_names = {}

    def __init__(self, file_names):
        self._file_names = file_names

    @classmethod
    def update_existing_files(cls, directory):
        file_list = []
        for file_name in os.listdir(directory):
            file_path = os.path.join(directory, file_name)

            if os.path.isfile(file_path) and ".tex" in file_path:
                with open(file_path, "r") as file:
                    data = file.read()
                    cls._file_names[data] = file_path.replace(".tex", "")

        return list(cls._file_names.values())


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

    _cls_container = Section("Dummy")
    cls_container = Section("Dummy")
    _container = Section("Dummy")
    cls_path = "."
    _caption = "Figure describes the numerical data"
    _label = "fig:markerIsMissing"
    _units = {}
    _autoreport = False
    # _frame = TimeDataFrame()
    _frame = None
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


def plots_no():
    num = 0
    while True:
        yield num
        num += 1


default_colors = ["red", "blue", "orange", "teal", "black", "green"]


class BaseIndexTransformer:
    """_summary_

    Returns
    -------
    _type_
        _description_

    maintenance by BCh
    """

    _coord_level = -1
    _name_level = "auto"

    def __init__(self, data):

        self._data = data

    @property
    def index(self):
        return self._data.index

    @property
    def columns(self):
        return self._data.columns

    def _spot_axis_label(self):

        cols = self.columns
        lvl_to_drop = self._coord_level

        if isinstance(cols, pd.MultiIndex):
            ylabel = list(cols.get_level_values(lvl_to_drop).unique())

        else:
            ylabel = [cols.name]

        return ylabel

    def _spot_labels(self):

        cols = self.columns
        lvl_to_drop = self._coord_level

        # print('raw input',cols)

        if isinstance(cols, pd.MultiIndex):

            if len(self._spot_axis_label()) == 1:
                cols = cols.droplevel(lvl_to_drop)

        if isinstance(cols, pd.MultiIndex):
            labels = list(cols.to_flat_index())

        #             #print('len',len(labels),labels,type(labels))

        #             if not isinstance(labels[0],tuple):
        #                 labels =[(label,)  for label in labels]

        else:

            labels = [(entry,) for entry in cols]

        # print('len',len(labels),labels,type(labels))
        return labels

    def get_columns_str(self):

        cols = self.columns
        lvl_to_drop = self._coord_level

        labels = self._spot_labels()

        return [", ".join([str(entry) for entry in label]) for label in labels]

    def get_axis_str(self):

        return ", ".join([str(entry) for entry in self._spot_axis_label()])


class DataAxis(Axis, ReportModule):
    """
    Args
    ----
    options: str, list or `~.Options`
    Options to format the axis environment.

    maintenance by BCh
    """

    _latex_name = "axis"
    _handle = "plotEntry"
    _at = None

    _x_axis_name = "x"
    _y_axis_name = "y"
    _legend_fontsize = r"\small "
    _default_colours = default_colors
    _default_columns_no = 2

    _height = NoEscape(r"5cm")
    _width = "0.9\\textwidth"
    _x_axis_empty = False
    _y_limit = None
    _line_style = None
    _legend_pos = "north east"
    _default_transformer = BaseIndexTransformer
    _reversedx = None
    _manual_at = None

    def __init__(
        self,
        plotdata=None,
        height=None,
        width=None,
        options=None,
        *,
        y_limit=None,
        x_axis_empty=None,
        colour=None,
        at=None,
        reversedx=None,
        handle=None,
        data=None,
    ):
        """
        Args
        ----
        options: str, list or `~.Options`
            Options to format the axis environment.
        """
        self._plotdata = plotdata

        if handle is not None:
            self._handle = handle

        if at is not None:
            self._at = at

        self._colour = colour

        if height is not None:
            self._height = height

        if width is not None:
            self._width = width

        if x_axis_empty is not None:
            self._x_axis_empty = x_axis_empty

        if y_limit is not None:
            self._y_limit = y_limit

        if reversedx is not None:
            self._reversedx = reversedx
            display("Reversed x is " + f"{self._reversedx}")

        if options is None:
            options = self._axis_options

        super().__init__(options=options, data=data)

        if plotdata is not None:
            for data in self._plots:
                self.append(data)

    @classmethod
    def set_at(self, at_val):
        self._manual_at = at_val

    @classmethod
    def set_default_colours(cls, colour=None):
        cls._default_colours = colour
        TikZPlot._default_colours = colour

        return cls

    @classmethod
    def set_width(cls, width=None):
        cls._width = width
        return cls

    @classmethod
    def set_height(cls, height=None):
        cls._height = height
        return cls

    @property
    def handle(self):
        return self._handle

    @property
    def at(self):

        if self._at is not None:
            return self._at
        else:
            return None  # return '(0,0)'

    @property
    def _axis_options(self):

        kwargs = {
            "height": NoEscape(self._height),
            "anchor": "north west",
            "ymajorgrids": "true",
            "xmajorgrids": "true",
            #'ylabel': self.y_axis_name,
        }

        if self.handle:
            kwargs["name"] = self.handle

        at_option = []

        if self.at is not None:
            at_option = [NoEscape(f"at={self.at}")]

        if self._manual_at is not None:
            display("MANUAL")
            display(self._manual_at)
            at_option = [
                NoEscape("legend style={font=\small,at=" + self._manual_at + "}")
            ]
        #         else:
        #             at_option = ['']
        ax_options = []

        if self._x_axis_empty:
            ax_options.append(
                NoEscape("xticklabels=\empty"),
            )
        else:
            ax_options.append(
                NoEscape(f"xlabel= {self.x_axis_name}"),
            )

        if self.__class__._y_limit is not None:
            ay_options = [
                NoEscape(f"ymin={(self.__class__._y_limit)[0]}"),
                NoEscape(f"ymax={(self.__class__._y_limit)[1]}"),
            ]
        else:
            ay_options = []

        if self._reversedx == True:
            ax_options.append("x dir=reverse")

        base_options = (
            [
                "grid style=dashed",
                f"legend pos={self._legend_pos}",
                NoEscape("legend style={font=\small}"),
                NoEscape(f"legend columns={self._default_columns_no}"),
                NoEscape("/pgf/number format/.cd,1000 sep={}"),
                # NoEscape('at={(0,0)}'),
                NoEscape(f"width={self._width}"),
                NoEscape(f"xmin={min(self.plotdata.index)}"),
                NoEscape(f"xmax={max(self.plotdata.index)}"),
                NoEscape(f"ylabel= {self.y_axis_name}"),
            ]
            + at_option
            + ax_options
            + ay_options
        )

        return Options(*base_options, **kwargs)

    @property
    def plotdata(self):
        """
        Returns data which is used for plot
        """

        return self._plotdata

    @property
    def height(self):
        """
        Returns dataframe index as numpy array related to data for plot
        """
        return self._height

    @property
    def width(self):
        """
        Returns dataframe index as numpy array related to data for plot
        """
        return self._width

    @property
    def _index(self):
        """
        Returns dataframe index as numpy array related to data for plot
        """
        return self._plotdata.index.to_numpy()

    @property
    def _index_limits(self):
        """
        Returns minimum and maximum index limits
        """
        if isinstance(self.plotdata, (pd.DataFrame, pd.Series)):
            xmin = min(self.plotdata.index)
            xmax = max(self.plotdata.index)
        else:
            xmin = None
            xmax = None

        return xmin, xmax

    @property
    def _values_limits(self):
        """
        Returns minimum and maximum values limits
        """
        if isinstance(self.plotdata, (pd.DataFrame, pd.Series)):
            ymin = -0.25  # min(self.plotdata.index)
            ymax = 0.25  # max(self.plotdata.index)
        else:
            ymin = None
            ymax = None

        return ymin, ymax

    @property
    def x_axis_name(self):
        """
        Returns name of plot "x" axis. If name is not given, returns default axis name.
        """
        if isinstance(self.plotdata, (pd.DataFrame, pd.Series)):
            name = self.plotdata.index.name
        else:
            name = None

        return "{" + name + "}" if name is not None else "{" + self._x_axis_name + "}"

    @property
    def _data_description(self):
        return self._default_transformer(self.plotdata)

    @property
    def y_axis_name(self):
        """
        Returns y axis name of the plot. If name is not defined, sets name of the axis to default "y"
        """

        if isinstance(self.plotdata, (pd.DataFrame)):
            cols_name = self.plotdata.columns.name

            if isinstance(self.plotdata.columns, pd.MultiIndex):

                name = self._data_description.get_axis_str()

            elif isinstance(
                self.plotdata.columns, pd.Index
            ):  # temporary solution for further improvement - must be rewritten
                name = ", ".join(self.plotdata.columns.unique())
            else:
                name = cols_name

        elif isinstance(self.plotdata, (pd.Series)):
            name = self.plotdata.name
        else:
            name = None

        return "{" + name + "}" if name is not None else "{" + self._y_axis_name + "}"

    def _preview(self):
        """
        Allows to create a preview of prepared plot
        """

        self.plotdata.plot(ylabel=self.y_axis_name)

    @property
    def _transformer(self):
        return self._default_transformer

    @property
    def labels(self):
        """
        Returns plot dataframe columns' labels.
        """

        trans_class = self._transformer

        if isinstance(self.plotdata, (pd.DataFrame)):

            labels = self._data_description.get_columns_str()

        elif isinstance(self.plotdata, (pd.Series)):
            labels = [self.plotdata.name]
        else:
            labels = ["Data"]

        return labels

    # Zaprogramować zabezpieczenie metody na data integrity!!!
    @property
    def _coordinates(self):
        """
        Returns data list which is prepared for plots.
        """
        data = self.plotdata
        if isinstance(data, pd.DataFrame):
            data_for_plot_list = [
                zip(self._index, plot_data) for label, plot_data in list(data.items())
            ]
        else:
            data_for_plot_list = [zip(self._index, data.to_numpy())]

        return data_for_plot_list  # [(x_elem,y_elem)  for x_elem,y_elem in data_for_plot_list if not np.isnan(y_elem) ]

    @property
    def legend_fontsize(self):
        """
        Returns size of the font used in plot's legend.
        """
        return self._legend_fontsize

    def _label_colour(self, no):
        """
        Returns colour of the label used in plot.
        """
        if self._colour:
            colour = self._colour
        else:
            indicator = no % len(self._default_colours)
            colour = self._default_colours[indicator]
        return colour

    @property
    def _plots(self):

        coords_pack = self._coordinates
        labels = self.labels
        plots_no = len(self._coordinates)

        def fix_nan(coords_list):

            return [
                (x_elem, y_elem)
                for x_elem, y_elem in coords_list
                if not np.isnan(y_elem)
            ]

        colours = [self._label_colour(no) for no, elem in enumerate(coords_pack)]

        if self._line_style is None:
            line_style_list = ["solid"] * plots_no

        elif type(self._line_style) == list:

            line_style_list = (self._line_style * plots_no)[0:plots_no]

        else:
            line_style_list = [self._line_style] * plots_no

        # return list(zip(coords_pack,labels,colours))
        if self.at is not None:
            plot_list = [
                Plot(
                    name=NoEscape(
                        NoEscape(self.legend_fontsize) + NoEscape(str(label))
                    ),
                    coordinates=fix_nan(list(coords)),
                    options=Options(f"{line_style}", color=colour),
                )
                for coords, label, colour, line_style in (
                    zip(coords_pack, labels, colours, line_style_list)
                )
            ]
        else:
            plot_list = [
                Plot(
                    name=NoEscape(
                        NoEscape(self.legend_fontsize) + NoEscape(str(label))
                    ),
                    coordinates=fix_nan(list(coords)),
                    options=Options(f"{line_style}", color=colour),
                )
                for coords, label, colour, line_style in (
                    zip(coords_pack, labels, colours, line_style_list)
                )
            ]

        return plot_list

    def _repr_markdown_(self):
        self._preview()
        return " "


class TikZPlot(TikZ, ReportModule):
    """
    Args
    ----
    options: str, list or `~.Options`
    Options to format the axis environment.

    maintenance by BCh
    """

    _latex_name = "tikzpicture"

    _x_axis_name = "x"
    _y_axis_name = "y"
    _legend_fontsize = r"\small "
    _default_colours = default_colors
    _default_columns_no = 1

    _subplots = False
    _height = 5
    _width = 0.9
    _ylim = None

    _subplots_gap = 0.0
    _subplots_horizontal_gap = 0.5
    _default_grid = None
    _grid = None

    _in_figure = False
    _figure_gen = lambda: Figure(position="htbp")
    _default_figure_env = None
    _image_parameters = {"width": None}

    _floats_no_gen = plots_no()
    _default_path = "./tikzplots"
    _filename = None
    _prefix = None
    _reverse_x = None

    _picture = True
    _caption = "Default caption"

    def __init__(
        self,
        plotdata=None,
        subplots=None,
        height=None,
        width=None,
        reverse_x=None,
        options=None,
        *,
        arguments=None,
        start_arguments=None,
        **kwargs,
    ):
        """
        Args
        ----
        options: str, list or `~.Options`
            Options to format the axis environment.
        """
        self._plotdata = plotdata
        self._subplots = subplots

        if options is None:
            options = self._axis_options

        if height is not None:
            self._height = height

        if width is not None:
            self._width = width

        if reverse_x is not None:
            self._reverse_x = reverse_x

        self._selected_colours = None

        super().__init__(
            options=options,
            arguments=arguments,
            start_arguments=start_arguments,
            **kwargs,
        )

        if plotdata is not None:
            if self.subplots:

                plots_no = len(plotdata.columns)
                cols, rows = self.grid
                empty_axis_list = self.grid_nodes[:-cols] + [False] * cols
                colours_list = self._default_colours * plots_no

                self._selected_colours = colours_list
                labels = self.ylabel_list

                for no, (column, pos, empty_axis) in enumerate(
                    zip(plotdata.columns, self.grid_nodes, empty_axis_list)
                ):

                    data = plotdata[[column]]
                    # data.columns.name = labels[no]
                    if reverse_x == True:
                        ax_data = self.axis_type(
                            data,
                            height=self.subplot_height,
                            width=self.subplot_width,
                            x_axis_empty=empty_axis,
                            colour=colours_list[no],
                            at=pos,
                            handle=f"subplot{no}",
                            reversedx=self._reverse_x,
                        )
                    else:
                        ax_data = self.axis_type(
                            data,
                            height=self.subplot_height,
                            width=self.subplot_width,
                            x_axis_empty=empty_axis,
                            colour=colours_list[no],
                            at=pos,
                            handle=f"subplot{no}",
                        )

                    self.append(ax_data)
            else:
                if reverse_x == True:
                    self.append(
                        self.axis_type(
                            plotdata,
                            height=self.height,
                            width=self.width,
                            reversedx=self._reverse_x,
                        )
                    )
                else:
                    self.append(
                        self.axis_type(plotdata, height=self.height, width=self.width)
                    )

    @classmethod
    def set_figure_environment(cls, env=None):

        from .report import Picture, StarredPicture

        if env == "figure*":
            # self._figure_gen = lambda: StarredPicture(caption=caption,width=NoEscape('\\textwidth'))
            cls._default_figure_env = StarredPicture
        elif env == "figure":
            # self._figure_gen = lambda: Figure(position='htbp')
            cls._default_figure_env = Picture

        return cls

    @property
    def _figure_env(self):
        from .report import Picture

        if self.__class__._default_figure_env is None:

            return Picture
        else:
            return self.__class__._default_figure_env

    @property
    def axis_type(self):
        return DataAxis

    @property
    def ylabel_list(self):

        plotdata = self._plotdata

        if plotdata is not None:

            if isinstance(plotdata.columns, pd.MultiIndex):
                return plotdata.columns.levels
            else:

                return plotdata.columns.unique().to_list()

        else:

            return []

    @property
    def subplots(self):
        if self._subplots is None:
            return self.__class__._subplots
        else:
            return self._subplots

    @property
    def _axis_options(self):

        return None

        if isinstance(self._ylim, list):
            return Options(
                "grid style=dashed",
                NoEscape("legend style={font=\small}"),
                NoEscape(f"legend columns={self._default_columns_no}"),
                NoEscape(f"width={self._width}"),
                NoEscape(f"xmin={min(self._plotdata.index)}"),
                NoEscape(f"xmax={max(self._plotdata.index)}"),
                NoEscape(f"ymin={min(self._ylim)}"),
                NoEscape(f"ymax={max(self._ylim)}"),
                height=NoEscape(self._height),
                anchor="north west",
                ymajorgrids="true",
                xmajorgrids="true",
                xlabel=self._x_axis_name,
                ylabel=self._y_axis_name,
            )
        elif self._ylim is None:
            return Options(
                "grid style=dashed",
                NoEscape("legend style={font=\small}"),
                NoEscape(f"legend columns={self._default_columns_no}"),
                NoEscape(f"width={self._width}"),
                NoEscape(f"xmin={min(self._plotdata.index)}"),
                NoEscape(f"xmax={max(self._plotdata.index)}"),
                height=NoEscape(self._height),
                anchor="north west",
                ymajorgrids="true",
                xmajorgrids="true",
                xlabel=self._x_axis_name,
                ylabel=self._y_axis_name,
            )

    @property
    def grid(self):

        if self._grid is None:
            grid = (1, len(self._plotdata.columns))
        else:
            grid = self._grid

        return grid

    #     @property
    #     def grid(self):
    #         if self._default_grid is None:
    #             return self.grid
    #         else:
    #             return self._default_grid

    #     @grid.setter
    #     def grid(self, size):
    #         self._grid_size=size

    @property
    def grid_nodes(self):

        cols_no, rows_no = self.grid

        coordinates = []

        for row_id in range(rows_no):
            y_coordinate = row_id * (self._height / rows_no)
            y_coord_round = round(y_coordinate, 2)

            for col_id in range(cols_no):
                if cols_no > 1:

                    x_coordinate = col_id * (self._width / cols_no)

                else:
                    x_coordinate = 0

                x_coord_round = round(x_coordinate, 2)
                position = "{" + f"({x_coord_round}cm , -{y_coord_round}cm)" + "}"
                coordinates.append(position)

        return coordinates

    @property
    def subplot_height(self):

        rows_no = self.grid[1]

        return NoEscape(f"{self._height/rows_no+1.5-self._subplots_gap}cm")

    @property
    def subplot_width(self):

        cols_no = self.grid[0]

        if cols_no > 1:

            return NoEscape(
                f"{self._width/cols_no-2-self._subplots_horizontal_gap-self._subplots_gap}cm"
            )

        else:
            return self.width

    @property
    def height(self):
        """
        Returns dataframe index as numpy array related to data for plot
        """
        if isinstance(self._height, float) or isinstance(self._height, int):
            return NoEscape(f"{self._height}cm")
        if isinstance(self._height, tuple):
            return NoEscape(f"{self._height[0]}" + self._height[1])
        if isinstance(self._height, str):
            return NoEscape(self._height)
        if isinstance(self._height, NoEscape):
            return self._height

    @property
    def width(self):
        """
        Returns dataframe index as numpy array related to data for plot
        """
        if isinstance(self._width, float) or isinstance(self._width, int):
            return f"{self._width}\\textwidth"
        if isinstance(self._width, tuple):
            return f"{self._width[0]}" + self._width[1]
        if isinstance(self._width, str):
            return self._width

    def in_figure(
        self,
        filename=None,
        position=None,
        caption=None,
        width=None,
        height=None,
        marker=None,
        **kwargs,
    ):

        from .documents import tikz

        ReportCache.update_existing_files(self.__class__._default_path)

        obj = copy.copy(self)
        obj._in_figure = True

        if caption is not None:
            obj._caption = caption

        standalone_plot = tikz.TikzStandalone()

        standalone_plot.append(self)

        filename = self.filename

        fig = self.__class__._figure_gen()
        fig.packages.append(Package("float"))

        img_params = self.__class__._image_parameters
        width = self.__class__._image_parameters["width"]

        if self._picture:
            from .report import Picture

            key = standalone_plot.dumps()

            if key in ReportCache._file_names:
                filename = ReportCache._file_names[key]

            else:

                standalone_plot.generate_pdf(
                    filename, clean_tex=False, compiler_args=["--lualatex"]
                )
                ReportCache._file_names[key] = filename

            PictureEnv = self._figure_env

            fig = PictureEnv(
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

    def _report(self):
        if self._in_figure is True:
            standalone_plot = tikz.TikzStandalone()
            standalone_plot.append(self)

            filename = self.filename

            fig = self.__class__._figure_gen()
            fig.packages.append(Package("float"))

            img_params = self.__class__._image_parameters

            if self._picture:
                from .report import Picture

                standalone_plot.generate_pdf(
                    filename, clean_tex=False, compiler_args=["--lualatex"]
                )

                width = self.__class__._image_parameters["width"]

                fig = self.__class__._figure_gen()
            else:
                standalone_plot.generate_tex(filename)

                fig.append(Command(command="input", arguments=filename))

            fig.add_caption(self._caption)

            self.cls_container.append(fig)
            display(fig)

        else:
            self.cls_container.append(self)
            self._plotdata.plot(
                subplots=self.subplots,
                # color=self._selected_colours
            )

    def _repr_markdown_(self):

        # self._plotdata.plot(subplots=self.subplots,color=self._selected_colours)
        self._report()

        return " "

    @classmethod
    def set_default_path(cls, path):
        cls._default_path = path


class DataMethods:

    _figure_gen = lambda: Figure(position="htbp")
    _image_parameters = {"width": NoEscape(r"0.9\textwidth")}
    _legend_fontsize = r" "
    _label_fontsize = r"\small "
    _template = Document(
        documentclass="standalone", geometry_options=None, document_options=["tikz"]
    )
    # _extra

    @classmethod
    def set_document_template(cls, template=None):
        if template is not None:
            cls._template = template
        return cls

    @classmethod
    def set_default_label_fontsize(cls, fontsize=None):
        if fontsize is not None:
            cls._label_fontsize = fontsize
        return cls

    @classmethod
    def set_default_legend_fontsize(cls, fontsize=None):
        if fontsize is not None:
            cls._legend_fontsize = fontsize
        return cls

    @classmethod
    def set_default_figure_generator(cls, figure_generator=None):
        if figure_generator is not None:
            cls._figure_gen = figure_generator
        return cls

    @classmethod
    def set_default_image_parameters(cls, image_parameters=None):
        if image_parameters is not None:
            cls._image_parameters = image_parameters
        return cls

    def _pylatex_tikz(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"7cm"),
        width=NoEscape("0.9\\textwidth"),
        x_axis_description=",xlabel={$t$},x unit=\\si{\\second},",
        y_axis_description=" ",
        subplots=False,
        extra_commands=None,
        options=None,
        smooth=False,
    ):

        if smooth:
            radius_str = NoEscape(f",rounded corners=1mm,")
        else:
            radius_str = ""

        labels_list = [label for label in self.columns]
        # labels_list=[(vlatex(label)) for label in self.columns]

        data_for_plot_list = [
            zip(self.index, plot_data) for label, plot_data in list(self.items())
        ]

        plots_no = len(data_for_plot_list)

        colors_multiplicator = np.ceil(plots_no / len(colors_list))

        # print(y_axis_description)
        # y_axis_description= ' '
        plot_options = (
            NoEscape(
                "anchor=north west,ymajorgrids=true,xmajorgrids=true,grid style=dashed,legend style={font=\\small},"
                + NoEscape(y_axis_description)
            )
            + NoEscape(", height=")
            + height
            + NoEscape(",width=")
            + width
            + NoEscape(f",xmin={min(self.index)},xmax={max(self.index)}")
        )

        print(plot_options)

        # with doc.create(Figure(position='!htb')) as fig:

        plot_options_list = [plot_options + NoEscape(",xticklabels=\empty,")] * (
            plots_no - 1
        )
        plot_options_list.append(plot_options + NoEscape(x_axis_description))

        tikzpicture = TikZ(options=options)

        legend_font_size = self.__class__._legend_fontsize

        if subplots == False:

            with tikzpicture.create(Axis(options=plot_options_list[-1])) as plot:

                for data_for_plot, label, color in zip(
                    data_for_plot_list,
                    labels_list,
                    colors_list * int(colors_multiplicator),
                ):
                    coordinates = data_for_plot

                    plot.append(
                        Plot(
                            name=NoEscape(NoEscape(legend_font_size + str(label))),
                            coordinates=coordinates,
                            options="color=" + color + ",solid" + radius_str,
                        )
                    )

        else:
            at_option = NoEscape("")
            for no, combined_plot_data in enumerate(
                zip(
                    data_for_plot_list,
                    labels_list,
                    colors_list * int(colors_multiplicator),
                )
            ):
                data_for_plot, label, color = combined_plot_data
                plot_name = NoEscape(",name=plot" + str(no) + ",")
                with tikzpicture.create(
                    Axis(options=plot_options_list[no] + plot_name + at_option)
                ) as plot:
                    coordinates = data_for_plot
                    plot.append(
                        Plot(
                            name=NoEscape(
                                NoEscape(legend_font_size) + NoEscape(str(label))
                            ),
                            coordinates=coordinates,
                            options="color=" + color + ",solid" + radius_str,
                        )
                    )
                    # at_option=NoEscape('at=(plot'+str(no)+'.below south west),')
                    at_option = NoEscape(
                        "at={(plot" + str(no) + ".south west)},yshift=-0.1cm,"
                    )

        if extra_commands is not None:
            tikzpicture.append(extra_commands)

        return tikzpicture

    def to_pylatex_plot(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"6cm"),
        width=NoEscape(r"0.49\textwidth"),
        x_axis_description=",xlabel={$t$},x unit=\si{\second},",
        y_axis_description="",
        subplots=False,
        extra_commands=None,
    ):

        tikz_pic = self._pylatex_tikz(
            filename,
            labels_list,
            colors_list,
            height,
            width,
            x_axis_description,
            y_axis_description,
            subplots,
            extra_commands=extra_commands,
        )

        return tikz_pic

    def to_standalone_plot(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"6cm"),
        width=NoEscape(r"0.49\textwidth"),
        x_axis_description=",xlabel={$t$},x unit=\si{\second},",
        y_axis_description="",
        subplots=False,
        legend_pos="north east",
        extra_commands=None,
        options=None,
        smooth=False,
        picture=True,
        template=None,
        *arg,
        **kwargs,
    ):

        geometry_options = {
            "margin": "0cm",
        }

        if template is None:

            doc = self.__class__._template

        doc = copy.deepcopy(self.__class__._template)
        # doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

        doc.packages.append(Package("siunitx"))
        doc.packages.append(Package("mathtools"))
        doc.packages.append(Package("float"))
        doc.packages.append(Package("tikz"))
        doc.packages.append(Package("pgfplots"))
        doc.append(Command("usepgfplotslibrary", arguments="units"))
        doc.append(Command("usetikzlibrary", arguments="spy"))

        doc.append(
            Command(
                "pgfplotsset",
                arguments=NoEscape(
                    r"compat=newest,label style={fontsize},legend pos=".format(
                        fontsize="{" + f"font={self.__class__._label_fontsize}" + "}"
                    )
                    + str(legend_pos)
                ),
            )
        )

        tikz_pic = self._pylatex_tikz(
            filename,
            labels_list,
            colors_list,
            height,
            width,
            x_axis_description,
            y_axis_description,
            subplots,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
        )

        doc.append(tikz_pic)

        if picture:
            doc.generate_pdf(filename, clean_tex=False, compiler_args=["--lualatex"])
        else:
            doc.generate_tex(filename)

        return doc.dumps()

    def to_tikz_plot(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"6cm"),
        width=NoEscape(r"0.49\textwidth"),
        x_axis_description=",xlabel={$t$},x unit=\si{\second},",
        y_axis_description="",
        subplots=False,
        legend_pos="north east",
        extra_commands=None,
        options=None,
        smooth=False,
        picture=True,
    ):

        return self.to_standalone_plot(
            filename,
            labels_list,
            colors_list,
            height,
            width,
            x_axis_description,
            y_axis_description,
            subplots,
            legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
            picture=picture,
        )

    def to_standalone_figure(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"6cm"),
        width=NoEscape(r"0.49\textwidth"),
        x_axis_description=",xlabel={$t$},x unit=\si{\second},",
        y_axis_description="",
        subplots=False,
        legend_pos="north east",
        extra_commands=None,
        options=None,
        smooth=False,
        picture=True,
    ):

        self.to_standalone_plot(
            filename,
            labels_list,
            colors_list,
            height,
            width,
            x_axis_description,
            y_axis_description,
            subplots,
            legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
            picture=picture,
        )

        fig = self.__class__._figure_gen()
        img_params = self.__class__._image_parameters

        if picture:

            width = self.__class__._image_parameters["width"]
            # width=NoEscape(r'0.8\textwidth')

            fig.add_image(filename, width=width)
        else:
            fig.append(Command(command="input", arguments=filename))

        return fig

    def to_pylatex_tikz(
        self, height=None, width=None, subplots=None, reverse_x=None, options=None
    ):
        return TikZPlot(
            LatexDataFrame.formatted(self),
            height=height,
            width=width,
            subplots=subplots,
            reverse_x=reverse_x,
            options=options,
        )

    def to_latex_dataframe(self):
        return LatexDataFrame.formatted(self)


class SpectralMethods(DataMethods):

    def is_uniformly_distributed(self):

        sample_length = max(self.index) - min(self.index) / len(self.index) - 1
        step_list = [
            self.index[i + 1] - self.index[i] for i in range(len(self.index) - 1)
        ]

        if all(np.round(element - step_list[0], 10) == 0 for element in step_list):
            return True
        else:
            return False

    #     def spectrum(self):

    #     spectrum={name:fft.fftshift(fft.fft(data)) for name,data in self.items()}
    #     f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

    #     return SpectrumFrame(data=spectrum,index=fft.fftshift(f_span))

    def to_time_domain(self):

        sampling_rate = self.index[1] - self.index[0]
        t_span = fft.fftshift(fft.fftfreq(len(self.index), d=sampling_rate))
        t_span = t_span[-1] + t_span

        timeseries = {
            # name: data._sort_for_ifft()
            name: fft.ifft(data._sort_for_ifft())
            for name, data in self.items()
        }

        return TimeDataFrame(data=(timeseries), index=t_span)

    def double_sided_rms(self):

        f_span_shifted_ds = fft.fftshift(self.index)
        spectrum_shifted_ds = fft.fftshift(abs(self) / len(self))

        return SpectrumSeries(
            data=spectrum_shifted_ds, index=f_span_shifted_ds, name=self.name
        )


#     def double_sided_spec(self):

#         f_span_shifted_ds=fft.fftshift(self.index)
#         spectrum_shifted_ds={name:fft.fftshift(abs(data)/len(data)) for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

#     def single_sided_spec(self):

#         f_span_shifted_ss=np.positive(fft.fftshift(self.index))
#         spectrum_shifted_ss={name:fft.fftshift(abs(data)/len(data))*np.heaviside([self.index],1)*2 for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)


class EntryWithUnit:
    _units = {}
    _latex_backend = lambda obj: obj if isinstance(obj, str) else vlatex(obj)

    _left_par = "["
    _right_par = "]"
    _unit_separator = " "

    @classmethod
    def set_default_surrounding(cls, left_or_both, right=None):

        cls._left_par = left_or_both[0]
        cls._right_par = left_or_both[1]

        return cls

    @classmethod
    def set_default_units(cls, units={}):

        cls._units = units
        return cls

    def __new__(cls, obj, units=None, latex_backend=None, **kwargs):
        obj_with_unit = super().__new__(cls)
        obj_with_unit._obj = obj

        if units is not None:
            obj_with_unit._units = units
        else:
            obj_with_unit._units = cls._units

        if isinstance(obj, Relational):
            obj_with_unit._left_par = ""
            obj_with_unit._right_par = ""
            obj_with_unit._quantity = obj_with_unit._obj.lhs
        else:
            obj_with_unit._quantity = obj_with_unit._obj

        has_unit = obj_with_unit._set_quantity_unit()

        if has_unit:

            obj_with_unit._obj = obj

            if latex_backend is not None:
                obj_with_unit._latex_backend = latex_backend
            else:
                obj_with_unit._latex_backend = cls._latex_backend

            return obj_with_unit

        else:

            return obj

    def _set_quantity_unit(self):

        if self._quantity in self._units:

            self._unit = self._units[self._quantity]
            return True
        else:

            self._unit = None

            return False

    def __str__(self):
        entry_str = self._obj.__str__()
        unit = self._unit
        left_par = self._left_par
        right_par = self._right_par

        sep = self._unit_separator

        if unit:
            return f"{entry_str}{sep}{left_par}{unit.__str__()}{right_par}"
        else:
            return f"{entry_str}"

    def __repr__(self):
        entry_str = self._obj.__repr__()
        unit = self._unit
        left_par = self._left_par
        right_par = self._right_par

        sep = self._unit_separator

        if unit:
            return f"{entry_str}{sep}{left_par}{unit.__repr__()}{right_par}"
        else:
            return f"{entry_str}"

    def _latex(self, *args):

        # print('abc')

        entry_str = self._latex_backend(self._obj)
        unit = self._unit
        left_par = self._left_par
        right_par = self._right_par

        sep = self._unit_separator
        # print('abc')
        # display(entry_str)
        # print(self._obj)
        # display(unit)

        if unit:
            if isinstance(unit, str):
                return f"{entry_str}{sep}{left_par}{unit}{right_par}"
            elif unit is None:
                return f"{entry_str}{sep}{left_par}-{right_par}"
            else:
                return f"{entry_str}{sep}{left_par}{latex(unit)}{right_par}"
        else:
            return f"{self._obj}"


class DataTable(Table, ReportModule):

    _position = "H"
    _latex_name = "table"
    packages = [
        Package("booktabs"),
        Package("longtable"),
        Package("caption"),
        Package("float"),
        Command("captionsetup[table]{skip=5pt}"),
    ]

    def __init__(self, numerical_data, position=None):

        if position is not None:
            self._position = position

        super().__init__(position=self._position)
        ##print(numerical_data)
        self._numerical_data = numerical_data
        self.position = position

    def add_table(
        self,
        numerical_data=None,
        index=False,
        longtable=False,
        multirow=True,
        column_format=None,
    ):
        self.append(NoEscape("\\centering"))
        self.append(NoEscape("%%%%%%%%%%%%%% Table %%%%%%%%%%%%%%%"))
        #         if numerical_data!=None:
        #             self._numerical_data=numerical_data

        tab = self._numerical_data

        #         latex_code=tab.style.to_latex(#index=index,
        #                                    #escape=False,
        #                              #longtable=longtable,
        # #                    multirow=multirow,
        #                     hrules=True,
        #                     column_format=column_format).replace(
        #                                  '\\toprule',
        #                                  '\\toprule \n \\midrule').replace(
        #                                      '\\bottomrule',
        #                                      '\\midrule \n \\bottomrule')
        latex_code = (
            tab.to_latex(
                index=index,
                # escape=False,
                # longtable=longtable,
                #                    multirow=multirow,
                # hrules=True,
                column_format=column_format,
            )
            .replace("\\toprule", "\\toprule \n \\midrule")
            .replace("\\bottomrule", "\\midrule \n \\bottomrule")
        )

        self.append(NoEscape(latex_code))

    @classmethod
    def set_default_position(cls, position=""):
        cls._position = position

    @classmethod
    def set_caption_space(cls, space):
        cls.packages.remove(Command("captionsetup[table]{skip=5pt}"))
        cls.packages.append(Command("captionsetup[table]{skip=" + str(space) + "pt}"))

        return cls


class BPASTSDataTable(DataTable):

    def add_table(
        self,
        numerical_data=None,
        index=False,
        longtable=False,
        multirow=True,
        column_format=None,
    ):
        self.append(NoEscape("\\centering"))
        self.append(NoEscape("%%%%%%%%%%%%%% Table %%%%%%%%%%%%%%%"))
        #         if numerical_data!=None:
        #             self._numerical_data=numerical_data

        tab = self._numerical_data

        #         latex_code=tab.style.to_latex(#index=index,
        #                                    #escape=False,
        #                              #longtable=longtable,
        # #                    multirow=multirow,
        #                     hrules=True,
        #                     column_format=column_format).replace(
        #                                  '\\toprule',
        #                                  '\\toprule \n \\midrule').replace(
        #                                      '\\bottomrule',
        #                                      '\\midrule \n \\bottomrule')

        format_string = "|l"
        col_format = len(tab.columns) * format_string + "|"
        if index is not False:
            col_format = col_format + "l|"
        latex_code = (
            tab.to_latex(
                index=index,
                # escape=False,
                # longtable=longtable,
                #                    multirow=multirow,
                # hrules=True,
                column_format=col_format,
            )
            .replace("\\\\", "\\\\ \\hline")
            .replace("\\toprule", "\\hline")
            .replace("\\midrule", "\\hline")
            .replace("\\bottomrule", "")
        )

        self.append(NoEscape(latex_code))


class MarkerRegistry(dict):

    _prefix = "automrk"
    _markers_dict = {}

    def __init__(self, prefix="automrk", sufix=None):
        super().__init__()
        self._prefix = prefix


class AutoMarker:
    _markers_dict = {}
    _prefix = "eq"
    _name = None
    _floats_no_gen = plots_no()
    _sufix = ""

    @classmethod
    def set_default_sufix(cls, sufix=""):
        cls._sufix = sufix

    @classmethod
    def add_marker(cls, elem, marker):

        from .report import Picture

        if isinstance(elem, (AdaptableDataFrame)):
            elem_id = elem._get_str_key()
            if elem._subplot is None:
                prefix = "tab"
            else:
                prefix = "fig"

        elif isinstance(elem, Picture):
            elem_id = elem._get_str_key()
            prefix = "fig"

        elif isinstance(elem, (pd.DataFrame, pd.Series)):
            elem_id = elem.style.to_latex()
            prefix = "fig"

        elif isinstance(elem, Matrix):
            elem_id = ImmutableMatrix(elem)
            prefix = "eq"
        else:
            elem_id = elem

        cls._markers_dict[elem_id] = marker

    #         return elem_id, prefix

    def __init__(self, elem, prefix="eq", name=None, sufix=None):

        from .report import Picture

        if isinstance(elem, (AdaptableDataFrame)):
            elem_id = elem._get_str_key()
            if elem._subplot is None:
                prefix = "tab"
            else:
                prefix = "fig"

        elif isinstance(elem, Picture):
            elem_id = elem._get_str_key()
            prefix = "fig"

        elif isinstance(elem, (pd.DataFrame, pd.Series)):
            elem_id = elem.style.to_latex()
            prefix = "fig"

        elif isinstance(elem, Matrix):
            elem_id = ImmutableMatrix(elem)
            prefix = "eq"
        else:
            elem_id = elem

        self._marker_name = elem.__class__.__name__
        self._elem_id = elem_id
        self._prefix = prefix
        self._sufix = sufix

        self._get_marker()

    def _get_marker(self, elem=None):

        if elem is None:
            elem = self._elem_id

        available_markers = self._markers_dict

        if elem in available_markers:
            marker = available_markers[elem]
        else:

            marker = Marker(
                f"Mrk{self._marker_name}{next(self._floats_no_gen)}{self.__class__._sufix}",
                prefix=self._prefix,
            )
            self._markers_dict[elem] = marker

        self._marker = marker
        return marker

    @property
    def marker(self):
        return self._marker

    def __repr__(self):
        return f"AutoMarker for {self.marker}"

    def __str__(self):
        return Ref(self.marker).dumps()


class BasicFormattingTools(DataMethods):

    _floats_no_gen = plots_no()

    _latex_backend = lambda obj: obj if isinstance(obj, str) else vlatex(obj)
    _label_formatter = lambda obj: (
        f"${_latex_backend(obj)}$"
        if isinstance(obj, (Expr, Eq, EntryWithUnit))
        else obj
    )
    _unit_selector = EntryWithUnit
    _domain = None
    _units = {}
    _applying_func = lambda x: x
    _init_ops = True

    _default_sep = ", "
    _container = []
    _default_path = "./tikzplots"
    _picture = True

    _subplot = False
    _caption = "Default caption"

    _default_width = NoEscape(r"0.9\textwidth")
    _default_height = NoEscape(r"6cm")
    _preview_mode = False

    @classmethod
    def set_default_width(cls, width=NoEscape(r"0.9\textwidth")):
        cls._default_width = width
        return cls

    @classmethod
    def set_preview_mode(cls, preview=False):
        cls._preview_mode = preview
        return cls

    @classmethod
    def set_default_height(cls, height=NoEscape(r"6cm")):
        cls._default_height = height

        return cls

    @classmethod
    def set_default_column_separator(cls, sep=", "):
        cls._default_sep = sep

        return cls

    @classmethod
    def set_picture_mode(cls, picture=False):

        cls._picture = picture
        return cls

    @classmethod
    def set_directory(cls, path="./"):

        cls._default_path = path
        return cls

    @classmethod
    def set_default_container(cls, container=[]):
        cls._container = container

        return cls

    @classmethod
    def set_container(cls, container=[]):
        cls._container = container

        return cls

    @classmethod
    def set_default_units(cls, units={}):

        cls._units = units
        return cls

    @classmethod
    def set_default_unit_selector(cls, selector=EntryWithUnit):

        cls._unit_selector = selector
        return cls

    @staticmethod
    def _label_formatter(obj):

        # print('type',type(obj))

        # print(isinstance(obj, (Expr, Eq, EntryWithUnit)))

        latex_backend = BasicFormattingTools._latex_backend

        if isinstance(obj, (Symbol, Function, Expr, Eq, EntryWithUnit)):
            formatted_obj = f"${latex_backend(obj)}$"
        else:
            formatted_obj = obj
        # print('effect',formatted_str)

        return formatted_obj

    def applying_method(self, data, func=None, **kwargs):
        if func:
            # print('func is used')
            ops_func = func

        elif self.__class__._applying_func is not None:
            # print('class func is used')

            ops_func = self.__class__._applying_func
        else:
            # print('identity is used')
            ops_func = lambda data: data

        return ops_func(data)

    def _modify_axis(self, func, axis=0):

        new_obj = self.copy()
        new_obj_idx = new_obj.axes[axis]
        # idx_frame = new_obj_idx.to_frame().map(func)

        # print('idx',new_obj_idx)
        # display('map',idx_frame)

        # if isinstance(new_obj_idx, pd.MultiIndex):
        #    new_obj_idx = pd.MultiIndex.from_frame(idx_frame)
        #    #new_obj_idx.names=map(func,new_obj_idx.names)
        # else:
        #    #new_obj_idx = pd.Index((idx_frame),name=new_obj_idx.name)
        #    new_obj_idx = new_obj_idx.map(func)

        return new_obj.rename(func, axis=axis)

    def _modify_axis_name(self, func, axis=0):
        new_obj = self.copy()
        new_obj_idx = new_obj.axes[axis]

        if isinstance(new_obj_idx, pd.MultiIndex):

            new_obj_idx.names = map(func, new_obj_idx.names)

        else:

            # new_obj_idx = pd.Index((idx_frame),name=new_obj_idx.name)
            new_obj_idx = new_obj_idx.map(func)
            new_obj_idx.name = func(new_obj_idx.name)

        return new_obj.set_axis(new_obj_idx, axis=axis)

    def set_multiindex_axis(self, axis=0):

        if axis == "index":
            axis = 0
        elif axis == "colums":
            axis = 1

        idx = self.axes[axis]

        if isinstance(idx, pd.MultiIndex):

            new_obj = self.copy()

        else:
            midx = pd.MultiIndex.from_tuples(idx)
            new_obj = self.copy().set_axis(midx, axis=axis)

        return new_obj

    def set_multiindex_columns(self):

        return self.set_multiindex_axis(axis=1)

    def set_flat_index_axis(self, axis=0):

        if axis == "index":
            axis = 0
        elif axis == "colums":
            axis = 1

        idx = self.axes[axis]

        if isinstance(idx, pd.MultiIndex):
            midx = idx.to_flat_index()
            new_obj = self.copy().set_axis(midx, axis=axis)

        else:

            new_obj = self.copy()

        return new_obj

    def set_flat_index_columns(self):

        return self.set_flat_index_axis(axis=1)

    def switch_axis_type(self, axis=0):

        if axis == "index":
            axis = 0
        elif axis == "colums":
            axis = 1

        idx = self.axes[axis]

        if isinstance(idx, pd.MultiIndex):

            new_obj = self.set_flat_index_axis(axis=axis)

        else:
            new_obj = self.copy().set_multiindex_axis(axis=axis)

        return new_obj

    def switch_index_type(self):

        return self.switch_axis_type(axis=0)

    def format_axis_names(self, formatter=None, axis=0):
        if formatter is None:

            # formatter = self.__class__._label_formatter
            formatter = self._label_formatter

        new_frame = self._modify_axis(formatter, axis=axis)._modify_axis_name(
            formatter, axis=axis
        )

        return new_frame

    def format_index_names(self, formatter=None):

        return self.format_axis_names(formatter=formatter, axis=0)

    def format_columns_names(self, formatter=None):

        return self.format_axis_names(formatter=formatter, axis=1)

    def format_axes_names(self, formatter=None):

        return self.format_axis_names(formatter=formatter, axis=1).format_axis_names(
            formatter=formatter, axis=0
        )

    def set_str_index_axis(self, axis=0, separator=", "):

        new_obj = self.copy()
        new_idx = new_obj.axes[axis]

        if isinstance(new_idx, pd.MultiIndex):
            new_obj = self.copy().set_flat_index_axis(axis=axis)

            new_obj = new_obj.set_axis(
                [separator.join(map(str, entry)) for entry in new_idx], axis=axis
            )

        return new_obj

    def set_str_index_columns(self, separator=", "):

        return self.set_str_index_axis(axis=1, separator=separator)

    def to_eq_axis_form(self, axis=0):
        new_obj = self.copy()

        return new_obj

    def to_named_axis_form(self, axis=0):
        new_obj = self.copy()

        idx = new_obj.axes[axis]

        if isinstance(idx, pd.MultiIndex):
            # print('MultiIndex modification has been not supported yet')
            numbered_obj = new_obj.format_axis_names(
                lambda name: float(name.rhs) if isinstance(name, Eq) else name,
                axis=axis,
            )
            # sym_wyn.loc[:,ix[:,:,dyn_sys.phi_1]].plot()

            new_obj_with_name = numbered_obj.format_columns_names(
                lambda name: float(name.rhs) if isinstance(name, Eq) else name
            )
            new_obj_with_name.axes[axis].names = [
                level_idx[0].lhs if isinstance(level_idx[0], Eq) else level_idx.name
                for level_idx in idx.levels
            ]

        else:
            name = idx.name

            if all([isinstance(entry, Relational) for entry in idx]):
                # TODO: add numerical types recognizing
                new_obj_with_name = new_obj._modify_axis(
                    lambda elem: (
                        float(elem.rhs) if isinstance(elem.rhs, Number) else elem.rhs
                    ),
                    axis=axis,
                )

                # print(new_obj_with_name.axes[axis].name)
                new_obj_with_name.axes[axis].name = idx[0].lhs

        return new_obj_with_name  # .rename(Symbol('l_w'),axis=axis)

    def to_named_index_form(self):
        return self.to_named_axis_form(axis=0)

    def to_named_columns_form(self):
        return self.to_named_axis_form(axis=1)

    def fit_units_to_axis(self, unit_selector=None, units=None, axis=0):
        if unit_selector is None:
            unit_selector = self.__class__._unit_selector

        if units is None:
            units = self.__class__._units

        units_fitter = unit_selector.set_default_units(units)

        new_frame = self._modify_axis(units_fitter, axis=axis)._modify_axis_name(
            units_fitter, axis=axis
        )

        return new_frame

    def fit_units_to_index(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(unit_selector=unit_selector, units=units, axis=0)

    def fit_units_to_columns(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(unit_selector=unit_selector, units=units, axis=1)

    def fit_units_to_axes(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(
            unit_selector=unit_selector, units=units, axis=1
        ).fit_units_to_axis(unit_selector=unit_selector, units=units, axis=0)

    def switch_columns_type(self):

        return self.switch_axis_type(axis=1)

    def _format_entry(self, obj, formatter=None):
        if formatter is not None:
            return formatter(obj)
        else:
            return self.__class__._label_formatter(obj)

    def plotted(
        self,
        filename=None,
        labels_list=None,
        colors_list=default_colors,
        height=None,
        width=None,
        x_axis_description=None,
        y_axis_description=None,
        subplots=False,
        legend_pos="north east",
        extra_commands=None,
        options=None,
        container=None,
        label=None,
        caption=None,
        smooth=False,
        picture=None,
        preview=None,
        *args,
        **kwargs,
    ):
        ",xlabel={$t$},x unit=\si{\second},"

        if height is None:
            height = self._default_height

        if width is None:
            width = self._default_width

        print("self++++++++++++++")
        print(self._ylabel)
        print("self++++++++++++++")
        self._raw_title = " "

        plotted_frame = self.copy()
        plotted_frame._ylabel = self._ylabel
        plotted_frame._raw_title = self._raw_title
        print("copy", plotted_frame._ylabel)

        col_idx = plotted_frame.columns

        latex_backend = self.__class__._latex_backend

        # title recognition
        # +++++++++++++++++++
        if isinstance(col_idx, pd.MultiIndex):
            if len(col_idx.get_level_values(0).unique()) == 1:

                plotted_frame._raw_title = col_idx.get_level_values(0).unique()[0]
                plotted_frame = plotted_frame.droplevel(0, axis=1)

        col_idx = plotted_frame.columns

        if y_axis_description is None and isinstance(col_idx, pd.MultiIndex):

            print("index2transform", col_idx.get_level_values(-1).unique())
            if len(col_idx.get_level_values(-1).unique()) == 1:

                print("tu jestem")
                label_raw = EntryWithUnit(col_idx.get_level_values(-1).unique()[0])
                if isinstance(label_raw, str):
                    ylabel = label_raw
                else:
                    ylabel = latex_backend(
                        NoEscape(
                            EntryWithUnit(
                                NoEscape(col_idx.get_level_values(-1).unique()[0])
                            )
                        )
                    )
                # print('detect')
                # print(ylabel)

                y_axis_description = "ylabel={$" + ylabel + "$},"
                y_axis_description = y_axis_description.replace("$$", "$")

                plotted_frame = plotted_frame.droplevel(-1, axis=1)
                plotted_frame._ylabel = ylabel
            else:
                print(" a tu jestem w else")
                ylabel_list = [
                    latex_backend(EntryWithUnit(label))
                    for label in col_idx.get_level_values(-1).unique()
                ]

                # print('ylabels',ylabel_list)

                ylabel = ", ".join(set(ylabel_list))

                # y_axis_description = 'ylabel={$' + ylabel + '$},'

                y_axis_description = "ylabel={$" + NoEscape(ylabel) + "$},"

                plotted_frame._ylabel = ylabel

        elif plotted_frame._ylabel is not None:
            ylabel = plotted_frame._ylabel

            # y_axis_description = 'ylabel={$' + ylabel + '$},'

            y_axis_description = "ylabel={$" + NoEscape(ylabel) + "$},"
            plotted_frame._ylabel = ylabel

        else:
            print("pure else")

            ylabel = plotted_frame._ylabel
            y_axis_description = ""

        if x_axis_description is None:

            _xlabel = plotted_frame.index.name
            if isinstance(_xlabel, str):
                xlabel_desc = _xlabel
            else:
                xlabel_desc = latex(_xlabel)

            x_axis_description = ",xlabel={$" + xlabel_desc + "$},"

        plotted_frame = plotted_frame.set_str_index_columns()
        plotted_frame = plotted_frame.rename_axis(None, axis=1)

        if filename is None:
            filename = f"{self.__class__._default_path}/plot{self.__class__.__name__}{next(self.__class__._floats_no_gen)}"

        if container is None:
            container = self.__class__._container

        if picture is None:
            picture = self.__class__._picture

        # print(ylabel)
        # print(type(ylabel))
        # print(y_axis_description)
        fig = plotted_frame.to_standalone_figure(
            filename,
            labels_list=labels_list,
            colors_list=colors_list,
            height=height,
            width=width,
            x_axis_description=x_axis_description.replace("$$", "$"),
            y_axis_description=y_axis_description.replace("$$", "$"),
            subplots=subplots,
            legend_pos=legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
            picture=picture,
        )

        ############################ to pack as method
        if caption is not None:
            fig.add_caption(NoEscape(caption))
            plotted_frame._caption = caption
        elif plotted_frame._caption is not None:
            fig.add_caption(NoEscape(plotted_frame._caption))
            plotted_frame._caption = plotted_frame._caption

        else:
            fig.add_caption(NoEscape(plotted_frame.__class__._caption))
            plotted_frame._caption = plotted_frame.__class__._caption

        caption = plotted_frame._caption
        #################################33 to as method

        ############################ to pack as method
        if subplots is not None:

            plotted_frame._subplot = subplots
        elif plotted_frame._subplot is not None:

            plotted_frame._subplot = plotted_frame._subplot

        else:
            plotted_frame._subplot = plotted_frame.__class__._subplot

        subplots = plotted_frame._subplot

        #################################33 to as method

        # ylabel = plotted_frame._ylabel

        if label is not None:
            AutoMarker.add_marker(plotted_frame._get_str_key(), label)
            fig.append(Label(label))
        else:
            auto_mrk = AutoMarker(plotted_frame).marker
            fig.append(Label(auto_mrk))

        if preview is None:
            preview = self.__class__._preview_mode

        if preview:
            plotted_frame.plot(ylabel=ylabel, subplots=subplots)
            print("==============")
            print(y_axis_description.replace("$$", "$"))
            print("==============")
            plt.ylabel((f"${ylabel}$").replace("$$", "$"))
            plt.title(plotted_frame._raw_title)
            plt.show()
            display(IPMarkdown(caption))
            container.append(fig)

        print("==============")
        print(ylabel)
        print("==============")

        plotted_frame._ylabel = ylabel
        return plotted_frame  # .plot(ylabel=ylabel,subplots=subplots)

    def reported(
        self,
        container=None,
        index=True,
        label=None,
        caption=None,
        multirow=True,
        column_format=None,
        longtable=None,
        style=None,
        *args,
        **kwargs,
    ):

        if container is None:
            container = self.__class__._container

        tab = DataTable(self)
        if style == "bpasts":
            tab = BPASTSDataTable(self)

        if caption is not None:
            tab.add_caption(NoEscape(caption))

        tab.add_table(
            index=index,
            multirow=multirow,
            column_format=column_format,
            longtable=longtable,
            **kwargs,
        )

        if label is not None:
            tab.append(Label(label))

        if label is not None:
            AutoMarker.add_marker(self.style.to_latex(), label)
            tab.append(Label(label))
        else:
            # old version
            # auto_mrk = AutoMarker(self.style.to_latex()).marker
            # new option
            auto_mrk = AutoMarker(self).marker
            tab.append(Label(auto_mrk))

        container.append(tab)

        return self.copy()


class AdaptableSeries(pd.Series, BasicFormattingTools):
    r"""
    Basic class for formatting data plots. It provides methods for setting options

    Arguments
    =========

    Methods
    =======

    Example
    =======

    """

    @property
    def _constructor(self):
        return AdaptableSeries

    @property
    def _constructor_expanddim(self):
        return AdaptableDataFrame

    def __init__(
        self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False
    ):

        super().__init__(
            data=data, index=index, dtype=dtype, name=name, copy=copy, fastpath=fastpath
        )
        self._reported = False


class AdaptableDataFrame(pd.DataFrame, BasicFormattingTools):

    @property
    def _constructor(self):
        return AdaptableDataFrame

    @property
    def _constructor_sliced(self):
        return AdaptableSeries

    @classmethod
    def _init_with_ops(
        cls, data=None, index=None, columns=None, dtype=None, copy=None, **kwargs
    ):

        raw_frame = cls(data=data, index=index, columns=columns, dtype=dtype, copy=copy)
        new_frame = raw_frame.applying_method(raw_frame, **kwargs)

        return cls(data=new_frame, index=index, columns=columns, dtype=dtype, copy=copy)

    @classmethod
    def formatted(cls, data=None, index=None, columns=None, dtype=None, copy=None):
        return cls._init_with_ops(
            data=data, index=index, columns=columns, dtype=dtype, copy=copy
        )

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=None):
        super().__init__(
            data=data, index=index, columns=columns, dtype=dtype, copy=copy
        )
        self._reported = False

        self._subplot = None
        self._caption = None
        self._prepared_fig = None
        self._ylabel = None
        self._raw_title = None

    def _get_str_key(self):
        # return self.to_latex()+f'subplot={self._subplot}, self._caption{self._caption} '
        return self.style.to_latex() + f"subplot={self._subplot}"


    def smooth_data(
        self,
        step= 2,
        interp_method= "cubic",
        window= 9,
        order= 6,
        interp_method_2= "linear"
    ):

        data_try = (
            self.iloc[::step]
            .reindex(self.index)
            .interpolate(method=interp_method, limit_direction="both")
            .rolling(window=window, center=True).mean()
            .interpolate(method=interp_method_2, limit_direction="both")
            .reset_index(drop=True)
        )

        y = data_try.to_numpy()
        peak_max, = argrelmax(y, order=order)
        peak_min, = argrelmin(y, order=order)
    
        return data_try#, peak_max, peak_min


class LatexDataFrame(AdaptableDataFrame):
    _applying_func = lambda obj: (obj).fit_units_to_axes().format_axes_names()

    @property
    def _constructor_sliced(self):
        return LatexSeries

    @property
    def _constructor(self):
        return LatexDataFrame


class LatexSeries(AdaptableSeries):

    @property
    def _constructor(self):
        return LatexSeries

    @property
    def _constructor_expanddim(self):
        return LatexDataFrame


class ComputationalErrorFrame(AdaptableDataFrame):
    # _applying_func = lambda obj: obj.join(((obj[obj.columns[1]]-obj[obj.columns[0]]).div(obj[obj.columns[0]],axis=0))).set_axis(list(obj.columns)+[Symbol('\\delta')],axis=1)
    _applying_func = lambda obj: obj.join(
        (obj[(obj.columns[1])] - obj[(obj.columns[0])])
        .div(obj[obj.columns[0]])
        .rename("Difference")
    )

    @property
    def _constructor(self):
        return ComputationalErrorFrame

    @property
    def _constructor_sliced(self):
        return ComputationalErrorSeries


class ComputationalErrorSeries(AdaptableSeries):

    @property
    def _constructor(self):
        return ComputationalErrorSeries

    @property
    def _constructor_expanddim(self):
        return ComputationalErrorFrame


class ParametersSummaryFrame(AdaptableDataFrame):
    _applying_func = (
        lambda frame: frame.abs()
        .max()
        .reset_index()
        .pivot(columns=["level_0", "level_2"], index=["level_1"])[0]
    )

    @property
    def _common_constructor_series(self):
        return ParameterSummarySeries


#     def _filtter(self, filtter=None):
#         if self.columns.nlevels == 2:
#             filtter = lambda frame: frame.abs().max().reset_index(
#                 level=1).pivot(columns=['level_1'])
#         else:
#             filtter =

#         return filtter


class ParameterSummarySeries(AdaptableSeries):

    @property
    def _common_constructor_series(self):
        return ParameterSummaryFrame


class NumericalAnalysisDataFrame(AdaptableDataFrame):
    """
    NumericalAnalysisDataFrame class allows you to present and format data in the form of a DataFrame.

    Example
        =======

            >>>from dynpy import *
            >>>from dynpy.solvers.linear import *
            >>>from sympy import *
            >>>from dynpy.utilities.adaptable import *

            >>>m, c, k, F = symbols('m c k F', positive = True)
            >>>t = Symbol('t')
            >>>x = Function('x')(t)

            >>>subs_data = {
                        m: 10,
                        c: 5,
                        k: 2,
                        }

            >>>tspan_1=np.linspace(0,10, 100)
            >>>ode_sys = ODESystem(odes=Matrix([k*x + F + c*x.diff(t)]), dvars=Matrix([x]), ode_order=1)

            >>>X = [x, x.diff(t)]
            >>>num_df = NumericalAnalysisDataFrame.from_model(ode_sys.solution.with_ics([1]),
                                              F,
                                              [1, 2, 5],
                                              reference_data=subs_data,
                                              coordinates=X,
                                              index=pd.Index(tspan_1,name=t),
                                              ics=None
                                              )

            >>>display(num_df)
    """

    _applying_func = None

    _metadata = ["_applying_func"]

    @classmethod
    def _init_with_ops(
        cls, data=None, index=None, columns=None, dtype=None, copy=None, **kwargs
    ):

        raw_frame = cls(
            data=data,
            index=index,
            columns=columns,
            dtype=dtype,
            copy=copy,
            func=lambda obj: obj,
        )

        new_frame = raw_frame.applying_method(raw_frame, **kwargs)

        return cls(
            data=new_frame,
            index=index,
            columns=columns,
            dtype=dtype,
            copy=copy,
            func=lambda obj: obj,
        )

    @classmethod
    def formatted(
        cls, data=None, index=None, columns=None, dtype=None, copy=None, **kwargs
    ):
        return cls._init_with_ops(
            data=data, index=index, columns=columns, dtype=dtype, copy=copy, **kwargs
        )

    @classmethod
    def from_model(
        cls,
        model,
        parameter,
        span,
        reference_data=None,
        coordinates=None,
        index=None,
        ics=None,
    ):
        """
        Class method that allows you to create and customize a DataFrame object based on the ODESystem object.
        The prepared DataFrame gives the possibility to perform numerical simulations depending on the defined parameter.

        Arguments
        =========

            model : ODESystem object

            parameter : Symbol
                Independent constant

            span : List or array
                Parameter values.

        Default Arguments
        =================

            reference_data : Dictionary
                The values of the constants of the equation.

            coordinates : List or array
                Column headers.

            index : pandas.Index
                Index values and index symbol.

            ics:
                Initial values.

        Example
        =======

            >>>from dynpy import *
            >>>from dynpy.solvers.linear import *
            >>>from sympy import *
            >>>from dynpy.utilities.adaptable import *

            >>>m, c, k, F = symbols('m c k F', positive = True)
            >>>t = Symbol('t')
            >>>x = Function('x')(t)

            >>>subs_data = {
                        m: 10,
                        c: 5,
                        k: 2,
                        }

            >>>tspan_1=np.linspace(0,10, 100)
            >>>ode_sys = ODESystem(odes=Matrix([k*x + F + c*x.diff(t)]), dvars=Matrix([x]), ode_order=1)

            >>>X = [x, x.diff(t)]
            >>>num_df = NumericalAnalysisDataFrame.from_model(ode_sys.solution.with_ics([1]),
                                              F,
                                              [1, 2, 5],
                                              reference_data=subs_data,
                                              coordinates=X,
                                              index=pd.Index(tspan_1,name=t),
                                              ics=None
                                              )

            >>>display(num_df)
        """

        if coordinates is None:
            coords = list(model.dvars)
        else:
            coords = coordinates

        params_list = [Eq(parameter, value, evaluate=False) for value in span]

        if isinstance(model, (list, tuple)):
            models_list = model
        else:
            models_list = [model]

        if reference_data is None:
            prepared_models = models_list
        else:
            prepared_models = [
                sim_model.subs({**reference_data, parameter: parameter}).doit()
                for sim_model in models_list
            ]

        num_cases = pd.MultiIndex.from_product([prepared_models, params_list, coords])

        return cls(
            data=ics,
            index=index,
            columns=num_cases,
        )

    def __init__(
        self,
        data=None,
        index=None,
        columns=None,
        model=None,
        ics=None,
        dtype=None,
        copy=None,
        **kwargs,
    ):
        # _try_evat='test'
        # print(f'custom init of {type(self)}')
        """
        Arguments
        =========

            data : List, DataFrame or DataSeries
                The data that will be in the final DataFrame.

            index : list or array
                Parameter values.

            columns : List
                Column headers.

        Example 1
        =========

            >>>from dynpy import *
            >>>from dynpy.solvers.linear import *
            >>>from sympy import *
            >>>from dynpy.utilities.adaptable import *

            >>>m = Symbol('m')
            >>>num_df = NumericalAnalysisDataFrame(data=[0, 1, 2, 3, 4, 5, 6], index=np.linspace(0, 6, 7), columns=[m])
            >>>display(num_df)

        Example 2
        =========

            >>>from dynpy import *
            >>>from dynpy.solvers.linear import *
            >>>from sympy import *
            >>>from dynpy.utilities.adaptable import *

            >>>m, c, k, F = symbols('m c k F', positive = True)
            >>>t = Symbol('t')
            >>>x = Function('x')(t)

            >>>tspan_1=np.linspace(0,10, 100)
            >>>subs_data_2 = {
                m: 10,
                c: 5,
                k: 2,
                F: 1
            }

            >>>ode_sys = ODESystem(odes=Matrix([k*x + F + c*x.diff(t)]), dvars=Matrix([x]), ode_order=1)
            >>>sol_num = ode_sys.subs(subs_data_2).solution.with_ics([1]).numerized().compute_solution(tspan_1)

            >>>num_df = NumericalAnalysisDataFrame(data=sol_num)
            >>>display(num_df)

        """
        super().__init__(
            data=data, index=index, columns=columns, dtype=dtype, copy=copy
        )
        #         self._numerical_model = model
        #         self._ics_list=ics
        self._comp_time = None

    @property
    def _constructor(self):
        return NumericalAnalysisDataFrame

    @property
    def _constructor_sliced(self):
        return NumericalAnalisysSeries

    def _spot_model(self, current_case):

        columns_index = self.columns

        if "model" in columns_index.names:
            current_models = columns_index.to_frame()[model][current_case]
        else:
            current_model = self._numerical_model

    def perform_simulations(
        self,
        model_level_name=0,
        coord_level_name=-1,
        ics=None,
        backend=None,
        dependencies=None,
        expand=False,
    ):
        """
        A method that allows you to perform numerical simulations on the created DataFrame.

        Example
        =======

        >>>sym_num = num_df.perform_simulations()

        >>>display(sym_num)
        >>>display(sym_num.plot())

        """

        # display(self.columns.droplevel(coord_level_name).unique())

        computed_data = self.copy()

        computed_data._comp_time = AdaptableDataFrame(
            columns=self.columns.droplevel(coord_level_name).unique()
        )

        for case_data in self.columns.droplevel(coord_level_name).unique():

            model = case_data[model_level_name]

            params_dict = {}

            for param_eq in case_data[1:]:

                params_dict[param_eq.lhs] = param_eq.rhs

            # numerized_model = model.numerized(params_dict, backend=backend,expand=expand)
            # try
            numerized_model = model.subs(params_dict).numerized(
                backend=backend, expand=expand
            )

            t_span = np.asarray((self.index))

            t0 = t_span[0]

            ics_series = self[case_data].T[t0]

            # print('xxxxxxxxxxxxxxxxxxx',ics_series)

            #             default_ics=model.default_ics()

            ics_list = []
            for coord in numerized_model.ics_dvars:
                given_value = ics_series[coord]

                # print('given',given_value,type(given_value))
                # print('cond: ',given_value is None,np.isnan( given_value))
                if given_value is None:
                    ics_val = model.default_ics()[coord]
                elif np.isnan(given_value):
                    ics_val = model.default_ics()[coord]
                else:
                    ics_val = given_value

                # print('selected value \n',ics_val)
                ics_list.append(ics_val)

            # print('ics list \n',ics_list)

            result = numerized_model.compute_solution(
                t_span, ics_list, params_values=params_dict
            )
            result_array = result.T.to_numpy()

            ######### extra dependencies handling - it should be reimplemented with new method #####
            if isinstance(dependencies, dict):

                for key, expr in dependencies.items():
                    lambdified_expr = lambdify(
                        (result.index.name, list(result.columns)),
                        expr.subs(params_dict).expand().doit(),
                        "numpy",
                    )

                    single_data = lambdified_expr(t_span, result_array)

                    result[key] = single_data
            #### END ##### extra dependencies handling

            computed_data[case_data] = result[computed_data[case_data].columns]

            sim_time = pd.Series(
                [result._get_comp_time()], index=["Computing time [s]"]
            )

            computed_data._comp_time[case_data] = sim_time

        return computed_data

    def compute_solution(
        self,
        t_span=None,
        ic_list=None,
        t_eval=None,
        params_values=None,
        method="RK45",
        derivatives=False,
        model_level_name=0,
        coord_level_name=-1,
        backend=None,
        dependencies=None,
        output=None,
        expand=False,
    ):

        computed_data = self.copy()

        if t_span is not None:
            computed_data = computed_data.reindex(t_span)

        result = computed_data.perform_simulations(
            model_level_name=model_level_name,
            coord_level_name=coord_level_name,
            ics=ic_list,
            backend=backend,
            dependencies=dependencies,
            expand=expand,
        )

        if output == "NA":
            return TimeDataFrame(result.droplevel(0))
        else:
            return TimeDataFrame(result)

    def with_ics(self, ics=None, ivar0=0, sol0=None):
        computed_data = self.copy()

        if isinstance(ics, (list, tuple)):
            ics_list = ics
            computed_data = self.insert_ics(ics_list=ics_list)
        elif isinstance(ics, dict):
            dvars_list = self.get_dvars_in_order(computed_data.columns)
            ics_list = [ics[coord] for coord in dvars_list]
            computed_data = self.insert_ics(ics_list=ics_list)
        else:
            raise TypeError(
                f"Expected list, tuple or dictionary, received: {type(ics)}"
            )

        return computed_data

    def insert_ics(self, ics_list):
        computed_data = self.copy()

        if len(ics_list) < len(computed_data.columns):
            param_span = len(computed_data.columns) // len(ics_list)
            ics_list = self.multiply_list(ics_list, param_span)
        computed_data.values[0] = ics_list
        return computed_data

    def multiply_list(self, ics_list, param_span):
        duplicated_list = []
        for _ in range(param_span):
            duplicated_list.extend(ics_list)
        return duplicated_list

    def get_dvars_in_order(self, columns):
        """This method is used to get dvars list in correct order"""
        seen = set()
        return [col[2] for col in columns if not (col[2] in seen or seen.add(col[2]))]


class NumericalAnalisysSeries(AdaptableSeries):

    @property
    def _constructor(self):
        return NumericalAnalisysSeries

    @property
    def _constructor_expanddim(self):
        return NumericalAnalysisDataFrame


# class LatexDataFrame(AdaptableDataFrame):
#     _applying_func = lambda obj: (obj).set_multiindex_columns().format_columns_names().set_multiindex_columns()

#     @classmethod
#     def _init_with_ops(cls,data=None, index=None, columns=None, dtype=None, copy=None,**kwargs):

#         raw_frame= cls(data=data, index=index, columns=columns, dtype=dtype, copy=copy)

#         print('_init_with_ops')
#         display(raw_frame)

#         new_frame=raw_frame.applying_method(raw_frame,**kwargs)

#         print('_init_without_ops')
#         display(new_frame)

#         return cls(data=new_frame, index=index, columns=columns, dtype=dtype, copy=copy,func=lambda obj: obj)

#     @classmethod
#     def formatted(cls,data=None, index=None, columns=None, dtype=None, copy=None,**kwargs):
#         return cls._init_with_ops(data=data, index=index, columns=columns, dtype=dtype, copy=copy)

#     @property
#     def _common_constructor_series(self):
#         return LatexSeries

# class LatexSeries(AdaptableSeries):
#     @property
#     def _common_constructor_series(self):
#         return LatexDataFrame


class TimeDomainMethods(DataMethods):

    def _set_comp_time(self, time):
        self._comp_time = time
        return None

    def _get_comp_time(self):

        try:
            obj = self._comp_time
        except:
            obj = None
        else:
            obj = self._comp_time

        return obj

    def gradient(self):

        data_gradient = np.gradient(self.to_numpy(), (self.index))

        return TimeSeries(data=data_gradient, index=self.index, name=self.name)

    def _pure_fft(self):

        spectrum = fft.fft(self.to_numpy())

        # f_span = fft.fftfreq(len(self.index), d=self.index[1] - self.index[0])

        spectrum = SpectrumSeries(data=spectrum)
        # spectrum.index.name = Symbol('f')

        return spectrum

    def is_uniformly_distributed(self):

        sample_length = max(self.index) - min(self.index) / len(self.index) - 1
        step_list = [
            self.index[i + 1] - self.index[i] for i in range(len(self.index) - 1)
        ]

        if all(np.round(element - step_list[0], 10) == 0 for element in step_list):
            return True
        else:
            return False

    def to_frequency_domain(self):

        spectrum = fft.fft(self.to_numpy())

        f_span = fft.fftfreq(len(self.index), d=self.index[1] - self.index[0])

        spectrum = SpectrumSeries(data=spectrum, index=(f_span), name=self.name)
        spectrum.index.name = Symbol("f")

        return spectrum.sort_index()


#     def to_frequency_domain(self):

#         spectrum={name:fft.fft(data) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#         return SpectrumFrame(data=spectrum,index=f_span)


class SpectrumSeries(AdaptableSeries, SpectralMethods):

    @property
    def _constructor(self):
        return SpectrumSeries

    @property
    def _constructor_expanddim(self):
        return SpectrumFrame

    #     def time_series(self):

    #         sampling_rate=self.index[1]-self.index[0]
    #         t_span=fft.fftfreq(len(self.index),d=sampling_rate)
    #         t_span=t_span[-1]+t_span

    #         timeseries={name:fft.ifftshift(fft.ifft( data ) )for name,data in self.items()}

    #         return TimeDataSeries(data=(timeseries),index=t_span)

    def _sort_for_ifft(self):

        return SpectrumSeries(fft.fftshift(self.to_numpy()))


class SpectrumFrame(AdaptableDataFrame, SpectralMethods):

    @property
    def _constructor(self):
        return SpectrumFrame

    @property
    def _constructor_sliced(self):
        return SpectrumSeries

    def to_tikz_plot(
        self,
        filename,
        labels_list=None,
        colors_list=["red", "blue", "orange"],
        x_axis_description=",xlabel={$f$},x unit=\si{\hertz},",
        y_axis_description=None,
        legend_pos="north east",
        extra_commands=None,
        smooth=False,
    ):

        if y_axis_description == None:
            y_axis_description = "ylabel={$" + self.name + "$},"

        return super().to_tikz_plot(
            filename=filename,
            labels_list=labels_list,
            colors_list=colors_list,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            legend_pos=legend_pos,
            extra_commands=extra_commands,
            smooth=smooth,
        )

    def to_standalone_figure(
        self,
        filename,
        labels_list=None,
        colors_list=default_colors,
        height=NoEscape(r"7cm"),
        width=NoEscape(r"0.9\textwidth"),
        x_axis_description=",xlabel={$f$},x unit=\si{\hertz},",
        y_axis_description="",
        subplots=False,
        legend_pos="north east",
        extra_commands=None,
        options=None,
        smooth=False,
    ):

        return super().to_standalone_figure(
            filename=filename,
            labels_list=labels_list,
            colors_list=colors_list,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            legend_pos=legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
        )

    def _sort_for_ifft(self):

        return SpectrumFrame(fft.fftshift(self), columns=self.columns, index=None)

    def double_sided_rms(self):

        spectrum_shifted_ds = {
            name: data.double_sided_rms() for name, data in self.items()
        }
        # f_span_shifted_ds = pd.Index(fft.fftshift(self.index),
        #                             name=self.index.name)
        f_span_shifted_ds = self.index

        return SpectrumFrame(data=spectrum_shifted_ds, index=f_span_shifted_ds)

    def single_sided_rms(self):

        spectrum_shifted_ss = {
            name: data.double_sided_rms()
            * np.heaviside(data.double_sided_rms().index, 0.5)
            * 2
            for name, data in self.items()
        }
        # f_span_shifted_ss = pd.Index(fft.fftshift(self.index),
        #                             name=self.index.name)
        f_span_shifted_ss = self.index

        return SpectrumFrame(data=spectrum_shifted_ss, index=f_span_shifted_ss)


#     def time_series(self):

#         sampling_rate=self.index[1]-self.index[0]
#         t_span=fft.fftshift(fft.fftfreq(len(self.index),d=sampling_rate))
#         t_span=t_span[-1]+t_span

#         timeseries={name:fft.ifftshift(fft.ifft( data ) )for name,data in self.items()}

#         return TimeDataFrame(data=(timeseries),index=t_span)

#     def is_uniformly_distributed(self):

#         sample_length=max(self.index)-min(self.index)/len(self.index)-1
#         step_list=[self.index[i+1] - self.index[i] for i in range(len(self.index)-1)]

#         if all(  np.round(element -  step_list[0],10) == 0 for element in step_list):
#             return True
#         else:
#             return False

#     def double_sided_shifted(self):

#         f_span=self.f_span

#         f_span_shifted_ds=fft.fftshift(f_span)
#         spectrum_shifted_ds={name:fft.fftshift(fft.fft(data)) for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

#     def single_sided_shifted(self):

#         f_span=self.f_span

#         f_span_shifted_ss=np.positive(fft.fftshift(f_span))
#         spectrum_shifted_ss={name:fft.fftshift(fft.fft(data))*np.heaviside([data],2) for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)


class TimeSeries(AdaptableSeries, TimeDomainMethods):

    @property
    def _constructor(self):
        return TimeSeries

    @property
    def _constructor_expanddim(self):
        return TimeDataFrame

    #     def spectrum(self):

    #         spectrum={name:(fft.fft(data)) for name,data in self.items()}
    #         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])
    #         return SpectrumSeries(data=spectrum,index=(f_span))

    def to_tikz_plot(
        self,
        filename,
        labels_list=None,
        colors_list=["red", "blue", "orange"],
        x_axis_description="xlabel={$t$},x unit=\si{\second},",
        y_axis_description=None,
        extra_commands=None,
        options=None,
        smooth=False,
    ):

        if y_axis_description == None:
            y_axis_description = "ylabel={$" + self.name + "$},"

        aux_DataFrame = TimeDataFrame(data=self, index=self.index)
        dumped_tex = aux_DataFrame.to_tikz_plot(
            filename=filename,
            labels_list=labels_list,
            colors_list=colors_list,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
        )
        return dumped_tex


class TimeDataFrame(AdaptableDataFrame, TimeDomainMethods):

    @property
    def _constructor(self):
        return TimeDataFrame

    @property
    def _constructor_sliced(self):
        return TimeSeries

    def _pure_fft(self):

        spectral_data = {name: data._pure_fft() for name, data in self.items()}

        return SpectrumFrame(data=spectral_data)

    def _resampled(self, sample_length, num_samples):

        from scipy import signal

        new_idx = np.linspace(0, sample_length, num_samples, endpoint=False)
        resampled_sig = {
            name: signal.resample(data.to_numpy(), num_samples)
            for name, data in self.items()
        }

        return TimeDataFrame(data=resampled_sig, index=new_idx)

    def to_frequency_domain(self):

        spectral_data = {
            name: data.to_frequency_domain() for name, data in self.items()
        }
        f_span = fft.fftfreq(len(self.index), d=self.index[1] - self.index[0])

        return SpectrumFrame(
            data=spectral_data, index=spectral_data[next(iter(spectral_data))].index
        ).sort_index()

    def gradient(self):

        data_gradient = {name: data.gradient() for name, data in self.items()}
        return TimeDataFrame(
            data=data_gradient, index=data_gradient[next(iter(data_gradient))].index
        )
