import numpy as np
from numpy import fft
from pandas import DataFrame, Series
from pylatex import (
    Alignat,
    Axis,
    Command,
    Document,
    Figure,
    Math,
    Package,
    Plot,
    Quantity,
    Section,
    Subsection,
    Tabular,
    TikZ,
)
from pylatex.utils import NoEscape, italic
from sympy.physics.mechanics import vlatex

from .adaptable import (
    DataMethods,
    SpectralMethods,
    SpectrumFrame,
    SpectrumSeries,
    TimeDataFrame,
    TimeDomainMethods,
    TimeSeries,
)

default_colors = ["red", "blue", "orange", "teal", "black", "green"]
