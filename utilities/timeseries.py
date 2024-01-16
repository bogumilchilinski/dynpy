from pandas import (Series,DataFrame)
from numpy import (fft)
import numpy as np
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure , Alignat, Package,Quantity, Command
from pylatex.utils import italic,NoEscape
from sympy.physics.mechanics import vlatex

from .adaptable import DataMethods, SpectralMethods, TimeDomainMethods, SpectrumSeries, SpectrumFrame, TimeSeries, TimeDataFrame

default_colors=['red','blue','orange','teal','black','green']