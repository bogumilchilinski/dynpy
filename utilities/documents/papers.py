import datetime
import os
import shutil
from typing import List, Optional, Union

from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    NewPage,
    Package,
    Tabularx,
)
from pylatex.base_classes import Environment

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape, dumps_list  # italic,

from ..report import (
    CurrentContainer,
    IPMarkdown,
    Markdown,
    ObjectCode,
    ReportText,
    display,
)

#1217
class Article(Document):
    pass