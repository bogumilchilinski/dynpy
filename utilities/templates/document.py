from pylatex import (Document, Package,
                     #Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
                    )
#from pylatex.section import Paragraph, Chapter
from pylatex.utils import (#italic, 
                           NoEscape)


class TikzCaSCStandalone(Document):
    latex_name = 'document'
    packages = [
        Package('natbib', options=['numbers']),

    ]

    def __init__(self,
                 default_filepath='default_filepath',
                 *,
                 documentclass='standalone',
                 document_options=NoEscape('tikz,class=cas-sc'),
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=False,
                 textcomp=True,
                 microtype=None,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

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


class TikzStandalone(Document):
    latex_name = 'document'
    packages = [
        #         Package('siunitx'),
        #         Package('amsmath'),
        #         Package('float'),
        #         Package('tikz'),
        #         Package('pgfplots'),
        #         Command('usepgfplotslibrary', arguments='units'),
        #         Command('usetikzlibrary', arguments='spy')
    ]

    def __init__(self,
                 default_filepath='default_filepath',
                 documentclass='standalone',
                 document_options=[NoEscape('tikz')],
                 fontenc='T1',
                 inputenc='utf8',
                 font_size='normalsize',
                 lmodern=True,
                 textcomp=True,
                 microtype=None,
                 page_numbers=True,
                 indent=None,
                 geometry_options=None,
                 data=None):

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

