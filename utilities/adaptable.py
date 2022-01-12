from numpy import (fft)
import numpy as np
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Alignat, Package, Quantity, Command, Label, Table, Marker,Ref
from pylatex.utils import italic, NoEscape

from sympy import Matrix, ImmutableMatrix, symbols, Symbol, Eq, Expr, latex, Float, Function, Number
from sympy.core.relational import Relational

from sympy.physics.mechanics import vlatex

from IPython.display import display, Markdown, Latex

import pandas as pd

import matplotlib.pyplot as plt
import copy
#from number impo


def plots_no():
    num = 0
    while True:
        yield num
        num += 1


default_colors = ['red', 'blue', 'orange', 'teal', 'black', 'green']


class DataMethods:
    
    _figure_gen = lambda: Figure(position='H')
    _image_parameters={'width':NoEscape('0.9\textwidth')}
    _legend_fontsize = r' '
    _label_fontsize = r'\small '
    _template = Document(documentclass='standalone',geometry_options=None,document_options=["tikz"])
    #_extra

    @classmethod
    def set_document_template(cls,template=None):
        if template is not None:
            cls._template = template
        return cls


    @classmethod
    def set_default_label_fontsize(cls,fontsize=None):
        if fontsize is not None:
            cls._label_fontsize = fontsize
        return cls
    
    @classmethod
    def set_default_legend_fontsize(cls,fontsize=None):
        if fontsize is not None:
            cls._legend_fontsize = fontsize
        return cls
    
    
    @classmethod
    def set_default_figure_generator(cls,figure_generator=None):
        if figure_generator is not None:
            cls._figure_gen = figure_generator
        return cls
    
    @classmethod
    def set_default_image_parameters(cls,image_parameters=None):
        if image_parameters is not None:
            cls._image_parameters = image_parameters
        return cls    
    
    
    
    
    def _pylatex_tikz(self,
                      filename,
                      labels_list=None,
                      colors_list=default_colors,
                      height=NoEscape(r'7cm'),
                      width=NoEscape(r'0.9\textwidth'),
                      x_axis_description=',xlabel={$t$},x unit=\si{\second},',
                      y_axis_description='',
                      subplots=False,
                      extra_commands=None,
                      options=None,
                      smooth=False):

        if smooth:
            radius_str = NoEscape(f',rounded corners=1mm,')
        else:
            radius_str = ''

        labels_list = [label for label in self.columns]
        #labels_list=[(vlatex(label)) for label in self.columns]

        data_for_plot_list = [
            zip(self.index, plot_data)
            for label, plot_data in list(self.items())
        ]

        plots_no = len(data_for_plot_list)

        colors_multiplicator = np.ceil(plots_no / len(colors_list))

        plot_options = NoEscape(
            'anchor=north west,ymajorgrids=true,xmajorgrids=true,grid style=dashed,legend style={font=\small},'
            + y_axis_description) + NoEscape(',height=') + height + NoEscape(
                ',width=') + width + NoEscape(
                    f',xmin={min(self.index)},xmax={max(self.index)}')

        #with doc.create(Figure(position='!htb')) as fig:

        plot_options_list = [plot_options + NoEscape(',xticklabels=\empty,')
                             ] * (plots_no - 1)
        plot_options_list.append(plot_options + NoEscape(x_axis_description))

        tikzpicture = TikZ(options=options)

        
        legend_font_size=self.__class__._legend_fontsize
        
        if subplots == False:

            with tikzpicture.create(
                    Axis(options=plot_options_list[-1])) as plot:

                for data_for_plot, label, color in zip(
                        data_for_plot_list, labels_list,
                        colors_list * int(colors_multiplicator)):
                    coordinates = data_for_plot

                    plot.append(
                        Plot(name=NoEscape(NoEscape(legend_font_size + str(label))),
                             coordinates=coordinates,
                             options='color=' + color + ',solid' + radius_str))

        else:
            at_option = NoEscape('')
            for no, combined_plot_data in enumerate(
                    zip(data_for_plot_list, labels_list,
                        colors_list * int(colors_multiplicator))):
                data_for_plot, label, color = combined_plot_data
                plot_name = NoEscape(',name=plot' + str(no) + ',')
                with tikzpicture.create(
                        Axis(options=plot_options_list[no] + plot_name +
                             at_option)) as plot:
                    coordinates = data_for_plot
                    plot.append(
                        Plot(name=NoEscape(
                            NoEscape(legend_font_size) + NoEscape(str(label))),
                             coordinates=coordinates,
                             options='color=' + color + ',solid' + radius_str))
                    #at_option=NoEscape('at=(plot'+str(no)+'.below south west),')
                    at_option = NoEscape('at={(plot' + str(no) +
                                         '.south west)},yshift=-0.1cm,')

        if extra_commands is not None:
            tikzpicture.append(extra_commands)

        return tikzpicture

    def to_pylatex_plot(
            self,
            filename,
            labels_list=None,
            colors_list=default_colors,
            height=NoEscape(r'6cm'),
            width=NoEscape(r'0.49\textwidth'),
            x_axis_description=',xlabel={$t$},x unit=\si{\second},',
            y_axis_description='',
            subplots=False,
            extra_commands=None):

        tikz_pic = self._pylatex_tikz(filename,
                                      labels_list,
                                      colors_list,
                                      height,
                                      width,
                                      x_axis_description,
                                      y_axis_description,
                                      subplots,
                                      extra_commands=extra_commands)

        return tikz_pic

    def to_standalone_plot(
            self,
            filename,
            labels_list=None,
            colors_list=default_colors,
            height=NoEscape(r'6cm'),
            width=NoEscape(r'0.49\textwidth'),
            x_axis_description=',xlabel={$t$},x unit=\si{\second},',
            y_axis_description='',
            subplots=False,
            legend_pos='north east',
            extra_commands=None,
            options=None,
            smooth=False,
            picture=True,
            template=None,
            *arg,
            **kwargs):

        geometry_options = {
            "margin": "0cm",
        }
        
        if template is None:
            
            doc = self.__class__._template
            
        doc = copy.deepcopy(self.__class__._template)
        #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

        doc.packages.append(Package('siunitx'))
        doc.packages.append(Package('mathtools'))
        doc.packages.append(Package('float'))
        doc.packages.append(Package('tikz'))
        doc.packages.append(Package('pgfplots'))
        doc.append(Command('usepgfplotslibrary', arguments='units'))
        doc.append(Command('usetikzlibrary', arguments='spy'))

        doc.append(
            Command(
                'pgfplotsset',
                arguments=NoEscape(
                    r'compat=newest,label style={fontsize},legend pos='.format(fontsize='{'+f'font={self.__class__._label_fontsize}' +'}') +
                    str(legend_pos))))

        tikz_pic = self._pylatex_tikz(filename,
                                      labels_list,
                                      colors_list,
                                      height,
                                      width,
                                      x_axis_description,
                                      y_axis_description,
                                      subplots,
                                      extra_commands=extra_commands,
                                      options=options,
                                      smooth=smooth)

        doc.append(tikz_pic)
        
        if picture:
            doc.generate_pdf(filename, clean_tex=False)
        else:
            doc.generate_tex(filename)            
        
        return doc.dumps()

    def to_tikz_plot(self,
                     filename,
                     labels_list=None,
                     colors_list=default_colors,
                     height=NoEscape(r'6cm'),
                     width=NoEscape(r'0.49\textwidth'),
                     x_axis_description=',xlabel={$t$},x unit=\si{\second},',
                     y_axis_description='',
                     subplots=False,
                     legend_pos='north east',
                     extra_commands=None,
                     options=None,
                     smooth=False,
                     picture=True,
                     ):

        return self.to_standalone_plot(filename,
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
                                       picture=picture)

    def to_standalone_figure(
            self,
            filename,
            labels_list=None,
            colors_list=default_colors,
            height=NoEscape(r'6cm'),
            width=NoEscape(r'0.49\textwidth'),
            x_axis_description=',xlabel={$t$},x unit=\si{\second},',
            y_axis_description='',
            subplots=False,
            legend_pos='north east',
            extra_commands=None,
            options=None,
            smooth=False,
            picture=True,
            ):

        self.to_standalone_plot(filename,
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
                                picture=picture)
       
        fig = self.__class__._figure_gen()
        img_params=self.__class__._image_parameters
        
        if picture:
            fig.add_image(filename,width=NoEscape(r'0.9\textwidth'))
        else:    
            fig.append(Command(command='input',arguments=filename))

        return fig


class SpectralMethods(DataMethods):
    def is_uniformly_distributed(self):

        sample_length = max(self.index) - min(self.index) / len(self.index) - 1
        step_list = [
            self.index[i + 1] - self.index[i]
            for i in range(len(self.index) - 1)
        ]

        if all(
                np.round(element - step_list[0], 10) == 0
                for element in step_list):
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
            name: fft.ifftshift(fft.ifft(data))
            for name, data in self.items()
        }

        return TimeDataFrame(data=(timeseries), index=t_span)

    def double_sided_rms(self):

        f_span_shifted_ds = (fft.fftshift(self.index))
        spectrum_shifted_ds = fft.fftshift(abs(self) / len(self))

        return SpectrumSeries(data=spectrum_shifted_ds,
                              index=f_span_shifted_ds,
                              name=self.name)


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
    _latex_backend = vlatex

    @classmethod
    def set_default_units(cls, units={}):

        cls._units = units
        return cls

    def __new__(cls, obj, units=None, latex_backend=None, **kwargs):
        obj_with_unit = super().__new__(cls)
        obj_with_unit._obj = obj
        obj_with_unit._left_par = '['
        obj_with_unit._right_par = ']'

        if units is not None:
            obj_with_unit._units = units
        else:
            obj_with_unit._units = cls._units

        if isinstance(obj, Relational):
            obj_with_unit._left_par = ''
            obj_with_unit._right_par = ''
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

        if unit:
            return f'{entry_str} {left_par}{unit.__str__()}{right_par}'
        else:
            return f'{entry_str}'

    def __repr__(self):
        entry_str = self._obj.__repr__()
        unit = self._unit
        left_par = self._left_par
        right_par = self._right_par

        if unit:
            return f'{entry_str} {left_par}{unit.__repr__()}{right_par}'
        else:
            return f'{entry_str}'

    def _latex(self, *args):

        entry_str = self._latex_backend(self._obj)
        unit = self._unit
        left_par = self._left_par
        right_par = self._right_par

        if unit:
            return f'{entry_str} {left_par}{unit:~L}{right_par}'
        else:
            return f'{self._obj}'

        
class DataTable(Table):
    _latex_name = 'table'

    def __init__(self, numerical_data, position=None):
        super().__init__(position=position)
        ##print(numerical_data)
        self._numerical_data = numerical_data
        self.position = position

    def add_table(self, numerical_data=None, index=False, longtable=False):
        self.append(NoEscape('\\centering'))
        self.append(NoEscape('%%%%%%%%%%%%%% Table %%%%%%%%%%%%%%%'))
        #         if numerical_data!=None:
        #             self._numerical_data=numerical_data

        tab = self._numerical_data
        self.append(
            NoEscape(tab.to_latex(index=index, escape=False, longtable=longtable).replace('\\toprule','\\toprule \n \\midrule').replace('\\bottomrule','\\midrule \n \\bottomrule'))
            )
        
        
class MarkerRegistry(dict):
    
    _prefix = 'automrk'
    _markers_dict={}
    
    def __init__(self,prefix='automrk',sufix=None):
        super().__init__()
        self._prefix = prefix


class AutoMarker:
    _markers_dict={}
    _prefix = 'eq'
    _name = None
    _floats_no_gen = plots_no()
    
    @classmethod
    def add_marker(cls,elem,marker):
        
        if isinstance(elem,(AdaptableDataFrame)):
            elem_id = elem._get_str_key()
            if elem._subplot is None:
                prefix = 'tab'
            else:
                prefix = 'fig'
        
        elif isinstance(elem,(pd.DataFrame,pd.Series)):
            elem_id = elem.to_latex()
            prefix = 'fig'
            
        elif isinstance(elem,Matrix):
            elem_id = ImmutableMatrix(elem)
            prefix='eq'
        else:
            elem_id = elem
        
        
        cls._markers_dict[elem_id]=marker
        
        return None
    
    
    def __init__(self,elem,prefix='eq',name=None,sufix=None):

        
        if isinstance(elem,(AdaptableDataFrame)):
            elem_id = elem._get_str_key()
            if elem._subplot is None:
                prefix = 'tab'
            else:
                prefix = 'fig'
        
        elif isinstance(elem,(pd.DataFrame,pd.Series)):
            elem_id = elem.to_latex()
            prefix = 'fig'
            
        elif isinstance(elem,Matrix):
            elem_id = ImmutableMatrix(elem)
            prefix='eq'
        else:
            elem_id = elem

            
        self._marker_name = elem.__class__.__name__
        self._elem_id = elem_id
        self._prefix = prefix
        self._sufix = sufix
            
        self._get_marker()

        

    def _get_marker(self,elem=None):
        if elem is None:
            elem = self._elem_id
            
        available_markers=self._markers_dict
        

            
        if elem in available_markers:
            marker = available_markers[elem]
        else:
            

            marker = Marker(f'Mrk{self._marker_name}{next(self._floats_no_gen)}',prefix=self._prefix)
            self._markers_dict[elem] = marker
        

        
        self._marker = marker
        return marker
    
    @property
    def marker(self):
        return self._marker
        
    def __repr__(self):
        return f'AutoMarker for {self._marker}'
    
    def __str__(self):
        return Ref(self._marker).dumps()
        
        

class BasicFormattingTools(DataMethods):

    _floats_no_gen = plots_no()

    _latex_backend = vlatex
    _label_formatter = lambda obj: f'${_latex_backend(obj)}$' if isinstance(
        obj, (Expr, Eq, EntryWithUnit)) else obj
    _unit_selector = EntryWithUnit
    _domain = None
    _units = {}
    _applying_func = lambda x: x
    _init_ops = True

    _default_sep = ', '
    _container = []
    _default_path = './tikzplots'
    _picture=False
    
    _subplot=False
    _caption = 'Default caption'
    
    _default_width =  NoEscape(r'0.9\textwidth')
    _default_height = NoEscape(r'6cm')
    _preview_mode = False
    
    @classmethod
    def set_default_width(cls, width=  NoEscape(r'0.9\textwidth')):
        cls._default_width = width
        return cls
    
    @classmethod
    def set_preview_mode(cls, preview=False):
        cls._preview_mode = preview
        return cls

    
    @classmethod
    def set_default_height(cls, height = NoEscape(r'6cm')):
        cls._default_height = height

        return cls


    
    @classmethod
    def set_default_column_separator(cls, sep=', '):
        cls._default_sep = sep

        return cls

    @classmethod
    def set_picture_mode(cls, picture=False):

        cls._picture=picture
        return cls
    
    @classmethod
    def set_directory(cls, path='./'):

        cls._default_path = path
        return cls
    
    @classmethod
    def set_default_container(cls, container=[]):
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

        #print('type',type(obj))

        #print(isinstance(obj, (Expr, Eq, EntryWithUnit)))
        
        latex_backend = BasicFormattingTools._latex_backend

        if isinstance(obj, (Symbol, Function, Expr, Eq, EntryWithUnit)):
            formatted_obj = f'${latex_backend(obj)}$'
        else:
            formatted_obj = obj
        #print('effect',formatted_str)

        return formatted_obj

    def applying_method(self, data, func=None, **kwargs):
        if func:
            #print('func is used')
            ops_func = func

        elif self.__class__._applying_func is not None:
            #print('class func is used')

            ops_func = self.__class__._applying_func
        else:
            #print('identity is used')
            ops_func = lambda data: data

        return ops_func(data)

    def _modify_axis(self, func, axis=0):
        
        

        new_obj = self.copy()
        new_obj_idx = new_obj.axes[axis]
        idx_frame = new_obj_idx.to_frame().applymap(func)

        #print('idx',new_obj_idx)
        
        #display('map',idx_frame)
        
        if isinstance(new_obj_idx, pd.MultiIndex):
            new_obj_idx = pd.MultiIndex.from_frame(idx_frame)
            #new_obj_idx.names=map(func,new_obj_idx.names)

        else:

            #new_obj_idx = pd.Index((idx_frame),name=new_obj_idx.name)
            new_obj_idx = new_obj_idx.map(func)

        return new_obj.set_axis(new_obj_idx, axis=axis)

    def _modify_axis_name(self, func, axis=0):
        new_obj = self.copy()
        new_obj_idx = new_obj.axes[axis]

        if isinstance(new_obj_idx, pd.MultiIndex):

            new_obj_idx.names = map(func, new_obj_idx.names)

        else:

            #new_obj_idx = pd.Index((idx_frame),name=new_obj_idx.name)
            new_obj_idx = new_obj_idx.map(func)
            new_obj_idx.name = func(new_obj_idx.name)

        return new_obj.set_axis(new_obj_idx, axis=axis)

    def set_multiindex_axis(self, axis=0):

        if axis == 'index':
            axis = 0
        elif axis == 'colums':
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

        if axis == 'index':
            axis = 0
        elif axis == 'colums':
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

        if axis == 'index':
            axis = 0
        elif axis == 'colums':
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

            #formatter = self.__class__._label_formatter
            formatter = self._label_formatter

        new_frame = self._modify_axis(formatter,
                                      axis=axis)._modify_axis_name(formatter,
                                                                   axis=axis)

        return new_frame

    def format_index_names(self, formatter=None):

        return self.format_axis_names(formatter=formatter, axis=0)

    def format_columns_names(self, formatter=None):

        return self.format_axis_names(formatter=formatter, axis=1)

    def format_axes_names(self, formatter=None):

        return self.format_axis_names(
            formatter=formatter, axis=1).format_axis_names(formatter=formatter,
                                                           axis=0)

    def set_str_index_axis(self, axis=0, separator=', '):

        new_obj = self.copy()
        new_idx = new_obj.axes[axis]

        if isinstance(new_idx, pd.MultiIndex):
            new_obj = self.copy().set_flat_index_axis(axis=axis)

            new_obj = new_obj.set_axis(
                [separator.join(map(str, entry)) for entry in new_idx],
                axis=axis)

        return new_obj

    def set_str_index_columns(self, separator=', '):

        return self.set_str_index_axis(axis=1, separator=separator)

    def to_eq_axis_form(self, axis=0):
        new_obj = self.copy()

        return new_obj

    def to_named_axis_form(self, axis=0):
        new_obj = self.copy()

        idx = new_obj.axes[axis]

        if isinstance(idx, pd.MultiIndex):
            #print('MultiIndex modification has been not supported yet')
            numbered_obj=new_obj.format_axis_names(lambda name: float(name.rhs) if isinstance(name,Eq) else name  ,axis=axis)
            #sym_wyn.loc[:,ix[:,:,dyn_sys.phi_1]].plot()

            new_obj_with_name=numbered_obj.format_columns_names(lambda name: float(name.rhs) if isinstance(name,Eq) else name  )
            new_obj_with_name.axes[axis].names=[ level_idx[0].lhs if isinstance(level_idx[0],Eq) else  level_idx.name  for level_idx  in idx.levels ]
            

        else:
            name = idx.name
            
            if all([isinstance(entry, Relational) for entry in idx]):
                #TODO: add numerical types recognizing
                new_obj_with_name = new_obj._modify_axis(
                    lambda elem: float(elem.rhs) if isinstance(elem.rhs,Number) else elem.rhs, axis=axis)

                #print(new_obj_with_name.axes[axis].name)
                new_obj_with_name.axes[axis].name = idx[0].lhs

        return new_obj_with_name  #.rename(Symbol('l_w'),axis=axis)

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

        new_frame = self._modify_axis(
            units_fitter, axis=axis)._modify_axis_name(units_fitter, axis=axis)

        return new_frame

    def fit_units_to_index(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(unit_selector=unit_selector,
                                      units=units,
                                      axis=0)

    def fit_units_to_columns(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(unit_selector=unit_selector,
                                      units=units,
                                      axis=1)

    def fit_units_to_axes(self, unit_selector=None, units=None):
        return self.fit_units_to_axis(unit_selector=unit_selector,
                                      units=units,
                                      axis=1).fit_units_to_axis(
                                          unit_selector=unit_selector,
                                          units=units,
                                          axis=0)

    def switch_columns_type(self):

        return self.switch_axis_type(axis=1)

    def _format_entry(self, obj, formatter=None):
        if formatter is not None:
            return formatter(obj)
        else:
            return self.__class__._label_formatter(obj)

    def plotted(self,
                filename=None,
                labels_list=None,
                colors_list=default_colors,
                height=None,
                width=None,
                x_axis_description=None,
                y_axis_description=None,
                subplots=False,
                legend_pos='north east',
                extra_commands=None,
                options=None,
                container=None,
                label=None,
                caption=None,
                smooth=False,
                picture=None,
                preview=None,
                *args,
                **kwargs):
        
        ',xlabel={$t$},x unit=\si{\second},'

        if height is None:
            height = self._default_height
            
        if width is None:
            width = self._default_width
        
        #print(self._ylabel)
        
        plotted_frame = self.copy()
        #print('copy',plotted_frame._ylabel)
        
        col_idx = plotted_frame.columns

        latex_backend = self.__class__._latex_backend

        if y_axis_description is None and isinstance(col_idx, pd.MultiIndex):

            #print('index2transform',col_idx.get_level_values(-1).unique())
            if len(col_idx.get_level_values(-1).unique())==1:
                ylabel = latex_backend(
                    EntryWithUnit(col_idx.get_level_values(-1).unique()[0]))

                y_axis_description = 'ylabel={' + ylabel + '},'


                plotted_frame = plotted_frame.droplevel(-1, axis=1)
            else:
                ylabel_list  = [latex_backend(EntryWithUnit(label))   for label in col_idx.get_level_values(-1).unique()]
                
                #print('ylabels',ylabel_list)
                
                ylabel = ', '.join(ylabel_list)

                y_axis_description = 'ylabel={' + ylabel + '},'                
                

            plotted_frame._ylabel=ylabel

            
            
            
        elif self._ylabel is not None:
            ylabel = self._ylabel
            y_axis_description = 'ylabel={' + ylabel + '},'
            

        else:
            
            y_axis_description = ''

        if x_axis_description is None:
            x_axis_description = ',xlabel={' + plotted_frame.index.name  +  '},'
        
            
            
        plotted_frame = plotted_frame.set_str_index_columns()
        plotted_frame = plotted_frame.rename_axis(None,axis=1)

        if filename is None:
            filename = f'{self.__class__._default_path}/plot{self.__class__.__name__}{next(self.__class__._floats_no_gen)}'

        if container is None:
            container = self.__class__._container

        if picture is None:
            picture = self.__class__._picture
            
            
        fig = plotted_frame.to_standalone_figure(
            filename,
            labels_list=labels_list,
            colors_list=colors_list,
            height=height,
            width=width,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            subplots=subplots,
            legend_pos=legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth,
            picture = picture)

        
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
        plotted_frame._ylabel=ylabel
 





        if label is not None:
            AutoMarker.add_marker(plotted_frame._get_str_key(),label)
            fig.append(Label(label))
        else:
            auto_mrk=AutoMarker(plotted_frame).marker
            fig.append(Label(auto_mrk))

        
        
        if preview is None:
            preview = self.__class__._preview_mode
        
        if preview:
            plotted_frame.plot(ylabel=ylabel,subplots=subplots)
            plt.show()
            display(Markdown(caption))
            container.append(fig)
        
        return plotted_frame#.plot(ylabel=ylabel,subplots=subplots)



        
        
    
    
    def reported(self,
                 container=None,
                 index=True,
                 label=None,
                 caption=None,
                 *args,
                 **kwargs):

        if container is None:
            container = self.__class__._container

        tab = DataTable(self)

        if caption is not None:
            tab.add_caption(NoEscape(caption))

        tab.add_table(index=index)

        if label is not None:
            tab.append(Label(label))

        if label is not None:
            AutoMarker.add_marker(self.to_latex(),label)
            tab.append(Label(label))
        else:
            auto_mrk=AutoMarker(self.to_latex()).marker
            tab.append(Label(auto_mrk))
            

        container.append(tab)
        
        

        return self.copy()


class AdaptableSeries(pd.Series, BasicFormattingTools):
    r'''
    Basic class for formatting data plots. It provides methods for setting options 
    
    Arguments
    =========

    Methods
    =======

    Example
    =======

    '''
    @property
    def _constructor(self):
        return AdaptableSeries

    @property
    def _constructor_expanddim(self):
        return AdaptableDataFrame

    def __init__(self,
                 data=None,
                 index=None,
                 dtype=None,
                 name=None,
                 copy=False,
                 fastpath=False):

        super().__init__(data=data,
                         index=index,
                         dtype=dtype,
                         name=name,
                         copy=copy,
                         fastpath=fastpath)
        self._reported = False


class AdaptableDataFrame(pd.DataFrame, BasicFormattingTools):
    @property
    def _constructor(self):
        return AdaptableDataFrame

    @property
    def _constructor_sliced(self):
        return AdaptableSeries

    @classmethod
    def _init_with_ops(cls,
                       data=None,
                       index=None,
                       columns=None,
                       dtype=None,
                       copy=None,
                       **kwargs):

        raw_frame = cls(data=data,
                        index=index,
                        columns=columns,
                        dtype=dtype,
                        copy=copy)
        new_frame = raw_frame.applying_method(raw_frame, **kwargs)

        return cls(data=new_frame,
                   index=index,
                   columns=columns,
                   dtype=dtype,
                   copy=copy)

    @classmethod
    def formatted(cls,
                  data=None,
                  index=None,
                  columns=None,
                  dtype=None,
                  copy=None):
        return cls._init_with_ops(data=data,
                                  index=index,
                                  columns=columns,
                                  dtype=dtype,
                                  copy=copy)

    def __init__(self,
                 data=None,
                 index=None,
                 columns=None,
                 dtype=None,
                 copy=None):
        super().__init__(data=data,
                         index=index,
                         columns=columns,
                         dtype=dtype,
                         copy=copy)
        self._reported = False
        
        self._subplot = None
        self._caption = None
        self._prepared_fig = None
        self._ylabel = None


    def _get_str_key(self):
        #return self.to_latex()+f'subplot={self._subplot}, self._caption{self._caption} '
        return self.to_latex()+f'subplot={self._subplot}'
    
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
    #_applying_func = lambda obj: obj.join(((obj[obj.columns[1]]-obj[obj.columns[0]]).div(obj[obj.columns[0]],axis=0))).set_axis(list(obj.columns)+[Symbol('\\delta')],axis=1)
    _applying_func = lambda obj: obj.join((obj[(obj.columns[1])] - obj[
        (obj.columns[0])]).div(obj[obj.columns[0]]).rename('Difference'))

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
    _applying_func = lambda frame: frame.abs().max().reset_index().pivot(
        columns=['level_0', 'level_2'], index=['level_1'])[0]

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
    _applying_func = None

    _metadata = ["_applying_func"]

    @classmethod
    def _init_with_ops(cls,
                       data=None,
                       index=None,
                       columns=None,
                       dtype=None,
                       copy=None,
                       **kwargs):

        raw_frame = cls(data=data,
                        index=index,
                        columns=columns,
                        dtype=dtype,
                        copy=copy,
                        func=lambda obj: obj)

        print('_init_with_ops')
        display(raw_frame)

        new_frame = raw_frame.applying_method(raw_frame, **kwargs)

        print('_init_without_ops')
        display(new_frame)

        return cls(data=new_frame,
                   index=index,
                   columns=columns,
                   dtype=dtype,
                   copy=copy,
                   func=lambda obj: obj)

    @classmethod
    def formatted(cls,
                  data=None,
                  index=None,
                  columns=None,
                  dtype=None,
                  copy=None,
                  **kwargs):
        return cls._init_with_ops(data=data,
                                  index=index,
                                  columns=columns,
                                  dtype=dtype,
                                  copy=copy,
                                  **kwargs)

    def __init__(self,
                 data=None,
                 index=None,
                 columns=None,
                 model=None,
                 ics=None,
                 dtype=None,
                 copy=None,
                 **kwargs):
        #_try_evat='test'
        #print(f'custom init of {type(self)}')

        super().__init__(data=data,
                         index=index,
                         columns=columns,
                         dtype=dtype,
                         copy=copy)
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

        if 'model' in columns_index.names:
            current_models = columns_index.to_frame()[model][current_case]
        else:
            current_model = self._numerical_model

    def perform_simulations(self,
                            model_level_name=0,
                            coord_level_name=-1,
                            ics=None,
                            backend=None):

        display(self.columns.droplevel(coord_level_name).unique())

        computed_data = self.copy()

        computed_data._comp_time = AdaptableDataFrame(
            columns=self.columns.droplevel(coord_level_name).unique())

        for case_data in self.columns.droplevel(coord_level_name).unique():

            model = case_data[model_level_name]

            params_dict = {}

            for param_eq in case_data[1:]:

                params_dict[param_eq.lhs] = param_eq.rhs

            numerized_model = model.numerized(params_dict, backend=backend)

            t_span = np.asarray((self.index))
            print(type(t_span))

            t0 = t_span[0]

            ics_series = (self[case_data].T[t0])

            ics_list = [
                np.float(ics_series[coord])
                for coord in numerized_model.ics_dvars
            ]

            print(ics_list)

            result = numerized_model.compute_solution(t_span, ics_list)

            computed_data[case_data] = result[computed_data[case_data].columns]

            sim_time = pd.Series([result._get_comp_time()],
                                 index=['Computing time [s]'])

            computed_data._comp_time[case_data] = sim_time

        return (computed_data)


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

    def is_uniformly_distributed(self):

        sample_length = max(self.index) - min(self.index) / len(self.index) - 1
        step_list = [
            self.index[i + 1] - self.index[i]
            for i in range(len(self.index) - 1)
        ]

        if all(
                np.round(element - step_list[0], 10) == 0
                for element in step_list):
            return True
        else:
            return False

    def to_frequency_domain(self):

        spectrum = (fft.fft(self.to_numpy()))

        f_span = fft.fftfreq(len(self.index), d=self.index[1] - self.index[0])

        spectrum = SpectrumSeries(data=spectrum,
                                  index=(f_span),
                                  name=self.name)
        spectrum.index.name = Symbol('f')

        return spectrum


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


class SpectrumFrame(AdaptableDataFrame, SpectralMethods):
    @property
    def _constructor(self):
        return SpectrumFrame

    @property
    def _constructor_sliced(self):
        return SpectrumSeries

    def to_tikz_plot(self,
                     filename,
                     labels_list=None,
                     colors_list=['red', 'blue', 'orange'],
                     x_axis_description=',xlabel={$f$},x unit=\si{\hertz},',
                     y_axis_description=None,
                     legend_pos='north east',
                     extra_commands=None,
                     smooth=False):

        if y_axis_description == None:
            y_axis_description = 'ylabel=$' + self.name + '$,'

        return super().to_tikz_plot(filename=filename,
                                    labels_list=labels_list,
                                    colors_list=colors_list,
                                    x_axis_description=x_axis_description,
                                    y_axis_description=y_axis_description,
                                    legend_pos=legend_pos,
                                    extra_commands=extra_commands,
                                    smooth=smooth)

    def to_standalone_figure(
            self,
            filename,
            labels_list=None,
            colors_list=default_colors,
            height=NoEscape(r'7cm'),
            width=NoEscape(r'0.9\textwidth'),
            x_axis_description=',xlabel={$f$},x unit=\si{\hertz},',
            y_axis_description='',
            subplots=False,
            legend_pos='north east',
            extra_commands=None,
            options=None,
            smooth=False):

        return super().to_standalone_figure(
            filename=filename,
            labels_list=labels_list,
            colors_list=colors_list,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            legend_pos=legend_pos,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth)

    def double_sided_rms(self):

        spectrum_shifted_ds = {
            name: data.double_sided_rms()
            for name, data in self.items()
        }
        f_span_shifted_ds = pd.Index(fft.fftshift(self.index),
                                     name=self.index.name)

        return SpectrumFrame(data=spectrum_shifted_ds, index=f_span_shifted_ds)

    def single_sided_rms(self):

        spectrum_shifted_ss = {
            name: data.double_sided_rms() *
            np.heaviside(data.double_sided_rms().index, 0.5) * 2
            for name, data in self.items()
        }
        f_span_shifted_ss = pd.Index(fft.fftshift(self.index),
                                     name=self.index.name)

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

    def to_tikz_plot(self,
                     filename,
                     labels_list=None,
                     colors_list=['red', 'blue', 'orange'],
                     x_axis_description='xlabel={$t$},x unit=\si{\second},',
                     y_axis_description=None,
                     extra_commands=None,
                     options=None,
                     smooth=False):

        if y_axis_description == None:
            y_axis_description = 'ylabel=$' + self.name + '$,'

        aux_DataFrame = TimeDataFrame(data=self, index=self.index)
        dumped_tex = aux_DataFrame.to_tikz_plot(
            filename=filename,
            labels_list=labels_list,
            colors_list=colors_list,
            x_axis_description=x_axis_description,
            y_axis_description=y_axis_description,
            extra_commands=extra_commands,
            options=options,
            smooth=smooth)
        return dumped_tex


class TimeDataFrame(AdaptableDataFrame, TimeDomainMethods):
    @property
    def _constructor(self):
        return TimeDataFrame

    @property
    def _constructor_sliced(self):
        return TimeSeries

    def to_frequency_domain(self):

        spectral_data = {
            name: data.to_frequency_domain()
            for name, data in self.items()
        }
        f_span = fft.fftfreq(len(self.index), d=self.index[1] - self.index[0])

        return SpectrumFrame(data=spectral_data,
                             index=spectral_data[next(
                                 iter(spectral_data))].index)

    def gradient(self):

        data_gradient = {name: data.gradient() for name, data in self.items()}
        return TimeDataFrame(data=data_gradient,
                             index=data_gradient[next(
                                 iter(data_gradient))].index)
