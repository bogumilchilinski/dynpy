from numpy import (fft)
import numpy as np
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure , Alignat, Package,Quantity, Command
from pylatex.utils import italic,NoEscape

from sympy import Matrix, symbols, Symbol, Eq, Expr,latex
from sympy.core.relational import Relational

from sympy.physics.mechanics import vlatex

import pandas as pd

default_colors=['red','blue','orange','teal','black','green']

class DataMethods:
    
    
    
    def _pylatex_tikz(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.9\textwidth'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False,extra_commands=None,options=None):

        
       

        
        labels_list=[label for label in self.columns]
        #labels_list=[(vlatex(label)) for label in self.columns]

        data_for_plot_list=[zip(self.index,plot_data)  for label,plot_data in list(self.items())]

        plots_no=len(data_for_plot_list)
        
        colors_multiplicator=np.ceil(plots_no/len(colors_list))

        plot_options = NoEscape('anchor=north west,ymajorgrids=true,xmajorgrids=true,grid style=dashed,legend style={font=\small},'+y_axis_description)+NoEscape(',height=')+height+NoEscape(',width=')+width+NoEscape(f',xmin={min(self.index)},xmax={max(self.index)}')
        
        #with doc.create(Figure(position='!htb')) as fig:

        plot_options_list=[plot_options+NoEscape(',xticklabels=\empty,')]*(plots_no-1)
        plot_options_list.append( plot_options+NoEscape(x_axis_description) )
        
        
        
        tikzpicture = TikZ(options=options)
        
        if subplots==False:

            with tikzpicture.create(Axis(options=plot_options_list[-1])) as plot:

                for data_for_plot,label,color in zip(data_for_plot_list,labels_list,colors_list*int(colors_multiplicator)) :
                    coordinates = data_for_plot

                    plot.append(Plot(name=NoEscape(NoEscape(r'\tiny '+label)), coordinates=coordinates,options='color='+color+',solid'+NoEscape(f',rounded corners=1mm,')))

        else:
            at_option=NoEscape('')
            for no,combined_plot_data in enumerate(zip(data_for_plot_list,labels_list,colors_list*int(colors_multiplicator))):
                data_for_plot,label,color = combined_plot_data
                plot_name=NoEscape(',name=plot'+str(no)+',')
                with tikzpicture.create(Axis(options=plot_options_list[no]+plot_name+at_option )) as plot:
                    coordinates = data_for_plot
                    plot.append(Plot(name=NoEscape(NoEscape(r'\tiny '+label)), coordinates=coordinates,options='color='+color+',solid'+NoEscape(f',rounded corners=1mm,')))
                    #at_option=NoEscape('at=(plot'+str(no)+'.below south west),')
                    at_option=NoEscape('at=(plot'+str(no)+'.south west),')
        
        if extra_commands is not None:
            tikzpicture.append(extra_commands)

        return tikzpicture

    def to_pylatex_plot(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.9\textwidth'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False,extra_commands=None):
        
        
        tikz_pic=self._pylatex_tikz(filename,labels_list,colors_list,height, width,x_axis_description,y_axis_description,subplots,extra_commands=extra_commands)
        
        return tikz_pic
    
    def to_standalone_plot(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.5\textwidth'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False,legend_pos='north east',extra_commands=None,options=None):
    
        geometry_options ={"margin": "0cm",}
        doc = Document(documentclass='standalone',geometry_options=None,document_options=["tikz"])
        #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))

        doc.packages.append(Package('siunitx'))
        doc.packages.append(Package('amsmath'))
        doc.packages.append(Package('float'))
        doc.packages.append(Package('tikz'))        
        doc.packages.append(Package('pgfplots'))        
        doc.append(Command('usepgfplotslibrary',arguments='units'))
        doc.append(Command('usetikzlibrary',arguments='spy'))
        
        doc.append(Command('pgfplotsset',arguments=NoEscape(r'compat=newest,label style={font=\small},legend pos='+str(legend_pos))))
        
        tikz_pic=self._pylatex_tikz(filename,labels_list,colors_list,height, width,x_axis_description,y_axis_description,subplots,extra_commands=extra_commands,options=options)
        
        doc.append(tikz_pic)
        
        doc.generate_pdf(filename,clean_tex=False)
        return doc.dumps()

    def to_tikz_plot(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.9\textwidth'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False,legend_pos='north east',extra_commands=None,options=None):


        return self.to_standalone_plot(filename,labels_list,colors_list,height, width,x_axis_description,y_axis_description,subplots,legend_pos,extra_commands=extra_commands,options=options)
    
    def to_standalone_figure(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.9\textwidth'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False,legend_pos='north east',extra_commands=None,options=None):


        self.to_standalone_plot(filename,labels_list,colors_list,height, width,x_axis_description,y_axis_description,subplots,legend_pos,extra_commands=extra_commands,options=options)
        fig = Figure(position='H')
        fig.add_image(filename,width=width)
        
        return fig
    
    
class SpectralMethods(DataMethods):

    def is_uniformly_distributed(self):

        sample_length=max(self.index)-min(self.index)/len(self.index)-1
        step_list=[self.index[i+1] - self.index[i] for i in range(len(self.index)-1)]

        if all(  np.round(element -  step_list[0],10) == 0 for element in step_list):
            return True
        else:
            return False

#     def spectrum(self):

#     spectrum={name:fft.fftshift(fft.fft(data)) for name,data in self.items()}
#     f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#     return SpectrumFrame(data=spectrum,index=fft.fftshift(f_span))

    def to_time_domain(self):

        sampling_rate=self.index[1]-self.index[0]
        t_span=fft.fftshift(fft.fftfreq(len(self.index),d=sampling_rate))
        t_span=t_span[-1]+t_span

        timeseries={name:fft.ifftshift(fft.ifft( data ) )for name,data in self.items()}

        return TimeDataFrame(data=(timeseries),index=t_span)

    def double_sided_rms(self):

        f_span_shifted_ds=(fft.fftshift(self.index))
        spectrum_shifted_ds=fft.fftshift(abs(self)/len(self))

        return SpectrumSeries(data=spectrum_shifted_ds,index=f_span_shifted_ds,name=self.name)

#     def double_sided_spec(self):

#         f_span_shifted_ds=fft.fftshift(self.index)
#         spectrum_shifted_ds={name:fft.fftshift(abs(data)/len(data)) for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

#     def single_sided_spec(self):

#         f_span_shifted_ss=np.positive(fft.fftshift(self.index))
#         spectrum_shifted_ss={name:fft.fftshift(abs(data)/len(data))*np.heaviside([self.index],1)*2 for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)

class EntryWithUnit:
    _units={}
    _latex_backend=latex

    
    @classmethod
    def set_default_units(cls, units={}):

        cls._units = units
        return cls
    
    def __new__(cls,obj,units=None,latex_backend=None,**kwargs):
        obj_with_unit= super().__new__(cls)

        if units is not None:
            obj_with_unit._units= units
        else:
            obj_with_unit._units= cls._units

        if isinstance(obj,Relational):
            obj_with_unit._left_par = ''
            obj_with_unit._right_par = ''
            obj_with_unit._quantity = obj_with_unit._obj.lhs
        else:
            obj_with_unit._quantity=obj_with_unit._obj
        
        has_unit = obj_with_unit._set_quantity_unit()

        if has_unit:

            obj_with_unit._obj =obj
            obj_with_unit._unit=None
            obj_with_unit._left_par = '['
            obj_with_unit._right_par = ']'
            
                
            if latex_backend is not None:
                obj_with_unit._latex_backend= latex_backend
            else:
                obj_with_unit._latex_backend= cls._latex_backend

            return obj_with_unit

        else:
            return obj


        
    def _set_quantity_unit(self):
        if self._quantity in self._units:
            self._unit = self._units[self._quantity]
            return True
        else:
            self._unit=None
            return False


            
    def __str__(self):
        entry_str=self._obj.__str__()
        unit = self._unit
        left_par=self._left_par
        right_par=self._right_par
        
        if unit:
            return f'{entry_str} {left_par}{unit.__str__()}{right_par}'
        else:
            return f'{entry_str}'

    def __repr__(self):
        entry_str=self._obj.__repr__()
        unit = self._unit
        left_par=self._left_par
        right_par=self._right_par
        
        if unit:
            return f'{entry_str} {left_par}{unit.__repr__()}{right_par}'
        else:
            return f'{entry_str}'
        
    def _latex(self,*args):
        entry_str=self._latex_backend(self._obj)
        unit = self._unit
        left_par=self._left_par
        right_par=self._right_par
        
        
        if unit:
            return f'{entry_str} {left_par}{unit:~L}{right_par}'
        else:
            return f'{self._obj}'


class BasicFormattingTools:
    
    
    _latex_backend=vlatex
    _label_formatter = lambda obj: f'${vlatex(obj)}$' if isinstance(obj,(Expr,Eq,EntryWithUnit)) else obj
    _unit_selector = EntryWithUnit
    _domain=None
    _units={}
    _applying_func = lambda x: x
    _init_ops=True

    _default_sep = ', '
    
    @classmethod
    def set_default_column_separator(cls, sep=', '):
        cls._default_sep = sep

        return cls

    @classmethod
    def set_default_units(cls, units={}):

        cls._units = units
        return cls

    @classmethod
    def set_default_unit_selector(cls, selector=EntryWithUnit):

        cls._unit_selector = selector
        return cls
    
    def set_multiindex_axis(self,axis=0):
        
        if axis == 'index':
            axis = 0
        elif axis == 'colums':
            axis = 1
            
        
        idx = self.axes[axis]
        
        if isinstance(idx,pd.MultiIndex):
       
            new_obj = self.copy()
            
        else:
            midx = pd.MultiIndex.from_tuples(idx)
            new_obj = self.copy().set_axis(midx,axis=axis)

        return new_obj

    
    def set_multiindex_columns(self):
        
        return self.set_multiindex_axis(axis=1)
    
    def set_flat_index_axis(self,axis=0):
        
        if axis == 'index':
            axis = 0
        elif axis == 'colums':
            axis = 1

        idx = self.axes[axis]
        
        if isinstance(idx,pd.MultiIndex):
            midx = idx.to_flat_index()
            new_obj = self.copy().set_axis(midx,axis=axis)

        else:

            new_obj = self.copy()

        return new_obj
    
    def set_flat_index_columns(self):
        
        return self.set_flat_index_axis(axis=1)    

    def switch_axis_type(self,axis=0):
        
        if axis == 'index':
            axis = 0
        elif axis == 'colums':
            axis = 1
            
        
        idx = self.axes[axis]
        
        if isinstance(idx,pd.MultiIndex):
            
            new_obj = self.set_flat_index_axis(axis=axis)
            
        else:
            new_obj = self.copy().set_multiindex_axis(axis=axis)

        return new_obj
    
    
    def switch_index_type(self):

        return self.switch_axis_type(axis=0)


    
    
    
    def applying_method(self,data,func=None,**kwargs):
        if func:
            #print('func is used')
            ops_func = func

        elif self.__class__._applying_func is not None:
            print('class func is used')

            ops_func = self.__class__._applying_func
        else:
            #print('identity is used')
            ops_func = lambda data: data

        
        return ops_func(data)

    def _modify_axis(self,func,axis=0):
        
        new_obj = self.copy()
        new_obj_idx= new_obj.axes[axis]
        idx_frame = new_obj_idx.to_frame().applymap(func)
        
        if isinstance(new_obj_idx,pd.MultiIndex):
            new_obj_idx = pd.MultiIndex.from_frame(idx_frame)
        else:
            new_obj_idx = pd.Index(idx_frame,name=new_obj_idx.name)
            


        return new_obj.set_axis(new_obj_idx ,axis=axis)
    

    def format_axis_names(self,formatter=None,axis=0):
        if formatter is None:
            
            formatter = self.__class__._label_formatter

        new_frame=self._modify_axis(formatter,axis=axis)
        
        return new_frame

    def format_index_names(self,formatter=None):

        return self.format_axis_names(formatter=formatter,axis=0)

    def format_columns_names(self,formatter=None):

        return self.format_axis_names(formatter=formatter,axis=1)

    def format_axes_names(self,formatter=None):

        return self.format_axis_names(formatter=formatter,axis=1).format_axis_names(formatter=formatter,axis=0)
    
    def fit_units_to_axis(self,unit_selector=None,axis=0):
        if unit_selector is None:
            unit_selector = self.__class__._unit_selector
            
        new_frame=self._modify_axis(unit_selector.set_default_units(self.__class__._units),axis=axis)
        
        return new_frame

    def fit_units_to_index(self,unit_selector=None):        
        return self.fit_units_to_axis(unit_selector=unit_selector,axis=0)


    def fit_units_to_columns(self,unit_selector=None):        
        return self.fit_units_to_axis(unit_selector=unit_selector,axis=1)

    def fit_units_to_axes(self,unit_selector=None):        
        return self.fit_units_to_axis(unit_selector=unit_selector,axis=1).fit_units_to_axis(unit_selector=unit_selector,axis=0)

    def switch_columns_type(self):


        return self.switch_axis_type(axis=1)
    
    def _format_entry(self, obj,formatter=None):
        if formatter is not None:
            return formatter(obj)
        else:
            return self.__class__._label_formatter(obj)




class AdaptableSeries(pd.Series,BasicFormattingTools):

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





class AdaptableDataFrame(pd.DataFrame,BasicFormattingTools):

    @property
    def _constructor(self):
        return AdaptableDataFrame

    @property
    def _constructor_sliced(self):
        return AdaptableSeries


    @classmethod
    def _init_with_ops(cls,data=None, index=None, columns=None, dtype=None, copy=None,**kwargs):

        raw_frame= cls(data=data, index=index, columns=columns, dtype=dtype, copy=copy)
        new_frame=raw_frame.applying_method(raw_frame,**kwargs)

       
        return cls(data=new_frame, index=index, columns=columns, dtype=dtype, copy=copy)

    @classmethod
    def formatted(cls,data=None, index=None, columns=None, dtype=None, copy=None):
        return cls._init_with_ops(data=data, index=index, columns=columns, dtype=dtype, copy=copy)
    
    

class LatexDataFrame(AdaptableDataFrame
                                 ,BasicFormattingTools):
    _applying_func = lambda obj: (obj).set_multiindex_columns().fit_units_to_axes().format_axes_names().set_multiindex_columns()

    
    @property
    def _constructor_sliced(self):
        return LatexSeries
    @property
    def _constructor(self):
        return LatexDataFrame  

class LatexSeries(AdaptableSeries
                              ,BasicFormattingTools):
    @property
    def _constructor(self):
        return LatexSeries
    @property
    def _constructor_expanddim(self):
        return LatexDataFrame 



class TimeDomainMethods(DataMethods):

    
    def _set_comp_time(self,time):
        self._comp_time=time
        return None
    
    def _get_comp_time(self):
        
        try:
            obj = self._comp_time
        except:
            obj=None
        else:
            obj = self._comp_time
            
        return obj
        
    def gradient(self):
        
        
 
        data_gradient=np.gradient(self.to_numpy(),(self.index) )

        return TimeSeries( data=data_gradient,index=self.index,name=self.name  )


    def is_uniformly_distributed(self):

        sample_length=max(self.index)-min(self.index)/len(self.index)-1
        step_list=[self.index[i+1] - self.index[i] for i in range(len(self.index)-1)]

        if all(  np.round(element -  step_list[0],10) == 0 for element in step_list):
            return True
        else:
            return False

    def to_frequency_domain(self):

        spectrum=(fft.fft(self.to_numpy()))

        f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])
        
        spectrum = SpectrumSeries(data=spectrum,index=(f_span),name=self.name)
        spectrum.index.name='f'
        
        return spectrum

#     def to_frequency_domain(self):

#         spectrum={name:fft.fft(data) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#         return SpectrumFrame(data=spectrum,index=f_span)






class SpectrumSeries(AdaptableSeries,SpectralMethods):

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

class SpectrumFrame(AdaptableDataFrame,SpectralMethods):

    @property
    def _constructor(self):
        return SpectrumFrame


    @property
    def _constructor_sliced(self):
        return SpectrumSeries

    
    def to_tikz_plot(self,filename,labels_list=None,colors_list=['red','blue','orange'],x_axis_description=',xlabel={$f$},x unit=\si{\hertz},',y_axis_description=None,legend_pos='north east',extra_commands=None):
        
        if y_axis_description == None:
            y_axis_description='ylabel=$'+self.name+'$,'
        
        

        return super().to_tikz_plot(filename=filename,labels_list=labels_list,colors_list=colors_list,x_axis_description=x_axis_description,y_axis_description=y_axis_description,legend_pos=legend_pos,extra_commands=extra_commands)

    def to_standalone_figure(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'7cm'), width=NoEscape(r'0.9\textwidth'),x_axis_description=',xlabel={$f$},x unit=\si{\hertz},',y_axis_description='',subplots=False,legend_pos='north east',extra_commands=None,options=None):

        return super().to_standalone_figure(filename=filename,labels_list=labels_list,colors_list=colors_list,x_axis_description=x_axis_description,y_axis_description=y_axis_description,legend_pos=legend_pos,extra_commands=extra_commands,options=options)
    

    def double_sided_rms(self):

        spectrum_shifted_ds={name:data.double_sided_rms() for name,data in self.items()}
        f_span_shifted_ds=(fft.fftshift(self.index))

        return SpectrumFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

    def single_sided_rms(self):

        spectrum_shifted_ss={name:data.double_sided_rms()*np.heaviside(data.double_sided_rms().index,0.5)*2 for name,data in self.items()}
        f_span_shifted_ss=fft.fftshift(self.index)

        return SpectrumFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)

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




class TimeSeries(AdaptableSeries,TimeDomainMethods):

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

    def to_tikz_plot(self,filename,labels_list=None,colors_list=['red','blue','orange'],x_axis_description='xlabel={$t$},x unit=\si{\second},',y_axis_description=None,extra_commands=None,options=None):
        
        if y_axis_description == None:
            y_axis_description='ylabel=$'+self.name+'$,'
        
        
        aux_DataFrame=TimeDataFrame(data=self,index=self.index)
        dumped_tex=aux_DataFrame.to_tikz_plot(filename=filename,labels_list=labels_list,colors_list=colors_list,x_axis_description=x_axis_description,y_axis_description=y_axis_description,extra_commands=extra_commands,options=options)
        return dumped_tex

class TimeDataFrame(AdaptableDataFrame,TimeDomainMethods):

    @property
    def _constructor(self):
        return TimeDataFrame


    @property
    def _constructor_sliced(self):
        return TimeSeries

    def to_frequency_domain(self):

        spectral_data={name:data.to_frequency_domain() for name,data in self.items()}
        f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

        return SpectrumFrame(data=spectral_data,index=spectral_data[next(iter(spectral_data))].index )

    def gradient(self):


        data_gradient={name:data.gradient() for name,data in self.items()}
        return TimeDataFrame(data=data_gradient,index=data_gradient[next(iter(data_gradient))].index )


