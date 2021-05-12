from pandas import (Series,DataFrame)
from numpy import (fft)
import numpy as np
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure , Alignat, Package,Quantity, Command
from pylatex.utils import italic,NoEscape
from sympy.physics.mechanics import vlatex


default_colors=['red','blue','orange','teal','black','green']

class DataMethods:
    def to_tikz_plot(self,filename,labels_list=None,colors_list=default_colors,height=NoEscape(r'6cm'), width=NoEscape(r'15cm'),x_axis_description=',xlabel={$t$},x unit=\si{\second},',y_axis_description='',subplots=False):

        
        geometry_options ={"margin": "0cm",}
        doc = Document(documentclass='standalone',geometry_options=None,document_options=["tikz"])
        #doc=Document(documentclass='subfiles',document_options=NoEscape('bch_r4.tex'))
        doc.append(Command('usepgfplotslibrary',arguments='units'))
        doc.packages.append(Package('siunitx'))
        doc.packages.append(Package('amsmath'))
        doc.packages.append(Package('float'))

        
        labels_list=['$'+vlatex(label)+'$' for label in self.columns]

        data_for_plot_list=[zip(self.index,plot_data)  for label,plot_data in list(self.items())]

        plots_no=len(data_for_plot_list)
        
        colors_multiplicator=np.ceil(plots_no/len(colors_list))

        plot_options = NoEscape('anchor=north west,ymajorgrids=true,xmajorgrids=true,grid style=dashed,legend style={font=\small},'+y_axis_description)+NoEscape(',height=')+height+NoEscape(',width=')+width
        
        #with doc.create(Figure(position='!htb')) as fig:

        plot_options_list=[plot_options+NoEscape(',xticklabels=\empty,')]*(plots_no-1)
        plot_options_list.append( plot_options+NoEscape(x_axis_description) )
        
        with doc.create(TikZ()) as tikzpicture:
            if subplots==False:

                with tikzpicture.create(Axis(options=plot_options_list[-1])) as plot:

                    for data_for_plot,label,color in zip(data_for_plot_list,labels_list,colors_list*int(colors_multiplicator)) :
                        coordinates = data_for_plot

                        plot.append(Plot(name=NoEscape(NoEscape(r'\tiny '+label)), coordinates=coordinates,options='color='+color+',solid'))

            else:
                at_option=NoEscape('')
                for no,combined_plot_data in enumerate(zip(data_for_plot_list,labels_list,colors_list*int(colors_multiplicator))):
                    data_for_plot,label,color = combined_plot_data
                    plot_name=NoEscape(',name=plot'+str(no)+',')
                    with tikzpicture.create(Axis(options=plot_options_list[no]+plot_name+at_option )) as plot:
                        coordinates = data_for_plot
                        plot.append(Plot(name=NoEscape(NoEscape(r'\tiny '+label)), coordinates=coordinates,options='color='+color+',solid'))
                        #at_option=NoEscape('at=(plot'+str(no)+'.below south west),')
                        at_option=NoEscape('at=(plot'+str(no)+'.south west),')
                        

        doc.generate_tex(filename)
        return doc.dumps()


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

        return TimeSeries(data=spectrum_shifted_ds,index=f_span_shifted_ds,name=self.name)

#     def double_sided_spec(self):

#         f_span_shifted_ds=fft.fftshift(self.index)
#         spectrum_shifted_ds={name:fft.fftshift(abs(data)/len(data)) for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

#     def single_sided_spec(self):

#         f_span_shifted_ss=np.positive(fft.fftshift(self.index))
#         spectrum_shifted_ss={name:fft.fftshift(abs(data)/len(data))*np.heaviside([self.index],1)*2 for name,data in self.items()}

#         return TimeDataFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)

class TimeDomainMethods(DataMethods):

    def gradient(self):
        
        data_gradient={name:np.gradient(data,self.index ) for name,data in self.items()}

        return type(self)( data=data_gradient,index=self.index,name=self.name  )

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
        return SpectrumSeries(data=spectrum,index=(f_span),name=self.name)

#     def to_frequency_domain(self):

#         spectrum={name:fft.fft(data) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#         return SpectrumFrame(data=spectrum,index=f_span)


class SpectrumSeries(Series,SpectralMethods):

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

class SpectrumFrame(DataFrame,SpectralMethods):

    @property
    def _constructor(self):
        return SpectrumFrame

    @property
    def _constructor_expanddim(self):
        return SpectrumFrame

    @property
    def _constructor_sliced(self):
        return SpectrumSeries


    def double_sided_rms(self):

        spectrum_shifted_ds={name:data.double_sided_rms() for name,data in self.items()}
        f_span_shifted_ds=(fft.fftshift(self.index))

        return TimeDataFrame(data=spectrum_shifted_ds,index=f_span_shifted_ds)

    def single_sided_rms(self):

        spectrum_shifted_ss={name:data.double_sided_rms()*np.heaviside(data.double_sided_rms().index,0.5)*2 for name,data in self.items()}
        f_span_shifted_ss=fft.fftshift(self.index)

        return TimeDataFrame(data=spectrum_shifted_ss,index=f_span_shifted_ss)

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


class TimeSeries(Series,TimeDomainMethods):

    @property
    def _constructor(self):
        return TimeSeries

    @property
    def _constructor_expanddim(self):
        return TimeDataFrame

    @property
    def _constructor_sliced(self):
        return TimeSeries

#     def spectrum(self):

#         spectrum={name:(fft.fft(data)) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])
#         return SpectrumSeries(data=spectrum,index=(f_span))

    def to_tikz_plot(self,filename,labels_list=None,colors_list=['red','blue','orange'],x_axis_description='xlabel={$t$},x unit=\si{\second},',y_axis_description=None):
        
        if y_axis_description == None:
            y_axis_description='ylabel=$'+self.name+'$,'
        
        
        aux_dataframe=TimeDataFrame(data=self,index=self.index)
        dumped_tex=aux_dataframe.to_tikz_plot(filename=filename,labels_list=labels_list,colors_list=colors_list,x_axis_description=x_axis_description,y_axis_description=y_axis_description)
        return dumped_tex

class TimeDataFrame(DataFrame,TimeDomainMethods):

    @property
    def _constructor(self):
        return TimeDataFrame

    @property
    def _constructor_expanddim(self):
        return TimeDataFrame

    @property
    def _constructor_sliced(self):
        return TimeSeries

    def to_frequency_domain(self):

        spectral_data={name:data.to_frequency_domain() for name,data in self.items()}
        f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

        return SpectrumFrame(data=spectral_data,index=spectral_data[next(iter(spectral_data))].index )

#     def spectrum(self):

#         spectrum={name:fft.fftshift(fft.fft(data)) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#         return SpectrumFrame(data=spectrum,index=fft.fftshift(f_span))

#     def spectrum(self):

#         spectrum={name:fft.fft(data) for name,data in self.items()}
#         f_span=fft.fftfreq(len(self.index),d=self.index[1]-self.index[0])

#         return SpectrumFrame(data=spectrum,index=f_span)

#     def is_uniformly_distributed(self):

#         sample_length=max(self.index)-min(self.index)/len(self.index)-1
#         step_list=[self.index[i+1] - self.index[i] for i in range(len(self.index)-1)]

#         if all(  np.round(element -  step_list[0],10) == 0 for element in step_list):
#             return True
#         else:
#             return False
