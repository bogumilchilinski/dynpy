from pandas import (Series,DataFrame)
from numpy import (fft)
import numpy as np

class SpectralMethods:

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

class TimeDomainMethods:

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
