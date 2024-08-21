# -*- coding: utf-8 -*-
from scipy.signal import welch
from scipy.signal.windows import tukey

def psd(strain,
        sampling_rate,
        segment_duration,
        window_function  = None,
        overlap_fraction = 0.5,
        nfft = None,
        return_onesided = False):
    """
    Computes the power spectral density of a time series using the Welch method
   
    :param strain: input time series
        :type strain: array
    :param sampling_rate: sampling rate of the input strain (Hz)
        :type sampling_rate: float
    :param segment_duration: length of each segment (in seconds) for the application of the Welch method
        :type segment_duration: float
    :param window_function: window type
        :type window_function: string, tuple or array-like
    :param overlap_fraction: overlap between segments (as a fraction of the segment duration)
        :type overlap_fraction: float
    :param nfft: total number of points for the fft. If None is equal to the length of the dataset
        :type nfft: int
    :param return_onesided: returns onesided (twoside) power spectral density if True (False) 
        :type return_onesided: boolean 
    
    :return: frequencies, corresponding power spectral density
    :rtype: array, array
    """
    
    # if no window has been passed
    # create a default Tukey window
    if window_function is None:
        padding = 0.4/segment_duration
        window_function  = tukey(int(sampling_rate*segment_duration),padding)
    
    # number of samples per segments
    N = int(segment_duration*sampling_rate)
    
    # compute the PSD and its frequencies
    frequencies, psd = welch(strain,
                             fs              = sampling_rate,
                             window          = window_function,
                             nperseg         = N,
                             noverlap        = int(overlap_fraction*N),
                             return_onesided = return_onesided,
                             scaling         = 'density',
                             nfft = nfft)
                             
    return frequencies, psd
