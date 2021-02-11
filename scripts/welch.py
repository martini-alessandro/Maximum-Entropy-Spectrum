# -*- coding: utf-8 -*-
from scipy.signal import tukey, welch

def psd(strain,
        sampling_rate,
        segment_duration,
        window_function  = None,
        overlap_fraction = 0.5):
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
                             return_onesided = True,
                             scaling         = 'density')
                             
    return frequencies, psd
