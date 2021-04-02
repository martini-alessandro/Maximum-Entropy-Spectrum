# -*- coding: utf-8 -*-
"""
Module that generates a random time series with a given power spectral density
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def generate_data(f,
                  psd,
                  T, 
                  sampling_rate = 1.,
                  fmin = None,
                  fmax = None,
                  asd = False):
    """
    Generate a time series with a given power spectral density 

    Parameters
    ----------
    f : 'np.ndarray'
        The frequencies over which the power spectral density is evaluated. (Shape (N,))
    psd : 'np.ndarray'
        The power spectral density (Shape (N,))
    T : 'np.float'
        The total time of the observation 
    sampling_rate : 'np.float', optional
        The sampling rate of the output series. The default is 1..
    fmin : 'np.float', optional
        The minimum frequency available. The default is None.
    fmax : 'np.float', optional
        Tha maximum frequency available. The default is None.
    asd : 'boolean', optional
        If True, takes the square of the input power spectral density. The default is False.

    Returns
    -------
    times: 'np.ndarray' 
        The sampling time vector. (Shape (N,))
    time_series: 'np.ndarray'
        The output time series  (Shape (N,))
    frequencies: 'np.ndarray'
        The sampling frequencies (Shape (N,))
    frequency_series: 'np.ndarray'
        The output series in frequency domain (Shape (N,))
    psd: 'np.ndarray'
        The frequencies interpolated power spectral density (Shape (N,))
    """
    # f, psd = np.loadtxt(psd_file, unpack=True)
    if asd is True : psd *= psd
    # generate an interpolant for the PSD
    psd_int = interp1d(f, psd, bounds_error=False, fill_value='extrapolate')
    df      = 1 / T
    N       = int(sampling_rate * T)
    times   = np.linspace(0, T, N) 
    if fmin == None: fmin = 0
    if fmax == None: fmax = (N / 2) / T
    # filter out the bad bits
    kmin = np.int(fmin/df)
    kmax = np.int(fmax/df) + 1
    

    # generate the FD noise
    frequencies = df * np.arange(kmin, kmax) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
    frequency_series = np.zeros(len(frequencies), dtype = np.complex128)


    sigma = np.sqrt(psd_int(frequencies) /  df * .5) 
    frequency_series = sigma * (np.random.normal(0, 1, len(sigma)) + 1j * np.random.normal(0, 1, len(sigma)))
      
    # inverse FFT to return the TD strain
    time_series = np.fft.irfft(frequency_series, n=N) * df * N
    return times, time_series, frequencies, frequency_series, psd_int(frequencies)


