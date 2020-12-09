# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:50:17 2020

@author: Workplace
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def generate_data(f,
                  psd,
                  T= 16.0,
                  starttime = 0,
                  sampling_rate = 100000.,
                  fmin = None,
                  fmax = None,
                  zero_noise = False,
                  asd = False):
    # f, psd = np.loadtxt(psd_file, unpack=True)
    if asd is True : psd *= psd
    # generate an interpolant for the PSD
    psd_int = interp1d(f, psd, bounds_error=False, fill_value= 'extrapolate')
    df      = 1 / T
    N       = int(sampling_rate * T)
    times   = np.linspace(starttime, starttime + T , N) 
    if fmin == None: fmin = 0
    if fmax == None: fmax = (N / 2) / T
    # filter out the bad bits
    kmin = np.int(fmin/df)
    kmax = np.int(fmax/df) + 1
    

    # generate the FD noise
    frequencies      = df * np.linspace(kmin, kmax, int(N / 2 + 1)) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
    frequency_series = np.zeros(len(frequencies), dtype = np.complex128)

    if zero_noise is False:
        sigma = np.sqrt(psd_int(frequencies) /  df * .5) 
        phi = np.random.uniform(0, 2 * np.pi, len(sigma))
        frequency_series = sigma * (np.random.normal(0, 1, len(sigma)) + 1j * np.random.normal(0, 1, len(sigma)))
        
      
    # inverse FFT to return the TD strain
    time_series = np.fft.irfft(frequency_series) * df * N 
    return times, time_series, frequencies, frequency_series, psd_int(frequencies)

    