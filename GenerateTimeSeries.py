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
                  sampling_rate = 1.,
                  fmin = None,
                  fmax = None,
                  zero_noise = False,
                  asd = False):
    # f, psd = np.loadtxt(psd_file, unpack=True)
    if asd is True : psd *= psd
    # generate an interpolant for the PSD
    #psd_int = interp1d(f, psd, bounds_error=False, fill_value= 'extrapolate')
    psd_int = interp1d(f, psd, bounds_error=False, fill_value= 0.) #proposed change
    df      = 1 / T
    N       = int(sampling_rate * T)
    times   = np.linspace(starttime, starttime + T , N) 
    if fmin == None: fmin = 0
    if fmax == None: fmax = (N / 2) / T
    # filter out the bad bits
    kmin = np.int(fmin/df)
    kmax = np.int(fmax/df) + 1
    

    # generate the FD noise
    frequencies = df * np.linspace(kmin, kmax, int(N / 2 + 1)) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
    frequency_series = np.zeros(len(frequencies), dtype = np.complex128)

    if zero_noise is False:
        sigma = np.sqrt(psd_int(frequencies) /  df * .5) 
        phi = np.random.uniform(0, 2 * np.pi, len(sigma))
        frequency_series = sigma * (np.random.normal(0, 1, len(sigma)) + 1j * np.random.normal(0, 1, len(sigma)))
      
    # inverse FFT to return the TD strain
    time_series = np.fft.irfft(frequency_series) * df * N 
    return times, time_series, frequencies, frequency_series, psd_int(frequencies)

def generate_noise_mesa(mesa_obj,
                  T= 16.0,
                  starttime = 0,
                  sampling_rate = 1.,
                  fmin = None,
                  fmax = None,
                  N_series = 1):
    """
    Generate some noise starting from a mesa object. The noise generated has the same features as the data given in input to the mesa.
        
        params:
            mesa_obj                    Instance of mesa class, to generate the noise from
            T: `float`                  Length (in seconds) of the signal to generate
            sampling_rate: `np.float`   Sampling rate for the time series to generate
            fmin: `float`               Minimum frequency in the signal (if None, is equal to zero)
            fmax: `float`               Maximum frequency in the signal (if None, Nyquist frequency is used: f_Ny = 0.5*sampling_rate)
            N_series `float`            Number of time series to generate
        
        return
            spectrum: `np.ndarray`   PSD of the model, including both positive and negative frequencies (shape (N,))
    """
    df      = 1 / T
    N       = int(sampling_rate * T)
    times   = np.linspace(starttime, starttime + T , N) 
    if fmin == None: fmin = 0
    if fmax == None: fmax = (N / 2) / T
    # filter out the bad bits
    kmin = np.int(fmin/df)
    kmax = np.int(fmax/df) + 1
    
    # generate the FD noise
    frequencies = df * np.linspace(kmin, kmax, int(N / 2 + 1)) #(D,) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
    psd = mesa_obj.spectrum(1/sampling_rate, frequencies)

    sigma = np.sqrt(psd /  df * .5) #(D,)
    phi = np.random.uniform(0, 2 * np.pi, len(sigma))
    frequency_series = np.einsum('ij,j -> ij',np.random.normal(0, 1, (N_series,len(sigma))) + 1j * np.random.normal(0, 1, (N_series,len(sigma)) ), sigma) #(N_series,D)
      
    # inverse FFT to return the TD strain
    time_series = np.fft.irfft(frequency_series) * df * N ##(N_series, N )
    return times, np.squeeze(time_series), frequencies, np.squeeze(frequency_series), psd



    
