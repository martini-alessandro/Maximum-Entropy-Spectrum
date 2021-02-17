#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:56:33 2021

@author: alessandro
"""

import sys
sys.path.insert(0,'..')
from memspectrum import MESA
import welch
import numpy as np
import matplotlib.pyplot as plt
from init_plotting import init_plotting
from GenerateTimeSeries import generate_data
from scipy.interpolate import interp1d 



def compute_spectra(time_series, dt):
    welch_frequency, welch_spectrum = [], []
    MESA_frequency, MESA_spectrum = [], [] 
    for i in range(len(time_series)): 
        None 

def relative_error(real, estimate):
    return np.abs(real - estimate) / real 



if __name__ == '__main__': 
    #Sampling variables 
    dt = 1. / 4096
    times = np.array([1, 5, 10, 100, 1000])
    N = times // dt 
    segment_duration = np.array([521, 1024, 2048, 8192, 32768]) 
    Ny = 0.5 / dt 
    
    #Generate normal PSD and import Ligo PSD 
    wnoise_frequency = np.arange(0, N[0] // 2 + 1) / times[0]
    wnoise_spectrum = np.repeat(1, wnoise_frequency.size) * dt 
    ligo_frequency, ligo_spectrum = np.loadtxt('LIGO-P1200087-v18-AdV_DESIGN_psd.dat', unpack=True)     
    ligo_interp = interp1d(ligo_frequency, ligo_spectrum)
    #Lists containing all simulated arrays, frequencies and psds 
    print('Generating white noise series')
    Wtime, Wtime_series, Wfrequency, Wfrequency_series, Wpsd = generate_data(wnoise_frequency,
                                                                        wnoise_spectrum,
                                                                        times[-1],
                                                                        1/ dt)
    print('Generating Ligo noise series')
    Ltime, Ltime_series, Lfrequency, Lfrequency_series, Lpsd = generate_data(ligo_frequency,
                                                                             ligo_spectrum,
                                                                             times[-1],
                                                                             1 / dt)
    
    
 