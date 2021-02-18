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



def relative_error(real, estimate):
    return np.abs(real - estimate) / real 



if __name__ == '__main__': 
    import os 
    save_dir = os.getcwd() + '/alessandro_fake_comparisons/'
    seglen_factor = 1
    save = False 
    white_noise = False 
    Ligo_noise = True 
    
    #Sampling variables 
    dt = 1. / 4096
    Ny = 0.5 / dt 
    times = np.array([1, 5, 10, 100, 1000])
    N = times // dt 
    segment_length = np.array([521, 1024, 2048, 8192, 32768]) 
    segment_length *= seglen_factor

    
    #Generate normal PSD and import Ligo PSD 
    if white_noise: 
        wnoise_frequency, wnoise_spectrum = np.arange(0, N[0] // 2 + 1) / times[0], np.repeat(1, N[0] // 2 + 1) * dt
        print('Generating white noise series')
        Wtime, Wtime_series, Wfrequency, Wfrequency_series, Wpsd = generate_data(wnoise_frequency,
                                                                        wnoise_spectrum,
                                                                        times[-1],
                                                                        1/ dt)
        w_interp = interp1d(wnoise_frequency, wnoise_spectrum) 
        
    #ligo_frequency, ligo_spectrum = np.loadtxt('LIGO-P1200087-v18-AdV_DESIGN_psd.dat', unpack=True)     
    if Ligo_noise: 
        ligo_frequency, ligo_spectrum, _ = np.loadtxt('GWTC1_GW150914_PSDs.dat', unpack = True)
    #Generating Noise 
   
        print('Generating Ligo noise series')
        Ltime, Ltime_series, Lfrequency, Lfrequency_series, Lpsd = generate_data(ligo_frequency,
                                                                             ligo_spectrum,
                                                                             times[-1],
                                                                             1 / dt)
        l_interp = interp1d(ligo_frequency, ligo_spectrum, fill_value = 'extrapolate')    
    
    M = MESA()
    init_plotting()
    
    for i, n in enumerate(segment_length):
        if white_noise: 
            w_welchFreq, w_welchSpectrum = welch.psd(Wtime_series[:int(N[i])], 1 / dt, n * dt)
            M.solve(Wtime_series[:int(N[i])])
            w_mesaSpectrum, w_mesaFreq = M.spectrum(dt)
            
            fig, ax = plt.subplots() 
            ax.loglog(w_welchFreq[ : n // 2], w_welchSpectrum[: n // 2], color = 'blue')
            ax.loglog(w_mesaFreq[ : int(N[i] // 2)], w_mesaSpectrum[: int(N[i] // 2)], color = 'r')
            ax.loglog(wnoise_frequency, wnoise_spectrum, '--', color = 'k')
            ax.set_title('White Noise reconstrcution with {} points'.format(N[i]))
            if save: 
                fig.savefig(save_dir + 'White Noise {} points ({}).pdf'.format(int(N[i]), seglen_factor))
                
        if Ligo_noise: 
            l_welchFreq, l_welchSpectrum = welch.psd(Ltime_series[:int(N[i])], 1 / dt, n * dt )
            M.solve(Ltime_series[:int(N[i])])
            l_mesaSpectrum, l_mesaFreq = M.spectrum(dt)
    
            fig2, ax2 = plt.subplots()
            ax2.loglog(l_welchFreq[: n // 2], l_welchSpectrum[: n // 2], color = 'blue')
            ax2.loglog(l_mesaFreq[: int(N[i] // 2)], l_mesaSpectrum[:int(N[i] // 2)], color = 'red')
            ax2.loglog(ligo_frequency, ligo_spectrum, '--', color = 'k')
            ax2.set_title('Ligo Noise reconstrcution with {} points'.format(N[i]))
            ax2.set_xlim(ligo_frequency.min(), ligo_frequency.max())
            if save: 
                fig2.savefig(save_dir +'GW150914 {} points ({} seglen).pdf'.format(int(N[i]), seglen_factor))
        

    
    
    
 
