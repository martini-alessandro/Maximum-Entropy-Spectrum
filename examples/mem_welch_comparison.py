#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:56:33 2021

@author: alessandro
"""

from memspectrum import MESA
from memspectrum import GenerateTimeSeries
from scipy.signal import welch, tukey
import scipy.stats 
import numpy as np
import matplotlib.pyplot as plt


def welch_estimate(noise, N, dt, overlap_fraction = .5):
    sampling_rate = 1 / dt
    segment_duration = N * dt #lunghezza sottoarray
    padding = 0.4 / segment_duration 
    window_function  = tukey(int(sampling_rate*segment_duration), padding)
    wf, ws = welch(noise, fs  = 1 / dt,
                              window = window_function,
                              nperseg         = N,
                              noverlap        = int(overlap_fraction*N),
                              return_onesided = False,
                              scaling         = 'density',
                              nfft = noise.size) #dopo N // 2 frequ negative
    return wf, ws 


def relative_error(real, estimate):
    return np.abs(real - estimate) / real 


if __name__ == '__main__': 
    
    #Generate the Time Series 
    f_ny = 5 
    dt = 1 / (2 * f_ny)
    N = 8192
    T = N * dt 
    frequency = np.arange(0, N // 2 + 1) / T 
    bnoise_psd = 1 / (frequency ** 2 + 1)
    normal_psd = scipy.stats.norm.pdf(frequency, 2.5, .5) 
    number_of_segments = 5 #Number of segments for Welch method
    
    #Generate Time series 
    time, normal_series, frequency, normal_frequency_series, normal_psd =\
        GenerateTimeSeries.generate_data(frequency, normal_psd, T, 1 / dt)
        
    
    time, bNoise_series, frequency, bnoise_frequency_series, bnoise_psd =\
        GenerateTimeSeries.generate_data(frequency, bnoise_psd, T, 1 / dt)
    
    #Estimate the spectrum for normal PSD 
    M = MESA(normal_series)
    M.solve() 
    normal_ms, normal_mf = M.spectrum(dt)
    normal_wf, normal_ws = welch_estimate(normal_series, 8192 // number_of_segments, dt)
    
    #Estimate the spectrum for brown noise 
    M = MESA(bNoise_series)
    M.solve() 
    bNoise_ms, bNoise_mf = M.spectrum(dt)
    bNoise_wf, bNoise_ws = welch_estimate(bNoise_series, 8192 // number_of_segments, dt)
    
    #ComputeErrors 
    normal_m_error = relative_error(normal_psd[:-1], normal_ms[:N // 2])
    normal_w_error = relative_error(normal_psd[:-1], normal_ws[:N // 2])
    bnoise_m_error = relative_error(bnoise_psd[:-1], bNoise_ms[:N // 2])
    bnoise_w_error = relative_error(bnoise_psd[:-1], bNoise_ws[:N // 2])
    
    #Normal_psd_plot
    fig, ax = plt.subplots(2)
    ax[0].loglog(normal_wf[: N // 2], normal_ws[: N // 2], color = 'green' , label = 'Welch')
    ax[0].loglog(frequency[:-1], normal_psd[:-1], '--', color = 'k')
    ax[0].loglog(normal_mf[: N // 2], normal_ms[: N // 2], color = 'r', label = 'MESA')
    ax[0].set_ylabel('Spectrum [a.u.]', fontsize = 17)
    ax[1].loglog(normal_wf[: N // 2], normal_w_error, '.', color = 'green', label = 'Welch')
    ax[1].loglog(normal_mf[: N // 2], normal_m_error, '.', color = 'r', label = 'MESA')
    ax[1].set_ylabel('Relative Error', fontsize = 17)
    ax[1].set_xlabel('Frequency [a.u]', fontsize = 17)
    ax[0].legend(fontsize = 17)
    ax[1].legend(fontsize = 17)
    
    #brown noise PSD plot 
    fig2, ax2 = plt.subplots(2)
    ax2[0].loglog(bNoise_wf[: N // 2], bNoise_ws[: N // 2], color = 'green', label = 'Welch')
    ax2[0].loglog(frequency[ : -1], bnoise_psd[ : -1], '--', color = 'k')
    ax2[0].loglog(bNoise_mf[: N // 2], bNoise_ms[: N // 2], color = 'r', label = 'MESA')
    ax2[0].set_ylabel('Spectrum [a.u.]', fontsize = 17)
    ax2[1].loglog(bNoise_wf[: N // 2], bnoise_w_error, '.', color = 'green', label = 'Welch')
    ax2[1].loglog(bNoise_mf[: N // 2], bnoise_m_error, '.', color = 'r', label = 'MESA')
    ax2[1].set_xlabel('Frequency [a.u.]', fontsize = 17)
    ax2[1].set_ylabel('Relative Error [a.u.]', fontsize = 17)
    ax2[0].legend(fontsize = 17)
    ax2[1].legend(fontsize = 17)
    