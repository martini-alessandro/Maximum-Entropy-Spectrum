#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:17:55 2021

@author: alessandro
"""
import sys
sys.path.insert(0, '..')
sys.path.insert(0, '../memspectrum')
from memspectrum import MESA 
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
from style_sheet import init_plotting
from scipy.interpolate import interp1d 
import matplotlib.ticker as ticker
import tqdm
import time
import os 


#Decide if loading (load = True) or simulating (load = False) 
load_noise = False 
load_estimates = True 
simulate_noise = False
save = False
plot = True 
save_plots = False 

#Initialise directories 
load_dir = 'psd_comparison/LIGO_FPE_VM/'
save_dir = 'psd_comparison/LIGO_FPE_VM/'
if save: 
    save_exists_bool = os.path.exists(save_dir)
    if not save_exists_bool:
        raise FileExistsError("Saving directory doesn't exist. Change saving directory and try again")
noise_file = 'psd_comparison/noise.npy'
psd_file = 'psd_comparison/psd.txt'
plot_savedir = None

#Initialise analysis methods and MESA object
methods = ['FPE', 'VM']
M = MESA() 

#Initialise plot infos 
colors = {'FPE': 'k', 'VM': 'k'}
fmts = {'FPE': '.', 'VM': 'o'}
linestyles = {'FPE': '--', 'VM': 'dashdot'}
alphas = {'FPE': .9, 'VM': 1}
errbar_colors = {'FPE': 'gray', 'VM': 'dimgray'}
lws = {'FPE': 1.3, 'VM': 1.3}
errbar_alphas = {'FPE': 0, 'VM': 0}

#Initialise simulation infos 
frequency, spectrum, _ = np.loadtxt('GWTC1_GW150914_PSDs.dat', unpack = True)
number_of_simulations = 1000
dt = 1./2048
number_of_points = 40960  
T = number_of_points * dt


#Initialize noise and spectrum arrays
spectra = {}
noise = [] 
for method in methods: 
    spectra[method] = []

#If load noise is True, load the noise, otherwhise generates noise vectors with
#characteristics given by #Initialise simulation infos
if load_noise:
    print('Loading Noise relisations\n')
    noise = np.load(load_dir + 'noise.npy')
    
if simulate_noise:     
    print('Simulating data\n')
    for i in tqdm.tqdm(range(number_of_simulations)):
        time.sleep(0.01)
        t, time_series, freq, frequency_series, psd =\
        generate_data(frequency, spectrum, T, 1 / dt)
        noise.append(time_series)

#Load psd from Load directory and the spectra with given methods list 
if load_estimates: 
    spectra = {}
    for method in methods: 
        print('Loading {} spectra\n'.format(method))
        freq, psd = np.loadtxt(load_dir + 'psd.txt')
        spectra[method] = np.load(load_dir + method + '_spectra.npy')
    
else: 
    print('Solving MESA for noise vectors\n')
    for i in tqdm.tqdm(range(number_of_simulations)):
        for method in methods:
            p, a_k, opt = M.solve(noise[i], optimisation_method = method, early_stop = False)
            f, sp = M.spectrum(dt)
            spectra[method].append(sp[: number_of_points // 2 + 1]) 
    spectra[method] = np.array(spectra[method])
#Saves noise, psd and estimats in save_dir 
if save:
    print('Saving noise and psd')
    np.save(save_dir + 'noise.npy', noise)
    np.savetxt(save_dir + 'psd.txt', np.array([freq, psd]))
    print('Saving PSD estimates')
    for method in methods: 
        np.save(save_dir + method + '_spectra.npy', spectra[method])
     
#Plots estimates with plot colors and styles given above. 
if plot: 
    print('Plotting')
    mean = {}
    median = {}
    lower_percentiles, lp = {}, 5
    upper_percentiles, up = {}, 95
    for method in methods: 
        mean[method] = spectra[method].mean(0)
        median[method] = np.median(spectra[method], axis = 0)
        lower_percentiles[method], upper_percentiles[method] =\
            np.percentile(spectra[method], (lp, up), axis = 0)
            
    fig, ax = plt.subplots() 
    axins = ax.inset_axes([.1,.65,.35,.3]) 
    axins2 = ax.inset_axes([.55,.5, .2, .45])

    for method in methods: 
        ax.loglog(freq,mean[method], color=colors[method], alpha=alphas[method],\
                  linestyle = linestyles[method], label = method, lw = lws[method])
        axins.loglog(freq, mean[method], color = colors[method], alpha = alphas[method],\
                     linestyle = linestyles[method], label = method, lw = lws[method])
        axins2.loglog(freq, mean[method], color = colors[method], alpha = alphas[method],\
                     linestyle = linestyles[method], label = method, lw = lws[method])
    ax.loglog(freq, psd, color = 'r', alpha = .5, label = 'True PSD', lw = 2)
    axins.loglog(freq, psd, color = 'r', alpha = .5, lw = 2)
    axins2.loglog(freq, psd, color = 'r', alpha = .5, lw = 2)
    
    axins.set_xticks([], minor = True), axins2.set_xticks([], minor = True)
    # axins.set_xticklabels([]), axins2.set_xticklabels([])
    axins.set_yticks([], minor = True), axins.set_yticklabels([])
    axins2.set_yticks([], minor = True), axins2.set_yticklabels([])

    #Set plot limits 
    axins.set_xlim(30.1, 45.4)
    axins.set_ylim(1.9e-46, 6.1e-42)
    axins2.set_xlim(495, 515.6)
    axins2.set_ylim(5.1e-47, 3.7e-42)
    ax.set_xlim(26,1024)
    ax.set_ylim(1.25e-47,1.26e-38)
    
    
    ax.indicate_inset_zoom(axins, edgecolor="black")
    ax.indicate_inset_zoom(axins2, edgecolor="black")
    ax.set_xlabel('Frequency [Hz]', fontsize = 25)
    ax.set_ylabel('Power Spectral Density $[Hz^{-1/2}]$', fontsize = 25)
    fig.legend(bbox_to_anchor=(0.47, 0.45, 0.5, 0.5), fontsize = 20) 
