#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:17:55 2021

@author: alessandro
"""
import sys
sys.path.insert(0, '..')
from memspectrum import MESA 
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
from init_plotting import init_plotting
from scipy.interpolate import interp1d 

def relative_error(real, estimate): 
    return np.abs(real - estimate) / real 

PSD = 'ligo' #normal or ligo 
save = False
savedir = '../paper/Images/optimisers_comparison/' + PSD + '/'
simulate_data = False 
methods = ['FPE', 'CAT', 'OBD']
M = MESA() 
init_plotting() 
colors = {'FPE': 'r', 'CAT': 'k', 'OBD': 'green'}

#Initialize single spectrum dictionaries 
spectra, optimisers, orders, errors = {}, {}, {}, {} 
for method in methods: 
    spectra[method], optimisers[method], orders[method], errors[method] = [], [], [], []

#ensemble dictionaries 
median, p5, p95, ensemble_error = {}, {}, {}, {}

#Generating noise 
if PSD.lower() == 'normal': 
        number_of_simulations, dt, number_of_points = 1000, 1./10, 3000 
        f = np.linspace(0, .5 / dt, number_of_points) 
        spectrum = scipy.stats.norm.pdf(f, 2.5, .5)

elif PSD.lower() == 'ligo':  
        number_of_simulations, dt, number_of_points = 500, 1./2048, 40960
        f, spectrum = np.loadtxt('LIGO-P1200087-v18-AdV_DESIGN_psd.dat', unpack = True)
        
T = number_of_points * dt 

if simulate_data: 
    for i in range(number_of_simulations): 
        time, time_series, frequency, frequency_series, psd =\
            generate_data(f, spectrum, T, 1 / dt)
        for method in methods:
            p, a_k, opt = M.solve(time_series, optimisation_method = method, early_stop = False)
            optimisers[method].append(opt)
            orders[method].append(M.a_k.size)
            spectra[method].append(M.spectrum(dt)[0][: number_of_points // 2])    

elif not simulate_data:
    frequency = np.linspace(0, .5 / dt, number_of_points // 2 + 1)
    psd_int = interp1d(f, spectrum, fill_value='extrapolate')
    psd = psd_int(frequency)
    
    if PSD.lower() == 'ligo': 
        s = np.load('ligo_spectra.npy')
        o = np.load('ligo_optimisers.npy')
        orde = np.load('ligo_orders.npy')
        
        for i, method in enumerate(methods):
            spectra[method] = s[i]; optimisers[method] = o[i]; orders[method] = orde[i]
#Set unity measures 
if PSD == 'ligo': 
    y = r"$PSD \left(\frac{1}{Hz} \right)$"

    
if PSD == 'normal': y = r"$PSD \left(a.u.\right)$"

#Plots for comparing methods 
ord_fig, ord_ax = plt.subplots()            #Error VS order for each method
err_fig, err_ax = plt.subplots(1, 3)        #Error hist for each method 
aro_fig, aro_ax = plt.subplots(1, 3)        #AR order estimate hist for each method
fig4, ax4 = plt.subplots()                  #Plot with every spectral estimate    

ord_ax.set_xlabel('Filter Length estimate'); ord_ax.set_ylabel('Frequency averaged error')
ax4.set_ylabel(y); ax4.set_xlabel(r"$f(Hz)$")
err_ax[0].set_ylabel('Percentage')
aro_ax[0].set_ylabel('Percentage')

#Error histogram plot 

for i, method in enumerate(methods):
    err_ax[i].set_xlabel('error ({})'.format(method))
    aro_ax[i].set_xlabel('AR order ({})'.format(method))



for i, method in enumerate(methods):
    #Singl and ensemble properties and errors 
    spectra[method],   optimisers[method]= np.array(spectra[method]), np.array(optimisers[method])
    errors[method] = relative_error(psd[:-1], spectra[method]).mean(1)
    p5[method], median[method], p95[method] = np.percentile(spectra[method],(5, 50, 95), axis = 0)
    ensemble_error[method] = relative_error(psd[:-1], median[method])

    
    #Plot optimizer behaviour
    fig, ax = plt.subplots() 
    ax.plot(optimisers[method].mean(0), '.', color = 'k')
    ax.set_xlabel('Filter length')
    ax.set_ylabel('{} optimizer value'.format(method))
    
    #Plot error VS order 
    fig2, ax2 = plt.subplots() 
    ax2.plot(orders[method], errors[method], '.', color = 'k')
    ax2.set_xlabel('Filter order estimate')
    ax2.set_ylabel('Frequency averaged errors ({})'.format(method))
    
    #plot order in one plot only
    ord_ax.plot(orders[method], errors[method], '.', color = colors[method], label = method)

    #Plot reproduction and frequency error 
    fig3, ax3 = plt.subplots(2)
    ax3[0].loglog(frequency[:-1], median[method], color = 'r')
    ax3[0].loglog(frequency, psd, color = 'k', ls = '--')
    if PSD == 'normal': ax3[0].set_xlim(ax3[0].get_xlim()[0], .5 / dt)
    elif PSD == 'ligo': ax3[0].set_xlim(f.min(), .5 / dt)
    ax3[0].fill_between(frequency[:-1], p5[method], p95[method], color = 'blue', alpha = .5)
    ax3[1].loglog(frequency[:-1], ensemble_error[method], '.', color = 'k')
    ax3[1].set_xlabel(r'$f(Hz)$'); ax3[0].set_ylabel(y); ax3[1].set_ylabel('Percentage Error')
    

    
    #Plot histogram for errors and orders 
    
    err_ax[i].hist(errors[method], bins = 70, density = True, color = 'k', histtype = 'step')
    aro_ax[i].hist(orders[method], bins = 70, density = True, color = 'k', histtype = 'step')

    ax4.loglog(frequency[:-1], median[method], color = colors[method], label = method)
    
    if save: 
        fig.savefig(savedir + '{} optimizer values.pdf'.format(method), bbox_inches = 'tight')
        fig2.savefig(savedir + '{} error VS order.pdf'.format(method), bbox_inches = 'tight')
        fig3.savefig(savedir + '{} spectrum estim.pdf'.format(method), bbox_inches = 'tight')

    print('{}; {}: r_ensemble: {}; sigma_r: {}\n'.format(PSD, method,
                                                     ensemble_error[method].mean(),
                                                     ensemble_error[method].std()))
    
ord_fig.legend(loc = 'lower right')
ax4.loglog(frequency, psd, '--', color = 'k') 
ax4.set_xlim(f.min(), .5 / dt)
fig4.legend(loc = 'upper right')
if save: 
    err_fig.savefig(savedir + 'errors hist.pdf', bbox_inches = 'tight')
    ord_fig.savefig(savedir + 'error VS order comparison.pdf', bbox_inches = 'tight')
    aro_fig.savefig(savedir + 'orders hist.pdf', bbox_inches = 'tight')
    fig4.savefig(savedir + 'compare estimates.pdf', bbox_inches = 'tight')

        
        


      