# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 12:53:23 2020

@author: Workplace
"""
import unittest
import numpy as np
import cpnest.model
from scipy import stats
from memspectrum import MESA
import os
from scipy.signal import welch, tukey
import matplotlib.pyplot as plt 

class Data(object):
    
    def __init__(self, x, y): 
        self.x = x 
        self.y = y
        
    def residuals(self, model, *parameters):
        model = Model(self, model)
        return self.y - model(*parameters)
    
    
class SinusoidModel(cpnest.model.Model):
    """
    An n-dimensional gaussian
    """
    def __init__(self, data, dt, model = 'sin2'): 
        self.names = ['Amplitude1', 'Amplitude2', 'Frequency']
        self.bounds = [[5, 10], [12, 18], [5,10]]
        self.data = data
        self.model = model
        self.dt = dt
        self.spectra = [] 
        self.N = len(self.data.x) 
        self.duration = dt * self.N 
        
    def log_likelihood(self, x):
        #Compute Residuals
        A1, A2, F = x['Amplitude1'], x['Amplitude2'], x['Frequency']
        residuals = self.data.residuals(self.model, A1, A2, F)
        
        psd = self.computeSpectrum(residuals) 
        #Compute Likelihood 
        
        Fresiduals =  np.fft.fft(residuals)  
        # Calculate the first term in Eq (10) of Veitch (2015)
        termA = (
            -2.0
            * np.vdot(Fresiduals, Fresiduals / psd)
            / self.duration
        )

        # Calculate the second term in Eq (10) of Veitch (2015)
        termB = -0.5 * np.sum(np.log(0.5 * np.pi * self.duration * psd))


        logL = termA + termB
        return logL
        
    def log_prior(self, p):
        logP = super(SinusoidModel,self).log_prior(p)
#        for n in p.names: logP += -0.5*(p[n])**2
        return logP
    
    def computeSpectrum(self, residuals):
        M = MESA()
        M.solve(residuals, optimisation_method="CAT")
        freq, spectrum = M.spectrum(self.dt)
        self.spectra.append(spectrum)
        return spectrum
    
    

class Model(object):
    
    def __init__(self, data, model):
        """Generate Models from a given set a data, following the chosen distribution.
        For information about order for the parameters check the documentation for 
        every available model. Model is not case sensitive"""
        #Initialize variables
        availModels = ["exponential", "powerlaw", "line", "cauchy", "sin", "sin2",\
                       "gauss"]
        self.data = data
        
        #Models should be in models
        if model.lower() in availModels: 
            self.model = model
        else: 
            raise ValueError("'{}' Model not available. Available models are {}"\
                             .format(model, ', '.join(availModels)))
        
    def __call__(self, *parameters):
        
        if self.model.lower() == 'exponential':
            return self.exponential(parameters[0], parameters[1], parameters[2])
        if self.model.lower() =='powerlaw':
            return self.powerLaw(parameters[0], parameters[1], parameters[2])
        if self.model.lower() == 'line':
            return self.line(parameters[0], parameters[1])
        if self.model.lower() == 'cauchy':
            return self.cauchy(parameters[0], parameters[1], parameters[2], 
                               parameters[3])
        if self.model.lower() == 'sin': 
            return self.sin(parameters[0], parameters[1], parameters[2])
        if self.model.lower() == 'sin2':
            return self.sin2(parameters[0], parameters[1], parameters[2])
        if self.model.lower() == 'gauss':
            return self.gauss(parameters[0], parameters[1], parameters[2])
        
    #Basic models 
    def expoential(self, base, height, floor): 
        """Generate the expected value if datas were exponentially distributed.
        """
        return height * base ** self.data.x + floor
     
    def powerLaw(self, exponent, height, floor):
        """Generate the expected values for power Law distibuted datas"""
        
        return height * (self.data.x ** exponent) + floor
    
    def line(self, slope, intercept):
        """Genereate expected values for linear prediction""" 
        return slope * self.data.x + intercept
    
    def cauchy(self, mode, scale, height, floor):
        """Return the image for cauchy - distributed datas height parameters
        is to be interpreted as height / scale """
        return height  / (1 + ((self.data.x - mode) / scale)) ** 2 + floor
    
    def sin(self, amplitude, frequency, phase):
        return amplitude * np.sin(frequency * self.data.x + phase)
    
    def sin2(self, amplitude1, amplitude2, frequency): 
        return amplitude1 * np.sin(frequency * self.data.x) +\
               amplitude2 * np.cos(frequency * self.data.x)
               
    def gauss(self, amplitude, mean, size): 
        return stats.norm.pdf(self.data.x, mean, size) * amplitude 
        

if __name__== '__main__':
    # amplitude, frequency, phase = 23.6, 7.3, 5.4
    amplitude1, amplitude2, frequency= 8, 14, 7.3
    #Generate Data 
    NumberOfData = 300
    dt = 1 / (4 * frequency)
    x = np.arange(0, NumberOfData) * dt
    T = NumberOfData  * dt
    error = np.random.normal(0, amplitude1, NumberOfData)
    m = MESA() 


    # y = amplitude * np.sin(frequency * x + phase)
    y = amplitude1 * np.sin(frequency * x) +\
        amplitude2 * np.cos(frequency * x) +\
        error
    D = Data(x, y)
    
    #Define Frequency arrays 
    Ny = 1 / (2 * frequency)
  
    #Compute spectrum 
    c = SinusoidModel(D, dt, 'sin2')
    work=cpnest.CPNest(c,   verbose      =  2,
                            poolsize     = 32,
                            nthreads     = 2,
                            nlive        = 1024,
                            maxmcmc      = 500,
                            output       = os.getcwd())

    work.run()
    print('log Evidence {0}'.format(work.NS.logZ))
    x = work.posterior_samples.ravel()
    np.savetxt('RealDatas.txt', y)
    np.savetxt('PosteriorSamples.txt', x)


 
