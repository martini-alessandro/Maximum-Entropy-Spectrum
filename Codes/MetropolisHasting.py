# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:05:38 2020

@author: Alessandro Martini

This class is meant to implement a Markov Chain Monte Carlo algorithm, to 
explore the posterior space parameter and infere posterior distribution shape.
"""
import numpy as np 
import MESAAlgorithm
import sys

class Data(object): 
    def __init__(self, x, y, e): 
        self.x = x 
        self.y = y 
        self.e = e
        
    def generateData(self, model, *parameters):
        m = Model(self, model)
        return m(*parameters)
    
    def residuals(self, model, *parameters):
        return self.x - self.generateData(model, *parameters)

class Prior(object):
        
    def __init__(self, prior):
        """Create a prior object that act on parameters""" 
        #Check that choesn prior is available
        priors = ['uniform','jeffreys']
        if prior.lower() in priors:
            self.prior = prior
        else: 
            raise ValueError("Invalid prior: possible choices are \
                             {}".format("'Uniform' and 'Jeffreys'"))
    def __call__(self, parameter, value = None, log = True): 
            
        #If no value passed, set value to be parameter's value
        if value == None: value = parameter.value
            
        #Compute prior if value is in parameter's domain 
        if parameter.min < value < parameter.max:
            if self.prior.lower() == 'uniform':
                #Returns value for uniform prior
                if log: return 0
                else: return 1
            elif self.prior.lower() == 'jeffreys':
                #Returns value for Jeffreys prior
                if log: return np.log(1 / value)
                else: return 1 / value
            
                
        #Compute prior if value is not in parameter's domain
        else: 
            if log: return -np.inf
            else: return 0
     
class Parameter(object):
    
    def __init__(self, _min, _max, prior, update, value = None):
        """Create a Parameter object whose domain is between specifiend _min
        and _max values, with given prior distribution""" 
        self.min = _min
        self.max = _max
        self.prior = prior
        self.update = update
        #consider adding prior in __init__ . Might be easier to implement in Metropolis, when calling posterior 
        if value == None: self._setValue()
        else: self.value = value
    
    def priorValue(self, value = None, log = True):
        """Compute the value of the Prior for the parameter. If No Value is 
        passed, compute value of prior on current value of the parameter, 
        otherwise, if value is 'proposal', compute prior on the update value.
        If some specific value is passed, compute prior on specified value""" 
        p = Prior(self.prior) 
        
        #If None value is passed, compute prior on current value
        if value == None:
            pass 
        #If Value is proposal, compute prior on the proposed value
        elif str(value).lower() == 'proposal':
            value = self.proposal
        return p(self, value, log)
    
    def _update(self, upd = None):
        """Update the value of the parameter and record it in self.newValue"""
        if upd == None: upd = self.update
        self.proposal = self.value + np.random.uniform(-upd, upd)
        return self.proposal
    
    def _accept(self):
        """Accept proposal and value to be equal to new Value""" 
        self.value = self.proposal
    
    
    def _setValue(self):
        """Compute random value for the parameter between max and min values"""
        self.value = np.random.uniform(self.min, self.max)
        return None 
    
class Model(object):
    
    
    def __init__(self, data, model):
        """Generate Models from a given set a data, following the chosen distribution.
        For information about order for the parameters check the documentation for 
        every available model. Model is not case sensitive"""
        #Initialize variables
        availModels = ["exponential", "powerlaw", "line", "cauchy", "sin"]
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


class Likelihood(object): 
    
    def __init__(self, data, model, distribution, log = True):
         """Create a Likelihood object for the chosen dataset. Model chosen
         to describe datas and functional form of the likelihood are required.
         both model and distribution are not case sensitive""" 
         #Initialize variables 
         availDistributions = ['exp', 'norm']
         self.data = data
         self.model = Model(data, model)
         
         #Log can only be True o False 
         if log in [True, False]: 
             self.log = log
         else:
             raise ValueError('Log can only be True o False')
         
         #Distribution should be of correct Type
         if distribution.lower() in availDistributions:
             self.distribution = distribution
         else: raise ValueError('{} distribution is not available.\
                               Available distributions are.'\
                                   .format(distribution, ', '.join(availDistributions)))
     
    def __call__(self, *parameters):
        """Compute value of the likelihood on data for the chosen model""" 
        
        #Compute Value of Normal distributed likelihood
        if self.distribution.lower() == 'norm':
            logL = -.5 * (self._residuals(*parameters) ** 2).sum()
            if   self.log == True: return logL
            elif self.log == False: return np.exp(logL)
        
        #Compute Value of exponential distributed likelihood
        elif self.distribution.lower() == 'exp':
            logL = - (self.data.y / self.model(*parameters)).sum()
            if   self.log == True: return logL
            elif self.log == False: return np.exp(logL)
            
    def _residuals(self, *parameters):
        """Compute residuals between obtained data and model predictions""" 
        return self.data.y - self.model(*parameters)
    
    
class Posterior(object):
    
    def __init__(self, data, likelihood_distribution, model, log = True):
        
        self.data = data
        self.likelihood = likelihood_distribution
        self.model = model
        if log in [True, False]:
           self.log = log
        else:
            raise ValueError('log can only be True o False')

    
    def value(self, *parameters):
        """Compute Value for the posterior distribution on proposed parameter
        values, for given likelihood and data""" 
        
        #Record Value for parameter and initialize likelihood
        par = self._LikelihoodParameters(0, *parameters)
        l = Likelihood(self.data, self.model, self.likelihood, self.log)

        
        #Compute value of the posterior for current parameter value
        posteriorValue = l(*par)
        for parameter in parameters:
            if type(parameter) == Parameter:
                posteriorValue += parameter.priorValue()
        return posteriorValue
    
    def proposal(self, *parameters):
        """Compute Value for the posterior distribution on proposed parameter 
        values, for given likelihood and data"""
        
        #Record Value for parameter object to compute Likelihood on values
        par = self._LikelihoodParameters(1, *parameters)
        l = Likelihood(self.data, self.model, self.likelihood, self.log)
        
        #Compute value of the posterior for proposed parameters's value
        posteriorValue = l(*par)
        for parameter in parameters:
            if type(parameter) == Parameter:
                posteriorValue += parameter.priorValue('proposal')
        return posteriorValue
        
    def _LikelihoodParameters(self, val, *parameters):
        """From mixed parameters return a tuple containing just values for every
        parameter. Compute current value if parameter Values is set to 'current'
        or 0, or proposed values if set to 'proposal'
        """
        #Create a tuple containing values for the parameters as float objects
        par = []
        for p in parameters:
            if type(p) == Parameter:
                if val == 0: par.append(p.value)
                else:        par.append(p.proposal)
            elif type(p) == float or type(p) == int:
                par.append(p)
        par = tuple(par)
        return par
    

class Metropolis(object): 
    def __init__(self, N, data, likelihood_distribution, model, *parameters):
        self.N = N 
        self.data = data
        self.posterior = Posterior(self.data, likelihood_distribution, model)
        self.parameters = parameters
        self.samples = []
    
    def solve(self): 
        for _ in range(self.N): 
            sys.stdout.write('\r acceptance percentage %f' % ((_ + 1) / self.N))
            #Updata parameters value to compute they proposal value
            self._updateParameters()
            
            #Compute posterior difference
            r = self.posterior.proposal(*self.parameters) -\
                self.posterior.value(*self.parameters)
            
            #accept values if condition met
            if r > np.log(np.random.uniform(0, 1)): 
        
                self._acceptNewValues()
    
            #update samples array
            self._updateSamples()
        
        return np.array(self.samples)
            
    def _updateParameters(self):
        #If parameter object, update its value, otherwise pass
        for p in self.parameters: 
            if type(p) == Parameter:
                p._update()
        return None
    
    def _acceptNewValues(self):
        #Accept proposed value for parameters object 
        for p in self.parameters:
            if type(p) == Parameter:
                p._accept()
        return None
    
    def _updateSamples(self):
        newSamples = []
        for p in self.parameters:
            if type(p) == Parameter: newSamples.append(p.value)
            else: newSamples.append(p)
        self.samples.append(newSamples)
        return None
    

        
        