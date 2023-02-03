# specify the module that needs to be imported relative and its path
from importlib.machinery import SourceFileLoader

module_path = '/home/alessandro/Documenti/PhD/MESApaper/Maximum-Entropy-Spectrum/memspectrum.py'
module_name = 'memspectrum'
mem = SourceFileLoader(module_name,module_path).load_module() 

import numpy as np 
import matplotlib.pyplot as plt

class Whitening(object):
    
    def __init__(self, data, sampling_interval, data_duration):
        """
        Class the implements whitening of the data using different techniques        

        Parameters
        ----------
        data : np.ndarray
            The data to be whitened
        sampling_interval: np.float
            The sampling interval of the data
        data_duration: np.float
            The time duration of the data's observation 
        spectral_interval 
        
        reflection_coefficients : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.data = data
        self.N = len(data)
        self.sampling_interval = sampling_interval
        self.data_duration = data_duration 
    
    def whiten(self, spectrum_interval, ar_order, method = 'lattice', lag = None):
        self.ar_order = ar_order
        self.spectrum_interval = spectrum_interval
        self.lag = self._lag(method, lag)
        self.method = method #do I really need that? 
        
        if method.lower() == 'lattice':
            #Put here the computation of reflection coefficients? 
            return self._lattice()
        
        
        
                
        ####Il ciclo sui segmenti per il whitening conviene farlo con un while\
        ####Che si interrompe quando si supera la "fine" del array (T_max > T[-1])
        ###Si può usare l'indice per farlo
    def _create_data_segments(self): 
        ####Il modo più facile è storare gli indici 
        #1 
        samples_per_interval = self.lag / self.sampling_interval
        N_max = samples_per_interval 
        i_max = N - samples_per_interval 
        while N_max < N:
        step_0
        step_1 = spectrum_interval up to spectrum_interval + samples_per_interval 
        step_2 = 2 * spectrum_interval up to  * (2 * spectrum_interval+ samples per interval) 
        step_i = from i * spectrum_interval to i * spectrum_interval + nsamples
       
        
        #min lages
        return 0 
    def _lattice(self):
        
        return FPE
    
    
    def data_segments(self): 
        return 0 
    
    def _compute_reflection_coefficients(self):
        return 0 
        
    def prediction_errors(self, ar_order = 300):
        self.ar_order = ar_order
        FPE, BPE = np.zeros((ar_order + 1, self.N)), np.zeros((ar_order + 1, self.N))
        FPE[0], BPE[0] = self.data, self.data 
        for i in range(ar_order): 
            FPE[i+1,i+1:] = FPE[i, i+1:] + self.reflection_coefficients[i] * BPE[i,i:-1]
            BPE[i+1,i+1:] = BPE[i,i:-1] + self.reflection_coefficients[i] * FPE[i,i+1:]
            
        return FPE, BPE
   
    def _lag(self, method, lag):
        if method.lower() == 'lattice' and lag == None:
            seconds_lost = self.ar_order * self.sampling_interval
            lag = .5* (self.spectrum_interval - seconds_lost)
            
        elif type(lag) != float: 
            raise TypeError('methods different from "lattice" need lag to be\
                                 specified as a float number')
        elif lag > self.duration or lag <= 0:
            raise ValueError('Lag must be greater than 0 and lower than the \
                             total signal duration')
                             
        return lag

if __name__ == '__main__':
    import_ = False
    if import_: 
        data = np.loadtxt('/home/alessandro/Documenti/PhD/MESApaper/Maximum-Entropy-Spectrum/examples/On_source_toymodel/Noise_series/noise_realisations.txt')
        dt = 1/4096
    compute_FPE = False
    if compute_FPE: 
        data1 = data[0,:]
        m = mem.MESA() 
        P, ak, opt, k, P_array = m.solve(data1, early_stop = False)
        w = Whitening(data1, k)
        FPE, BPE = w.prediction_errors(k.size)
    i = 130
    P2, ak2, opt2, k, P_array = m.solve(FPE[i])
    f, sp = m.spectrum(dt)
    plt.loglog(f, sp)
    
    
